import tempfile
import tkinter as tk
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock
from tkinter import ttk

from ophiuchus.theme import COLORS, FONT_STACK
from ophiuchus.app import (
    OphiuchusApp,
    assess_vesta_preflight,
    build_library_cache_from_folder,
    build_output_dir,
    default_cache_path,
    default_candidate_dir,
    default_app_state_path,
    default_library_path,
    inspect_peak_from_app,
    initial_window_geometry,
    import_local_cifs_to_library,
    import_target_cif_to_library,
    load_database_api_config,
    load_app_state,
    library_manager_rows,
    library_target_phase_options,
    load_vesta_config,
    apply_vesta_env_config,
    materials_project_harvest_from_app,
    resolve_target_phase_selection,
    repair_app_state_paths,
    run_library_analysis_from_app,
    save_database_api_config,
    save_app_state,
    save_mp_api_key_to_env,
    save_formula_vesta_reference,
    save_vesta_config,
    set_library_entry_enabled,
    vesta_preflight_lines,
    validate_analysis_inputs,
    workbench_sections,
)
from ophiuchus.xrd.models import Candidate, CandidateScore
from ophiuchus.session_storage import TransientAnalysisStore


class AppExperienceTests(unittest.TestCase):
    def test_initial_window_geometry_fits_common_laptop_work_area(self):
        geometry = initial_window_geometry(1536, 864)
        width, height = (int(value) for value in geometry.split("+", 1)[0].split("x"))

        self.assertEqual(width, 1240)
        self.assertLessEqual(height, 700)
        self.assertGreaterEqual(height, 620)
        self.assertRegex(geometry, r"^\d+x\d+\+\d+\+\d+$")

    def test_main_workbench_keeps_sidebar_compact_and_result_tabs_visible(self):
        app = OphiuchusApp()
        try:
            app.geometry("1240x673+0+0")
            app.update()

            self.assertLessEqual(app.sidebar.winfo_width(), 155)
            self.assertGreaterEqual(app.tabs.winfo_width(), 500)
            visible_tabs = set()
            for x in range(0, app.tabs.winfo_width(), 4):
                try:
                    visible_tabs.add(app.tabs.index(f"@{x},10"))
                except tk.TclError:
                    pass
            self.assertIn(len(app.tabs.tabs()) - 1, visible_tabs)
        finally:
            app.destroy()

    def test_advanced_inputs_are_collapsed_until_requested(self):
        app = OphiuchusApp()
        try:
            app.update()
            self.assertFalse(app.advanced_inputs.winfo_ismapped())
            self.assertEqual(app.advanced_toggle_var.get(), "显示数据与保存设置")

            app.advanced_toggle_button.invoke()
            app.update()

            self.assertTrue(app.advanced_inputs.winfo_ismapped())
            self.assertEqual(app.advanced_toggle_var.get(), "收起数据与保存设置")
        finally:
            app.destroy()

    def test_core_inputs_fit_without_scrolling_at_default_height(self):
        app = OphiuchusApp()
        try:
            app.geometry("1240x673+0+0")
            app.update()
            scrollregion = [int(float(value)) for value in app.workflow_canvas.cget("scrollregion").split()]
            content_height = scrollregion[3] - scrollregion[1]

            self.assertLessEqual(content_height, app.workflow_canvas.winfo_height())
        finally:
            app.destroy()

    def test_target_phase_toolbar_keeps_both_commands_visible(self):
        app = OphiuchusApp()
        try:
            app.geometry("1240x673+0+0")
            app.update()
            row = app.target_phase_combo.master

            self.assertGreaterEqual(app.refresh_target_button.winfo_width(), 50)
            self.assertGreaterEqual(app.import_target_button.winfo_width(), 65)
            self.assertLessEqual(
                app.import_target_button.winfo_x() + app.import_target_button.winfo_width(),
                row.winfo_width(),
            )
        finally:
            app.destroy()

    def test_main_commands_are_grouped_by_workflow_stage(self):
        app = OphiuchusApp()
        try:
            self.assertEqual(int(app.run_button.grid_info()["row"]), 0)
            self.assertEqual(int(app.library_run_button.grid_info()["row"]), 0)
            self.assertIs(app.database_button.master, app.setup_actions)
            self.assertIs(app.vesta_button.master, app.setup_actions)
            self.assertIs(app.save_analysis_button.master, app.result_actions)
            self.assertIs(app.phase_stripping_button.master, app.result_actions)
            self.assertIs(app.refinement_button.master, app.result_actions)
        finally:
            app.destroy()

    def test_app_destroy_cancels_pending_result_poll(self):
        app = OphiuchusApp()
        callback_id = app._poll_after_id

        app.destroy()

        self.assertNotIn(callback_id, app.tk.call("after", "info"))

    def test_database_dialog_exposes_background_common_oxide_supplement(self):
        app = OphiuchusApp()
        app.elements_var.set("Zr Fe Ge")
        try:
            with mock.patch("ophiuchus.app.threading.Thread") as thread:
                app._open_database_dialog()
                app.update_idletasks()
                dialogs = [item for item in app.winfo_children() if isinstance(item, tk.Toplevel)]
                dialog = next(item for item in dialogs if item.title() == "连接 API / 数据库")

                def descendants(widget):
                    rows = list(widget.winfo_children())
                    for child in list(rows):
                        rows.extend(descendants(child))
                    return rows

                buttons = [item for item in descendants(dialog) if isinstance(item, ttk.Button)]
                supplement = next(item for item in buttons if str(item.cget("text")) == "补齐常见氧化物")
                supplement.invoke()

            self.assertTrue(app.running)
            thread.assert_called_once()
            self.assertIn("氧化物", app.status_var.get())
        finally:
            for item in app.winfo_children():
                if isinstance(item, tk.Toplevel):
                    item.destroy()
            app.destroy()

    def test_oxide_supplement_worker_uses_current_main_elements_and_reports_queue_result(self):
        app = OphiuchusApp()
        app.elements_var.set("Fe Ge")
        app.library_var.set("library.sqlite")
        try:
            with mock.patch(
                "ophiuchus.app.supplement_common_oxide_library",
                return_value={"imported": 1, "requested_formulas": ["FeO"], "missing_after": []},
            ) as supplement:
                app._oxide_supplement_worker()
            kind, payload = app.result_queue.get_nowait()

            self.assertEqual(kind, "oxide_supplement_ok")
            self.assertEqual(payload["imported"], 1)
            self.assertEqual(supplement.call_args.args[1], {"Fe", "Ge"})
        finally:
            app.destroy()
    def test_default_candidate_dir_prefers_existing_categorized_desktop_structure_folder(self):
        with tempfile.TemporaryDirectory() as tmp:
            desktop = Path(tmp) / "Desktop"
            expected = desktop / "03_实验数据与分析" / "结构与CIF" / "结构"
            expected.mkdir(parents=True)
            with mock.patch("ophiuchus.app.desktop_dir", return_value=desktop):
                self.assertEqual(default_candidate_dir(Path(tmp) / "project"), expected)

    def test_repair_app_state_paths_remaps_moved_project_files_and_clears_missing_xrd(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "Ophiuchus"
            (root / "data").mkdir(parents=True)
            (root / "results" / "980").mkdir(parents=True)
            old_root = Path(tmp) / "old" / "Ophiuchus" / "Ophiuchus"
            state = {
                "xrd_file": str(Path(tmp) / "missing.asc"),
                "out_dir": str(old_root / "results" / "980"),
                "cache_path": str(old_root / "data" / "ophi_xrd_cache.sqlite"),
                "library_path": str(old_root / "data" / "ophi_library.sqlite"),
                "vesta_exe": str(old_root / "tools" / "VESTA.exe"),
                "rietan_exe": str(old_root / "tools" / "RIETAN.exe"),
            }

            repaired = repair_app_state_paths(state, root=root)

            self.assertEqual(repaired["xrd_file"], "")
            self.assertEqual(Path(repaired["out_dir"]), root / "results" / "980")
            self.assertEqual(Path(repaired["cache_path"]), root / "data" / "ophi_xrd_cache.sqlite")
            self.assertEqual(Path(repaired["library_path"]), root / "data" / "ophi_library.sqlite")
            self.assertNotIn("vesta_exe", repaired)
            self.assertNotIn("rietan_exe", repaired)

    def test_repair_app_state_paths_does_not_duplicate_results_folder_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "Ophiuchus"
            (root / "results").mkdir(parents=True)
            state = {"out_dir": str(Path(tmp) / "old" / "Ophiuchus" / "results")}

            repaired = repair_app_state_paths(state, root=root)

            self.assertEqual(Path(repaired["out_dir"]), root / "results")

    def test_build_output_dir_uses_xrd_stem_under_project_results(self):
        out = build_output_dir(Path(r"C:\data\ZrFe6Ge4 900.asc"), Path(r"C:\project"))
        self.assertEqual(out, Path(r"C:\project\results\ZrFe6Ge4_900"))

    def test_default_cache_path_lives_under_project_data(self):
        cache = default_cache_path(Path(r"C:\project"))
        self.assertEqual(cache, Path(r"C:\project\data\ophi_xrd_cache.sqlite"))

    def test_default_library_path_lives_under_project_data(self):
        library = default_library_path(Path(r"C:\project"))
        self.assertEqual(library, Path(r"C:\project\data\ophi_library.sqlite"))

    def test_default_app_state_path_lives_under_project_data(self):
        state = default_app_state_path(Path(r"C:\project"))
        self.assertEqual(state, Path(r"C:\project\data\ophi_app_state.json"))

    def test_load_app_state_accepts_windows_utf8_bom(self):
        with tempfile.TemporaryDirectory() as tmp:
            state_path = Path(tmp) / "state.json"
            state_path.write_bytes(b"\xef\xbb\xbf{\"elements\": \"Zr Fe Ge\"}")

            state = load_app_state(state_path)

        self.assertEqual(state["elements"], "Zr Fe Ge")

    def test_workbench_sections_match_phase2_navigation(self):
        labels = [section["label"] for section in workbench_sections()]
        self.assertEqual(labels, ["项目总览", "样品与输入", "结构库", "XRD 分析", "物相证据", "设置"])

    def test_theme_supports_light_chinese_ui(self):
        self.assertEqual(COLORS["background"], "#eef3f8")
        self.assertIn("Microsoft YaHei UI", FONT_STACK)
        self.assertIn("PingFang SC", FONT_STACK)
        self.assertNotEqual(COLORS["background"], "#0b1020")

    def test_save_mp_key_and_dry_run_harvest_preview(self):
        with tempfile.TemporaryDirectory() as tmp:
            env_path = Path(tmp) / ".env"
            save_mp_api_key_to_env(env_path, "abc123")
            text = env_path.read_text(encoding="utf-8")
            preview = materials_project_harvest_from_app(
                library_db=Path(tmp) / "library.sqlite",
                elements="Zr Fe Ge",
                impurities="O",
                mode="normal",
                api_key="",
                env_path=env_path,
                dry_run=True,
            )
        self.assertIn("MP_API_KEY=abc123", text)
        self.assertIn("Fe-Ge-Zr", preview["chemical_systems"])
        self.assertIn("Zr-O", preview["chemical_systems"])

    def test_database_api_config_saves_mp_provider_without_exposing_key(self):
        with tempfile.TemporaryDirectory() as tmp:
            env_path = Path(tmp) / ".env"
            saved = save_database_api_config(env_path, "Materials Project", "secret-key")
            text = env_path.read_text(encoding="utf-8")
        self.assertTrue(saved["configured"])
        self.assertEqual(saved["provider"], "materials_project")
        self.assertIn("OPHI_STRUCTURE_DATABASE=materials_project", text)
        self.assertIn("MP_API_KEY=secret-key", text)

    def test_database_api_config_loads_saved_key_and_does_not_overwrite_with_blank(self):
        with tempfile.TemporaryDirectory() as tmp:
            env_path = Path(tmp) / ".env"
            save_database_api_config(env_path, "Materials Project", "secret-key", endpoint="https://api.materialsproject.org")
            save_database_api_config(env_path, "Materials Project", "", endpoint="")
            loaded = load_database_api_config(env_path)
        self.assertEqual(loaded["provider"], "materials_project")
        self.assertEqual(loaded["api_key"], "secret-key")
        self.assertEqual(loaded["endpoint"], "https://api.materialsproject.org")

    def test_app_state_round_trips_last_main_screen_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            state_path = Path(tmp) / "state.json"
            save_app_state(
                state_path,
                {
                    "xrd_file": r"C:\data\last.asc",
                    "candidate_dir": r"C:\data\cifs",
                    "elements": "Zr V Ge Sn",
                    "extra_elements": "",
                    "target_phase_label": "Zr3V3GeSn4 | local | abc",
                    "mp_api_key": "must-not-be-saved-here",
                },
            )
            loaded = load_app_state(state_path)
        self.assertEqual(loaded["xrd_file"], r"C:\data\last.asc")
        self.assertEqual(loaded["elements"], "Zr V Ge Sn")
        self.assertEqual(loaded["extra_elements"], "")
        self.assertNotIn("mp_api_key", loaded)

    def test_app_constructs_three_column_workbench(self):
        from ophiuchus.app import OphiuchusApp

        app = OphiuchusApp()
        try:
            self.assertTrue(hasattr(app, "sidebar"))
            self.assertTrue(hasattr(app, "workspace_panel"))
            self.assertTrue(hasattr(app, "workflow_canvas"))
            self.assertTrue(hasattr(app, "workflow_scrollbar"))
            self.assertTrue(hasattr(app, "workflow_content"))
            self.assertTrue(hasattr(app, "tabs"))
            self.assertIn("library", app.nav_buttons)
            self.assertTrue(hasattr(app, "database_provider_var"))
            self.assertTrue(hasattr(app, "vesta_button"))
            self.assertTrue(hasattr(app, "rietan_exe_var"))
            self.assertTrue(hasattr(app, "target_phase_combo"))
            self.assertTrue(hasattr(app, "periodic_table_button"))
            self.assertGreaterEqual(int(app.target_phase_combo.cget("height")), 8)
            self.assertTrue(hasattr(app, "save_plot_button"))
            self.assertTrue(hasattr(app, "save_analysis_button"))
            self.assertTrue(hasattr(app, "phase_stripping_button"))
            self.assertTrue(hasattr(app, "refinement_button"))
            self.assertTrue(app.save_analysis_button.instate(["disabled"]))
            self.assertTrue(app.phase_stripping_button.instate(["disabled"]))
            self.assertTrue(app.refinement_button.instate(["disabled"]))
            self.assertTrue(hasattr(app, "candidate_phase_text"))
            self.assertTrue(hasattr(app, "used_structures_text"))
            app._select_section("library")
            self.assertEqual(app.status_var.get(), "当前工作区：结构库")
        finally:
            app.destroy()

    def test_sidebar_help_button_opens_about_dialog_with_contact_actions(self):
        app = OphiuchusApp()
        try:
            app.help_button.invoke()
            app.update()

            dialog = next(
                child
                for child in app.winfo_children()
                if isinstance(child, tk.Toplevel) and child.title() == "帮助与关于"
            )

            def descendants(widget):
                children = list(widget.winfo_children())
                for child in list(children):
                    children.extend(descendants(child))
                return children

            widgets = descendants(dialog)
            label_text = "\n".join(
                str(widget.cget("text")) for widget in widgets if isinstance(widget, ttk.Label)
            )
            buttons = {
                str(widget.cget("text")): widget for widget in widgets if isinstance(widget, ttk.Button)
            }

            self.assertIn("wanyc@issp.u-tokyo.ac.jp", label_text)
            self.assertIn("打开操作手册", buttons)
            self.assertIn("复制邮箱", buttons)
            self.assertIn("发送邮件", buttons)

            buttons["复制邮箱"].invoke()
            self.assertEqual(app.clipboard_get(), "wanyc@issp.u-tokyo.ac.jp")
        finally:
            for child in app.winfo_children():
                if isinstance(child, tk.Toplevel):
                    child.destroy()
            app.destroy()

    def test_main_analyses_wait_until_open_refinement_session_is_closed(self):
        app = OphiuchusApp()
        try:
            with (
                mock.patch.object(app, "_refinement_is_open", return_value=True),
                mock.patch("ophiuchus.app.messagebox.showwarning") as warning,
            ):
                app._run()
                app._run_library()

            self.assertEqual(warning.call_count, 2)
            self.assertIn("精修", warning.call_args_list[0].args[0])
            self.assertIn("精修", warning.call_args_list[1].args[0])
            self.assertFalse(app.running)
        finally:
            app.destroy()

    def test_vesta_dialog_stays_non_modal_visible_and_reuses_existing_window(self):
        app = OphiuchusApp()
        app.update()
        try:
            app.vesta_button.invoke()
            app.update()
            first = app._vesta_dialog

            self.assertIsNotNone(first)
            self.assertTrue(first.winfo_viewable())
            self.assertIsNone(app.grab_current())
            self.assertNotEqual(first.geometry().split("+", 1)[1], "0+0")
            self.assertEqual(str(app.vesta_button.cget("text")), "VESTA")
            self.assertEqual(first.title(), "VESTA / RIETAN 设置")
            self.assertTrue(bool(first.attributes("-topmost")))

            app.vesta_button.invoke()
            app.update()
            self.assertIs(app._vesta_dialog, first)
        finally:
            if getattr(app, "_vesta_dialog", None) is not None:
                app._vesta_dialog.destroy()
            app.destroy()

    def test_vesta_config_saves_local_executable_and_reference_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe = root / "VESTA.exe"
            rietan = root / "RIETAN.exe"
            refs = root / "refs"
            exe.write_text("", encoding="utf-8")
            rietan.write_text("", encoding="utf-8")
            refs.mkdir()
            env_path = root / ".env"
            saved = save_vesta_config(env_path, str(exe), str(refs), str(rietan))
            loaded = load_vesta_config(env_path)
        self.assertTrue(saved["vesta_exe_exists"])
        self.assertTrue(saved["rietan_exe_exists"])
        self.assertTrue(saved["reference_dir_exists"])
        self.assertEqual(loaded["vesta_exe"], str(exe))
        self.assertEqual(loaded["rietan_exe"], str(rietan))
        self.assertEqual(loaded["reference_dir"], str(refs))

    def test_vesta_config_ignores_reference_directory_from_another_computer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fallback = root / "XRD"
            fallback.mkdir()
            env_path = root / ".env"
            env_path.write_text(
                "OPHI_VESTA_REFERENCE_DIR=C:\\Users\\someone_else\\Desktop\\XRD\n",
                encoding="utf-8",
            )

            with mock.patch("ophiuchus.app.default_xrd_dir", return_value=fallback):
                loaded = load_vesta_config(env_path)

        self.assertEqual(loaded["reference_dir"], str(fallback))

    def test_vesta_config_ignores_missing_executable_from_another_computer(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            discovered = root / "VESTA.exe"
            discovered.write_text("", encoding="utf-8")
            env_path = root / ".env"
            env_path.write_text(
                "OPHI_VESTA_EXE=C:\\Users\\someone_else\\Desktop\\VESTA.exe\n",
                encoding="utf-8",
            )

            with mock.patch("ophiuchus.app.find_vesta_executable", return_value=str(discovered)):
                loaded = load_vesta_config(env_path)

        self.assertEqual(loaded["vesta_exe"], str(discovered))

    def test_simulation_validation_summary_lines_show_failed_vesta_candidates(self):
        failed = Candidate("bad", "ZrFe6Ge4", "library:local", "", ["Zr", "Fe", "Ge"], None)
        failed.simulation_validation = {"status": "failed"}
        passed = Candidate("good", "Fe", "library:local", "", ["Fe"], None)
        passed.simulation_validation = {"status": "passed"}
        result = type("Result", (), {"top_scores": [CandidateScore(failed, 0.2, [], [], [], []), CandidateScore(passed, 0.8, [], [], [], [])]})()
        lines = OphiuchusApp._simulation_validation_summary_lines(None, result)
        joined = "\n".join(lines)
        self.assertIn("failed: 1", joined)
        self.assertIn("passed: 1", joined)
        self.assertIn("ZrFe6Ge4", joined)

    def test_vesta_preflight_reports_available_and_missing_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe = root / "VESTA.exe"
            refs = root / "refs"
            exe.write_text("", encoding="utf-8")
            refs.mkdir()
            (refs / "ZrFe6Ge4 VESTA.int").write_text("20 100\n30 80\n", encoding="utf-8")
            status = assess_vesta_preflight(str(exe), str(refs), ["ZrFe6Ge4", "Zr3V3GeSn4"])
        self.assertTrue(status["vesta_exe_exists"])
        self.assertTrue(status["reference_dir_exists"])
        self.assertEqual(status["formulas_checked"], 2)
        self.assertEqual(status["references_found"], 1)
        self.assertEqual(status["missing_references"], ["Zr3V3GeSn4"])
        self.assertIn("ZrFe6Ge4", status["reference_paths"])
        self.assertEqual(status["reference_candidate_counts"]["ZrFe6Ge4"], 1)

    def test_save_formula_vesta_reference_writes_stable_env_key(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / "Zr3V3GeSn4 VESTA.int"
            ref.write_text("20 100\n", encoding="utf-8")
            env_path = root / ".env"
            saved = save_formula_vesta_reference(env_path, "Zr3V3GeSn4", str(ref))
            values = env_path.read_text(encoding="utf-8")
        self.assertTrue(saved["reference_exists"])
        self.assertIn("OPHI_VESTA_REFERENCE_ZR3V3GESN4=", values)

    def test_apply_vesta_env_config_restores_formula_overrides_for_new_session(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / "Zr3V3GeSn4 VESTA.int"
            ref.write_text("20 100\n", encoding="utf-8")
            env_path = root / ".env"
            save_formula_vesta_reference(env_path, "Zr3V3GeSn4", str(ref))
            import os

            os.environ.pop("OPHI_VESTA_REFERENCE_ZR3V3GESN4", None)
            apply_vesta_env_config(env_path)
            restored = os.environ.get("OPHI_VESTA_REFERENCE_ZR3V3GESN4")
        self.assertEqual(restored, str(ref))

    def test_vesta_preflight_lines_are_actionable(self):
        lines = vesta_preflight_lines(
            {
                "vesta_exe_exists": True,
                "reference_dir_exists": True,
                "formulas_checked": 2,
                "references_found": 1,
                "missing_references": ["Zr3V3GeSn4"],
                "reference_paths": {"ZrFe6Ge4": "ref.int"},
                "reference_candidate_counts": {"ZrFe6Ge4": 2},
                "reference_candidates": {"ZrFe6Ge4": ["ref_a.int", "ref_b.int"]},
            }
        )
        joined = "\n".join(lines)
        self.assertIn("VESTA 预检", joined)
        self.assertIn("1/2", joined)
        self.assertIn("Zr3V3GeSn4", joined)
        self.assertIn("备选", joined)
        self.assertIn("ref_a.int", joined)

    def test_import_and_cache_helpers_support_local_cifs(self):
        cif_text = """
data_Fe
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Fe1 1.0 0.0 0.0 0.0 Fe
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "Fe.cif").write_text(cif_text, encoding="utf-8")
            db = root / "library.sqlite"
            imported = import_local_cifs_to_library(cifs, db)
            cached = build_library_cache_from_folder(db)
        self.assertEqual(imported["imported"], 1)
        self.assertEqual(cached["simulated"], 1)

    def test_import_target_cif_imports_one_file_and_returns_target_label(self):
        cif_text = """
data_ZrVGeSn
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Zr1 1.0 0.0 0.0 0.0 Zr
V1 1.0 0.5 0.5 0.5 V
Ge1 1.0 0.25 0.25 0.25 Ge
Sn1 1.0 0.75 0.75 0.75 Sn
"""
        other_cif = """
data_Sn
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Sn1 1.0 0.0 0.0 0.0 Sn
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            target = cifs / "ZrVGeSn.cif"
            target.write_text(cif_text, encoding="utf-8")
            (cifs / "Sn.cif").write_text(other_cif, encoding="utf-8")
            db = root / "library.sqlite"
            imported = import_target_cif_to_library(target, db)
            rows = library_target_phase_options(db, "Zr V Ge Sn")
        self.assertEqual(imported["imported"], 1)
        self.assertIn("ZrVGeSn", imported["label"])
        self.assertEqual([row["label"] for row in rows], [imported["label"]])

    def test_library_manager_rows_and_enabled_toggle(self):
        cif_text = """
data_Fe
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Fe1 1.0 0.0 0.0 0.0 Fe
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "Fe.cif").write_text(cif_text, encoding="utf-8")
            db = root / "library.sqlite"
            import_local_cifs_to_library(cifs, db)
            build_library_cache_from_folder(db)
            rows = library_manager_rows(db)
            set_library_entry_enabled(db, rows[0]["internal_id"], False)
            disabled_rows = library_manager_rows(db)
        self.assertEqual(rows[0]["formula"], "Fe")
        self.assertEqual(rows[0]["xrd_cache"], "validated")
        self.assertEqual(rows[0]["backend_name"], "Ophi Validated pymatgen XRD")
        self.assertTrue(rows[0]["enabled"])
        self.assertFalse(disabled_rows[0]["enabled"])

    def test_target_phase_options_prefer_local_exact_main_element_system(self):
        main_cif = """
data_ZrFeGe
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Zr1 1.0 0.0 0.0 0.0 Zr
Fe1 1.0 0.5 0.5 0.5 Fe
Ge1 1.0 0.25 0.25 0.25 Ge
"""
        fe_cif = """
data_Fe
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Fe1 1.0 0.0 0.0 0.0 Fe
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "ZrFeGe.cif").write_text(main_cif, encoding="utf-8")
            (cifs / "Fe.cif").write_text(fe_cif, encoding="utf-8")
            db = root / "library.sqlite"
            import_local_cifs_to_library(cifs, db)
            options = library_target_phase_options(db, "Zr Fe Ge")
        self.assertIn("ZrFeGe", options[0]["label"])

    def test_target_phase_options_exclude_subsystem_impurities(self):
        target_cif = """
data_ZrVGeSn
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Zr1 1.0 0.0 0.0 0.0 Zr
V1 1.0 0.5 0.5 0.5 V
Ge1 1.0 0.25 0.25 0.25 Ge
Sn1 1.0 0.75 0.75 0.75 Sn
"""
        sn_cif = """
data_Sn
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Sn1 1.0 0.0 0.0 0.0 Sn
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "ZrVGeSn.cif").write_text(target_cif, encoding="utf-8")
            (cifs / "Sn.cif").write_text(sn_cif, encoding="utf-8")
            db = root / "library.sqlite"
            import_local_cifs_to_library(cifs, db)
            options = library_target_phase_options(db, "Zr V Ge Sn")
        labels = [row["label"] for row in options]
        self.assertTrue(any("ZrVGeSn" in label for label in labels))
        self.assertFalse(any(label.startswith("Sn |") for label in labels))

    def test_target_phase_selection_drops_stale_old_system_label(self):
        rows = [
            {"label": "Zr3V3GeSn4 | local | new", "id": "new"},
            {"label": "SnGe | materials_project | impurity", "id": "impurity"},
        ]
        label, internal_id, options = resolve_target_phase_selection(rows, "ZrFe6Ge4 | local | old")
        self.assertEqual(label, "")
        self.assertEqual(internal_id, "")
        self.assertNotIn("ZrFe6Ge4 | local | old", options)

    def test_inspect_peak_from_app_uses_library_cache(self):
        cif_text = """
data_Fe
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Fe1 1.0 0.0 0.0 0.0 Fe
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "Fe.cif").write_text(cif_text, encoding="utf-8")
            db = root / "library.sqlite"
            import_local_cifs_to_library(cifs, db)
            build_library_cache_from_folder(db)
            evidence = inspect_peak_from_app(db, 22.2, tolerance_deg=0.3)
        self.assertEqual(evidence["status"], "matched")
        self.assertEqual(evidence["nearby_peaks"][0]["formula"], "Fe")

    def test_run_library_analysis_helper_supports_cached_library(self):
        cif_text = """
data_Fe
_cell_length_a 4.000
_cell_length_b 4.000
_cell_length_c 4.000
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_occupancy
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_type_symbol
Fe1 1.0 0.0 0.0 0.0 Fe
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cifs = root / "cifs"
            cifs.mkdir()
            (cifs / "Fe.cif").write_text(cif_text, encoding="utf-8")
            xrd = root / "sample.xy"
            rows = []
            for i in range(2000):
                x = 10.0 + i * 0.03
                y = 1.0 + 100.0 * pow(2.718281828, -((x - 22.2) ** 2) / (2 * 0.06**2))
                rows.append(f"{x:.3f} {y:.5f}")
            xrd.write_text("\n".join(rows), encoding="utf-8")
            db = root / "library.sqlite"
            import_local_cifs_to_library(cifs, db)
            result = run_library_analysis_from_app(
                xrd, db, "Fe", "", root / "out",
                scientific_safe_mode=False,
                force_recompute=True,
            )
        self.assertTrue(result.top_scores)
        self.assertEqual(result.top_scores[0].candidate.source, "library:local")
        self.assertFalse(result.scientific_runtime["scientific_safe_mode"])
        self.assertTrue(result.scientific_runtime["force_recompute"])

    def test_validate_inputs_reports_actionable_missing_values(self):
        errors = validate_analysis_inputs("", "", "", "")
        joined = "\n".join(errors)
        self.assertIn("请选择实验 XRD 文件", joined)
        self.assertIn("请选择候选结构/模拟峰文件夹", joined)
        self.assertIn("请填写主元素", joined)
        self.assertNotIn("请选择输出文件夹", joined)

    def test_validate_inputs_checks_real_paths_and_candidate_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            xrd = root / "sample.asc"
            xrd.write_text("*START = 10\n*STEP = 1\n*COUNT = 1\n1", encoding="utf-8")
            empty = root / "empty"
            empty.mkdir()
            errors = validate_analysis_inputs(str(xrd), str(empty), "Zr Fe Ge", str(root / "out"))
        self.assertIn("候选文件夹里没有找到 CIF 或峰表文件", "\n".join(errors))

    def test_validate_inputs_blocks_elements_missing_from_inferred_formula_but_allows_extras(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            xrd = root / "Zr3V3GeSn4.asc"
            xrd.write_text("10 1\n11 2", encoding="utf-8")
            candidates = root / "candidates"
            candidates.mkdir()
            (candidates / "phase.cif").write_text("data_test", encoding="utf-8")

            missing = validate_analysis_inputs(str(xrd), str(candidates), "Zr Ge Sn", str(root / "out"))
            with_extra = validate_analysis_inputs(str(xrd), str(candidates), "Zr V Ge Sn O", str(root / "out"))

        self.assertIn("V", "\n".join(missing))
        self.assertFalse(any("实验谱名称推断" in error for error in with_extra))

    def test_gui_worker_uses_single_transient_session_instead_of_user_results_folder(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            permanent = root / "results"
            app = OphiuchusApp()
            try:
                app.analysis_store = TransientAnalysisStore(root / "analysis-session")
                app._analysis_run_id = 1
                snapshot = {
                    "xrd_file": str(root / "sample.xy"),
                    "candidate_dir": str(root / "candidates"),
                    "elements": ("Fe",),
                    "extra_elements": (),
                    "cache_path": "",
                }
                called = {}

                def fake_run_analysis(**kwargs):
                    called.update(kwargs)
                    out = Path(kwargs["out_dir"])
                    out.mkdir(parents=True, exist_ok=True)
                    report = out / "results.json"
                    report.write_text("{}", encoding="utf-8")
                    return SimpleNamespace(outputs={"json": str(report)})

                with mock.patch("ophiuchus.app.run_analysis", side_effect=fake_run_analysis):
                    app._worker(1, snapshot)
                kind, payload = app.result_queue.get_nowait()

                self.assertEqual(kind, "analysis_ok")
                self.assertEqual(payload["run_id"], 1)
                self.assertTrue(Path(payload["result"].outputs["json"]).is_file())
                self.assertEqual(Path(called["out_dir"]), app.analysis_store.pending_path)
                self.assertTrue(called["use_rietan_display"])
                self.assertFalse(permanent.exists())
            finally:
                app.destroy()


if __name__ == "__main__":
    unittest.main()
