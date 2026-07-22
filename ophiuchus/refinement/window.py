from __future__ import annotations

import csv
from dataclasses import asdict
import json
from pathlib import Path
import queue
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Callable, Iterable

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from ophiuchus.phase_stripping.models import AnalysisContext
from ophiuchus.refinement.conclusion import (
    SampleConclusion,
    build_sample_conclusion,
    discover_sibling_peak_references,
)
from ophiuchus.refinement.oxide_candidates import load_controlled_oxide_candidates
from ophiuchus.refinement.evidence_window import CandidateEvidenceDialog
from ophiuchus.theme import COLORS, SPACING, configure_matplotlib_fonts
from ophiuchus.xrd.models import Candidate
from ophiuchus.xrd.multiphase_models import (
    MultiphaseRefinementResult,
    MultiphaseRefinementSettings,
    PhaseRefinementInput,
)
from ophiuchus.xrd.multiphase_refinement import RietanMultiphaseBackend
from ophiuchus.xrd.quantification import (
    QuantificationAssessment,
    assess_quantification,
    initial_scale_variants,
)
from ophiuchus.xrd.refinement import (
    RefinementSettings,
    RietanRefinementBackend,
    RietanRefinementResult,
)


MODE_SETTINGS = {
    "背景 + 尺度": (False, False),
    "加零点偏移": (True, False),
    "加峰宽 U/V/W": (True, True),
}


class SampleConclusionDialog(tk.Toplevel):
    def __init__(
        self,
        parent: tk.Misc,
        conclusion: SampleConclusion,
        *,
        on_exclude_formula: Callable[[str], None] | None = None,
    ) -> None:
        super().__init__(parent)
        self.title("Ophi 样品结论")
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        width = min(820, screen_width - 100)
        height = min(720, screen_height - 100)
        left = max(0, (screen_width - width) // 2)
        top = max(0, (screen_height - height) // 2)
        self.geometry(f"{width}x{height}+{left}+{top}")
        self.minsize(620, 520)
        self.configure(bg=COLORS["background"])
        self.transient(parent)
        self.conclusion = conclusion
        self._on_exclude_formula = on_exclude_formula
        self.evidence_dialog: CandidateEvidenceDialog | None = None

        shell = ttk.Frame(self, style="Shell.TFrame", padding=SPACING["window_pad"])
        shell.pack(fill="both", expand=True)
        ttk.Label(shell, text="样品结论", style="Title.TLabel").pack(anchor="w", pady=(0, 10))
        text_shell = ttk.Frame(shell, style="Panel.TFrame", padding=1)
        text_shell.pack(fill="both", expand=True)
        self.report_text = tk.Text(
            text_shell,
            wrap="word",
            relief="flat",
            borderwidth=0,
            padx=16,
            pady=14,
            bg=COLORS["panel"],
            fg=COLORS["text"],
            font=("Microsoft YaHei UI", 10),
            spacing1=2,
            spacing3=5,
        )
        scrollbar = ttk.Scrollbar(text_shell, orient="vertical", command=self.report_text.yview)
        self.report_text.configure(yscrollcommand=scrollbar.set)
        self.report_text.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.report_text.insert("1.0", conclusion.to_chinese_text())
        self.report_text.configure(state="disabled")
        actions = ttk.Frame(shell, style="Shell.TFrame")
        actions.pack(fill="x", pady=(10, 0))
        self.evidence_button = ttk.Button(actions, text="候选证据", command=self._open_evidence)
        self.evidence_button.pack(side="left")
        if not conclusion.impurity_estimates:
            self.evidence_button.state(["disabled"])
        self.copy_button = ttk.Button(actions, text="复制结论", command=self._copy_report)
        self.copy_button.pack(side="left", padx=(8, 0))
        ttk.Button(actions, text="关闭", command=self.destroy).pack(side="right")
        self.lift()
        self.focus_force()

    def _copy_report(self) -> None:
        self.clipboard_clear()
        self.clipboard_append(self.conclusion.to_chinese_text())
        self.update_idletasks()
        self.copy_button.configure(text="已复制")

    def _open_evidence(self) -> None:
        if self.evidence_dialog is not None:
            try:
                if self.evidence_dialog.winfo_exists():
                    self.evidence_dialog.lift()
                    return
            except tk.TclError:
                pass
        self.evidence_dialog = CandidateEvidenceDialog(
            self,
            self.conclusion,
            on_exclude_formula=self._on_exclude_formula,
        )


class RefinementWindow(tk.Toplevel):
    def __init__(
        self,
        parent: tk.Misc,
        context: AnalysisContext,
        candidates: Iterable[Candidate],
        *,
        backend: RietanRefinementBackend | None = None,
        multiphase_backend: RietanMultiphaseBackend | None = None,
        supporting_candidates: Iterable[Candidate] | None = None,
        oxide_library_path: str | Path | None = None,
        on_closed: Callable[[], None] | None = None,
    ) -> None:
        super().__init__(parent)
        self.title("RIETAN-FP 受约束精修")
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        width = min(1440, max(1100, screen_width - 120))
        height = min(860, max(700, screen_height - 120))
        left = max(0, (screen_width - width) // 2)
        top = max(0, (screen_height - height) // 2)
        self.geometry(f"{width}x{height}+{left}+{top}")
        self.minsize(980, 640)
        self.configure(bg=COLORS["background"])
        self.context = context
        supplied_candidates = list(candidates)
        self.candidates = [
            candidate for candidate in supplied_candidates if Path(candidate.source_path).suffix.lower() == ".cif"
        ]
        self.supporting_candidates = list(supporting_candidates or supplied_candidates)
        self.oxide_library_path = str(oxide_library_path).strip() if oxide_library_path else ""
        self.backend = backend or RietanRefinementBackend()
        self.multiphase_backend = multiphase_backend or RietanMultiphaseBackend()
        self._on_closed = on_closed
        self._queue: queue.Queue[tuple[str, object]] = queue.Queue()
        self._running = False
        self._closed = False
        self._poll_after_id: str | None = None
        self.result: RietanRefinementResult | None = None
        self.multiphase_result: MultiphaseRefinementResult | None = None
        self.quantification_assessment: QuantificationAssessment | None = None
        self.stability_results: tuple[MultiphaseRefinementResult, ...] = ()
        self.conclusion: SampleConclusion | None = None
        self.conclusion_dialog: SampleConclusionDialog | None = None
        self._active_candidate: Candidate | None = None
        self._conclusion_candidates: list[Candidate] = []
        self._oxide_formulas: tuple[str, ...] = ()
        self._oxide_warnings: tuple[str, ...] = ()
        self._excluded_formulas: set[str] = set()

        self._candidate_by_label = {
            f"{candidate.formula_pretty} | {candidate.candidate_id}": candidate for candidate in self.candidates
        }
        self.target_var = tk.StringVar(value=next(iter(self._candidate_by_label), ""))
        self.workflow_var = tk.StringVar(value="单相确认")
        self.impurity_vars: dict[str, tk.BooleanVar] = {}
        self.mode_var = tk.StringVar(value="背景 + 尺度")
        self.background_terms_var = tk.IntVar(value=6)
        self.refine_lattice_var = tk.BooleanVar(value=False)
        self.status_var = tk.StringVar(value="选择目标相和精修层级后开始；原始实验谱不会被覆盖。")
        self.workflow_note_var = tk.StringVar(value="单相确认，不输出 wt%")
        self.metrics_var = tk.StringVar(value="尚未运行")
        self.locked_parameters_var = tk.StringVar(
            value="锁定：原子坐标、占位率、热参数、择优取向、各向异性展宽"
        )
        self.parameter_var = tk.StringVar(value="")

        self._build_layout()
        self.protocol("WM_DELETE_WINDOW", self._request_close)
        self._schedule_poll()

    def _build_layout(self) -> None:
        shell = ttk.Frame(self, style="Shell.TFrame", padding=SPACING["window_pad"])
        shell.pack(fill="both", expand=True)
        header = ttk.Frame(shell, style="Shell.TFrame")
        header.pack(fill="x", pady=(0, 12))
        ttk.Label(header, text="RIETAN-FP 受约束精修", style="Title.TLabel").pack(side="left")
        ttk.Label(header, textvariable=self.workflow_note_var, style="Muted.TLabel").pack(side="right")

        body = ttk.Panedwindow(shell, orient="horizontal")
        body.pack(fill="both", expand=True)
        control_shell = ttk.Frame(body, style="Panel.TFrame")
        control_shell.rowconfigure(0, weight=1)
        control_shell.columnconfigure(0, weight=1)
        self.control_canvas = tk.Canvas(
            control_shell,
            width=340,
            bg=COLORS["panel"],
            highlightthickness=0,
            borderwidth=0,
        )
        self.control_scrollbar = ttk.Scrollbar(
            control_shell,
            orient="vertical",
            command=self.control_canvas.yview,
        )
        self.control_canvas.configure(yscrollcommand=self.control_scrollbar.set)
        self.control_canvas.grid(row=0, column=0, sticky="nsew")
        self.control_scrollbar.grid(row=0, column=1, sticky="ns")
        controls = ttk.Frame(self.control_canvas, style="Panel.TFrame", padding=16)
        self.controls_frame = controls
        controls_window = self.control_canvas.create_window((0, 0), window=controls, anchor="nw")
        controls.bind(
            "<Configure>",
            lambda _event: self.control_canvas.configure(scrollregion=self.control_canvas.bbox("all")),
        )
        self.control_canvas.bind(
            "<Configure>",
            lambda event: self.control_canvas.itemconfigure(controls_window, width=event.width),
        )
        plot_panel = ttk.Frame(body, style="Panel.TFrame", padding=12)
        body.add(control_shell, weight=2)
        body.add(plot_panel, weight=6)
        self._build_controls(controls)
        self._bind_control_wheel(controls)
        self._build_plot(plot_panel)

    def _bind_control_wheel(self, widget: tk.Misc) -> None:
        widget.bind("<MouseWheel>", self._scroll_controls, add="+")
        for child in widget.winfo_children():
            self._bind_control_wheel(child)

    def _scroll_controls(self, event: tk.Event) -> None:
        delta = int(-event.delta / 120) if event.delta else 0
        if delta:
            self.control_canvas.yview_scroll(delta, "units")

    def _build_controls(self, pane: ttk.Frame) -> None:
        pane.columnconfigure(0, weight=1)
        ttk.Label(pane, text="工作流", style="Section.TLabel").grid(row=0, column=0, sticky="w")
        workflow = ttk.Frame(pane, style="Panel.TFrame")
        workflow.grid(row=1, column=0, sticky="ew", pady=(6, 14))
        ttk.Radiobutton(
            workflow,
            text="单相确认",
            value="单相确认",
            variable=self.workflow_var,
            command=self._on_workflow_changed,
            style="Toolbutton",
        ).pack(side="left", fill="x", expand=True)
        ttk.Radiobutton(
            workflow,
            text="多相实验定量",
            value="多相实验定量",
            variable=self.workflow_var,
            command=self._on_workflow_changed,
            style="Toolbutton",
        ).pack(side="left", fill="x", expand=True, padx=(6, 0))

        ttk.Label(pane, text="目标相", style="Section.TLabel").grid(row=2, column=0, sticky="w")
        self.target_combo = ttk.Combobox(
            pane,
            textvariable=self.target_var,
            values=tuple(self._candidate_by_label),
            state="readonly",
            height=10,
        )
        self.target_combo.grid(row=3, column=0, sticky="ew", pady=(6, 12))
        self.target_combo.bind("<<ComboboxSelected>>", self._on_target_changed)

        self.impurity_label = ttk.Label(pane, text="杂质相（最多三个）", style="Section.TLabel")
        self.impurity_label.grid(row=4, column=0, sticky="w")
        self.impurity_frame = ttk.Frame(pane, style="Panel.TFrame")
        self.impurity_frame.grid(row=5, column=0, sticky="ew", pady=(6, 14))
        self._rebuild_impurity_controls(reset=True)

        ttk.Label(pane, text="精修层级", style="Section.TLabel").grid(row=6, column=0, sticky="w")
        self.mode_combo = ttk.Combobox(
            pane,
            textvariable=self.mode_var,
            values=tuple(MODE_SETTINGS),
            state="readonly",
        )
        self.mode_combo.grid(row=7, column=0, sticky="ew", pady=(6, 10))
        ttk.Label(
            pane,
            text="建议先跑背景 + 尺度，再逐级开放零点和峰宽。峰宽结果需要仪器标准样校准。",
            style="PanelMuted.TLabel",
            wraplength=300,
            justify="left",
        ).grid(row=8, column=0, sticky="ew", pady=(0, 12))

        terms = ttk.Frame(pane, style="Panel.TFrame")
        terms.grid(row=9, column=0, sticky="ew")
        ttk.Label(terms, text="背景项数", style="Panel.TLabel").pack(side="left")
        ttk.Spinbox(terms, from_=2, to=10, textvariable=self.background_terms_var, width=5).pack(side="right")
        ttk.Checkbutton(
            pane,
            text="开放晶胞参数（谨慎）",
            variable=self.refine_lattice_var,
            style="Panel.TCheckbutton",
        ).grid(row=10, column=0, sticky="w", pady=(10, 4))
        ttk.Label(
            pane,
            textvariable=self.locked_parameters_var,
            style="PanelMuted.TLabel",
            wraplength=300,
            justify="left",
        ).grid(row=11, column=0, sticky="ew", pady=(0, 14))

        self.start_button = ttk.Button(pane, text="开始受约束精修", style="Primary.TButton", command=self._start)
        self.start_button.grid(row=12, column=0, sticky="ew")
        if not self.candidates or not self.backend.available:
            self.start_button.state(["disabled"])
        self.export_button = ttk.Button(pane, text="导出本次精修", command=self._export)
        self.export_button.grid(row=13, column=0, sticky="ew", pady=(8, 0))
        self.export_button.state(["disabled"])
        self.conclusion_button = ttk.Button(pane, text="查看样品结论", command=self._open_conclusion)
        self.conclusion_button.grid(row=14, column=0, sticky="ew", pady=(8, 0))
        self.conclusion_button.state(["disabled"])

        ttk.Separator(pane).grid(row=15, column=0, sticky="ew", pady=14)
        ttk.Label(pane, text="结果", style="Section.TLabel").grid(row=16, column=0, sticky="w")
        ttk.Label(pane, textvariable=self.metrics_var, style="Panel.TLabel", wraplength=300, justify="left").grid(
            row=17, column=0, sticky="ew", pady=(6, 4)
        )
        ttk.Label(pane, textvariable=self.parameter_var, style="PanelMuted.TLabel", wraplength=300, justify="left").grid(
            row=18, column=0, sticky="ew"
        )
        ttk.Label(pane, textvariable=self.status_var, style="PanelMuted.TLabel", wraplength=300, justify="left").grid(
            row=19, column=0, sticky="ew", pady=(14, 0)
        )
        self._on_workflow_changed()

    def _rebuild_impurity_controls(self, *, reset: bool = False) -> None:
        for child in self.impurity_frame.winfo_children():
            child.destroy()
        target_label = self.target_var.get()
        old_values = {label: variable.get() for label, variable in self.impurity_vars.items()}
        self.impurity_vars = {}
        selected_count = 0
        for label in self._candidate_by_label:
            if label == target_label:
                continue
            selected = old_values.get(label, False)
            if reset and selected_count < 3:
                selected = True
            if selected:
                selected_count += 1
            variable = tk.BooleanVar(value=selected and selected_count <= 3)
            self.impurity_vars[label] = variable
            ttk.Checkbutton(
                self.impurity_frame,
                text=label,
                variable=variable,
                command=self._limit_impurity_selection,
                style="Panel.TCheckbutton",
            ).pack(anchor="w", fill="x", pady=2)

    def _limit_impurity_selection(self) -> None:
        selected = [label for label, variable in self.impurity_vars.items() if variable.get()]
        if len(selected) <= 3:
            return
        self.impurity_vars[selected[-1]].set(False)
        messagebox.showwarning("杂质相过多", "首版联合精修最多选择三个杂质相。", parent=self)

    def _on_target_changed(self, _event: tk.Event | None = None) -> None:
        self._rebuild_impurity_controls(reset=True)

    def _on_workflow_changed(self) -> None:
        multiphase = self.workflow_var.get() == "多相实验定量"
        if multiphase:
            self.impurity_label.grid()
            self.impurity_frame.grid()
            self.start_button.configure(text="运行三组联合精修")
            self.workflow_note_var.set("多相联合精修；未过闸门不显示 wt%")
            self.conclusion_button.state(["disabled"])
            available = bool(self.candidates) and self.multiphase_backend.available
        else:
            self.impurity_label.grid_remove()
            self.impurity_frame.grid_remove()
            self.start_button.configure(text="开始受约束精修")
            self.workflow_note_var.set("单相确认，不输出 wt%")
            available = bool(self.candidates) and self.backend.available
        self.start_button.state(["!disabled"] if available and not self._running else ["disabled"])

    def selected_phases(self) -> list[Candidate]:
        target = self._candidate_by_label.get(self.target_var.get())
        selected = [] if target is None else [target]
        selected.extend(
            self._candidate_by_label[label]
            for label, variable in self.impurity_vars.items()
            if variable.get() and label in self._candidate_by_label
        )
        return selected[:4]

    def _build_multiphase_inputs(
        self,
        candidates: Iterable[Candidate] | None = None,
    ) -> tuple[PhaseRefinementInput, ...]:
        selected = list(candidates) if candidates is not None else self.selected_phases()
        if not 2 <= len(selected) <= 4:
            raise ValueError("多相实验定量需要一个目标相和一至三个杂质相。")
        return tuple(
            PhaseRefinementInput(
                candidate.candidate_id,
                candidate.formula_pretty,
                Path(candidate.source_path),
                "target" if index == 0 else "impurity",
            )
            for index, candidate in enumerate(selected)
        )

    def _build_plot(self, pane: ttk.Frame) -> None:
        pane.rowconfigure(0, weight=1)
        pane.columnconfigure(0, weight=1)
        configure_matplotlib_fonts()
        self.figure = Figure(figsize=(9.0, 6.5), dpi=100, facecolor=COLORS["panel"])
        self.pattern_axis, self.residual_axis = self.figure.subplots(
            2,
            1,
            sharex=True,
            gridspec_kw={"height_ratios": (3.3, 1.25), "hspace": 0.06},
        )
        self.canvas = FigureCanvasTkAgg(self.figure, master=pane)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        toolbar_frame = ttk.Frame(pane, style="Panel.TFrame")
        toolbar_frame.grid(row=1, column=0, sticky="ew")
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side="left", fill="x")
        self._show_initial_pattern()

    def _show_initial_pattern(self) -> None:
        self.pattern_axis.clear()
        self.residual_axis.clear()
        self.pattern_axis.plot(self.context.x, self.context.intensity, color="#111827", linewidth=0.85, label="实验谱")
        self.pattern_axis.set_ylabel("强度 (counts)")
        self.pattern_axis.legend(frameon=False, loc="upper right")
        self.residual_axis.axhline(0.0, color="#7b8798", linewidth=0.8)
        self.residual_axis.set_ylabel("差值")
        self.residual_axis.set_xlabel("2theta (degree)")
        self._style_axes()
        self.canvas.draw_idle()

    def _style_axes(self) -> None:
        for axis in (self.pattern_axis, self.residual_axis):
            axis.grid(color=COLORS["border_soft"], linewidth=0.5, alpha=0.65)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)

    def _start(self) -> None:
        if self._running:
            return
        if self.workflow_var.get() == "多相实验定量":
            self._start_multiphase()
            return
        candidate = self._candidate_by_label.get(self.target_var.get())
        if candidate is None:
            messagebox.showwarning("没有目标相", "请选择一个 CIF 目标相。", parent=self)
            return
        try:
            background_terms = int(self.background_terms_var.get())
            refine_zero, refine_profile = MODE_SETTINGS[self.mode_var.get()]
            settings = RefinementSettings(
                two_theta_min=self.context.two_theta_range[0],
                two_theta_max=self.context.two_theta_range[1],
                background_terms=background_terms,
                refine_zero_shift=refine_zero,
                refine_profile=refine_profile,
                refine_lattice=bool(self.refine_lattice_var.get()),
                radiation=self.context.radiation,
            )
        except (KeyError, TypeError, ValueError) as exc:
            messagebox.showerror("精修设置无效", str(exc), parent=self)
            return
        self._running = True
        self._active_candidate = candidate
        self._excluded_formulas.clear()
        self._conclusion_candidates = []
        self.start_button.state(["disabled"])
        self.export_button.state(["disabled"])
        self.conclusion_button.state(["disabled"])
        self.status_var.set("RIETAN-FP 正在精修；随后会受控检查常见氧化物。")
        threading.Thread(target=self._worker, args=(candidate, settings), daemon=True).start()

    def _start_multiphase(self) -> None:
        selected_candidates = self.selected_phases()
        try:
            phases = self._build_multiphase_inputs(selected_candidates)
            refine_zero, refine_profile = MODE_SETTINGS[self.mode_var.get()]
            settings = MultiphaseRefinementSettings(
                two_theta_min=self.context.two_theta_range[0],
                two_theta_max=self.context.two_theta_range[1],
                background_terms=int(self.background_terms_var.get()),
                refine_zero_shift=refine_zero,
                refine_profile=refine_profile,
                refine_lattice=bool(self.refine_lattice_var.get()),
                radiation=self.context.radiation,
            )
        except (KeyError, TypeError, ValueError, FileNotFoundError) as exc:
            messagebox.showerror("联合精修设置无效", str(exc), parent=self)
            return
        validation_states = {
            phase.phase_id: str(candidate.simulation_validation.get("status", "missing"))
            for phase, candidate in zip(phases, selected_candidates)
        }
        self._running = True
        self._active_candidate = selected_candidates[0]
        self.result = None
        self.conclusion = None
        self.multiphase_result = None
        self.quantification_assessment = None
        self.stability_results = ()
        self.start_button.state(["disabled"])
        self.export_button.state(["disabled"])
        self.conclusion_button.state(["disabled"])
        self.status_var.set("正在从三组初始尺度运行 RIETAN-FP 联合精修；未过闸门不会显示 wt%。")
        threading.Thread(
            target=self._multiphase_worker,
            args=(phases, settings, validation_states),
            daemon=True,
        ).start()

    def _multiphase_worker(
        self,
        phases: tuple[PhaseRefinementInput, ...],
        settings: MultiphaseRefinementSettings,
        validation_states: dict[str, str],
    ) -> None:
        try:
            results = tuple(
                self.multiphase_backend.refine(
                    variant,
                    self.context.x,
                    self.context.intensity,
                    settings,
                )
                for variant in initial_scale_variants(phases)
            )
            assessment = assess_quantification(
                results[0],
                results[1:],
                validation_states,
                stability_spread_limit_percent=settings.stability_spread_limit_percent,
            )
            self._queue.put(("multiphase_ok", (results[0], assessment, results)))
        except Exception as exc:
            self._queue.put(("error", exc))

    def _worker(self, candidate: Candidate, settings: RefinementSettings) -> None:
        try:
            result = self.backend.refine(
                candidate.source_path,
                self.context.x,
                self.context.intensity,
                settings,
            )
            allowed_elements = set(candidate.elements) - {"O"}
            references = discover_sibling_peak_references(
                self.context.source_path,
                allowed_elements=allowed_elements,
                target_formula=candidate.formula_pretty,
            )
            oxide_candidates: list[Candidate] = []
            oxide_formulas: tuple[str, ...] = ()
            oxide_warnings: tuple[str, ...] = ()
            if self.oxide_library_path:
                try:
                    oxide_result = load_controlled_oxide_candidates(
                        self.oxide_library_path,
                        set(candidate.elements) - {"O"},
                        radiation=self.context.radiation,
                        two_theta_range=self.context.two_theta_range,
                    )
                    oxide_candidates = list(oxide_result.candidates)
                    oxide_formulas = oxide_result.formulas_checked
                    oxide_warnings = oxide_result.warnings
                except Exception as exc:
                    oxide_warnings = (f"氧化物二次筛选未完成：{exc}",)
            conclusion = build_sample_conclusion(
                result,
                candidate,
                [*self.supporting_candidates, *references, *oxide_candidates],
                oxide_formulas_checked=oxide_formulas,
                oxide_screening_warnings=oxide_warnings,
            )
            self._conclusion_candidates = [*self.supporting_candidates, *references, *oxide_candidates]
            self._oxide_formulas = oxide_formulas
            self._oxide_warnings = oxide_warnings
            self._queue.put(("ok", (result, conclusion)))
        except Exception as exc:
            self._queue.put(("error", exc))

    def _poll_queue(self) -> None:
        self._poll_after_id = None
        if self._closed:
            return
        try:
            kind, payload = self._queue.get_nowait()
        except queue.Empty:
            self._schedule_poll()
            return
        if kind == "conclusion_ok":
            self._running = False
            self._on_workflow_changed()
            self.conclusion = payload
            self.conclusion_button.state(["!disabled"])
            self.status_var.set("已按人工排除条件重算杂质结论；RIETAN 目标相结果未重新运行。")
            if self.conclusion_dialog is not None:
                try:
                    self.conclusion_dialog.destroy()
                except tk.TclError:
                    pass
                self.conclusion_dialog = None
            self._open_conclusion()
            self._schedule_poll()
            return
        if kind == "conclusion_error":
            self._running = False
            self._on_workflow_changed()
            self.conclusion_button.state(["!disabled"])
            self.status_var.set(f"杂质结论重算失败：{payload}")
            self._schedule_poll()
            return
        self._running = False
        self._on_workflow_changed()
        if kind == "multiphase_ok":
            self._show_multiphase_result(payload[0], payload[1], payload[2])
        elif kind == "ok":
            if isinstance(payload, tuple) and len(payload) == 2:
                self._show_result(payload[0], payload[1])
            else:
                self._show_result(payload)
        else:
            self.status_var.set(f"精修失败：{payload}")
            messagebox.showerror("RIETAN-FP 精修失败", str(payload), parent=self)
        self._schedule_poll()

    def _schedule_poll(self) -> None:
        if not self._closed and self.winfo_exists():
            self._poll_after_id = self.after(100, self._poll_queue)

    def _show_result(
        self,
        result: RietanRefinementResult,
        conclusion: SampleConclusion | None = None,
    ) -> None:
        self.result = result
        self.conclusion = conclusion
        self.pattern_axis.clear()
        self.residual_axis.clear()
        self.pattern_axis.plot(
            result.two_theta_deg,
            result.observed_intensity,
            color="#111827",
            linewidth=0.75,
            label="实验谱",
        )
        self.pattern_axis.plot(
            result.two_theta_deg,
            result.calculated_intensity,
            color="#d93025",
            linewidth=1.0,
            label="计算谱",
        )
        self.pattern_axis.plot(
            result.two_theta_deg,
            result.background_intensity,
            color="#3974d8",
            linewidth=0.9,
            linestyle="--",
            label="背景",
        )
        tick_height = max(float(result.observed_intensity.max()) * 0.035, 1.0)
        self.pattern_axis.vlines(
            result.reflection_two_theta_deg,
            0.0,
            tick_height,
            color="#15a34a",
            linewidth=0.55,
            alpha=0.7,
            label="Bragg 峰位",
        )
        self.residual_axis.plot(
            result.two_theta_deg,
            result.residual_intensity,
            color="#7c3aed",
            linewidth=0.8,
            label="差值",
        )
        self.residual_axis.axhline(0.0, color="#7b8798", linewidth=0.8)
        self.pattern_axis.set_ylabel("强度 (counts)")
        self.residual_axis.set_ylabel("实验 - 计算")
        self.residual_axis.set_xlabel("2theta (degree)")
        self.pattern_axis.legend(frameon=False, loc="upper right", ncols=4, fontsize=8)
        self.residual_axis.legend(frameon=False, loc="upper right", fontsize=8)
        self._style_axes()
        self.canvas.draw_idle()

        s_text = f"S {result.s_value:.4f}\n" if result.s_value is not None else ""
        self.metrics_var.set(
            f"Rwp {result.rwp_percent:.3f}%\nRp {result.rp_percent:.3f}%\n"
            f"{s_text}GofF {result.goodness_of_fit:.4f}"
        )
        labels = {
            "zero_shift": "零点偏移",
            "scale": "尺度",
            "fwhm_u": "U",
            "fwhm_v": "V",
            "fwhm_w": "W",
            "cell_a": "a",
            "cell_b": "b",
            "cell_c": "c",
        }
        self.parameter_var.set(
            "\n".join(f"{labels[key]}: {value:.7g}" for key, value in result.parameters.items() if key in labels)
        )
        if result.warnings:
            self.status_var.set("需要人工检查：\n" + "\n".join(f"- {warning}" for warning in result.warnings))
        else:
            self.status_var.set("精修完成；当前参数未触发 Ophi 的基本可信度警告。")
        self.export_button.state(["!disabled"])
        if conclusion is not None:
            self.conclusion_button.state(["!disabled"])
        self.after_idle(lambda: self.control_canvas.yview_moveto(1.0))

    def _show_multiphase_result(
        self,
        result: MultiphaseRefinementResult,
        assessment: QuantificationAssessment,
        stability_results: tuple[MultiphaseRefinementResult, ...],
    ) -> None:
        self.result = None
        self.conclusion = None
        self.multiphase_result = result
        self.quantification_assessment = assessment
        self.stability_results = tuple(stability_results)
        self.pattern_axis.clear()
        self.residual_axis.clear()
        self.pattern_axis.plot(
            result.two_theta_deg,
            result.observed_intensity,
            color="#111827",
            linewidth=0.72,
            label="实验谱",
        )
        self.pattern_axis.plot(
            result.two_theta_deg,
            result.calculated_intensity,
            color="#d93025",
            linewidth=1.0,
            label="RIETAN 总计算谱",
        )
        self.pattern_axis.plot(
            result.two_theta_deg,
            result.background_intensity,
            color="#3974d8",
            linewidth=0.85,
            linestyle="--",
            label="背景",
        )
        colors = ("#15803d", "#d97706", "#0891b2", "#be123c")
        observed_min = float(np.min(result.observed_intensity))
        observed_span = max(float(np.ptp(result.observed_intensity)), 1.0)
        for index, phase in enumerate(result.phases):
            baseline = observed_min - observed_span * (0.045 + index * 0.025)
            self.pattern_axis.vlines(
                phase.reflection_two_theta_deg,
                baseline,
                baseline + observed_span * 0.016,
                color=colors[index % len(colors)],
                linewidth=0.75,
                alpha=0.9,
                label=f"{phase.formula} 反射",
            )
        self.residual_axis.plot(
            result.two_theta_deg,
            result.residual_intensity,
            color="#6d28d9",
            linewidth=0.78,
            label="实验 - 计算",
        )
        self.residual_axis.axhline(0.0, color="#7b8798", linewidth=0.8)
        self.pattern_axis.set_ylabel("强度 (counts)")
        self.residual_axis.set_ylabel("差值")
        self.residual_axis.set_xlabel("2theta (degree)")
        self.pattern_axis.legend(frameon=False, loc="upper right", ncols=3, fontsize=8)
        self.residual_axis.legend(frameon=False, loc="upper right", fontsize=8)
        self._style_axes()
        self.canvas.draw_idle()

        self.metrics_var.set(
            f"{assessment.label}\nRwp {result.rwp_percent:.3f}%\n"
            f"Rp {result.rp_percent:.3f}%\nGofF {result.goodness_of_fit:.4f}"
        )
        if assessment.allow_weight_percent:
            lines = []
            for phase in result.phases:
                stability = assessment.stability_ranges.get(phase.phase_id)
                range_text = "" if stability is None else f"  [{stability[0]:.3f}, {stability[1]:.3f}]"
                lines.append(f"{phase.formula}: {phase.weight_percent:.3f} wt%{range_text}")
            self.parameter_var.set("\n".join(lines))
        else:
            self.parameter_var.set("各相精确 wt% 已隐藏。请先解决下列可信度问题。")
        if assessment.findings:
            self.status_var.set("可信度检查：\n" + "\n".join(f"- {item.message}" for item in assessment.findings))
        else:
            self.status_var.set("三组初始尺度得到稳定结果；当前结果标记为实验性定量，尚未经过实测标样认证。")
        self.export_button.state(["!disabled"])
        self.conclusion_button.state(["disabled"])
        self.after_idle(lambda: self.control_canvas.yview_moveto(1.0))

    def _open_conclusion(self) -> None:
        if self.conclusion is None:
            return
        if self.conclusion_dialog is not None:
            try:
                if self.conclusion_dialog.winfo_exists():
                    self.conclusion_dialog.lift()
                    self.conclusion_dialog.focus_force()
                    return
            except tk.TclError:
                pass
        self.conclusion_dialog = SampleConclusionDialog(
            self,
            self.conclusion,
            on_exclude_formula=self._exclude_candidate_formula,
        )

    def _exclude_candidate_formula(self, formula: str) -> None:
        if self._running or self.result is None or self._active_candidate is None:
            return
        clean = formula.strip()
        if not clean or clean.lower() == self._active_candidate.formula_pretty.strip().lower():
            return
        self._excluded_formulas.add(clean)
        self._running = True
        self.start_button.state(["disabled"])
        self.conclusion_button.state(["disabled"])
        self.status_var.set(f"正在排除 {clean} 并重算杂质结论；不会重新运行 RIETAN。")
        threading.Thread(target=self._conclusion_worker, daemon=True).start()

    def _conclusion_worker(self) -> None:
        try:
            assert self.result is not None and self._active_candidate is not None
            excluded_keys = {item.lower() for item in self._excluded_formulas}
            candidates = [
                item
                for item in self._conclusion_candidates
                if item.formula_pretty.strip().lower() not in excluded_keys
            ]
            conclusion = build_sample_conclusion(
                self.result,
                self._active_candidate,
                candidates,
                oxide_formulas_checked=self._oxide_formulas,
                oxide_screening_warnings=self._oxide_warnings,
                excluded_formulas=sorted(self._excluded_formulas),
            )
            self._queue.put(("conclusion_ok", conclusion))
        except Exception as exc:
            self._queue.put(("conclusion_error", exc))

    def _export(self) -> None:
        if self.result is None and self.multiphase_result is None:
            return
        target = filedialog.askdirectory(title="选择精修导出文件夹", parent=self)
        if not target:
            return
        root = Path(target)
        root.mkdir(parents=True, exist_ok=True)
        if self.multiphase_result is not None:
            self._export_multiphase(root)
            return
        self.figure.savefig(root / "rietan_refinement.png", dpi=300, bbox_inches="tight")
        with (root / "rietan_refinement_profile.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("two_theta", "observed", "calculated", "residual", "background"))
            writer.writerows(
                zip(
                    self.result.two_theta_deg,
                    self.result.observed_intensity,
                    self.result.calculated_intensity,
                    self.result.residual_intensity,
                    self.result.background_intensity,
                )
            )
        payload = {
            "rwp_percent": self.result.rwp_percent,
            "rp_percent": self.result.rp_percent,
            "s_value": self.result.s_value,
            "goodness_of_fit": self.result.goodness_of_fit,
            "parameters": self.result.parameters,
            "warnings": list(self.result.warnings),
            "provenance": self.result.provenance,
            "sample_conclusion": None
            if self.conclusion is None
            else {
                "target_formula": self.conclusion.target_formula,
                "target_evidence_label": self.conclusion.target_evidence_label,
                "competing_models": self.conclusion.competing_models,
                "oxide_formulas_checked": list(self.conclusion.oxide_formulas_checked),
                "oxide_screening_warnings": list(self.conclusion.oxide_screening_warnings),
                "excluded_formulas": list(self.conclusion.excluded_formulas),
                "impurity_estimates": [asdict(item) for item in self.conclusion.impurity_estimates],
                "report_text": self.conclusion.to_chinese_text(),
            },
        }
        (root / "rietan_refinement.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        self.status_var.set(f"精修结果已导出：{root}")

    def _export_multiphase(self, root: Path) -> None:
        result = self.multiphase_result
        assessment = self.quantification_assessment
        if result is None or assessment is None:
            return
        self.figure.savefig(root / "rietan_multiphase_refinement.png", dpi=300, bbox_inches="tight")
        with (root / "rietan_multiphase_profile.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("two_theta", "observed", "calculated", "residual", "background"))
            writer.writerows(
                zip(
                    result.two_theta_deg,
                    result.observed_intensity,
                    result.calculated_intensity,
                    result.residual_intensity,
                    result.background_intensity,
                )
            )
        with (root / "rietan_multiphase_reflections.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("phase_id", "formula", "two_theta"))
            for phase in result.phases:
                writer.writerows((phase.phase_id, phase.formula, value) for value in phase.reflection_two_theta_deg)
        payload = {
            "quantification_status": assessment.label,
            "allow_weight_percent": assessment.allow_weight_percent,
            "rwp_percent": result.rwp_percent,
            "rp_percent": result.rp_percent,
            "goodness_of_fit": result.goodness_of_fit,
            "phases": [
                {
                    "phase_id": phase.phase_id,
                    "formula": phase.formula,
                    "scale": phase.scale,
                    "z": phase.z,
                    "molar_mass": phase.molar_mass,
                    "volume_angstrom3": phase.volume_angstrom3,
                    "weight_percent": phase.weight_percent if assessment.allow_weight_percent else None,
                    "stability_range_wt_percent": assessment.stability_ranges.get(phase.phase_id),
                    "reflection_two_theta_deg": list(phase.reflection_two_theta_deg),
                }
                for phase in result.phases
            ],
            "findings": [
                {"code": item.code, "level": item.level.value, "message": item.message}
                for item in assessment.findings
            ],
            "warnings": list(result.warnings),
            "provenance": result.provenance,
            "stability_runs": [
                {
                    "rwp_percent": run.rwp_percent,
                    "weights": {phase.phase_id: phase.weight_percent for phase in run.phases},
                }
                for run in self.stability_results
            ],
        }
        (root / "rietan_multiphase_refinement.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        self.status_var.set(f"多相精修证据包已导出：{root}")

    def _request_close(self) -> None:
        self._closed = True
        if self.conclusion_dialog is not None:
            try:
                self.conclusion_dialog.destroy()
            except tk.TclError:
                pass
        self.destroy()
        if self._on_closed is not None:
            self._on_closed()

    def destroy(self) -> None:
        self._closed = True
        canvas = getattr(self, "canvas", None)
        draw_id = getattr(canvas, "_idle_draw_id", None)
        if draw_id is not None:
            try:
                self.after_cancel(draw_id)
            except tk.TclError:
                pass
            canvas._idle_draw_id = None
        if self._poll_after_id is not None:
            try:
                self.after_cancel(self._poll_after_id)
            except tk.TclError:
                pass
            self._poll_after_id = None
        try:
            super().destroy()
        except tk.TclError:
            pass


def open_or_raise_refinement_window(
    parent: tk.Misc,
    existing: RefinementWindow | None,
    *,
    context: AnalysisContext,
    candidates: Iterable[Candidate],
    backend: RietanRefinementBackend | None = None,
    multiphase_backend: RietanMultiphaseBackend | None = None,
    supporting_candidates: Iterable[Candidate] | None = None,
    oxide_library_path: str | Path | None = None,
    on_closed: Callable[[], None] | None = None,
) -> RefinementWindow:
    if existing is not None:
        try:
            if existing.winfo_exists():
                existing.deiconify()
                existing.lift()
                existing.focus_force()
                return existing
        except tk.TclError:
            pass
    return RefinementWindow(
        parent,
        context,
        candidates,
        backend=backend,
        multiphase_backend=multiphase_backend,
        supporting_candidates=supporting_candidates,
        oxide_library_path=oxide_library_path,
        on_closed=on_closed,
    )
