import tempfile
import unittest
from os import environ
from pathlib import Path
from unittest.mock import patch

from ophiuchus.library.models import StructureEntry
from ophiuchus.library.database import StructureLibrary
from ophiuchus.library.local_provider import LocalFolderProvider
from ophiuchus.xrd.config import XRDConfig
from ophiuchus.xrd.matching import prioritize_impurities_after_main, score_candidate
from ophiuchus.xrd.models import Candidate, Peak
from ophiuchus.xrd.models import PatternPoint, XrdPattern


ROOT = Path(__file__).parents[1]
REAL_ZR_CIF = ROOT / "xrd_sandbox" / "inputs" / "Zr3V3GeSn4.cif"
REAL_TI_CIF = ROOT / "xrd_sandbox" / "inputs" / "Ti3V3GeSn4.cif"


class DirectXrdBackendTests(unittest.TestCase):
    def test_explicit_low_score_target_stays_first(self):
        experimental = [Peak("e1", 20.0, 100.0), Peak("e2", 30.0, 80.0)]
        target = Candidate("target", "Target", "library:local", "target.cif", ["A", "B"], "target-hash")
        target.theory_peaks = [Peak("t1", 20.0, 100.0), Peak("t2", 60.0, 100.0)]
        tin = Candidate("tin", "Sn", "library:materials_project", "sn.cif", ["Sn"], "sn-hash")
        tin.theory_peaks = [Peak("s1", 20.0, 100.0), Peak("s2", 30.0, 80.0)]
        scores = [score_candidate(experimental, target), score_candidate(experimental, tin)]

        ordered = prioritize_impurities_after_main(
            scores,
            experimental,
            main_elements=["A", "B"],
            target_candidate_id="target",
        )

        self.assertLess(scores[0].score, 0.75)
        self.assertEqual(ordered[0].candidate.candidate_id, "target")

    def test_validated_backend_preserves_full_real_zr_pattern(self):
        from ophiuchus.xrd.backend import SimulationContext, ValidatedXRDBackend

        config = XRDConfig(
            radiation_source="CuKalpha12",
            wavelength_angstrom=1.54056,
            two_theta_min=10.0,
            two_theta_max=90.0,
            intensity_threshold=1.0,
        )
        pattern = ValidatedXRDBackend().simulate_cif(
            REAL_ZR_CIF,
            config,
            SimulationContext(structure_id="real-zr", source="local"),
        )

        self.assertEqual(len(pattern.two_theta_deg), 146)
        self.assertEqual(len(pattern.raw_intensity), 146)
        self.assertEqual(len(pattern.normalized_intensity), 146)
        self.assertEqual(max(pattern.normalized_intensity), 100.0)
        self.assertTrue(any(value < 1.0 for value in pattern.normalized_intensity))
        self.assertEqual(pattern.cif_sha256, pattern.calculate_cif_sha256(REAL_ZR_CIF))

    def test_real_ti_p1_near_degenerate_reflections_are_merged_once_before_normalization(self):
        from ophiuchus.xrd.backend import SimulationContext, ValidatedXRDBackend

        pattern = ValidatedXRDBackend().simulate_cif(
            REAL_TI_CIF,
            XRDConfig(two_theta_min=2.5, two_theta_max=90.5, peak_merge_tolerance_deg=0.03),
            SimulationContext(structure_id="real-ti", source="local", formula="Ti3V3GeSn4"),
        )
        kalpha1 = [
            (two_theta, intensity)
            for two_theta, intensity, component in zip(
                pattern.two_theta_deg, pattern.normalized_intensity, pattern.line_component
            )
            if component == "Kalpha1"
        ]
        expected_vesta = {
            28.01: 60.05339,
            30.08: 24.96655,
            32.92: 21.56566,
            34.40: 44.37757,
            38.12: 59.62114,
            39.73: 73.00159,
            40.32: 100.0,
            44.79: 15.30658,
        }

        self.assertEqual(sum(1 for x, _ in kalpha1 if 39.70 <= x <= 39.75), 1)
        for position, expected_intensity in expected_vesta.items():
            actual_position, actual_intensity = min(kalpha1, key=lambda row: abs(row[0] - position))
            self.assertLessEqual(abs(actual_position - position), 0.03)
            self.assertAlmostEqual(actual_intensity, expected_intensity, delta=15.0)

    def test_first_three_impurity_slots_prefer_distinct_formulas(self):
        experimental = [Peak(f"e{i}", 20.0 + i * 5.0, 100.0 - i * 5.0) for i in range(6)]
        target = Candidate("target", "Target", "library:local", "target.cif", ["A"], "target")
        target.theory_peaks = [Peak("t", 20.0, 100.0)]
        candidates = [target]
        for index, formula in enumerate(["Ti", "Ti", "Ti", "Ge", "Sn"], 1):
            candidate = Candidate(f"c{index}", formula, "library:test", f"{index}.cif", [formula], str(index))
            candidate.theory_peaks = [Peak(f"p{index}", 25.0 + min(index, 4) * 5.0, 100.0)]
            candidates.append(candidate)
        ordered = prioritize_impurities_after_main(
            [score_candidate(experimental, candidate, tolerance_deg=0.2) for candidate in candidates],
            experimental,
            target_candidate_id="target",
        )

        self.assertEqual(ordered[0].candidate.candidate_id, "target")
        self.assertEqual(len({score.candidate.formula_pretty for score in ordered[1:4]}), 3)

    def test_curated_vesta_top_peak_table_validates_coverage_without_claiming_exhaustiveness(self):
        from ophiuchus.xrd.backend import SimulationContext, ValidatedXRDBackend
        from ophiuchus.xrd.simulation_trust import apply_vesta_trust_check

        pattern = ValidatedXRDBackend().simulate_cif(
            REAL_TI_CIF,
            XRDConfig(two_theta_min=2.5, two_theta_max=90.5),
            SimulationContext("ti", "local", "Ti3V3GeSn4"),
        )
        candidate = Candidate("ti", "Ti3V3GeSn4", "library:local", str(REAL_TI_CIF), ["Ti", "V", "Ge", "Sn"], "ti")
        candidate.simulated_pattern = pattern
        candidate.theory_peaks = pattern.to_peaks()
        with tempfile.TemporaryDirectory() as tmp:
            phase_dir = Path(tmp) / "Ti3V3GeSn4"
            phase_dir.mkdir()
            (phase_dir / "MONI.int").write_text("28.01 60\n30.08 25\n32.92 22\n34.40 44\n38.12 60\n39.73 73\n40.32 100\n44.79 15\n", encoding="utf-8")
            (phase_dir / "simulated_peak_positions.csv").write_text(
                "two_theta,intensity_norm\n28.01,60\n30.08,25\n32.92,22\n34.40,44\n38.12,60\n39.73,73\n40.32,100\n44.79,15\n",
                encoding="utf-8",
            )
            with patch.dict(environ, {"OPHI_VESTA_REFERENCE_DIR": tmp}, clear=False):
                status = apply_vesta_trust_check(candidate, 2.5, 90.5)

        self.assertEqual(status["status"], "passed")
        self.assertEqual(status["reference_scope"], "curated_major_peaks")
        self.assertGreater(status["extra_strong_count"], 0)
        self.assertFalse(status["extra_strong_used_as_failure"])

    def test_public_simulation_route_does_not_replace_cif_with_vesta_reference(self):
        from ophiuchus.xrd.validation import simulate_cif_with_config

        with tempfile.TemporaryDirectory() as tmp:
            reference = Path(tmp) / "Zr3V3GeSn4 VESTA.int"
            reference.write_text("20.0 100\n", encoding="utf-8")
            config = XRDConfig(two_theta_min=10.0, two_theta_max=90.0)
            with patch.dict(environ, {"OPHI_VESTA_REFERENCE_DIR": tmp}, clear=False):
                peaks = simulate_cif_with_config(REAL_ZR_CIF, config)

        self.assertEqual(len(peaks), 146)
        self.assertFalse(any(abs(peak.two_theta - 20.0) < 1e-10 for peak in peaks))

    def test_candidate_route_uses_canonical_backend_without_dropping_weak_peaks(self):
        from ophiuchus.xrd.candidates import simulate_or_load_peaks

        candidate = Candidate(
            "real-zr",
            "Zr3V3GeSn4",
            "local_cif",
            str(REAL_ZR_CIF),
            ["Zr", "V", "Ge", "Sn"],
            "fixture-hash",
        )
        peaks = simulate_or_load_peaks(candidate, radiation="CuKalpha12", two_theta_range=(10.0, 90.0))

        self.assertEqual(len(peaks), 146)
        self.assertIsNotNone(candidate.simulated_pattern)
        self.assertEqual(
            tuple(peak.two_theta for peak in peaks),
            candidate.simulated_pattern.two_theta_deg,
        )

    def test_phase_grouping_merges_duplicate_records_but_keeps_polymorphs(self):
        from ophiuchus.library.phase_grouping import group_phase_candidates
        from ophiuchus.xrd.backend import SimulationContext, ValidatedXRDBackend

        config = XRDConfig(two_theta_min=10.0, two_theta_max=90.0)
        pattern = ValidatedXRDBackend().simulate_cif(
            REAL_ZR_CIF,
            config,
            SimulationContext(structure_id="a", source="local"),
        )
        duplicate = StructureEntry(
            "a", "local", "a.cif", "Zr3V3GeSn4", "Zr3V3GeSn4",
            ["Zr", "V", "Ge", "Sn"], str(REAL_ZR_CIF), "same-hash", space_group_number=1,
        )
        same_phase = StructureEntry(
            "b", "materials_project", "mp-a", "Zr3V3GeSn4", "Zr3V3GeSn4",
            ["Zr", "V", "Ge", "Sn"], str(REAL_ZR_CIF), "same-hash", space_group_number=1,
        )
        polymorph = StructureEntry(
            "c", "materials_project", "mp-b", "Zr3V3GeSn4", "Zr3V3GeSn4",
            ["Zr", "V", "Ge", "Sn"], str(REAL_ZR_CIF), "other-hash", space_group_number=2,
        )

        groups = group_phase_candidates(
            [duplicate, same_phase, polymorph],
            {"a": pattern, "b": pattern, "c": pattern},
            target_structure_id="a",
        )

        self.assertEqual(len(groups), 2)
        target_group = next(group for group in groups if group.contains_structure("a"))
        self.assertEqual({entry.internal_id for entry in target_group.entries}, {"a", "b"})
        self.assertEqual(target_group.representative.internal_id, "a")

    def test_scientific_safe_mode_bypasses_cache_and_v2_cache_roundtrips_losslessly(self):
        from ophiuchus.library.xrd_cache import build_library_xrd_cache, load_validated_pattern

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cifs"
            source.mkdir()
            target = source / REAL_ZR_CIF.name
            target.write_bytes(REAL_ZR_CIF.read_bytes())
            library = StructureLibrary(root / "library.sqlite")
            LocalFolderProvider().import_folder(source, library)
            entry = library.list_structures()[0]

            first = build_library_xrd_cache(
                library,
                structure_ids=[entry.internal_id],
                radiation="CuKalpha12",
                two_theta_range=(10.0, 90.0),
                scientific_safe_mode=True,
            )
            second = build_library_xrd_cache(
                library,
                structure_ids=[entry.internal_id],
                radiation="CuKalpha12",
                two_theta_range=(10.0, 90.0),
                scientific_safe_mode=True,
            )
            normal = build_library_xrd_cache(
                library,
                structure_ids=[entry.internal_id],
                radiation="CuKalpha12",
                two_theta_range=(10.0, 90.0),
                scientific_safe_mode=False,
            )
            restored = load_validated_pattern(
                library,
                entry,
                radiation="CuKalpha12",
                two_theta_range=(10.0, 90.0),
            )
            from ophiuchus.xrd.backend import SimulationContext, ValidatedXRDBackend

            expected = ValidatedXRDBackend().simulate_cif(
                target,
                XRDConfig(two_theta_min=10.0, two_theta_max=90.0),
                SimulationContext(entry.internal_id, entry.source, entry.reduced_formula),
            )

        self.assertEqual(first["freshly_simulated"], 1)
        self.assertEqual(second["freshly_simulated"], 1)
        self.assertEqual(second["validated_cache_hits"], 0)
        self.assertEqual(normal["validated_cache_hits"], 1)
        self.assertIsNotNone(restored)
        self.assertEqual(len(restored.two_theta_deg), 146)
        self.assertEqual(restored.two_theta_deg, expected.two_theta_deg)
        self.assertEqual(restored.raw_intensity, expected.raw_intensity)
        self.assertEqual(restored.normalized_intensity, expected.normalized_intensity)

    def test_v2_cache_invalidates_when_range_changes(self):
        from ophiuchus.library.xrd_cache import build_library_xrd_cache, load_validated_pattern

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cifs"
            source.mkdir()
            target = source / REAL_ZR_CIF.name
            target.write_bytes(REAL_ZR_CIF.read_bytes())
            library = StructureLibrary(root / "library.sqlite")
            LocalFolderProvider().import_folder(source, library)
            entry = library.list_structures()[0]
            build_library_xrd_cache(
                library,
                structure_ids=[entry.internal_id],
                radiation="CuKalpha12",
                two_theta_range=(10.0, 90.0),
                scientific_safe_mode=False,
            )

            wrong_range = load_validated_pattern(
                library,
                entry,
                radiation="CuKalpha12",
                two_theta_range=(20.0, 80.0),
            )

        self.assertIsNone(wrong_range)

    def test_v2_cache_invalidates_for_changed_cif_radiation_or_backend(self):
        from ophiuchus.library.xrd_cache import build_library_xrd_cache, load_validated_pattern
        from ophiuchus.xrd.backend import ValidatedXRDBackend

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cifs"
            source.mkdir()
            target = source / REAL_ZR_CIF.name
            target.write_bytes(REAL_ZR_CIF.read_bytes())
            library = StructureLibrary(root / "library.sqlite")
            LocalFolderProvider().import_folder(source, library)
            entry = library.list_structures()[0]
            build_library_xrd_cache(
                library,
                structure_ids=[entry.internal_id],
                radiation="CuKalpha12",
                two_theta_range=(10.0, 90.0),
                scientific_safe_mode=False,
            )

            wrong_radiation = load_validated_pattern(
                library, entry, radiation="CoKa", two_theta_range=(10.0, 90.0)
            )
            with patch.object(ValidatedXRDBackend, "engine_version", "future-backend"):
                wrong_backend = load_validated_pattern(
                    library, entry, radiation="CuKalpha12", two_theta_range=(10.0, 90.0)
                )
            cached_cif = library.path.parent / entry.cached_file_path
            cached_cif.write_bytes(cached_cif.read_bytes() + b"\n# changed\n")
            changed_cif = load_validated_pattern(
                library, entry, radiation="CuKalpha12", two_theta_range=(10.0, 90.0)
            )

        self.assertIsNone(wrong_radiation)
        self.assertIsNone(wrong_backend)
        self.assertIsNone(changed_cif)

    def test_materials_project_export_route_matches_direct_backend(self):
        from pymatgen.core import Lattice, Structure

        from ophiuchus.library.mp_provider import MaterialsProjectProvider
        from ophiuchus.library.xrd_cache import build_library_xrd_cache, load_validated_pattern
        from ophiuchus.xrd.backend import SimulationContext, ValidatedXRDBackend

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            library = StructureLibrary(root / "library.sqlite")
            structure = Structure(Lattice.cubic(4.0), ["Fe"], [[0, 0, 0]])
            provider = MaterialsProjectProvider(api_key="test")
            entry = provider._save_record(
                library,
                {
                    "material_id": "mp-test",
                    "formula_pretty": "Fe",
                    "elements": ["Fe"],
                    "chemsys": "Fe",
                    "structure": structure,
                    "symmetry": {"symbol": "Pm-3m", "number": 221, "crystal_system": "cubic"},
                },
            )
            library.upsert_structure(entry)
            build_library_xrd_cache(
                library,
                structure_ids=[entry.internal_id],
                radiation="CuKalpha12",
                two_theta_range=(10.0, 90.0),
                scientific_safe_mode=True,
            )
            cached = load_validated_pattern(
                library, entry, radiation="CuKalpha12", two_theta_range=(10.0, 90.0)
            )
            cif_path = library.path.parent / entry.cached_file_path
            direct = ValidatedXRDBackend().simulate_cif(
                cif_path,
                XRDConfig(two_theta_min=10.0, two_theta_max=90.0),
                SimulationContext(entry.internal_id, entry.source, entry.reduced_formula, entry.space_group_number),
            )

        self.assertIsNotNone(cached)
        self.assertEqual(cached.cif_sha256, direct.cif_sha256)
        self.assertEqual(cached.settings_fingerprint, direct.settings_fingerprint)
        self.assertEqual(cached.two_theta_deg, direct.two_theta_deg)
        self.assertEqual(cached.raw_intensity, direct.raw_intensity)
        self.assertEqual(cached.normalized_intensity, direct.normalized_intensity)

    def test_structure_library_analysis_keeps_selected_target_pattern_and_groups_duplicate_entries(self):
        from ophiuchus.library.analysis import run_library_analysis
        from ophiuchus.xrd.backend import SimulationContext, ValidatedXRDBackend

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cifs"
            source.mkdir()
            target_cif = source / REAL_ZR_CIF.name
            target_cif.write_bytes(REAL_ZR_CIF.read_bytes())
            library = StructureLibrary(root / "library.sqlite")
            provider = LocalFolderProvider()
            provider.import_folder(source, library, source="local")
            provider.import_folder(source, library, source="target")
            target_entry = next(entry for entry in library.list_structures() if entry.source == "target")
            config = XRDConfig(two_theta_min=10.0, two_theta_max=90.0)
            direct = ValidatedXRDBackend().simulate_cif(
                target_cif,
                config,
                SimulationContext(target_entry.internal_id, target_entry.source, target_entry.reduced_formula),
            )
            xrd = root / "sample.xy"
            rows = []
            strongest = sorted(
                zip(direct.two_theta_deg, direct.normalized_intensity),
                key=lambda item: item[1],
                reverse=True,
            )[:12]
            for index in range(4001):
                x = 10.0 + index * 0.02
                y = 1.0
                for center, intensity in strongest:
                    y += intensity * pow(2.718281828, -((x - center) ** 2) / (2 * 0.06**2))
                rows.append(f"{x:.4f} {y:.8f}")
            xrd.write_text("\n".join(rows), encoding="utf-8")
            reference_dir = root / "vesta"
            reference_dir.mkdir()
            (reference_dir / "Zr3V3GeSn4 VESTA.int").write_text("20.0 100\n", encoding="utf-8")

            with patch.dict(environ, {"OPHI_VESTA_REFERENCE_DIR": str(reference_dir)}, clear=False):
                result = run_library_analysis(
                    xrd_file=xrd,
                    library_db=library.path,
                    elements=["Zr", "V", "Ge", "Sn"],
                    out_dir=root / "out",
                    radiation="CuKalpha12",
                    two_theta_min=10.0,
                    two_theta_max=90.0,
                    target_candidate_id=target_entry.internal_id,
                    scientific_safe_mode=True,
                )

        self.assertEqual(len(result.phase_candidates), 1)
        self.assertEqual(len(result.phase_candidates[0].entries), 2)
        self.assertIsNotNone(result.target_score)
        self.assertEqual(result.top_scores[0].candidate.candidate_id, result.target_score.candidate.candidate_id)
        self.assertEqual(
            result.target_score.candidate.simulated_pattern.two_theta_deg,
            direct.two_theta_deg,
        )

    def test_matplotlib_receives_canonical_target_two_theta_without_vesta_replacement(self):
        from ophiuchus.xrd.backend import SimulationContext, ValidatedXRDBackend
        from ophiuchus.xrd.plotting import build_xrd_figure
        from ophiuchus.xrd.models import CandidateScore

        direct = ValidatedXRDBackend().simulate_cif(
            REAL_ZR_CIF,
            XRDConfig(two_theta_min=10.0, two_theta_max=90.0),
            SimulationContext("target", "local", "Zr3V3GeSn4"),
        )
        candidate = Candidate(
            "target", "Zr3V3GeSn4", "phase:local", str(REAL_ZR_CIF),
            ["Zr", "V", "Ge", "Sn"], direct.cif_sha256,
            simulated_pattern=direct,
        )
        candidate.theory_peaks = direct.to_peaks()
        score = CandidateScore(candidate, 0.5, [], [], [], [])
        experimental = XrdPattern(
            [PatternPoint(10.0, 1.0), PatternPoint(50.0, 2.0), PatternPoint(90.0, 1.0)]
        )
        with tempfile.TemporaryDirectory() as tmp:
            reference = Path(tmp) / "Zr3V3GeSn4 VESTA.int"
            reference.write_text("20.0 100\n", encoding="utf-8")
            with patch.dict(environ, {"OPHI_VESTA_REFERENCE_DIR": tmp}, clear=False):
                figure = build_xrd_figure(experimental, [], [score], max_candidates=4)
        try:
            target_collection = next(
                collection
                for collection in figure.axes[1].collections
                if collection.get_label() == "_ophi_target_pattern"
            )
            plotted_x = tuple(float(segment[0][0]) for segment in target_collection.get_segments())
        finally:
            import matplotlib.pyplot as plt

            plt.close(figure)

        self.assertEqual(plotted_x, direct.two_theta_deg)
        self.assertNotIn(20.0, plotted_x)

    def test_report_serializes_full_canonical_pattern_for_replotting(self):
        import json

        from ophiuchus.xrd.backend import SimulationContext, ValidatedXRDBackend
        from ophiuchus.xrd.models import CandidateScore, MultiPhaseExplanation
        from ophiuchus.xrd.report import write_reports

        direct = ValidatedXRDBackend().simulate_cif(
            REAL_ZR_CIF,
            XRDConfig(two_theta_min=10.0, two_theta_max=90.0),
            SimulationContext("target", "local", "Zr3V3GeSn4"),
        )
        candidate = Candidate(
            "target", "Zr3V3GeSn4", "phase:local", str(REAL_ZR_CIF),
            ["Zr", "V", "Ge", "Sn"], direct.cif_sha256,
            simulated_pattern=direct,
            phase_entry_ids=["target", "duplicate"],
            space_group_symbol="P1",
            space_group_number=1,
        )
        candidate.theory_peaks = direct.to_peaks()
        score = CandidateScore(candidate, 0.5, [], [], [], [])
        explanation = MultiPhaseExplanation([score], [], [])
        with tempfile.TemporaryDirectory() as tmp:
            outputs = write_reports(
                tmp,
                {"two_theta_range": [10.0, 90.0]},
                [],
                [score],
                explanation,
            )
            payload = json.loads(Path(outputs["json"]).read_text(encoding="utf-8"))

        serialized = payload["top_candidates"][0]["simulated_pattern"]
        self.assertEqual(tuple(serialized["two_theta_deg"]), direct.two_theta_deg)
        self.assertEqual(tuple(serialized["normalized_intensity"]), direct.normalized_intensity)
        self.assertEqual(payload["top_candidates"][0]["phase_entry_ids"], ["target", "duplicate"])
        self.assertEqual(payload["top_candidates"][0]["space_group_number"], 1)


if __name__ == "__main__":
    unittest.main()
