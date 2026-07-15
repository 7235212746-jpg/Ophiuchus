import csv
import json
import tempfile
import unittest
from pathlib import Path

from ophiuchus.library.database import StructureLibrary
from ophiuchus.library.exports import export_analysis_bundle, next_export_path, write_analysis_json
from ophiuchus.xrd.models import Candidate, CandidateScore, MultiPhaseExplanation, Peak, PeakMatch
from ophiuchus.xrd.pipeline import AnalysisResult
from ophiuchus.library.inspector import inspect_peak
from ophiuchus.library.models import LibraryPeak, StructureEntry


class LibraryExportsInspectorTests(unittest.TestCase):
    def test_next_export_path_auto_increments(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            first = next_export_path(folder, "ZrFe6Ge4", "ZFG_980C72h", "20260627", "XRD_match", "png")
            first.write_text("old", encoding="utf-8")
            second = next_export_path(folder, "ZrFe6Ge4", "ZFG_980C72h", "20260627", "XRD_match", "png")
        self.assertTrue(str(first).endswith("ZrFe6Ge4_ZFG_980C72h_20260627_XRD_match_v01.png"))
        self.assertTrue(str(second).endswith("ZrFe6Ge4_ZFG_980C72h_20260627_XRD_match_v02.png"))

    def test_write_analysis_json_keeps_traceability_and_caveat(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "analysis.json"
            write_analysis_json(
                out,
                project="P",
                sample_id="S",
                settings={"radiation": "CuKa"},
                selected_library_entry_ids=["local_1"],
                candidate_phases=[{"formula": "Fe"}],
            )
            payload = json.loads(out.read_text(encoding="utf-8"))
        self.assertEqual(payload["project"], "P")
        self.assertEqual(payload["selected_library_entry_ids"], ["local_1"])
        self.assertIn("not Rietveld", payload["scientific_caveat"])

    def test_export_analysis_bundle_writes_named_png_csv_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_png = root / "source.png"
            source_png.write_bytes(b"fakepng")
            source_diag = root / "diagnostic.png"
            source_diag.write_bytes(b"fakediag")
            candidate = Candidate("local_fe", "Fe", "library:local", "Fe.cif", ["Fe"], "hash")
            theory = Peak("t1", 22.2, 100.0)
            exp = Peak("e1", 22.25, 90.0)
            score = CandidateScore(candidate, 1.0, [PeakMatch(theory, exp, 0.05)], [], [exp], [])
            result = AnalysisResult(
                experimental_peaks=[exp],
                candidates=[candidate],
                top_scores=[score],
                explanation=MultiPhaseExplanation([score], [exp], []),
                outputs={"xrd_presentation_plot": str(source_png), "xrd_diagnostic_plot": str(source_diag)},
                warnings=[],
            )
            outputs = export_analysis_bundle(result, root / "exports", "FeProject", "FeSample", date_text="20260629")
            with Path(outputs["peak_table"]).open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            payload = json.loads(Path(outputs["analysis_report"]).read_text(encoding="utf-8"))
        self.assertTrue(outputs["xrd_presentation"].endswith("FeProject_FeSample_20260629_xrd_candidate_screening_presentation_v01.png"))
        self.assertTrue(outputs["xrd_diagnostic"].endswith("FeProject_FeSample_20260629_xrd_candidate_screening_diagnostic_v01.png"))
        self.assertTrue(outputs["peak_table"].endswith("FeProject_FeSample_20260629_peak_table_v01.csv"))
        self.assertTrue(outputs["analysis_report"].endswith("FeProject_FeSample_20260629_analysis_report_v01.json"))
        self.assertEqual(rows[0]["candidate_phase"], "Fe")
        self.assertEqual(payload["selected_library_entry_ids"], ["local_fe"])

    def test_inspect_peak_returns_nearby_local_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            library = StructureLibrary(Path(tmp) / "library.sqlite")
            entry = StructureEntry(
                internal_id="local_fe",
                source="local",
                source_id="Fe.cif",
                formula="Fe",
                reduced_formula="Fe",
                elements=["Fe"],
                cached_file_path="library/structures/local/local_fe.cif",
                structure_hash="abc",
            )
            library.upsert_structure(entry)
            library.store_xrd_peaks(
                "local_fe",
                [LibraryPeak("local_fe", 22.20, 100.0, "CuKa", "settings1", peak_id="p1", hkl="100")],
                "CuKa",
                "settings1",
                10.0,
                90.0,
            )
            evidence = inspect_peak(library, 22.25, tolerance_deg=0.1, radiation="CuKa", settings_hash="settings1")
        self.assertEqual(evidence["status"], "matched")
        self.assertEqual(evidence["nearby_peaks"][0]["formula"], "Fe")
        self.assertAlmostEqual(evidence["nearby_peaks"][0]["delta"], 0.05)


if __name__ == "__main__":
    unittest.main()
