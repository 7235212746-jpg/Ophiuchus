import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ophiuchus.xrd.config import XRDConfig
from ophiuchus.xrd.validation import simulate_cif_with_config
from ophiuchus.xrd.vesta_reference import find_local_vesta_reference, find_local_vesta_references, load_local_vesta_reference_peaks


class VestaReferenceTests(unittest.TestCase):
    def test_simulation_keeps_vesta_reference_as_validation_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            xrd_root = root / "XRD"
            target = xrd_root / "ZrFe6Ge4 sum"
            target.mkdir(parents=True)
            (target / "MONI.int").write_text("35.36 100\n74.80 33\n", encoding="utf-8")
            (target / "reference_simulated_peak_positions.csv").write_text(
                "two_theta,intensity_norm\n35.36,66.36725\n74.80,22.19231\n",
                encoding="utf-8",
            )
            cif = root / "ZrFe6Ge4.cif"
            cif.write_text(
                "data_ZrFe6Ge4\n"
                "_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n"
                "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
                "loop_\n"
                "_atom_site_label\n_atom_site_type_symbol\n"
                "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
                "Zr1 Zr 0 0 0\nFe1 Fe 0.5 0.5 0.5\nGe1 Ge 0.25 0.25 0.25\n",
                encoding="utf-8",
            )

            with patch.dict("os.environ", {"OPHI_VESTA_REFERENCE_DIR": str(xrd_root), "OPHI_REQUIRE_VESTA_REFERENCE": "1"}):
                peaks = simulate_cif_with_config(cif, XRDConfig.validation_default())

        self.assertNotEqual([(round(p.two_theta, 2), round(p.intensity, 5)) for p in peaks], [(35.36, 66.36725), (74.8, 22.19231)])
        self.assertTrue(all(p.peak_id.startswith("canonical_") for p in peaks))

    def test_reference_matching_does_not_match_short_element_inside_formula(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "XRD"
            target = root / "ZrFe6Ge4 sum"
            target.mkdir(parents=True)
            (target / "MONI.int").write_text("35.36 100\n", encoding="utf-8")
            with patch.dict("os.environ", {"OPHI_VESTA_REFERENCE_DIR": str(root)}):
                self.assertIsNotNone(find_local_vesta_reference("ZrFe6Ge4"))
                self.assertIsNone(find_local_vesta_reference("Zr"))
                self.assertIsNone(find_local_vesta_reference("Fe"))
                self.assertIsNone(find_local_vesta_reference("Ge"))

    def test_reference_matching_accepts_exact_formula_int_filename(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "structure"
            root.mkdir(parents=True)
            reference = root / "Zr3V3GeSn4 VESTA.int"
            reference.write_text("19.24 100\n37.39 62\n", encoding="utf-8")
            with patch.dict("os.environ", {"OPHI_VESTA_REFERENCE_DIR": str(root)}):
                matched = find_local_vesta_reference("Zr3V3GeSn4")
                self.assertIsNotNone(matched)
                self.assertEqual(matched["pattern"], reference)
                self.assertIsNone(find_local_vesta_reference("Zr"))
                self.assertIsNone(find_local_vesta_reference("ZrVGeSn"))

    def test_reference_matching_does_not_use_nonmatching_child_int_in_formula_folder(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "XRD"
            target = root / "Zr3V3GeSn4 400-900 100h"
            target.mkdir(parents=True)
            (target / "Zr.int").write_text("20 100\n", encoding="utf-8")
            (target / "simulated_peak_positions.csv").write_text("two_theta,intensity_norm\n30,100\n", encoding="utf-8")
            with patch.dict("os.environ", {"OPHI_VESTA_REFERENCE_DIR": str(root)}):
                matched = find_local_vesta_reference("Zr3V3GeSn4")
        self.assertIsNotNone(matched)
        self.assertEqual(matched["pattern"].name, "simulated_peak_positions.csv")

    def test_reference_matching_exposes_ordered_alternatives(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "XRD"
            a = root / "Zr3V3GeSn4 900 100h"
            b = root / "Zr3V3GeSn4 950 100h"
            a.mkdir(parents=True)
            b.mkdir(parents=True)
            (a / "simulated_peak_positions.csv").write_text("two_theta,intensity_norm\n20,100\n", encoding="utf-8")
            (b / "simulated_peak_positions.csv").write_text("two_theta,intensity_norm\n30,100\n", encoding="utf-8")
            exact = root / "Zr3V3GeSn4 VESTA.int"
            exact.write_text("40 100\n", encoding="utf-8")
            with patch.dict("os.environ", {"OPHI_VESTA_REFERENCE_DIR": str(root)}):
                refs = find_local_vesta_references("Zr3V3GeSn4")
        self.assertGreaterEqual(len(refs), 3)
        self.assertEqual(refs[0]["pattern"].name, "Zr3V3GeSn4 VESTA.int")
        self.assertTrue(all("rank_reason" in item for item in refs))

    def test_reference_matching_honors_formula_specific_override(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "XRD"
            auto = root / "Zr3V3GeSn4 950 100h"
            manual = root / "manual"
            auto.mkdir(parents=True)
            manual.mkdir()
            auto_ref = auto / "simulated_peak_positions.csv"
            manual_ref = manual / "Zr3V3GeSn4 custom VESTA.int"
            auto_ref.write_text("two_theta,intensity_norm\n30,100\n", encoding="utf-8")
            manual_ref.write_text("40 100\n", encoding="utf-8")
            with patch.dict(
                "os.environ",
                {
                    "OPHI_VESTA_REFERENCE_DIR": str(root),
                    "OPHI_VESTA_REFERENCE_ZR3V3GESN4": str(manual_ref),
                },
            ):
                matched = find_local_vesta_reference("Zr3V3GeSn4")
        self.assertEqual(matched["pattern"], manual_ref)
        self.assertEqual(matched["rank_reason"], "formula-specific override")

    def test_reference_peaks_load_from_standard_peak_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "XRD"
            target = root / "ZrFe6Ge4 sum"
            target.mkdir(parents=True)
            (target / "MONI.int").write_text("35.36 100\n", encoding="utf-8")
            (target / "reference_simulated_peak_positions.csv").write_text(
                "two_theta,intensity_norm\n35.36,100\n74.80,20\n",
                encoding="utf-8",
            )
            with patch.dict("os.environ", {"OPHI_VESTA_REFERENCE_DIR": str(root)}):
                peaks, path = load_local_vesta_reference_peaks("ZrFe6Ge4", 10, 90)
        self.assertEqual([round(peak.two_theta, 2) for peak in peaks], [35.36, 74.8])
        self.assertIsNotNone(path)


if __name__ == "__main__":
    unittest.main()
