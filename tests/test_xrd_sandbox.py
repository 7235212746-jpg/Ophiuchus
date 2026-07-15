import csv
import tempfile
import unittest
from pathlib import Path

from xrd_sandbox.reproduce_xrd import run_xrd_reproduction, simulate_with_pymatgen


class XrdSandboxTests(unittest.TestCase):
    def test_run_xrd_reproduction_writes_debug_tables_and_plot(self):
        cif_text = """
data_Fe
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a 2.866
_cell_length_b 2.866
_cell_length_c 2.866
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0 0 0 1
Fe2 Fe 0.5 0.5 0.5 1
"""
        reference_csv = "two_theta,intensity_norm,hkl\n44.67,100,110\n65.02,20,200\n82.33,30,211\n"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cif = root / "Fe.cif"
            ref = root / "Fe_reference.csv"
            cif.write_text(cif_text, encoding="utf-8")
            ref.write_text(reference_csv, encoding="utf-8")
            outputs = run_xrd_reproduction(cif, ref, root / "out", two_theta_min=30, two_theta_max=90)
            for key in ["python_peaks_csv", "reference_peaks_csv", "comparison_csv", "summary_json", "plot_png"]:
                self.assertTrue(Path(outputs[key]).is_file(), key)
            self.assertGreater(Path(outputs["plot_png"]).stat().st_size, 1000)
            with Path(outputs["comparison_csv"]).open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        statuses = {row["status"] for row in rows}
        self.assertIn("matched", statuses)

    def test_run_xrd_reproduction_writes_profile_comparison_for_continuous_reference(self):
        cif_text = """
data_Fe
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a 2.866
_cell_length_b 2.866
_cell_length_c 2.866
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0 0 0 1
Fe2 Fe 0.5 0.5 0.5 1
"""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cif = root / "Fe.cif"
            ref = root / "Fe_reference.int"
            cif.write_text(cif_text, encoding="utf-8")
            rows = ["GENERAL$", "7000"]
            for i in range(3000, 9001):
                x = i / 100.0
                y = 100.0 if abs(x - 44.67) < 0.02 else 0.0
                rows.append(f"{x:.2f} {y:.5f} 0.0")
            ref.write_text("\n".join(rows), encoding="utf-8")
            outputs = run_xrd_reproduction(cif, ref, root / "out", two_theta_min=30, two_theta_max=90)
            self.assertTrue(Path(outputs["profile_comparison_csv"]).is_file())
            self.assertTrue(Path(outputs["profile_plot_png"]).is_file())

    def test_cu_kalpha12_line_model_adds_high_angle_doublet(self):
        cif_text = """
data_Fe
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a 2.866
_cell_length_b 2.866
_cell_length_c 2.866
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Fe1 Fe 0 0 0 1
Fe2 Fe 0.5 0.5 0.5 1
"""
        with tempfile.TemporaryDirectory() as tmp:
            cif = Path(tmp) / "Fe.cif"
            cif.write_text(cif_text, encoding="utf-8")
            kalpha1 = simulate_with_pymatgen(cif, (30, 90), wavelength=1.54056, line_model="kalpha1")
            doublet = simulate_with_pymatgen(cif, (30, 90), wavelength=1.54056, line_model="cu_kalpha12")
        self.assertGreater(len(doublet), len(kalpha1))
        high_angle = [peak.two_theta for peak in doublet if 82.0 <= peak.two_theta <= 82.7]
        self.assertGreaterEqual(len(high_angle), 2)


if __name__ == "__main__":
    unittest.main()
