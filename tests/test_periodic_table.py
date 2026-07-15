from pathlib import Path
import unittest

from ophiuchus.periodic_table import (
    ELEMENTS,
    element_scope_mismatch,
    infer_elements_from_xrd_path,
    parse_element_symbols,
)


class PeriodicTableTests(unittest.TestCase):
    def test_contains_all_118_elements_in_atomic_order_with_display_positions(self):
        self.assertEqual(len(ELEMENTS), 118)
        self.assertEqual([element.atomic_number for element in ELEMENTS], list(range(1, 119)))
        self.assertEqual(ELEMENTS[0].symbol, "H")
        self.assertEqual(ELEMENTS[-1].symbol, "Og")
        self.assertTrue(all(element.display_row >= 1 and element.display_column >= 1 for element in ELEMENTS))

    def test_parse_normalizes_case_spacing_and_rejects_invalid_symbols(self):
        self.assertEqual(parse_element_symbols(" zr, FE   ge sn "), ("Zr", "Fe", "Ge", "Sn"))
        self.assertEqual(parse_element_symbols("Ge Zr Ge"), ("Ge", "Zr"))
        with self.assertRaisesRegex(ValueError, "invalid element"):
            parse_element_symbols("Zr Xx Ge")

    def test_formula_inference_prefers_file_then_parent_and_rejects_generic_names(self):
        self.assertEqual(
            infer_elements_from_xrd_path(Path("C:/XRD/other/Zr3V3GeSn4.asc")),
            ("Zr", "V", "Ge", "Sn"),
        )
        self.assertEqual(
            infer_elements_from_xrd_path(Path("C:/XRD/ZrFe6Ge4 in different Temperature/980.asc")),
            ("Zr", "Fe", "Ge"),
        )
        self.assertEqual(infer_elements_from_xrd_path(Path("C:/XRD/run/ARC.asc")), ())
        self.assertEqual(infer_elements_from_xrd_path(Path("C:/XRD/run/sample.asc")), ())

    def test_scope_check_blocks_only_missing_inferred_elements(self):
        missing, extra = element_scope_mismatch(("Zr", "Ge"), ("Zr", "V", "Ge"))
        self.assertEqual(missing, ("V",))
        self.assertEqual(extra, ())

        missing, extra = element_scope_mismatch(("Zr", "V", "Ge", "Sn"), ("Zr", "V", "Ge"))
        self.assertEqual(missing, ())
        self.assertEqual(extra, ("Sn",))


if __name__ == "__main__":
    unittest.main()
