import tempfile
import unittest
from pathlib import Path

from ophiuchus.library.database import StructureLibrary
from ophiuchus.library.local_provider import LocalFolderProvider
from ophiuchus.library.providers import ProviderError, ProviderRegistry, StubProvider


SIMPLE_CIF = """
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


class LibraryProviderTests(unittest.TestCase):
    def test_provider_registry_lists_registered_provider(self):
        registry = ProviderRegistry()
        provider = StubProvider("COD", "public_api", "COD integration is not configured yet.")
        registry.register(provider)
        self.assertEqual(registry.names(), ["COD"])
        self.assertIs(registry.get("COD"), provider)

    def test_stub_provider_raises_clear_error(self):
        provider = StubProvider("ICSD_manual", "restricted_manual", "Manual import only; no scraping.")
        with self.assertRaises(ProviderError) as ctx:
            provider.search_by_elements(["Fe"], max_elements=1, include_subsystems=True, filters={})
        self.assertIn("Manual import only", str(ctx.exception))

    def test_local_folder_provider_imports_cif_into_library(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cifs"
            source.mkdir()
            (source / "Fe.cif").write_text(SIMPLE_CIF, encoding="utf-8")
            library = StructureLibrary(root / "library.sqlite")
            summary = LocalFolderProvider().import_folder(source, library)
            entries = library.list_structures()
            cached = root / entries[0].cached_file_path
            self.assertEqual(summary["imported"], 1)
            self.assertEqual(entries[0].formula, "Fe")
            self.assertEqual(entries[0].source, "local")
            self.assertEqual(entries[0].space_group_number, 221)
            self.assertEqual(entries[0].crystal_system, "cubic")
            self.assertAlmostEqual(entries[0].lattice_parameters["a"], 4.0)
            self.assertTrue(cached.exists())
            self.assertIn("User-provided", entries[0].license_or_access_note)


if __name__ == "__main__":
    unittest.main()
