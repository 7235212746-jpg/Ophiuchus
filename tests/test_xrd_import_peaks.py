import random
import tempfile
import unittest
from pathlib import Path

from ophiuchus.xrd.importers import infer_xrd_range, load_xrd_file, normalize_pattern
from ophiuchus.xrd.peaks import find_peaks


class XrdImportPeakTests(unittest.TestCase):
    def test_reads_rigaku_asc_metadata_and_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.asc"
            path.write_text(
                "\n".join(
                    [
                        "*SAMPLE = ZrFe6Ge4",
                        "*DATE = 09-Jun-26 17:54",
                        "*START = 3",
                        "*STOP = 3.06",
                        "*STEP = 0.02",
                        "*COUNT = 4",
                        "10, 20, 15, 5",
                    ]
                ),
                encoding="utf-8",
            )
            pattern = load_xrd_file(path)
        self.assertEqual(pattern.metadata["sample"], "ZrFe6Ge4")
        self.assertEqual(len(pattern.points), 4)
        self.assertAlmostEqual(pattern.points[-1].two_theta, 3.06)

    def test_reads_two_column_text_with_comments_and_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "pattern.xy"
            path.write_text(
                "# comment\n2theta intensity\n10 5\n10.5 15\n11 10\n",
                encoding="utf-8",
            )
            pattern = load_xrd_file(path, normalize=True)
        self.assertEqual([p.two_theta for p in pattern.points], [10.0, 10.5, 11.0])
        self.assertAlmostEqual(max(p.intensity for p in pattern.points), 100.0)

    def test_infer_xrd_range_from_folder_combines_readable_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "a.xy").write_text("10 1\n20 2\n", encoding="utf-8")
            (root / "b.csv").write_text("30,1\n45,2\n", encoding="utf-8")
            (root / "bad.xy").write_text("not xrd", encoding="utf-8")
            info = infer_xrd_range(root, padding_deg=0.25)
        self.assertAlmostEqual(info["two_theta_min"], 9.75)
        self.assertAlmostEqual(info["two_theta_max"], 45.25)
        self.assertEqual(len(info["files_read"]), 2)
        self.assertTrue(info["warnings"])

    def test_normalize_handles_zero_pattern(self):
        points = normalize_pattern([(10.0, 0.0), (11.0, 0.0)])
        self.assertEqual([p.intensity for p in points], [0.0, 0.0])

    def test_randomized_peak_recovery_without_scipy(self):
        rng = random.Random(20260627)
        for _ in range(12):
            centers = sorted(rng.sample([20 + i for i in range(45)], 4))
            points = []
            for i in range(3000):
                x = 10.0 + i * 0.02
                y = rng.uniform(0.0, 0.6)
                for center in centers:
                    y += 100.0 * pow(2.718281828, -((x - center) ** 2) / (2 * 0.05**2))
                points.append((x, y))
            peaks = find_peaks(points, min_height=10.0, min_distance_deg=0.5, max_peaks=20)
            recovered = [p.two_theta for p in peaks]
            for center in centers:
                self.assertTrue(
                    any(abs(x - center) <= 0.08 for x in recovered),
                    f"missing generated peak near {center}; recovered={recovered}",
                )


if __name__ == "__main__":
    unittest.main()
