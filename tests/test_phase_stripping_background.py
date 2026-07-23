import unittest

import numpy as np

from ophiuchus.phase_stripping.background import estimate_xrd_background


class PhaseStrippingBackgroundTests(unittest.TestCase):
    def test_asls_recovers_smooth_background_without_absorbing_narrow_peaks(self):
        x = np.arange(10.0, 80.001, 0.02)
        expected_background = 180.0 + 0.8 * (x - 10.0)
        expected_background += 55.0 * np.exp(-0.5 * ((x - 18.0) / 5.0) ** 2)
        peaks = 900.0 * np.exp(-0.5 * ((x - 32.0) / 0.07) ** 2)
        peaks += 500.0 * np.exp(-0.5 * ((x - 61.0) / 0.10) ** 2)

        estimate = estimate_xrd_background(x, expected_background + peaks)

        self.assertEqual(estimate.method, "asls")
        self.assertEqual(estimate.values.shape, x.shape)
        self.assertFalse(estimate.values.flags.writeable)
        self.assertLess(float(np.sqrt(np.mean((estimate.values - expected_background) ** 2))), 12.0)
        corrected = expected_background + peaks - estimate.values
        self.assertGreater(float(corrected[np.argmin(abs(x - 32.0))]), 850.0)

    def test_estimator_rejects_invalid_or_nonfinite_inputs(self):
        with self.assertRaisesRegex(ValueError, "same shape"):
            estimate_xrd_background(np.array([1.0, 2.0]), np.array([1.0]))
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            estimate_xrd_background(np.array([1.0, 1.0]), np.array([2.0, 3.0]))
        with self.assertRaisesRegex(ValueError, "finite"):
            estimate_xrd_background(np.array([1.0, 2.0]), np.array([2.0, np.nan]))


if __name__ == "__main__":
    unittest.main()
