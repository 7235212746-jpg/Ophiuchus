# VESTA Validation

Date: 2026-06-30

## Reference Case

- Structure: `data\library\structures\local\local_362b46a32755d09f.cif`
- Reference pattern: `C:\path\to\MONI.int`
- Radiation: Cu Kalpha1
- Wavelength: 1.54056 Angstrom
- 2theta range: 10-90 degree

## Debug Outputs

- B=0 Ophi peaks: `results\validation_reports\ZrFe6Ge4_debug_B0_peaks.csv`
- B=0 comparison: `results\validation_reports\ZrFe6Ge4_debug_B0_vs_MONI.csv`
- B=4 Ophi peaks: `results\validation_reports\ZrFe6Ge4_debug_B4_peaks.csv`
- B=4 comparison: `results\validation_reports\ZrFe6Ge4_debug_B4_vs_MONI.csv`

## Result

The peak positions agree well. The major failure was relative intensity, especially at high angle.

Without Debye-Waller correction:

- 70.07 degree: VESTA 13.3666, Ophi 25.1127, ratio 1.879
- 74.80 degree: VESTA 18.4305, Ophi 39.8468, ratio 2.162
- strong intensity mismatch rows: 3

With a uniform Debye-Waller B factor of 4.0:

- 70.07 degree: VESTA 13.3666, Ophi 13.4861, ratio 1.009
- 74.80 degree: VESTA 18.4305, Ophi 18.7408, ratio 1.017
- strong intensity mismatch rows: 0

## Interpretation

The earlier high-angle intensity problem was caused mainly by comparing an un-damped pymatgen stick pattern against a VESTA continuous/profile pattern. Adding a uniform Debye-Waller factor suppresses high-angle intensity and brings the main high-angle peaks close to the VESTA reference.

This does not prove B=4.0 is a universal physical value. It is a documented calibration setting for this reference case. Future trusted use should expose B/profile settings to the user and keep peak-position validation separate from intensity-shape validation.
