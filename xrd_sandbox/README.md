# XRD Sandbox

This folder is an Ophi-independent XRD reproduction workbench.

It is for debugging the scientific calculation before the logic is moved back into Ophi.

## Run

Put CIF files here:

```text
.\xrd_sandbox\inputs
```

Example target file:

```text
.\xrd_sandbox\inputs\Zr3V3GeSn4.cif
```

```powershell
python -m xrd_sandbox.reproduce_xrd `
  --cif ".\xrd_sandbox\inputs\Zr3V3GeSn4.cif" `
  --reference "C:\path\to\Zr3V3GeSn4 VESTA.int" `
  --out ".\xrd_sandbox\runs\Zr3V3GeSn4" `
  --two-theta-min 10 `
  --two-theta-max 90 `
  --wavelength 1.54056 `
  --tolerance 0.10
```

## Outputs

- `python_simulated_peaks.csv`: direct pymatgen CIF simulation.
- `vesta_reference_peaks.csv`: parsed reference peaks from the VESTA/reference file.
- `peak_comparison.csv`: matched, missing, extra peaks, and intensity ratios.
- `summary.json`: configuration and summary metrics.
- `comparison_plot.png`: reference peaks and Python simulated peaks in a two-panel figure.

## Rule

Do not treat this as an Ophi feature yet. First use it to compare peak positions and intensities against trusted VESTA output, then port only the validated correction back into Ophi.
