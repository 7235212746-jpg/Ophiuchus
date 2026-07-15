from __future__ import annotations

from pathlib import Path

from .config import XRDConfig
from .models import Peak


def simulate_cif_with_pymatgen(cif_file: str | Path, config: XRDConfig) -> list[Peak]:
    from .backend import SimulationContext, ValidatedXRDBackend

    path = Path(cif_file)
    return ValidatedXRDBackend().simulate_cif(
        path,
        config,
        SimulationContext(structure_id=path.stem, source="pymatgen_compatibility"),
    ).to_peaks()
