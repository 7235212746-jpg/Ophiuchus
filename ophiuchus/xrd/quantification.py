from __future__ import annotations

from collections.abc import Iterable

import numpy as np


def weight_fractions_from_zmv(
    rows: Iterable[tuple[str, float, float, float, float]],
) -> dict[str, float]:
    weighted: dict[str, float] = {}
    for phase_id, scale, z, molar_mass, volume in rows:
        values = np.asarray([scale, z, molar_mass, volume], dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError("Every quantitative phase requires positive finite scale, Z, M, and V values.")
        if phase_id in weighted:
            raise ValueError(f"Duplicate quantitative phase id: {phase_id}")
        weighted[phase_id] = float(np.prod(values))
    if not weighted:
        raise ValueError("At least one quantitative phase is required.")
    total = sum(weighted.values())
    return {phase_id: value / total * 100.0 for phase_id, value in weighted.items()}
