from __future__ import annotations

from .models import PatternPoint, Peak


def find_peaks(
    points: list[PatternPoint] | list[tuple[float, float]],
    two_theta_min: float = 10.0,
    two_theta_max: float = 90.0,
    min_height: float = 1.0,
    min_distance_deg: float = 0.15,
    max_peaks: int = 80,
    smooth_window: int = 1,
    min_prominence: float = 0.0,
) -> list[Peak]:
    rows = [
        (p.two_theta, p.intensity) if isinstance(p, PatternPoint) else (float(p[0]), float(p[1]))
        for p in points
    ]
    rows = [(x, y) for x, y in rows if two_theta_min <= x <= two_theta_max]
    if len(rows) < 3:
        return []
    xs = [x for x, _ in rows]
    ys = _smooth([y for _, y in rows], smooth_window)
    selected: list[int]
    prominence_by_index: dict[int, float] = {}
    try:
        import numpy as np
        from scipy.signal import find_peaks as scipy_find_peaks

        steps = [xs[index + 1] - xs[index] for index in range(len(xs) - 1) if xs[index + 1] > xs[index]]
        median_step = float(np.median(steps)) if steps else min_distance_deg
        distance_points = max(1, int(round(min_distance_deg / max(median_step, 1e-12))))
        indices, properties = scipy_find_peaks(
            np.asarray(ys, dtype=float),
            height=min_height,
            prominence=max(0.0, min_prominence),
            distance=distance_points,
        )
        prominence_by_index = {
            int(index): float(prominence)
            for index, prominence in zip(indices, properties.get("prominences", []))
        }
        selected = sorted((int(index) for index in indices), key=lambda index: ys[index], reverse=True)[:max_peaks]
    except Exception:
        candidates: list[int] = []
        for i in range(1, len(rows) - 1):
            local_prominence = max(0.0, ys[i] - max(ys[i - 1], ys[i + 1]))
            if (
                ys[i] >= min_height
                and local_prominence >= min_prominence
                and ys[i] > ys[i - 1]
                and ys[i] >= ys[i + 1]
            ):
                candidates.append(i)
                prominence_by_index[i] = local_prominence
        candidates.sort(key=lambda i: ys[i], reverse=True)
        selected = []
        for idx in candidates:
            if all(abs(xs[idx] - xs[old]) >= min_distance_deg for old in selected):
                selected.append(idx)
            if len(selected) >= max_peaks:
                break
    selected.sort(key=lambda i: xs[i])
    peaks: list[Peak] = []
    for serial, idx in enumerate(selected, 1):
        prominence = prominence_by_index.get(idx)
        if prominence is None:
            left = ys[idx - 1] if idx > 0 else ys[idx]
            right = ys[idx + 1] if idx < len(ys) - 1 else ys[idx]
            prominence = max(0.0, ys[idx] - max(left, right))
        peaks.append(Peak(f"exp_{serial}", xs[idx], rows[idx][1], prominence=prominence))
    return peaks


def _smooth(values: list[float], window: int) -> list[float]:
    if window <= 1 or len(values) < 3:
        return values[:]
    if window % 2 == 0:
        window += 1
    radius = window // 2
    smoothed: list[float] = []
    for i in range(len(values)):
        lo = max(0, i - radius)
        hi = min(len(values), i + radius + 1)
        smoothed.append(sum(values[lo:hi]) / (hi - lo))
    return smoothed
