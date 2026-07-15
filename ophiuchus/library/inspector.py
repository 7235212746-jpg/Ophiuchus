from __future__ import annotations

import csv
import json
from pathlib import Path

from .database import StructureLibrary


def inspect_peak(
    library: StructureLibrary,
    experimental_two_theta: float,
    tolerance_deg: float = 0.2,
    enabled_only: bool = True,
    radiation: str = "CuKa",
    settings_hash: str = "",
) -> dict[str, object]:
    nearby: list[dict[str, object]] = []
    for entry in library.list_structures(enabled_only=enabled_only):
        peaks = library.load_xrd_peaks(entry.internal_id, radiation, settings_hash)
        strong_peaks = sorted([peak for peak in peaks if peak.relative_intensity >= 10.0], key=lambda peak: peak.relative_intensity, reverse=True)
        for peak in peaks:
            delta = abs(peak.two_theta - experimental_two_theta)
            if delta <= tolerance_deg:
                nearby.append(
                    {
                        "structure_internal_id": entry.internal_id,
                        "formula": entry.reduced_formula or entry.formula,
                        "source": entry.source,
                        "source_id": entry.source_id,
                        "structure_hash": entry.structure_hash,
                        "enabled": entry.enabled_for_matching,
                        "cif_path": str(Path(library.path.parent) / entry.cached_file_path),
                        "calculated_two_theta": peak.two_theta,
                        "experimental_two_theta": experimental_two_theta,
                        "delta": round(delta, 4),
                        "hkl": peak.hkl,
                        "theoretical_relative_intensity": peak.relative_intensity,
                        "strong_theoretical_peak_count": len(strong_peaks),
                        "top_strong_theoretical_peaks": [
                            {"two_theta": round(item.two_theta, 4), "intensity": round(item.relative_intensity, 3), "hkl": item.hkl or ""}
                            for item in strong_peaks[:8]
                        ],
                    }
                )
    nearby.sort(key=lambda item: (item["delta"], -float(item["theoretical_relative_intensity"])))
    warnings = []
    if len({item["structure_internal_id"] for item in nearby}) > 1:
        warnings.append("overlap warning: multiple library phases have nearby simulated peaks")
    return {
        "status": "matched" if nearby else "unresolved",
        "experimental_two_theta": experimental_two_theta,
        "tolerance_deg": tolerance_deg,
        "nearby_peaks": nearby,
        "warnings": warnings,
        "scientific_note": "Peak evidence is local cached screening evidence, not a definitive phase assignment.",
    }


def export_peak_inspection(evidence: dict[str, object], out_base: str | Path) -> dict[str, str]:
    base = Path(out_base)
    base.parent.mkdir(parents=True, exist_ok=True)
    json_path = base.with_suffix(".json")
    csv_path = base.with_suffix(".csv")
    json_path.write_text(json.dumps(evidence, indent=2, ensure_ascii=False), encoding="utf-8")
    rows = list(evidence.get("nearby_peaks") or [])
    fieldnames = [
        "structure_internal_id",
        "formula",
        "source",
        "source_id",
        "structure_hash",
        "enabled",
        "cif_path",
        "experimental_two_theta",
        "calculated_two_theta",
        "delta",
        "hkl",
        "theoretical_relative_intensity",
        "strong_theoretical_peak_count",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    return {"json": str(json_path), "csv": str(csv_path)}
