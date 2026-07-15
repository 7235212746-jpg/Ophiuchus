from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class XRDConfig:
    radiation_source: str = "CuKalpha12"
    wavelength_angstrom: float = 1.54056
    two_theta_min: float = 10.0
    two_theta_max: float = 90.0
    peak_merge_tolerance_deg: float = 0.03
    intensity_threshold: float = 1.0
    line_model: str = "cu_kalpha12"
    normalization: str = "max_100"
    angle_unit: str = "degree_2theta"
    debye_waller_b: float = 0.0

    @classmethod
    def validation_default(cls) -> "XRDConfig":
        return cls()

    def two_theta_range(self) -> tuple[float, float]:
        return (self.two_theta_min, self.two_theta_max)

    def config_id(self) -> str:
        payload = json.dumps(asdict(self), sort_keys=True)
        digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:10]
        return f"{self.radiation_source}_{self.wavelength_angstrom:.5f}_{digest}"

    def to_dict(self) -> dict[str, object]:
        return asdict(self)
