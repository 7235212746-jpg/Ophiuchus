from __future__ import annotations

from itertools import combinations


def _clean(elements: list[str] | None) -> list[str]:
    return sorted({element.strip() for element in elements or [] if element and element.strip()})


def _system(items: tuple[str, ...] | list[str]) -> str:
    return "-".join(sorted(items))


def _append_unique(target: list[str], value: str, max_systems: int) -> None:
    if value not in target and len(target) < max_systems:
        target.append(value)


def generate_chemical_systems(
    target_elements: list[str],
    impurity_elements: list[str] | None = None,
    mode: str = "normal",
    max_systems: int = 80,
) -> list[str]:
    targets = _clean(target_elements)
    impurities = _clean(impurity_elements)
    if mode not in {"conservative", "normal", "broad"}:
        raise ValueError("mode must be conservative, normal, or broad")
    systems: list[str] = []
    for element in targets:
        _append_unique(systems, element, max_systems)
    for size in (2, 3):
        for combo in combinations(targets, size):
            _append_unique(systems, _system(combo), max_systems)
    if mode in {"normal", "broad"}:
        for target in targets:
            for impurity in impurities:
                _append_unique(systems, f"{target}-{impurity}", max_systems)
    if mode == "broad":
        for impurity in impurities:
            for pair in combinations(targets, 2):
                _append_unique(systems, _system((*pair, impurity)), max_systems)
    return systems
