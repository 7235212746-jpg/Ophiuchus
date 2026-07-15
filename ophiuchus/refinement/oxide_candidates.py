from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from pymatgen.core import Composition

from ophiuchus.library.database import StructureLibrary
from ophiuchus.library.models import StructureEntry
from ophiuchus.library.mp_provider import MaterialsProjectProvider
from ophiuchus.library.xrd_cache import (
    build_library_xrd_cache,
    library_entries_to_candidates,
    load_validated_pattern,
)
from ophiuchus.xrd.models import Candidate


COMMON_OXIDES_BY_ELEMENT = {
    "Fe": ("FeO", "Fe2O3", "Fe3O4"),
    "Zr": ("ZrO2",),
    "Ge": ("GeO2",),
}

COMMON_OXIDE_SPACE_GROUP_NUMBERS = {
    "FeO": (225,),
    "Fe2O3": (167, 206),
    "Fe3O4": (227,),
    "ZrO2": (14, 137, 225),
    "GeO2": (136, 152, 154),
}

# MP's corrected 0 K hull strongly penalizes rocksalt FeO, while wustite is a
# common experimental impurity. The exception remains prototype-restricted.
COMMON_OXIDE_MAX_HULL = {"FeO": 0.50}


@dataclass(frozen=True)
class ControlledOxideLoadResult:
    candidates: tuple[Candidate, ...]
    formulas_checked: tuple[str, ...]
    warnings: tuple[str, ...]


def select_controlled_oxide_entries(
    entries: Iterable[StructureEntry],
    main_elements: set[str],
    *,
    max_polymorphs_per_formula: int = 3,
    maximum_energy_above_hull: float = 0.12,
) -> list[StructureEntry]:
    main = {element for element in main_elements if element and element != "O"}
    formula_order = list(expected_common_oxide_formulas(main))
    allowed_formulas = set(formula_order)
    grouped: dict[str, list[StructureEntry]] = {formula: [] for formula in formula_order}
    for item in entries:
        item_elements = set(item.elements)
        if not item.enabled_for_matching or "O" not in item_elements:
            continue
        if not item_elements.issubset(main | {"O"}) or len(item_elements - {"O"}) != 1:
            continue
        formula = _reduced_formula(item.reduced_formula or item.formula)
        if formula not in allowed_formulas:
            continue
        preferred_space_groups = COMMON_OXIDE_SPACE_GROUP_NUMBERS[formula]
        if item.space_group_number not in preferred_space_groups and item.is_experimental is not True:
            continue
        allowed_hull = max(maximum_energy_above_hull, COMMON_OXIDE_MAX_HULL.get(formula, maximum_energy_above_hull))
        if item.energy_above_hull is not None and item.energy_above_hull > allowed_hull:
            continue
        grouped[formula].append(item)

    selected: list[StructureEntry] = []
    limit = max(1, int(max_polymorphs_per_formula))
    for formula in formula_order:
        preferred_space_groups = COMMON_OXIDE_SPACE_GROUP_NUMBERS[formula]
        ordered = sorted(grouped[formula], key=lambda item: _entry_priority(item, preferred_space_groups))
        seen_space_groups: set[str] = set()
        for item in ordered:
            space_group = (item.space_group_symbol or f"number:{item.space_group_number}" or item.internal_id).strip()
            if space_group in seen_space_groups:
                continue
            seen_space_groups.add(space_group)
            selected.append(item)
            if len(seen_space_groups) >= limit:
                break
    return selected


def load_controlled_oxide_candidates(
    library_path: str | Path,
    main_elements: set[str],
    *,
    radiation: str,
    two_theta_range: tuple[float, float],
) -> ControlledOxideLoadResult:
    library = StructureLibrary(library_path)
    entries = select_controlled_oxide_entries(library.list_structures(enabled_only=True), main_elements)
    selected_formulas = tuple(dict.fromkeys(_reduced_formula(item.reduced_formula or item.formula) for item in entries))
    missing = [
        item.internal_id
        for item in entries
        if load_validated_pattern(
            library,
            item,
            radiation=radiation,
            two_theta_range=two_theta_range,
        )
        is None
    ]
    warnings = [
        f"{formula}：结构库缺少常见晶型，未参与匹配。"
        for formula in expected_common_oxide_formulas(set(main_elements))
        if formula not in selected_formulas
    ]
    if missing:
        summary = build_library_xrd_cache(
            library,
            structure_ids=missing,
            radiation=radiation,
            two_theta_range=two_theta_range,
            scientific_safe_mode=True,
        )
        warnings.extend(str(item) for item in summary.get("warnings", []))
    ids = {item.internal_id for item in entries}
    candidates = library_entries_to_candidates(
        library,
        list(main_elements | {"O"}),
        radiation=radiation,
        two_theta_range=two_theta_range,
        candidate_ids=ids,
    )
    entries_by_id = {item.internal_id: item for item in entries}
    for candidate in candidates:
        item = entries_by_id[candidate.candidate_id]
        candidate.space_group_symbol = item.space_group_symbol or ""
        candidate.space_group_number = item.space_group_number
        candidate.simulation_state = "validated"
    loaded_ids = {candidate.candidate_id for candidate in candidates}
    for item in entries:
        if item.internal_id not in loaded_ids:
            warnings.append(f"{item.reduced_formula or item.formula} ({item.internal_id}) 未生成可验证模拟谱。")
    formulas_checked = tuple(dict.fromkeys(candidate.formula_pretty for candidate in candidates))
    return ControlledOxideLoadResult(tuple(candidates), formulas_checked, tuple(dict.fromkeys(warnings)))


def expected_common_oxide_formulas(main_elements: set[str]) -> tuple[str, ...]:
    main = {element for element in main_elements if element and element != "O"}
    return tuple(
        formula
        for element in ("Fe", "Zr", "Ge")
        if element in main
        for formula in COMMON_OXIDES_BY_ELEMENT[element]
    )


def missing_common_oxide_requirements(
    library: StructureLibrary,
    main_elements: set[str],
) -> dict[str, tuple[int, ...]]:
    entries = library.list_structures(enabled_only=True)
    requirements: dict[str, tuple[int, ...]] = {}
    for formula in expected_common_oxide_formulas(main_elements):
        allowed_space_groups = COMMON_OXIDE_SPACE_GROUP_NUMBERS[formula]
        present = any(
            _reduced_formula(item.reduced_formula or item.formula) == formula
            and item.space_group_number in allowed_space_groups
            for item in entries
        )
        if not present:
            requirements[formula] = allowed_space_groups
    return requirements


def supplement_common_oxide_library(
    library_path: str | Path,
    main_elements: set[str],
    *,
    provider: MaterialsProjectProvider,
    radiation: str = "CuKa",
    two_theta_range: tuple[float, float] = (10.0, 90.0),
    progress_callback=None,
) -> dict[str, object]:
    library = StructureLibrary(library_path)
    requirements = missing_common_oxide_requirements(library, main_elements)
    if requirements:
        harvest = provider.harvest_formulas_to_library(
            library,
            requirements,
            maximum_energy_above_hull_by_formula={
                formula: max(0.12, COMMON_OXIDE_MAX_HULL.get(formula, 0.12))
                for formula in requirements
            },
            progress_callback=progress_callback,
        )
    else:
        harvest = {
            "provider": "materials_project",
            "retrieved": 0,
            "imported": 0,
            "imported_ids": [],
            "skipped_duplicates": 0,
            "skipped_prototype": 0,
            "failed": 0,
            "warnings": [],
        }
    imported_ids = [str(item) for item in harvest.get("imported_ids", [])]
    cache_summary: dict[str, object] = {"checked": 0, "simulated": 0, "failed": 0, "warnings": []}
    if imported_ids:
        cache_summary = build_library_xrd_cache(
            library,
            structure_ids=imported_ids,
            radiation=radiation,
            two_theta_range=two_theta_range,
            scientific_safe_mode=True,
        )
    missing_after = missing_common_oxide_requirements(library, main_elements)
    return {
        **harvest,
        "requested_formulas": list(requirements),
        "missing_after": list(missing_after),
        "xrd_cache": cache_summary,
    }


def _reduced_formula(formula: str) -> str:
    try:
        return Composition(formula).reduced_formula
    except Exception:
        return formula.replace(" ", "")


def _entry_priority(item: StructureEntry, preferred_space_groups: tuple[int, ...]) -> tuple[float, float, float, str]:
    try:
        prototype_priority = float(preferred_space_groups.index(int(item.space_group_number)))
    except (TypeError, ValueError):
        prototype_priority = float(len(preferred_space_groups))
    experimental_priority = 0.0 if item.is_experimental is True else 1.0
    hull = float(item.energy_above_hull) if item.energy_above_hull is not None else 0.08
    return prototype_priority, experimental_priority, hull, item.internal_id
