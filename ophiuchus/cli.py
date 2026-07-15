from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .app import launch_app
from .library.analysis import run_library_analysis
from .library.database import StructureLibrary
from .library.exports import export_analysis_bundle
from .library.local_provider import LocalFolderProvider
from .library.mp_provider import MaterialsProjectProvider
from .library.providers import ProviderError
from .library.systems import generate_chemical_systems
from .library.xrd_cache import build_library_xrd_cache
from .xrd.cache import CandidateCache, build_candidate_cache
from .xrd.config import XRDConfig
from .xrd.debug import compare_debug_peak_tables, simulate_cif_debug_table
from .xrd.pipeline import run_analysis
from .xrd.validation import ValidationRunner, simulate_cif_with_config


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="ophiuchus", description="Ophiuchus local research workflow tools")
    sub = parser.add_subparsers(dest="command")

    analyze = sub.add_parser("analyze", help="Run Ophi XRD Candidate Screener")
    analyze.add_argument("--xrd", required=True, help="Experimental XRD file")
    analyze.add_argument("--cif-dir", action="append", required=True, help="Local CIF or peak-list directory; repeatable")
    analyze.add_argument("--elements", nargs="+", required=True, help="Main elements, such as Zr Fe Ge")
    analyze.add_argument("--extra-elements", nargs="*", default=[], help="Optional impurity/contamination elements")
    analyze.add_argument("--radiation", default="CuKa")
    analyze.add_argument("--two-theta-min", type=float, default=10.0)
    analyze.add_argument("--two-theta-max", type=float, default=90.0)
    analyze.add_argument("--tolerance", type=float, default=0.20)
    analyze.add_argument("--max-candidates", type=int, default=20)
    analyze.add_argument("--max-phases", type=int, default=3)
    analyze.add_argument("--cache-path", help="Optional local SQLite cache for simulated candidate patterns")
    analyze.add_argument("--no-cache", action="store_true", help="Do not read existing cached patterns")
    analyze.add_argument("--out-dir", required=True)

    cache_build = sub.add_parser("cache-build", help="Build or update the local XRD candidate cache")
    cache_build.add_argument("--cif-dir", action="append", required=True, help="Local CIF or peak-list directory; repeatable")
    cache_build.add_argument("--elements", nargs="+", required=True, help="Allowed elements for candidate filtering")
    cache_build.add_argument("--radiation", default="CuKa")
    cache_build.add_argument("--two-theta-min", type=float, default=10.0)
    cache_build.add_argument("--two-theta-max", type=float, default=90.0)
    cache_build.add_argument("--cache-path", required=True, help="SQLite cache path")

    library_import = sub.add_parser("library-import", help="Import local CIF files into the Ophi structure library")
    library_import.add_argument("--folder", required=True, help="Folder containing legal local CIF files")
    library_import.add_argument("--library-db", required=True, help="Ophi structure library SQLite path")

    library_list = sub.add_parser("library-list", help="List local structure library entries")
    library_list.add_argument("--library-db", required=True, help="Ophi structure library SQLite path")

    library_cache = sub.add_parser("library-cache-xrd", help="Build simulated XRD cache for library entries")
    library_cache.add_argument("--library-db", required=True, help="Ophi structure library SQLite path")
    library_cache.add_argument("--radiation", default="CuKa")
    library_cache.add_argument("--two-theta-min", type=float, default=10.0)
    library_cache.add_argument("--two-theta-max", type=float, default=90.0)

    library_analyze = sub.add_parser("library-analyze", help="Analyze experimental XRD using enabled local structure library entries")
    library_analyze.add_argument("--xrd", required=True, help="Experimental XRD file")
    library_analyze.add_argument("--library-db", required=True, help="Ophi structure library SQLite path")
    library_analyze.add_argument("--elements", nargs="+", required=True, help="Main elements, such as Zr Fe Ge")
    library_analyze.add_argument("--extra-elements", nargs="*", default=[], help="Optional impurity/contamination elements")
    library_analyze.add_argument("--radiation", default="CuKa")
    library_analyze.add_argument("--two-theta-min", type=float, default=10.0)
    library_analyze.add_argument("--two-theta-max", type=float, default=90.0)
    library_analyze.add_argument("--tolerance", type=float, default=0.20)
    library_analyze.add_argument("--max-candidates", type=int, default=20)
    library_analyze.add_argument("--max-phases", type=int, default=3)
    library_analyze.add_argument("--out-dir", required=True)
    library_analyze.add_argument("--project", default="OphiProject")
    library_analyze.add_argument("--sample-id", default="Sample")
    library_analyze.add_argument("--target-candidate-id", help="Explicit target phase structure id; all other phases are treated as impurities")
    library_analyze.add_argument("--target-formula", help="Explicit target phase formula if a structure id is unavailable")
    library_analyze.add_argument("--export-folder", help="Optional deterministic export folder for PNG/CSV/JSON bundle")
    library_analyze.add_argument("--export-date", help="Optional YYYYMMDD export date override")

    mp_test = sub.add_parser("mp-test", help="Test Materials Project API key and mp-api availability")
    mp_test.add_argument("--api-key", help="Materials Project API key; otherwise MP_API_KEY/.env is used")
    mp_test.add_argument("--env-path", help="Optional .env path")

    mp_harvest = sub.add_parser("mp-harvest", help="Harvest structures from Materials Project into the local library")
    mp_harvest.add_argument("--library-db", required=True, help="Ophi structure library SQLite path")
    mp_harvest.add_argument("--elements", nargs="+", required=True, help="Main target elements")
    mp_harvest.add_argument("--impurities", nargs="*", default=[], help="Optional impurity/substitution elements")
    mp_harvest.add_argument("--mode", choices=["conservative", "normal", "broad"], default="normal")
    mp_harvest.add_argument("--max-systems", type=int, default=80)
    mp_harvest.add_argument("--max-entries-per-system", type=int, default=50)
    mp_harvest.add_argument("--api-key", help="Materials Project API key; otherwise MP_API_KEY/.env is used")
    mp_harvest.add_argument("--env-path", help="Optional .env path")
    mp_harvest.add_argument("--dry-run", action="store_true", help="Only show generated chemical systems")

    xrd_validate = sub.add_parser("xrd-validate", help="Validate Ophi CIF XRD simulation against a VESTA/reference peak table")
    xrd_validate.add_argument("--cif", required=True, help="CIF file to simulate")
    xrd_validate.add_argument("--reference", required=True, help="VESTA/reference peak table CSV/TXT")
    xrd_validate.add_argument("--phase", help="Optional phase name filter for all_peak_positions.csv")
    xrd_validate.add_argument("--out-dir", required=True, help="Validation report output folder")
    xrd_validate.add_argument("--label", help="Report label")
    xrd_validate.add_argument("--two-theta-min", type=float, default=10.0)
    xrd_validate.add_argument("--two-theta-max", type=float, default=90.0)
    xrd_validate.add_argument("--wavelength", type=float, default=1.54056, help="Validation wavelength in Angstrom; default Cu Kalpha1/Kalpha2")
    xrd_validate.add_argument("--debye-waller-b", type=float, default=0.0, help="Uniform Debye-Waller B factor for debug validation")

    xrd_debug_simulate = sub.add_parser("xrd-debug-simulate", help="Export a transparent single-CIF simulated XRD peak table")
    xrd_debug_simulate.add_argument("--cif", required=True)
    xrd_debug_simulate.add_argument("--out", required=True)
    xrd_debug_simulate.add_argument("--radiation", default="CuKalpha12")
    xrd_debug_simulate.add_argument("--wavelength", type=float, default=1.54056)
    xrd_debug_simulate.add_argument("--two-theta-min", type=float, default=5.0)
    xrd_debug_simulate.add_argument("--two-theta-max", type=float, default=90.0)
    xrd_debug_simulate.add_argument("--debye-waller-b", type=float, default=0.0)

    xrd_debug_compare = sub.add_parser("xrd-debug-compare", help="Compare an Ophi debug peak table against a VESTA/reference table")
    xrd_debug_compare.add_argument("--ophi", required=True)
    xrd_debug_compare.add_argument("--reference", required=True)
    xrd_debug_compare.add_argument("--out", required=True)
    xrd_debug_compare.add_argument("--position-tolerance", type=float, default=0.05)
    xrd_debug_compare.add_argument("--phase")

    app = sub.add_parser("app", help="Open the local desktop window")

    args = parser.parse_args(argv)
    if args.command == "analyze":
        result = run_analysis(
            xrd_file=args.xrd,
            candidate_dirs=args.cif_dir,
            elements=args.elements,
            extra_elements=args.extra_elements,
            out_dir=args.out_dir,
            radiation=args.radiation,
            two_theta_min=args.two_theta_min,
            two_theta_max=args.two_theta_max,
            tolerance_deg=args.tolerance,
            max_candidates=args.max_candidates,
            max_phases=args.max_phases,
            cache_path=args.cache_path,
            use_cache=not args.no_cache,
        )
        print("Ophi XRD analysis finished.")
        for name, path in result.outputs.items():
            print(f"{name}: {path}")
        for warning in result.warnings:
            print(f"warning: {warning}", file=sys.stderr)
        return 0
    if args.command == "cache-build":
        cache = CandidateCache(args.cache_path)
        summary = build_candidate_cache(
            cache,
            args.cif_dir,
            set(args.elements),
            radiation=args.radiation,
            two_theta_range=(args.two_theta_min, args.two_theta_max),
        )
        payload = {"cache_path": args.cache_path, **summary, "stats": cache.stats()}
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 0
    if args.command == "library-import":
        library = StructureLibrary(args.library_db)
        summary = LocalFolderProvider().import_folder(args.folder, library)
        print(json.dumps({"library_db": args.library_db, **summary}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "library-list":
        library = StructureLibrary(args.library_db)
        rows = [
            {
                "internal_id": entry.internal_id,
                "formula": entry.formula,
                "chemical_system": entry.chemical_system,
                "source": entry.source,
                "enabled_for_matching": entry.enabled_for_matching,
            }
            for entry in library.list_structures()
        ]
        print(json.dumps({"library_db": args.library_db, "structures": rows}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "library-cache-xrd":
        library = StructureLibrary(args.library_db)
        summary = build_library_xrd_cache(
            library,
            radiation=args.radiation,
            two_theta_range=(args.two_theta_min, args.two_theta_max),
        )
        print(json.dumps({"library_db": args.library_db, **summary}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "library-analyze":
        result = run_library_analysis(
            xrd_file=args.xrd,
            library_db=args.library_db,
            elements=args.elements,
            extra_elements=args.extra_elements,
            out_dir=args.out_dir,
            project=args.project,
            sample_id=args.sample_id,
            radiation=args.radiation,
            two_theta_min=args.two_theta_min,
            two_theta_max=args.two_theta_max,
            tolerance_deg=args.tolerance,
            max_candidates=args.max_candidates,
            max_phases=args.max_phases,
            target_candidate_id=args.target_candidate_id,
            target_formula=args.target_formula,
        )
        payload = {"outputs": result.outputs, "warnings": result.warnings}
        if args.export_folder:
            payload["exports"] = export_analysis_bundle(
                result,
                args.export_folder,
                project=args.project,
                sample_id=args.sample_id,
                date_text=args.export_date,
            )
        print("Ophi library analysis finished.")
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 0
    if args.command == "mp-test":
        provider = MaterialsProjectProvider(api_key=args.api_key, env_path=args.env_path)
        try:
            payload = provider.test_connection()
        except ProviderError as exc:
            print(str(exc), file=sys.stderr)
            return 2
        print(json.dumps(payload, indent=2, ensure_ascii=False))
        return 0
    if args.command == "mp-harvest":
        systems = generate_chemical_systems(args.elements, args.impurities, mode=args.mode, max_systems=args.max_systems)
        if args.dry_run:
            print(json.dumps({"chemical_systems": systems, "count": len(systems)}, indent=2, ensure_ascii=False))
            return 0
        provider = MaterialsProjectProvider(api_key=args.api_key, env_path=args.env_path)
        library = StructureLibrary(args.library_db)
        try:
            summary = provider.harvest_to_library(library, systems, max_entries_per_system=args.max_entries_per_system)
            cache_summary = build_library_xrd_cache(library)
        except ProviderError as exc:
            print(str(exc), file=sys.stderr)
            return 2
        print(json.dumps({"library_db": args.library_db, **summary, "xrd_cache": cache_summary}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "xrd-validate":
        config = XRDConfig.validation_default()
        config = XRDConfig(
            radiation_source="CuKalpha12" if abs(args.wavelength - 1.54056) < 1e-6 else f"lambda_{args.wavelength:.5f}",
            wavelength_angstrom=args.wavelength,
            two_theta_min=args.two_theta_min,
            two_theta_max=args.two_theta_max,
            debye_waller_b=args.debye_waller_b,
        )
        outputs = ValidationRunner(simulator=simulate_cif_with_config).run(
            cif_file=args.cif,
            reference_file=args.reference,
            out_dir=args.out_dir,
            config=config,
            label=args.label or Path(args.cif).stem,
            phase=args.phase,
        )
        print("Ophi XRD validation finished.")
        print(json.dumps(outputs, indent=2, ensure_ascii=False))
        return 0
    if args.command == "xrd-debug-simulate":
        config = XRDConfig(
            radiation_source=args.radiation,
            wavelength_angstrom=args.wavelength,
            two_theta_min=args.two_theta_min,
            two_theta_max=args.two_theta_max,
            debye_waller_b=args.debye_waller_b,
        )
        out = simulate_cif_debug_table(args.cif, args.out, config)
        print(json.dumps({"debug_peaks": out, "config": config.to_dict()}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "xrd-debug-compare":
        out = compare_debug_peak_tables(
            args.ophi,
            args.reference,
            args.out,
            position_tolerance=args.position_tolerance,
            phase=args.phase,
        )
        print(json.dumps({"comparison_csv": out}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "app":
        launch_app()
        return 0
    app.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
