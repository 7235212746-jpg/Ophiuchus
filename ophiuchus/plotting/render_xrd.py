from __future__ import annotations

import argparse
import json

from .xrd_plotter import render_xrd_from_analysis_json


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Render an Ophiuchus XRD plot from analysis JSON")
    parser.add_argument("--analysis-json", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--svg")
    args = parser.parse_args(argv)
    outputs = render_xrd_from_analysis_json(args.analysis_json, args.out, out_svg=args.svg)
    print(json.dumps(outputs, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
