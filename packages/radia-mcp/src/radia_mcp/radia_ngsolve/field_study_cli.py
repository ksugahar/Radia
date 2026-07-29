"""Headless bridge from the Radia Simulink Field Study block to NGSolve."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

from .vol2d_dynamics import analyze_vol2d_dynamics
from .vol2d_scalar import analyze_vol2d_scalar


SCALAR_PHYSICS = {"electrostatic", "current_flow", "steady_heat"}


def _content(exports: Mapping[str, Any], name: str) -> str:
    row = exports.get(name)
    if not isinstance(row, Mapping) or not isinstance(row.get("content"), str):
        raise ValueError(f"field study did not produce {name}")
    return str(row["content"])


def _exports(result: Mapping[str, Any], physics: str) -> Mapping[str, Any]:
    if physics in SCALAR_PHYSICS:
        exports = result.get("exports")
    else:
        nested = result.get("result")
        exports = nested.get("exports") if isinstance(nested, Mapping) else None
    if not isinstance(exports, Mapping):
        raise ValueError("field study did not produce an export bundle")
    return exports


def run(request: Mapping[str, Any], *, msh_output: Path) -> dict[str, Any]:
    physics = str(request.get("physics", ""))
    prepared = dict(request)
    prepared["export_basename"] = msh_output.stem
    if physics in SCALAR_PHYSICS:
        prepared["operation"] = "solve"
        result = analyze_vol2d_scalar(prepared)
    elif physics == "harmonic_eddy":
        prepared["operation"] = "harmonic"
        result = analyze_vol2d_dynamics(prepared)
    else:
        raise ValueError(
            "physics must be electrostatic, current_flow, steady_heat, or harmonic_eddy"
        )

    exports = _exports(result, physics)
    msh_output.parent.mkdir(parents=True, exist_ok=True)
    msh_output.write_text(_content(exports, "gmsh_msh"), encoding="utf-8")
    geo = msh_output.with_suffix(".geo")
    geo.write_text(_content(exports, "gmsh_geo"), encoding="utf-8")
    Path(str(geo) + ".opt").write_text(
        _content(exports, "gmsh_geo_opt"), encoding="utf-8"
    )
    Path(str(msh_output) + ".opt").write_text(
        _content(exports, "gmsh_msh_opt"), encoding="utf-8"
    )
    return dict(result)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--msh-output", required=True)
    args = parser.parse_args(argv)

    request_path = Path(args.request).resolve()
    output_path = Path(args.output).resolve()
    msh_output = Path(args.msh_output).resolve()
    request = json.loads(request_path.read_text(encoding="utf-8-sig"))
    if not isinstance(request, dict):
        raise ValueError("field study request must be a JSON object")
    result = run(request, msh_output=msh_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
