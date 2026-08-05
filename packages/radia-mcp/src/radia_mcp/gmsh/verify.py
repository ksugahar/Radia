"""One-call verification gate (run_matlab_test_file twin).

Runs every applicable check on a GMSH artifact and returns a
structured pass/fail report -- the "test runner" position in the
matlab-mcp-core-server verb set.  A ``.msh`` gets the structural +
Jacobian + field-finiteness gates plus its sibling ``.geo`` (if any);
a ``.geo`` gets the deep launch check plus every merged ``.msh``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .msh_inspect import validate_geo, validate_msh


def verify_artifact(path: str | Path,
                    check_jacobians: bool = True,
                    check_options: bool = False) -> dict[str, Any]:
    """Run all applicable gates on a .msh or .geo artifact."""
    p = Path(path)
    if not p.is_file():
        return {"ok": False, "artifact": str(p),
                "error": f"file not found: {p}"}

    gates: list[dict[str, Any]] = []

    def _run(name: str, result: dict[str, Any]) -> dict[str, Any]:
        gates.append({
            "gate": name,
            "status": result.get("status", "needs_attention"),
            "failed_checks": [k for k, v in result.get("checks", {}).items()
                              if not v],
            "errors": result.get("errors", [])[:3],
        })
        return result

    suffix = p.suffix.lower()
    if suffix == ".geo":
        geo = _run(f"geo:{p.name}",
                   validate_geo(p, deep=True, check_options=check_options))
        for target in geo.get("merge_targets", []):
            resolved = str(target.get("resolved", ""))
            if target.get("exists") and resolved.lower().endswith(".msh"):
                _run(f"msh:{Path(resolved).name}",
                     validate_msh(Path(resolved),
                                  check_jacobians=check_jacobians))
    elif suffix == ".msh":
        _run(f"msh:{p.name}",
             validate_msh(p, check_jacobians=check_jacobians))
        sibling_geo = p.with_suffix(".geo")
        if sibling_geo.is_file():
            _run(f"geo:{sibling_geo.name}",
                 validate_geo(sibling_geo, deep=True,
                              check_options=check_options))
    else:
        return {"ok": False, "artifact": str(p),
                "error": f"unsupported artifact type {p.suffix!r} "
                         f"(expected .msh or .geo)"}

    failed = [g["gate"] for g in gates if g["status"] != "ok"]
    return {
        "ok": not failed,
        "artifact": str(p),
        "passed": [g["gate"] for g in gates if g["status"] == "ok"],
        "failed": failed,
        "gates": gates,
        "jacobians_checked": bool(check_jacobians),
        "options_checked": bool(check_options),
    }
