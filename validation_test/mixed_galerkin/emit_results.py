"""Write the one numerical artifact for the mixed Galerkin study.

`validation_test/` owns the numbers. Documentation notebooks and talk material
read this file; they do not recompute it. That rule exists because the
alternative was tried: the IGTE deck recomputed the cylinder and sphere errors
itself and reported 0.0525 % where this directory reports 0.0639 %, on the same
frequency grid, with the same Y_mixed_galerkin. The physics had one
implementation and the *measurement* had two, and the second one differenced
magnitudes instead of complex values, which silently discards the phase error.

So each case exposes summary(), which runs its sweep once and returns the
numbers; the script's printed table and this file are the same sweep. Every
case records its metric as a string, because the metric is exactly what was
ambiguous.

Run:
    python emit_results.py                 # analytic cases (seconds)
    python emit_results.py --out other.json

The NGSolve cases (cube3d/06_ngsolve_ground_truth.py, lshape3d_ngsolve_mellin.py)
are not called from here yet: they are FEM solves rather than quadrature, and
belong on a quiet compute host under the lab's benchmark rule. Their numbers
are recorded in `pending_cases` so this file states what it does not cover
rather than quietly omitting it.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent

# (module path, label) -- each module must expose summary() -> dict
CASES = [
    ("cylinder/01_no_d_baseline.py", "cylinder_planar_sibc"),
    ("cylinder/02_senior_tower_truncation.py", "cylinder_senior_tower"),
    ("sphere/01_no_d_baseline.py", "sphere_planar_sibc"),
    ("sphere/02_hoibc_gamma1.py", "sphere_hoibc_gamma1"),
    ("square2d/01_corner_envelope.py", "square2d_corner_envelope"),
]

PENDING = {
    "cube3d_ngsolve_ground_truth": {
        "script": "cube3d/06_ngsolve_ground_truth.py",
        "why": "NGSolve FEM solve; belongs on a quiet compute host",
        "known_result": "0.33 % at rank 20 with the closed K_ss, "
                        "measured against the FEM solution",
    },
    "lshape3d_mellin": {
        "script": "lshape3d_ngsolve_mellin.py",
        "why": "NGSolve FEM solve; belongs on a quiet compute host",
        "known_result": "K_SIBC = S sqrt(sigma/mu) holds to 1.5 % "
                        "on a body with one concave edge",
    },
}


def load(rel: str):
    """Import a numbered script by path; '01_...' is not an identifier."""
    path = HERE / rel
    name = "mg_" + rel.replace("/", "_").replace(".py", "")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def versions() -> dict:
    out = {"python": platform.python_version(), "platform": platform.platform()}
    for name in ("numpy", "scipy"):
        try:
            out[name] = __import__(name).__version__
        except ImportError:
            out[name] = None
    return out


def build() -> dict:
    cases = {}
    for rel, label in CASES:
        mod = load(rel)
        if not hasattr(mod, "summary"):
            raise AttributeError(f"{rel} has no summary(); the emitter needs one")
        r = mod.summary()
        if r.get("case") != label:
            raise ValueError(f"{rel} reports case {r.get('case')!r}, expected {label!r}")
        r["script"] = rel
        cases[label] = r
        print(f"  {label:28s} <- {rel}")
    return {
        "schema": "radia.validation.mixed_galerkin.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "versions": versions(),
        "owner": "validation_test/mixed_galerkin -- documentation and talk "
                 "material read this file and must not recompute it",
        "cases": cases,
        "pending_cases": PENDING,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(HERE / "results" / "mixed_galerkin_results.json"))
    args = ap.parse_args()

    print("running the analytic mixed Galerkin cases ...")
    data = build()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    for label, r in data["cases"].items():
        headline = (r.get("max_error_pct")
                    or r.get("planar_max_error_pct")
                    or (r.get("by_n_dof", {}).get("4", {}).get("max_error_pct")))
        print(f"  {label:28s} {headline:.5f} %" if headline is not None
              else f"  {label:28s} (see file)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
