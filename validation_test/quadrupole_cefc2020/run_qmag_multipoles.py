"""Allowed multipoles of the CEFC 2020 quadrupole from an HDiv-MMM solve.

Solves the linear ``mu_r`` case on ``--mesh`` (HEX or TET, any curve order) with the HDiv order
``--hdiv-order``, evaluates the field on ``--n-phi`` points of the circle ``r = --radius`` in the
midplane and expands ``B_y + i B_x = sum_n C_n (z / r)^(n - 1)`` by FFT.  Reported: the main
quadrupole ``B_2``, the allowed harmonics ``b_6, b_10, b_14`` (and skew ``a_n``) in units of
``1e-4 B_2``, and the quadrupole-forbidden ``b_3, b_4, b_5`` as a mesh-symmetry check.  This is
the observable that separates curved from straight pole faces and BDM1 from BDM2 (README,
"Multipole convergence"); ``B_perp(15 mm)`` alone does not.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import ngsolve as ng  # noqa: E402
import numpy as np  # noqa: E402
import radia as rad  # noqa: E402

import qmag_case as Q  # noqa: E402
from radia import vim  # noqa: E402
from radia.vim import mesh_conformity_report  # noqa: E402

SCHEMA = "radia.qmag-cefc2020-multipoles.v1"


def circle_points(radius: float, n_phi: int) -> np.ndarray:
    phi = 2.0 * np.pi * np.arange(int(n_phi)) / int(n_phi)
    return np.column_stack([radius * np.cos(phi), radius * np.sin(phi), np.zeros(len(phi))])


def multipoles(field: np.ndarray, orders=(6, 10, 14), forbidden=(3, 4, 5)) -> dict:
    """FFT of ``B_y + i B_x`` on equispaced circle points; ``C[k]`` multiplies ``(z/r)^k`` (n = k + 1)."""
    field = np.asarray(field, dtype=float)
    n_phi = len(field)
    coefficients = np.fft.fft(field[:, 1] + 1j * field[:, 0]) / n_phi
    main = coefficients[1]
    if abs(main) == 0.0:
        raise RuntimeError("zero quadrupole component")
    result = {"B2_T": float(abs(main)), "n_phi": int(n_phi)}
    for n in orders:
        ratio = coefficients[n - 1] / main
        result[f"b{n}_units"] = float(ratio.real * 1.0e4)
        result[f"a{n}_units"] = float(ratio.imag * 1.0e4)
    for n in forbidden:
        result[f"forbidden_b{n}_units"] = float(abs(coefficients[n - 1] / main) * 1.0e4)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mu-r", type=float, default=1000.0)
    parser.add_argument("--hdiv-order", type=int, choices=(1, 2), default=1)
    parser.add_argument("--gram-eps", type=float, default=1.0e-10)
    parser.add_argument("--radius", type=float, default=0.015, help="circle radius in metres (pole tip at 0.020)")
    parser.add_argument("--n-phi", type=int, default=128)
    parser.add_argument("--current-density", type=float, default=3.0, help="A/mm^2")
    options = parser.parse_args(argv)
    mesh_path = options.mesh.resolve()
    mesh = ng.Mesh(str(mesh_path))
    conformity = mesh_conformity_report(mesh)
    if not conformity["conforming"]:
        raise RuntimeError(f"non-conforming iron mesh: {conformity}")
    points = circle_points(options.radius, options.n_phi)
    rad.UtiDelAll()
    coil, coil_manifest = Q.build_qmag_coils(options.current_density)
    source_h = rad.RadiaField(coil, "h")
    started = time.perf_counter()
    with ng.TaskManager():
        result = vim.Solve(mesh, H_ext=source_h, order=options.hdiv_order, gram_eps=options.gram_eps, tol=1e-8,
                           maxit=12000, preconditioner="auto", mu_r=float(options.mu_r))
        demag_h = np.asarray(vim.FieldFromSolution(result, points, algorithm="direct"), dtype=float)
    field = Q.MU0 * (demag_h + np.asarray(rad.Fld(coil, "h", points), dtype=float))
    wall = time.perf_counter() - started
    families = {}
    for el in mesh.Elements(ng.VOL):
        key = str(el.type).rsplit(".", 1)[-1]
        families[key] = families.get(key, 0) + 1
    report = {
        "schema": SCHEMA, "generated_at_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
        "radia_version": getattr(rad, "__version__", None), "radia_file": rad.__file__,
        "mesh": {"path": str(mesh_path), "sha256": hashlib.sha256(mesh_path.read_bytes()).hexdigest(),
                 "ne": int(mesh.ne), "nv": int(mesh.nv), "families": families, "conformity": conformity},
        "hdiv_order": int(options.hdiv_order), "mu_r": float(options.mu_r), "gram_eps": float(options.gram_eps),
        "ndof": int(result["ndof"]), "linear_iterations": result.get("linear_iterations"),
        "coils": coil_manifest, "radius_m": float(options.radius),
        "multipoles": multipoles(field), "field_T": field.tolist(), "wall_s": wall,
    }
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(report, indent=1), encoding="utf-8")
    m = report["multipoles"]
    print(f"{mesh_path.name} order {options.hdiv_order}: ndof {report['ndof']}, B2 {m['B2_T']:.5f} T, "
          f"b6 {m['b6_units']:+.3f} b10 {m['b10_units']:+.3f} b14 {m['b14_units']:+.3f} units, "
          f"forbidden b3 {m['forbidden_b3_units']:.3f}", flush=True)
    print("wrote", options.output, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
