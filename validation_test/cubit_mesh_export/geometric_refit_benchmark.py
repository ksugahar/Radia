"""Netgen L2 coefficients versus the exporter's geometric refit (no Cubit GUI).

The same Cubit TET and HEX unit-sphere meshes are exported at orders 1-5
twice: once with Netgen's coefficients (CUBIT_MESH_EXPORT_GEOMETRIC_REFIT=0)
and once with the geometric refit.  Geometry error is measured on the
boundary against the exact sphere; the field error uses the electrostatic
manufactured problem of paper_sphere_benchmark.py.  check-vol runs on every
file.  Raw runs stay outside the checkout; the result JSON is written next to
this script.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import os
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
from ngsolve import (
    BND, H1, VOL, BilinearForm, CoefficientFunction, GridFunction, IntegrationRule,
    Integrate, LinearForm, Mesh, TaskManager, cos, dx, grad, sin, sqrt, x, y, z,
)
from cubit_mesh_export.check import check_consistency

ORDERS = (1, 2, 3, 4, 5)
FIELD_ORDERS = (1, 2, 3, 4)
MESHES = {"tet": ("tetmesh", "0.65"), "hex": ("sphere", "0.4")}


def journal(run_dir: Path) -> str:
    parts = []
    for kind, (scheme, size) in MESHES.items():
        parts.append(
            "reset\ncreate sphere radius 1\n"
            f"volume all scheme {scheme}\nvolume all size {size}\n"
            "mesh volume all\nblock 1 volume all\nblock 1 name \"body\"\n"
            "sideset 1 surface all\nsideset 1 name \"outer\"\n")
        for p in ORDERS:
            vol = (run_dir / f"{kind}_o{p}.vol").as_posix()
            parts.append(f'export netgen "{vol}" order {p} overwrite\n')
    return "".join(parts) + "exit 0\n"


def run_cubit(exe: Path, plugin_dir: Path, run_dir: Path, refit: bool, timeout: int) -> float:
    run_dir.mkdir(parents=True, exist_ok=True)
    jou = run_dir / "driver.jou"
    jou.write_text(journal(run_dir), encoding="utf-8")
    env = dict(os.environ, CUBIT_MESH_EXPORT_GEOMETRIC_REFIT="1" if refit else "0")
    started = time.monotonic()
    proc = subprocess.run(
        [str(exe), "-batch", "-nographics", "-nojournal", "-noinitfile",
         "-commandplugindir", str(plugin_dir), str(jou)],
        cwd=run_dir, capture_output=True, text=True, encoding="utf-8",
        errors="replace", timeout=timeout, env=env)
    elapsed = time.monotonic() - started
    (run_dir / "cubit.log").write_text(
        f"exit={proc.returncode}\n=== STDOUT ===\n{proc.stdout}\n=== STDERR ===\n{proc.stderr}",
        encoding="utf-8")
    marker = "geometric refit disabled" if not refit else "geometric refit - "
    if marker not in proc.stdout:
        raise RuntimeError(f"refit={refit}: plugin did not report '{marker}'; see {run_dir / 'cubit.log'}")
    for kind in MESHES:
        for p in ORDERS:
            if not (run_dir / f"{kind}_o{p}.vol").is_file():
                raise RuntimeError(f"missing {kind}_o{p}.vol; see {run_dir / 'cubit.log'}")
    return elapsed


def solve(mesh: Mesh, order: int) -> dict:
    """Manufactured -Delta(phi)=3*pi^2*phi of paper_sphere_benchmark.py.

    The error norms are integrated at order 2p+8: NGSolve's default
    Integrate order (5) under-integrates them from p=3 on (the hex p=3
    field error reads 0.23 instead of 0.36).
    """
    phi = sin(math.pi * x) * cos(math.pi * y) * cos(math.pi * z)
    grad_phi = CoefficientFunction((
        math.pi * cos(math.pi * x) * cos(math.pi * y) * cos(math.pi * z),
        -math.pi * sin(math.pi * x) * sin(math.pi * y) * cos(math.pi * z),
        -math.pi * sin(math.pi * x) * cos(math.pi * y) * sin(math.pi * z),
    ))
    fes = H1(mesh, order=order, dirichlet=".*")
    gfu = GridFunction(fes)
    gfu.Set(phi, BND)
    u, v = fes.TnT()
    a = BilinearForm(grad(u) * grad(v) * dx).Assemble()
    rhs = LinearForm((3 * math.pi ** 2 * phi) * v * dx).Assemble()
    residual = rhs.vec - a.mat * gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * residual
    dphi = grad(gfu) - grad_phi
    q = 2 * order + 8
    return {"ndof": fes.ndof,
            "error_integration_order": q,
            "phi_l2_error": float(sqrt(Integrate((gfu - phi) ** 2, mesh, order=q))),
            "electric_field_l2_error": float(sqrt(Integrate(dphi * dphi, mesh, order=q)))}


def geometry(mesh: Mesh, order: int) -> dict:
    r = sqrt(x * x + y * y + z * z)
    q = 2 * order + 6
    vol = float(Integrate(CoefficientFunction(1), mesh, VOL, order=q))
    area = float(Integrate(CoefficientFunction(1), mesh, BND, order=q))
    radial_l2 = math.sqrt(float(Integrate((r - 1) ** 2, mesh, BND, order=q)) / area)
    s = np.linspace(0.05, 0.95, 19)
    pts = [(a, b, 0) for a in s for b in s if a + b < 1.0] + [(a, 0.0, 0) for a in s]
    ir = IntegrationRule(points=pts, weights=[1] * len(pts))
    xyz = np.array(CoefficientFunction((x, y, z))(mesh.MapToAllElements(ir, BND)))
    return {"volume_error": vol / (4 * math.pi / 3) - 1,
            "area_error": area / (4 * math.pi) - 1,
            "radial_l2": radial_l2,
            "radial_max_sampled": float(np.abs(np.linalg.norm(xyz, axis=1) - 1).max())}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cubit-exe", type=Path,
                        default=Path(r"C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.exe"))
    parser.add_argument("--plugin-dir", type=Path, required=True,
                        help="directory holding the cubit_mesh_export.ccm under test")
    parser.add_argument("--output", type=Path, default=Path(r"C:\temp\cubit-geometric-refit-benchmark"))
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()

    rows, timing = [], {}
    with TaskManager():
        for refit in (False, True):
            label = "geometric_refit" if refit else "netgen_l2"
            run_dir = args.output / label
            timing[label] = run_cubit(args.cubit_exe, args.plugin_dir, run_dir, refit, args.timeout)
            for kind in MESHES:
                for p in ORDERS:
                    vol = run_dir / f"{kind}_o{p}.vol"
                    mesh = Mesh(str(vol))
                    row = {"coefficients": label, "kind": kind, "order": p, **geometry(mesh, p)}
                    if p in FIELD_ORDERS:
                        fe = solve(mesh, p)
                        row.update(fe)
                    gate = check_consistency(vol, strict_labels=True,
                                             required_boundaries=("outer",),
                                             required_materials=("body",))
                    row.update(vol_check_passed=gate["passed"],
                               invalid_jacobian_samples=gate["quality"]["invalid_jacobian_sample_count"])
                    rows.append(row)
                    print(f"{label:16s} {kind} p={p}: radial_l2 {row['radial_l2']:.3e}  "
                          f"vol {row['volume_error']:+.3e}", flush=True)

    def pick(label, kind, p):
        return next(r for r in rows if (r["coefficients"], r["kind"], r["order"]) == (label, kind, p))

    # Acceptance: the refit never worsens geometry, keeps valid Jacobians,
    # removes the odd-order stall (p=3 beats p=2 by 10x) and beats Netgen at
    # p=5 by 100x, and p=5 beats p=4 (it did not while the .vol writer
    # rounded coefficients to 8 digits, before 2.1.1).
    for kind in MESHES:
        for p in ORDERS[1:]:
            new, old = pick("geometric_refit", kind, p), pick("netgen_l2", kind, p)
            if new["radial_l2"] > old["radial_l2"]:
                raise AssertionError(f"refit worsened {kind} p={p}")
            if not new["vol_check_passed"] or new["invalid_jacobian_samples"]:
                raise AssertionError(f"refit {kind} p={p} failed check-vol")
        if pick("geometric_refit", kind, 3)["radial_l2"] * 10 >= pick("geometric_refit", kind, 2)["radial_l2"]:
            raise AssertionError(f"refit {kind}: p=3 is not 10x better than p=2")
        if pick("geometric_refit", kind, 5)["radial_l2"] * 100 >= pick("netgen_l2", kind, 5)["radial_l2"]:
            raise AssertionError(f"refit {kind}: p=5 is not 100x better than Netgen p=5")
        if pick("geometric_refit", kind, 5)["radial_l2"] >= pick("geometric_refit", kind, 4)["radial_l2"]:
            raise AssertionError(f"refit {kind}: p=5 does not improve on p=4")

    result = {
        "protocol": "unit sphere; Cubit tet size 0.65 / hex sphere scheme size 0.4; "
                    "export netgen orders 1-5 with Netgen L2 coefficients and with the geometric refit; "
                    "field problem of paper_sphere_benchmark.py for orders 1-4",
        "environment": {"platform_class": platform.system(), "python": platform.python_version(),
                        "platform": platform.platform(),
                        "ngsolve": importlib.metadata.version("ngsolve"),
                        "plugin_dir": str(args.plugin_dir), "cubit_exe": str(args.cubit_exe)},
        "cubit_seconds": timing,
        "rows": rows,
    }
    target = Path(__file__).with_name("geometric_refit_benchmark_results.json")
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {target}")


if __name__ == "__main__":
    main()
