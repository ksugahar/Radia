"""Physical spherical-capacitor benchmark for Cubit-exported curved TET meshes.

The concentric electrodes have radii a=0.3 and b=1, potentials 1 and 0,
and uniform relative permittivity 1. Exact capacitance is 4*pi*a*b/(b-a).
Cubit and Netgen/OCC use independent meshing paths with nearby H1 DOF counts.
All Cubit execution is batch/nographics; no GUI is started.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import platform
import subprocess
import time
from pathlib import Path

from cubit_mesh_export.check import check_consistency
from netgen.occ import Pnt, Sphere
from ngsolve import (
    BND, VOL, BilinearForm, CoefficientFunction, GridFunction, H1, Integrate,
    LinearForm, Mesh, TaskManager, sqrt, x, y, z, grad, dx,
)

A, B = 0.3, 1.0
C_EXACT = 4 * math.pi * A * B / (B-A)
ORDERS = (1, 2, 3)


def solve(mesh: Mesh, order: int) -> dict:
    names = set(mesh.GetBoundaries())
    materials = set(mesh.GetMaterials())
    if names != {"inner", "outer"} or materials != {"dielectric"}:
        raise AssertionError(f"unexpected labels: {names}, {materials}")
    fes = H1(mesh, order=order, dirichlet="inner|outer")
    gfu = GridFunction(fes)
    gfu.Set(mesh.BoundaryCF({"inner": 1, "outer": 0}), BND)
    u, v = fes.TnT()
    stiffness = BilinearForm(grad(u)*grad(v)*dx).Assemble()
    rhs = LinearForm(fes).Assemble()
    residual = rhs.vec - stiffness.mat*gfu.vec
    gfu.vec.data += stiffness.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")*residual
    radial = sqrt(x*x+y*y+z*z)
    factor = A*B/(B-A)
    phi_exact = factor*(1/radial - 1/B)
    grad_exact = -factor*CoefficientFunction((x, y, z))/(radial**3)
    grad_error = grad(gfu)-grad_exact
    energy = float(Integrate(grad(gfu)*grad(gfu), mesh, VOL))
    areas = {name: float(Integrate(CoefficientFunction(1), mesh, BND,
                                  definedon=mesh.Boundaries(name)))
             for name in ("inner", "outer")}
    kinds: dict[str, int] = {}
    for el in mesh.Elements(VOL):
        kind = str(el.type)
        kinds[kind] = kinds.get(kind, 0)+1
    return {
        "order": order, "elements": kinds, "ndof": fes.ndof,
        "boundary_areas": areas,
        "boundary_label_preservation": sum(area > 0 for area in areas.values())/2,
        "material_label_preservation": 1.0,
        "capacitance": energy, "exact_capacitance": C_EXACT,
        "capacitance_error_pct": 100*(energy/C_EXACT-1),
        "potential_l2_error": float(sqrt(Integrate((gfu-phi_exact)**2, mesh))),
        "electric_field_l2_error": float(sqrt(Integrate(grad_error*grad_error, mesh))),
    }


def journal(run_dir: Path) -> str:
    exports = "\n".join(
        f'export netgen "{(run_dir / f"capacitor_o{p}.vol").as_posix()}" order {p} overwrite'
        for p in ORDERS
    )
    return (
        "reset\ncreate sphere radius 1\ncreate sphere radius 0.3\n"
        "subtract volume 2 from volume 1\n"
        "volume all scheme tetmesh\nvolume all size 0.22\nmesh volume all\n"
        "block 1 volume all\nblock 1 name \"dielectric\"\n"
        "sideset 1 surface 1\nsideset 1 name \"outer\"\n"
        "sideset 2 surface 3\nsideset 2 name \"inner\"\n"
        f"{exports}\nexit 0\n"
    )


def cubit_run(exe: Path, run_dir: Path, timeout: int) -> list[dict]:
    run_dir.mkdir(parents=True, exist_ok=True)
    driver = run_dir / "driver.jou"
    driver.write_text(journal(run_dir), encoding="utf-8")
    started = time.monotonic()
    proc = subprocess.run(
        [str(exe), "-batch", "-nographics", "-nojournal",
         "-commandplugindir", str(exe.parent / "plugins"), str(driver)],
        cwd=run_dir, capture_output=True, text=True, encoding="utf-8",
        errors="replace", timeout=timeout,
    )
    elapsed = time.monotonic()-started
    (run_dir / "cubit.log").write_text(
        f"exit={proc.returncode}\n=== STDOUT ===\n{proc.stdout}\n=== STDERR ===\n{proc.stderr}",
        encoding="utf-8",
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Cubit exited {proc.returncode}; see {run_dir / 'cubit.log'}")
    rows = []
    for order in ORDERS:
        vol = run_dir / f"capacitor_o{order}.vol"
        if not vol.is_file():
            raise RuntimeError(f"Cubit exit {proc.returncode}; missing {vol}")
        gate = check_consistency(vol, strict_labels=True,
                                 required_boundaries=("inner", "outer"),
                                 required_materials=("dielectric",))
        row = solve(Mesh(str(vol)), order)
        row.update(route="cubit_export_netgen", run_dir=str(run_dir),
                   cubit_exit=proc.returncode, mesh_seconds=elapsed,
                   vol_check_passed=gate["passed"],
                   vol_check_warnings=gate["warnings"],
                   invalid_jacobian_samples=gate["quality"]["invalid_jacobian_sample_count"])
        rows.append(row)
    return rows


def occ_run(run_dir: Path) -> list[dict]:
    run_dir.mkdir(parents=True, exist_ok=True)
    outer = Sphere(Pnt(0, 0, 0), B)
    inner = Sphere(Pnt(0, 0, 0), A)
    outer.faces.name = "outer"
    inner.faces.name = "inner"
    shell = outer-inner
    shell.mat("dielectric")
    started = time.monotonic()
    mesh = shell.GenerateMesh(maxh=0.16)
    elapsed = time.monotonic()-started
    rows = []
    for order in ORDERS:
        mesh.Curve(order)
        row = solve(mesh, order)
        row.update(route="netgen_occ_native", run_dir=str(run_dir),
                   mesh_seconds=elapsed)
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cubit-exe", type=Path,
                        default=Path(r"C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.exe"))
    parser.add_argument("--output", type=Path,
                        default=Path(r"C:\temp\cubit-paper-capacitor-benchmark"))
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()
    if args.repeats < 2:
        parser.error("at least two independent runs are required")
    if not args.cubit_exe.is_file():
        parser.error(f"Cubit executable not found: {args.cubit_exe}")
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    with TaskManager():
        for repeat in range(args.repeats):
            print(f"Cubit capacitor repeat {repeat+1}/{args.repeats}", flush=True)
            rows.extend({**row, "repeat": repeat+1} for row in cubit_run(
                args.cubit_exe, args.output / f"cubit_{repeat+1}", args.timeout))
            print(f"Netgen/OCC capacitor repeat {repeat+1}/{args.repeats}", flush=True)
            rows.extend({**row, "repeat": repeat+1} for row in occ_run(
                args.output / f"occ_{repeat+1}"))
    result = {
        "protocol": "concentric spherical capacitor a=0.3 b=1 Vinner=1 Vouter=0 eps=1; Cubit TET size=0.22; Netgen OCC maxh=0.16; fixed topology within each p-series",
        "environment": {"python": platform.python_version(),
                        "platform": platform.platform(),
                        "ngsolve": importlib.metadata.version("ngsolve"),
                        "cubit_mesh_export": importlib.metadata.version("cubit-mesh-export")},
        "cubit_exe": str(args.cubit_exe), "rows": rows,
    }
    for route in ("cubit_export_netgen", "netgen_occ_native"):
        for repeat in range(1, args.repeats+1):
            series = [r for r in rows if (r["route"], r["repeat"]) == (route, repeat)]
            series.sort(key=lambda r: r["order"])
            if len(series) != len(ORDERS):
                raise AssertionError(f"incomplete series: {route} repeat {repeat}")
            if not all(series[i+1]["electric_field_l2_error"] < series[i]["electric_field_l2_error"]
                       and abs(series[i+1]["capacitance_error_pct"])
                       < abs(series[i]["capacitance_error_pct"])
                       for i in range(len(series)-1)):
                raise AssertionError(f"capacitor convergence failed: {route}")
            if route == "cubit_export_netgen" and not all(
                    r["vol_check_passed"] and r["invalid_jacobian_samples"] == 0
                    for r in series[1:]):
                raise AssertionError("high-order capacitor quality gate failed")
            if not all(r["boundary_label_preservation"] == 1 for r in series):
                raise AssertionError("capacitor electrode label missing")
    max_repeat_delta = 0.0
    for row in rows:
        if row["repeat"] == 1:
            continue
        first = next(r for r in rows if (r["route"], r["order"], r["repeat"])
                     == (row["route"], row["order"], 1))
        if row["elements"] != first["elements"] or row["ndof"] != first["ndof"]:
            raise AssertionError("repeat topology/DOF drift")
        max_repeat_delta = max(max_repeat_delta, *(
            abs(row[k]-first[k]) for k in
            ("capacitance", "potential_l2_error", "electric_field_l2_error")))
    result["max_repeat_metric_delta"] = max_repeat_delta
    if max_repeat_delta > 1e-10:
        raise AssertionError(f"repeat metric drift: {max_repeat_delta}")
    path = args.output / "results.json"
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
