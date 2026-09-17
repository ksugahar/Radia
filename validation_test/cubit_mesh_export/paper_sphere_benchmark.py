"""Reproducible Cubit-to-NGSolve curved-mesh benchmark (no Cubit GUI).

Two independent Cubit batch runs produce TET and HEX sphere meshes at orders
1--3. Netgen/OCC TET is the independent conversion/meshing comparator. The
electrostatic manufactured problem is -Delta(phi)=3*pi**2*phi with
phi=sin(pi*x)cos(pi*y)cos(pi*z) and exact Dirichlet data on the whole sphere.
Raw runs, logs and a machine-readable result are retained outside the checkout.
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

from ngsolve import (
    BND, VOL, BilinearForm, CoefficientFunction, GridFunction, H1, Integrate,
    LinearForm, Mesh, TaskManager, cos, dx, grad, sin, sqrt, x, y, z,
)
from netgen.occ import Pnt, Sphere
from cubit_mesh_export.check import check_consistency


EXPECTED_VOLUME = 4 * math.pi / 3
EXPECTED_BOUNDARIES = {"outer"}
EXPECTED_MATERIALS = {"body"}
SAMPLE_BOUNDARIES = {"source", "sink", "sibc", "coil_surface", "outer"}
SAMPLE_MATERIALS = {"coil", "workpiece", "air"}
ORDERS = (1, 2, 3)


def solve(mesh: Mesh, order: int) -> dict:
    phi = sin(math.pi*x) * cos(math.pi*y) * cos(math.pi*z)
    grad_phi = CoefficientFunction((
        math.pi*cos(math.pi*x)*cos(math.pi*y)*cos(math.pi*z),
        -math.pi*sin(math.pi*x)*sin(math.pi*y)*cos(math.pi*z),
        -math.pi*sin(math.pi*x)*cos(math.pi*y)*sin(math.pi*z),
    ))
    fes = H1(mesh, order=order, dirichlet=".*")
    gfu = GridFunction(fes)
    gfu.Set(phi, BND)
    u, v = fes.TnT()
    a = BilinearForm(grad(u)*grad(v)*dx).Assemble()
    rhs = LinearForm((3*math.pi**2*phi)*v*dx).Assemble()
    residual = rhs.vec - a.mat*gfu.vec
    gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * residual
    dphi = grad(gfu) - grad_phi
    volume = float(Integrate(CoefficientFunction(1), mesh, VOL))
    boundaries = sorted(set(mesh.GetBoundaries()))
    materials = sorted(set(mesh.GetMaterials()))
    outer_area = (float(Integrate(CoefficientFunction(1), mesh, BND,
                         definedon=mesh.Boundaries("outer")))
                  if "outer" in boundaries else 0.0)
    kinds: dict[str, int] = {}
    for el in mesh.Elements(VOL):
        kind = str(el.type)
        kinds[kind] = kinds.get(kind, 0) + 1
    return {
        "order": order,
        "elements": kinds,
        "ndof": fes.ndof,
        "volume": volume,
        "volume_error_pct": 100*(volume/EXPECTED_VOLUME - 1),
        "boundaries": boundaries,
        "materials": materials,
        "boundary_label_preservation": float(outer_area > 0),
        "material_label_preservation": len(EXPECTED_MATERIALS & set(materials))/len(EXPECTED_MATERIALS),
        "phi_l2_error": float(sqrt(Integrate((gfu-phi)**2, mesh))),
        "electric_field_l2_error": float(sqrt(Integrate(dphi*dphi, mesh))),
    }


def cubit_journal(kind: str, run_dir: Path) -> str:
    scheme = "tetmesh" if kind == "tet" else "sphere"
    size = "0.65" if kind == "tet" else "0.4"
    exports = "\n".join(
        f'export netgen "{(run_dir / f"{kind}_o{order}.vol").as_posix()}" order {order} overwrite'
        for order in ORDERS
    )
    return (
        "reset\ncreate sphere radius 1\n"
        f"volume all scheme {scheme}\nvolume all size {size}\n"
        "mesh volume all\nblock 1 volume all\nblock 1 name \"body\"\n"
        "sideset 1 surface all\nsideset 1 name \"outer\"\n"
        f"{exports}\nexit 0\n"
    )


def run_cubit(exe: Path, kind: str, run_dir: Path, timeout: int) -> list[dict]:
    run_dir.mkdir(parents=True, exist_ok=True)
    journal = run_dir / "driver.jou"
    journal.write_text(cubit_journal(kind, run_dir), encoding="utf-8")
    started = time.monotonic()
    proc = subprocess.run(
        [str(exe), "-batch", "-nographics", "-nojournal", str(journal)],
        cwd=run_dir, capture_output=True, text=True, encoding="utf-8",
        errors="replace", timeout=timeout,
    )
    elapsed = time.monotonic()-started
    (run_dir / "cubit.log").write_text(
        f"exit={proc.returncode}\n=== STDOUT ===\n{proc.stdout}\n=== STDERR ===\n{proc.stderr}",
        encoding="utf-8",
    )
    rows = []
    for order in ORDERS:
        vol = run_dir / f"{kind}_o{order}.vol"
        if not vol.is_file():
            raise RuntimeError(f"Cubit exit {proc.returncode}; missing {vol}; see {run_dir / 'cubit.log'}")
        row = solve(Mesh(str(vol)), order)
        gate = check_consistency(vol, strict_labels=True,
                                 required_boundaries=("outer",),
                                 required_materials=("body",))
        row.update(route="cubit_export_netgen", kind=kind, run_dir=str(run_dir),
                   cubit_exit=proc.returncode, mesh_seconds=elapsed,
                   vol_check_passed=gate["passed"],
                   invalid_jacobian_samples=gate["quality"]["invalid_jacobian_sample_count"],
                   vol_check_warnings=gate["warnings"])
        rows.append(row)
    return rows


def run_occ(run_dir: Path) -> list[dict]:
    run_dir.mkdir(parents=True, exist_ok=True)
    shape = Sphere(Pnt(0, 0, 0), 1)
    shape.faces.name = "outer"
    shape.mat("body")
    rows = []
    for order in ORDERS:
        started = time.monotonic()
        mesh = shape.GenerateMesh(maxh=0.65)
        mesh.Curve(order)
        mesh_seconds = time.monotonic()-started
        row = solve(mesh, order)
        row.update(route="netgen_occ_native", kind="tet", run_dir=str(run_dir),
                   mesh_seconds=mesh_seconds)
        rows.append(row)
    return rows


def run_label_case(exe: Path, run_dir: Path, timeout: int) -> dict:
    """Check nontrivial sideset/block labels on the bundled EM sample."""
    run_dir.mkdir(parents=True, exist_ok=True)
    sample = (Path(__file__).resolve().parents[2] / "packages" /
              "cubit-mesh-export" / "src" / "cubit_mesh_export" /
              "cubit_gui" / "solver_ready_sample.jou")
    if not sample.is_file():
        raise FileNotFoundError(sample)
    vol = run_dir / "em_sample_o2.vol"
    driver = run_dir / "driver.jou"
    driver.write_text(
        f'play "{sample.as_posix()}"\n'
        f'export netgen "{vol.as_posix()}" order 2 overwrite\nexit 0\n',
        encoding="utf-8",
    )
    proc = subprocess.run(
        [str(exe), "-batch", "-nographics", "-nojournal", str(driver)],
        cwd=run_dir, capture_output=True, text=True, encoding="utf-8",
        errors="replace", timeout=timeout,
    )
    (run_dir / "cubit.log").write_text(
        f"exit={proc.returncode}\n=== STDOUT ===\n{proc.stdout}\n=== STDERR ===\n{proc.stderr}",
        encoding="utf-8",
    )
    if not vol.is_file():
        raise RuntimeError(f"Cubit exit {proc.returncode}; missing {vol}")
    mesh = Mesh(str(vol))
    boundaries = sorted(set(mesh.GetBoundaries()))
    materials = sorted(set(mesh.GetMaterials()))
    areas = {name: float(Integrate(CoefficientFunction(1), mesh, BND,
                                  definedon=mesh.Boundaries(name))) if name in boundaries else 0.0
             for name in SAMPLE_BOUNDARIES}
    positive = {name for name, area in areas.items() if area > 0}
    gate = check_consistency(vol, strict_labels=True,
                             required_boundaries=tuple(sorted(SAMPLE_BOUNDARIES)),
                             required_materials=tuple(sorted(SAMPLE_MATERIALS)))
    return {
        "route": "cubit_export_netgen", "case": "gapped_coil_workpiece_air",
        "order": 2, "boundaries": boundaries, "materials": materials,
        "boundary_areas": areas,
        "boundary_label_preservation": len(positive & SAMPLE_BOUNDARIES)/len(SAMPLE_BOUNDARIES),
        "material_label_preservation": len(set(materials) & SAMPLE_MATERIALS)/len(SAMPLE_MATERIALS),
        "vol_check_passed": gate["passed"],
        "invalid_jacobian_samples": gate["quality"]["invalid_jacobian_sample_count"],
        "cubit_exit": proc.returncode, "run_dir": str(run_dir),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cubit-exe", type=Path,
                        default=Path(r"C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.exe"))
    parser.add_argument("--output", type=Path,
                        default=Path(r"C:\temp\cubit-paper-sphere-benchmark"))
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()
    if args.repeats < 2:
        parser.error("at least two independent runs are required for reproducibility")
    if not args.cubit_exe.is_file():
        parser.error(f"Cubit executable not found: {args.cubit_exe}")
    args.output.mkdir(parents=True, exist_ok=True)
    rows = []
    label_rows = []
    with TaskManager():
        for repeat in range(args.repeats):
            for kind in ("tet", "hex"):
                print(f"Cubit {kind} repeat {repeat+1}/{args.repeats}", flush=True)
                rows.extend({**row, "repeat": repeat+1} for row in run_cubit(
                    args.cubit_exe, kind, args.output / f"cubit_{kind}_{repeat+1}", args.timeout))
            print(f"Netgen/OCC tet repeat {repeat+1}/{args.repeats}", flush=True)
            rows.extend({**row, "repeat": repeat+1} for row in run_occ(
                args.output / f"occ_tet_{repeat+1}"))
            print(f"Cubit multi-label EM sample repeat {repeat+1}/{args.repeats}", flush=True)
            label_rows.append({**run_label_case(args.cubit_exe,
                args.output / f"em_labels_{repeat+1}", args.timeout), "repeat": repeat+1})
    result = {
        "protocol": "sphere radius=1; Cubit tet max size=0.65, hex sphere size=0.4; Netgen OCC maxh=0.65; p=1..3; repeated independent mesh generations",
        "field_problem": "electrostatic -Delta(phi)=3*pi^2*phi; phi=sin(pi*x)cos(pi*y)cos(pi*z); exact full-boundary Dirichlet; E=-grad(phi)",
        "environment": {"python": platform.python_version(),
                        "platform": platform.platform(),
                        "ngsolve": importlib.metadata.version("ngsolve"),
                        "cubit_mesh_export": importlib.metadata.version("cubit-mesh-export")},
        "cubit_exe": str(args.cubit_exe), "rows": rows,
        "multi_label_rows": label_rows,
    }
    target = args.output / "results.json"
    target.write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Wrote {target}")


if __name__ == "__main__":
    main()
