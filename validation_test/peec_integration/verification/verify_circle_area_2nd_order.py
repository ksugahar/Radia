"""Verify curved Netgen/NGSolve area convergence on a circular disk.

Historical versions compared GMSH-generated Tri6 elements against Netgen.  The
current Radia policy is stricter: GMSH is visualization-only, while analysis
meshes come from Netgen/OCC or Cubit/Coreform and are passed to NGSolve as
native meshes or .vol files.

This script keeps the same validation goal: show that curving the mesh greatly
improves the geometric area error for a circle.
"""

from __future__ import annotations

import math


def netgen_circle_area(radius: float, mesh_size: float,
                       curve_order: int = 1) -> dict:
    """Create a Netgen/OCC disk mesh and integrate its area in NGSolve."""
    from netgen.occ import WorkPlane, OCCGeometry
    from netgen.meshing import MeshingParameters
    from ngsolve import CF, Integrate, Mesh, TaskManager

    wp = WorkPlane()
    wp.MoveTo(radius, 0).Arc(radius, 180).Arc(radius, 180)
    face = wp.Face()
    geo = OCCGeometry(face, dim=2)
    mp = MeshingParameters(maxh=mesh_size)

    with TaskManager():
        ngmesh = geo.GenerateMesh(mp)
        mesh = Mesh(ngmesh)
        if curve_order > 1:
            mesh.Curve(curve_order)
        area = Integrate(CF(1.0), mesh)

    exact_area = math.pi * radius**2
    error = abs(area - exact_area) / exact_area
    return {
        "area": float(area),
        "exact": float(exact_area),
        "error": float(error),
        "n_elements": int(mesh.ne),
        "order": int(curve_order),
        "method": f"Netgen Curve({curve_order})",
    }


def convergence_rate(results: list[dict]) -> float:
    """Estimate convergence rate from first and last rows using h ~ N^-1/2."""
    if len(results) < 2:
        return float("nan")
    first = results[0]
    last = results[-1]
    if first["error"] <= 0 or last["error"] <= 0:
        return float("nan")
    h_ratio = math.sqrt(first["n_elements"] / last["n_elements"])
    if h_ratio <= 0 or h_ratio == 1:
        return float("nan")
    return math.log(first["error"] / last["error"]) / math.log(1.0 / h_ratio)


def run_sweep(radius: float = 1.0) -> dict:
    mesh_sizes = [0.8, 0.5, 0.3, 0.2, 0.15, 0.1]
    orders = [1, 2, 3]
    groups: dict[int, list[dict]] = {order: [] for order in orders}

    for h in mesh_sizes:
        for order in orders:
            groups[order].append(netgen_circle_area(radius, h, order))

    return {
        "radius": radius,
        "exact_area": math.pi * radius**2,
        "groups": groups,
        "rates": {order: convergence_rate(rows)
                  for order, rows in groups.items()},
    }


def print_results(result: dict) -> None:
    print()
    print("Circle area accuracy: Netgen Curve(p)")
    print("=" * 70)
    print(f"  Radius = {result['radius']}, exact area = {result['exact_area']:.10f}")

    for order, rows in result["groups"].items():
        print()
        print(f"  Netgen Curve({order}):")
        print(f"    {'N_elem':>7s}  {'Area':>14s}  {'Rel error':>12s}")
        for row in rows:
            print(f"    {row['n_elements']:7d}  "
                  f"{row['area']:14.10f}  {row['error']:12.3e}")

    print()
    print("Convergence summary")
    print("=" * 70)
    print(f"  {'Method':<18s} {'N_min':>7s} {'Err_min':>12s} "
          f"{'N_max':>7s} {'Err_max':>12s} {'Rate':>8s}")
    for order, rows in result["groups"].items():
        first = rows[0]
        last = rows[-1]
        print(f"  {f'Netgen Curve({order})':<18s} "
              f"{first['n_elements']:7d} {first['error']:12.3e} "
              f"{last['n_elements']:7d} {last['error']:12.3e} "
              f"{result['rates'][order]:8.2f}")


def main() -> None:
    result = run_sweep(radius=1.0)
    print_results(result)

    final_linear = result["groups"][1][-1]["error"]
    final_curved = result["groups"][2][-1]["error"]
    if final_curved >= final_linear:
        raise SystemExit(
            "FAIL: Curve(2) should improve the final area error over Curve(1)"
        )
    print()
    print("Overall: PASS")


if __name__ == "__main__":
    main()
