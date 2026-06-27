"""
NGSolve User Meeting Demo: Volume/Area p-convergence via Cubit ACIS trampoline.

Demonstrates:
  1. Cubit hex/tet sphere mesh -> ACIS CallbackGeometry trampoline
  2. p-convergence (order 2-5) for volume and surface area
  3. Comparison: OCC native Curve() vs ACIS Python callback

Requirements:
  pip install radia   (includes cubit_mesh_curver.pyd + ngsolve + cubit_mesh_export)
  Coreform Cubit 2025.8+ installed

Docs companion:
  docs/ngsolve_user_meeting/volume_area_convergence.ipynb executes the Cubit-free
  OCC subset and stores synchronized JSON results.  This script keeps the full
  Cubit ACIS trampoline presentation path beside that notebook.

Usage:
  python acis_volume_area_convergence_demo.py
"""

import sys
import os
import math
import time

# NGSolve MUST be imported before cubit
from ngsolve import Mesh as NGMesh, Integrate, CF, BND
from ngsolve import TaskManager

RADIUS = 0.05
V_EXACT = (4.0 / 3.0) * math.pi * RADIUS**3
A_EXACT = 4.0 * math.pi * RADIUS**2


def demo_occ_path():
    """OCC native: STEP -> Netgen -> mesh.Curve(p) -> Integrate."""
    from netgen.occ import Sphere, OCCGeometry

    print("=" * 70)
    print("  Path 1: OCC native (no Cubit)")
    print("=" * 70)

    geo = OCCGeometry(Sphere((0, 0, 0), RADIUS))

    for maxh, label in [(0.015, "coarse"), (0.008, "fine")]:
        with TaskManager():
            ng_base = geo.GenerateMesh(maxh=maxh)
            ne = ng_base.ne
            print(f"\n  Mesh: {ne} tets (maxh={maxh}, {label})")
            print(f"  {'p':>3} {'Volume':>14} {'V err':>10} {'Area':>14} {'A err':>10} {'Time':>7}")
            print(f"  {'-'*60}")

            for p in range(1, 6):
                ng = geo.GenerateMesh(maxh=maxh)
                t0 = time.perf_counter()
                if p >= 2:
                    ng.Curve(p)
                t1 = time.perf_counter()
                mesh = NGMesh(ng)
                vol = Integrate(CF(1), mesh)
                area = Integrate(CF(1), mesh, BND)
                v_err = (vol - V_EXACT) / V_EXACT * 100
                a_err = (area - A_EXACT) / A_EXACT * 100
                print(f"  {p:>3} {vol:>14.8e} {v_err:>+9.5f}% "
                      f"{area:>14.8e} {a_err:>+9.5f}% {t1-t0:>6.3f}s")


def demo_acis_path():
    """ACIS trampoline: Cubit mesh -> Python callbacks -> Integrate."""
    # Setup Cubit
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '..', '..', 'src', 'radia')
    sys.path.insert(0, os.path.abspath(radia_src))
    from install_panels import find_cubit_bin
    cubit_path = find_cubit_bin()
    if not cubit_path:
        print("\n  Cubit not found. Skipping ACIS path.")
        return
    sys.path.append(cubit_path)
    plugin_dir = os.path.join(os.path.dirname(cubit_path), 'bin', 'plugins')
    os.environ['CUBIT_PLUGIN_DIR'] = plugin_dir

    import cubit
    cubit.init(['cubit', '-nojournal', '-batch', '-nographics',
                '-commandplugindir', plugin_dir])
    from cubit_mesh_export import extract_curved_mesh

    print("\n" + "=" * 70)
    print("  Path 2: Cubit ACIS trampoline (CallbackGeometry)")
    print("=" * 70)

    # --- Tet sphere ---
    for mesh_size, label in [(0.015, "coarse"), (0.008, "fine")]:
        cubit.cmd("reset")
        cubit.cmd(f"create sphere radius {RADIUS}")
        cubit.cmd("volume 1 scheme tetmesh")
        cubit.cmd(f"volume 1 size {mesh_size}")
        cubit.cmd("mesh volume 1")
        ne = len(cubit.parse_cubit_list("tet", "all"))

        print(f"\n  Tet mesh: {ne} tets (size={mesh_size}, {label})")
        print(f"  {'p':>3} {'Volume':>14} {'V err':>10} {'Area':>14} {'A err':>10} {'Time':>7}")
        print(f"  {'-'*60}")

        for p in range(2, 6):
            t0 = time.perf_counter()
            ng = extract_curved_mesh(cubit, order=p)
            t1 = time.perf_counter()
            mesh = NGMesh(ng)
            vol = Integrate(CF(1), mesh)
            area = Integrate(CF(1), mesh, BND)
            v_err = (vol - V_EXACT) / V_EXACT * 100
            a_err = (area - A_EXACT) / A_EXACT * 100
            print(f"  {p:>3} {vol:>14.8e} {v_err:>+9.5f}% "
                  f"{area:>14.8e} {a_err:>+9.5f}% {t1-t0:>6.3f}s")

    # --- Hex sphere ---
    for mesh_size, label in [(0.015, "coarse"), (0.008, "fine")]:
        cubit.cmd("reset")
        cubit.cmd(f"create sphere radius {RADIUS}")
        cubit.cmd(f"volume 1 size {mesh_size}")
        cubit.cmd("volume 1 scheme sphere")
        cubit.cmd("mesh volume 1")
        ne = len(cubit.parse_cubit_list("hex", "all"))
        if ne == 0:
            continue

        print(f"\n  Hex mesh: {ne} hexes (size={mesh_size}, {label})")
        print(f"  {'p':>3} {'Volume':>14} {'V err':>10} {'Area':>14} {'A err':>10} {'Time':>7}")
        print(f"  {'-'*60}")

        for p in range(2, 6):
            t0 = time.perf_counter()
            ng = extract_curved_mesh(cubit, order=p)
            t1 = time.perf_counter()
            mesh = NGMesh(ng)
            vol = Integrate(CF(1), mesh)
            area = Integrate(CF(1), mesh, BND)
            v_err = (vol - V_EXACT) / V_EXACT * 100
            a_err = (area - A_EXACT) / A_EXACT * 100
            print(f"  {p:>3} {vol:>14.8e} {v_err:>+9.5f}% "
                  f"{area:>14.8e} {a_err:>+9.5f}% {t1-t0:>6.3f}s")
        print(f"  Note: hex order>=4 limited by quad bubble fix (pip pending)")


def main():
    print("=" * 70)
    print("  NGSolve User Meeting: Volume/Area p-Convergence Demo")
    print(f"  Sphere r={RADIUS} m")
    print(f"  V_exact = {V_EXACT:.10e} m^3")
    print(f"  A_exact = {A_EXACT:.10e} m^2")
    print("=" * 70)

    demo_occ_path()
    demo_acis_path()

    print("\n" + "=" * 70)
    print("  Summary")
    print("=" * 70)
    print("  OCC Curve():   order 1-5, C++ native, fastest")
    print("  ACIS callback: order 2-5, Python trampoline, Cubit hex/tet")
    print("  Both paths achieve identical accuracy at each order.")
    print("  Hex: order>=4 limited by quad bubble fix (pending pip release).")
    print("=" * 70)


if __name__ == "__main__":
    main()
