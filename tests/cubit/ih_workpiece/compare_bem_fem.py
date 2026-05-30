"""Compare BEM-SIBC vs FEM-SIBC on the same IH geometry.

Generates a gapped-torus coil + cylinder workpiece + air sphere via OCC,
saves as .vol, then runs both solvers and compares P_total.

No Cubit required. No Kelvin transform (truncated air sphere + Dirichlet).
"""

import json
import os
import subprocess
import sys
import tempfile
import time

_this_dir = os.path.dirname(os.path.abspath(__file__))
_panels_dir = os.path.abspath(os.path.join(_this_dir, '..', '..', '..', 'src', 'radia', 'panels'))
sys.path.insert(0, _panels_dir)


def generate_mesh(vol_path, R_coil=0.030, a_coil=0.003, gap_deg=5,
                  R_wp=0.010, H_wp=0.020, R_air=0.080,
                  maxh_coil=0.005, maxh_air=0.012, order=2):
    """Generate OCC mesh with coil+workpiece+air."""
    from calc_common import build_occ_ih_mesh_simple
    mesh, info = build_occ_ih_mesh_simple(
        R_coil=R_coil, a_coil=a_coil, gap_deg=gap_deg,
        R_wp=R_wp, H_wp=H_wp, R_air=R_air,
        maxh_coil=maxh_coil, maxh_air=maxh_air, order=order)
    mesh.ngmesh.Save(vol_path)

    from ngsolve import VOL, BND
    print(f"Mesh: {mesh.GetNE(VOL)} vol, {mesh.GetNE(BND)} surf")
    print(f"Materials: {list(mesh.GetMaterials())}")
    print(f"Boundaries: {list(set(mesh.GetBoundaries()))}")
    return info


def run_bem(vol_path, frequency, sigma, mu_r, R_wp, tmpdir):
    """Run BEM-SIBC with LINEAR Z_s (same as FEM) for fair comparison.

    Uses bem_sibc_solver directly (not subprocess) to control Z_s exactly.
    Coil = Biot-Savart filament (same torus geometry, zero cross-section).
    """
    import math
    import numpy as np

    sys.path.insert(0, os.path.join(_panels_dir, '..'))
    from calc_common import MU_0, write_surface_only_vol

    print(f"\n=== BEM-SIBC (linear Z_s, Biot-Savart coil) ===")
    t0 = time.perf_counter()

    # Extract workpiece surface
    wp_vol = os.path.join(tmpdir, 'workpiece_only.vol')
    _extract_workpiece_vol(vol_path, wp_vol)

    from ngsolve import Mesh, Integrate, CF, BND
    mesh = Mesh(wp_vol)

    # Linear SIBC Z_s (same formula as FEM)
    omega = 2 * math.pi * frequency
    rho = 1.0 / sigma
    mu_eff = MU_0 * mu_r
    delta = math.sqrt(2 * rho / (omega * mu_eff))
    Z_s = complex(1, 1) * rho / delta
    print(f"Z_s = {Z_s:.6e}, |Z_s| = {abs(Z_s):.6e}, delta = {delta*1e3:.4f}mm")

    from calc_common import setup_paths
    setup_paths()
    from radia.bem_sibc_solver import ScalarBIESIBCSolver, compute_phi_inc_from_loop

    # BEM assembly
    solver = ScalarBIESIBCSolver(mesh, order=1)
    print(f"BEM: {solver.ndof} DOFs, assembly {solver.t_assembly:.1f}s")

    # phi_inc from Biot-Savart loop (same geometry as OCC torus)
    node_coords = np.array([[mesh.vertices[i].point[j] for j in range(3)]
                            for i in range(mesh.nv)])
    phi_inc = compute_phi_inc_from_loop(
        node_coords, loop_center=[0, 0, 0], loop_radius=0.030,
        current=1.0, n_quad=30, gap_deg=5)

    # Single solve with fixed linear Z_s (no Karl iteration)
    result = solver.solve(phi_inc, Z_s=Z_s, omega=omega)
    t_bem = time.perf_counter() - t0

    H_t_rms = result['H_t_rms']
    P_density = result['P_density']
    area = result['area']
    P_total = P_density * area

    out = {
        'P_total_W': float(P_total),
        'H_t_rms_Am': float(H_t_rms),
        'Z_s_abs_Ohm': float(abs(Z_s)),
        'Z_s_str': str(Z_s),
        'ndof': solver.ndof,
        'delta_mm': float(delta * 1e3),
        'area_m2': float(area),
        't_wall': t_bem,
    }
    print(f"P_total = {P_total:.6e} W")
    print(f"H_t_rms = {H_t_rms:.4f} A/m")
    print(f"Area = {area:.6e} m^2")
    print(f"time = {t_bem:.1f}s")
    return out


def _extract_workpiece_vol(full_vol, wp_vol):
    """Extract workpiece from full volume mesh as its own mesh."""
    # For BEM, we just need the workpiece surface.
    # Load the full mesh, extract workpiece boundary elements.
    from ngsolve import Mesh, BND
    from ngsolve import TaskManager
    import netgen.meshing as ngm

    mesh = Mesh(full_vol)

    # Create a new mesh with only wp_surface boundary elements
    ngmesh_new = ngm.Mesh(dim=3)

    bnd_labels = list(mesh.GetBoundaries())
    old_to_new = {}

    # Add face descriptor for wp_surface
    fd = ngm.FaceDescriptor(bc=1)
    ngmesh_new.Add(fd)
    ngmesh_new.SetBCName(0, "wp_surface")

    for el in mesh.Elements(BND):
        lbl = bnd_labels[el.index]
        if lbl != "wp_surface":
            continue
        new_verts = []
        for v in el.vertices:
            if v.nr not in old_to_new:
                pt = mesh.vertices[v.nr].point
                old_to_new[v.nr] = ngmesh_new.Add(
                    ngm.MeshPoint(ngm.Pnt(pt[0], pt[1], pt[2])))
            new_verts.append(old_to_new[v.nr])
        se = ngm.Element2D(1, new_verts)
        ngmesh_new.Add(se)

    wp_mesh = Mesh(ngmesh_new)
    with TaskManager():
        wp_mesh.Curve(2)
        wp_mesh.ngmesh.Save(wp_vol)
        print(f"Extracted workpiece: {wp_mesh.GetNE(BND)} elements, "
              f"{wp_mesh.nv} vertices")


def run_fem(vol_path, frequency, sigma, material, mu_r, R_wp, tmpdir):
    """Run FEM-SIBC (T0 source/sink + Robin BC on wp_surface)."""
    calc_script = os.path.join(_panels_dir, 'calc_fem_kelvin.py')

    cmd = [
        sys.executable, calc_script,
        "--vol", vol_path,
        "--fes-order", "1",
        "--frequency", str(frequency),
        "--sigma", str(sigma),
        "--impedance", "sibc",
        "--mu-r", str(mu_r),
        "--material", material,
        "--current", "1.0",
        "--a-coil", "0.003",
        "--r-wp", str(R_wp),
        "--max-iter", "1",
        "--solver", "pardiso",
        "--msh-output", os.path.join(tmpdir, "fem_result.msh"),
    ]
    print(f"\n=== FEM-SIBC (no Kelvin) ===")
    print(f"CMD: {' '.join(cmd[-16:])}")
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    t_fem = time.perf_counter() - t0

    if proc.returncode != 0:
        print(f"FEM FAILED:")
        print(f"stderr: {proc.stderr[-1000:]}")
        print(f"stdout: {proc.stdout[-500:]}")
        return None

    result = json.loads(proc.stdout.strip().split("\n")[-1])
    if "error" in result:
        print(f"FEM ERROR: {result['error']}")
        return None

    result['t_wall'] = t_fem
    print(f"P_total = {result['P_total']:.6e} W")
    print(f"H_t_rms = {result['H_t_rms']:.4f} A/m")
    print(f"Z_s = {result['Z_s']}")
    print(f"L = {result['L']*1e9:.2f} nH")
    print(f"ndof = {result['ndof']}, iter = {result['iterations']}, "
          f"time = {t_fem:.1f}s")
    return result


def main():
    tmpdir = tempfile.mkdtemp(prefix="ih_compare_")
    vol_path = os.path.join(tmpdir, "ih_model.vol")

    # Common parameters
    frequency = 7000
    sigma = 2e6
    material = "steel"
    mu_r = 100.0
    R_wp = 0.010

    print("=" * 60)
    print("BEM-SIBC vs FEM-SIBC comparison")
    print("=" * 60)
    print(f"Frequency: {frequency} Hz")
    print(f"Material: {material}, sigma={sigma:.1e} S/m, mu_r={mu_r}")
    print(f"Workpiece: R={R_wp*1e3}mm cylinder")
    print(f"Coil: R=30mm torus, a=3mm, gap=5deg, I=1A")
    print(f"FEM: truncated air R=80mm (no Kelvin)")

    # Step 1: Generate mesh
    print("\n=== Mesh generation ===")
    info = generate_mesh(vol_path, R_wp=R_wp, maxh_coil=0.004, maxh_air=0.010)

    # Step 2: Run BEM-SIBC (linear Z_s, Biot-Savart coil)
    bem_result = run_bem(vol_path, frequency, sigma, mu_r, R_wp, tmpdir)

    # Step 3: Run FEM-SIBC (linear Z_s, T0 source/sink coil)
    fem_result = run_fem(vol_path, frequency, sigma, material, mu_r, R_wp,
                         tmpdir)

    # Step 4: Compare
    print("\n" + "=" * 60)
    print("COMPARISON")
    print("=" * 60)

    if bem_result and fem_result:
        P_bem = bem_result['P_total_W']
        P_fem = fem_result['P_total']
        H_bem = bem_result['H_t_rms_Am']
        H_fem = fem_result['H_t_rms']

        diff_P = (P_bem - P_fem) / max(abs(P_fem), 1e-30) * 100
        diff_H = (H_bem - H_fem) / max(abs(H_fem), 1e-30) * 100

        print(f"{'':20s} {'BEM-SIBC':>14s} {'FEM-SIBC':>14s} {'diff':>10s}")
        print(f"{'P_total [W]':20s} {P_bem:14.6e} {P_fem:14.6e} {diff_P:+9.1f}%")
        print(f"{'H_t_rms [A/m]':20s} {H_bem:14.4f} {H_fem:14.4f} {diff_H:+9.1f}%")
        print(f"Z_s BEM  = {bem_result['Z_s_str']}")
        print(f"Z_s FEM  = {fem_result['Z_s']}")
        print(f"delta = {bem_result['delta_mm']:.4f} mm")
        print(f"t_wall BEM = {bem_result['t_wall']:.1f}s, "
              f"FEM = {fem_result['t_wall']:.1f}s")

        # Expected: same Z_s (both linear), difference from coil model
        # BEM: Biot-Savart filament (zero cross-section)
        # FEM: T0 source/sink (finite cross-section a=3mm)
        if abs(diff_P) > 50:
            print("\nWARNING: P_total differs by >50%. Check coil model.")
    else:
        if not bem_result:
            print("BEM failed!")
        if not fem_result:
            print("FEM failed!")

    print(f"\nTemp dir: {tmpdir}")
    return 0 if (bem_result and fem_result) else 1


if __name__ == "__main__":
    sys.exit(main())
