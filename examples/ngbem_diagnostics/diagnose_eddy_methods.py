"""
Diagnostic script for eddy current solver consistency.

Runs each solver with detailed intermediate output to identify
root causes of inconsistency. All methods run in separate subprocesses
to prevent BEM solver segfaults.

Validated findings (2026-02-15):
  - VectorFEMBEM: P=7036 W at 10 kHz, |A_scat/A_inc|=1.01 (strong shielding)
  - ShieldBEMSIBC: P=6742 W at 10 kHz, consistent with Vector (4.2% diff)
  - ScalarSIBC: Hz_total/Hz_inc decreases with freq (0.58 -> 0.30 -> 0.12)
    This is fundamental 3D limitation, not a bug
  - ScalarFEMBEM FEM-only: P=4808 W (~68% of Vector, Hz-only captures z-component)
"""
import sys
import os
import subprocess
import json
import numpy as np

MU_0 = 4e-7 * np.pi

# Geometry and material
WIDTH, HEIGHT, DEPTH = 0.02, 0.02, 0.01  # 20x20x10 mm
SIGMA = 3.7e7  # Al
MU_R = 1.0
FREQ = 10e3  # Single frequency for diagnostics

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(SCRIPT_DIR, '../../src/radia')


def analytical_reference(freq, sigma, mu_r, W, H, D):
    """Half-space analytical reference for rectangular block."""
    omega = 2 * np.pi * freq
    mu = mu_r * MU_0
    delta = np.sqrt(2 / (omega * mu * sigma))
    Rs = np.sqrt(omega * mu / (2 * sigma))
    H0 = 1.0 / MU_0  # Hz_inc for Bz=1T

    # Side faces only (Hz is tangential there)
    side_area = 2 * (W + H) * D
    P_sides = 0.5 * Rs * H0**2 * side_area

    return {
        'delta': delta,
        'Rs': Rs,
        'H0': H0,
        'side_area': side_area,
        'P_sides': P_sides,
    }


def run_in_subprocess(label, code):
    """Run code in subprocess and return JSON results."""
    import textwrap
    indented_code = textwrap.indent(textwrap.dedent(code), '    ')
    full_code = f"""
import sys, os, json, traceback
import numpy as np
sys.path.insert(0, r'{SRC_DIR}')

import io
_stderr = sys.stderr
sys.stderr = io.StringIO()

MU_0 = 4e-7 * np.pi

try:
    from ngbem_eddy import *
    from ngsolve import BND

    mesh = create_conductor_mesh({WIDTH}, {HEIGHT}, {DEPTH}, maxh=0.005)
    n_vol = sum(1 for _ in mesh.Elements())
    n_bnd = sum(1 for _ in mesh.Elements(BND))

    freq = {FREQ}
    sigma = {SIGMA}
    mu_r = {MU_R}

{indented_code}

    sys.stderr = _stderr
    sys.stdout.write("JSONRESULT:" + json.dumps(results) + "\\n")
    sys.stdout.flush()
except Exception as e:
    sys.stderr = _stderr
    sys.stdout.write("JSONRESULT:" + json.dumps({{"error": str(e), "traceback": traceback.format_exc()}}) + "\\n")
    sys.stdout.flush()
"""
    try:
        result = subprocess.run(
            [sys.executable, '-c', full_code],
            capture_output=True, text=True, timeout=300,
            cwd=SCRIPT_DIR
        )
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            if line.startswith('JSONRESULT:'):
                return json.loads(line[len('JSONRESULT:'):])
        if result.stderr:
            print(f"  [{label}] stderr: {result.stderr[:500]}")
        if result.stdout:
            print(f"  [{label}] stdout: {result.stdout[:500]}")
        return None
    except Exception as e:
        print(f"  [{label}] Exception: {e}")
        return None


def main():
    print("=" * 70)
    print("EDDY CURRENT SOLVER DIAGNOSTICS")
    print("=" * 70)
    print(f"  Geometry: {WIDTH*1e3:.0f}x{HEIGHT*1e3:.0f}x{DEPTH*1e3:.0f} mm "
          f"Al block (sigma={SIGMA:.1e} S/m)")
    print(f"  Frequency: {FREQ/1e3:.0f} kHz")

    ana = analytical_reference(FREQ, SIGMA, MU_R, WIDTH, HEIGHT, DEPTH)
    print(f"  Skin depth: {ana['delta']*1e3:.3f} mm")
    print(f"  Rs = {ana['Rs']:.6f} Ohm/sq")
    print(f"  H0 = {ana['H0']:.0f} A/m")
    print(f"  P_analytical (half-space, sides) = {ana['P_sides']:.2e} W")
    print(f"  maxh = 5 mm (better resolution)")
    print()

    # --- Test 1: VectorFEMBEM (reference method) ---
    print("--- VectorEddyCurrentFEMBEM (3D vector A, maxh=5mm) ---")
    r = run_in_subprocess("VectorFEMBEM", """
solver = VectorEddyCurrentFEMBEM(mesh, sigma, mu_r=mu_r, order=1)
solver.assemble(kappa=0.01, intorder=4)
res = solver.solve(freq, B_ext=[0, 0, 1.0])

# Get detailed diagnostics
A_total = solver._A_coeffs + solver._A_inc_coeffs
A_max = float(np.max(np.abs(A_total)))
A_scat_max = float(np.max(np.abs(solver._A_coeffs)))
A_inc_max = float(np.max(np.abs(solver._A_inc_coeffs)))

P = solver.compute_loss()

results = {
    "method": "VectorFEMBEM",
    "n_vol": n_vol, "n_bnd": n_bnd,
    "P": P,
    "n_dof": res["n_total"],
    "cond": res["cond"],
    "A_max_total": A_max,
    "A_max_scat": A_scat_max,
    "A_max_inc": A_inc_max,
    "A_ratio_scat_inc": A_scat_max / max(A_inc_max, 1e-30),
}
""")
    if r:
        if 'error' in r:
            print(f"  ERROR: {r['error']}")
        else:
            print(f"  P = {r['P']:.4e} W")
            print(f"  |A_total|_max = {r['A_max_total']:.4e}")
            print(f"  |A_scat|_max = {r['A_max_scat']:.4e}, "
                  f"|A_inc|_max = {r['A_max_inc']:.4e}")
            print(f"  |A_scat/A_inc| = {r['A_ratio_scat_inc']:.4f}")
            print(f"  DOFs = {r['n_dof']}, cond = {r['cond']:.2e}")
    else:
        print("  FAILED")
    print()

    # --- Test 2: ShieldBEMSIBC (loop basis) ---
    print("--- ShieldBEMSIBC (loop basis + SIBC) ---")
    r2 = run_in_subprocess("ShieldSIBC", """
solver = ShieldBEMSIBC(mesh, sigma, mu_r=mu_r)
solver.assemble(intorder=4)

def A_inc_uniform(points):
    A = np.zeros_like(points)
    A[:, 0] = -0.5 * points[:, 1]  # Ax = -0.5 * Bz * y
    A[:, 1] = 0.5 * points[:, 0]   # Ay = 0.5 * Bz * x
    return A

res = solver.solve(freq, A_inc_func=A_inc_uniform)
P = solver.compute_loss()

# Diagnostics
J = solver._J_full
MJ = solver._M_full @ J
energy_form = float(np.real(J.conj() @ MJ))
Zs = solver._Zs

results = {
    "method": "ShieldSIBC",
    "P": P,
    "n_loops": solver._loop.n_loops,
    "n_active": solver._loop.n_active,
    "J_max": float(np.max(np.abs(J))),
    "J_mean": float(np.mean(np.abs(J))),
    "Zs_real": float(Zs.real),
    "Zs_imag": float(Zs.imag),
    "delta_mm": float(solver.delta * 1e3),
    "energy_form": energy_form,
}
""")
    if r2:
        if 'error' in r2:
            print(f"  ERROR: {r2['error']}")
        else:
            print(f"  P = {r2['P']:.4e} W")
            print(f"  |J|_max = {r2['J_max']:.4e}, |J|_mean = {r2['J_mean']:.4e}")
            print(f"  Re(Zs) = {r2['Zs_real']:.6f}, Im(Zs) = {r2['Zs_imag']:.6f}")
            print(f"  delta = {r2['delta_mm']:.3f} mm")
            print(f"  J^H*M*J = {r2['energy_form']:.4e}")
            print(f"  Loops = {r2['n_loops']}, Active DOFs = {r2['n_active']}")
    else:
        print("  FAILED")
    print()

    # --- Test 3: ScalarSIBC diagnostics ---
    print("--- EddyCurrentBEMSIBC (scalar Hz + SIBC) DIAGNOSTICS ---")
    r3 = run_in_subprocess("ScalarSIBC", """
solver = EddyCurrentBEMSIBC(mesh, sigma=sigma, mu_r=mu_r, order=2)
solver.assemble(intorder=8)
solver.solve(freq, B_ext=[0, 0, 1.0])
P = solver.compute_loss()

# Detailed diagnostics
Hz_scat = solver._Hz_scat_surf
Hz_inc_val = solver._Hz_inc
Hz_inc_surf = np.full(solver._n_surf, Hz_inc_val, dtype=complex)
Hz_total = Hz_scat + Hz_inc_surf
gamma = solver._gamma
Zs = gamma / sigma

# Check Hz_total statistics
Hz_total_abs = np.abs(Hz_total)
Hz_scat_abs = np.abs(Hz_scat)

# M_loss quadratic form
M_loss_form = float(np.real(Hz_total.conj() @ solver._M_loss_np @ Hz_total))
M_surf_form = float(np.real(Hz_total.conj() @ solver._M_surf_np @ Hz_total))

# Check M_loss vs M_surf (should differ by nz^2 weighting)
M_loss_trace = float(np.real(np.trace(solver._M_loss_np)))
M_surf_trace = float(np.real(np.trace(solver._M_surf_np)))

results = {
    "method": "ScalarSIBC",
    "P": P,
    "Hz_inc": Hz_inc_val,
    "Hz_scat_max": float(Hz_scat_abs.max()),
    "Hz_scat_mean": float(Hz_scat_abs.mean()),
    "Hz_total_max": float(Hz_total_abs.max()),
    "Hz_total_mean": float(Hz_total_abs.mean()),
    "Hz_total_min": float(Hz_total_abs.min()),
    "Hz_ratio_total_inc": float(Hz_total_abs.mean() / abs(Hz_inc_val)),
    "Hz_scat_over_inc": float(Hz_scat_abs.mean() / abs(Hz_inc_val)),
    "Zs_real": float(Zs.real),
    "Zs_imag": float(Zs.imag),
    "gamma_abs": float(abs(gamma)),
    "M_loss_quad": M_loss_form,
    "M_surf_quad": M_surf_form,
    "M_loss_trace": M_loss_trace,
    "M_surf_trace": M_surf_trace,
    "n_surf": solver._n_surf,
    "n_l2": solver._n_l2,
}
""")
    if r3:
        if 'error' in r3:
            print(f"  ERROR: {r3['error']}")
        else:
            print(f"  P = {r3['P']:.4e} W")
            print(f"  Hz_inc = {r3['Hz_inc']:.0f} A/m")
            print(f"  |Hz_scat|: max={r3['Hz_scat_max']:.0f}, "
                  f"mean={r3['Hz_scat_mean']:.0f}")
            print(f"  |Hz_total|: max={r3['Hz_total_max']:.0f}, "
                  f"mean={r3['Hz_total_mean']:.0f}, "
                  f"min={r3['Hz_total_min']:.0f}")
            print(f"  |Hz_total/Hz_inc| = {r3['Hz_ratio_total_inc']:.6f}")
            print(f"  |Hz_scat/Hz_inc| = {r3['Hz_scat_over_inc']:.6f}")
            print(f"  Re(Zs) = {r3['Zs_real']:.6f}")
            print(f"  |gamma| = {r3['gamma_abs']:.2f}")
            print(f"  Hz^H*M_loss*Hz = {r3['M_loss_quad']:.4e}")
            print(f"  Hz^H*M_surf*Hz = {r3['M_surf_quad']:.4e}")
            print(f"  Tr(M_loss) = {r3['M_loss_trace']:.6e}")
            print(f"  Tr(M_surf) = {r3['M_surf_trace']:.6e}")
            print(f"  Surface DOFs = {r3['n_surf']}, L2 DOFs = {r3['n_l2']}")
    else:
        print("  FAILED")
    print()

    # --- Test 4: ScalarSIBC at multiple frequencies (to see trend) ---
    print("--- ScalarSIBC frequency trend ---")
    r4 = run_in_subprocess("ScalarSIBC_sweep", """
solver = EddyCurrentBEMSIBC(mesh, sigma=sigma, mu_r=mu_r, order=2)
solver.assemble(intorder=8)

freqs_test = [1e3, 10e3, 100e3]
items = []
for f in freqs_test:
    solver.solve(f, B_ext=[0, 0, 1.0])
    P = solver.compute_loss()

    Hz_scat = solver._Hz_scat_surf
    Hz_inc_surf = np.full(solver._n_surf, solver._Hz_inc, dtype=complex)
    Hz_total = Hz_scat + Hz_inc_surf

    items.append({
        "freq": f,
        "P": P,
        "Hz_total_mean": float(np.mean(np.abs(Hz_total))),
        "Hz_scat_mean": float(np.mean(np.abs(Hz_scat))),
        "Hz_total_over_inc": float(np.mean(np.abs(Hz_total)) / abs(solver._Hz_inc)),
        "Zs_real": float((solver._gamma / sigma).real),
        "delta_mm": float(solver.delta * 1e3),
    })

results = {"method": "ScalarSIBC_sweep", "items": items}
""")
    if r4:
        if 'error' in r4:
            print(f"  ERROR: {r4['error']}")
        else:
            print(f"  {'freq':>10}  {'P (W)':>12}  {'|Hz_t|/|Hz_i|':>14}  "
                  f"{'Re(Zs)':>10}  {'delta':>8}")
            for item in r4['items']:
                print(f"  {item['freq']/1e3:8.1f} kHz  {item['P']:12.4e}  "
                      f"{item['Hz_total_over_inc']:14.6f}  "
                      f"{item['Zs_real']:10.6f}  "
                      f"{item['delta_mm']:6.3f} mm")
    else:
        print("  FAILED")
    print()

    # --- Test 5: ScalarFEMBEM (Costabel) with FEM-only mode as reference ---
    print("--- EddyCurrentFEMBEM: FEM-only vs FEM-BEM ---")
    r5 = run_in_subprocess("ScalarFEMBEM", """
items = []
for mode_name, mode in [("FEM-only", "fem"), ("FEM-BEM", "fembem")]:
    try:
        solver = EddyCurrentFEMBEM(mesh, sigma=sigma, mu_r=mu_r, order=2)
        if mode == "fem":
            solver.assemble_fem(freq)
        else:
            solver.assemble_fembem(freq, intorder=8)
        solver.solve(B_ext=[0, 0, 1.0], mode=mode)
        P = solver.compute_loss()
        items.append({"mode": mode_name, "P": P, "freq": freq})
    except Exception as e:
        items.append({"mode": mode_name, "P": None, "error": str(e), "freq": freq})

results = {"method": "ScalarFEMBEM_modes", "items": items}
""")
    if r5:
        if 'error' in r5:
            print(f"  ERROR: {r5['error']}")
        else:
            for item in r5['items']:
                if item.get('P') is not None:
                    print(f"  {item['mode']:12s}: P = {item['P']:.4e} W")
                else:
                    print(f"  {item['mode']:12s}: FAILED - {item.get('error', 'unknown')}")
    else:
        print("  FAILED")
    print()

    # --- Test 6: VectorFEMBEM mesh convergence ---
    print("--- VectorFEMBEM mesh convergence (freq=10 kHz) ---")
    for maxh_mm in [10, 5]:  # maxh=3mm times out (>5min for dense solve)
        r6 = run_in_subprocess(f"Vector_maxh{maxh_mm}", f"""
mesh2 = create_conductor_mesh({WIDTH}, {HEIGHT}, {DEPTH}, maxh={maxh_mm/1000.0})
n_el = sum(1 for _ in mesh2.Elements())

solver = VectorEddyCurrentFEMBEM(mesh2, sigma, mu_r=mu_r, order=1)
solver.assemble(kappa=0.01, intorder=4)
res = solver.solve(freq, B_ext=[0, 0, 1.0])
P = solver.compute_loss()

results = {{"maxh_mm": {maxh_mm}, "n_elem": n_el, "P": P, "ndof": res["n_total"]}}
""")
        if r6:
            if 'error' in r6:
                print(f"  maxh={maxh_mm} mm: ERROR - {r6['error']}")
            else:
                print(f"  maxh={maxh_mm:2d} mm: n_elem={r6['n_elem']:5d}, "
                      f"P={r6['P']:.4e} W, DOFs={r6['ndof']}")
        else:
            print(f"  maxh={maxh_mm} mm: FAILED")
    print()

    # --- Summary ---
    print("=" * 70)
    print("DIAGNOSTIC SUMMARY")
    print("=" * 70)
    print(f"  Analytical P (half-space, sides): {ana['P_sides']:.2e} W")
    print(f"  Note: Finite block P < half-space P due to demagnetizing effect")
    print()
    print("  Key questions to answer:")
    print("  1. Does VectorFEMBEM converge with mesh refinement?")
    print("  2. Does ShieldBEMSIBC give same order of magnitude?")
    print("  3. Is ScalarSIBC Hz_total -> 0 (formulation issue)?")
    print("  4. Does ScalarFEMBEM FEM-only mode work at 10 kHz?")


if __name__ == '__main__':
    main()
