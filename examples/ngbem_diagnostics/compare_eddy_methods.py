"""
Cross-method comparison: VectorFEMBEM vs ShieldBEMSIBC eddy current solvers.

Runs each solver independently (in subprocess) and compares loss results.
Subprocess isolation prevents segfaults from multiple BEM solver instances.

Validated results (2026-02-15):
  - VectorFEMBEM and ShieldBEMSIBC agree within 4.2% at 10 kHz
  - Each method has a complementary validity range:
    * VectorFEMBEM: accurate when maxh << skin depth (low-medium freq)
    * ShieldBEMSIBC: accurate when skin depth << conductor size (medium-high freq)
  - ScalarFEM (Hz only) gives ~60-70% of Vector loss (Hz-only limitation)
  - The crossover frequency (both valid) depends on mesh size and geometry
"""
import sys
import os
import subprocess
import json
import textwrap
import numpy as np

MU_0 = 4e-7 * np.pi

# Common parameters
WIDTH, HEIGHT, DEPTH = 0.02, 0.02, 0.01  # 20x20x10 mm
SIGMA = 3.7e7  # Al
MU_R = 1.0
MAXH = 0.005  # 5mm mesh (good balance speed/accuracy)
FREQS = [1e3, 3e3, 10e3, 30e3, 100e3]

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(SCRIPT_DIR, '../../src/radia')


def run_method(method_name, code):
    """Run a method in a subprocess and return results."""
    indented = textwrap.indent(textwrap.dedent(code), '    ')
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

    mesh = create_conductor_mesh({WIDTH}, {HEIGHT}, {DEPTH}, maxh={MAXH})
    n_vol = sum(1 for _ in mesh.Elements())
    n_bnd = sum(1 for _ in mesh.Elements(BND))

    freqs = {FREQS}
    sigma = {SIGMA}
    mu_r = {MU_R}

{indented}

    sys.stderr = _stderr
    sys.stdout.write("JSONRESULT:" + json.dumps(results) + "\\n")
    sys.stdout.flush()
except Exception as e:
    sys.stderr = _stderr
    sys.stdout.write("JSONRESULT:" + json.dumps({{"error": str(e), "tb": traceback.format_exc()}}) + "\\n")
    sys.stdout.flush()
"""
    try:
        result = subprocess.run(
            [sys.executable, '-c', full_code],
            capture_output=True, text=True, timeout=600,
            cwd=SCRIPT_DIR
        )
        # Find JSON line by marker
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            if line.startswith('JSONRESULT:'):
                return json.loads(line[len('JSONRESULT:'):])
        # Print stderr for debugging
        if result.stderr:
            print(f"  [{method_name}] stderr: {result.stderr[:500]}")
        if result.stdout:
            print(f"  [{method_name}] stdout: {result.stdout[:500]}")
        return None
    except subprocess.TimeoutExpired:
        print(f"  [{method_name}] Timeout (600s)")
        return None
    except Exception as e:
        print(f"  [{method_name}] Exception: {e}")
        return None


def main():
    print("=" * 70)
    print("Cross-Method Comparison: Eddy Current Solvers")
    print("=" * 70)
    print(f"  Geometry: {WIDTH*1e3:.0f}x{HEIGHT*1e3:.0f}x{DEPTH*1e3:.0f} mm "
          f"aluminum block")
    print(f"  sigma = {SIGMA:.1e} S/m, mu_r = {MU_R}")
    print(f"  maxh = {MAXH*1e3:.0f} mm")
    print(f"  Frequencies: {[f'{f/1e3:.0f} kHz' for f in FREQS]}")
    print()

    all_results = {}

    # --- Method 1: VectorEddyCurrentFEMBEM ---
    print("--- Method 1: VectorEddyCurrentFEMBEM (3D vector A) ---")
    r = run_method("VectorFEMBEM", """
        solver = VectorEddyCurrentFEMBEM(mesh, sigma, mu_r=mu_r, order=1)
        solver.assemble(kappa=0.01, intorder=4)

        results = {"method": "VectorFEMBEM", "losses": [],
                   "n_vol": n_vol, "n_bnd": n_bnd}
        for f in freqs:
            try:
                res = solver.solve(f, B_ext=[0, 0, 1.0])
                P = solver.compute_loss()
                results["losses"].append({
                    "freq": f, "P": P,
                    "cond": res["cond"],
                    "n_dof": res["n_total"],
                })
            except Exception as e:
                results["losses"].append({"freq": f, "P": None,
                                          "error": str(e)})
    """)
    if r:
        if 'error' in r:
            print(f"  ERROR: {r['error']}")
            print(f"  Traceback: {r.get('tb', '')[:300]}")
        else:
            all_results["VectorFEMBEM"] = r
            for item in r["losses"]:
                if item.get("P") is not None:
                    print(f"  f={item['freq']/1e3:8.1f} kHz  P={item['P']:.4e} W  "
                          f"cond={item['cond']:.2e}  DOFs={item['n_dof']}")
                else:
                    print(f"  f={item['freq']/1e3:8.1f} kHz  FAILED: "
                          f"{item.get('error', 'unknown')}")
    else:
        print("  FAILED")

    # --- Method 2: ShieldBEMSIBC (loop basis) ---
    print("\n--- Method 2: ShieldBEMSIBC (loop basis + SIBC) ---")
    r = run_method("ShieldSIBC", """
        solver = ShieldBEMSIBC(mesh, sigma, mu_r=mu_r)
        solver.assemble(intorder=4)

        n_loops = solver._loop.n_loops
        results = {"method": "ShieldSIBC", "losses": [],
                   "n_loops": n_loops}
        for f in freqs:
            try:
                def A_inc_uniform(points):
                    # A for uniform Bz=1T: A = 0.5 * B x r
                    A = np.zeros_like(points)
                    A[:, 0] = -0.5 * points[:, 1]  # Ax = -0.5 * Bz * y
                    A[:, 1] = 0.5 * points[:, 0]   # Ay = 0.5 * Bz * x
                    return A
                solver.solve(f, A_inc_func=A_inc_uniform)
                P = solver.compute_loss()
                results["losses"].append({"freq": f, "P": P})
            except Exception as e:
                results["losses"].append({"freq": f, "P": None,
                                          "error": str(e)})
    """)
    if r:
        if 'error' in r:
            print(f"  ERROR: {r['error']}")
        else:
            all_results["ShieldSIBC"] = r
            for item in r["losses"]:
                if item.get("P") is not None:
                    print(f"  f={item['freq']/1e3:8.1f} kHz  P={item['P']:.4e} W")
                else:
                    print(f"  f={item['freq']/1e3:8.1f} kHz  FAILED: "
                          f"{item.get('error', 'unknown')}")
            print(f"  Loops: {r.get('n_loops', 'N/A')}")
    else:
        print("  FAILED")

    # --- Method 3: ScalarFEMBEM FEM-only (for reference) ---
    print("\n--- Method 3: EddyCurrentFEMBEM FEM-only (scalar Hz, Dirichlet) ---")
    r = run_method("ScalarFEM", """
        results = {"method": "ScalarFEM_only", "losses": []}
        for f in freqs:
            try:
                solver = EddyCurrentFEMBEM(mesh, sigma=sigma,
                                            mu_r=mu_r, order=2)
                solver.assemble_fem(f)
                solver.solve(B_ext=[0, 0, 1.0], mode='fem')
                P = solver.compute_loss()
                results["losses"].append({"freq": f, "P": P})
            except Exception as e:
                results["losses"].append({"freq": f, "P": None,
                                          "error": str(e)})
    """)
    if r:
        if 'error' in r:
            print(f"  ERROR: {r['error']}")
        else:
            all_results["ScalarFEM"] = r
            for item in r["losses"]:
                if item.get("P") is not None:
                    print(f"  f={item['freq']/1e3:8.1f} kHz  P={item['P']:.4e} W")
                else:
                    print(f"  f={item['freq']/1e3:8.1f} kHz  FAILED: "
                          f"{item.get('error', 'unknown')}")
    else:
        print("  FAILED")

    # --- Comparison Table ---
    print("\n" + "=" * 70)
    print("CROSS-METHOD COMPARISON TABLE")
    print("=" * 70)

    methods = list(all_results.keys())
    header = f"{'f (kHz)':>10}"
    for m in methods:
        header += f"  {m:>16}"
    header += f"  {'Analytical':>16}"
    print(header)

    for i, f in enumerate(FREQS):
        delta = np.sqrt(2 / (2*np.pi*f * MU_0 * SIGMA * MU_R))
        Rs = np.sqrt(2*np.pi*f * MU_0 * MU_R / (2 * SIGMA))
        H0 = 1.0 / MU_0
        side_area = 2 * (WIDTH + HEIGHT) * DEPTH
        P_ana = 0.5 * Rs * H0**2 * side_area

        row = f"{f/1e3:10.1f}"
        for m in methods:
            losses = all_results[m]["losses"]
            if i < len(losses) and losses[i].get("P") is not None:
                row += f"  {losses[i]['P']:16.4e}"
            else:
                row += f"  {'N/A':>16}"
        row += f"  {P_ana:16.4e}"
        print(row)

    # Ratios vs VectorFEMBEM
    ref_method = "VectorFEMBEM" if "VectorFEMBEM" in all_results else None
    if ref_method is None:
        ref_method = "ShieldSIBC" if "ShieldSIBC" in all_results else None

    if ref_method:
        ref = all_results[ref_method]["losses"]
        print(f"\nRatios (vs {ref_method}):")
        for m in methods:
            if m == ref_method:
                continue
            losses = all_results[m]["losses"]
            for i, f in enumerate(FREQS):
                if (i < len(losses) and i < len(ref)
                        and losses[i].get("P") is not None
                        and ref[i].get("P") is not None
                        and ref[i]["P"] > 0):
                    ratio = losses[i]["P"] / ref[i]["P"]
                    delta = np.sqrt(2 / (2*np.pi*f * MU_0 * SIGMA * MU_R))
                    print(f"  {m:>16} @ {f/1e3:.0f} kHz "
                          f"(delta={delta*1e3:.2f} mm): "
                          f"ratio = {ratio:.4f}")

    # Analytical reference
    print(f"\nAnalytical reference (half-space, side faces only):")
    print(f"  B_ext = 1 T uniform, H0 = B/mu_0 = {1.0/MU_0:.0f} A/m")
    for f in FREQS:
        delta = np.sqrt(2 / (2*np.pi*f * MU_0 * SIGMA * MU_R))
        Rs = np.sqrt(2*np.pi*f * MU_0 * MU_R / (2 * SIGMA))
        H0 = 1.0 / MU_0
        side_area = 2 * (WIDTH + HEIGHT) * DEPTH
        P_ana = 0.5 * Rs * H0**2 * side_area
        print(f"  f={f/1e3:.0f} kHz: delta={delta*1e3:.3f} mm, "
              f"R_s={Rs:.6f} Ohm/sq, P_sides={P_ana:.2e} W")

    # Summary verdict
    print("\n" + "=" * 70)
    print("CONSISTENCY ASSESSMENT")
    print("=" * 70)

    if "VectorFEMBEM" in all_results and "ShieldSIBC" in all_results:
        v = all_results["VectorFEMBEM"]["losses"]
        s = all_results["ShieldSIBC"]["losses"]
        print(f"\n  VectorFEMBEM vs ShieldBEMSIBC (primary comparison):")
        for i, f in enumerate(FREQS):
            if (i < len(v) and i < len(s)
                    and v[i].get("P") is not None
                    and s[i].get("P") is not None):
                pv, ps = v[i]["P"], s[i]["P"]
                delta = np.sqrt(2 / (2*np.pi*f * MU_0 * SIGMA * MU_R))
                diff_pct = 100.0 * abs(pv - ps) / max(abs(pv), abs(ps))
                verdict = "OK" if diff_pct < 20 else "CHECK"
                print(f"    f={f/1e3:6.0f} kHz  "
                      f"Vector={pv:.2e}  Shield={ps:.2e}  "
                      f"diff={diff_pct:.1f}%  "
                      f"delta={delta*1e3:.2f} mm  [{verdict}]")

    if "ScalarFEM" in all_results:
        print(f"\n  ScalarFEM FEM-only (expected ~60-70% of VectorFEMBEM):")
        print(f"  (Scalar Hz captures only z-component; misses Hx,Hy in 3D)")
        if "VectorFEMBEM" in all_results:
            v = all_results["VectorFEMBEM"]["losses"]
            s = all_results["ScalarFEM"]["losses"]
            for i, f in enumerate(FREQS):
                if (i < len(v) and i < len(s)
                        and v[i].get("P") is not None
                        and s[i].get("P") is not None
                        and v[i]["P"] > 0):
                    ratio = s[i]["P"] / v[i]["P"]
                    print(f"    f={f/1e3:6.0f} kHz  "
                          f"Scalar/Vector = {ratio:.3f}")

    print("\n  Method validity notes:")
    print("  - VectorFEMBEM: Valid for all 3D, requires mesh << skin depth")
    print("  - ShieldBEMSIBC: Valid when delta << conductor size (SIBC)")
    print("  - ScalarFEM (Hz): Underestimates for 3D (only z-component)")
    print("  - ScalarSIBC: NOT suitable for 3D blocks (Hz_total -> 0)")


if __name__ == '__main__':
    main()
