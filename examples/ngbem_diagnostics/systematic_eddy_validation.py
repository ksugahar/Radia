"""
Systematic eddy current solver validation.

Sweeps geometry (plate vs cube) x frequency to map out
where each method is valid and WHY they disagree.

Physics:
  For a block in uniform Bz, tangential H on each face is:
    - Side faces (x,y normals): Hz is tangential -> captured by scalar Hz
    - Top/bottom (z normal):    Hz is NORMAL -> Hx,Hy are tangential
  Scalar Hz formulation ONLY captures Hz-driven eddy currents,
  missing the Hx,Hy contributions entirely.

  SIBC assumes delta << conductor dimension (half-space approximation).
  When delta ~ thickness, SIBC breaks down.

  FEM requires mesh << delta. When delta << maxh, FEM saturates.
"""
import sys
import os
import subprocess
import json
import textwrap
import numpy as np

MU_0 = 4e-7 * np.pi

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(SCRIPT_DIR, '../../src/radia')

SIGMA = 3.7e7  # Al
MU_R = 1.0
FREQS = [1e3, 3e3, 10e3, 30e3, 100e3]


def skin_depth(freq):
    return np.sqrt(2 / (2 * np.pi * freq * MU_0 * SIGMA * MU_R))


def Rs(freq):
    return np.sqrt(2 * np.pi * freq * MU_0 * MU_R / (2 * SIGMA))


def analytical_halfspace(freq, W, H, D):
    """Half-space loss: side faces only (Hz tangential)."""
    H0 = 1.0 / MU_0
    side_area = 2 * (W + H) * D
    return 0.5 * Rs(freq) * H0**2 * side_area


def analytical_all_faces(freq, W, H, D):
    """Half-space loss: ALL 6 faces (upper bound, all faces contribute)."""
    H0 = 1.0 / MU_0
    total_area = 2 * (W * H + W * D + H * D)
    return 0.5 * Rs(freq) * H0**2 * total_area


def run_method(label, code, geometry, maxh):
    """Run solver in subprocess."""
    W, H, D = geometry
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

    mesh = create_conductor_mesh({W}, {H}, {D}, maxh={maxh})
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
        for line in result.stdout.strip().split('\n'):
            line = line.strip()
            if line.startswith('JSONRESULT:'):
                return json.loads(line[len('JSONRESULT:'):])
        if result.stderr:
            print(f"  [{label}] stderr: {result.stderr[:300]}")
        return None
    except subprocess.TimeoutExpired:
        print(f"  [{label}] Timeout")
        return None
    except Exception as e:
        print(f"  [{label}] Exception: {e}")
        return None


def run_vector_fembem(geometry, maxh):
    """VectorEddyCurrentFEMBEM: full 3D vector A."""
    return run_method("VectorFEMBEM", """
solver = VectorEddyCurrentFEMBEM(mesh, sigma, mu_r=mu_r, order=1)
solver.assemble(kappa=0.01, intorder=4)

results = {"losses": []}
for f in freqs:
    try:
        res = solver.solve(f, B_ext=[0, 0, 1.0])
        P = solver.compute_loss()
        results["losses"].append({"freq": f, "P": P, "ndof": res["n_total"]})
    except Exception as e:
        results["losses"].append({"freq": f, "P": None, "error": str(e)})
""", geometry, maxh)


def run_shield_sibc(geometry, maxh):
    """ShieldBEMSIBC: loop-basis BEM + SIBC."""
    return run_method("ShieldSIBC", """
solver = ShieldBEMSIBC(mesh, sigma, mu_r=mu_r)
solver.assemble(intorder=4)

results = {"losses": [], "n_loops": solver._loop.n_loops}
for f in freqs:
    try:
        def A_inc_uniform(points):
            A = np.zeros_like(points)
            A[:, 0] = -0.5 * points[:, 1]
            A[:, 1] = 0.5 * points[:, 0]
            return A
        solver.solve(f, A_inc_func=A_inc_uniform)
        P = solver.compute_loss()
        results["losses"].append({"freq": f, "P": P})
    except Exception as e:
        results["losses"].append({"freq": f, "P": None, "error": str(e)})
""", geometry, maxh)


def run_scalar_fem(geometry, maxh):
    """EddyCurrentFEMBEM FEM-only: scalar Hz."""
    return run_method("ScalarFEM", """
results = {"losses": []}
for f in freqs:
    try:
        solver = EddyCurrentFEMBEM(mesh, sigma=sigma, mu_r=mu_r, order=2)
        solver.assemble_fem(f)
        solver.solve(B_ext=[0, 0, 1.0], mode='fem')
        P = solver.compute_loss()
        results["losses"].append({"freq": f, "P": P})
    except Exception as e:
        results["losses"].append({"freq": f, "P": None, "error": str(e)})
""", geometry, maxh)


def run_scalar_sibc(geometry, maxh):
    """EddyCurrentBEMSIBC: scalar Hz + SIBC."""
    return run_method("ScalarSIBC", """
solver = EddyCurrentBEMSIBC(mesh, sigma=sigma, mu_r=mu_r, order=2)
solver.assemble(intorder=8)

results = {"losses": []}
for f in freqs:
    try:
        solver.solve(f, B_ext=[0, 0, 1.0])
        P = solver.compute_loss()
        Hz_scat = solver._Hz_scat_surf
        Hz_inc_surf = np.full(solver._n_surf, solver._Hz_inc, dtype=complex)
        Hz_total = Hz_scat + Hz_inc_surf
        ratio = float(np.mean(np.abs(Hz_total)) / abs(solver._Hz_inc))
        results["losses"].append({"freq": f, "P": P, "Hz_ratio": ratio})
    except Exception as e:
        results["losses"].append({"freq": f, "P": None, "error": str(e)})
""", geometry, maxh)


def extract_P(result, freq_idx):
    """Extract loss from result dict."""
    if result is None:
        return None
    if 'error' in result:
        return None
    losses = result.get("losses", [])
    if freq_idx < len(losses) and losses[freq_idx].get("P") is not None:
        return losses[freq_idx]["P"]
    return None


def main():
    print("=" * 78)
    print("SYSTEMATIC EDDY CURRENT VALIDATION")
    print("Geometry x Frequency sweep: method validity mapping")
    print("=" * 78)

    # Two geometries: same footprint, different thickness
    # Key variable: side_face_fraction = A_sides / A_total
    #   Hz is tangential on side faces only (for Bz excitation)
    #   Scalar Hz only captures loss from Hz -> correlates with side fraction
    geometries = {
        "Thin 20x20x5mm": {
            "dims": (0.02, 0.02, 0.005),  # thin block
            "maxh": 0.004,
            "note": "Side fraction=0.33, Hz tangential on 1/3 of area",
        },
        "Block 20x20x10mm": {
            "dims": (0.02, 0.02, 0.01),  # standard test case
            "maxh": 0.005,
            "note": "Side fraction=0.50, Hz tangential on 1/2 of area",
        },
    }

    methods = {
        "VectorFEMBEM": run_vector_fembem,
        "ShieldSIBC": run_shield_sibc,
        "ScalarFEM": run_scalar_fem,
        "ScalarSIBC": run_scalar_sibc,
    }

    all_data = {}

    for geo_name, geo_info in geometries.items():
        W, H, D = geo_info["dims"]
        maxh = geo_info["maxh"]
        min_dim = min(W, H, D)

        print(f"\n{'='*78}")
        print(f"Geometry: {geo_name}")
        print(f"  Dimensions: {W*1e3:.1f} x {H*1e3:.1f} x {D*1e3:.1f} mm")
        print(f"  Volume: {W*H*D*1e9:.0f} mm^3")
        print(f"  Min dimension: {min_dim*1e3:.1f} mm")
        print(f"  maxh: {maxh*1e3:.1f} mm")
        print(f"  Note: {geo_info['note']}")

        # Face area analysis
        A_top_bot = 2 * W * H
        A_sides = 2 * (W + H) * D
        A_total = A_top_bot + A_sides
        print(f"  Face areas: top+bot={A_top_bot*1e6:.0f} mm^2, "
              f"sides={A_sides*1e6:.0f} mm^2, "
              f"ratio(top_bot/total)={A_top_bot/A_total:.2f}")

        # Frequency info
        for f in FREQS:
            d = skin_depth(f)
            print(f"    f={f/1e3:6.0f} kHz: delta={d*1e3:.3f} mm, "
                  f"delta/min_dim={d/min_dim:.3f}, "
                  f"maxh/delta={maxh/d:.1f}")

        geo_data = {}

        # Run each method
        for method_name, method_func in methods.items():
            print(f"\n  Running {method_name}...", end=" ", flush=True)
            result = method_func(geo_info["dims"], maxh)
            geo_data[method_name] = result

            if result is None:
                print("FAILED")
            elif 'error' in result:
                print(f"ERROR: {result['error'][:100]}")
            else:
                losses = result.get("losses", [])
                ok_count = sum(1 for l in losses if l.get("P") is not None)
                print(f"OK ({ok_count}/{len(FREQS)} freqs)")

        all_data[geo_name] = geo_data

        # Comparison table for this geometry
        print(f"\n  {'f(kHz)':>8} {'delta':>8} {'d/Dmin':>7} "
              f"{'Vector':>12} {'Shield':>12} {'ScalFEM':>12} "
              f"{'ScalSIBC':>12} {'Ana(side)':>12} {'Ana(all)':>12}")
        print(f"  {'-'*8} {'-'*8} {'-'*7} "
              f"{'-'*12} {'-'*12} {'-'*12} "
              f"{'-'*12} {'-'*12} {'-'*12}")

        for i, f in enumerate(FREQS):
            d = skin_depth(f)
            P_ana_side = analytical_halfspace(f, W, H, D)
            P_ana_all = analytical_all_faces(f, W, H, D)

            row = f"  {f/1e3:8.0f} {d*1e3:7.3f}m {d/min_dim:7.3f}"

            for mn in ["VectorFEMBEM", "ShieldSIBC", "ScalarFEM", "ScalarSIBC"]:
                P = extract_P(geo_data.get(mn), i)
                if P is not None:
                    row += f" {P:12.1f}"
                else:
                    row += f" {'N/A':>12}"

            row += f" {P_ana_side:12.1f} {P_ana_all:12.1f}"
            print(row)

        # Ratio analysis: each method vs VectorFEMBEM
        ref = geo_data.get("VectorFEMBEM")
        if ref and 'error' not in ref:
            print(f"\n  Ratios vs VectorFEMBEM:")
            print(f"  {'f(kHz)':>8} {'d/Dmin':>7} "
                  f"{'Shield':>8} {'ScalFEM':>8} {'ScalSIBC':>8} "
                  f"{'Validity':>30}")
            for i, f in enumerate(FREQS):
                d = skin_depth(f)
                P_ref = extract_P(ref, i)
                if P_ref is None or P_ref == 0:
                    continue

                row = f"  {f/1e3:8.0f} {d/min_dim:7.3f}"
                notes = []

                for mn in ["ShieldSIBC", "ScalarFEM", "ScalarSIBC"]:
                    P = extract_P(geo_data.get(mn), i)
                    if P is not None:
                        ratio = P / P_ref
                        row += f" {ratio:8.3f}"
                    else:
                        row += f" {'N/A':>8}"

                # Validity assessment
                if d / min_dim > 0.5:
                    notes.append("SIBC invalid(d~Dmin)")
                if maxh / d > 3:
                    notes.append("FEM coarse(maxh>>d)")
                if d / min_dim < 0.1:
                    notes.append("SIBC good")
                if maxh / d < 1:
                    notes.append("FEM good")
                row += f" {', '.join(notes):>30}" if notes else ""
                print(row)

    # === GLOBAL ANALYSIS ===
    print(f"\n\n{'='*78}")
    print("SYSTEMATIC ANALYSIS")
    print("=" * 78)

    print("""
WHY DO METHODS DISAGREE?

1. Scalar Hz (ScalarFEM, ScalarSIBC) -- fundamental 3D limitation
   ---------------------------------------------------------------
   For block in Bz, tangential H field on each face:
     Side faces (x,y normals): Hz IS tangential -> captured
     Top/bottom (z normal):    Hz is NORMAL -> Hx,Hy tangential, MISSED

   ScalarFEM loss formula: P = 1/(2*sigma) * integral(|dHz/dx|^2+|dHz/dy|^2)
   This only counts J induced by Hz component.
   Missing: Jz from dHy/dx - dHx/dy

   Expected ratio (ScalarFEM / VectorFEMBEM):
     - Plate (top/bot >> sides): ~0.0 to 0.3  (Hz mostly normal on big faces)
     - Cube (balanced):          ~0.5 to 0.7  (partial coverage)
     - Long rod (sides >> top):  ~0.7 to 0.9  (Hz tangential on big faces)

2. ScalarSIBC -- scalar Hz + SIBC compound error
   -----------------------------------------------
   Has BOTH the scalar Hz limitation AND the SIBC approximation.
   Additionally: Hz_total -> 0 as frequency increases (conductor shields Hz)
   But the REAL tangential field (Hx,Hy on top/bottom) does NOT go to zero.
   Result: loss DECREASES with frequency (physically wrong).

3. ShieldBEMSIBC -- SIBC validity range
   ------------------------------------
   Uses full vector surface currents (correct 3D).
   SIBC valid when delta << min(W,H,D):
     delta/Dmin < 0.1: excellent agreement with VectorFEMBEM
     delta/Dmin ~ 0.3: moderate agreement
     delta/Dmin > 0.5: SIBC breaks down (over- or under-estimates)

4. VectorFEMBEM -- mesh resolution limit
   -------------------------------------
   Full 3D correct formulation but needs mesh << delta.
   When maxh/delta > 3: loss saturates (mesh cannot resolve skin layer)
   Solution: finer mesh or use ShieldSIBC for high frequency.

VALIDITY MAP:
                    Low freq              Mid freq             High freq
                  (delta>>Dmin)        (delta~Dmin)          (delta<<Dmin)
  VectorFEMBEM:  [VALID if mesh OK] [VALID if mesh OK]  [INVALID: mesh too coarse]
  ShieldSIBC:    [INVALID: SIBC]    [MARGINAL]          [VALID: SIBC correct]
  ScalarFEM:     [PARTIAL: Hz only] [PARTIAL: Hz only]  [PARTIAL: Hz only]
  ScalarSIBC:    [WRONG: compound]  [WRONG: compound]   [WRONG: Hz->0]

RECOMMENDATION:
  - Low-medium freq (delta > maxh):    VectorEddyCurrentFEMBEM
  - Medium-high freq (delta << Dmin):  ShieldBEMSIBC
  - Crossover region:                  Both agree -> use either
  - NEVER use ScalarSIBC for 3D blocks (fundamentally broken)
  - ScalarFEM: quick estimate only (~60-70% of true loss)
""")

    # Quantitative summary
    print("QUANTITATIVE SUMMARY BY GEOMETRY:")
    for geo_name, geo_info in geometries.items():
        W, H, D = geo_info["dims"]
        min_dim = min(W, H, D)
        maxh = geo_info["maxh"]
        geo_data = all_data.get(geo_name, {})

        print(f"\n  {geo_name} (min_dim={min_dim*1e3:.1f}mm, maxh={maxh*1e3:.1f}mm):")

        # Find the "sweet spot" frequency where both methods are valid
        ref = geo_data.get("VectorFEMBEM")
        shield = geo_data.get("ShieldSIBC")
        if ref and shield and 'error' not in ref and 'error' not in shield:
            best_diff = float('inf')
            best_freq = None
            for i, f in enumerate(FREQS):
                d = skin_depth(f)
                P_v = extract_P(ref, i)
                P_s = extract_P(shield, i)
                if P_v and P_s and P_v > 0:
                    diff = abs(P_v - P_s) / P_v
                    if diff < best_diff:
                        best_diff = diff
                        best_freq = f
            if best_freq:
                d = skin_depth(best_freq)
                print(f"    Best agreement at f={best_freq/1e3:.0f} kHz "
                      f"(delta={d*1e3:.2f}mm, delta/Dmin={d/min_dim:.2f}): "
                      f"diff={best_diff*100:.1f}%")

        # ScalarFEM ratio
        scalar = geo_data.get("ScalarFEM")
        if ref and scalar and 'error' not in ref and 'error' not in scalar:
            ratios = []
            for i in range(len(FREQS)):
                P_v = extract_P(ref, i)
                P_s = extract_P(scalar, i)
                if P_v and P_s and P_v > 0:
                    ratios.append(P_s / P_v)
            if ratios:
                print(f"    ScalarFEM/VectorFEMBEM ratio: "
                      f"{min(ratios):.2f} - {max(ratios):.2f} "
                      f"(mean={np.mean(ratios):.2f})")
                A_top_bot = 2 * W * H
                A_total = A_top_bot + 2 * (W + H) * D
                face_ratio = 1 - A_top_bot / A_total
                print(f"    Side face fraction: {face_ratio:.2f} "
                      f"(Hz tangential on sides only)")


if __name__ == '__main__':
    main()
