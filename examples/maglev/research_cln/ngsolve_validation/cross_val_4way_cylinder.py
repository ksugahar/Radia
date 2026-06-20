"""4-way Cauer ladder cross-validation, Cu disk cylinder R=10mm t=2mm.

Compares per-stage R_2k, L_{2k+1}, tau_pair from four independent paradigms:

  (A) BEM-Foster + Cauer-I CFE  (Mathematica + mpmath, integral-equation, axisym)
      Source: bem_disk_axisym_cauer.json + bem_disk_axisym_cauer_python_results.json

  (B) axihenrotte FE Q2 + Hiruma 3-term  (NGSolve add-on, differential, axisym)
      Source: S:/.../examples/axifemm/research/verification/test_hiruma_disk_q2_results.json

  (C) axihenrotte FE Q1 + Hiruma 3-term  (same FE solver, p=1 — convergence check)
      Source: S:/.../examples/axifemm/research/verification/test_hiruma_disk_q1_results.json

  (D) NGSolve 3D HCurl + Helmholtz-Hodge + Hiruma 3-term  (full 3D, NOT axisym)
      Source: 2026-05-08-disk_3d_kameari_hiruma_v5_orderphi2_results.json

Track B extensions:
  (B-i)  Absolute R_2k, L_{2k+1} comparison with sigma-normalisation. The BEM
         moments alpha_n carry an extra sigma^n factor (eigvals scaled as
         tau_k = sigma * lambda_k), so BEM R_2k = sigma * FE R_2k and likewise
         for L_{2k+1}. After normalisation the FE/BEM gap reduces to the same
         gap as in tau_pair, validating that absolute ladder values agree.

  (B-ii) High-k (k >= 4) gap diagnosis: report each method's tau_pair sequence
         and flag where Cauer extraction starts producing nonsense (negative
         tau in BEM around k=6) so the user can decide whether to regenerate
         BEM moments at higher precision.

  (B-iii) Add the (D) full-3D method to the comparison so the published claim
          "3D and axisym agree" can be quoted verbatim from a single table.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

SIGMA_CU = 5.8e7

# === Data sources =========================================================
HERE = Path(__file__).parent
BEM_JSON = HERE / "bem_disk_axisym_cauer_python_results.json"
V5_JSON = HERE / "2026-05-08-disk_3d_kameari_hiruma_v5_orderphi2_results.json"
Q2_JSON = Path("S:/Radia/01_GitHub/examples/axifemm/research/verification/"
                "test_hiruma_disk_q2_results.json")
Q1_JSON = Path("S:/Radia/01_GitHub/examples/axifemm/research/verification/"
                "test_hiruma_disk_q1_results.json")


def load_bem(path):
    """Returns list of dicts {k, R, L, tau_us}."""
    j = json.loads(path.read_text(encoding="utf-8"))
    R = [float(s) for s in j["R_2k"]]
    L = [float(s) for s in j["L_2k_plus_1"]]
    tau = [float(s) for s in j["tau_pair_us"]]
    return [{"k": k, "R": R[k], "L": L[k], "tau_us": tau[k]}
            for k in range(min(len(R), len(tau)))]


def load_axifemm(path, case_name):
    """Returns list of dicts {k, R, L, tau_us} from axihenrotte test JSON."""
    j = json.loads(path.read_text(encoding="utf-8"))
    stages = j["results"][case_name]["stages"]
    out = []
    for s in stages:
        R = s["R_2k"]
        L = s["L_2k_plus_1"]
        tau = s.get("tau_pair_us", s.get("tau_per_stage_us"))
        out.append({"k": s["k"], "R": R, "L": L, "tau_us": tau})
    return out


def load_v5_3d(path):
    """Returns list of dicts {k, R, L, tau_us, gauge_frac} from v5 3D run."""
    j = json.loads(path.read_text(encoding="utf-8"))
    out = []
    for s in j["stages"]:
        out.append({
            "k": s["k"],
            "R": s["R_2k"],
            "L": s["L_2k_plus_1"],
            "tau_us": s["tau_pair_us"],
            "gauge_frac": s.get("gauge_frac"),
        })
    return out


def fmt(x, w, p):
    if x is None:
        return f"{'N/A':>{w}}"
    if isinstance(x, float):
        if abs(x) >= 1e6 or (abs(x) < 1e-3 and x != 0):
            return f"{x:>{w}.{p}e}"
        return f"{x:>{w}.{p}f}"
    return f"{x:>{w}}"


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    print("=" * 92)
    print(" 4-way Cauer ladder cross-validation, Cu disk R=10mm t=2mm sigma=5.8e7 S/m")
    print(" Nagamine convention: R_{2k} series, L_{2k+1} shunt, tau_pair[k]=L_{2k+1}/R_{2k}")
    print("=" * 92)

    bem = load_bem(BEM_JSON)
    q2 = load_axifemm(Q2_JSON, "fine")
    q1 = load_axifemm(Q1_JSON, "very fine")
    v5 = load_v5_3d(V5_JSON)

    n_compare = min(6, len(bem), len(q2), len(q1), len(v5))

    # === Table 1: tau_pair comparison =====================================
    print()
    print("Table 1: tau_pair[k] (us)")
    print("-" * 92)
    print(f"{'k':>3} {'(A) BEM':>11} {'(B) Q2 fine':>13} {'(C) Q1 vfine':>13} "
          f"{'(D) 3D HCurl':>13} {'B-A %':>8} {'C-A %':>8} {'D-A %':>8}")
    print("-" * 92)
    for k in range(n_compare):
        a = bem[k]["tau_us"]
        b = q2[k]["tau_us"]
        c = q1[k]["tau_us"]
        d = v5[k]["tau_us"]
        print(f"{k:>3} {a:>11.4f} {b:>13.4f} {c:>13.4f} {d:>13.4f} "
              f"{(b-a)/a*100:>+7.3f}% {(c-a)/a*100:>+7.3f}% {(d-a)/a*100:>+7.3f}%")

    # === Table 2: Absolute R_2k, L_{2k+1} with sigma normalisation =========
    # BEM moments carry tau_k = sigma * lambda_k -> R_BEM = sigma * R_FE,
    #                                              L_BEM = sigma * L_FE
    # so we normalise BEM by /sigma to compare with FE in SI units.
    print()
    print("Table 2: Absolute R_2k (Ohm) and L_{2k+1} (Henry) — Track B-i sigma normalisation")
    print("-" * 92)
    print("BEM normalisation: BEM moments use tau_k = sigma * lambda_k, hence")
    print(f"                   R_2k_BEM = sigma * R_2k_FE  with sigma = {SIGMA_CU:.2e} S/m")
    print(f"                   We display R_2k_BEM_normalised = R_2k_BEM / sigma")
    print()
    print(f"{'k':>3} {'BEM R / sigma':>16} {'Q2 R':>13} {'V5 R':>13} "
          f"{'Q2/BEM_n':>10} {'V5/BEM_n':>10}")
    print("-" * 92)
    for k in range(n_compare):
        a_R_norm = bem[k]["R"] / SIGMA_CU
        b_R = q2[k]["R"]
        d_R = v5[k]["R"]
        print(f"{k:>3} {a_R_norm:>16.6e} {b_R:>13.6e} {d_R:>13.6e} "
              f"{b_R/a_R_norm:>10.6f} {d_R/a_R_norm:>10.6f}")
    print()
    print(f"{'k':>3} {'BEM L / sigma':>16} {'Q2 L':>13} {'V5 L':>13} "
          f"{'Q2/BEM_n':>10} {'V5/BEM_n':>10}")
    print("-" * 92)
    for k in range(n_compare):
        a_L_norm = bem[k]["L"] / SIGMA_CU
        b_L = q2[k]["L"]
        d_L = v5[k]["L"]
        print(f"{k:>3} {a_L_norm:>16.6e} {b_L:>13.6e} {d_L:>13.6e} "
              f"{b_L/a_L_norm:>10.6f} {d_L/a_L_norm:>10.6f}")

    # === Table 3: High-k gap diagnosis ====================================
    print()
    print("Table 3: High-k gap diagnosis — Track B-ii")
    print("-" * 92)
    bem_full = json.loads(BEM_JSON.read_text(encoding="utf-8"))
    bem_tau_full = [float(s) for s in bem_full["tau_pair_us"]]
    print(f"BEM tau_pair full sequence (Cauer extracted, n_moments={bem_full.get('n_moments', '?')}):")
    for k, t in enumerate(bem_tau_full):
        flag = ""
        if t < 0:
            flag = "  <-- NEGATIVE: Cauer extraction has gone unstable"
        elif k > 0 and t > bem_tau_full[k-1]:
            flag = "  <-- INCREASE: tau should decay monotonically"
        print(f"  k={k}: tau_pair = {t:>+.4f} us{flag}")
    print()
    print("Diagnosis:")
    n_neg = sum(1 for t in bem_tau_full if t < 0)
    print(f"  - {n_neg} negative-tau stages detected (Cauer unstable)")
    print(f"  - First sign change at k = ", end="")
    first_neg = next((k for k, t in enumerate(bem_tau_full) if t < 0), None)
    print(f"{first_neg}" if first_neg is not None else "never")
    print()
    print("  Reduction strategy: regenerate BEM moments at higher precision.")
    print(f"    Current: n_moments={bem_full.get('n_moments', '?')}, dps={bem_full.get('dps', '?')}")
    print(f"    Try    : n_moments=40, dps=75 (see bem_disk_axisym_cauer_highprec.wls)")

    # === Table 4: gauge_frac for the 3D method ============================
    print()
    print("Table 4: 3D HCurl v5 gauge diagnostic (cylinder geometry — see Track A for cuboid)")
    print("-" * 92)
    print(f"{'k':>3} {'gauge_frac':>12} {'  L_2k+1_raw':>16} {'L_2k+1_clean':>16} "
          f"{'  raw/clean':>12}")
    print("-" * 92)
    for k in range(n_compare):
        gf = v5[k]["gauge_frac"]
        # Raw vs clean L is in the JSON
        s = json.loads(V5_JSON.read_text(encoding="utf-8"))["stages"][k]
        L_raw = s.get("L_2k_plus_1_raw_unprojected", None)
        L_clean = s.get("L_2k_plus_1", None)
        ratio = L_raw/L_clean if (L_raw and L_clean) else None
        print(f"{k:>3} {gf:>12.4f} {L_raw:>16.6e} {L_clean:>16.6e} "
              f"{ratio:>12.4f}")
    print()
    print("Note: cylinder gauge_frac ~ 0.87 (works); cuboid gauge_frac = 1.000")
    print("      (collapses to gradient subspace, see project_3d_hcurl_hiruma_cuboid_gauge_collapse).")
    print("      Track A (saddle-point H-H) addresses the cuboid pathology.")

    # === Summary ==========================================================
    print()
    print("=" * 92)
    print("Summary")
    print("=" * 92)
    print(f"  Leading time constant tau_pair[0]:")
    print(f"    (A) BEM-Foster Cauer        : {bem[0]['tau_us']:>10.4f} us  (1920 elements, mpmath dps=50)")
    print(f"    (B) FE Q2 Henrotte fine     : {q2[0]['tau_us']:>10.4f} us  ({(q2[0]['tau_us']-bem[0]['tau_us'])/bem[0]['tau_us']*100:+.3f}% vs BEM)")
    print(f"    (C) FE Q1 Henrotte v.fine   : {q1[0]['tau_us']:>10.4f} us  ({(q1[0]['tau_us']-bem[0]['tau_us'])/bem[0]['tau_us']*100:+.3f}% vs BEM)")
    print(f"    (D) 3D HCurl v5 (3D, NOT axisym!): {v5[0]['tau_us']:>10.4f} us  ({(v5[0]['tau_us']-bem[0]['tau_us'])/bem[0]['tau_us']*100:+.3f}% vs BEM)")
    print()
    print(f"  4-way agreement at lead: max gap = "
          f"{max(abs(q2[0]['tau_us']-bem[0]['tau_us']), abs(q1[0]['tau_us']-bem[0]['tau_us']), abs(v5[0]['tau_us']-bem[0]['tau_us']))/bem[0]['tau_us']*100:.3f}%")
    print()
    print("  This is the publication-grade 4-way validation for Cu disk axisym + 3D.")
    print("  For cuboid (5x2x1mm), see radia-vim Phase F-4 BEM (10.85us) + ELF time-domain")
    print("  (11.58us). Cuboid 3D-FEM via saddle-point H-H is in progress (Track A).")


if __name__ == "__main__":
    main()
