"""
hdiv_vim_physics.py -- Phase 5 physics validation: does the hand-built HDiv operator reproduce the
PHYSICAL demag factors (the basis-INVARIANT quantity), resolving the "3x raw-spectrum / 23-vs-8
iters" normalization question?

The raw N spectrum is scale/basis-dependent (mine [-0.095,5.56] vs NGSolve [-0.357,1.515]) -> NOT a
meaningful comparison.  The PHYSICS is the generalized eigenproblem (N, M_mass): its eigenvalues are
the DEMAG FACTORS, bounded in [0,1] and basis-invariant.  If my hand-built (N, M_mass) gives the
SAME demag factors as NGSolve's (N, M_mass), then my +-1/V,1/area normalization is a VALID CHOICE
(the 3x raw difference is pure basis/scale).  If they differ / leave [0,1], there's a missing factor.

Compares: hand-built generalized eig (design-spec N + RT0 mass) vs NGSolve generalized eig (prototype
build() N + mass_hdiv).  Both use the same centroid-monopole G + cube/subdivided self, so the demag
factors should match.  Decisive resolution of the normalization question.
"""
import json, os
import numpy as np
import scipy.linalg as sla

HERE = os.path.dirname(os.path.abspath(__file__))
import sys; sys.path.insert(0, HERE)
from hdiv_vim_structured import (build_structured_rt0, charge_map, coulomb_gram, assemble_mass)

def gen_demag_factors(N, M):
    """generalized eig (N, M); M SPD -> real eigenvalues = demag factors (loops at 0)."""
    w = sla.eigh(N, M, eigvals_only=True)
    return np.sort(w.real)

print("=== Phase 5 physics: generalized eig (N, M_mass) = demag factors (basis-invariant) ===")

# --- hand-built (== C++ rad_hdiv_vim), regular 3x3x3 ---
faces, n_cell, cell_c, cell_V = build_structured_rt0(3, 3, 3, h=1.0)
B, _, _ = charge_map(faces, n_cell, cell_V)
G = coulomb_gram(faces, n_cell, cell_c, cell_V)
N_h = B.T @ G @ B
M_h = assemble_mass(faces, n_cell, h=1.0)
spd = np.all(np.linalg.eigvalsh(M_h) > 0)
df_h = gen_demag_factors(N_h, M_h)
print(f"\n[hand-built] M_mass SPD={spd}; raw N spectrum [{np.linalg.eigvalsh(N_h).min():.3e}, "
      f"{np.linalg.eigvalsh(N_h).max():.3e}]")
print(f"[hand-built] DEMAG FACTORS (N,M_mass): [{df_h[0]:.4e}, {df_h[-1]:.4e}], "
      f"top-3 {np.round(df_h[-3:],4)}, n in (0,1]={int(np.sum((df_h>1e-9)&(df_h<=1.0+1e-6)))}, "
      f"n>1+eps={int(np.sum(df_h>1.0+1e-3))}")

# --- NGSolve prototype, regular 3x3x3 ---
try:
    from hdiv_demag_quad_self import build as proto_build
    d = proto_build(0.0, 6, 8, n=3)
    N_ng, S = d["N"], d["Sstar"]
    # NGSolve mass_hdiv is available on the build dict
    M_ng = d["mass_hdiv"]
    df_ng = gen_demag_factors(N_ng, M_ng)
    print(f"\n[NGSolve]    DEMAG FACTORS (N,M_mass): [{df_ng[0]:.4e}, {df_ng[-1]:.4e}], "
          f"top-3 {np.round(df_ng[-3:],4)}")
    # compare the sorted POSITIVE demag factors (drop the n_loop near-zero + any tiny negatives)
    pos_h = df_h[df_h > 1e-6]; pos_ng = df_ng[df_ng > 1e-6]
    k = min(len(pos_h), len(pos_ng))
    match = np.max(np.abs(pos_h[-k:] - pos_ng[-k:])) if k else float("nan")
    print(f"\n  positive demag-factor match (top {k}): max|d(factor)| = {match:.3e}")
    verdict = (df_h[-1] < 1.05 and match < 1e-2)
    print(f"  => {'PASS: hand-built reproduces NGSolve demag factors -> +-1/V,1/area is a VALID '
          'normalization (raw 3x diff was pure basis/scale)' if verdict else 'MISMATCH: investigate '
          '(possible missing geometric factor)'}")
    out = {"handbuilt_demag_range": [float(df_h[0]), float(df_h[-1])],
           "ngsolve_demag_range": [float(df_ng[0]), float(df_ng[-1])],
           "match_max_abs": float(match), "verdict_valid_normalization": bool(verdict)}
except Exception as e:
    print(f"\n[NGSolve comparison skipped: {e}]")
    out = {"handbuilt_demag_range": [float(df_h[0]), float(df_h[-1])]}

with open(os.path.join(HERE, "hdiv_vim_physics.json"), "w") as f:
    json.dump(out, f, indent=2)
print("\nsaved", os.path.join(HERE, "hdiv_vim_physics.json"))
