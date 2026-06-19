# -*- coding: utf-8 -*-
r"""
demo_budget_dofcost.py  (Track A -- kelvin branch)
==================================================
The COMMITTED, asserted backing for the manuscript error-budget table (Tab.3 /
table:budget) and the DoF-cost figure (fig_dofcost).  Both were live-reproducible
from committed production functions but had NO test asserting the specific values;
this demo pins them (and the BEM open-BC number that had no committed assertion).

On a SHARED interior shell mesh (R_in=0.5, R_out=1.0, maxh=0.4, order 2, dipole),
swap ONLY the Gamma operator -- exact-DtN Robin, Kelvin-DtN Robin, BEM-DtN Schur --
so the interior FEM error CANCELS and the residual is each method's pure open-BC
error.  The DoF-cost sweep refines only the exterior (Kelvin ball) with the interior
fixed, reading the open-BC error vs the exterior DoF increment.

PROVEN HERE (all asserted; every 'ok' gated on a computed number):
  [A] BUDGET: interior FEM error 5.34%; Kelvin open-BC error 0.118% (= 1/45 of it);
      BEM open-BC error 0.172% (= 1/31).  Both closures are >=1 order below the
      interior FEM error -> the open boundary is NOT the bottleneck.
  [B] DOF-COST: even the coarsest Kelvin ball (DeltaDoF~58, ~Gamma scale) gives
      open-BC error 1.2e-3 = 1/45 of the interior error; refining the ball
      (DeltaDoF 58->301) drops it 1.2e-3 -> 7.5e-5 but that is over-resolving a
      non-bottleneck.  The error CONVERGES and stays BELOW the interior FEM error
      at every exterior mesh.

Uses the committed production functions kelvin_vs_exact_open_bc_error (with BEM) and
kelvin_openbc_error_vs_exterior_mesh.  Slow (a dense BEM Schur solve + a 4-mesh
exterior sweep): minutes.
"""
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
from radia_mcp.radia_ngsolve.fem_bem_coupling import (
    kelvin_vs_exact_open_bc_error, kelvin_openbc_error_vs_exterior_mesh)

N_FAIL = 0
def check(name, cond, detail=""):
    global N_FAIL
    if not cond: N_FAIL += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}{('  -- ' + detail) if detail else ''}")

print("=" * 76)
print(" demo_budget_dofcost: error budget (Tab.3) + DoF-cost (fig_dofcost)")
print("=" * 76)

# ===========================================================================
# [A] ERROR BUDGET -- interior FEM vs isolated Kelvin / BEM open-BC error
# ===========================================================================
print("\n[A] Error budget (shell dipole, R_in=0.5/R_out=1.0, maxh=0.4, order 2):")
b = kelvin_vs_exact_open_bc_error(R_inner=0.5, R_outer=1.0, maxh=0.4, order=2, degree=1, include_bem=True)
fem = b["interior_fem_error"]; kel = b["kelvin_openbc_error"]; bem = b["bem_openbc_error"]
print(f"    interior FEM error (open BC exact) = {fem*100:.2f}%")
print(f"    Kelvin open-BC error               = {kel*100:.3f}%   (interior/Kelvin = {fem/kel:.1f})")
print(f"    BEM    open-BC error               = {bem*100:.3f}%   (interior/BEM    = {fem/bem:.1f})")
check("interior FEM error ~ 5.34% (in [4.5%, 6.0%])", 0.045 < fem < 0.060, f"{fem*100:.2f}%")
check("Kelvin open-BC error ~ 0.12% (< 0.2%)", kel < 2e-3, f"{kel*100:.3f}%")
check("BEM open-BC error ~ 0.17% (< 0.25%)", bem < 2.5e-3, f"{bem*100:.3f}%")
check("Kelvin open-BC error is ~1/45 of interior FEM error (ratio > 35)", fem / kel > 35, f"ratio = {fem/kel:.1f}")
check("BEM open-BC error is ~1/31 of interior FEM error (ratio > 28)", fem / bem > 28, f"ratio = {fem/bem:.1f}")

# ===========================================================================
# [B] DoF-COST -- open-BC error vs exterior (Kelvin ball) DoF increment
# ===========================================================================
print("\n[B] DoF-cost: refine only the exterior (Kelvin ball), interior fixed:")
d = kelvin_openbc_error_vs_exterior_mesh(R_inner=0.5, R_outer=1.0, maxh=0.4, order=2, degree=1)
fem2 = d["interior_fem_error"]
for m in d["per_mesh"]:
    print(f"    kmaxh={m['kelvin_maxh']:.2f}  DeltaDoF={m['kelvin_ndof']:4d}  open-BC={m['kelvin_openbc_error']:.2e}  (interior/open-BC = {m['ratio_fem_over_openbc']:.0f})")
coarsest = d["per_mesh"][0]; finest = d["per_mesh"][-1]
check("coarsest Kelvin ball is Gamma-scale: DeltaDoF ~ 58 (in [40,80])",
      40 <= coarsest["kelvin_ndof"] <= 80, f"DeltaDoF = {coarsest['kelvin_ndof']}")
check("coarsest ball already gives open-BC error ~1.2e-3 = 1/45 of interior (ratio > 35)",
      coarsest["kelvin_openbc_error"] < 2e-3 and coarsest["ratio_fem_over_openbc"] > 35,
      f"open-BC = {coarsest['kelvin_openbc_error']:.2e}, ratio = {coarsest['ratio_fem_over_openbc']:.0f}")
check("refining the ball (DeltaDoF -> ~301) drops open-BC error to ~7.5e-5 (< 1.5e-4)",
      finest["kelvin_ndof"] > 250 and finest["kelvin_openbc_error"] < 1.5e-4,
      f"DeltaDoF = {finest['kelvin_ndof']}, open-BC = {finest['kelvin_openbc_error']:.2e}")
check("open-BC error converges under exterior refinement AND stays below interior FEM error",
      d["converges"] and d["always_below_fem"],
      f"converges={d['converges']}, always_below_fem={d['always_below_fem']}")

print("\n" + "=" * 76)
if N_FAIL == 0:
    print(" ALL CHECKS PASSED -- the open boundary is not the bottleneck: Kelvin closes it")
    print(" with a Gamma-scale ball (DeltaDoF~58) at 1/45 of the interior FEM error; BEM 1/31.")
else:
    print(f" {N_FAIL} CHECK(S) FAILED")
print("=" * 76)
assert N_FAIL == 0, f"{N_FAIL} checks failed"
