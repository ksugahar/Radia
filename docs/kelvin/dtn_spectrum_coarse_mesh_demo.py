"""DtN spectrum vs mesh -- why open boundaries are accurate on COARSE meshes.

A spectral (DtN-matrix) explanation of the coarse-mesh accuracy that Kameari
demonstrated empirically (by mesh refinement) for the Kelvin transformation.
This does NOT reimplement Kelvin (the lab's radia.kelvin_* toolkit already does
that).  It MEASURES the effective exterior DtN two ways and
shows both give the spectral fact that explains the coarse-mesh accuracy:
  Part A  the BEM matrix  Lambda_h = V^-1(-1/2 M + K)  (boundary only), and
  Part B  the KELVIN closure realised by volume FEM on the inverted ball.
Both are the SAME continuous operator Lambda_ext; the two discretisations expose
the coarse-mesh accuracy through complementary mechanisms (Part B note below).

Idea
----
Truncating the unbounded exterior at a surface Gamma and imposing the exact
transparent condition there means imposing the exterior Dirichlet-to-Neumann
(Steklov-Poincare) operator  d u/d n = Lambda_ext u  on Gamma.  Every open-BC
closure (Kelvin / BEM / PML / infinite elements) is a discretization of this one
operator.  On a sphere of radius R it is diagonalised by the spherical harmonics
with the mesh-independent eigenvalue ladder

    Lambda_ext Y_n = -(n+1)/R Y_n     (degree n, multiplicity 2n+1)

We assemble the DISCRETE DtN matrix  Lambda_h = V^-1(-1/2 M + K)  (ngsolve.bem)
on a sphere at a few mesh sizes, take its eigenvalues, and match them to the
analytic ladder.  Three facts come out (and are printed below):

  1. the LOW modes (n<=2) are already accurate on the COARSEST mesh
     (dipole ~0.07 %, quadrupole ~0.2 %) -- the coarse-mesh accuracy itself,
     read straight off the matrix, before solving any field problem;
  2. at EVERY mesh the per-degree error grows with degree n (the spectral
     signature: accurate at the bottom, degraded toward the resolution limit);
  3. refinement lowers every mode steeply (~O(h^4) order-1 superconvergence) and
     WIDENS the accurate band (n<=2 -> n<=4) -- it buys higher modes, not
     materially better low ones.

Because a compact source's exterior field is dominated by the low multipoles,
the coarse mesh already resolves everything that matters -- which matches
Kameari's result, stated as a property of Lambda_h instead of a refinement curve.

Part B (Kelvin, measured): the Kelvin inversion maps the exterior mode r^-(n+1)Y_n
to the solid harmonic rho^n Y_n on the ball -- a degree-n polynomial -- so order-p
FEM reproduces mode n's effective DtN EXACTLY iff p >= n. The dominant dipole
inverts to a LINEAR field, so even order-1 FEM nails it on a coarse mesh (residual
is geometry-limited and converges fast); the quadrupole needs order >= 2. So the
Kelvin closure is coarse-mesh accurate for exactly the low modes a compact source
radiates -- the SAME conclusion as the BEM side, via a complementary mechanism
(boundary smoothness for BEM; polynomial order-threshold for Kelvin).

Part C (accuracy vs EXTERIOR mesh size): in the Kelvin transform the unbounded
exterior IS the Kelvin ball, so the "exterior region mesh size" is kelvin_maxh.
On a shared shell INTERIOR mesh (fixed), refine ONLY the exterior and swap the Γ
operator (exact-DtN Robin vs Kelvin-DtN Robin); the interior FEM error CANCELS, so
the residual is the pure exterior (open-BC) error. It CONVERGES with exterior
refinement (~x4/level) but stays far below the FIXED interior FEM error at every
level (45-700x). This is Kameari's exterior-refinement accuracy check, ISOLATED
from the interior FEM error -- a coarse exterior mesh already suffices.

Companion: bem_integral.{exterior_dtn_spectrum, dtn_spectrum_vs_mesh} (BEM side),
           fem_bem_coupling.{kelvin_dtn_eigenvalue,                  (Kelvin DtN),
                             kelvin_vs_exact_open_bc_error,          (isolated error),
                             kelvin_openbc_error_vs_exterior_mesh}   (vs exterior mesh)
Knowledge:  MCP tool  dtn_coarse_mesh  (overview / numerics / api / applications)
Gate:       validation_test/radia_mcp/test_dtn_spectrum_coarse.py

NB the dense BEM densification (Part A) is O(ndof^2); the default 2-mesh sweep
(ndof 336 + 564) takes a few minutes. The Kelvin volume solves (Part B) are fast.
A sphere surface mesh FLOORS at ndof=336, so every maxh >= 0.5 gives the same
coarsest mesh; use maxh <= 0.5 for distinct levels.
Run:  python docs/kelvin/dtn_spectrum_coarse_mesh_demo.py
"""
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
_RADIA_MCP_SRC = _REPO_ROOT / "packages" / "radia-mcp" / "src"
if str(_RADIA_MCP_SRC) not in sys.path:
    sys.path.insert(0, str(_RADIA_MCP_SRC))

from radia_mcp.radia_ngsolve.bem_integral import exterior_dtn_spectrum  # noqa: E402
from radia_mcp.radia_ngsolve.fem_bem_coupling import (  # noqa: E402
    kelvin_dtn_eigenvalue, kelvin_openbc_error_vs_exterior_mesh,
)

R = 1.0
ORDER = 1
INTORDER = 10
NMAX = 4
MAXH_SWEEP = (0.5, 0.4)        # distinct levels (ndof 336, 564); add 0.3 for ndof 924
BAND_TOL = 0.005               # "accurate" = eigenvalue within 0.5 %
LOW_NMAX = 2                   # monopole..quadrupole = the engineering-dominant block


def _accurate_band(modes, tol):
    """Highest degree n whose eigenvalue is within tol (stops at first failure)."""
    band = -1
    for m in modes:
        if m["rel_err"] < tol:
            band = m["n"]
        else:
            break
    return band


def main():
    print(__doc__.splitlines()[0])
    print(f"sphere R={R}, BEM SurfaceL2 order={ORDER}, intorder={INTORDER}")
    print(f"analytic ladder: lambda_n = -(n+1)/R, multiplicity 2n+1\n")

    per_mesh = []
    for maxh in MAXH_SWEEP:
        print(f"... assembling DtN matrix at maxh={maxh} (dense BEM, please wait)")
        res = exterior_dtn_spectrum(R=R, maxh=maxh, order=ORDER,
                                    intorder=INTORDER, nmax=NMAX)
        per_mesh.append(res)
        band = _accurate_band(res["modes"], BAND_TOL)
        print(f"--- maxh={maxh}  ndof={res['ndof']}  "
              f"accurate band (rel_err<{BAND_TOL:.1%}): n<={band} ---")
        print("   n  name        lam_exact  lam_mean   rel_err")
        names = {0: "monopole", 1: "dipole", 2: "quadrupole", 3: "octupole", 4: "n=4"}
        for m in res["modes"]:
            print(f"   {m['n']}  {names.get(m['n'],''):10s} "
                  f"{m['lambda_exact']:+8.4f}  {m['lambda_mean']:+8.4f}  "
                  f"{m['rel_err']:.3e}")
        print()

    # --- the three spectral facts, as numbers ---
    coarse = per_mesh[0]
    coarse_low = max(m["rel_err"] for m in coarse["modes"] if m["n"] <= LOW_NMAX)
    monotonic = all(
        all(res["modes"][i]["rel_err"] < res["modes"][i + 1]["rel_err"]
            for i in range(len(res["modes"]) - 1))
        for res in per_mesh)
    bands = [_accurate_band(res["modes"], BAND_TOL) for res in per_mesh]
    dipole0 = next(m for m in coarse["modes"] if m["n"] == 1)

    # per-degree convergence ratio coarse -> finest (if a sweep was run)
    if len(per_mesh) >= 2:
        fine = per_mesh[-1]
        ratios = []
        for mc in coarse["modes"]:
            mf = next((x for x in fine["modes"] if x["n"] == mc["n"]), None)
            if mf and mf["rel_err"] > 0:
                ratios.append(mc["rel_err"] / mf["rel_err"])
        drop = sum(ratios) / len(ratios) if ratios else float("nan")
    else:
        drop = float("nan")

    print("=== spectral facts (this is the coarse-mesh accuracy, read off Lambda_h) ===")
    print(f"  1. coarsest mesh (ndof={coarse['ndof']}): low modes n<={LOW_NMAX} "
          f"accurate to {coarse_low*100:.2f}% ; dipole on -2/R to {dipole0['rel_err']*100:.3f}%")
    print(f"  2. error ordered by degree at every mesh (monotone in n): {monotonic}")
    print(f"  3. accurate band per mesh {list(MAXH_SWEEP)} -> {bands} "
          f"(refinement widens it); mean error drop coarse->fine ~ x{drop:.1f}")
    print()
    print("  => low multipoles (which dominate a compact source's far field) are")
    print("     already resolved on the coarsest mesh; refining only widens the band.")
    print()

    # === Part B: the KELVIN closure, measured the same way (volume FEM) ========
    # Kelvin inverts mode n to a degree-n polynomial, so order-p FEM is exact
    # iff p >= n. The dominant dipole inverts to a LINEAR field -> order-1 coarse
    # accurate; the quadrupole needs order >= 2. This is the Kelvin realisation's
    # effective DtN, MEASURED (no longer argued by analogy).
    print("=== Kelvin closure effective DtN (volume FEM), MEASURED vs -(n+1)/R ===")
    print("   Kelvin maps mode n -> degree-n polynomial; order-p FEM exact iff p>=n.\n")
    print("   mode        order  maxh   ndof   lam_eff    rel_err")
    names = {1: "dipole", 2: "quadrupole", 3: "octupole"}
    kel = {}
    # dipole at order 1 across meshes: inverts to LINEAR -> coarse accurate, converges
    for maxh in (0.6, 0.4, 0.25):
        r = kelvin_dtn_eigenvalue(R=R, degree=1, maxh=maxh, order=1)
        kel[("dip", maxh)] = r
        print(f"   {names[1]:10s}  {r['order']}      {maxh:.2f}  {r['ndof']:5d}  "
              f"{r['lam']:+.5f}  {r['rel_err']:.2e}")
    # quadrupole order 1 vs 2 at one mesh: the polynomial-order THRESHOLD
    for order in (1, 2):
        r = kelvin_dtn_eigenvalue(R=R, degree=2, maxh=0.4, order=order)
        kel[("quad", order)] = r
        print(f"   {names[2]:10s}  {r['order']}      0.40  {r['ndof']:5d}  "
              f"{r['lam']:+.5f}  {r['rel_err']:.2e}")
    print()
    print("  => dipole (the dominant mode) inverts to a LINEAR field -> accurate")
    print("     already at order 1 on the coarsest mesh, converging fast (geometry-")
    print("     limited). The quadrupole needs order >= 2. So the Kelvin closure IS")
    print("     coarse-mesh accurate for the modes a compact source actually radiates")
    print("     -- now MEASURED, matching the BEM picture above (shared operator).")
    print()

    # === Part C: open-BC accuracy vs the EXTERIOR (Kelvin ball) mesh size ========
    # Same shell INTERIOR mesh (fixed), refine only the EXTERIOR (Kelvin) region.
    # The interior FEM error cancels, so the residual is the pure exterior (open-BC)
    # error -- Kameari's exterior refinement, ISOLATED. It converges with exterior
    # refinement but stays far below the FIXED interior FEM error at every level.
    print("=== open-BC accuracy vs EXTERIOR (Kelvin ball) mesh size, interior FIXED ===")
    sw = kelvin_openbc_error_vs_exterior_mesh(
        kelvin_maxh_list=(0.7, 0.5, 0.35, 0.25), R_inner=0.5, R_outer=1.0,
        maxh=0.4, order=2, degree=1, kelvin_order=1)
    print(f"   interior FEM error (FIXED, open BC exact)  : {sw['interior_fem_error']*100:6.2f}%")
    print("   ext_maxh  ext_ndof   Kelvin open-BC err   FEM/open-BC")
    for m in sw["per_mesh"]:
        print(f"     {m['kelvin_maxh']:.2f}     {m['kelvin_ndof']:5d}      "
              f"{m['kelvin_openbc_error']*100:7.4f}%        {m['ratio_fem_over_openbc']:5.0f}x")
    print(f"   converges with exterior refinement : {sw['converges']}")
    print(f"   open-BC error < interior FEM error at every exterior mesh : {sw['always_below_fem']}")
    print()
    print("  => refining the EXTERIOR (Kelvin) mesh cuts the open-BC error (~x4/level),")
    print("     but the interior FEM error is FIXED and dominates (45-700x larger).")
    print("     So a COARSE exterior mesh already suffices -- Kameari's coarse-mesh")
    print("     accuracy, verified on the EXTERIOR mesh size as a separated error budget.")
    print()

    # verdict: BEM spectral facts + Kelvin DtN coarse-accuracy + exterior-mesh error budget
    kel_dip_coarse = kel[("dip", 0.6)]["rel_err"]
    kel_dip_fine = kel[("dip", 0.25)]["rel_err"]
    kel_quad_o2 = kel[("quad", 2)]["rel_err"]
    ok = (coarse_low < 0.005 and monotonic and bands[0] >= 2 and bands[-1] >= bands[0]
          and dipole0["rel_err"] < 0.002                       # BEM dipole coarse-accurate
          and kel_dip_coarse < 0.02 and kel_dip_fine < kel_dip_coarse  # Kelvin dipole coarse-acc + converges
          and kel_quad_o2 < 0.01                               # Kelvin quad exact at order 2
          and sw["always_below_fem"] and sw["converges"])      # exterior-mesh error budget
    print("OK" if ok else "OFF -- a measured fact is not as expected, inspect tables above")


if __name__ == "__main__":
    main()
