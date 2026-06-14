r"""Scaling-FFAG super-ferric proton-gantry pole: achromatic design via the hodograph,
certified by COMPLEMENTARY (A vs phi) bounds.

WHY THIS EXAMPLE
----------------
A fixed-field (scaling-FFAG) gantry delivers a RANGE of beam momenta WITHOUT
re-exciting the magnets -- the win for a compact proton therapy gantry (fast
energy switching, no per-energy ramp).  The enabling property is ACHROMATICITY
= zero chromaticity = the betatron tune is momentum-independent.  For a scaling
field that is exactly

    B_y(r) = B0 (r/r0)^k          (k = field index),

orbits of different momenta are GEOMETRICALLY SIMILAR (p ~ r^{k+1}), so the tune
is the same for every momentum.  Achromaticity is therefore EXACTLY the single
condition

    k(r) = d log B_y / d log r = const            across the radial aperture.

THE HODOGRAPH CONNECTION (why log-r is the natural chart)
--------------------------------------------------------
Momentum enters as a SCALING of the orbit (r ~ p^{1/(k+1)}).  In u = log r the
scaling r -> lam r becomes a TRANSLATION u -> u + log lam, and the scaling field
B ~ r^k = e^{k u} is translation-covariant: log B is a STRAIGHT LINE in u with
slope k.  So
    achromatic  <=>  log B vs log r is a straight line (k = const = slope).
The pole gap g(r) ~ r^{-k} (B ~ 1/g) is, in u, log g linear in u -- a "uniform"
pole in the log chart.  This is the 2-variable conformal (log) chart; the
non-linear (saturation) reshape in Steps 2-3 drops to the single-variable von
Mises chart, because the full 2-variable swap FOLDS once mu = mu(q) (see
chaplygin_inverse_vonmises_2d.py).

COMPLEMENTARY (A vs phi) CERTIFICATION
--------------------------------------
The achromaticity of a FINITE pole is a numerical claim, so we bracket it.  The
SAME gap is solved two complementary ways (the energy <-> co-energy Legendre
pair; A and phi are the conjugate (flux-function, scalar-potential) pair of the
hodograph):

  * phi-formulation : scalar potential, Dirichlet phi on the poles
                      (the high-mu equipotential), B = -grad phi.
  * A-formulation   : flux function A_z, Dirichlet A on the flux walls
                      (the conjugate BCs), B = (A_y, -A_x).

The two solutions converge to the same field from discretisation-COMPLEMENTARY
sides; the global magnetic energy is rigorously bracketed (Synge hypercircle /
Rikabi-Bryant-Freeman, monotone BH => convex energy => the bracket survives into
saturation).  The LOCAL field index k(r) is not a strict two-sided theorem, but
k_A(r) and k_phi(r) converge from opposite sides, so their GAP certifies that a
flat k(r) is physics, not a mesh artefact -- and (Steps 2-3) drives the reshape.

NOTE the field-index k is SCALE-INVARIANT (d log B / d log r ignores the overall
amplitude), so phi (unit potential) and A (unit flux) need NOT share a
normalisation to be compared -- only the field SHAPE matters.

SCOPE (this file = Step 1, the linear foundation)
-------------------------------------------------
  Step 1 (HERE): linear high-mu scaling pole; measure k_phi(r), k_A(r); show the
                 A/phi bracket certifies how flat the naive g ~ r^{-k} pole is.
  Step 2 (next): Froehlich mu(B) saturation -> k droops at the high-r (high-B)
                 edge = the super-ferric operating wall.
  Step 3 (next): reshape the pole (von Mises / log chart, no remesh) to restore
                 k(r) = const into saturation.

run:  python scaling_ffag_pole_2d.py            # complementary bracket (fast)
      python scaling_ffag_pole_2d.py --fig       # + figure
"""
import argparse
import math
import os
import sys

import numpy as np

# --------------------------------------------------------------------------- #
# proton gantry parameters (the radial aperture from the momentum range)
# --------------------------------------------------------------------------- #
M_PROTON = 938.272            # MeV
K_INDEX = 5.0                 # scaling field index B ~ r^k
T_LO, T_HI = 70.0, 250.0      # proton kinetic-energy range (MeV) -> momentum band
R0 = 1.0                      # reference radius (normalised)
G0 = 0.10                     # gap at r0 (normalised); B ~ 1/g


def _rigidity_ratio(T_lo=T_LO, T_hi=T_HI):
    """p_hi / p_lo for kinetic energies T_hi, T_lo (relativistic protons)."""
    p = lambda T: math.sqrt(T * (T + 2 * M_PROTON))        # pc in MeV  # noqa: E731
    return p(T_hi) / p(T_lo)


def aperture_radii(k=K_INDEX, r0=R0, t_lo=T_LO, t_hi=T_HI):
    """Radial aperture [r_min, r_max] (centred on r0) spanned by the momentum
    band, using the scaling relation r ~ p^{1/(k+1)}."""
    ratio = _rigidity_ratio(t_lo, t_hi) ** (1.0 / (k + 1.0))   # r_hi / r_lo
    return r0 / math.sqrt(ratio), r0 * math.sqrt(ratio)


def scaling_gap(r, k=K_INDEX, g0=G0, r0=R0):
    """Naive scaling pole gap g(r) = g0 (r/r0)^{-k}  (0th order: B ~ 1/g ~ r^k)."""
    return g0 * (np.asarray(r, dtype=float) / r0) ** (-k)


# --------------------------------------------------------------------------- #
# the new field-quality metric: field index k(r) = d log B_y / d log r
# --------------------------------------------------------------------------- #
def field_index(rs, By, r_fit=None):
    """Local field index k(r) = d(log B_y)/d(log r), by a centred log-log slope.

    rs, By : 1-D arrays (By > 0 on the aperture).  Returns k at the interior
    sample points (length len(rs) - 2 by central difference), with rs_mid."""
    rs = np.asarray(rs, dtype=float)
    By = np.asarray(By, dtype=float)
    lr, lB = np.log(rs), np.log(np.abs(By))
    k = (lB[2:] - lB[:-2]) / (lr[2:] - lr[:-2])              # central difference
    return rs[1:-1], k


# --------------------------------------------------------------------------- #
# complementary 2-D solves on the scaling gap (upper half, high-mu limit)
# --------------------------------------------------------------------------- #
def _gap_geometry(r_min, r_max, k, g0, r0, n_face=60):
    """Upper-half air gap between the median (y=0) and the scaling pole face
    y = g(r)/2.  Boundaries: 'median' (y=0), 'pole' (the face), 'rmin'/'rmax'
    (the radial walls = flux walls)."""
    from netgen.occ import WorkPlane, OCCGeometry
    rs = np.linspace(r_min, r_max, n_face)
    yf = scaling_gap(rs, k, g0, r0) / 2.0
    wp = WorkPlane().MoveTo(r_min, 0.0)
    wp.LineTo(r_max, 0.0)                                    # median  (y=0)
    wp.LineTo(r_max, float(yf[-1]))                         # rmax wall
    for ri, yi in zip(rs[-2::-1], yf[-2::-1]):              # pole face (back)
        wp.LineTo(float(ri), float(yi))
    wp.LineTo(r_min, 0.0)                                    # rmin wall
    face = wp.Face()
    face.faces.name = "gap"
    # name edges by midpoint geometry (no hardcoded ids)
    for e in face.edges:
        c = e.center
        if abs(c.y) < 1e-9:
            e.name = "median"
        elif abs(c.x - r_min) < 1e-6:
            e.name = "rmin"
        elif abs(c.x - r_max) < 1e-6:
            e.name = "rmax"
        else:
            e.name = "pole"
    return OCCGeometry(face, dim=2)


def _median_By(gfu, kind, rs_eval, mesh, y_probe):
    """Sample B_y on the median plane.  phi: B_y = -dphi/dy ; A: B_y = -dA/dx."""
    from ngsolve import grad
    g = grad(gfu)
    out = []
    for r in rs_eval:
        mp = mesh(float(r), float(y_probe))
        gv = g(mp)
        out.append(-gv[1] if kind == "phi" else -gv[0])
    return np.array(out)


def solve_complementary(r_min, r_max, k=K_INDEX, g0=G0, r0=R0,
                        order=4, maxh=0.02, n_eval=41, buffer=1.25):
    """Solve the scaling gap with BOTH complementary formulations and return the
    median field index k_phi(r), k_A(r).

    The Laplace problem is solved on a radial domain BUFFERED beyond the
    measurement aperture [r_min, r_max] (by the factor `buffer` each side) so the
    radial-wall fringing stays OUT of the measurement window -- the field index
    is then read in the clean bulk.

    phi : Dirichlet phi=0 on median, phi=1 on pole (equipotential), natural on
          the radial walls.            B = -grad phi.
    A   : Dirichlet A=0 on rmin, A=1 on rmax (flux specified), natural on median
          and pole (the dual BCs).     B = (A_y, -A_x)  ->  median B_y = -A_x.
    """
    from ngsolve import (H1, BilinearForm, GridFunction, grad, dx, TaskManager,
                         Mesh)
    r_lo, r_hi = r_min / buffer, r_max * buffer            # solve domain
    geo = _gap_geometry(r_lo, r_hi, k, g0, r0)
    rs_eval = np.linspace(r_min, r_max, n_eval)            # measurement window
    y_probe = 0.02 * scaling_gap(r0, k, g0, r0)            # just off the median
    out = {}
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        # ---- phi-formulation (Dirichlet on poles) ----
        fes_p = H1(mesh, order=order, dirichlet="median|pole")
        u, v = fes_p.TnT()
        a = BilinearForm(fes_p)
        a += grad(u) * grad(v) * dx
        a.Assemble()
        gfp = GridFunction(fes_p)
        gfp.Set(mesh.BoundaryCF({"pole": 1.0, "median": 0.0}, default=0.0),
                definedon=mesh.Boundaries("median|pole"))
        r = gfp.vec.CreateVector()
        r.data = -a.mat * gfp.vec
        gfp.vec.data += a.mat.Inverse(fes_p.FreeDofs(),
                                      inverse="sparsecholesky") * r
        By_p = _median_By(gfp, "phi", rs_eval, mesh, y_probe)
        # ---- A-formulation (Dirichlet on flux walls = the dual BCs) ----
        fes_a = H1(mesh, order=order, dirichlet="rmin|rmax")
        u, v = fes_a.TnT()
        a2 = BilinearForm(fes_a)
        a2 += grad(u) * grad(v) * dx
        a2.Assemble()
        gfa = GridFunction(fes_a)
        gfa.Set(mesh.BoundaryCF({"rmax": 1.0, "rmin": 0.0}, default=0.0),
                definedon=mesh.Boundaries("rmin|rmax"))
        r2 = gfa.vec.CreateVector()
        r2.data = -a2.mat * gfa.vec
        gfa.vec.data += a2.mat.Inverse(fes_a.FreeDofs(),
                                       inverse="sparsecholesky") * r2
        By_a = _median_By(gfa, "A", rs_eval, mesh, y_probe)
    rmid_p, k_p = field_index(rs_eval, By_p)
    rmid_a, k_a = field_index(rs_eval, By_a)
    out["rs_eval"] = rs_eval
    out["By_phi"] = By_p
    out["By_A"] = By_a
    out["r_index"] = rmid_p
    out["k_phi"] = k_p
    out["k_A"] = k_a
    out["ndof"] = (fes_p.ndof, fes_a.ndof)
    return out


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def run(k=K_INDEX, g0=G0, r0=R0, order=4, maxh=0.008):
    """Step 1: linear scaling pole; analytic check + complementary k(r) bracket."""
    # (i) analytic check: field_index on B ~ r^k returns k.
    rs = np.linspace(0.5, 2.0, 200)
    k_meas = field_index(rs, (rs / r0) ** k)[1]
    analytic_err = float(np.max(np.abs(k_meas - k)))

    # (ii) physical aperture from the proton momentum band.
    r_min, r_max = aperture_radii(k=k, r0=r0)
    sol = solve_complementary(r_min, r_max, k=k, g0=g0, r0=r0,
                              order=order, maxh=maxh)

    # certification: where the A/phi bracket is tightest, k is best resolved.
    k_lo = np.minimum(sol["k_phi"], sol["k_A"])
    k_hi = np.maximum(sol["k_phi"], sol["k_A"])
    k_mid = 0.5 * (sol["k_phi"] + sol["k_A"])
    gap = float(np.max(k_hi - k_lo))
    # interior window (drop the 2 edge points each side -- wall fringing)
    m = slice(2, -2)
    k_flat_dev = float(np.max(np.abs(k_mid[m] - k)))         # how flat (vs k)
    return {
        "k_design": float(k),
        "aperture": (float(r_min), float(r_max)),
        "rigidity_ratio": float(_rigidity_ratio()),
        "analytic_index_err": analytic_err,
        "bracket_gap_max": gap,
        "k_interior_dev_from_design": k_flat_dev,
        "k_interior_mid_mean": float(np.mean(k_mid[m])),
        "ndof": sol["ndof"],
        "_sol": sol,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fig", action="store_true")
    ap.add_argument("--order", type=int, default=4)
    ap.add_argument("--maxh", type=float, default=0.008)
    args = ap.parse_args()

    print("=" * 74)
    print("Scaling-FFAG proton-gantry pole -- Step 1: linear scaling + A/phi bracket")
    print("=" * 74)
    res = run(order=args.order, maxh=args.maxh)
    rmin, rmax = res["aperture"]
    print(f"field index k_design        : {res['k_design']:.3f}")
    print(f"momentum band p_hi/p_lo     : {res['rigidity_ratio']:.3f}"
          f"  (T {T_LO:.0f}-{T_HI:.0f} MeV)")
    print(f"radial aperture [r_min,r_max]: [{rmin:.4f}, {rmax:.4f}]"
          f"  (ratio {rmax / rmin:.3f})")
    print(f"analytic field_index error  : {res['analytic_index_err']:.2e}"
          f"  (B~r^k -> k)")
    print(f"ndof (phi, A)               : {res['ndof']}")
    print(f"A/phi bracket gap (max)     : {res['bracket_gap_max']:.3e}")
    print(f"k interior mean (mid)       : {res['k_interior_mid_mean']:.3f}")
    print(f"k interior dev from design  : {res['k_interior_dev_from_design']:.3e}")

    if args.fig:
        _figure(res)


def _figure(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sol = res["_sol"]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    rs = sol["rs_eval"]
    ax[0].loglog(rs, np.abs(sol["By_phi"]), "o-", ms=3, label="phi (B=-grad phi)")
    ax[0].loglog(rs, np.abs(sol["By_A"]), "s-", ms=3, label="A (B=curl A)")
    ax[0].loglog(rs, (rs / res["k_design"] * 0 + 1) * np.abs(sol["By_phi"][0])
                 * (rs / rs[0]) ** res["k_design"], "k--", lw=1,
                 label=f"ideal r^{res['k_design']:.0f}")
    ax[0].set_xlabel("r"); ax[0].set_ylabel("|B_y| (median)")
    ax[0].legend(fontsize=8)
    ri = sol["r_index"]
    ax[1].plot(ri, sol["k_phi"], "o-", ms=3, label="k_phi")
    ax[1].plot(ri, sol["k_A"], "s-", ms=3, label="k_A")
    ax[1].axhline(res["k_design"], color="k", ls="--", lw=1,
                  label=f"k_design={res['k_design']:.0f}")
    ax[1].fill_between(ri, np.minimum(sol["k_phi"], sol["k_A"]),
                       np.maximum(sol["k_phi"], sol["k_A"]),
                       alpha=0.2, label="A/phi bracket")
    ax[1].set_xlabel("r"); ax[1].set_ylabel("field index k(r)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "scaling_ffag_pole_2d.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
