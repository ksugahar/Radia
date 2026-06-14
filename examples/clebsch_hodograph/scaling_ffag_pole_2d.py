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


def scaling_gap(r, k=K_INDEX, g0=G0, r0=R0, gamma=0.0, gamma2=0.0):
    """Scaling pole gap.  Naive (gamma=gamma2=0): g(r) = g0 (r/r0)^{-k} (B ~ r^k).

    The Step-3 RESHAPE adds a log-r polynomial correction (u = log(r/r0)):
        g(r) = g0 * exp(-k u - gamma/2 u^2 - gamma2/6 u^3),
    so the LOCAL geometric index k_geom = -d log g / d log r tilts as a quadratic
        k_geom(u) = k + gamma*u + gamma2/2 * u^2.
    gamma cancels the LINEAR tilt of k(r), gamma2 the CURVATURE -- the 2-parameter
    reshape that flattens the SATURATED field index (cf. the 2-param Newton in
    clebsch_pole_shape_optimization_2d that nulls b3 AND b5)."""
    u = np.log(np.asarray(r, dtype=float) / r0)
    return g0 * np.exp(-k * u - 0.5 * gamma * u * u - gamma2 * u ** 3 / 6.0)


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
# Step 2: SATURATION -- a super-ferric iron pole (Froehlich mu(B)) droops k(r)
#         at the high-r (high-B) edge as the drive rises.
# --------------------------------------------------------------------------- #
MU0 = 4.0e-7 * math.pi
MUR0_IRON = 2000.0          # Froehlich mu_r at B=0
BK_IRON = 1.2               # Froehlich knee (T)


def mu_r_froehlich(B, mur0=MUR0_IRON, Bk=BK_IRON):
    """Froehlich saturation mu_r(|B|): mur0 at B=0 -> 1 as B -> infinity."""
    return 1.0 + (mur0 - 1.0) / (1.0 + (B / Bk) ** 2)


def _gap_iron_geometry(r_lo, r_hi, k, g0, r0, y_top, n_face=60, gamma=0.0,
                       gamma2=0.0):
    """Upper half: an AIR gap (median y=0 to the scaling pole face y=g(r)/2)
    glued to an IRON pole (pole face up to a flat back at y_top).  The pole
    face is now an interior gap/iron interface (not a BC).  Boundaries:
    'median' (y=0), 'irontop' (the driven back), 'rmin'/'rmax'.  `gamma`,
    `gamma2` are the Step-3 reshape parameters (see scaling_gap)."""
    from netgen.occ import WorkPlane, OCCGeometry, Glue
    rs = np.linspace(r_lo, r_hi, n_face)
    yf = scaling_gap(rs, k, g0, r0, gamma, gamma2) / 2.0
    wp = WorkPlane().MoveTo(r_lo, 0.0)
    wp.LineTo(r_hi, 0.0)                                     # median
    wp.LineTo(r_hi, float(yf[-1]))                          # rmax (gap)
    for ri, yi in zip(rs[-2::-1], yf[-2::-1]):              # pole face (gap top)
        wp.LineTo(float(ri), float(yi))
    wp.Close()
    gap = wp.Face(); gap.faces.name = "gap"
    wp2 = WorkPlane().MoveTo(r_lo, float(yf[0]))
    for ri, yi in zip(rs[1:], yf[1:]):                      # pole face (iron bottom)
        wp2.LineTo(float(ri), float(yi))
    wp2.LineTo(r_hi, y_top)                                 # rmax (iron)
    wp2.LineTo(r_lo, y_top)                                 # irontop
    wp2.Close()
    iron = wp2.Face(); iron.faces.name = "iron"
    shape = Glue([gap, iron])
    for e in shape.edges:
        c = e.center
        if abs(c.y) < 1e-9:
            e.name = "median"
        elif abs(c.y - y_top) < 1e-9:
            e.name = "irontop"
        elif abs(c.x - r_lo) < 1e-6:
            e.name = "rmin"
        elif abs(c.x - r_hi) < 1e-6:
            e.name = "rmax"
        else:
            e.name = "poleface"                            # interior interface
    return OCCGeometry(shape, dim=2)


def solve_saturated(mmf, r_min, r_max, k=K_INDEX, g0=G0, r0=R0, pole_t=0.03,
                    order=3, maxh=0.006, buffer=1.25, n_eval=41,
                    relax=0.5, tol=1e-4, maxit=60, gamma=0.0, gamma2=0.0):
    """Nonlinear scalar-potential solve with a Froehlich iron pole: the magnetic
    scalar potential phi (phi=mmf on the iron back, phi=0 on the median) drives
    flux through the air gap; B = -mu0 mu_r(|B|) grad phi.  Picard on mu_r(B).
    Returns the median field index k(r) and the saturation diagnostics.  `gamma`,
    `gamma2` are the Step-3 pole reshape (0,0 = naive scaling pole)."""
    from ngsolve import (H1, L2, BilinearForm, GridFunction, grad, dx, CF, Norm,
                         TaskManager, Mesh)
    r_lo, r_hi = r_min / buffer, r_max * buffer
    y_top = float(scaling_gap(r_lo, k, g0, r0, gamma, gamma2) / 2.0 + pole_t)
    geo = _gap_iron_geometry(r_lo, r_hi, k, g0, r0, y_top, gamma=gamma,
                             gamma2=gamma2)
    rs_eval = np.linspace(r_min, r_max, n_eval)
    y_probe = 0.02 * float(scaling_gap(r0, k, g0, r0))
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        fes = H1(mesh, order=order, dirichlet="median|irontop")
        u, v = fes.TnT()
        gfu = GridFunction(fes)
        bccf = mesh.BoundaryCF({"irontop": mmf, "median": 0.0}, default=0.0)
        fes_mu = L2(mesh, order=0)
        mu_gf = GridFunction(fes_mu)
        mu_gf.Set(mesh.MaterialCF({"iron": MUR0_IRON}, default=1.0))   # start
        iron_ind = mesh.MaterialCF({"iron": 1.0}, default=0.0)
        resid = 1.0
        it = 0
        for it in range(1, maxit + 1):
            a = BilinearForm(fes)
            a += mu_gf * grad(u) * grad(v) * dx
            a.Assemble()
            gfu.Set(bccf, definedon=mesh.Boundaries("median|irontop"))
            r = gfu.vec.CreateVector()
            r.data = -a.mat * gfu.vec
            gfu.vec.data += a.mat.Inverse(fes.FreeDofs(),
                                          inverse="sparsecholesky") * r
            B = MU0 * mu_gf * Norm(grad(gfu))               # |B| with current mu
            froh = 1.0 + (MUR0_IRON - 1.0) / (1.0 + (B / BK_IRON) ** 2)
            mu_target = (1.0 - iron_ind) * CF(1.0) + iron_ind * froh
            mu_new = GridFunction(fes_mu)
            mu_new.Set(mu_target)
            d = mu_new.vec.CreateVector()
            d.data = mu_new.vec - mu_gf.vec
            denom = (mu_gf.vec.Norm() or 1.0)
            resid = d.Norm() / denom
            mu_gf.vec.data += relax * d
            if resid < tol:
                break
        # median B_y = -mu0 dphi/dy (gap, mu_r=1)
        g = grad(gfu)
        By = np.array([-MU0 * g(mesh(float(r), y_probe))[1] for r in rs_eval])
    rmid, kk = field_index(rs_eval, By)
    return {
        "mmf": float(mmf),
        "rs_eval": rs_eval, "By": By,
        "r_index": rmid, "k": kk,
        "B_gap_min": float(np.min(np.abs(By))),
        "B_gap_max": float(np.max(np.abs(By))),
        "iters": int(it), "resid": float(resid),
        "ndof": fes.ndof,
    }


def run_step2(b_targets=(0.6, 1.4, 2.0, 2.6), k=K_INDEX, g0=G0, r0=R0,
              order=3, maxh=0.006):
    """Sweep the drive (parametrised by the target B_gap at r0); the field index
    k(r) is referenced to the LOWEST (unsaturated) drive, so the saturation
    signal Dk(r) = k(drive) - k(reference) is isolated from the geometric
    baseline.  As the high-r edge (highest B) crosses the iron knee, Dk develops
    a negative dip there: the achromaticity degrades at the high-energy edge of
    the momentum acceptance -- the super-ferric operating wall."""
    r_min, r_max = aperture_radii(k=k, r0=r0)
    mmfs = [bt * g0 / (2.0 * MU0) for bt in b_targets]
    sols = [solve_saturated(mmf, r_min, r_max, k=k, g0=g0, r0=r0,
                            order=order, maxh=maxh) for mmf in mmfs]
    m = slice(2, -2)                                        # interior window
    k_ref = sols[0]["k"]                                    # unsaturated baseline
    levels = []
    for s in sols:
        dk = s["k"] - k_ref
        levels.append({
            "B_target": s["mmf"] * 2.0 * MU0 / g0,
            "B_gap_min": s["B_gap_min"], "B_gap_max": s["B_gap_max"],
            "k_hi_edge": float(s["k"][m][-1]),
            "dk_hi": float(dk[m][-1]),                      # high-r index loss
            "dk_lo": float(dk[m][0]),                       # low-r index loss
            "dk_tilt": float(dk[m][-1] - dk[m][0]),         # differential droop
            "iters": s["iters"], "resid": s["resid"],
            "_s": s,
        })
    return {"k_design": float(k), "aperture": (float(r_min), float(r_max)),
            "Bk_iron": BK_IRON, "k_ref_interior": k_ref[m].tolist(),
            "levels": levels}


# --------------------------------------------------------------------------- #
# Step 3: RESHAPE -- a 1-parameter pole correction (gamma) restores k(r)=const
#         into saturation (the von Mises / log-chart reshape, single-valued).
# --------------------------------------------------------------------------- #
def _index_residual(mmf, gamma, gamma2, r_min, r_max, k, g0, r0, order, maxh):
    """Saturated solve at reshape (gamma, gamma2); return the field-index
    flatness residual (tilt, curvature), the peak-to-peak, and the solve.
      tilt = k(r_hi) - k(r_lo)        (the linear slope of k(r))
      curv = k(r_mid) - mean(edges)   (the bow / quadratic part)."""
    s = solve_saturated(mmf, r_min, r_max, k=k, g0=g0, r0=r0,
                        order=order, maxh=maxh, gamma=gamma, gamma2=gamma2)
    ki = s["k"][2:-2]
    mid = len(ki) // 2
    tilt = float(ki[-1] - ki[0])
    curv = float(ki[mid] - 0.5 * (ki[-1] + ki[0]))
    return tilt, curv, float(ki.max() - ki.min()), s


def run_step3(B_design=1.8, k=K_INDEX, g0=G0, r0=R0, order=3, maxh=0.008,
              tol=3e-3, newton_steps=3, dg=0.4):
    """At the design (saturated) excitation, find the 2-parameter pole reshape
    (gamma, gamma2) that drives BOTH the field-index tilt AND curvature to zero
    -- restoring a flat (achromatic) k(r) into saturation.  The local geometric
    index k_geom(u) = k + gamma*u + gamma2/2*u^2 (u=log r/r0) has 2 knobs; a 2-D
    Newton on (tilt, curv)=0 (finite-difference Jacobian) flattens the saturated
    index, the single-valued von Mises/log-chart reshape.  Returns naive vs
    reshaped flatness."""
    r_min, r_max = aperture_radii(k=k, r0=r0)
    mmf = B_design * g0 / (2.0 * MU0)

    def resid(g, g2):
        return _index_residual(mmf, g, g2, r_min, r_max, k, g0, r0, order, maxh)

    t0, c0, p0, s0 = resid(0.0, 0.0)                        # naive pole
    g, g2 = 0.0, 0.0
    t, c, p, s = t0, c0, p0, s0
    evals = [{"gamma": 0.0, "gamma2": 0.0, "tilt": t0, "curv": c0, "ptp": p0}]
    best = {"gamma": 0.0, "gamma2": 0.0, "tilt": t0, "curv": c0, "ptp": p0,
            "_s": s0}
    for _ in range(newton_steps):
        if max(abs(t), abs(c)) < tol:
            break
        # finite-difference 2x2 Jacobian d(tilt,curv)/d(gamma,gamma2)
        tg, cg, _, _ = resid(g + dg, g2)
        tg2, cg2, _, _ = resid(g, g2 + dg)
        J = np.array([[(tg - t) / dg, (tg2 - t) / dg],
                      [(cg - c) / dg, (cg2 - c) / dg]])
        try:
            step = np.linalg.solve(J, np.array([t, c]))
        except np.linalg.LinAlgError:
            break
        g, g2 = g - step[0], g2 - step[1]
        t, c, p, s = resid(g, g2)
        evals.append({"gamma": g, "gamma2": g2, "tilt": t, "curv": c, "ptp": p})
        if p < best["ptp"]:
            best = {"gamma": g, "gamma2": g2, "tilt": t, "curv": c, "ptp": p,
                    "_s": s}
    return {
        "k_design": float(k), "B_design": float(B_design),
        "aperture": (float(r_min), float(r_max)),
        "naive": {"gamma": 0.0, "gamma2": 0.0, "tilt": t0, "curv": c0,
                  "ptp": p0, "_s": s0},
        "reshaped": best,
        "ptp_improvement": float(p0 / best["ptp"]) if best["ptp"] > 0 else None,
        "evals": evals,
    }


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
    ap.add_argument("--step2", action="store_true",
                    help="also run the saturation sweep (Froehlich iron pole)")
    ap.add_argument("--step3", action="store_true",
                    help="also run the 2-parameter pole reshape (restore flat k)")
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

    if args.step2:
        print("\n" + "=" * 74)
        print("Step 2: SATURATION -- Froehlich iron pole droops k at the high-r edge")
        print("=" * 74)
        s2 = run_step2()
        print(f"iron knee Bk                : {s2['Bk_iron']:.2f} T")
        print(f"{'B_tgt(T)':<9}{'B_gap[min,max]':<16}{'k_hi':<8}"
              f"{'dk_hi':<9}{'dk_tilt':<9}{'iters':<6}")
        for L in s2["levels"]:
            print(f"{L['B_target']:<9.1f}[{L['B_gap_min']:.2f}, {L['B_gap_max']:.2f}]"
                  f"   {L['k_hi_edge']:<8.3f}{L['dk_hi']:<+9.3f}"
                  f"{L['dk_tilt']:<+9.3f}{L['iters']:<6}")
        if args.fig:
            _figure_step2(s2)

    if args.step3:
        print("\n" + "=" * 74)
        print("Step 3: RESHAPE -- 2-param pole correction restores flat k(r) in saturation")
        print("=" * 74)
        s3 = run_step3()
        n, o = s3["naive"], s3["reshaped"]
        print(f"design excitation B@r0      : {s3['B_design']:.1f} T")
        print(f"naive pole   k(r) ptp       : {n['ptp']:.4f}"
              f"  (tilt {n['tilt']:+.4f}, curv {n['curv']:+.4f})")
        print(f"reshaped     k(r) ptp       : {o['ptp']:.4f}"
              f"  (tilt {o['tilt']:+.4f}, curv {o['curv']:+.4f})")
        print(f"reshape (gamma, gamma2)     : ({o['gamma']:+.3f}, {o['gamma2']:+.3f})")
        print(f"flatness improvement        : {s3['ptp_improvement']:.1f}x"
              f"  ({len(s3['evals'])-1} Newton step(s))")
        if args.fig:
            _figure_step3(s3)


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


def _figure_step2(s2):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    k_ref = np.array(s2["k_ref_interior"])
    for L in s2["levels"]:
        s = L["_s"]
        m = slice(2, -2)
        ri = s["r_index"][m]
        dk = s["k"][m] - k_ref
        ax[0].plot(ri, dk, "o-", ms=3,
                   label=f"B@r0={L['B_target']:.1f}T (max {L['B_gap_max']:.1f})")
        ax[1].plot(s["rs_eval"], np.abs(s["By"]), "-",
                   label=f"B@r0={L['B_target']:.1f}T")
    ax[0].axhline(0, color="k", lw=0.8)
    ax[0].set_xlabel("r"); ax[0].set_ylabel("dk(r) = k - k(unsaturated)")
    ax[0].set_title("field-index droop (high-r edge saturates)")
    ax[0].legend(fontsize=7)
    ax[1].axhline(s2["Bk_iron"], color="r", ls="--", lw=1, label="iron knee Bk")
    ax[1].set_xlabel("r"); ax[1].set_ylabel("|B_gap|(r)  [T]")
    ax[1].set_title("gap field: high-r end crosses the knee")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "scaling_ffag_pole_2d_saturation.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


def _figure_step3(s3):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n, o = s3["naive"], s3["reshaped"]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for tag, d, style in (("naive", n, "o-"), ("reshaped", o, "s-")):
        s = d["_s"]
        m = slice(2, -2)
        ax[0].plot(s["r_index"][m], s["k"][m], style, ms=3,
                   label=f"{tag} (ptp {d['ptp']:.3f})")
    ax[0].axhline(s3["k_design"], color="k", ls="--", lw=1,
                  label=f"k_design={s3['k_design']:.0f}")
    ax[0].set_xlabel("r"); ax[0].set_ylabel("field index k(r)")
    ax[0].set_title(f"saturated k(r) flattened {s3['ptp_improvement']:.1f}x"
                    f" @ B@r0={s3['B_design']:.1f}T")
    ax[0].legend(fontsize=8)
    # pole gap profiles (the shape change)
    r_min, r_max = s3["aperture"]
    rs = np.linspace(r_min, r_max, 60)
    ax[1].plot(rs, scaling_gap(rs, K_INDEX, G0, R0) / 2.0, "o-", ms=3,
               label="naive g(r)/2")
    ax[1].plot(rs, scaling_gap(rs, K_INDEX, G0, R0, o["gamma"], o["gamma2"]) / 2.0,
               "s-", ms=3, label="reshaped g(r)/2")
    ax[1].set_xlabel("r"); ax[1].set_ylabel("pole-face height g(r)/2")
    ax[1].set_title("the pole reshape (log-chart 2-param)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "scaling_ffag_pole_2d_reshape.png")
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
