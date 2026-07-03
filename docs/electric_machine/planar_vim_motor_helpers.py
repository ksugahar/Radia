"""Reference and coupling helpers for the planar HDiv-VIM motor showcase notebook.

The SOLVER lives in radia.vim (PlanarDemagBody / solve_planar_demag / maxwell_torque_circle --
the promoted 2D layer); this module holds the APPLICATION-side pieces the showcase compares
against and couples with:

* Case A (salient-bar reluctance motor): a 6-wire sin-distribution stator + an ALL-IN-ONE
  exact-Newton nonlinear A_z FEM reference (closed-form nu(B) inversion of the saturating law,
  coil disks curved + measured-area J normalization -- an inscribed-polygon disk under-carries
  its current by ~5%).
* Case B (induction physics): the boundary-matching Bessel closed form for a solid conducting
  cylinder in a uniform ROTATING field, and the reduced complex A_z solve of the same problem.
* Case C (mini cage induction machine): iron core (VIM) + conducting bar ring (reduced complex
  A_z FEM) weak coupling -- including the SINGLE-VALUED polar construction of the iron's
  conjugate potential (the plain atan2 formula has per-charge branch cuts that a SURROUNDING
  bar ring crosses) -- and the all-in-one frozen-secant FEM reference.

TaskManager: helpers do NOT wrap; the notebook cells wrap (caller-wraps policy).
"""
from __future__ import annotations

import numpy as np
from scipy.special import jv, jvp
from scipy.interpolate import RegularGridInterpolator

import ngsolve as ng
from ngsolve import solvers as ngsolvers
from netgen.occ import WorkPlane, OCCGeometry, Glue

MU0 = 4e-7 * np.pi

# =====================================================================================
# saturating law (shared by every nonlinear case in the showcase)
# =====================================================================================
CHI0, MSAT = 1000.0, 1.2e6
K_LAW = CHI0 / MSAT


def bh_table(n=400):
    """The saturating law M(H) = chi0 H / (1 + chi0 H / Msat) as a [[H, B]] table."""
    H = np.logspace(0, 7.5, n)
    M = CHI0 * H / (1.0 + CHI0 * H / MSAT)
    return np.stack([H, MU0 * (H + M)], axis=1).tolist()


def M_of_h(h):
    return CHI0 * np.asarray(h, float) / (1.0 + CHI0 * np.abs(h) / MSAT)


def chi_sec(h):
    return CHI0 / (1.0 + CHI0 * np.asarray(h, float) / MSAT)


def fixed_point_disk(H0, D=0.5):
    """Analytic uniform-body root M = Mof(H0 - D*M) (odd law; bracket [0, H0/D])."""
    lo, hi = 0.0, H0 / D
    f = lambda M: M - float(M_of_h(H0 - D * M))
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _nu_rotor_cf(gradu):
    """Exact reluctivity CF of the saturating law: B(H) inverts in closed form
    (k H^2 + (1 + chi0 - k B/mu0) H - B/mu0 = 0)."""
    b = ng.sqrt(gradu * gradu + 1e-18)
    beta = (1.0 + CHI0) - K_LAW * b / MU0
    Hb = (-beta + ng.sqrt(beta * beta + 4.0 * K_LAW * b / MU0)) / (2.0 * K_LAW)
    return Hb / b, b


# =====================================================================================
# Case A: salient-bar reluctance motor (6-wire sin stator)
# =====================================================================================
A_BAR, B_BAR = 0.24, 0.12          # rotor bar (salient, no closed form)
RS, RC_A = 0.28, 0.20              # wire radius / Maxwell gap circle
NW, K0 = 6, 8.0e5                  # wires; sheet amplitude [A/m] -> H_center ~ 4e5 (knee)
WIRE_R = 0.012
PHI_K = np.radians(30.0 + 60.0 * np.arange(NW))
I_K = K0 * np.sin(PHI_K) * RS * (2 * np.pi / NW)


def wire_pos(theta_deg):
    """Wire positions in the ROTOR frame (rotor rotated +theta == wires rotated -theta)."""
    a = PHI_K - np.radians(theta_deg)
    return np.stack([RS * np.cos(a), RS * np.sin(a)], axis=1)


def H_wires_np(P, W):
    H = np.zeros_like(P)
    for (cx, cy), I in zip(W, I_K):
        dx = P[:, 0] - cx
        dy = P[:, 1] - cy
        r2 = dx * dx + dy * dy
        H[:, 0] += I / (2 * np.pi) * (-dy / r2)
        H[:, 1] += I / (2 * np.pi) * (dx / r2)
    return H


def H_wires_cf(W):
    terms = [ng.CoefficientFunction((I / (2 * np.pi) * (-(ng.y - cy)),
                                     I / (2 * np.pi) * (ng.x - cx)))
             / ((ng.x - cx) ** 2 + (ng.y - cy) ** 2) for (cx, cy), I in zip(W, I_K)]
    cf = terms[0]
    for t in terms[1:]:
        cf = cf + t
    return cf


def _bar_case_setup(theta_deg, maxh_box=0.5, rotor_h=0.02, airin_h=0.04, order=2):
    """Shared Case-A geometry + space + normalized coil source.
    Lessons baked in: coil disks are CURVED and J is normalized by the MEASURED disk area
    (an inscribed-polygon disk under-carries I by ~5% -> torque -10%)."""
    W = wire_pos(theta_deg)
    box = WorkPlane().RectangleC(8, 8).Face()
    box.edges.name = "outer"
    rotor = WorkPlane().RectangleC(A_BAR, B_BAR).Face()
    rotor.faces.name = "rotor"
    rotor.faces.maxh = rotor_h
    ring = WorkPlane().Circle(0, 0, 0.5).Face()
    coils = []
    for k, ((cx, cy), I) in enumerate(zip(W, I_K)):
        c = WorkPlane().Circle(cx, cy, WIRE_R).Face()
        c.faces.name = f"coil{k}"
        c.faces.maxh = 0.008
        coils.append(c)
    airin = ring - rotor
    for c in coils:
        airin = airin - c
    airin.faces.name = "airin"
    airin.faces.maxh = airin_h
    airout = box - ring
    airout.faces.name = "airout"
    geo = Glue([airout, airin, rotor] + coils)
    mesh = ng.Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=maxh_box))
    mesh.Curve(3)
    fes = ng.H1(mesh, order=order, dirichlet="outer")
    coil_area = {k: ng.Integrate(ng.CoefficientFunction(1.0), mesh,
                                 definedon=mesh.Materials(f"coil{k}")) for k in range(NW)}
    Jcf = mesh.MaterialCF({f"coil{k}": I_K[k] / coil_area[k] for k in range(NW)}, default=0.0)
    return mesh, fes, Jcf


def _maxwell_torque_gradA(mesh, gfA, Rc=None, nphi=1440):
    Rc = RC_A if Rc is None else Rc
    acc = 0.0
    for p in np.linspace(0, 2 * np.pi, nphi, endpoint=False):
        x, y = Rc * np.cos(p), Rc * np.sin(p)
        g = ng.grad(gfA)(mesh(x, y))
        Bx, By = g[1], -g[0]
        Br = Bx * np.cos(p) + By * np.sin(p)
        Bp = -Bx * np.sin(p) + By * np.cos(p)
        acc += Br * Bp
    return Rc * Rc / MU0 * (2 * np.pi / nphi) * acc


def fem_reference_bar(theta_deg, maxh_box=0.5, rotor_h=0.02, airin_h=0.04, order=2,
                      nphi=1440):
    """ALL-IN-ONE exact-Newton nonlinear A_z reference for Case A; returns (torque, ndof, iters).
    The closed-form nu(B) inversion makes this a true Newton (quadratic, 6-9 iters cold-start
    at deep saturation); torque via the Maxwell circle at RC_A from point-sampled grad A."""
    mesh, fes, Jcf = _bar_case_setup(theta_deg, maxh_box, rotor_h, airin_h, order)
    u, v = fes.TnT()
    nu_rot, _ = _nu_rotor_cf(ng.grad(u))
    ind = mesh.MaterialCF({"rotor": 1.0}, default=0.0)
    nu_tot = ind * nu_rot + (1.0 - ind) * (1.0 / MU0)
    a = ng.BilinearForm(fes)
    a += (nu_tot * ng.grad(u) * ng.grad(v) - Jcf * v) * ng.dx
    gfA = ng.GridFunction(fes)
    ret = ngsolvers.Newton(a, gfA, maxit=100, maxerr=1e-9, dampfactor=0.6,
                           printing=False, inverse="sparsecholesky")
    status, iters = (ret if isinstance(ret, tuple) else (ret, -1))
    if status != 0:
        raise RuntimeError(f"Case-A FEM Newton not converged at theta={theta_deg}")
    return _maxwell_torque_gradA(mesh, gfA, nphi=nphi), fes.ndof, iters


def fem_reference_bar_secant(theta_deg, nouter=60, relax=0.3, maxh_box=0.5, rotor_h=0.02,
                             airin_h=0.04, order=2):
    """The DELIBERATELY-KEPT failure mode for the reference-audit demonstration: the
    per-element secant-nu Picard on the same Case-A problem.  At knee-level drive the corner
    elements swing across the BH knee and the iteration PLATEAUS (dA ~ 0.1) -- returns the
    residual history instead of raising, so the plateau can be shown next to the exact Newton.
    Do NOT use this as a reference; see bug pattern `reference-secant-picard-oscillation`."""
    mesh, fes, Jcf = _bar_case_setup(theta_deg, maxh_box, rotor_h, airin_h, order)
    u, v = fes.TnT()
    fes0 = ng.L2(mesh, order=0)
    gfnu = ng.GridFunction(fes0)
    areas = ng.Integrate(ng.CoefficientFunction(1.0), mesh, element_wise=True)
    areas = np.array([areas[k] for k in range(mesh.ne)])
    el2dof = np.zeros(mesh.ne, dtype=int)
    is_rotor = np.zeros(mesh.ne, dtype=bool)
    for el in mesh.Elements(ng.VOL):
        el2dof[el.nr] = fes0.GetDofNrs(el)[0]
        is_rotor[el.nr] = mesh[el].mat == "rotor"
    H_TAB = np.concatenate([[0.0], np.logspace(0, 7.5, 400)])
    B_TAB = MU0 * (H_TAB + CHI0 * H_TAB / (1.0 + CHI0 * H_TAB / MSAT))
    NU_TAB = np.empty_like(H_TAB)
    NU_TAB[1:] = H_TAB[1:] / B_TAB[1:]
    NU_TAB[0] = 1.0 / (MU0 * (1.0 + CHI0))
    f = ng.LinearForm(fes)
    f += Jcf * v * ng.dx
    f.Assemble()
    nu_e = np.full(mesh.ne, 1.0 / MU0)
    nu_e[is_rotor] = NU_TAB[0]
    gfA = ng.GridFunction(fes)
    hist = []
    for _ in range(nouter):
        gfnu.vec.FV().NumPy()[el2dof] = nu_e
        a = ng.BilinearForm(fes)
        a += gfnu * ng.grad(u) * ng.grad(v) * ng.dx
        a.Assemble()
        Aold = gfA.vec.FV().NumPy().copy()
        gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
        dA = np.linalg.norm(gfA.vec.FV().NumPy() - Aold) / max(
            np.linalg.norm(gfA.vec.FV().NumPy()), 1e-300)
        hist.append(dA)
        Bint = ng.Integrate(ng.Norm(ng.grad(gfA)), mesh, element_wise=True)
        Bel = np.array([Bint[k] for k in range(mesh.ne)]) / areas
        nu_new = nu_e.copy()
        nu_new[is_rotor] = np.interp(Bel[is_rotor], B_TAB, NU_TAB)
        nu_e = np.exp((1 - relax) * np.log(nu_e) + relax * np.log(nu_new))
    return np.array(hist), fes.ndof


# =====================================================================================
# Case B: solid conducting cylinder in a uniform ROTATING field (induction physics)
# =====================================================================================
A_CYL, SIGMA_CYL, B0_CYL = 0.10, 3.5e7, 0.05
H0_CYL = B0_CYL / MU0


def cylinder_analytic_T(omega, mu_r=1.0, a=A_CYL, sigma=SIGMA_CYL):
    """Boundary-matching Bessel closed form; time-averaged torque T = 2 pi H0 Re(beta)
    (beta = the exterior dipole coefficient; R-independence is analytic)."""
    mu = MU0 * mu_r
    k = np.sqrt(-1j * omega * mu * sigma)
    ka = k * a
    c = B0_CYL
    M = np.array([[jv(1, ka), -1.0 / a],
                  [(k / mu) * jvp(1, ka), 1.0 / (MU0 * a * a)]], dtype=complex)
    alpha, beta = np.linalg.solve(M, np.array([1j * c * a, 1j * c / MU0], dtype=complex))
    return 2.0 * np.pi * H0_CYL * float(np.real(beta))


def cylinder_reduced_T(omega, mu_r=1.0, L=2.0, maxh_cond=0.004, Rc=0.15, nphi=1440):
    """Reduced complex A_z solve of the rotating-field cylinder (A_s = B0 (y + j x)); the
    conductor's mu_r goes into the FEM reluctivity, and for mu_r != 1 the split adds the weak
    residual -(nu - nu0) grad(A_s).grad(v).  Time-averaged Maxwell-circle torque."""
    cond = WorkPlane().Circle(0, 0, A_CYL).Face()
    cond.faces.name = "cond"
    cond.faces.maxh = maxh_cond
    box = WorkPlane().RectangleC(2 * L, 2 * L).Face()
    box.edges.name = "router"
    air = box - cond
    air.faces.name = "cair"
    air.faces.maxh = L / 8
    geo = Glue([air, cond])
    mesh = ng.Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=L / 8))
    mesh.Curve(3)
    fes = ng.H1(mesh, order=2, dirichlet="router", complex=True)
    u, v = fes.TnT()
    sig = mesh.MaterialCF({"cond": SIGMA_CYL}, default=0.0)
    nu = mesh.MaterialCF({"cond": 1.0 / (MU0 * mu_r)}, default=1.0 / MU0)
    As_cf = B0_CYL * (ng.y + 1j * ng.x)
    a = ng.BilinearForm(fes)
    a += nu * ng.grad(u) * ng.grad(v) * ng.dx + 1j * omega * sig * u * v * ng.dx
    a.Assemble()
    f = ng.LinearForm(fes)
    f += -1j * omega * sig * As_cf * v * ng.dx
    if mu_r != 1.0:
        gAs = ng.CoefficientFunction((1j * B0_CYL, B0_CYL + 0j))
        f += -(nu - 1.0 / MU0) * gAs * ng.grad(v) * ng.dx
    f.Assemble()
    gfA = ng.GridFunction(fes)
    gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    phi = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    acc = 0.0
    for p in phi:
        x, y = Rc * np.cos(p), Rc * np.sin(p)
        g = ng.grad(gfA)(mesh(x, y))
        Bx = complex(g[1]) + B0_CYL
        By = -complex(g[0]) - 1j * B0_CYL
        Br = Bx * np.cos(p) + By * np.sin(p)
        Bp = -Bx * np.sin(p) + By * np.cos(p)
        acc += (Br * np.conj(Bp)).real
    return Rc * Rc / (2 * MU0) * (2 * np.pi / nphi) * acc


# =====================================================================================
# Case C: mini CAGE induction machine (VIM iron core + reduced-FEM bar ring)
# =====================================================================================
R_CORE = 0.08
N_BARIM, R_BARIM, A_BARIM = 8, 0.095, 0.010
SIGMA_BAR = 3.5e7
B0_IM = 0.05
H0_IM = B0_IM / MU0
RC_IM = 0.13
COLLAR_IM = (-3.2, -3.2, 3.2, 3.2)
BAR_POS = [(R_BARIM * np.cos(2 * np.pi * k / N_BARIM), R_BARIM * np.sin(2 * np.pi * k / N_BARIM))
           for k in range(N_BARIM)]
H_ROT_RE = ng.CoefficientFunction((H0_IM, 0.0))
H_ROT_IM = ng.CoefficientFunction((0.0, -H0_IM))


def _voxel_complex(func, bbox, nx, ny):
    x0, y0, x1, y1 = bbox
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(y0, y1, ny)
    X, Y = np.meshgrid(xs, ys)
    V = func(np.stack([X.ravel(), Y.ravel()], axis=1)).reshape(ny, nx).astype(complex)
    return ng.VoxelCoefficient((x0, y0), (x1, y1), np.ascontiguousarray(V), linear=True)


def Az_iron_polar_voxel(body, m_re, m_im, W, nx=340, ny=340, Nr=64, Nphi=1440):
    """SINGLE-VALUED iron A_z on a window that SURROUNDS the body, branch-cut-free.

    The direct conjugate-potential formula (PlanarDemagBody.Az_at) has a branch cut along the
    -x ray of every charge; a bar ring around the core crosses them.  Construction: integrate
    dA/dphi = mu0 r H_r (H_at is branch-free) around circles, anchored on the +x axis where the
    direct formula IS cut-free; closure over 2 pi is exact because the total charge is zero
    (Gauss) -- asserted."""
    r_in = R_CORE + 0.002
    r_out = float(np.hypot(W[2], W[3])) + 0.005
    rs = np.linspace(r_in, r_out, Nr)
    ph = np.linspace(0.0, 2 * np.pi, Nphi, endpoint=False)
    R, PH = np.meshgrid(rs, ph, indexing="ij")
    P = np.stack([(R * np.cos(PH)).ravel(), (R * np.sin(PH)).ravel()], axis=1)

    def one(m):
        H = np.empty_like(P)
        for k0 in range(0, len(P), 8000):
            H[k0:k0 + 8000] = body.H_at(P[k0:k0 + 8000], m)
        Hr = (H[:, 0] * np.cos(PH.ravel()) + H[:, 1] * np.sin(PH.ravel())).reshape(Nr, Nphi)
        A0 = body.Az_at(np.stack([rs, np.zeros(Nr)], axis=1), m)
        dphi = 2 * np.pi / Nphi
        mid = 0.5 * (Hr[:, 1:] + Hr[:, :-1])
        cum = np.concatenate([np.zeros((Nr, 1)), np.cumsum(mid, axis=1) * dphi], axis=1)
        A = A0[:, None] + MU0 * rs[:, None] * cum
        closure = MU0 * rs * (cum[:, -1] + 0.5 * (Hr[:, -1] + Hr[:, 0]) * dphi)
        scale = np.abs(A).max() + 1e-300
        if np.abs(closure).max() >= 5e-3 * scale:
            raise RuntimeError("Az polar closure failed: %.2e vs scale %.2e"
                               % (np.abs(closure).max(), scale))
        return A

    A_re, A_im = one(m_re), one(m_im)
    ph_w = np.append(ph, 2 * np.pi)
    itp_re = RegularGridInterpolator((rs, ph_w), np.concatenate([A_re, A_re[:, :1]], axis=1),
                                     bounds_error=False, fill_value=0.0)
    itp_im = RegularGridInterpolator((rs, ph_w), np.concatenate([A_im, A_im[:, :1]], axis=1),
                                     bounds_error=False, fill_value=0.0)
    xs = np.linspace(W[0], W[2], nx)
    ys = np.linspace(W[1], W[3], ny)
    X, Y = np.meshgrid(xs, ys)
    pts = np.stack([np.hypot(X, Y).ravel(), np.mod(np.arctan2(Y, X), 2 * np.pi).ravel()], axis=1)
    V = (itp_re(pts) + 1j * itp_im(pts)).reshape(ny, nx)
    return ng.VoxelCoefficient((W[0], W[1]), (W[2], W[3]), np.ascontiguousarray(V), linear=True)


def _bar_geo(bar_h=0.0025):
    box = WorkPlane().MoveTo(COLLAR_IM[0], COLLAR_IM[1]).Rectangle(
        COLLAR_IM[2] - COLLAR_IM[0], COLLAR_IM[3] - COLLAR_IM[1]).Face()
    box.edges.name = "router"
    bars = []
    for cx, cy in BAR_POS:
        b = WorkPlane().Circle(cx, cy, A_BARIM).Face()
        b.faces.name = "bar"
        b.faces.maxh = bar_h
        bars.append(b)
    air = box
    for b in bars:
        air = air - b
    air.faces.name = "cair"
    return Glue([air] + bars)


def reduced_bars_step(body, m_re, m_im, omega):
    """One complex eddy solve on the bar ring; returns (bar loss, current-quadrature Pk, Ik)."""
    geo = _bar_geo()
    mesh = ng.Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=0.35))
    mesh.Curve(3)
    Wwin = (-(R_BARIM + 1.6 * A_BARIM), -(R_BARIM + 1.6 * A_BARIM),
            (R_BARIM + 1.6 * A_BARIM), (R_BARIM + 1.6 * A_BARIM))
    As_cf = B0_IM * (ng.y + 1j * ng.x) + Az_iron_polar_voxel(body, m_re, m_im, Wwin)
    fes = ng.H1(mesh, order=2, dirichlet="router", complex=True)
    u, v = fes.TnT()
    sig = mesh.MaterialCF({"bar": SIGMA_BAR}, default=0.0)
    a = ng.BilinearForm(fes)
    a += (1.0 / MU0) * ng.grad(u) * ng.grad(v) * ng.dx + 1j * omega * sig * u * v * ng.dx
    a.Assemble()
    f = ng.LinearForm(fes)
    f += -1j * omega * sig * As_cf * v * ng.dx
    f.Assemble()
    gfA = ng.GridFunction(fes)
    gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    Atot = gfA + As_cf
    loss = 0.5 * SIGMA_BAR * omega ** 2 * ng.Integrate(
        Atot * ng.Conj(Atot), mesh, definedon=mesh.Materials("bar")).real
    segs = []
    ir = ng.IntegrationRule(ng.TRIG, 2)
    for el in mesh.Elements(ng.VOL):
        if mesh[el].mat != "bar":
            continue
        tr = mesh.GetTrafo(el)
        for ip in ir:
            mip = tr(ip)
            wq = ip.weight * abs(np.linalg.det(np.asarray(mip.jacobi)))
            segs.append((mip.point[0], mip.point[1],
                         -1j * omega * SIGMA_BAR * complex(Atot(mip)) * wq))
    Pk = np.array([(s[0], s[1]) for s in segs])
    Ik = np.array([s[2] for s in segs])
    return loss, Pk, Ik


def H_of_currents(P, Pk, Ik):
    dx = P[:, 0, None] - Pk[None, :, 0]
    dy = P[:, 1, None] - Pk[None, :, 1]
    r2 = dx * dx + dy * dy
    Hx = (Ik[None, :] / (2 * np.pi) * (-dy / r2)).sum(axis=1)
    Hy = (Ik[None, :] / (2 * np.pi) * (dx / r2)).sum(axis=1)
    return np.stack([Hx, Hy], axis=1)


def _H_voxels_complex(Pk, Ik, bbox, nx=110, ny=110):
    x0, y0, x1, y1 = bbox
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(y0, y1, ny)
    X, Y = np.meshgrid(xs, ys)
    H = H_of_currents(np.stack([X.ravel(), Y.ravel()], axis=1), Pk, Ik)

    def two(part):
        return ng.CoefficientFunction((
            ng.VoxelCoefficient((x0, y0), (x1, y1),
                                np.ascontiguousarray(part(H[:, 0]).reshape(ny, nx)), linear=True),
            ng.VoxelCoefficient((x0, y0), (x1, y1),
                                np.ascontiguousarray(part(H[:, 1]).reshape(ny, nx)), linear=True)))
    return two(np.real), two(np.imag)


def cage_coupled_slip(body, omega, max_stag=25, tol=1e-3, relax=0.7):
    """Stagger { bars <-> iron } until dM < tol (at HIGH slip the bar reaction grows with the
    bar currents, so a fixed stagger count stops converging -- iterate, mildly under-relaxed).
    Returns (T_avg, history)."""
    ib = (-1.05 * R_CORE, -1.05 * R_CORE, 1.05 * R_CORE, 1.05 * R_CORE)
    m_re, chi_e, _, _ = body.solve_nonlinear(M_of_h, chi_sec, body.project(H_ROT_RE))
    m_im = np.linalg.solve(body.weighted_mass(1.0 / chi_e) + body.N,
                           body.Md @ body.project(H_ROT_IM))
    hist = []
    Pk = Ik = None
    for it in range(max_stag):
        loss, Pk, Ik = reduced_bars_step(body, m_re, m_im, omega)
        Hre_cf, Him_cf = _H_voxels_complex(Pk, Ik, ib)
        m_re_new, chi_e, itn, _ = body.solve_nonlinear(
            M_of_h, chi_sec, body.project(H_ROT_RE + Hre_cf))
        m_im_new = np.linalg.solve(body.weighted_mass(1.0 / chi_e) + body.N,
                                   body.Md @ body.project(H_ROT_IM + Him_cf))
        dM = (np.linalg.norm(m_re_new - m_re) + np.linalg.norm(m_im_new - m_im)) \
            / max(np.linalg.norm(m_re_new), 1e-300)
        m_re = (1 - relax) * m_re + relax * m_re_new
        m_im = (1 - relax) * m_im + relax * m_im_new
        hist.append({"stagger": it, "loss": loss, "dM": dM, "picard": itn})
        if dM < tol:
            break
    else:
        raise RuntimeError(f"cage stagger NOT converged at slip {omega}: dM={dM:.2e}")
    nphi = 1440
    phi = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    P = np.stack([RC_IM * np.cos(phi), RC_IM * np.sin(phi)], axis=1)
    Hi = body.H_at(P, m_re) + 1j * body.H_at(P, m_im)
    Hb = H_of_currents(P, Pk, Ik)
    Bx = MU0 * (Hi[:, 0] + Hb[:, 0]) + B0_IM
    By = MU0 * (Hi[:, 1] + Hb[:, 1]) - 1j * B0_IM
    Br = Bx * np.cos(phi) + By * np.sin(phi)
    Bp = -Bx * np.sin(phi) + By * np.cos(phi)
    T = RC_IM * RC_IM / (2 * MU0) * (2 * np.pi / nphi) * float((Br * np.conj(Bp)).real.sum())
    return T, hist


def cage_allinone_slip(omega):
    """All-in-one reference: box 32, complex Dirichlet A = B0 (y + j x), frozen-secant iron nu
    from its own static real-part Newton, bars sigma.  Returns (T_avg, bar loss, ndof)."""
    box = WorkPlane().RectangleC(32, 32).Face()
    box.edges.name = "outer"
    core = WorkPlane().Circle(0, 0, R_CORE).Face()
    core.faces.name = "iron"
    core.faces.maxh = 0.012
    bars = []
    for cx, cy in BAR_POS:
        b = WorkPlane().Circle(cx, cy, A_BARIM).Face()
        b.faces.name = "bar"
        b.faces.maxh = 0.0025
        bars.append(b)
    ring = WorkPlane().Circle(0, 0, 0.6).Face()
    airin = ring - core
    for b in bars:
        airin = airin - b
    airin.faces.name = "airin"
    airin.faces.maxh = 0.03
    airout = box - ring
    airout.faces.name = "airout"
    geo = Glue([airout, airin, core] + bars)
    mesh = ng.Mesh(OCCGeometry(geo, dim=2).GenerateMesh(maxh=2.0))
    mesh.Curve(3)
    # (1) static real Newton (drive = the real part of the rotating field)
    fesR = ng.H1(mesh, order=2, dirichlet="outer")
    uR, vR = fesR.TnT()
    nu_rot, _ = _nu_rotor_cf(ng.grad(uR))
    ind = mesh.MaterialCF({"iron": 1.0}, default=0.0)
    nu_tot = ind * nu_rot + (1.0 - ind) * (1.0 / MU0)
    aR = ng.BilinearForm(fesR)
    aR += (nu_tot * ng.grad(uR) * ng.grad(vR)) * ng.dx
    gfA0 = ng.GridFunction(fesR)
    gfA0.Set(B0_IM * ng.y, definedon=mesh.Boundaries("outer"))
    ret = ngsolvers.Newton(aR, gfA0, maxit=100, maxerr=1e-9, dampfactor=0.6,
                           printing=False, inverse="sparsecholesky")
    status = ret[0] if isinstance(ret, tuple) else ret
    if status != 0:
        raise RuntimeError("cage all-in-one static Newton not converged")
    # (2) freeze per-element secant nu
    fes0 = ng.L2(mesh, order=0)
    gfnu = ng.GridFunction(fes0)
    areas = ng.Integrate(ng.CoefficientFunction(1.0), mesh, element_wise=True)
    areas = np.array([areas[k] for k in range(mesh.ne)])
    el2dof = np.zeros(mesh.ne, dtype=int)
    is_iron = np.zeros(mesh.ne, dtype=bool)
    for el in mesh.Elements(ng.VOL):
        el2dof[el.nr] = fes0.GetDofNrs(el)[0]
        is_iron[el.nr] = mesh[el].mat == "iron"
    H_TAB = np.concatenate([[0.0], np.logspace(0, 7.5, 400)])
    B_TAB = MU0 * (H_TAB + CHI0 * H_TAB / (1.0 + CHI0 * H_TAB / MSAT))
    NU_TAB = np.empty_like(H_TAB)
    NU_TAB[1:] = H_TAB[1:] / B_TAB[1:]
    NU_TAB[0] = 1.0 / (MU0 * (1.0 + CHI0))
    Bint = ng.Integrate(ng.Norm(ng.grad(gfA0)), mesh, element_wise=True)
    Bel = np.array([Bint[k] for k in range(mesh.ne)]) / areas
    nu_e = np.full(mesh.ne, 1.0 / MU0)
    nu_e[is_iron] = np.interp(Bel[is_iron], B_TAB, NU_TAB)
    gfnu.vec.FV().NumPy()[el2dof] = nu_e
    # (3) one complex solve with complex Dirichlet A = B0 (y + j x)
    fes = ng.H1(mesh, order=2, dirichlet="outer", complex=True)
    u, v = fes.TnT()
    sig = mesh.MaterialCF({"bar": SIGMA_BAR}, default=0.0)
    a = ng.BilinearForm(fes)
    a += gfnu * ng.grad(u) * ng.grad(v) * ng.dx + 1j * omega * sig * u * v * ng.dx
    a.Assemble()
    gfA = ng.GridFunction(fes)
    gfA.Set(B0_IM * (ng.y + 1j * ng.x), definedon=mesh.Boundaries("outer"))
    r = (-1.0) * (a.mat * gfA.vec)
    upd = gfA.vec.CreateVector()
    upd.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r
    gfA.vec.data += upd
    nphi = 1440
    acc = 0.0
    for p in np.linspace(0, 2 * np.pi, nphi, endpoint=False):
        x, y = RC_IM * np.cos(p), RC_IM * np.sin(p)
        g = ng.grad(gfA)(mesh(x, y))
        Bx, By = complex(g[1]), -complex(g[0])
        Br = Bx * np.cos(p) + By * np.sin(p)
        Bp = -Bx * np.sin(p) + By * np.cos(p)
        acc += (Br * np.conj(Bp)).real
    T = RC_IM * RC_IM / (2 * MU0) * (2 * np.pi / nphi) * acc
    loss = 0.5 * SIGMA_BAR * omega ** 2 * ng.Integrate(
        gfA * ng.Conj(gfA), mesh, definedon=mesh.Materials("bar")).real
    return T, loss, fes.ndof
