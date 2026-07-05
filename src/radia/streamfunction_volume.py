"""Volume (3D) stream-function coil design -- foliated current-Clebsch.

A 3D VOLUME stream function on a solid conductor: fix a foliation ``mu`` (here the
nested cylinders ``mu = r = sqrt(x^2 + y^2)``) and put a single scalar ``lambda``
on the volume.  The current is

    J = grad(lambda) x grad(mu)        (div J = 0 identically, sympy-verified)

With ``mu = r`` this is the azimuthal current ``J = (d lambda/dz) phi_hat``, so
each ``mu = const`` cylinder carries a surface stream function ``lambda`` -- a
STACK of surface stream functions = a true 3D bulk winding.  Wires are the
level-set intersections ``{lambda = const} INT {mu = const}`` = ``lambda``
contours on each cylinder, each carrying ``I = d_lambda * d_mu`` (equal current
by construction).

Because ``mu`` is FIXED the inverse design ``A lambda = B_target`` is a LINEAR
least-norm problem, solved by the existing kernel-agnostic
``radia.stream_function`` ACA+TSVD engine UNCHANGED -- only the per-DOF basis
current is the volume ``grad(phi_j) x grad(mu)`` instead of a surface
``n x grad``.

This module is the SHIPPED home of the pipeline validated by
``validation_test/feec/vim_legacy/foliated_clebsch_solenoid.py`` (Stage A) and
``validation_test/feec/vim_legacy/foliated_solenoid_wires.py`` (Stage B).  The headless panel backend
``radia.panels.calc_streamfunction_volume`` wraps ``design_volume_coil`` here.

Domain of validity (the honest frontier; see
``docs/stream_function/examples_catalog.ipynb`` and the VIM retirement ledger): clean equal-current
extraction needs the zero-helicity, fixed-(radial)-foliation, closed-streamline
regime.  A general target / arbitrary foliation hits the F1 (non-convex
foliation choice), F2 (ill-posed Clebsch recovery) or F3 (nonzero-helicity
obstruction) walls -- not handled here by design.

Caller wraps NGSolve work in ``with TaskManager():`` (Caller-Wraps policy);
this module never opens its own TaskManager.
"""
import math

import numpy as np

from ngsolve import (
    Mesh, H1, HCurl, GridFunction, CF, grad, curl, InnerProduct, dx, sqrt,
    x, y, z, LinearForm, Integrate, Cross,
)

from radia.stream_function import aca_tsvd, pseudo_inverse_solve

MU0 = 4e-7 * math.pi
PI = math.pi

__all__ = [
    "tube_extent", "radial_grad_mu_cf", "assemble_response", "solve_lambda",
    "helicity_relative", "extract_wires", "bz_segments", "bz_radia",
    "write_wires_msh", "design_volume_coil", "MU0",
]


# ----------------------------------------------------------------------
# geometry probe (derive the cylindrical tube extent from the loaded mesh)
# ----------------------------------------------------------------------
def tube_extent(mesh):
    """Return (r0, r1, zl): inner radius, outer radius, half-length [m].

    Derived from the mesh vertices -- the radial foliation ``mu = r`` only makes
    physical sense for an annular (hollow cylindrical) conductor about the z
    axis.  Raises if the conductor reaches the axis (r0 ~ 0): radial foliation
    is undefined there, and silently proceeding would produce wrong wires (No
    Fallbacks policy)."""
    r_min, r_max, z_max = math.inf, 0.0, 0.0
    for v in mesh.vertices:
        px, py, pz = v.point
        rr = math.hypot(px, py)
        r_min = min(r_min, rr)
        r_max = max(r_max, rr)
        z_max = max(z_max, abs(pz))
    if not (r_max > 0 and z_max > 0):
        raise ValueError("degenerate mesh extent (r_max or z_max is 0).")
    if r_min < 0.02 * r_max:
        raise ValueError(
            "radial foliation needs a HOLLOW cylindrical tube about the z axis; "
            f"the conductor reaches the axis (r_min={r_min:.4g} << "
            f"r_max={r_max:.4g}).  Mesh an annular winding region, or use a "
            "different foliation.")
    return r_min, r_max, z_max


def radial_grad_mu_cf():
    """grad(mu) for the radial foliation mu = r = sqrt(x^2 + y^2)."""
    r = sqrt(x * x + y * y)
    return CF((x / r, y / r, 0.0))


# ----------------------------------------------------------------------
# Stage A: linear inverse design (ACA+TSVD) for the foliated lambda
# ----------------------------------------------------------------------
def assemble_response(mesh, fes, targets):
    """Dense M x N Biot-Savart response A[i,j] = Bz at target i from basis j.

    Since grad(mu)_z = 0,
        [(grad v x grad mu) x K]_z = -(dv/dz)(grad mu . K),  K = (p_i-x)/|p_i-x|^3
    so each response row is a single LinearForm in (dv/dz)."""
    gmu = radial_grad_mu_cf()
    _, v = fes.TnT()
    X = CF((x, y, z))
    M, N = len(targets), fes.ndof
    A = np.zeros((M, N))
    for i, (px, py, pz) in enumerate(targets):
        P = CF((px, py, pz))
        d = P - X
        r3 = (InnerProduct(d, d)) ** 1.5
        K = d / r3
        lf = LinearForm(fes)
        lf += -(MU0 / (4 * PI)) * grad(v)[2] * InnerProduct(gmu, K) * dx
        lf.Assemble()
        A[i, :] = lf.vec.FV().NumPy()
    return A


def solve_lambda(mesh, fes, targets, b0, aca_eps=1e-10):
    """Least-norm foliated lambda for uniform axial Bz=b0 at the targets.

    Returns (gf_lambda, k_aca, fit_res)."""
    M, N = len(targets), fes.ndof
    A = assemble_response(mesh, fes, targets)
    res = aca_tsvd(M, N, lambda i, j: A[i, j],
                   modes=M, kmax=M, aca_eps=aca_eps)
    lam = pseudo_inverse_solve(res, np.full(M, b0))
    fit_res = float(np.linalg.norm(A @ lam - b0) / (abs(b0) * math.sqrt(M)))
    gf = GridFunction(fes)
    gf.vec.FV().NumPy()[:] = lam
    return gf, int(res.k_aca), fit_res


def helicity_relative(mesh, T_cf, order=3):
    """H_rel = |int T . curl T| / (||T|| ||curl T||) in [0,1] for a CF T.

    ~0 => Clebsch-representable (clean level-set wire extraction possible);
    ~1 => maximally helical (linked current lines, no global Clebsch)."""
    gfT = GridFunction(HCurl(mesh, order=order))
    gfT.Set(T_cf)
    J = curl(gfT)
    H = Integrate(InnerProduct(gfT, J) * dx, mesh)
    tn = math.sqrt(Integrate(InnerProduct(gfT, gfT) * dx, mesh))
    jn = math.sqrt(Integrate(InnerProduct(J, J) * dx, mesh))
    return abs(H) / (tn * jn + 1e-30)


# ----------------------------------------------------------------------
# Stage B: equal-current wire extraction + two-codebase Biot-Savart check
# ----------------------------------------------------------------------
def _sample_cylinder(mesh, gf, r_k, zl, nphi, nz):
    phis = np.linspace(0, 2 * PI, nphi, endpoint=True)
    zs = np.linspace(-0.95 * zl, 0.95 * zl, nz)
    LAM = np.zeros((nz, nphi))
    for j, ph in enumerate(phis):
        cx, sy = r_k * math.cos(ph), r_k * math.sin(ph)
        for i, zz in enumerate(zs):
            LAM[i, j] = gf(mesh(cx, sy, zz))
    return phis, zs, LAM


def extract_wires(mesh, gf, r0, r1, zl, n_leaves, rings_per_span,
                  nphi=49, nz=41):
    """Contour lambda on each mu-leaf cylinder at GLOBAL equal-Delta-lambda.

    Returns (wires, currents): wires = list of (n,3) polylines, currents =
    list of I = d_lambda * d_mu (equal-current by construction)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    K = int(n_leaves)
    radii = np.linspace(r0 + (r1 - r0) / (2 * K), r1 - (r1 - r0) / (2 * K), K)
    dmu = (r1 - r0) / K
    lam_vals = gf.vec.FV().NumPy()
    lam_span = float(lam_vals.max() - lam_vals.min())
    dlam = lam_span / float(rings_per_span)

    wires, currents = [], []
    lo, hi = +1e30, -1e30
    samples = []
    for r_k in radii:
        phis, zs, LAM = _sample_cylinder(mesh, gf, float(r_k), zl, nphi, nz)
        samples.append((float(r_k), phis, zs, LAM))
        lo, hi = min(lo, LAM.min()), max(hi, LAM.max())
    levels = np.arange(math.ceil(lo / dlam) * dlam, hi, dlam)
    for (r_k, phis, zs, LAM) in samples:
        PHI, Z = np.meshgrid(phis, zs)
        cs = plt.contour(PHI, Z, LAM, levels=levels)
        for segs in cs.allsegs:
            for seg in segs:
                if len(seg) < 2:
                    continue
                ph, zz = seg[:, 0], seg[:, 1]
                pts = np.column_stack([r_k * np.cos(ph), r_k * np.sin(ph), zz])
                if (ph.max() - ph.min()) > 1.8 * PI:        # close a full ring
                    pts = np.vstack([pts, pts[0]])
                wires.append(pts)
                currents.append(dlam * dmu)
        plt.close("all")
    return wires, currents


def bz_segments(wires, currents, obs):
    """Bz at obs from straight current segments (finite-segment Biot-Savart).

    B = (mu0 I/4pi)(a x b)(|a|+|b|)/(|a||b|(|a||b|+a.b)),  a=P0-P1, b=P0-P2."""
    out = []
    for p in obs:
        P0 = np.asarray(p, dtype=float)
        Bz = 0.0
        for W, I in zip(wires, currents):
            P1, P2 = W[:-1], W[1:]
            a, b = P0 - P1, P0 - P2
            am = np.linalg.norm(a, axis=1)
            bm = np.linalg.norm(b, axis=1)
            cr = np.cross(a, b)
            denom = am * bm * (am * bm + np.einsum("ij,ij->i", a, b))
            good = denom > 1e-30
            f = np.zeros_like(am)
            f[good] = (am[good] + bm[good]) / denom[good]
            Bz += (MU0 * I / (4 * PI)) * np.sum(f * cr[:, 2])
        out.append(Bz)
    return np.array(out)


def bz_radia(wires, currents, obs):
    """Independent cross-check: Bz from Radia rad.ObjFlmCur + rad.Fld."""
    import radia as rad
    rad.UtiDelAll()
    objs = []
    for W, I in zip(wires, currents):
        pts = [[float(px), float(py), float(pz)] for (px, py, pz) in W]
        objs.append(rad.ObjFlmCur(pts, float(I)))
    cnt = rad.ObjCnt(objs)
    out = [float(rad.Fld(cnt, "b", list(p))[2]) for p in obs]
    rad.UtiDelAll()
    return np.array(out)


# ----------------------------------------------------------------------
# wire -> GMSH .msh v4.1 (1D line elements; lab-wide v4.1-only policy)
# ----------------------------------------------------------------------
def write_wires_msh(wires, path):
    """Write the extracted wires as 1D line elements to a GMSH .msh v4.1 file.

    One curve entity holding all wire polylines (a viewer overlay of the coil).
    """
    nodes = []          # (x, y, z)
    lines = []          # (n1, n2) 1-based node tags
    bb = [math.inf, math.inf, math.inf, -math.inf, -math.inf, -math.inf]
    for W in wires:
        base = len(nodes)
        for (px, py, pz) in W:
            nodes.append((float(px), float(py), float(pz)))
            bb[0], bb[1], bb[2] = min(bb[0], px), min(bb[1], py), min(bb[2], pz)
            bb[3], bb[4], bb[5] = max(bb[3], px), max(bb[4], py), max(bb[5], pz)
        for k in range(len(W) - 1):
            lines.append((base + k + 1, base + k + 2))   # 1-based
    nN, nE = len(nodes), len(lines)
    if nN == 0:
        bb = [0, 0, 0, 0, 0, 0]
    out = []
    out.append("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
    out.append("$Entities\n0 1 0 0\n")
    out.append(f"1 {bb[0]:.9g} {bb[1]:.9g} {bb[2]:.9g} "
               f"{bb[3]:.9g} {bb[4]:.9g} {bb[5]:.9g} 0 0\n")
    out.append("$EndEntities\n")
    out.append("$Nodes\n")
    out.append(f"1 {nN} 1 {nN}\n")
    out.append(f"1 1 0 {nN}\n")                 # dim=1 tag=1 parametric=0 count
    for i in range(nN):
        out.append(f"{i + 1}\n")
    for (px, py, pz) in nodes:
        out.append(f"{px:.9g} {py:.9g} {pz:.9g}\n")
    out.append("$EndNodes\n")
    out.append("$Elements\n")
    out.append(f"1 {nE} 1 {nE}\n")
    out.append(f"1 1 1 {nE}\n")                 # dim=1 tag=1 type=1(line) count
    for e, (n1, n2) in enumerate(lines):
        out.append(f"{e + 1} {n1} {n2}\n")
    out.append("$EndElements\n")
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("".join(out))
    return nN, nE


# ----------------------------------------------------------------------
# orchestrator (the panel backend entry point)
# ----------------------------------------------------------------------
def design_volume_coil(mesh, target_bz=1.0e-3, n_leaves=3, fes_order=2,
                       aca_eps=1e-10, rings_per_span=14, n_targets=9,
                       nphi=49, nz=41):
    """Foliated-Clebsch volume stream-function solenoid design.

    Loads geometry from ``mesh`` (a hollow cylindrical conductor), solves the
    foliated lambda for a uniform axial Bz=target_bz, extracts equal-current
    wires, and cross-checks the wire Biot-Savart against Radia.

    Returns a dict of design metrics + the wire geometry.  Caller must be inside
    ``with TaskManager():``."""
    r0, r1, zl = tube_extent(mesh)
    fes = H1(mesh, order=int(fes_order))
    zt = np.linspace(-0.6 * zl, 0.6 * zl, int(n_targets))
    targets = [(0.0, 0.0, float(zz)) for zz in zt]

    gf, k_aca, fit_res = solve_lambda(mesh, fes, targets, target_bz, aca_eps)

    # helicity gate (radial foliation T = lambda*grad(mu) is Clebsch -> ~0)
    T_cf = gf * radial_grad_mu_cf()
    h_rel = helicity_relative(mesh, T_cf)

    wires, currents = extract_wires(mesh, gf, r0, r1, zl, n_leaves,
                                    rings_per_span, nphi, nz)
    bz_wire = bz_segments(wires, currents, targets)
    bz_rad = bz_radia(wires, currents, targets)

    field_res = float(np.linalg.norm(bz_wire - target_bz)
                      / (abs(target_bz) * math.sqrt(len(targets))))
    twocode = float(np.linalg.norm(bz_wire - bz_rad)
                    / (np.linalg.norm(bz_wire) + 1e-30))
    I_arr = np.array(currents) if currents else np.array([0.0])
    equal_cur = float((I_arr.max() - I_arr.min()) / (I_arr.mean() + 1e-30))

    return {
        "r0_m": float(r0), "r1_m": float(r1), "zl_m": float(zl),
        "ne": int(mesh.ne), "ndof": int(fes.ndof),
        "n_targets": int(len(targets)), "k_aca": int(k_aca),
        "fit_res": fit_res, "h_rel": float(h_rel),
        "n_wires": int(len(wires)),
        "current_per_wire_A": float(currents[0]) if currents else 0.0,
        "field_res": field_res, "two_codebase": twocode,
        "equal_current_spread": equal_cur,
        "mean_bz_T": float(bz_wire.mean()) if len(bz_wire) else 0.0,
        "_wires": wires, "_currents": currents,
    }
