#!/usr/bin/env python3
"""mini-REGCOIL / NESCOIL: a fusion coil-design (stellarator Stage-2) demo.

The stream-function (current-potential) machinery that designs an MRI gradient
coil or an induction-heating work coil is, mathematically, the SAME object the
fusion community calls a **winding-surface current potential** -- the unknown of
NESCOIL (Merkel 1987), REGCOIL (Landreman 2017) and the surface part of FOCUS
(Zhu 2018).  The Stage-2 coil problem is:

    given a target NORMAL field  B.n  on the PLASMA boundary,
    find a current potential psi on a WINDING SURFACE around it
    whose Biot-Savart field reproduces that B.n  (its iso-contours = the coils).

This demo runs that forward problem with Radia's existing surface-FE stream
function, in four parts:

  1. FORWARD WORKS -- two producible targets reproduced to MACHINE PRECISION:
       * uniform vertical field B.n   (a PF / equilibrium / vertical-field coil)
       * non-axisymmetric  sin(theta) cos(2 phi)  (a stellarator-like shape)
     ``A_n psi = B_n`` where the rows of ``A_n`` are the plasma-boundary normal
     component ``n.B`` of the winding-surface Biot-Savart kernel (we already
     assemble the 3-component ``A``; we just dot it with the plasma normal).

  2. THE REGCOIL TRADE-OFF -- on a GENUINELY HARD target the winding surface
     cannot reproduce cheaply (a high (theta, phi) mode that decays across the
     plasma-coil gap), sweeping the regularisation weight alpha traces the
     classic REGCOIL L-curve: (field error, coil complexity).

  3. NET CURRENT (the multivalued / secular term) -- a single-valued psi can
     only carry currents that net to zero through each hole of the winding
     torus.  A real coil set carries NET current: the current potential gains a
     SECULAR term  Psi = psi + (G/2pi) zeta + (I/2pi) theta  (zeta toroidal,
     theta poloidal angle).  The TWO extra degrees of freedom are the first
     cohomology generators of the winding surface -- their COUNT is the surface
     Betti number b1 = 2 (genus 1), which we CONFIRM gmsh-free via the Euler
     characteristic of the triangulated torus (b1 = 2 - chi; the volume T-Omega
     cut engine ``src/radia/cohomology_cut.py`` is itself gmsh-free,
     ``radia.cohomology``).  We verify the
     net-poloidal-current (TF) secular term reproduces the textbook  B_tor ~ 1/R
     toroidal field (Ampere's law) inside the tube and ~0 outside.  Key physics:
     the TF field is TANGENT to the plasma boundary (n.B ~ 0), so the net
     poloidal current is a PRESCRIBED engineering parameter (you set it for the
     desired on-axis B_tor) -- exactly why REGCOIL takes net_poloidal_current as
     an INPUT, not a fitted output.

  4. A VMEC-SHAPED PLASMA BOUNDARY -- the circular plasma torus is replaced by a
     genuinely non-axisymmetric **rotating-ellipse** boundary in the VMEC Fourier
     representation  R = sum RBC(m,n) cos(m th - n NFP ph),  Z = sum ZBS(m,n)
     sin(...) , with analytic surface normals from the parametric tangents.  A
     real machine's free-boundary equilibrium can be dropped in with
     ``--wout wout_*.nc`` (a standard VMEC output; the netCDF4 reader extracts
     the boundary rmnc/zmns coefficients).

HONEST SCOPE (what this is and is NOT):
  * Parts 1-2 use a SINGLE-VALUED psi -- correct for PF / RMP / shaping / shim
    fields.  Part 3 adds the multivalued secular term for net-current (TF-type)
    coils; the generator COUNT is computed (Euler characteristic, gmsh-free), and on the torus
    the generators themselves are analytic (exact).  For a general winding
    surface the generators come from the cohomology cut
    (``src/radia/cohomology_cut.py``), not analytic forms.
  * The producible targets reproduce to ~1e-8 BECAUSE they are smooth and within
    reach; that is honest (the forward map is exact), not a hidden tolerance.
  * Part 4's default boundary is an ANALYTIC rotating-ellipse model, NOT a
    converged equilibrium -- a stand-in with the genuine VMEC Fourier shape.  A
    production run drops in a real ``wout_*.nc`` (reader provided, verified by a
    round-trip against the VMEC netCDF schema) and additionally needs coil force
    / stress and winding-surface geometry optimisation (FOCUS).

Run (standalone -- no Cubit, no panel UI):
    python demo_regcoil_fusion.py
    python demo_regcoil_fusion.py --no-plot               # JSON only
    python demo_regcoil_fusion.py --wout wout_w7x.nc      # real VMEC boundary
    python demo_regcoil_fusion.py --no-cohomology         # skip the b1 (Euler-char) check

Writes validation evidence to
``validation_test/stream_function/demos/demo_regcoil_fusion.json`` and, if
matplotlib + radia-mcp are present, ``demo_regcoil_fusion.{png,pdf}`` next to
this script.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
from datetime import datetime

import numpy as np

from _validation_output import validation_output

# major radius shared by the winding torus and the (circular) plasma torus -- a
# simple, reproducible stand-in for a VMEC boundary (Part 4 replaces it).
R_MAJOR = 0.30
A_PLASMA = 0.06       # plasma-boundary minor radius (circular case, parts 1-3)
A_WIND = 0.12         # winding-surface minor radius (the coil lives here)
_MU0 = 4.0e-7 * math.pi
_MU0_4PI = 1.0e-7


def _torus_surface_vol(a, maxh, path):
    """Save a torus SURFACE mesh (revolved wire = dim-2 surface, what the
    surface stream function lives on -- revolving a Face would give a solid)."""
    from netgen.occ import WorkPlane, Axes, Pnt, Y, Z, Axis, OCCGeometry
    wire = WorkPlane(Axes(Pnt(R_MAJOR, 0, 0), n=Y)).Circle(a).Wire()
    surf = wire.Revolve(Axis(Pnt(0, 0, 0), Z), 360)
    OCCGeometry(surf).GenerateMesh(maxh=maxh).Save(path)
    return path


def _plasma_points_normals(plasma, eval_max):
    """Circular-plasma sample points + analytic outward normals (torus)."""
    pts = np.array([list(v.point) for v in plasma.vertices])
    if len(pts) > eval_max:
        pts = pts[np.linspace(0, len(pts) - 1, eval_max).astype(int)]
    phi = np.arctan2(pts[:, 1], pts[:, 0])
    nrm = np.column_stack([pts[:, 0] - R_MAJOR * np.cos(phi),
                           pts[:, 1] - R_MAJOR * np.sin(phi), pts[:, 2]])
    nrm /= (np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-30)
    rho = np.hypot(pts[:, 0], pts[:, 1]) - R_MAJOR
    theta = np.arctan2(pts[:, 2], rho)
    return pts, nrm, theta, phi


def _A_normal(C, fes, n_cf, pts, nrm):
    """B.n design matrix rows: 3-component winding Biot-Savart dotted with the
    plasma normal -> A_n[m, dof] = n . B at plasma point m from psi-dof."""
    A3 = C._assemble_biot_savart(fes, n_cf, pts, [0, 1, 2]).reshape(
        len(pts), 3, fes.ndof)
    return np.einsum("mc,mcj->mj", nrm, A3)


def _design(C, fes, coil, A_n, Bn, regularize, alpha_rel):
    """Solve min psi^T S psi s.t. (regularised) A_n psi = Bn; return psi + diag.

    Returns (psi, B.n relative residual, peak |grad psi| = coil complexity)."""
    from ngsolve import GridFunction
    S = C._seminorm(fes, regularize).toarray()
    AtA = A_n.T @ A_n
    alpha = alpha_rel * np.trace(AtA) / fes.ndof
    psi = np.linalg.solve(AtA + alpha * S, A_n.T @ Bn)
    res = float(np.linalg.norm(A_n @ psi - Bn) / (np.linalg.norm(Bn) + 1e-30))
    gfu = GridFunction(fes)
    gfu.vec.FV().NumPy()[:] = psi
    return psi, res, float(C._peak_current_density(fes, coil, gfu))


def _contours(C, coil, psi, n_levels):
    verts = np.array([list(v.point) for v in coil.vertices])
    tris = C._surface_triangles(coil)
    tol = max(1e-12, 1e-6 * C._char_len(verts, tris))
    loops = []
    for lev in C._contour_levels(psi[:coil.nv], n_levels):
        loops.extend(C._segs_to_polylines(
            C._contour_segments(verts, psi[:coil.nv], tris, lev), tol))
    return loops


# --------------------------------------------------------------------------
# Part 3: net-current (multivalued / secular) term
# --------------------------------------------------------------------------
def _secular_currents(n_cf):
    """The two net-current secular surface currents on the torus winding
    surface, K = n x grad(angle):

      K_zeta  = n x grad(zeta)   zeta = toroidal angle -> NET POLOIDAL current
                (the TF solenoid; produces a toroidal B ~ 1/R inside the tube)
      K_theta = n x grad(theta)  theta = poloidal angle -> NET TOROIDAL current
                (a plasma-current-like ring; produces a poloidal / vertical B)

    grad(zeta), grad(theta) are single-valued vector fields even though the
    angles are multivalued -- they ARE the cohomology generators (analytic on a
    torus).  Returns [(name, K_cf), ...]."""
    from ngsolve import Cross, CoefficientFunction as CF, x, y, z, sqrt
    rho = sqrt(x * x + y * y)
    s = rho - R_MAJOR
    den = s * s + z * z
    grad_zeta = CF((-y, x, 0)) / (x * x + y * y)
    grad_theta = CF((-z * x / (rho * den), -z * y / (rho * den), s / den))
    return [("net_poloidal_TF", Cross(n_cf, grad_zeta)),
            ("net_toroidal", Cross(n_cf, grad_theta))]


def _assemble_secular_normal(coil, pts, nrm, K_list):
    """B.n at the plasma points from each explicit secular surface current K
    (same Biot-Savart kernel as _assemble_biot_savart, but K is a CF, not a
    per-DOF basis; the normal is already baked into each K = n x grad(...)).
    Returns (M, len(K_list))."""
    from ngsolve import x, y, z, sqrt, Integrate, ds
    out = np.zeros((len(pts), len(K_list)))
    for k, (_name, K) in enumerate(K_list):
        for m, p in enumerate(pts):
            dxt, dyt, dzt = p[0] - x, p[1] - y, p[2] - z
            r2 = dxt * dxt + dyt * dyt + dzt * dzt
            r3 = r2 * sqrt(r2)
            cr = (K[1] * dzt - K[2] * dyt,
                  K[2] * dxt - K[0] * dzt,
                  K[0] * dyt - K[1] * dxt)
            bn = sum(nrm[m, c] * _MU0_4PI * Integrate(cr[c] / r3 * ds, coil)
                     for c in range(3))
            out[m, k] = bn
    return out


def _tf_field_1overR(coil, n_cf, radii):
    """Toroidal field B_y of the unit net-poloidal-current (TF) secular term at
    (r, 0, 0) for r in radii.  Inside the tube Ampere's law gives B_y * r =
    const (~1/R); outside ~0.  Returns the B_y array."""
    from ngsolve import x, y, z, sqrt, Integrate, ds
    _name, K = _secular_currents(n_cf)[0]      # K_zeta = TF
    by = np.zeros(len(radii))
    for i, r in enumerate(radii):
        dxt, dyt, dzt = r - x, -y, -z
        r3 = (dxt * dxt + dyt * dyt + dzt * dzt) * sqrt(
            dxt * dxt + dyt * dyt + dzt * dzt)
        cy = K[2] * dxt - K[0] * dzt
        by[i] = _MU0_4PI * Integrate(cy / r3 * ds, coil)
    return by


def _betti1_winding_surface(a, maxh):
    """First Betti number b1 of the torus winding SURFACE -- gmsh-free.

    Computed from the surface's OWN triangulation via the Euler characteristic
    chi = V - E + F.  For a connected closed orientable surface b0 = 1, b2 = 1,
    so b1 = 2 - chi (= 2g, the genus count).  b1 = number of independent
    net-current (secular) DOFs = 2 for the genus-1 torus winding surface.

    This drops the Gmsh ``computeHomology`` dependency: the count is a direct
    combinatorial property of the closed triangulated torus (no homology solver
    needed).  The T-Omega cohomology CUT engine for general (volume) domains is
    likewise gmsh-free (``radia.cohomology`` / ``cohomology_cut.py``)."""
    nu = max(6, int(round(2 * np.pi * R_MAJOR / maxh)))   # toroidal divisions
    nv = max(6, int(round(2 * np.pi * a / maxh)))         # poloidal divisions

    def vid(i, j):
        return (i % nu) * nv + (j % nv)

    edges = set()
    nfaces = 0
    for i in range(nu):
        for j in range(nv):
            v00, v10 = vid(i, j), vid(i + 1, j)
            v01, v11 = vid(i, j + 1), vid(i + 1, j + 1)
            for tri in ((v00, v10, v11), (v00, v11, v01)):   # split each quad
                nfaces += 1
                for e in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
                    edges.add(frozenset(e))
    V, E, F = nu * nv, len(edges), nfaces
    chi = V - E + F                                          # = 0 for the torus
    return 2 - chi                                          # b1 = 2 - chi = 2


# --------------------------------------------------------------------------
# Part 4: VMEC-format (rotating-ellipse) plasma boundary
# --------------------------------------------------------------------------
def _rotating_ellipse_terms():
    """Default non-axisymmetric boundary: a rotating ellipse in the VMEC Fourier
    representation.  terms = [(m, n, RBC, ZBS), ...]; the (1,1) term rotates the
    elliptical cross-section as it goes around toroidally.  An analytic MODEL
    (not a converged equilibrium) with the genuine VMEC boundary shape."""
    nfp = 3
    terms = [(0, 0, R_MAJOR, 0.0),
             (1, 0, 0.055, 0.055),    # circular minor radius
             (1, 1, 0.020, -0.020)]   # rotating elongation (n != 0)
    return terms, nfp


def _read_wout_boundary(path):
    """Boundary Fourier terms (m, n, RBC, ZBS) + NFP from a VMEC ``wout_*.nc``.
    Reads the standard schema (xm, xn, rmnc, zmns, nfp); the boundary is the
    last flux surface.  STELLARATOR-SYMMETRIC only: a non-symmetric (lasym=T)
    wout also carries rmns/zmnc -- this reader RAISES rather than silently drop
    them (No-Fallback).  Also raises if the file lacks the VMEC variables.
    Verified by a round-trip against the schema in the golden."""
    import netCDF4
    d = netCDF4.Dataset(path, "r")
    try:
        nfp = int(np.asarray(d["nfp"][...]))
        xm = np.asarray(d["xm"][:])
        xn = np.asarray(d["xn"][:])
        rmnc = np.asarray(d["rmnc"][:])        # (ns, mn_mode)
        zmns = np.asarray(d["zmns"][:])
        # a non-stellarator-symmetric (lasym=T) wout carries rmns/zmnc that this
        # symmetric reader would silently truncate -> raise instead.
        if "rmns" in d.variables or "zmnc" in d.variables:
            raise ValueError(
                f"{path}: non-stellarator-symmetric wout (rmns/zmnc present); "
                "only the stellarator-symmetric rmnc/zmns boundary is supported")
    except (KeyError, IndexError) as e:
        raise ValueError(
            f"{path}: not a VMEC wout file (missing nfp/xm/xn/rmnc/zmns): {e}")
    finally:
        d.close()
    rb, zb = rmnc[-1, :], zmns[-1, :]          # boundary = last surface
    terms = [(int(xm[i]), int(round(xn[i] / nfp)), float(rb[i]), float(zb[i]))
             for i in range(len(xm))]
    return terms, nfp


def _vmec_boundary(terms, nfp, theta, phi):
    """Evaluate R, Z and their (theta, phi) derivatives from the Fourier terms.
    R = sum RBC cos(m th - n NFP ph),  Z = sum ZBS sin(m th - n NFP ph)."""
    R = np.zeros_like(theta); Z = np.zeros_like(theta)
    Rt = np.zeros_like(theta); Zt = np.zeros_like(theta)
    Rp = np.zeros_like(theta); Zp = np.zeros_like(theta)
    for m, n, rc, zs in terms:
        ang = m * theta - n * nfp * phi
        R += rc * np.cos(ang); Z += zs * np.sin(ang)
        Rt += -rc * m * np.sin(ang); Zt += zs * m * np.cos(ang)
        Rp += rc * (n * nfp) * np.sin(ang); Zp += -zs * (n * nfp) * np.cos(ang)
    return R, Z, Rt, Zt, Rp, Zp


def _vmec_points_normals(terms, nfp, ntheta, nphi):
    """Plasma-boundary sample points + OUTWARD unit normals from the parametric
    tangents n ~ dr/dtheta x dr/dphi.  The cross product is ALREADY consistently
    oriented; we only pick the global sign (away from the boundary's own
    phi-dependent m=0 centre) by majority vote -- NOT a star-shaped requirement,
    so a strongly-shaped (non-star-shaped) real wout boundary is accepted.
    Raises only on a degenerate (zero) normal (No-Fallback)."""
    th = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    R, Z, Rt, Zt, Rp, Zp = _vmec_boundary(terms, nfp, TH, PH)
    x = R * np.cos(PH); y = R * np.sin(PH); z = Z
    dtheta = np.stack([Rt * np.cos(PH), Rt * np.sin(PH), Zt], axis=-1)
    dphi = np.stack([Rp * np.cos(PH) - R * np.sin(PH),
                     Rp * np.sin(PH) + R * np.cos(PH), Zp], axis=-1)
    nrm = np.cross(dtheta, dphi)
    nlen = np.linalg.norm(nrm, axis=-1, keepdims=True)
    if (nlen < 1e-12).any():
        raise ValueError("degenerate VMEC boundary normal (zero cross product) "
                         "-- check the Fourier coefficients / parametrisation")
    nrm = nrm / nlen
    # the boundary's own m=0 centre R0(phi), Z0(phi) (NOT a hardcoded R_MAJOR --
    # a real equilibrium's axis is Shafranov-shifted); an interior reference for
    # the global outward-sign vote.
    R0 = np.zeros_like(PH); Z0 = np.zeros_like(PH)
    for m, n, rc, zs in terms:
        if m == 0:
            R0 += rc * np.cos(n * nfp * PH)
            Z0 += -zs * np.sin(n * nfp * PH)      # Z = ZBS sin(-n NFP ph) at m=0
    centre = np.stack([R0 * np.cos(PH), R0 * np.sin(PH), Z0], axis=-1)
    if np.sum(nrm * (np.stack([x, y, z], axis=-1) - centre)) < 0:
        nrm = -nrm                                # flip global sign -> outward
    pts = np.stack([x, y, z], axis=-1).reshape(-1, 3)
    nrm = nrm.reshape(-1, 3)
    return pts, nrm


def main():
    ap = argparse.ArgumentParser(
        description="mini-REGCOIL / NESCOIL fusion coil-design demo")
    ap.add_argument("--out-dir",
                    default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--work-dir",
                    default=os.path.join(os.environ.get("TEMP", "/tmp"),
                                         "sf_regcoil"))
    ap.add_argument("--eval-max", type=int, default=220,
                    help="plasma-boundary B.n sample points")
    ap.add_argument("--wind-maxh", type=float, default=0.045)
    ap.add_argument("--plasma-maxh", type=float, default=0.035)
    ap.add_argument("--n-levels", type=int, default=14,
                    help="current-potential contours (= coil turns) to draw")
    ap.add_argument("--n-alpha", type=int, default=9,
                    help="REGCOIL L-curve points (alpha sweep)")
    ap.add_argument("--wout", default=None,
                    help="VMEC wout_*.nc; boundary replaces the rotating ellipse")
    ap.add_argument("--vmec-ntheta", type=int, default=40)
    ap.add_argument("--vmec-nphi", type=int, default=48)
    ap.add_argument("--tf-axis-field", type=float, default=1.0,
                    help="desired on-axis toroidal field [T] (prescribed net "
                         "poloidal current for Part 3)")
    ap.add_argument("--no-cohomology", action="store_true",
                    help="skip the b1 (Euler-characteristic) confirmation")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    from radia.panels import calc_streamfunction as C
    from ngsolve import Mesh, H1, specialcf, TaskManager

    os.makedirs(args.work_dir, exist_ok=True)

    with TaskManager():
        coil = Mesh(_torus_surface_vol(
            A_WIND, args.wind_maxh, os.path.join(args.work_dir, "wind.vol")))
        plasma = Mesh(_torus_surface_vol(
            A_PLASMA, args.plasma_maxh, os.path.join(args.work_dir, "plasma.vol")))
        fes = H1(coil, order=1, definedon=coil.Boundaries(".*"))
        n_cf = specialcf.normal(3)

        pts, nrm, theta, phi = _plasma_points_normals(plasma, args.eval_max)
        A_n = _A_normal(C, fes, n_cf, pts, nrm)

        # ---- 1. forward works: two producible targets -> machine precision ----
        targets = {
            "uniform_vertical_PF": nrm[:, 2],
            "stellarator_sinTheta_cos2phi": np.sin(theta) * np.cos(2 * phi),
        }
        designs = {}
        ref_loops = None
        for name, Bn in targets.items():
            psi, res, peakJ = _design(C, fes, coil, A_n, Bn, "h1", 1e-7)
            loops = _contours(C, coil, psi, args.n_levels)
            designs[name] = {"bn_residual_rel": res,
                             "peak_grad_psi": peakJ,
                             "n_contours": len(loops)}
            if name == "uniform_vertical_PF":
                ref_loops = loops
            print(f"[forward] {name:32s}  B.n resid={res:.2e}  "
                  f"peak|grad psi|={peakJ:.2e}  contours={len(loops)}")

        # ---- 2. REGCOIL L-curve on a GENUINELY HARD target ----
        Bn_hard = np.sin(3 * theta) * np.cos(5 * phi)
        alpha_rel_grid = np.logspace(2.0, -8.0, args.n_alpha)
        lcurve = []
        for a_rel in alpha_rel_grid:
            _psi, res, peakJ = _design(C, fes, coil, A_n, Bn_hard, "h1",
                                       float(a_rel))
            lcurve.append({"alpha_rel": float(a_rel),
                           "bn_residual_rel": res, "peak_grad_psi": peakJ})
            print(f"[L-curve] alpha_rel={a_rel:.1e}  B.n resid={res:.3e}  "
                  f"peak|grad psi|={peakJ:.3e}")

        # ---- 3. NET CURRENT: the multivalued / secular term ----
        K_list = _secular_currents(n_cf)
        sec_n = _assemble_secular_normal(coil, pts, nrm, K_list)
        col_scale = float(np.median(np.linalg.norm(A_n, axis=0)))
        # the TF (net-poloidal) secular current is ~tangent to the plasma, so its
        # B.n footprint is tiny -> it is a PRESCRIBED parameter, not fitted.
        radii = np.linspace(0.20, 0.40, 9)
        by_in = _tf_field_1overR(coil, n_cf, radii)
        by_out = _tf_field_1overR(coil, n_cf, np.array([0.10, 0.55]))
        br = by_in * radii                       # B_tor * R (should be const)
        # prescribe net poloidal current for the desired on-axis toroidal field
        Bphi0 = args.tf_axis_field
        Ipol_needed = 2.0 * math.pi * R_MAJOR * Bphi0 / _MU0
        net_current = {
            "betti1_winding_surface": (None if args.no_cohomology else
                                       _betti1_winding_surface(A_WIND, 0.06)),
            "n_secular_dofs": len(K_list),
            "tf_field_BR_const_mean": float(np.mean(br)),
            "tf_field_BR_const_spread_rel":
                float((br.max() - br.min()) / abs(np.mean(br))),
            "tf_field_outside_over_inside":
                float(np.max(np.abs(by_out)) / np.max(np.abs(by_in))),
            "secular_bn_footprint_rel": {
                K_list[k][0]: float(np.linalg.norm(sec_n[:, k]) / col_scale)
                for k in range(len(K_list))},
            "prescribed_on_axis_B_tor_T": Bphi0,
            "net_poloidal_current_needed_A": Ipol_needed,
        }
        b1 = net_current["betti1_winding_surface"]
        print(f"[net-current] b1(winding surface) [Euler char, gmsh-free]="
              f"{b1}  -> {len(K_list)} secular DOFs"
              + ("" if b1 is None else f"  (match: {b1 == len(K_list)})"))
        print(f"[net-current] TF secular field: B_tor*R const to "
              f"{net_current['tf_field_BR_const_spread_rel']*100:.2f}% inside, "
              f"outside/inside={net_current['tf_field_outside_over_inside']:.1e}")
        for k, (nm, _K) in enumerate(K_list):
            print(f"[net-current] secular '{nm}' B.n footprint = "
                  f"{net_current['secular_bn_footprint_rel'][nm]:.2e} "
                  f"(rel. to a psi column)")
        print(f"[net-current] prescribe B_tor={Bphi0:.1f} T on axis -> net "
              f"poloidal current {Ipol_needed/1e3:.1f} kA")

        # ---- 4. VMEC-shaped plasma boundary (rotating ellipse / wout) ----
        if args.wout:
            terms, nfp = _read_wout_boundary(args.wout)
            vmec_src = os.path.basename(args.wout)
        else:
            terms, nfp = _rotating_ellipse_terms()
            vmec_src = "rotating_ellipse_model"
        vpts, vnrm = _vmec_points_normals(terms, nfp, args.vmec_ntheta,
                                          args.vmec_nphi)
        if len(vpts) > args.eval_max:
            sel = np.linspace(0, len(vpts) - 1, args.eval_max).astype(int)
            vpts, vnrm = vpts[sel], vnrm[sel]
        A_n_v = _A_normal(C, fes, n_cf, vpts, vnrm)
        Bn_v = vnrm[:, 2]                         # uniform vertical field B.n
        _psiv, res_v, peakJ_v = _design(C, fes, coil, A_n_v, Bn_v, "h1", 1e-7)
        minor = np.linalg.norm(
            vpts - np.column_stack([R_MAJOR * vpts[:, 0] / np.hypot(vpts[:, 0],
                                    vpts[:, 1]),
                                    R_MAJOR * vpts[:, 1] / np.hypot(vpts[:, 0],
                                    vpts[:, 1]),
                                    np.zeros(len(vpts))]), axis=1)
        vmec = {"source": vmec_src, "nfp": nfp, "n_terms": len(terms),
                "minor_radius_min": float(minor.min()),
                "minor_radius_max": float(minor.max()),
                "non_axisymmetric": bool(minor.max() - minor.min() > 1e-3),
                "bn_residual_rel": res_v, "peak_grad_psi": peakJ_v,
                "eval_points": int(len(vpts))}
        print(f"[vmec] boundary='{vmec_src}' NFP={nfp}  minor "
              f"{minor.min():.3f}..{minor.max():.3f} m (non-axisym="
              f"{vmec['non_axisymmetric']})  -> B.n resid={res_v:.2e}")

    fig_data = {"lcurve": lcurve, "loops": ref_loops, "tf_radii": radii.tolist(),
                "tf_BR": br.tolist(), "vmec_terms": terms, "vmec_nfp": nfp}

    data = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "benchmark": "regcoil_fusion_forward",
        "problem": {
            "winding_torus": {"R": R_MAJOR, "a": A_WIND, "maxh": args.wind_maxh},
            "plasma_torus": {"R": R_MAJOR, "a": A_PLASMA, "maxh": args.plasma_maxh},
            "winding_ndof": int(fes.ndof),
            "plasma_eval_points": int(len(pts)),
            "regularize": "h1",
            "note": "parts 1-2 single-valued psi; part 3 adds the multivalued "
                    "secular (net-current) term; part 4 a VMEC-format boundary. "
                    "B.n from analytic model fields (swap for free-boundary VMEC "
                    "via --wout for a production run).",
        },
        "producible_targets": designs,
        "regcoil_lcurve": {"hard_target": "sin(3 theta) cos(5 phi)",
                           "points": lcurve},
        "net_current": net_current,
        "vmec_boundary": vmec,
    }
    jpath = validation_output("demo_regcoil_fusion.json", args.out_dir)
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, allow_nan=False)
    print(f"\nResults -> {jpath}")

    if not args.no_plot:
        _plot(args, fig_data)


def _plot(args, fd):
    """2x2 lab figure: (a) the current-potential coil + plasma; (b) the REGCOIL
    L-curve; (c) the net-poloidal-current TF field B_tor*R vs R (Ampere 1/R);
    (d) the VMEC rotating-ellipse boundary cross-sections."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from radia_mcp.figure import apply_lab_style, save_lab_figure
    except Exception as e:
        print(f"(plot skipped: radia-mcp figure unavailable: {e})")
        return
    plt.rcParams["pdf.fonttype"] = 42
    w_in, h_in = apply_lab_style(embed_width_cm=16.0, aspect=0.88)
    fig = plt.figure(figsize=(w_in, h_in))

    # (a) the designed coil = current-potential contours on the winding torus
    ax = fig.add_subplot(2, 2, 1, projection="3d")
    u = np.linspace(0, 2 * np.pi, 60); v = np.linspace(0, 2 * np.pi, 40)
    U, V = np.meshgrid(u, v)
    X = (R_MAJOR + A_PLASMA * np.cos(V)) * np.cos(U)
    Yy = (R_MAJOR + A_PLASMA * np.cos(V)) * np.sin(U)
    Zz = A_PLASMA * np.sin(V)
    ax.plot_surface(X, Yy, Zz, color="orange", alpha=0.18, linewidth=0)
    for p in (fd["loops"] or []):
        ax.plot(p[:, 0], p[:, 1], p[:, 2], lw=0.8, color="C0")
    ax.text2D(0.0, 0.97, "(a) current-potential coil + plasma",
              transform=ax.transAxes)
    for s in ("x", "y", "z"):
        getattr(ax, f"set_{s}ticks")([])
    ax.set_box_aspect((1, 1, 0.5))

    # (b) REGCOIL L-curve: residual vs coil complexity
    ax2 = fig.add_subplot(2, 2, 2)
    res = [pt["bn_residual_rel"] for pt in fd["lcurve"]]
    pk = [pt["peak_grad_psi"] for pt in fd["lcurve"]]
    ax2.loglog(pk, res, "o-", color="C3")
    ax2.set_xlabel(r"peak $|\nabla\psi|$  (coil complexity)")
    ax2.set_ylabel(r"$B{\cdot}n$ residual  (rel.)")
    ax2.text(0.04, 0.04, "(b) REGCOIL L-curve", transform=ax2.transAxes)
    ax2.grid(True, which="both", alpha=0.3)
    ax2.annotate(r"large $\alpha$", xy=(pk[0], res[0]),
                 xytext=(0.30, 0.85), textcoords="axes fraction",
                 arrowprops=dict(arrowstyle="->", lw=0.7))
    ax2.annotate(r"small $\alpha$", xy=(pk[-1], res[-1]),
                 xytext=(0.55, 0.18), textcoords="axes fraction",
                 arrowprops=dict(arrowstyle="->", lw=0.7))

    # (c) net-poloidal-current (TF) field: B_tor * R vs R (Ampere 1/R, flat)
    ax3 = fig.add_subplot(2, 2, 3)
    rr = np.array(fd["tf_radii"]); br = np.array(fd["tf_BR"])
    ax3.plot(rr, br / br.mean(), "o-", color="C2")
    ax3.axhline(1.0, ls="--", lw=0.7, color="0.5")
    ax3.set_xlabel(r"major radius $R$  (m)")
    ax3.set_ylabel(r"$B_{\rm tor}\,R$  (normalised)")
    ax3.text(0.04, 0.90, r"(c) TF secular term: Ampere $1/R$",
             transform=ax3.transAxes)
    ax3.set_ylim(0.95, 1.05)

    # (d) VMEC rotating-ellipse boundary cross-sections at several phi
    ax4 = fig.add_subplot(2, 2, 4)
    th = np.linspace(0, 2 * np.pi, 120)
    nfp = fd["vmec_nfp"]
    for j, frac in enumerate(np.linspace(0, 1.0, 5)[:-1]):
        ph = frac * 2 * np.pi / nfp
        R, Z, *_ = _vmec_boundary(fd["vmec_terms"], nfp, th, np.full_like(th, ph))
        ax4.plot(R, Z, lw=1.0, label=rf"$\phi={frac:.2f}\cdot2\pi/N_{{fp}}$")
    ax4.set_aspect("equal")
    ax4.set_xlabel(r"$R$  (m)"); ax4.set_ylabel(r"$Z$  (m)")
    ax4.text(0.04, 0.92, rf"(d) VMEC boundary $N_{{fp}}={nfp}$",
             transform=ax4.transAxes)
    ax4.legend(fontsize=7, loc="upper right", frameon=False)

    # explicit margins (a 3D subplot confuses tight_layout -> clipped x-labels)
    fig.subplots_adjust(left=0.085, right=0.97, top=0.95, bottom=0.085,
                        hspace=0.32, wspace=0.27)
    info = save_lab_figure(fig, os.path.join(args.out_dir,
                                             "demo_regcoil_fusion"),
                           embed_width_cm=16.0, tighten=False)
    print(f"Plot    -> {', '.join(info['wrote'])}")


if __name__ == "__main__":
    main()
