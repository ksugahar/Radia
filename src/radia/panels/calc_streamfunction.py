#!/usr/bin/env python3
"""calc_streamfunction.py -- Stream-Function (SF) coil design (Layer 4, headless).

Design a surface coil on a STANDALONE 2D surface mesh (.vol) that produces a
target field over an evaluation region (.vol, surface OR volume).  The target
and the computed field are NGSolve CoefficientFunctions (no .sol / .csv).

I/O (per the 2026-06 design decisions):
  --coil-vol   coil surface mesh (.vol, 2D-in-3D).  psi = H1 on it (FE-direct).
  --eval-vol   evaluation region (.vol).  Surface OR volume; the homogeneity
               adapts to the mesh dimension automatically.
  --target-cf  target field as a CoefficientFunction EXPRESSION of x, y, z.
               SCALAR expr -> Bz target (e.g. "x" = Gx, "1" = uniform Bz,
               "x*x-y*y" = C2 ellipse shim).  3-VECTOR expr "(Bx,By,Bz)" ->
               full vector-B target.
  --target-harmonic  target as a SPHERICAL HARMONIC (MRI gradient/shim basis):
               a name (Z0,X,Y,Z,Z2,ZX,ZY,C2,S2,Z3,..) or 'l=L,m=M', optionally
               weighted and summed, e.g. 'Z2' (pure 2nd-order shim), 'X' (Gx),
               'Z2:1.0,Z:0.1'.  Generates the solid-harmonic --target-cf
               polynomial.  Give --target-cf OR --target-harmonic.

In design mode the ACHIEVED Bz over the DSV is decomposed into solid
spherical harmonics (``result["harmonics"]``): the per-(l,m) field-RMS
spectrum, the LSQ residual, and -- when a --target-harmonic is given --
the PURITY (target-harmonic field fraction) and the largest CONTAMINANT.
These are the canonical gradient-coil quality metrics.

THREE MODES (--method):
  design       target -> A psi = B (folded-Tikhonov RegularizedTSVD) -> psi,
               field homogeneity over the eval region, peak surface current.
  pareto       (homogeneity, peak current density) Pareto front via a chosen
               --pareto-lever:
                 alpha     Tikhonov L-curve (one factorisation, swept alpha)
                 linf      minimax IRLS trajectory (peak current down at
                           ~constant homogeneity; re-uses the ONE ACA factor)
                 geometry  eval-region (DSV) scale sweep (dual of former size)
  manufacture  psi iso-contours (marching triangles) -> nearest-neighbour
               single-stroke wire -> optional sheet-metal wire distortion
               (--distort) -> CAD STEP (--step-output) -> PEEC L,R (--peec).
               Reports the discretised-wire field homogeneity (design loop
               closed).  GENERAL on any loaded surface mesh.

Surface-FE convention (verified): a standalone surface .vol loads with ne=0
(surface elements as boundary), so psi uses
``H1(coil, definedon=coil.Boundaries('.*'))`` with ``grad(v).Trace()`` and
``* ds`` -- the sphere-demo machinery generalised to a loaded mesh.
"""
import os
import sys
import argparse
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                                   # calc_common
sys.path.insert(0, os.path.join(_HERE, "..", ".."))         # src/radia package root
sys.path.insert(0, os.path.join(_HERE, ".."))               # src/radia

from calc_common import calc_main  # noqa: E402

_MU0_4PI = 1.0e-7


def _parse_target_cf(expr):
    """Eval a CoefficientFunction expression of x, y, z.  Returns (cf, dim)."""
    from ngsolve import x, y, z, CoefficientFunction, sin, cos, exp, sqrt, log
    import math
    scope = {"x": x, "y": y, "z": z, "CF": CoefficientFunction,
             "sin": sin, "cos": cos, "exp": exp, "sqrt": sqrt, "log": log,
             "pi": math.pi}
    val = eval(expr, {"__builtins__": {}}, scope)            # noqa: S307
    cf = CoefficientFunction(val)
    return cf, cf.dim


def _eval_target_np(expr, pts, vector_b):
    """Evaluate the target field EXPRESSION numerically at arbitrary points
    (no mesh containment needed -- so the geometry lever can scale the eval
    region beyond the original mesh).  Returns a flat array (1 or 3 / point)."""
    import math

    def _cf(*a):
        return tuple(a) if len(a) > 1 else a[0]

    out = []
    for p in pts:
        scope = {"x": float(p[0]), "y": float(p[1]), "z": float(p[2]),
                 "CF": _cf, "sin": np.sin, "cos": np.cos, "exp": np.exp,
                 "sqrt": np.sqrt, "log": np.log, "pi": math.pi}
        v = eval(expr, {"__builtins__": {}}, scope)            # noqa: S307
        out.extend([float(c) for c in v] if vector_b else [float(v)])
    return np.array(out, dtype=float)


# ==========================================================================
# Spherical-harmonic target + field decomposition (MRI gradient / shim basis).
#
# In a current-free region Bz is harmonic, so it expands in REAL REGULAR SOLID
# HARMONICS R_l^m(x,y,z) -- homogeneous harmonic polynomials (Laplace nabla^2
# R = 0).  These ARE the named MRI gradients/shims (Z0, X, Y, Z, Z2, ZX, ..).
# We store each as a Cartesian polynomial STRING and use the SAME table for
# BOTH the target (--target-harmonic) and the analysis of the achieved field
# (purity / contamination), so the round-trip is exact.  m>0 = cos(m.phi) (C),
# m<0 = sin(|m|.phi) (S); the leading (unnormalised) Cartesian form is used
# consistently throughout, so normalisation cancels in every reported ratio.
# ==========================================================================
_SOLID_HARMONICS = {
    (0, 0):  ("Z0",  "(1.0 + 0.0*x)"),
    (1, 0):  ("Z",   "(z)"),
    (1, 1):  ("X",   "(x)"),
    (1, -1): ("Y",   "(y)"),
    (2, 0):  ("Z2",  "(z*z - 0.5*(x*x + y*y))"),
    (2, 1):  ("ZX",  "(x*z)"),
    (2, -1): ("ZY",  "(y*z)"),
    (2, 2):  ("C2",  "(x*x - y*y)"),
    (2, -2): ("S2",  "(x*y)"),
    (3, 0):  ("Z3",  "(z*z*z - 1.5*z*(x*x + y*y))"),
    (3, 1):  ("Z2X", "(x*(4.0*z*z - x*x - y*y))"),
    (3, -1): ("Z2Y", "(y*(4.0*z*z - x*x - y*y))"),
    (3, 2):  ("ZC2", "(z*(x*x - y*y))"),
    (3, -2): ("ZS2", "(x*y*z)"),
    (3, 3):  ("C3",  "(x*(x*x - 3.0*y*y))"),
    (3, -3): ("S3",  "(y*(3.0*x*x - y*y))"),
    (4, 0):  ("Z4",   "(z*z*z*z - 3.0*z*z*(x*x + y*y) "
                      "+ 0.375*(x*x + y*y)*(x*x + y*y))"),
    (4, 1):  ("Z3X",  "(x*z*(4.0*z*z - 3.0*x*x - 3.0*y*y))"),
    (4, -1): ("Z3Y",  "(y*z*(4.0*z*z - 3.0*x*x - 3.0*y*y))"),
    (4, 2):  ("Z2C2", "((x*x - y*y)*(6.0*z*z - x*x - y*y))"),
    (4, -2): ("Z2S2", "(x*y*(6.0*z*z - x*x - y*y))"),
    (4, 3):  ("ZC3",  "(z*x*(x*x - 3.0*y*y))"),
    (4, -3): ("ZS3",  "(z*y*(3.0*x*x - y*y))"),
    (4, 4):  ("C4",   "(x*x*x*x - 6.0*x*x*y*y + y*y*y*y)"),
    (4, -4): ("S4",   "(x*y*(x*x - y*y))"),
}
_HARMONIC_BY_NAME = {nm.upper(): lm for lm, (nm, _) in _SOLID_HARMONICS.items()}
_HARMONIC_BY_NAME.update({"X2Y2": (2, 2), "XY": (2, -2), "GX": (1, 1),
                          "GY": (1, -1), "GZ": (1, 0), "Z1": (1, 0)})


def _harmonic_lookup(token):
    """Resolve one harmonic token -> (l, m).  Accepts a name (Z2, ZX, ..), an
    ``l=L,m=M`` pair, or a ``(L,M)`` tuple string.  Raises (No-Fallback) with
    the available names on an unknown token."""
    t = token.strip()
    tu = t.upper()
    if tu in _HARMONIC_BY_NAME:
        return _HARMONIC_BY_NAME[tu]
    s = t.replace("(", "").replace(")", "").replace(" ", "")
    bad = ValueError(                       # one clear, actionable message
        f"unknown harmonic '{token}'. Names: "
        f"{sorted({nm for nm, _ in _SOLID_HARMONICS.values()})}; "
        f"or 'l=L,m=M' / '(L,M)' with integer L,M.")
    # Fail LOUD but CLEARLY: a malformed pair ('l=4' with no m, '(4,)',
    # 'l=x,m=y') must give the actionable message above, NOT a raw KeyError
    # or int() ValueError.  (No-Fallback: this converts the parse exception
    # into one clear user error -- it does not try another path.)
    if s.lower().startswith("l="):
        try:
            d = dict(kv.split("=") for kv in s.split(","))
            lm = (int(d["l"]), int(d["m"]))
        except (KeyError, ValueError, IndexError):
            raise bad from None
    elif "," in s:
        try:
            a, b = s.split(",")[:2]
            lm = (int(a), int(b))
        except (ValueError, IndexError):
            raise bad from None
    else:
        raise bad
    if lm not in _SOLID_HARMONICS:
        raise ValueError(f"harmonic {lm} not in the table (l<=4 supported).")
    return lm


def _split_harmonic_terms(spec):
    """Split a --target-harmonic spec into its terms.  Commas separate terms,
    BUT a comma also appears INSIDE the ``l=L,m=M`` and ``(L,M)`` single-
    harmonic forms; this merges those back so e.g. ``l=4,m=-4`` and
    ``(2,0),Z:0.1`` parse correctly (names like ``Z2`` never contain a comma)."""
    raw = spec.split(",")
    out, i = [], 0
    while i < len(raw):
        frag = raw[i]
        while frag.count("(") > frag.count(")") and i + 1 < len(raw):
            i += 1
            frag += "," + raw[i]                 # rejoin a split (L,M) pair
        low = frag.lower().replace(" ", "")
        if low.startswith("l=") and "m=" not in low and i + 1 < len(raw):
            i += 1
            frag += "," + raw[i]                 # rejoin a split l=L,m=M pair
        if frag.strip():
            out.append(frag.strip())
        i += 1
    return out


def _harmonic_target_expr(spec):
    """Build the target-cf polynomial STRING for a spherical-harmonic spec.

    ``spec`` = comma-separated terms ``NAME`` or ``NAME:coeff`` (or
    ``l=L,m=M:coeff`` / ``(L,M):coeff``), e.g. ``Z2`` (pure Z2 shim),
    ``X`` (Gx gradient), ``Z2:1.0,Z:0.1`` (Z2 with a Z offset).  Returns
    ``(expr_string, [(l, m, coeff), ...])``."""
    terms, parts = [], []
    for chunk in _split_harmonic_terms(spec):
        chunk = chunk.strip()
        if not chunk:
            continue
        # a coeff is the part after the LAST ':' only if it parses as a float
        coeff = 1.0
        token = chunk
        if ":" in chunk:
            head, tail = chunk.rsplit(":", 1)
            try:
                coeff = float(tail)
                token = head
            except ValueError:
                token = chunk            # ':' belonged to the token (none today)
        lm = _harmonic_lookup(token)
        terms.append((lm[0], lm[1], coeff))
        parts.append(f"({coeff})*{_SOLID_HARMONICS[lm][1]}")
    if not terms:
        raise ValueError(f"empty --target-harmonic spec '{spec}'")
    return " + ".join(parts), terms


def _solid_harmonic_matrix(pts, lmax):
    """(H, labels): H[i, j] = R_lj(pts[i]); labels[j] = (l, m, name).  Columns
    are every table harmonic with l <= lmax."""
    xs, ys, zs = pts[:, 0], pts[:, 1], pts[:, 2]
    cols, labels = [], []
    for (l, m), (nm, expr) in sorted(_SOLID_HARMONICS.items()):
        if l > lmax:
            continue
        val = eval(expr, {"__builtins__": {}},                 # noqa: S307
                   {"x": xs, "y": ys, "z": zs})
        cols.append(np.broadcast_to(np.asarray(val, float), xs.shape))
        labels.append((l, m, nm))
    return np.column_stack(cols), labels


def _harmonic_decompose(pts, bz, lmax, target_terms=None):
    """Least-squares decomposition of an achieved Bz sample into solid
    harmonics over the eval points.  Returns the per-term field RMS spectrum,
    the LSQ residual fraction, and (if a target is given) the purity and the
    largest contaminant -- the canonical gradient-coil quality metrics.

    purity = RMS of the target harmonic's field / RMS of the full fit field
    (an energy fraction; the solid harmonics are orthogonal over a sphere, so
    over a roughly-spherical DSV this is the standard practical measure)."""
    H, labels = _solid_harmonic_matrix(pts, lmax)
    coeff, *_ = np.linalg.lstsq(H, bz, rcond=None)
    fit = H @ coeff
    total = float(np.sqrt(np.mean(fit * fit))) + 1e-300
    bz_rms = float(np.sqrt(np.mean(bz * bz))) + 1e-300
    term_rms = np.array([float(np.sqrt(np.mean((H[:, j] * coeff[j]) ** 2)))
                         for j in range(len(labels))])
    spectrum = []
    for j, (l, m, nm) in enumerate(labels):
        spectrum.append({"l": l, "m": m, "name": nm,
                         "coeff": float(coeff[j]),
                         "field_rms_T": term_rms[j],
                         "fraction": float(term_rms[j] / total)})
    spectrum.sort(key=lambda d: -d["field_rms_T"])
    out = {"lmax": lmax,
           "residual_fraction": float(np.linalg.norm(bz - fit)
                                      / (np.linalg.norm(bz) + 1e-300)),
           "spectrum": spectrum,
           "dominant": {"name": spectrum[0]["name"], "l": spectrum[0]["l"],
                        "m": spectrum[0]["m"], "coeff": spectrum[0]["coeff"]}}
    if target_terms:
        tset = {(l, m) for (l, m, _c) in target_terms}
        e_tot = float(np.sum(term_rms ** 2))
        e_tgt = float(np.sum([term_rms[j] ** 2 for j, (l, m, _n)
                              in enumerate(labels) if (l, m) in tset]))
        contam = [(labels[j][2], term_rms[j] / total) for j in range(len(labels))
                  if (labels[j][0], labels[j][1]) not in tset]
        contam.sort(key=lambda t: -t[1])
        out["purity"] = float(np.sqrt(e_tgt / (e_tot + 1e-300)))
        out["max_contaminant"] = ({"name": contam[0][0], "fraction": contam[0][1]}
                                  if contam else None)
    return out


def _resolve_target_harmonic(args):
    """If ``--target-harmonic`` is given, fill ``args.target_cf`` with the
    generated polynomial (idempotent).  Stores the parsed (l,m,coeff) terms on
    ``args._harmonic_terms`` for the analysis stage.  Exactly one of
    ``--target-cf`` / ``--target-harmonic`` must be set."""
    spec = getattr(args, "target_harmonic", "") or ""
    if getattr(args, "_harmonic_terms", None) is not None:
        return                                                  # already resolved
    if spec:
        if args.target_cf:
            raise ValueError("give --target-cf OR --target-harmonic, not both.")
        expr, terms = _harmonic_target_expr(spec)
        args.target_cf = expr
        args._harmonic_terms = terms
    else:
        args._harmonic_terms = []
        if not args.target_cf:
            raise ValueError("need --target-cf or --target-harmonic.")


def _assemble_biot_savart(fes, n, pts, comps):
    """A[row, dof]: field component(s) at each eval point from the FE-direct
    surface current K = n x grad_s phi (Setup B: grad(v).Trace(), * ds)."""
    from ngsolve import grad, ds, LinearForm, sqrt as ng_sqrt, Cross, x, y, z
    ndof = fes.ndof
    u, v = fes.TnT()
    A = np.zeros((len(pts) * len(comps), ndof))
    row = 0
    for p in pts:
        dxt, dyt, dzt = p[0] - x, p[1] - y, p[2] - z
        r2 = dxt * dxt + dyt * dyt + dzt * dzt
        r3 = r2 * ng_sqrt(r2)
        Kv = Cross(n, grad(v).Trace())
        cross = (Kv[1] * dzt - Kv[2] * dyt,
                 Kv[2] * dxt - Kv[0] * dzt,
                 Kv[0] * dyt - Kv[1] * dxt)
        for c in comps:
            Lf = LinearForm(fes)
            Lf += (_MU0_4PI / r3) * cross[c] * ds
            Lf.Assemble()
            A[row, :] = Lf.vec.FV().NumPy()
            row += 1
    return A


def _assemble_iron_reaction(fes, n, obs_pts, comps, iron_vol, mu_r,
                            iron_mat="iron", fem_order=2, ndof_cap=20000,
                            exact_source=False, quad_order=None):
    """A_react[row, dof]: the MATERIAL (iron) reaction field at each obs point
    from the FE-direct surface current of each psi-DOF, via the reduced scalar
    potential ``H = H_s - grad(Omega)`` on a Kelvin open-boundary FE mesh.

    The iron geometry is supplied as a ``.vol`` (``iron_vol``) carrying an iron
    material (``iron_mat``), the Kelvin ball material ``kelvin`` and the periodic
    ``kelvin_int``/``kelvin_ext`` identification + a ``GND`` vertex -- the same
    Kelvin pipeline ``calc_fem_kelvin`` consumes.  Returned A_react is in the
    SAME units (Tesla) and row order (``pt*len(comps)+comp``) as
    ``_assemble_biot_savart`` so the caller can do ``Ac += A_react``.

    Verified mechanism: ``validation_test/kelvin_dtn_spectrum/act8_03_general_iron_design.py``
    (M = M_free + M_react; material-aware design HITS, free-space design
    MISSES).  Source quality: the DEFAULT is the P1-nodal
    injection of B_s on the iron vertices (cheap); ``exact_source=True``
    evaluates B_s at the iron QUADRATURE points (analytic Biot-Savart, no P1
    interpolation -> removes the iron-vertex plateau, ~14% on a coarse shell) at
    the cost of one ``_assemble_biot_savart`` per quad point (much heavier --
    opt-in).  Cost: one Kelvin factorisation + one back-substitution PER OBS ROW
    (the obs-adjoint, NOT per psi-DOF); the adjoint identity is
    ``A_react[i,j] = (mu_r-1) int_iron B_s^j . grad(X_i)``, X_i = Ainv G^T.
    Verified == the per-psi-DOF forward to 1e-14.  ``ndof_cap`` bounds the dense
    source matrix memory."""
    from ngsolve import (Mesh, H1, Periodic, VectorH1, BilinearForm,
                         GridFunction, grad, dx, InnerProduct, ElementId, VOL,
                         IntegrationRule, x, y, z)
    sys.path.insert(0, os.path.dirname(__file__))
    from calc_common import detect_kelvin_offset

    if fes.ndof > ndof_cap:
        raise ValueError(
            f"--iron-vol material-aware kernel forms the dense source matrix "
            f"B_iron (3*n_iron x N_psi); N_psi={fes.ndof} > cap {ndof_cap}. "
            f"Coarsen the coil mesh / lower --order, or raise the cap if memory "
            f"allows (the back-substitution count is the obs rows, not N_psi).")

    mesh = Mesh(iron_vol)
    mats = set(mesh.GetMaterials())
    if iron_mat not in mats:
        raise ValueError(f"--iron-vol has no '{iron_mat}' material; "
                         f"available: {sorted(mats)}")
    if "kelvin" not in mats:
        raise ValueError(f"--iron-vol has no 'kelvin' material (open-boundary "
                         f"ball); available: {sorted(mats)}")
    if mesh.ngmesh.GetNrIdentifications() == 0:
        raise ValueError("--iron-vol has 'kelvin' material but no periodic "
                         "identification (kelvin_int<->kelvin_ext). Re-export "
                         "with the Kelvin pipeline.")

    # Kelvin ball geometry from the mesh: centre via detect_kelvin_offset,
    # radius = max kelvin-vertex distance from that centre (same as calc_fem_kelvin).
    kc = np.array(detect_kelvin_offset(mesh))
    kverts = sorted({v.nr for el in mesh.Elements() if el.mat == "kelvin"
                     for v in el.vertices})
    if not kverts:
        raise ValueError("--iron-vol lists 'kelvin' material but no volume "
                         "elements are tagged 'kelvin'.")
    Vc = np.array([list(mesh.vertices[i].point) for i in range(mesh.nv)])
    a_kelvin = float(np.max(np.linalg.norm(Vc[kverts] - kc, axis=1)))

    rp2 = (x - kc[0]) ** 2 + (y - kc[1]) ** 2 + (z - kc[2]) ** 2 + 1e-20
    mu = mesh.MaterialCF({iron_mat: float(mu_r), "kelvin": a_kelvin ** 2 / rp2},
                         default=1.0)
    fesO = Periodic(H1(mesh, order=fem_order, dirichlet="GND"))
    u, w = fesO.TnT()
    A = BilinearForm(mu * grad(u) * grad(w) * dx(bonus_intorder=2 * fem_order + 2))
    A.Assemble()
    Ainv = A.mat.Inverse(fesO.FreeDofs(), inverse="sparsecholesky")

    # OBS-ADJOINT shared covector G[row] = "-grad(.)(obs)[comp]" on fesO (element-
    # local basis gradients at the point).  A_react[i,j] = (mu_r-1) int_iron
    # B_s^j . grad(X_i), X_i = Ainv G[i] (one back-sub per obs row, not per
    # psi-DOF).  Verified == the per-psi-DOF forward to 1e-14.
    n_rows = len(obs_pts) * len(comps)
    G = np.zeros((n_rows, fesO.ndof))
    gt = GridFunction(fesO); arrgt = gt.vec.FV().NumPy()
    row = 0
    for p in obs_pts:
        mp = mesh(*p)
        dofs = fesO.GetDofNrs(ElementId(VOL, mp.nr))
        for c in comps:
            for gd in dofs:
                if gd < 0:
                    continue
                arrgt[:] = 0.0; arrgt[gd] = 1.0
                G[row, gd] = -grad(gt)(mp)[c]
            row += 1

    if exact_source:
        # EXACT source: B_s^j at the iron QUADRATURE points (analytic Biot-Savart,
        # no P1 interpolation -> removes the iron-vertex plateau).  Heavier: one
        # _assemble_biot_savart column block per quad point.
        qorder = int(quad_order) if quad_order else 2 * fem_order
        qpts, qw = [], []
        for el in mesh.Elements(VOL):
            if el.mat != iron_mat:
                continue
            trafo = mesh.GetTrafo(ElementId(VOL, el.nr))
            for ip in IntegrationRule(el.type, qorder):
                mip = trafo(ip)
                qpts.append([mip.point[0], mip.point[1], mip.point[2]])
                qw.append(ip.weight * mip.measure)
        qpts = np.array(qpts); qw = np.array(qw)
        B_q = _assemble_biot_savart(fes, n, qpts, [0, 1, 2]).reshape(
            len(qpts), 3, fes.ndof)
        mps = [mesh(*p) for p in qpts]
        gXi = GridFunction(fesO); src = gXi.vec.CreateVector()
        A_react = np.zeros((n_rows, fes.ndof))
        for i in range(n_rows):
            src.FV().NumPy()[:] = G[i]
            gXi.vec.data = Ainv * src          # X_i = Ainv G[i]
            GX = np.array([[grad(gXi)(mp)[c] for c in range(3)] for mp in mps])
            A_react[i, :] = (float(mu_r) - 1.0) * np.einsum(
                'qc,qcj->j', GX * qw[:, None], B_q)
        return A_react

    # P1 source (default, cheap): B_s sampled at the iron VERTICES (VectorH1
    # order 1), P1-interpolated in the coupling.  mu0 cancels (inject B_s, read
    # -grad(Omega) -> B_react in Tesla).  A_react = einsum(E_iron, B_iron),
    # E_iron = (D^T Ainv G^T) read on the iron-vertex source dofs.
    Vsrc = VectorH1(mesh, order=1)
    vt = Vsrc.TestFunction()
    DT = BilinearForm(trialspace=fesO, testspace=Vsrc)   # D^T : fesO -> Vsrc
    DT += (float(mu_r) - 1.0) * InnerProduct(vt, grad(u)) * \
        dx(definedon=mesh.Materials(iron_mat), bonus_intorder=4)
    DT.Assemble()
    iron = np.array(sorted({v.nr for el in mesh.Elements() if el.mat == iron_mat
                            for v in el.vertices}), dtype=int)
    B_iron = _assemble_biot_savart(fes, n, Vc[iron], [0, 1, 2]).reshape(
        len(iron), 3, fes.ndof)
    gX = GridFunction(fesO); gE = GridFunction(Vsrc)
    tmp = gX.vec.CreateVector()
    E_iron = np.zeros((n_rows, 3, len(iron)))
    for i in range(n_rows):
        gX.vec.FV().NumPy()[:] = G[i]
        tmp.data = Ainv * gX.vec               # X_i = Ainv G[i]   (A symmetric)
        gE.vec.data = DT.mat * tmp             # e_i = D^T X_i  (Vsrc)
        for k in range(3):
            E_iron[i, k, :] = gE.components[k].vec.FV().NumPy()[iron]
    return np.einsum('ikv,vkj->ij', E_iron, B_iron)


def _seminorm(fes, regularize):
    """SPD seminorm S over ALL DOFs as a SPARSE (scipy CSR) matrix: l2 =
    identity, h1 = surface H1 stiffness.  The caller reduces it to the
    independent-DOF space via the R matrix (``S_ind = R^T S R``).

    Sparse end-to-end on purpose: the FE stiffness is ~2 % dense, so the old
    ``np.array(a.mat.ToDense())`` was an O(N^2) memory/time wall (1.9 GB and
    ~13 s of dense folding at N=15k DOF -- see bench_sf_scaling).  CSR keeps it
    O(N); RegularizedTSVD.from_stiffness already has a scipy-sparse path."""
    import scipy.sparse as sp
    if regularize == "l2":
        return sp.identity(fes.ndof, format="csr")
    from ngsolve import BilinearForm, grad, ds, TaskManager
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u).Trace() * grad(v).Trace() * ds
    a += 1.0e-10 * u * v * ds
    with TaskManager():
        a.Assemble()
    r, c, val = a.mat.COO()
    return sp.csr_matrix((np.asarray(val), (np.asarray(r), np.asarray(c))),
                         shape=(a.mat.height, a.mat.width))


def _inductance_seminorm(coil, fes, order_J=0, ridge=1e-9):
    """Magnetic self-inductance quadratic form ``L_psi`` over the psi DOFs:
    ``psi^T L_psi psi = mu0 K^T SL K``, K = n x grad psi the surface current,
    SL the BEM single-layer (Laplace) operator (``ngsolve.bem.LaplaceSL``).

    ``min 1/2 psi^T L_psi psi  s.t.  A psi = B`` is the canonical MINIMUM-
    STORED-ENERGY gradient-coil objective (Turner / Forbes target-field
    method): low stored energy = fast slew rate.  Used as the regularisation
    seminorm S, so the design is ``RegularizedTSVD.from_stiffness(base,
    L_psi)`` -- the SAME folded-TSVD machinery as l2/h1, with a physically
    meaningful S instead of a smoothness proxy.

    Construction: K = n x grad psi maps H1 -> HDivSurface via the commuting
    interpolation (dual ``Set``), i.e. C is the discrete surface rot; then
    ``L_psi = mu0 C^T SL C``.  Validated -- torus L matches analytic to
    -0.6 %; psi=z cylinder gives a solenoid whose energy matches the Nagaoka
    coefficient (0.78 at 2R/L = 0.6).

    DENSE by nature: the inductance is a fully-coupled integral operator
    (every current element couples via 1/r), so L_psi is dense N x N -- unlike
    the sparse l2/h1 seminorms -- which caps min-inductance design at moderate
    N (the BEM SL assembly + the dense C^T SL C).  For very large meshes use
    h1 (a sparse smoothness proxy) or a future H-matrix SL.

    Caller holds the TaskManager context (panel calc, not a helper module).
    A tiny ``ridge`` lifts the harmless constant-psi null space (both A and
    L_psi annihilate an additive constant in psi); confine abe/on also grounds
    it via R."""
    import math
    from ngsolve import (HDivSurface, GridFunction, specialcf, Cross, grad,
                         ds)
    from ngsolve.bem import LaplaceSL
    from scipy.sparse import coo_matrix
    MU0 = 4.0e-7 * math.pi
    n = specialcf.normal(3)
    fJ = HDivSurface(coil, order=order_J)
    u, v = fJ.TnT()
    # BEM single-layer (dense, 100 % fill); COO extract (ToDense ~2500x slower)
    op = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
    rr, cc, vv = op.mat.COO()
    SL = coo_matrix((np.asarray(vv), (np.asarray(rr), np.asarray(cc))),
                    shape=(op.mat.height, op.mat.width)).toarray()
    # C: K = n x grad(psi_i) interpolated into the RT current space
    gp = GridFunction(fes)
    gJ = GridFunction(fJ)
    C = np.zeros((fJ.ndof, fes.ndof))
    for i in range(fes.ndof):
        gp.vec[:] = 0.0
        gp.vec[i] = 1.0
        gJ.Set(Cross(n, grad(gp).Trace()),
               definedon=coil.Boundaries(".*"), dual=True)
        C[:, i] = gJ.vec.FV().NumPy()
    L = MU0 * (C.T @ SL @ C)
    L = 0.5 * (L + L.T)                                   # symmetrise
    L += ridge * (np.trace(L) / L.shape[0]) * np.eye(L.shape[0])
    return L


def _build_abe_R(coil, fes, order):
    """Abe edge-equipotential node-current-potential constraint (Abe, IEEE
    Trans. Magn., 'A Design Tool/Technique for MRI GC using DUCAS', Appendix:
    eq.6 ``T = R T_IN`` with constraints A-1 / A-3).

    Each PHYSICAL boundary edge component of the coil patch is tied to ONE
    free constant (A-1: no current crosses the edge -> the contours close);
    one component (or one node, on a closed surface) is grounded (A-3: only
    grad(T) carries current, so the overall additive constant is fixed).
    Higher-order DOFs ON the physical boundary are zeroed (T constant along
    the edge).  A 'seam' edge of a periodic CAD face is NOT a physical
    boundary -- it is excluded by the element-adjacency test (a physical
    boundary mesh-edge belongs to exactly ONE surface element; a seam to two).

    Returns (R, n_components)."""
    from collections import defaultdict
    from ngsolve import BND
    nv, ndof = coil.nv, fes.ndof
    ecount = defaultdict(int)
    for el in coil.Elements(BND):
        vs = [v.nr for v in el.vertices]
        m = len(vs)
        for a in range(m):
            i, j = vs[a], vs[(a + 1) % m]
            ecount[(min(i, j), max(i, j))] += 1
    bnd_edges = [e for e, c in ecount.items() if c == 1]   # physical boundary
    adj = defaultdict(list)
    bverts = set()
    for i, j in bnd_edges:
        adj[i].append(j); adj[j].append(i); bverts |= {i, j}
    comp = {}
    cid = 0
    for v in bverts:
        if v in comp:
            continue
        stack = [v]; comp[v] = cid
        while stack:
            u_ = stack.pop()
            for w in adj[u_]:
                if w not in comp:
                    comp[w] = cid; stack.append(w)
        cid += 1
    ncomp = cid

    cols = []
    comp_col = {}
    # one free constant per boundary component, grounding component 0
    for c in range(ncomp):
        comp_col[c] = None if c == 0 else len(cols)
        if c != 0:
            cols.append([])
    grounded_vertex = None
    for vtx in range(nv):
        if vtx in comp:
            cc = comp_col[comp[vtx]]
            if cc is not None:
                cols[cc].append(vtx)
            # else: grounded boundary component -> psi = 0 there (no column)
        else:
            cols.append([vtx])                       # interior vertex
    if ncomp == 0:        # closed surface (no boundary): ground one node
        grounded_vertex = 0
        cols = [c for c in cols if c != [0]]
    # higher-order DOFs: interior keep, physical-boundary zero
    if order > 1:
        from ngsolve import H1
        fd = H1(coil, order=order, definedon=coil.Boundaries(".*"),
                dirichlet_bbnd=".*")
        bmask = ~np.array(fd.FreeDofs())
        for d in range(nv, ndof):
            if not bmask[d]:
                cols.append([d])                     # interior higher-order
            # else: boundary higher-order -> zeroed (no column)
    import scipy.sparse as sp
    rows = [d for ds_ in cols for d in ds_]
    cidx = [k for k, ds_ in enumerate(cols) for _ in ds_]
    R = sp.csr_matrix((np.ones(len(rows)), (rows, cidx)),
                      shape=(ndof, len(cols)))
    return R, ncomp


def _dof_reduction_R(args, coil, fes):
    """Independent-DOF -> all-DOF map R (Abe eq.6).  off/on: R = I[:, free] so
    the result is numerically identical to plain free-DOF selection.  abe:
    Abe's edge-equipotential constraint (handles gradient AND solenoid)."""
    from ngsolve import H1                                   # noqa: F401 (abe path)
    if getattr(args, "confine", "off") == "abe":
        R, _ = _build_abe_R(coil, fes, args.order)
        return R
    import scipy.sparse as sp
    fi = np.where(np.array(fes.FreeDofs()))[0]
    R = sp.csr_matrix((np.ones(len(fi)), (fi, np.arange(len(fi)))),
                      shape=(fes.ndof, len(fi)))
    return R


def _element_centroids(mesh):
    """Interior (non-vertex) measure points: element centroids of the eval mesh."""
    vpts = np.array([list(p.point) for p in mesh.vertices])
    cents = []
    for el in mesh.Elements():
        vs = [v.nr for v in el.vertices]
        cents.append(vpts[vs].mean(axis=0))
    return np.array(cents) if cents else vpts


def _peak_current_density(fes, coil, gfu):
    """max |K| = max |grad_s psi|, sampled via an L2 projection of the |grad|
    CF onto the boundary-H1 nodal values (the documented vertex-DOF workaround
    for boundary-defined H1 point eval)."""
    from ngsolve import GridFunction, grad, sqrt as ng_sqrt, InnerProduct
    Kmag = ng_sqrt(InnerProduct(grad(gfu).Trace(), grad(gfu).Trace()))
    gK = GridFunction(fes)
    gK.Set(Kmag, definedon=coil.Boundaries(".*"))
    return float(np.max(np.abs(gK.vec.FV().NumPy())))


def _build_problem(args, eval_scale=1.0):
    """Shared setup: coil FES, eval constraint + measure points, target, A, S.

    ``eval_scale`` (geometry lever) scales the eval-region points about their
    common centroid -- a larger homogeneous volume relative to the fixed coil
    is the dual of the classic former-size lever, and needs no coil re-mesh."""
    from ngsolve import Mesh, H1, specialcf
    from radia.stream_function import aca_tsvd, RegularizedTSVD

    _resolve_target_harmonic(args)          # --target-harmonic -> args.target_cf
    coil = Mesh(args.coil_vol)
    # confine = the stream-function (current-potential) boundary condition.
    #   off : none (current may cross the former edge -> contours can run off
    #         the ends).
    #   on  : psi = 0 on every boundary edge (simple gradient/shim BC; BREAKS
    #         solenoid-type targets where the two ends need different psi).
    #   abe : Abe edge-equipotential -- each physical boundary component gets
    #         ONE free constant + one ground (eq.6 R).  Closes the contours
    #         AND works for gradient and solenoid; strictly the best.
    # All three are expressed as one DOF-reduction matrix R (Abe eq.6
    # ``T = R T_IN``); off/on reduce to plain free-DOF column selection.
    confine = getattr(args, "confine", "off")
    diri = {"dirichlet_bbnd": ".*"} if confine == "on" else {}
    fes = H1(coil, order=args.order, definedon=coil.Boundaries(".*"), **diri)
    n = specialcf.normal(3)
    R = _dof_reduction_R(args, coil, fes)

    target_cf, tdim = _parse_target_cf(args.target_cf)
    if tdim not in (1, 3):
        raise ValueError(f"target-cf dim {tdim}; expected 1 (Bz) or 3 (B)")
    vector_b = (tdim == 3)
    comps = [0, 1, 2] if vector_b else [2]

    evalm = Mesh(args.eval_vol)
    cpts = np.array([list(p.point) for p in evalm.vertices])    # constraints
    if args.eval_max and len(cpts) > args.eval_max:
        cpts = cpts[np.linspace(0, len(cpts) - 1, args.eval_max).astype(int)]
    mpts = _element_centroids(evalm)                            # interior measure
    if args.eval_max and len(mpts) > args.eval_max:
        mpts = mpts[np.linspace(0, len(mpts) - 1, args.eval_max).astype(int)]

    if eval_scale != 1.0:
        c0 = cpts.mean(axis=0)
        cpts = c0 + (cpts - c0) * eval_scale
        mpts = c0 + (mpts - c0) * eval_scale

    Bc = _eval_target_np(args.target_cf, cpts, vector_b)
    Bm = _eval_target_np(args.target_cf, mpts, vector_b)
    Ac = _assemble_biot_savart(fes, n, cpts, comps)
    Am = _assemble_biot_savart(fes, n, mpts, comps)
    if getattr(args, "iron_vol", None):
        # MATERIAL-AWARE kernel: add the iron reaction field (reduced scalar
        # potential on a Kelvin open-boundary FE mesh) to the free-space Biot-
        # Savart operator, so the WHOLE downstream design (TSVD/ridge/pareto/
        # manufacture/distort/chain) sees the iron yoke/shield/core.  One Kelvin
        # factorisation serves both constraint (cpts) and measure (mpts) points.
        obs_all = np.vstack([cpts, mpts])
        A_react = _assemble_iron_reaction(
            fes, n, obs_all, comps, args.iron_vol, float(args.mu_r),
            iron_mat=getattr(args, "iron_mat", "iron"),
            exact_source=getattr(args, "iron_exact_source", False),
            quad_order=getattr(args, "iron_quad_order", None))
        nrc = len(cpts) * len(comps)
        Ac = Ac + A_react[:nrc]
        Am = Am + A_react[nrc:]
    if args.regularize == "inductance":
        # DENSE BEM (SL ~ (3N)^2, plus an N-iteration Set loop): O(N^2) memory.
        # No-Fallback: cap with a clear message instead of OOM/hanging on a big
        # mesh (the sparse l2/h1 seminorms scale to N~15k; this does not).
        if fes.ndof > 5000:
            raise ValueError(
                f"--regularize inductance is DENSE (BEM self-inductance, "
                f"O(N^2) memory + an N-column assembly loop); N={fes.ndof} "
                f"surface DOFs is too large (cap 5000). Use --regularize h1 "
                f"(sparse, the smoothness proxy) for large meshes, or coarsen "
                f"the mesh / lower --order.")
        # physical min-stored-energy objective: DENSE BEM self-inductance L_psi
        Lf = _inductance_seminorm(coil, fes)    # dense N x N (mu0 C^T SL C)
        # fold to the independent space WITHOUT a dense N x N intermediate of
        # R: (R^T Lf) is sparse@dense=dense (n_free x N); (.T) avoids dense@
        # sparse; final sparse@dense -> dense (n_free x n_free) SPD core.
        S = np.asarray(R.T @ (R.T @ Lf).T)      # = R^T Lf R (Lf symmetric)
    else:
        Sf = _seminorm(fes, args.regularize)    # full-DOF SPARSE seminorm
        S = (R.T @ Sf @ R).tocsr()              # reduce to independent (sparse)
    # Ac is dense (M x N, M small); R is sparse (N x n_free).  (R^T Ac^T)^T
    # keeps R sparse and yields the dense M x n_free design matrix.
    Af = np.asarray((R.T @ Ac.T).T)             # rows x n_independent (dense)
    Am_R = np.asarray((R.T @ Am.T).T)
    n_free = R.shape[1]
    base = aca_tsvd(Ac.shape[0], n_free, lambda i, j: float(Af[i, j]),
                    modes=Ac.shape[0], kmax=min(Ac.shape[0], n_free),
                    aca_eps=1e-10)
    reg = RegularizedTSVD.from_stiffness(base, S)
    # mean diag(Af^T Af) = mean of column norms^2 -- WITHOUT forming the
    # n_free x n_free dense Af^T Af (the other O(N^2) wall).
    md = float(np.mean(np.sum(Af * Af, axis=0)))
    return dict(coil=coil, fes=fes, R=R, n_free=n_free, n=n, vector_b=vector_b,
                comps=comps, Ac=Ac, Af=Af, Am=Am_R, Bc=Bc, Bm=Bm, reg=reg,
                base=base, target_cf=target_cf, evalm=evalm, cpts=cpts,
                mpts=mpts, eval_scale=eval_scale, confine=confine, S=S,
                md=md, n_constraint=len(cpts), n_measure=len(mpts))


# ==========================================================================
# ACTIVE SHIELDING: primary (inner) + shield (outer) coil design.
#
# Turner (1986) / Mansfield & Chapman (1986): two concentric cylindrical
# surfaces designed JOINTLY so that
#   (1) primary + shield together produce the target field inside the DSV
#   (2) primary + shield together cancel the stray field OUTSIDE the shield
#
# Implementation: fold both constraints into one stacked LS system
#
#   A_stack = [ A_dsv_p  | A_dsv_s  ]   (M_c rows -- DSV target)
#             [ w*A_ext_p| w*A_ext_s ]   (M_e rows -- stray null)
#   B_stack = [B_target; zeros(M_e)]
#
# psi_f = [psi_primary_free | psi_shield_free] (concatenated free DOFs).
# RegularizedTSVD handles the stacked system identically to the single-coil
# case; after solving, psi_f is split back into primary and shield parts.
#
# HONEST SCOPE (measured): this is a POINT-SAMPLED external nulling, NOT the
# analytic Fourier-mode cancellation of a textbook shielded gradient coil.
# The reported ``stray_rms`` is the combined field at INDEPENDENT external
# measure points (``e_mpts``, the external-mesh element centroids -- NOT the
# constraint vertices); the ``stray_fit_rms`` at the constraint points is the
# circular fit residual and is NOT the shielding quality.
#
# COVERAGE IS EVERYTHING: the shield-eval-vol must COVER the exterior the user
# cares about, not a thin slice.  With a LARGE external region (a full-length
# annular shell spanning the coil) the shielding is genuine and strong --
# measured ~86x (39 dB) overall and 100-1700x in the well-covered far region,
# DSV homogeneity preserved.  With a TOO-SMALL region (a thin mid-plane shell)
# the point-sampled nulling OVERFITS that slice (~4x locally) and can make the
# field WORSE at larger z.  The un-sampled GAP between the shield and the null
# region is essentially unshielded (~2-3x).  An analytic external-multipole-
# moment constraint would remove the sampling dependence -- a documented next
# step, not done here.
# ==========================================================================

def _build_shielded_problem(args, eval_scale=1.0):
    """Active-shielding design (primary + shield coil, Turner method).

    Requires ``args.shield_vol`` (outer shield surface mesh) and
    ``args.shield_eval_vol`` (external region where stray field = 0).
    Optional: ``args.shield_weight`` (default 1.0) -- weight on the
    stray-field null-constraint rows relative to the DSV target rows.

    Returns a dict compatible with ``_solve_and_metrics``.  Extra keys
    ``is_shielded``, ``n_primary``, ``n_shield``, etc. are used by
    ``run_design`` for split reporting."""
    import copy
    import scipy.sparse as sp
    from ngsolve import Mesh, H1, specialcf
    from radia.stream_function import aca_tsvd, RegularizedTSVD

    _resolve_target_harmonic(args)
    target_cf, tdim = _parse_target_cf(args.target_cf)
    if tdim not in (1, 3):
        raise ValueError(f"target-cf dim {tdim}; expected 1 (Bz) or 3 (B)")
    vector_b = (tdim == 3)
    comps = [0, 1, 2] if vector_b else [2]
    confine = getattr(args, "confine", "off")
    diri = {"dirichlet_bbnd": ".*"} if confine == "on" else {}

    if args.regularize == "inductance":
        raise ValueError("--regularize inductance is not supported with "
                         "--shield-vol. Use h1 (smooth) or l2 (min-norm).")
    if not getattr(args, "shield_eval_vol", None):
        raise ValueError("--shield-vol requires --shield-eval-vol (the "
                         "external region where the stray field must vanish).")

    # --- primary coil (inner surface) ---
    coil_p = Mesh(args.coil_vol)
    fes_p = H1(coil_p, order=args.order,
               definedon=coil_p.Boundaries(".*"), **diri)
    n_p = specialcf.normal(3)
    R_p = _dof_reduction_R(args, coil_p, fes_p)

    # --- shield coil (outer surface) ---
    coil_s = Mesh(args.shield_vol)
    fes_s = H1(coil_s, order=args.order,
               definedon=coil_s.Boundaries(".*"), **diri)
    n_s = specialcf.normal(3)
    # _dof_reduction_R uses args.confine and args.order, not args.coil_vol
    R_s = _dof_reduction_R(args, coil_s, fes_s)

    # --- DSV evaluation points (inside primary; TARGET region) ---
    evalm = Mesh(args.eval_vol)
    cpts = np.array([list(p.point) for p in evalm.vertices])
    if args.eval_max and len(cpts) > args.eval_max:
        cpts = cpts[np.linspace(0, len(cpts) - 1, args.eval_max).astype(int)]
    mpts = _element_centroids(evalm)
    if args.eval_max and len(mpts) > args.eval_max:
        mpts = mpts[np.linspace(0, len(mpts) - 1, args.eval_max).astype(int)]
    if eval_scale != 1.0:
        c0 = cpts.mean(axis=0)
        cpts = c0 + (cpts - c0) * eval_scale
        mpts = c0 + (mpts - c0) * eval_scale

    # --- external points (outside shield; STRAY FIELD = 0) ---
    # epts (vertices)   = the null CONSTRAINTS that drive the solve.
    # e_mpts (centroids)= INDEPENDENT measure points (NOT in the constraint
    #   set) used to report the HONEST generalising stray field -- the exact
    #   analogue of cpts (DSV constraints) vs mpts (DSV measure).  Reporting
    #   the stray at the constraint points themselves would be circular (an
    #   exact-fit underdetermined solve nulls its own constraints to ~0).
    extm = Mesh(args.shield_eval_vol)
    epts = np.array([list(p.point) for p in extm.vertices])
    if args.eval_max and len(epts) > args.eval_max:
        epts = epts[np.linspace(0, len(epts) - 1, args.eval_max).astype(int)]
    e_mpts = _element_centroids(extm)
    if args.eval_max and len(e_mpts) > args.eval_max:
        e_mpts = e_mpts[np.linspace(0, len(e_mpts) - 1,
                                    args.eval_max).astype(int)]

    Bc = _eval_target_np(args.target_cf, cpts, vector_b)
    Bm = _eval_target_np(args.target_cf, mpts, vector_b)

    # --- Biot-Savart matrices (primary/shield at 4 point sets) ---
    Ac_p = _assemble_biot_savart(fes_p, n_p, cpts, comps)
    Am_p = _assemble_biot_savart(fes_p, n_p, mpts, comps)
    Ae_p = _assemble_biot_savart(fes_p, n_p, epts, comps)
    Aem_p = _assemble_biot_savart(fes_p, n_p, e_mpts, comps)   # ext measure
    Ac_s = _assemble_biot_savart(fes_s, n_s, cpts, comps)
    Am_s = _assemble_biot_savart(fes_s, n_s, mpts, comps)
    Ae_s = _assemble_biot_savart(fes_s, n_s, epts, comps)
    Aem_s = _assemble_biot_savart(fes_s, n_s, e_mpts, comps)   # ext measure

    # R-reduce: (R^T A^T)^T -- keeps R sparse, result is dense M x n_free
    def _rR(A, R):
        return np.asarray((R.T @ A.T).T)

    n_fp, n_fs = R_p.shape[1], R_s.shape[1]
    Ac_pR = _rR(Ac_p, R_p)    # (M_c, n_fp)
    Am_pR = _rR(Am_p, R_p)    # (M_m, n_fp)
    Ae_pR = _rR(Ae_p, R_p)    # (M_e, n_fp)
    Aem_pR = _rR(Aem_p, R_p)  # (M_em, n_fp)
    Ac_sR = _rR(Ac_s, R_s)    # (M_c, n_fs)
    Am_sR = _rR(Am_s, R_s)    # (M_m, n_fs)
    Ae_sR = _rR(Ae_s, R_s)    # (M_e, n_fs)
    Aem_sR = _rR(Aem_s, R_s)  # (M_em, n_fs)

    # Stacked combined matrices (primary | shield columns)
    Ac_f = np.hstack([Ac_pR, Ac_sR])      # (M_c, n_fp+n_fs)
    Am_f = np.hstack([Am_pR, Am_sR])      # (M_m, n_fp+n_fs)  -- DSV homogeneity
    Ae_f = np.hstack([Ae_pR, Ae_sR])      # (M_e, n_fp+n_fs)  -- stray CONSTRAINT
    Aem_f = np.hstack([Aem_pR, Aem_sR])   # (M_em,n_fp+n_fs)  -- stray MEASURE

    w = float(getattr(args, "shield_weight", 1.0))
    n_free = n_fp + n_fs
    A_stack = np.vstack([Ac_f, w * Ae_f])            # (M_c+M_e, n_free)
    B_stack = np.concatenate([Bc, np.zeros(Ae_f.shape[0])])

    # Block-diagonal seminorm: each coil regularized independently
    Sf_p = _seminorm(fes_p, args.regularize)
    Sf_s = _seminorm(fes_s, args.regularize)
    S_p_ind = (R_p.T @ Sf_p @ R_p).tocsr()
    S_s_ind = (R_s.T @ Sf_s @ R_s).tocsr()
    S_block = sp.block_diag([S_p_ind, S_s_ind], format="csr")

    md = float(np.mean(np.sum(A_stack * A_stack, axis=0)))
    base = aca_tsvd(A_stack.shape[0], n_free,
                    lambda i, j: float(A_stack[i, j]),
                    modes=A_stack.shape[0],
                    kmax=min(A_stack.shape[0], n_free),
                    aca_eps=1e-10)
    reg = RegularizedTSVD.from_stiffness(base, S_block)

    return dict(
        # _solve_and_metrics-compatible:
        Af=A_stack, Am=Am_f, Bc=B_stack, Bm=Bm, reg=reg, md=md,
        # shield-specific (run_design uses these for split reporting):
        is_shielded=True,
        n_primary=n_fp, n_shield=n_fs,
        fes=fes_p, shield_fes=fes_s,
        coil=coil_p, shield_coil=coil_s,
        n=n_p, shield_n=n_s,
        R=sp.block_diag([R_p, R_s], format="csr"),  # block R (compat.)
        R_p=R_p, R_s=R_s,
        n_free=n_free, vector_b=vector_b, comps=comps,
        Ac_f=Ac_f, Ae_f=Ae_f, Aem_f=Aem_f, Bc_dsv=Bc,
        shield_weight=w,
        target_cf=args.target_cf, evalm=evalm,
        cpts=cpts, mpts=mpts, epts=epts, e_mpts=e_mpts,
        eval_scale=eval_scale, confine=confine, S=S_block,
        n_constraint=len(cpts), n_measure=len(mpts),
        n_external=len(epts), n_external_measure=len(e_mpts),
    )


def _solve_and_metrics(P, alpha):
    """Solve psi at a given alpha; return (psi_full, fit_rms, homogeneity, peak)."""
    Bc, Bm, Af, Am, reg, md = (P["Bc"], P["Bm"], P["Af"], P["Am"],
                               P["reg"], P["md"])
    psi_f = reg.solve(Bc, alpha=alpha * md) if alpha > 0 else reg.solve(Bc)
    fit_rms = float(np.linalg.norm(Af @ psi_f - Bc) / (np.linalg.norm(Bc) + 1e-30))
    # true homogeneity = field error at the INTERIOR measure points (not
    # constraints) -> captures the between-constraint deviation.  Af/Am are
    # already reduced to the independent-DOF space (Af = Ac R, Am = Am_full R).
    homo = float(np.linalg.norm(Am @ psi_f - Bm)
                 / (np.linalg.norm(Bm) + 1e-30))
    return psi_f, fit_rms, homo


# ==========================================================================
# Pareto lever bodies (design stage): linf (L-inf minimax) + geometry sweep.
# (alpha = the Tikhonov L-curve in _solve_and_metrics; sheet-metal is a
# manufacture-stage WIRE operation -- see run_manufacture --distort.)
# ==========================================================================
def _per_vertex_grad_norm(fes, coil, gfu):
    """Per-vertex |grad_s psi| via the documented boundary-H1 nodal Set
    workaround (point eval on a boundary-defined H1 returns 0)."""
    from ngsolve import GridFunction, grad, sqrt as ng_sqrt, InnerProduct
    Kmag = ng_sqrt(InnerProduct(grad(gfu).Trace(), grad(gfu).Trace()))
    gK = GridFunction(fes)
    gK.Set(Kmag, definedon=coil.Boundaries(".*"))
    return np.abs(gK.vec.FV().NumPy()[:coil.nv].copy())


def _weighted_surface_stiffness(fes, coil, w_vertex):
    """``int w(x) grad_s u . grad_s v ds`` over ALL DOFs; w from vertex values.
    Caller reduces to the independent space via R.  Caller holds the
    TaskManager context (this is a panel calc, not a helper module)."""
    from ngsolve import BilinearForm, GridFunction, H1, grad, ds
    fes_w = H1(coil, order=1, definedon=coil.Boundaries(".*"))
    gw = GridFunction(fes_w)
    gw.vec.FV().NumPy()[:coil.nv] = w_vertex
    import scipy.sparse as sp
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += gw * grad(u).Trace() * grad(v).Trace() * ds
    a += 1.0e-12 * u * v * ds
    a.Assemble()
    r, c, val = a.mat.COO()
    return sp.csr_matrix((np.asarray(val), (np.asarray(r), np.asarray(c))),
                         shape=(a.mat.height, a.mat.width))


def _linf_irls_front(P, alpha, n_iter, weight_ratio_max=20.0, damping=0.5,
                     hot_quantile=0.85, hot_boost=3.0):
    """L-inf (minimax peak current) lever via bounded-ratio active-set IRLS.

    Re-uses the ONE ACA+TSVD factorisation (``P['base']``); each iter rebuilds
    only the weighted surface-H1 stiffness and refolds RegularizedTSVD.  The
    returned 'front' is the IRLS trajectory: peak |grad psi| drops over
    iterations while the homogeneity stays ~constant (exact data fit), i.e.
    the lever pushes the design point DOWN the peak-current axis."""
    from ngsolve import GridFunction
    from radia.stream_function import RegularizedTSVD
    fes, coil, R = P["fes"], P["coil"], P["R"]
    Bc, Bm, Am, md, base = P["Bc"], P["Bm"], P["Am"], P["md"], P["base"]
    nv = coil.nv
    w = np.ones(nv)
    front = []
    for k in range(n_iter):
        Sw = R.T @ _weighted_surface_stiffness(fes, coil, w) @ R
        reg = RegularizedTSVD.from_stiffness(base, Sw)
        psi_f = reg.solve(Bc, alpha=alpha * md) if alpha > 0 else reg.solve(Bc)
        psi = R @ psi_f
        gfu = GridFunction(fes); gfu.vec.FV().NumPy()[:] = psi
        gnorm = _per_vertex_grad_norm(fes, coil, gfu)
        peak = float(np.max(gnorm))
        homo = float(np.linalg.norm(Am @ psi_f - Bm)
                     / (np.linalg.norm(Bm) + 1e-30))
        front.append({"iter": int(k), "homogeneity_rms": homo,
                      "peak_J": peak})
        # bounded-ratio reweight + active-set push + renormalise + damping
        rel = np.clip(gnorm / (peak + 1e-30), 1.0 / weight_ratio_max, 1.0)
        w_new = rel ** 2.0
        thresh = float(np.quantile(gnorm, hot_quantile))
        w_new = np.where(gnorm >= thresh, w_new * hot_boost, w_new)
        w_new = w_new / (np.mean(w_new) + 1e-30)
        w = (1.0 - damping) * w + damping * w_new
    return front


# ==========================================================================
# Manufacture stage: psi iso-contours -> single-stroke wire -> (sheet-metal)
# -> CAD(STEP) -> PEEC.  All GENERAL on a loaded surface mesh (marching
# triangles), not tied to a cylinder / plane parametrisation.
# ==========================================================================
def _surface_triangles(coil):
    """Vertex-index triangles of the coil surface (BND elements; fan-split
    any quads)."""
    from ngsolve import BND
    tris = []
    for el in coil.Elements(BND):
        vs = [v.nr for v in el.vertices]
        if len(vs) == 3:
            tris.append((vs[0], vs[1], vs[2]))
        elif len(vs) == 4:
            tris.append((vs[0], vs[1], vs[2]))
            tris.append((vs[0], vs[2], vs[3]))
    return tris


def _char_len(verts, tris):
    ls = [float(np.linalg.norm(verts[a] - verts[b]))
          for (a, b, c) in tris[:4000]]
    return float(np.median(ls)) if ls else 1.0


def _contour_levels(psi_v, n):
    lo, hi = float(np.min(psi_v)), float(np.max(psi_v))
    if hi - lo < 1e-30:
        return []
    return list(np.linspace(lo, hi, n + 2)[1:-1])


def _contour_segments(verts, fv, tris, level):
    """Marching-triangles iso-segments of (fv - level) == 0 on the surface."""
    segs = []
    for (a, b, c) in tris:
        fa, fb, fc = fv[a] - level, fv[b] - level, fv[c] - level
        pts = []
        for i, j, fi_, fj_ in ((a, b, fa, fb), (b, c, fb, fc), (c, a, fc, fa)):
            if (fi_ < 0.0) != (fj_ < 0.0):
                t = fi_ / (fi_ - fj_)
                pts.append(verts[i] + t * (verts[j] - verts[i]))
        if len(pts) == 2:
            segs.append((pts[0], pts[1]))
    return segs


def _ref_subtriangulation(sub):
    """Barycentric sub-lattice of the unit triangle + its micro-triangles."""
    idx = {}
    refpts = []
    for i in range(sub + 1):
        for j in range(sub + 1 - i):
            idx[(i, j)] = len(refpts)
            refpts.append((i / sub, j / sub))
    subtris = []
    for i in range(sub):
        for j in range(sub - i):
            subtris.append((idx[(i, j)], idx[(i + 1, j)], idx[(i, j + 1)]))
            if i + j < sub - 1:
                subtris.append((idx[(i + 1, j)], idx[(i + 1, j + 1)],
                                idx[(i, j + 1)]))
    return refpts, subtris


def _contour_segments_hp(coil, gfu, levels, sub):
    """Order-p ('high-p') iso-contours: subdivide each surface triangle into
    sub^2 micro-triangles in reference coordinates and evaluate the FULL-order
    psi GridFunction at every micro-vertex via ``GetTrafo`` (so the contour
    follows the curved order-2/3 psi WITHIN each element, which vertex-linear
    marching misses -- this is the FE analogue of the analytical flux-line
    trace, Hirahatake/Noguchi/Igarashi/Yamashita).  sub=1 reduces exactly to
    vertex-linear marching.  Returns {level: [segments]}."""
    from ngsolve import BND, IntegrationRule
    refpts, subtris = _ref_subtriangulation(sub)
    out = {L: [] for L in levels}
    n_nontri = 0
    for el in coil.Elements(BND):
        if len([v for v in el.vertices]) != 3:
            n_nontri += 1                    # quad/other: counted, raised below
            continue
        tr = coil.GetTrafo(el)
        P = np.empty((len(refpts), 3))
        F = np.empty(len(refpts))
        for k, (u, v) in enumerate(refpts):
            mip = tr(IntegrationRule([(u, v)], [1.0])[0])
            P[k] = [mip.point[0], mip.point[1], mip.point[2]]
            F[k] = float(gfu(mip))
        for L in levels:
            for (a, b, c) in subtris:
                fa, fb, fc = F[a] - L, F[b] - L, F[c] - L
                pts = []
                for i, j, fi_, fj_ in ((a, b, fa, fb), (b, c, fb, fc),
                                       (c, a, fc, fa)):
                    if (fi_ < 0.0) != (fj_ < 0.0):
                        t = fi_ / (fi_ - fj_)
                        pts.append(P[i] + t * (P[j] - P[i]))
                if len(pts) == 2:
                    out[L].append((pts[0], pts[1]))
    if n_nontri:                             # No-Fallback: do not silently drop
        raise ValueError(
            f"--contour-sub > 1 (high-order contour) needs a TRIANGULAR surface "
            f"mesh; found {n_nontri} non-triangular (quad) BND element(s). Use "
            f"--contour-sub 1 (vertex-linear marching handles quads), or remesh "
            f"the coil surface with triangles.")
    return out


def _segs_to_polylines(segs, tol):
    """Chain iso-segments into polylines (loops/arcs) by snapping shared
    endpoints (adjacent triangles emit identical crossing points on a shared
    edge, so exact-match snapping is robust)."""
    coords = []
    idx = {}

    def gid(p):
        key = (int(round(p[0] / tol)), int(round(p[1] / tol)),
               int(round(p[2] / tol)))
        if key not in idx:
            idx[key] = len(coords); coords.append(np.asarray(p, float))
        return idx[key]

    adj = {}
    for a, b in segs:
        ia, ib = gid(a), gid(b)
        if ia == ib:
            continue
        adj.setdefault(ia, []).append(ib)
        adj.setdefault(ib, []).append(ia)
    coords = np.array(coords) if coords else np.zeros((0, 3))
    visited = set()
    polylines = []
    for start in list(adj.keys()):
        for nb in adj[start]:
            e0 = frozenset((start, nb))
            if e0 in visited:
                continue
            visited.add(e0)
            poly = [start, nb]
            cur = nb
            while True:
                nxts = [x for x in adj.get(cur, [])
                        if frozenset((cur, x)) not in visited]
                if not nxts:
                    break
                nxt = nxts[0]
                visited.add(frozenset((cur, nxt)))
                poly.append(nxt); cur = nxt
                if cur == start:
                    break
            if len(poly) >= 2:
                polylines.append(coords[poly])
    return polylines


def _nn_order(cents):
    n = len(cents)
    used = [False] * n
    order = [0]; used[0] = True
    for _ in range(n - 1):
        last = cents[order[-1]]
        best, bd = -1, 1e30
        for j in range(n):
            if used[j]:
                continue
            d = float(np.linalg.norm(cents[j] - last))
            if d < bd:
                bd, best = d, j
        order.append(best); used[best] = True
    return order


def _loop_gap_matrix(loops, sub=14):
    """Pairwise minimum vertex-to-vertex distance between loops -- a cheap
    proxy for the shortest connector that can bridge loop i to loop j (the
    cut is placed at the nearest points).  Each loop is subsampled to <=
    ``sub`` vertices so the O(n_loops^2 * sub^2) cost stays small."""
    pts = []
    for p in loops:
        body = p[:-1] if (len(p) > 1 and
                          np.linalg.norm(p[0] - p[-1]) < 1e-12) else p
        if len(body) > sub:
            body = body[np.linspace(0, len(body) - 1, sub).astype(int)]
        pts.append(np.asarray(body, float))
    n = len(loops)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            g = float(np.min(np.linalg.norm(
                pts[i][:, None, :] - pts[j][None, :, :], axis=2)))
            D[i, j] = D[j, i] = g
    return D


def _two_opt_order(D, order, max_pass=8):
    """2-opt local search shortening the OPEN tour
    ``sum_k D[order[k], order[k+1]]`` (an open single stroke -- no wrap edge).
    Untangles the long inter-loop crossings a greedy nearest-neighbour visit
    leaves behind: starting from the NN order it repeatedly reverses a
    sub-tour whenever doing so shortens the two boundary connectors, so the
    long 'jump across the coil to a far lobe and back' rungs collapse into
    short local hops.  Cheap (n_loops ~ tens) and never lengthens the tour."""
    order = list(order)
    n = len(order)
    if n < 4:
        return order
    improved, npass = True, 0
    while improved and npass < max_pass:
        improved, npass = False, npass + 1
        for i in range(1, n - 1):
            a = order[i - 1]
            for j in range(i + 1, n):
                b, c = order[i], order[j]
                d = order[j + 1] if j + 1 < n else None
                before = D[a, b] + (D[c, d] if d is not None else 0.0)
                after = D[a, c] + (D[b, d] if d is not None else 0.0)
                if after + 1e-12 < before:
                    order[i:j + 1] = order[i:j + 1][::-1]
                    improved = True
    return order


def _single_stroke_chain(polylines):
    """Nearest-neighbour single stroke: visit contour loops in NN order, each
    closed into a full turn, opened at the point nearest the previous loop's
    end (the field-aware Kuijpers chain is the planar/cylinder demo
    refinement).  Returns (chain_xyz, n_connectors, total_connector_len)."""
    if not polylines:
        return np.zeros((0, 3)), 0, 0.0
    cents = [p.mean(axis=0) for p in polylines]
    order = _nn_order(cents)
    chain, conns, conn_len, last = [], 0, 0.0, None
    for oi in order:
        poly = polylines[oi]
        if last is not None:
            d = np.linalg.norm(poly - last, axis=1)
            s = int(np.argmin(d))
            poly = np.roll(poly, -s, axis=0)
            conn_len += float(np.linalg.norm(poly[0] - last)); conns += 1
        poly = np.vstack([poly, poly[0]])          # close the turn
        chain.append(poly); last = poly[-1]
    return np.vstack(chain), conns, conn_len


def _dedupe(path, tol=1e-10):
    if len(path) == 0:
        return path
    keep = [path[0]]
    for p in path[1:]:
        if np.linalg.norm(p - keep[-1]) > tol:
            keep.append(p)
    return np.array(keep)


def _field_aware_chain(loops, signs, mpts, target_flat, vector_b,
                       n_cut=None, passes=None):
    """Field-aware single-stroke chain (geometry-independent): visit the
    sign-aligned contour loops in a good order, then choose each loop's
    ENTRY/EXIT point (the cut) by coordinate descent to minimise the full
    one-current WIRE field error ``min_I ||I*(loops + connectors) - B||``.

    Cut-opt RESOLUTION (``n_cut`` candidate cuts/loop, ``passes`` descent
    rounds) is AUTO-scaled with the loop count when not given: more loops ->
    more connectors -> the joint cut-opt needs more candidates/rounds to drive
    their net field to zero.  This matters most for CLUSTERED levels (e.g.
    --optimize-levels, which packs contours to sharpen the loop-set field): at
    the under-resolved old default (n_cut=24, passes=4) a 45-loop clustered Z2
    wire delivered ~0.20 rms vs a 0.0098 loop-set floor; auto-scaling to
    (~50, ~8) recovers it to ~0.067 (3x).  More ``passes`` at a fixed ``n_cut``
    monotonically lower the cut-opt objective (coordinate descent re-scans the
    same grid); a higher ``n_cut`` samples a FINER but not strictly-nested cut
    grid, so it is a STRONG NET improvement (large on clustered/multi-region
    coils) but NOT a per-geometry guarantee -- a marginal case can tie or shift
    a few %.  Override with explicit n_cut/passes (CLI --chain-ncut/--chain-passes).

    The loop field is fixed; only the inter-loop connectors (the 'rungs')
    depend on the cuts.  Minimising the TOTAL error (not the connectors in
    isolation) is what makes this correct on BOTH closed contours (residual
    ~0 -> it drives the connectors' net field to zero, reaching the
    separate-turns floor with NO --distort) and open contours (large residual
    from the rim-closing chords -> it at least uses the connectors to cancel
    part of it).

    Visit ORDER: the same wire-error objective is minimised over a small set
    of candidate orders -- a nearest-neighbour seed and a 2-opt-shortened
    variant (``_two_opt_order``) -- keeping whichever the cut-opt drives
    lowest.  The 2-opt minimises connector LENGTH, which gives the cut-opt a
    tidier, shorter-rung start and HELPS most cases (LAB Gx: +9..+70 %), but
    a length-optimal reorder can break the rungs' symmetric stray-field
    cancellation and HURT others (LAB Gx abe nl=16: -78 %) -- exactly the
    documented 'shorter rungs != better field' trap (single_stroke.md).
    Selecting the lower-wire-error order makes the 2-opt strictly an upside:
    GUARANTEED never worse than nearest-neighbour, captures its gains where
    they are real.  (The pure-NN ``--chain nn`` baseline keeps NN order with
    no cut-opt, so it never sees the 2-opt -- which would only hurt it.)

    Returns (chain_xyz, n_connectors, connector_length)."""
    if not loops:
        return np.zeros((0, 3)), 0, 0.0
    n = len(loops)
    if n_cut is None:               # auto: ~1.1 cuts per loop, capped [24, 64]
        n_cut = int(min(64, max(24, round(1.1 * n))))
    if passes is None:              # auto: grow descent rounds with loop count
        passes = int(min(8, max(4, n // 6 + 4)))
    L0 = [_close_loop(p if s >= 0.0 else p[::-1])
          for p, s in zip(loops, signs)]
    ncomp = len(target_flat)
    den = float(np.linalg.norm(target_flat)) + 1e-30

    def fld(a, b):
        B = _segment_field_B(np.array([a, b]), mpts, 1.0)
        return B.reshape(-1) if vector_b else B[:, 2]

    def loop_fld(p):
        B = _segment_field_B(p, mpts, 1.0)
        return B.reshape(-1) if vector_b else B[:, 2]

    def rot(p, c):
        body = p[:-1]
        r = np.roll(body, -c, axis=0)
        return np.vstack([r, r[0]])

    lf_cache = [loop_fld(p) for p in L0]   # per-loop field (order-independent)

    def optimize(order):
        """Cut coordinate-descent for ONE visit order; returns
        (chain, connector_length, final_wire_err)."""
        L = [L0[i] for i in order]
        loops_field = np.sum([lf_cache[i] for i in order], axis=0)
        cuts = [0] * len(L)
        conn_f = [np.zeros(ncomp) for _ in range(max(0, len(L) - 1))]

        def set_conn(i):
            if 0 <= i < len(L) - 1:
                conn_f[i] = fld(rot(L[i], cuts[i])[0],
                                rot(L[i + 1], cuts[i + 1])[0])

        def wire_err():
            ct = loops_field + (np.sum(conn_f, axis=0) if conn_f else 0.0)
            I = _best_fit_current(ct, target_flat)
            return float(np.linalg.norm(I * ct - target_flat) / den)

        for i in range(len(L) - 1):
            set_conn(i)
        for _ in range(passes):
            for i in range(len(L)):
                nb = len(L[i]) - 1
                if nb < 2:
                    continue
                step = max(1, nb // n_cut)
                best, beste = cuts[i], None
                for c in range(0, nb, step):
                    cuts[i] = c
                    set_conn(i - 1)
                    set_conn(i)
                    e = wire_err()
                    if beste is None or e < beste:
                        beste, best = e, c
                cuts[i] = best
                set_conn(i - 1)
                set_conn(i)
        chain = _dedupe(np.vstack([rot(L[i], cuts[i])
                                   for i in range(len(L))]))
        clen = sum(float(np.linalg.norm(rot(L[i + 1], cuts[i + 1])[0]
                                        - rot(L[i], cuts[i])[0]))
                   for i in range(len(L) - 1))
        return chain, clen, wire_err()

    nn = _nn_order([p.mean(axis=0) for p in L0])
    cands = [nn]
    two = _two_opt_order(_loop_gap_matrix(L0), list(nn))
    if two != nn:
        cands.append(two)
    # (A lobe-aware order -- group loops by current sense, cross at the saddle --
    # was tried to shorten the long inter-lobe rung; it does NOT help: the saddle
    # sits at the DSV centre, so a SHORT crossover there contaminates the target
    # MORE than a long one routed away from it.  The "shorter rungs != better
    # field" trap; the cut-opt already finds the field-best of {nn, 2-opt}.)
    best = None
    for order in cands:
        chain, clen, err = optimize(order)
        if best is None or err < best[2]:
            best = (chain, clen, err)
    chain, clen, _ = best
    return chain, max(0, len(L0) - 1), clen


def _segment_field_B(path, obs, current=1.0):
    """Biot-Savart B (T) of a polyline of straight segments at obs points.
    General finite-wire formula (verified vs infinite-wire mu0 I/2 pi d)."""
    obs = np.asarray(obs, float)
    B = np.zeros((len(obs), 3))
    for k in range(len(path) - 1):
        p1, p2 = path[k], path[k + 1]
        e = p2 - p1
        L = float(np.linalg.norm(e))
        if L < 1e-12:
            continue
        e = e / L
        r1 = obs - p1; r2 = obs - p2
        n1 = np.linalg.norm(r1, axis=1); n2 = np.linalg.norm(r2, axis=1)
        c1 = (r1 @ e) / np.maximum(n1, 1e-30)
        c2 = (r2 @ e) / np.maximum(n2, 1e-30)
        rp = r1 - np.outer(r1 @ e, e)
        d = np.linalg.norm(rp, axis=1)
        phi = np.cross(np.tile(e, (len(obs), 1)), rp)
        pn = np.linalg.norm(phi, axis=1)
        good = (d > 1e-12) & (pn > 1e-30)
        phi_hat = np.zeros_like(phi)
        phi_hat[good] = phi[good] / pn[good, None]
        Bmag = np.zeros(len(obs))
        Bmag[good] = _MU0_4PI * current * (c1[good] - c2[good]) / d[good]
        B += Bmag[:, None] * phi_hat
    return B


def _best_fit_current(field_flat, target_flat):
    d = float(field_flat @ field_flat)
    return float(field_flat @ target_flat / d) if d > 0.0 else 0.0


def _wire_field_flat(path, mpts, vector_b, current=1.0):
    Bm = _segment_field_B(path, mpts, current)
    return Bm.reshape(-1) if vector_b else Bm[:, 2]


def _close_loop(p):
    if len(p) >= 2 and np.linalg.norm(p[0] - p[-1]) > 1e-12:
        return np.vstack([p, p[0]])
    return p


def _tri_grad_K(verts, psi_v, tris):
    """Per-triangle centroid + surface-current direction K = n_hat x grad_s(psi).
    The constant in-plane gradient of the linear psi is solved in an orthonormal
    tangent basis (robust on a curved surface in 3D)."""
    cen = np.empty((len(tris), 3)); K = np.empty((len(tris), 3)); m = 0
    for (a, b, c) in tris:
        pa = verts[a]; e1 = verts[b] - pa; e2 = verts[c] - pa
        nrm = np.cross(e1, e2); A2 = float(np.linalg.norm(nrm))
        l1 = float(np.linalg.norm(e1))
        if A2 < 1e-30 or l1 < 1e-30:
            continue
        t1 = e1 / l1; nhat = nrm / A2; t2 = np.cross(nhat, t1)
        M = np.array([[e1 @ t1, e1 @ t2], [e2 @ t1, e2 @ t2]])
        try:
            ab = np.linalg.solve(M, [psi_v[b] - psi_v[a], psi_v[c] - psi_v[a]])
        except np.linalg.LinAlgError:
            continue
        cen[m] = (pa + verts[b] + verts[c]) / 3.0
        K[m] = np.cross(nhat, ab[0] * t1 + ab[1] * t2)
        m += 1
    return cen[:m], K[:m]


def _build_gradpsi_kdt(verts, psi_v, tris):
    """KDTree of triangle centroids + their K = n_hat x grad_s(psi), built ONCE
    per psi so contour loops can be oriented by the surface-current winding."""
    from scipy.spatial import cKDTree
    cen, K = _tri_grad_K(verts, psi_v, tris)
    return cKDTree(cen), K


def _gradpsi_signs(loops, kdt, tri_K, mpts, vector_b, target_flat,
                   global_align=True):
    """Per-loop unit-current fields + the sign orienting each contour loop by the
    consistent surface-current winding K = n_hat x grad_s(psi) -- NOT per-loop
    field alignment.

    The marching-triangles polyline order is arbitrary, so its winding does NOT
    carry the current sign.  The old s_k = sign(f_k . B) ('flip every loop to
    help the target') recovers it ONLY for a MONOTONE psi (l=1 gradient): there
    every loop wants the same sense.  For a SADDLE psi (l>=2 shim) the natural
    winding has MIXED senses (the single wire reverses through the saddle), and
    sign(f.B) FLATTENS them so the one common current cancels itself (Z2 single-
    wire residual ~88%).  Orienting each loop so its traversal matches K -- a
    KDTree majority vote of (segment direction . K) over the loop -- recovers the
    physical stream-function winding (Z2 ~5.6%) and reduces to sign(f.B) for l=1.
    The K-vote signs are globally consistent BY CONSTRUCTION (K is one global
    vector field), so they may be summed across contour levels directly.  With
    ``global_align`` the OVERALL sense is then flipped once to align with the
    target (cosmetic -- the best-fit current absorbs it -- and it reproduces the
    old l=1 convention).  ``global_align`` MUST be False when this is called
    PER LEVEL (e.g. _optimize_contour_levels): a per-level overall flip would
    re-flatten the inter-level saddle reversals the K-vote just recovered.
    Returns (fields[list], signs[ndarray]); sign = +1 keep / -1 reverse the
    polyline."""
    if not loops:
        return [], np.zeros(0)
    signs = np.empty(len(loops))
    for li, p in enumerate(loops):
        mid = 0.5 * (p[:-1] + p[1:]); d = p[1:] - p[:-1]
        _, j = kdt.query(mid)
        signs[li] = (1.0 if float(np.einsum('ij,ij->i', d, tri_K[j]).sum()) >= 0.0
                     else -1.0)
    fields = [_wire_field_flat(_close_loop(p), mpts, vector_b, 1.0) for p in loops]
    if global_align and target_flat is not None:
        G = np.sum([s * f for s, f in zip(signs, fields)], axis=0)
        if float(G @ target_flat) < 0.0:
            signs = -signs                        # align global sense to target
    return fields, signs


def _optimize_contour_levels(psi_v, verts, tris, nlevels, mpts, target_flat,
                             vector_b, tol, n_rounds=3):
    """Optimize the N contour LEVELS (vs the uniform equal-deltaI default) to
    minimise the EQUAL-CURRENT (single series current) discrete-loop field error
    -- the few-turn quantisation fix.  Coordinate descent on the ORDERED levels
    (each level 1-D bounded between its neighbours), caching the per-level signed
    unit-current field so only the moved level's contour is re-extracted.  Keeps
    clean iso-contours (manufacturable) -- complements --distort's wire bend.
    Returns (levels, rms_uniform, rms_optimised)."""
    from scipy.optimize import minimize_scalar
    lo, hi = float(np.min(psi_v)), float(np.max(psi_v))
    if hi - lo < 1e-30 or nlevels < 1:
        return _contour_levels(psi_v, nlevels), float("nan"), float("nan")
    den = float(np.linalg.norm(target_flat)) + 1e-30
    kdt, tri_K = _build_gradpsi_kdt(verts, psi_v, tris)

    def level_field(lev):
        loops = _segs_to_polylines(_contour_segments(verts, psi_v, tris, lev), tol)
        if not loops:
            return np.zeros(len(target_flat))
        fields, signs = _gradpsi_signs(loops, kdt, tri_K, mpts, vector_b,
                                       target_flat, global_align=False)
        return np.sum([s * f for f, s in zip(fields, signs)], axis=0)

    def resid(F):
        d = float(F @ F)
        if d <= 0.0:
            return 1.0e30
        Iw = float(F @ target_flat) / d          # one common (series) current
        return float(np.linalg.norm(Iw * F - target_flat) / den)

    levels = list(_contour_levels(psi_v, nlevels))
    fk = [level_field(lv) for lv in levels]
    rms0 = resid(np.sum(fk, axis=0))
    eps = 1.0e-6 * (hi - lo)
    for _ in range(int(n_rounds)):
        for k in range(nlevels):
            F_other = (np.sum([fk[j] for j in range(nlevels) if j != k], axis=0)
                       if nlevels > 1 else np.zeros(len(target_flat)))
            a = levels[k - 1] if k > 0 else lo
            b = levels[k + 1] if k < nlevels - 1 else hi
            if b - a < 4.0 * eps:
                continue
            r = minimize_scalar(
                lambda lv: resid(F_other + level_field(lv)),
                bounds=(a + eps, b - eps), method="bounded",
                options={"maxiter": 25})
            levels[k] = float(r.x); fk[k] = level_field(levels[k])
    return levels, rms0, resid(np.sum(fk, axis=0))


def _greedy_coil_turns(cand_fields, target_flat, max_turns, target_rms=None,
                       cand_centroids=None, connector_weight=0.0, coil_diag=1.0):
    """LOW-TURN STRATEGY -- greedy constructive coil ("put up a pin, lay the next
    turn where it helps most").  Matching pursuit with a COMMON series current:
    keep adding the candidate loop (with sign) that best ALIGNS the unit combined
    field G = sum_j s_j f_j with the target (residual = min_I ||I G - target|| =
    the angle between G and target).  MONOTONE by construction -- stop when no
    candidate improves, ``max_turns`` reached, or residual <= ``target_rms``.
    Allows repeats (a re-picked loop = a stacked turn -> more local current).
    Returns (selected [(cand_idx, sign)], trace [{n_turns, rms}]).  This is the
    antidote to the NON-MONOTONIC contour-quantise turn-count search.

    CONNECTOR-AWARE selection (``connector_weight`` > 0, needs
    ``cand_centroids``): bias the SELECTION toward loops geometrically NEAR the
    loops already laid -- score = field_residual + connector_weight * (gap /
    coil_diag), gap = nearest distance from the candidate's centroid to an
    already-selected loop (0 for the first turn or a re-picked/stacked loop).
    This trades a little field optimality for a SHORT single-stroke connector
    ('rung') path -> a manufacturable wire, the fix for the isolated-pin chaining
    penalty (long rungs wreck the single-current homogeneity).  Only field-
    IMPROVING turns are accepted, so the rms-vs-turns trace stays MONOTONE
    regardless of the weight; ``connector_weight=0`` is the pure-field default."""
    den = float(np.linalg.norm(target_flat)) + 1e-30
    F = np.asarray(cand_fields, float)             # (D, M) unit-current fields
    if F.ndim != 2 or F.shape[0] == 0:
        return [], []
    cents = (np.asarray(cand_centroids, float)
             if (cand_centroids is not None and connector_weight > 0.0)
             else None)
    diag = float(coil_diag) + 1e-30
    G = np.zeros(F.shape[1]); selected, trace, sel_cents = [], [], []
    prev = 1.0                                      # empty coil: ||0 - t||/|t| = 1
    for _ in range(int(max_turns)):
        if cents is not None and sel_cents:         # connector gap per candidate
            sc = np.asarray(sel_cents)
            gap = np.min(np.linalg.norm(cents[:, None, :] - sc[None, :, :],
                                        axis=2), axis=1) / diag
        else:
            gap = np.zeros(F.shape[0])
        best = None                                 # (score, field_r, idx, sign)
        for s in (1.0, -1.0):
            Gp = G[None, :] + s * F                 # (D, M)
            d = np.einsum('dm,dm->d', Gp, Gp) + 1e-300
            Iw = (Gp @ target_flat) / d
            r = np.linalg.norm(Iw[:, None] * Gp - target_flat[None, :],
                               axis=1) / den
            score = r + connector_weight * gap
            improving = np.where(r < prev - 1e-12)[0]   # keep trace monotone
            if improving.size == 0:
                continue
            k = int(improving[int(np.argmin(score[improving]))])
            if best is None or score[k] < best[0]:
                best = (float(score[k]), float(r[k]), k, s)
        if best is None:                            # no candidate improves -> stop
            break
        _, fr, k, s = best
        G = G + s * F[k]; selected.append((k, s)); prev = fr
        if cents is not None:
            sel_cents.append(cents[k])
        trace.append({"n_turns": len(selected), "rms": fr})
        if target_rms is not None and fr <= target_rms:
            break
    return selected, trace


def _contour_loop_candidates(psi_v, verts, tris, mpts, vector_b, tol, n_levels=40):
    """Greedy dictionary (A): individual psi iso-contour loops at a DENSE set of
    levels, each a unit-current closed loop + its Biot-Savart field at mpts."""
    loops = []
    for lev in _contour_levels(psi_v, n_levels):
        loops.extend(_segs_to_polylines(
            _contour_segments(verts, psi_v, tris, lev), tol))
    loops = [p for p in loops if len(p) >= 3]
    fields = [_wire_field_flat(_close_loop(p), mpts, vector_b, 1.0) for p in loops]
    return loops, fields


def _pin_loop_candidates(verts, tris, mpts, vector_b, n_pins=120,
                         krings=(1, 2, 3)):
    """Greedy dictionary (B): ON-SURFACE k-RING loops at each PIN (a surface
    vertex) -- "wind a loop on a pin", multi-scale.  For pin v and ring k the
    loop is the mesh vertices at graph distance EXACTLY k from v, ordered by
    angle in v's tangent plane (k=1 = the natural 1-ring; larger k = bigger
    loops -> the "adaptive size" the small flat tangent circles lacked, and the
    loops follow the curved surface, not a flat chord).  Returns (loops, fields)."""
    nv = len(verts)
    adj = [set() for _ in range(nv)]
    for (a, b, c) in tris:
        adj[a].update((b, c)); adj[b].update((a, c)); adj[c].update((a, b))
    vn = np.zeros_like(verts, dtype=float)
    for (a, b, c) in tris:
        nrm = np.cross(verts[b] - verts[a], verts[c] - verts[a])
        vn[a] += nrm; vn[b] += nrm; vn[c] += nrm
    ln = np.linalg.norm(vn, axis=1); ln[ln == 0.0] = 1.0
    vn = vn / ln[:, None]
    pins = np.unique(np.linspace(0, nv - 1, min(n_pins, nv)).astype(int))
    kmax = max(krings)
    loops, fields = [], []
    for v in pins:
        # BFS graph distance from v up to kmax
        dist = {int(v): 0}; frontier = [int(v)]
        for d in range(1, kmax + 1):
            nf = []
            for u in frontier:
                for w in adj[u]:
                    if w not in dist:
                        dist[w] = d; nf.append(w)
            frontier = nf
        nrm = vn[v]
        t1 = np.cross(nrm, [1.0, 0.0, 0.0])
        if np.linalg.norm(t1) < 1e-6:
            t1 = np.cross(nrm, [0.0, 1.0, 0.0])
        t1 = t1 / np.linalg.norm(t1); t2 = np.cross(nrm, t1)
        for k in krings:
            shell = [w for w, dd in dist.items() if dd == k]
            if len(shell) < 3:
                continue
            rel = verts[shell] - verts[v]
            order = np.argsort(np.arctan2(rel @ t2, rel @ t1))
            loop = verts[[shell[i] for i in order]]
            loops.append(loop)
            fields.append(_wire_field_flat(_close_loop(loop), mpts, vector_b, 1.0))
    return loops, fields


def _surface_grad_mag(psi_v, verts, tris):
    """Per-vertex |grad_s psi| (the surface-current magnitude |K| = |n x grad
    psi| up to the constant) -- area-weighted average of each incident triangle's
    constant in-plane gradient.  Solves the 2x2 system for the gradient in an
    orthonormal tangent basis (robust on a curved surface in 3D).  Used to place
    bubble pins with line DENSITY proportional to |K| (the stream-function
    equal-flux packing law)."""
    nv = len(verts)
    acc = np.zeros(nv); wsum = np.zeros(nv)
    for (a, b, c) in tris:
        pa = verts[a]; e1 = verts[b] - pa; e2 = verts[c] - pa
        n = np.cross(e1, e2); area2 = float(np.linalg.norm(n))
        if area2 < 1e-30:
            continue
        ln1 = float(np.linalg.norm(e1))
        if ln1 < 1e-30:
            continue
        t1 = e1 / ln1; t2 = np.cross(n / area2, t1)   # orthonormal in-plane basis
        M = np.array([[e1 @ t1, e1 @ t2], [e2 @ t1, e2 @ t2]])
        rhs = np.array([psi_v[b] - psi_v[a], psi_v[c] - psi_v[a]])
        try:
            ab = np.linalg.solve(M, rhs)
        except np.linalg.LinAlgError:
            continue
        gmag = float(np.hypot(ab[0], ab[1])); area = 0.5 * area2
        for v in (a, b, c):
            acc[v] += gmag * area; wsum[v] += area
    wsum[wsum == 0.0] = 1.0
    return acc / wsum


def _bubble_pin_loops(psi_v, verts, tris, mpts, vector_b, n_target=80,
                      frac=0.7, kmax=6):
    """DENSE BUBBLE-TILING pin coil -- the bubble system (Hirahatake/Noguchi/
    Igarashi/Yamashita) applied to coil DISCRETISATION: drop small circulating
    loops over the surface with DENSITY proportional to |K| = |grad_s psi|
    (clearance radius r_i = k / sqrt(|K|_i), so loops pack where the current is
    strong), and size each loop to its bubble (grow the on-surface k-ring until
    the ring radius ~ frac * r_i).  A many-turn, MANUFACTURABLE alternative to
    long iso-contours: equal-current loops that TILE the surface, dense where
    needed, each a short local loop with a short hop to its neighbour (no long
    cross-coil rungs).  k is bisected so ~``n_target`` loops land.  Returns
    (loops, fields)."""
    nv = len(verts)
    gmag = _surface_grad_mag(psi_v, verts, tris)
    pos = gmag > 0.0
    if not np.any(pos):
        return [], []
    adj = [set() for _ in range(nv)]               # vertex adjacency
    for (a, b, c) in tris:
        adj[a].update((b, c)); adj[b].update((a, c)); adj[c].update((a, b))
    vn = np.zeros_like(verts, dtype=float)         # vertex normals
    for (a, b, c) in tris:
        nrm = np.cross(verts[b] - verts[a], verts[c] - verts[a])
        vn[a] += nrm; vn[b] += nrm; vn[c] += nrm
    ln = np.linalg.norm(vn, axis=1); ln[ln == 0.0] = 1.0
    vn = vn / ln[:, None]

    order = np.argsort(-gmag)                       # strongest |K| first
    diag = float(np.linalg.norm(verts.max(0) - verts.min(0))) + 1e-30

    def pack(k_const):
        seeds = []
        for idx in order:
            b = gmag[idx]
            if b <= 0.0:
                continue
            r = k_const / np.sqrt(b)
            p = verts[idx]
            if all(np.linalg.norm(p - verts[s]) > r for s in seeds):
                seeds.append(int(idx))
        return seeds

    # bisect k so ~n_target seeds land (monotone: bigger k -> fewer seeds)
    bref = float(np.median(gmag[pos]))
    klo = 1e-4 * diag * np.sqrt(bref); khi = 4.0 * diag * np.sqrt(bref)
    k_used = 0.5 * (klo + khi); pins = pack(k_used)
    for _ in range(18):
        k_used = 0.5 * (klo + khi); pins = pack(k_used)
        if len(pins) > n_target:
            klo = k_used
        elif len(pins) < n_target:
            khi = k_used
        else:
            break

    loops, fields = [], []
    for v in pins:
        target_r = frac * k_used / np.sqrt(gmag[v] + 1e-30)
        dist = {v: 0}; frontier = [v]; chosen = None
        for d in range(1, kmax + 1):               # grow rings to the bubble size
            nf = []
            for u in frontier:
                for w_ in adj[u]:
                    if w_ not in dist:
                        dist[w_] = d; nf.append(w_)
            frontier = nf
            if not nf:
                break
            chosen = d
            if float(np.mean(np.linalg.norm(verts[nf] - verts[v],
                                            axis=1))) >= target_r:
                break
        if chosen is None:
            continue
        shell = [w_ for w_, dd in dist.items() if dd == chosen]
        if len(shell) < 3:
            continue
        nrm = vn[v]
        t1 = np.cross(nrm, [1.0, 0.0, 0.0])
        if np.linalg.norm(t1) < 1e-6:
            t1 = np.cross(nrm, [0.0, 1.0, 0.0])
        t1 = t1 / np.linalg.norm(t1); t2 = np.cross(nrm, t1)
        rel = verts[shell] - verts[v]
        ang = np.argsort(np.arctan2(rel @ t2, rel @ t1))
        loop = _close_loop(verts[[shell[i] for i in ang]])   # closed k-ring loop
        loops.append(loop)
        fields.append(_wire_field_flat(loop, mpts, vector_b, 1.0))
    return loops, fields


def _greedy_trace_plot(trace, target, out_path):
    """Turn-budget UI: the greedy rms-vs-turns curve (monotone) -- read off the
    fewest turns meeting a spec.  Marks the --greedy-target line + the turn it is
    first met.  No in-figure title (lab rule -- caption carries it)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ns = [t["n_turns"] for t in trace]; rms = [t["rms"] for t in trace]
    fig, ax = plt.subplots(figsize=(4.0, 3.0))
    ax.plot(ns, rms, "o-", color="C0", label="greedy")
    if target is not None:
        ax.axhline(target, ls="--", color="C3", label="spec")
        met = [n for n, r in zip(ns, rms) if r <= target]
        if met:
            ax.axvline(min(met), ls=":", color="C2",
                       label="min turns = %d" % min(met))
    ax.set_xlabel("turns"); ax.set_ylabel("field residual rms")
    ax.set_ylim(bottom=0.0); ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out_path, dpi=150); plt.close(fig)


def _distort_wire(chain, mpts, target_flat, vector_b, n_grid=3, n_iter=5,
                  lam_rel=0.3, step=0.8, fd=1.0e-5):
    """Sheet-metal (bankin-ho) single-current wire distortion: bend the ONE
    series wire with a smooth 3D trilinear control grid (out-of-surface lever)
    to better hit the target, contour LEVELS fixed.  Gauss-Newton with a
    displacement penalty relative to mean(diag(JtJ)), I_w re-fit each step.
    Returns (chain_best, I_w, rms0, rms_best, max_disp)."""
    chain = np.asarray(chain, float)
    bb_min = chain.min(0); span = np.maximum(chain.max(0) - bb_min, 1e-9)
    g = int(n_grid)
    u = (chain - bb_min) / span * (g - 1)
    i0 = np.clip(np.floor(u).astype(int), 0, g - 2)
    f = u - i0

    def deform(dofs):
        C = dofs.reshape(g, g, g, 3)
        d = np.zeros_like(chain)
        for cx in (0, 1):
            for cy in (0, 1):
                for cz in (0, 1):
                    wgt = ((f[:, 0] if cx else 1 - f[:, 0])
                           * (f[:, 1] if cy else 1 - f[:, 1])
                           * (f[:, 2] if cz else 1 - f[:, 2]))
                    d += wgt[:, None] * C[i0[:, 0] + cx, i0[:, 1] + cy,
                                          i0[:, 2] + cz]
        return chain + d

    def fitI(fu):
        d = float(fu @ fu)
        return float(fu @ target_flat / d) if d else 0.0

    den = float(np.linalg.norm(target_flat)) + 1e-30
    ndof = g * g * g * 3
    p0 = deform(np.zeros(ndof))
    f0 = _wire_field_flat(p0, mpts, vector_b); I0 = fitI(f0)
    rms0 = float(np.linalg.norm(I0 * f0 - target_flat) / den)
    dofs = np.zeros(ndof)
    best = (rms0, p0.copy(), I0)
    for _ in range(int(n_iter)):
        w = deform(dofs)
        fb = _wire_field_flat(w, mpts, vector_b); Ib = fitI(fb)
        base = Ib * fb; r = target_flat - base
        J = np.empty((len(base), ndof))
        for dd in range(ndof):
            d2 = dofs.copy(); d2[dd] += fd
            fu = _wire_field_flat(deform(d2), mpts, vector_b)
            J[:, dd] = (fitI(fu) * fu - base) / fd
        JtJ = J.T @ J
        md = float(np.mean(np.diag(JtJ))) + 1e-30
        dofs = dofs + step * np.linalg.solve(
            JtJ + lam_rel * md * np.eye(ndof),
            J.T @ r - lam_rel * md * dofs)
        w = deform(dofs)
        fu = _wire_field_flat(w, mpts, vector_b); Iu = fitI(fu)
        rms = float(np.linalg.norm(Iu * fu - target_flat) / den)
        if rms < best[0]:
            best = (rms, w.copy(), Iu)
    disp = float(np.max(np.linalg.norm(best[1] - chain, axis=1)))
    return best[1], best[2], rms0, best[0], disp


def _write_step_polylines(polylines, filename):
    """Write contour polylines as OCC wire edges to a STEP file (general)."""
    from netgen.occ import Pnt, Segment, Glue
    edges = []
    for poly in polylines:
        for k in range(len(poly) - 1):
            p, q = poly[k], poly[k + 1]
            if np.linalg.norm(q - p) < 1e-9:
                continue
            edges.append(Segment(Pnt(float(p[0]), float(p[1]), float(p[2])),
                                  Pnt(float(q[0]), float(q[1]), float(q[2]))))
    if not edges:
        raise ValueError("no wire segments to export to STEP")
    Glue(edges).WriteStep(filename)


def _decimate_polyline(poly, step):
    """Keep the endpoints + every point at least ``step`` (m) of straight-line
    distance from the last kept point -- resamples a dense wire polyline to a
    coarser channel spine (fewer CAD primitives) WITHOUT changing its shape.
    ``step <= 0`` returns the polyline unchanged."""
    poly = np.asarray(poly, float)
    if step <= 0.0 or len(poly) < 3:
        return poly
    keep = [0]
    last = poly[0]
    for i in range(1, len(poly) - 1):
        if np.linalg.norm(poly[i] - last) >= step:
            keep.append(i)
            last = poly[i]
    keep.append(len(poly) - 1)
    return poly[keep]


def _write_former_step(chain, wire_diam, clearance, margin, wall, decimate,
                       filename, max_segments=300):
    """Printable FORMER for the delivered wire: a base plate with the wire
    CHANNEL cut out (plate - swept tube), so a water-soluble-support 3D print
    leaves the out-of-surface groove the single-stroke wire threads through.

    The channel = a true boolean UNION (``occ.Fuse``) of per-segment cylinders
    (radius = wire_diam/2 + clearance) + joint spheres, subtracted from the
    plate in ONE cut.  ``Fuse`` (not ``Glue``) is REQUIRED: the single-stroke
    connectors cross other turns, so the tube self-overlaps, and ``Glue`` (which
    assumes non-overlapping shapes) silently returns the plate uncut; a
    per-primitive sequential fuse (``body = body + s``) is O(N^2) and unusable
    (>2 min at ~800 primitives).  ``decimate`` (m) resamples the wire spine to
    bound the primitive count, and ``max_segments`` hard-caps the Fuse cost
    (coarsening the step further if needed -- reported, not silent).  Returns
    geometry provenance; raises (No-Fallback) if the cut yields no solid or
    removes no material (the channel missed the plate)."""
    import netgen.occ as occ
    chain0 = np.asarray(chain, float)
    chain = _decimate_polyline(chain0, decimate)
    decimate_used = float(decimate)
    if len(chain) - 1 > max_segments:                  # cap the Fuse cost
        seglen = np.linalg.norm(np.diff(chain0, axis=0), axis=1)
        decimate_used = float(seglen.sum()) / max_segments
        chain = _decimate_polyline(chain0, decimate_used)
    if len(chain) < 2:
        raise ValueError("former: wire chain has < 2 points")
    r = 0.5 * float(wire_diam) + float(clearance)
    tools = [occ.Sphere(occ.Pnt(*map(float, chain[0])), r)]
    for k in range(len(chain) - 1):
        p, q = chain[k], chain[k + 1]
        d = q - p
        L = float(np.linalg.norm(d))
        if L < 1e-12:
            continue
        tools.append(occ.Cylinder(occ.Pnt(*map(float, p)),
                                   occ.Dir(*map(float, d / L)), r=r, h=L))
        tools.append(occ.Sphere(occ.Pnt(*map(float, q)), r))
    mn = chain.min(axis=0) - np.array([margin, margin, r + wall])
    mx = chain.max(axis=0) + np.array([margin, margin, r + wall])
    plate = occ.Box(occ.Pnt(*map(float, mn)), occ.Pnt(*map(float, mx)))
    plate_vol = float(plate.mass)
    former = plate - occ.Fuse(tools)
    vol = float(former.mass)
    if not (vol > 0.0):
        raise ValueError("former: plate - channel yielded no solid (check "
                         "--wire-diam / --former-wall / --former-margin)")
    if vol >= plate_vol * (1.0 - 1e-9):
        raise ValueError("former: the channel removed no material (the wire "
                         "tube did not intersect the plate)")
    former.WriteStep(filename)
    return {
        "former_step": filename,
        "channel_segments": int(len(chain) - 1),
        "channel_radius_m": r,
        "decimate_step_m": decimate_used,
        # STEP coords are in METRES (Radia convention) -- scale x1000 in the
        # slicer if it interprets STEP as mm.
        "former_size_m": [float(mx[i] - mn[i]) for i in range(3)],
        "former_volume_m3": vol,
        "plate_volume_m3": plate_vol,
        "channel_volume_m3": plate_vol - vol,
    }


def _bubble_seeds(pts2d, bmag, k_const, n_max=240):
    """Bubble-system seed placement (Hirahatake/Noguchi/Igarashi/Yamashita):
    greedily drop seeds, biggest |B| first, keeping a clearance radius
    ``r = k_const / sqrt(|B|)`` -- so the resulting line DENSITY is
    proportional to |B| (the bubble is small where the field is strong, so
    more lines pack in).  Returns the seed coordinates (M, 2)."""
    order = np.argsort(-bmag)
    seeds = []
    for idx in order:
        b = bmag[idx]
        if b <= 0.0:
            continue
        p = pts2d[idx]
        r = k_const / np.sqrt(b)
        if all(np.hypot(p[0] - s[0], p[1] - s[1]) > r for s in seeds):
            seeds.append((p[0], p[1]))
            if len(seeds) >= n_max:
                break
    return np.array(seeds) if seeds else np.zeros((0, 2))


def _flux_line_plot(chain, current, out_path, plane="y", half=None,
                    ngrid=64, n_lines=46):
    """Bubble-system flux-line visualisation of the DESIGNED coil's field on a
    cut-plane (default y=0, the x-z plane).  Traces the in-plane B direction
    (matplotlib streamplot) seeded by the bubble system so the line density
    reflects |B| (definition (i): density ~ |B|; (ii): tangent ~ B).  The same
    'equal-flux-per-tube' principle as the stream-function contour wires."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    chain = np.asarray(chain, float)
    if half is None:
        half = 1.2 * float(np.max(np.abs(chain)))
    axn = {"x": 0, "y": 1, "z": 2}[plane]
    ia, ib = [c for c in (0, 1, 2) if c != axn]          # in-plane axes
    g = np.linspace(-half, half, ngrid)
    A, B = np.meshgrid(g, g)                              # in-plane grid
    obs = np.zeros((A.size, 3))
    obs[:, ia] = A.ravel(); obs[:, ib] = B.ravel()
    Bvec = _segment_field_B(chain, obs, current)
    Ba = Bvec[:, ia].reshape(ngrid, ngrid)
    Bb = Bvec[:, ib].reshape(ngrid, ngrid)
    bmag = np.hypot(Ba, Bb)
    # bubble seeds: tune k so ~n_lines seeds land
    bref = float(np.median(bmag[bmag > 0])) + 1e-30
    k = 2.0 * half / np.sqrt(n_lines) * np.sqrt(bref)
    seeds = _bubble_seeds(np.column_stack([A.ravel(), B.ravel()]),
                          bmag.ravel(), k)
    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    ax.imshow(bmag, extent=[-half, half, -half, half], origin="lower",
              cmap="magma", alpha=0.35, aspect="equal")
    lw = 0.6 + 2.4 * (bmag / (bmag.max() + 1e-30))
    try:
        ax.streamplot(g, g, Ba, Bb, color="white", linewidth=lw,
                      density=2.5, arrowsize=0.6,
                      start_points=seeds if len(seeds) else None)
    except Exception:                              # start_points off-grid edge
        ax.streamplot(g, g, Ba, Bb, color="white", linewidth=lw, density=2.0)
    lbl = "xyz"
    # lab rule 9 (units in PARENTHESES, not brackets) + rule 3 (italic math
    # symbol, upright unit): "x (m)" with x italic, not "x [m]".
    ax.set_xlabel(f"${lbl[ia]}$ (m)"); ax.set_ylabel(f"${lbl[ib]}$ (m)")
    # Lab figure convention: NO in-figure title -- the cut-plane and line count
    # are returned (the panel / docs caption carries them, not the PNG).
    plt.rcParams["pdf.fonttype"] = 42                  # TrueType (no Type-3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return {"flux_plot": out_path, "plane": f"{lbl[axn]}=0",
            "n_flux_lines": int(len(seeds))}


def _scan_map_3component(chain, current, axn, plane_point, half, n):
    """3-component (Bx,By,Bz) field of the single-stroke wire on a planar raster
    scan -- the SIMULATION TWIN of a 3-axis magnetometer (TMR / fluxgate)
    scanning a plane a standoff above the coil.

    ``axn`` is the plane-normal axis (0/1/2); the field is sampled on the two
    in-plane axes over a square (u, v) grid centred at ``plane_point``, half-
    width ``half`` (m), ``n`` points per axis, at the physical drive ``current``
    (A) (the field is linear in current, so this scales the whole map to the
    bench current the sensor sees).  Returns the grid, the three Cartesian
    components, and a homogeneity / transverse-leakage summary.

    ``Bnormal`` is the component along the plane normal (the wanted uniform
    field for a normal-aligned target); ``Bt`` is the in-plane transverse
    magnitude a 3-axis sensor captures (single-stroke connector stray) that a
    single-axis probe would miss."""
    ia, ib = [c for c in (0, 1, 2) if c != axn]
    u = np.linspace(-half, half, int(n))
    v = np.linspace(-half, half, int(n))
    U, V = np.meshgrid(u, v, indexing="ij")            # U[i,j]=u[i], V[i,j]=v[j]
    obs = np.empty((U.size, 3))
    obs[:, axn] = float(plane_point[axn])
    obs[:, ia] = float(plane_point[ia]) + U.ravel()
    obs[:, ib] = float(plane_point[ib]) + V.ravel()
    B = _segment_field_B(np.asarray(chain, float), obs, float(current))
    nn = int(n)
    Bx = B[:, 0].reshape(nn, nn); By = B[:, 1].reshape(nn, nn)
    Bz = B[:, 2].reshape(nn, nn)
    Bn = B[:, axn].reshape(nn, nn)                      # along-normal (main)
    Bt = np.hypot(B[:, ia], B[:, ib]).reshape(nn, nn)   # in-plane transverse
    Bmag = np.linalg.norm(B, axis=1).reshape(nn, nn)
    bn_mean = float(Bn.mean())
    den = abs(bn_mean) + 1e-30
    lbl = "xyz"
    return {
        "plane_normal": lbl[axn], "axis_u": lbl[ia], "axis_v": lbl[ib],
        "plane_point_m": [float(p) for p in plane_point],
        "half_m": float(half), "n": nn, "current_A": float(current),
        # full grid + components (durable measurement-comparison artifact)
        "grid_u_m": [float(x) for x in u],
        "grid_v_m": [float(x) for x in v],
        "Bx_T": Bx.tolist(), "By_T": By.tolist(), "Bz_T": Bz.tolist(),
        # summary (rendered in the panel Output; full arrays are in scan_output)
        "Bnormal_mean_T": bn_mean,
        "Bnormal_min_T": float(Bn.min()), "Bnormal_max_T": float(Bn.max()),
        "Bnormal_homogeneity_rms": float(np.sqrt(np.mean((Bn - bn_mean) ** 2))
                                         / den),
        "Bnormal_ptp_rel": float((Bn.max() - Bn.min()) / den),
        "Bt_mean_T": float(Bt.mean()), "Bt_max_T": float(Bt.max()),
        "Bt_over_Bn_max": float(Bt.max() / den),
        "Bmag_mean_T": float(Bmag.mean()),
        "Bmag_min_T": float(Bmag.min()), "Bmag_max_T": float(Bmag.max()),
    }


def _write_scan_map_json(path, args, scan):
    """Persist the full 3-component scan map + provenance to committed JSON
    (Data Persistence Policy: figure/result-backing data always saved next to
    its driver, never left transient).  This is the file a measured 3-axis
    scan (CSV, current-reversal averaged) is compared against."""
    import json
    import platform
    from datetime import datetime, timezone
    try:
        import radia as _r
        rv = getattr(_r, "__version__", None)
    except Exception:                                  # noqa: BLE001
        rv = None
    doc = {
        "artifact": "sf_scan_map_3component",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "radia_version": rv,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "coil_vol": os.path.basename(args.coil_vol),
        "eval_vol": os.path.basename(args.eval_vol),
        "target_cf": args.target_cf,
        "distort": bool(getattr(args, "distort", False)),
    }
    doc.update(scan)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(doc, f, indent=2, ensure_ascii=False)


def _equal_3d(ax, pts):
    """Equal aspect box for a 3D axis around a point cloud."""
    pts = np.asarray(pts, float)
    c = pts.mean(axis=0)
    r = float(np.max(np.abs(pts - c))) + 1e-9
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:                              # older matplotlib
        pass
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])


def _tube_mesh(path, radius, n_circ=12):
    """Sweep a circle of given radius along the polyline with a twist-free
    (parallel-transport) frame -> (X, Y, Z) for plot_surface.  This is the
    'with thickness' rendering of the single conductor."""
    P = np.asarray(path, float)
    if len(P) < 2:
        return None
    T = np.gradient(P, axis=0)
    T /= (np.linalg.norm(T, axis=1, keepdims=True) + 1e-30)
    a = np.array([0.0, 0.0, 1.0]) if abs(T[0, 2]) < 0.9 else np.array([1.0, 0, 0])
    Nprev = a - (a @ T[0]) * T[0]
    Nprev /= (np.linalg.norm(Nprev) + 1e-30)
    Ns = np.zeros_like(P)
    Ns[0] = Nprev
    for i in range(1, len(P)):
        v = np.cross(T[i - 1], T[i])
        s = float(np.linalg.norm(v))
        Ni = Nprev
        if s > 1e-9:                               # rotate the frame
            ax_ = v / s
            ang = np.arctan2(s, float(T[i - 1] @ T[i]))
            Ni = (Nprev * np.cos(ang) + np.cross(ax_, Nprev) * np.sin(ang)
                  + ax_ * (ax_ @ Nprev) * (1 - np.cos(ang)))
        Ni = Ni - (Ni @ T[i]) * T[i]
        Ni /= (np.linalg.norm(Ni) + 1e-30)
        Ns[i] = Ni
        Nprev = Ni
    Bb = np.cross(T, Ns)
    th = np.linspace(0.0, 2.0 * np.pi, n_circ)
    ct, st = np.cos(th), np.sin(th)
    X = np.empty((len(P), n_circ)); Y = np.empty_like(X); Z = np.empty_like(X)
    for i in range(len(P)):
        ring = (P[i][None, :] + radius * (ct[:, None] * Ns[i][None, :]
                                          + st[:, None] * Bb[i][None, :]))
        X[i], Y[i], Z[i] = ring[:, 0], ring[:, 1], ring[:, 2]
    return X, Y, Z


def _steps_plot(loops, chain, dist_chain, diam, out_path, target_cf=""):
    """Per-step manufacturing visualisation -- a 2x2 3D figure of the coil at
    each stage: (a) equal-current iso-contours (N = nlevels turns), (b) the
    single-stroke wire, (c) the sheet-metal (bankin-ho) distorted wire, (d) the
    wire WITH conductor thickness + distortion.

    Lab figure convention (applied inline -- the shipped panel calc cannot
    import radia_mcp.figure): NO in-figure title / suptitle; panels carry small
    (a)-(d) CORNER LABELS (not titles, so the _check_no_in_figure_title gate
    passes), TrueType fonts, English text; the full stage descriptions are
    returned in ``stages`` for the panel / docs caption."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d)
    plt.rcParams["pdf.fonttype"] = 42                  # TrueType (no Type-3)
    final = dist_chain if dist_chain is not None else chain
    distorted = dist_chain is not None
    allpts = np.vstack([chain, final])
    fig = plt.figure(figsize=(11.0, 9.0))

    def _tag(ax, s):
        ax.text2D(0.02, 0.96, s, transform=ax.transAxes, fontsize=10)

    ax = fig.add_subplot(2, 2, 1, projection="3d")
    for p in loops:
        ax.plot(p[:, 0], p[:, 1], p[:, 2], lw=0.7)
    _tag(ax, f"(a) contours (N={len(loops)})")
    _equal_3d(ax, allpts)

    ax = fig.add_subplot(2, 2, 2, projection="3d")
    ax.plot(chain[:, 0], chain[:, 1], chain[:, 2], lw=0.7, color="C0")
    _tag(ax, "(b) single-stroke wire")
    _equal_3d(ax, allpts)

    ax = fig.add_subplot(2, 2, 3, projection="3d")
    ax.plot(final[:, 0], final[:, 1], final[:, 2], lw=0.7, color="C3")
    _tag(ax, "(c) sheet-metal distorted" + ("" if distorted else " (= b)"))
    _equal_3d(ax, allpts)

    ax = fig.add_subplot(2, 2, 4, projection="3d")
    tube = _tube_mesh(final, max(diam, 1e-6) / 2.0)
    if tube is not None:
        ax.plot_surface(tube[0], tube[1], tube[2], color="C1",
                        alpha=0.9, linewidth=0, antialiased=True)
    _tag(ax, f"(d) thickness d={diam*1e3:.1f} mm")
    _equal_3d(ax, allpts)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return {"steps_plot": out_path, "target_cf": target_cf, "stages": [
        f"(a) equal-current iso-contours, N={len(loops)} turns",
        "(b) single-stroke (ikkaki) wire",
        "(c) sheet-metal (bankin-ho) distorted wire"
        + ("" if distorted else " -- none applied, equals (b)"),
        f"(d) with conductor thickness d={diam*1e3:.1f} mm + distortion",
    ]}


def _peec_inductance(chain, diam, sigma, freq):
    """PEEC L, R of the single-stroke wire (one port across the two ends)."""
    from radia.peec_matrices import PEECBuilder
    from radia.peec_topology import PEECCircuitSolver
    chain = _dedupe(chain, tol=1e-9)
    if len(chain) < 2:
        raise ValueError("chain too short for PEEC")
    b = PEECBuilder()
    ids = [b.add_node_at(float(p[0]), float(p[1]), float(p[2])) for p in chain]
    for k in range(len(chain) - 1):
        b.add_connected_segment(ids[k], ids[k + 1], diam, diam, sigma)
    b.add_port(ids[0], ids[-1])
    solver = PEECCircuitSolver(b.build_topology())
    Z = complex(solver.compute_port_impedance(freq))
    omega = 2.0 * np.pi * freq
    return {"freq_Hz": freq, "Z_re_ohm": Z.real, "Z_im_ohm": Z.imag,
            "L_H": (Z.imag / omega) if omega > 0 else 0.0, "R_ohm": Z.real,
            "n_nodes": int(len(chain))}


def _coil_L_at_nlevels(P, psi_v, verts, tris, tol, args, nlevels):
    """Coil inductance at a given turn count: vertex-linear contours -> nearest-
    neighbour single stroke -> PEEC L_coil.  The inductance-target search uses
    nn for speed (L_coil is dominated by the turns, not the connector cuts); the
    final manufacture re-chains with args.chain and reports the accurate L."""
    levels = _contour_levels(psi_v, nlevels)
    loops = []
    for lev in levels:
        loops.extend(_segs_to_polylines(
            _contour_segments(verts, psi_v, tris, lev), tol))
    if not loops:
        return None
    kdt, tri_K = _build_gradpsi_kdt(verts, psi_v, tris)
    fields, signs = _gradpsi_signs(loops, kdt, tri_K, P["mpts"], P["vector_b"],
                                   P["Bm"])
    loops_signed = [p if s >= 0.0 else p[::-1] for p, s in zip(loops, signs)]
    chain, _nc, _cl = _single_stroke_chain(loops_signed)
    chain = _dedupe(chain)
    if len(chain) < 2:
        return None
    pe = _peec_inductance(chain, args.wire_diam, 5.8e7, args.peec_freq)
    return float(pe["L_H"]), int(len(loops))


def _search_nlevels_for_L(P, psi_v, verts, tris, tol, args, L_target,
                          nl_min=4, nl_max=60):
    """Find the turn count (nlevels) whose single-stroke coil inductance is
    closest to ``L_target`` -- the IH resonance condition L_target = 1/(omega^2
    C).  L_coil is monotone increasing in nlevels (~N^1.5..2) so bisect on the
    integer turn count (memoised PEEC evals).  Returns the search result with
    the achievable [L_min, L_max] range and an ``in_range`` flag: an
    out-of-range target is REPORTED (No-Fallback), so the user fixes the
    geometry / tank capacitor rather than getting a silently-clamped coil."""
    trace, memo = [], {}

    def L_of(n):
        n = int(n)
        if n in memo:
            return memo[n]
        r = _coil_L_at_nlevels(P, psi_v, verts, tris, tol, args, n)
        L = r[0] if r else float("nan")
        memo[n] = L
        trace.append({"nlevels": n, "L_H": L})
        return L

    L_lo, L_hi = L_of(nl_min), L_of(nl_max)
    if np.isnan(L_lo) or np.isnan(L_hi):       # No-Fallback: nan poisons the
        raise ValueError(                       # bisection -> raise, don't drift
            f"cannot compute coil inductance at a search endpoint "
            f"(nlevels={nl_min} -> {L_lo}, nlevels={nl_max} -> {L_hi}): the "
            f"single-stroke chain is empty / too short at that turn count. "
            f"Check the design / contours or adjust --nlevels-max.")
    in_range = (min(L_lo, L_hi) <= L_target <= max(L_lo, L_hi))
    lo, hi = nl_min, nl_max
    best_n, best_L = nl_min, L_lo
    best_err = abs(L_lo - L_target)
    while hi - lo > 1:                              # bisection (L monotone up)
        mid = (lo + hi) // 2
        Lm = L_of(mid)
        if abs(Lm - L_target) < best_err:
            best_n, best_L, best_err = mid, Lm, abs(Lm - L_target)
        if Lm < L_target:
            lo = mid
        else:
            hi = mid
    for n in (lo, hi):
        Ln = L_of(n)
        if abs(Ln - L_target) < best_err:
            best_n, best_L, best_err = n, Ln, abs(Ln - L_target)
    return {"nlevels": int(best_n), "L_coil_search_H": float(best_L),
            "L_target_H": float(L_target),
            "rel_error_search": float(best_err / (abs(L_target) + 1e-30)),
            "L_range_H": [float(min(L_lo, L_hi)), float(max(L_lo, L_hi))],
            "nlevels_range": [int(nl_min), int(nl_max)],
            "in_range": bool(in_range), "n_peec_evals": len(memo)}


def run_design(args):
    from ngsolve import GridFunction, TaskManager
    t0 = time.perf_counter()
    with TaskManager():
        if getattr(args, "shield_vol", None) and getattr(args, "iron_vol", None):
            raise ValueError("--iron-vol (material-aware kernel) and "
                             "--shield-vol (active-shield coil) are not yet "
                             "combinable; choose one.")
        if getattr(args, "shield_vol", None):
            P = _build_shielded_problem(args)
        else:
            P = _build_problem(args)
        t1 = time.perf_counter()
        psi_f, fit_rms, homo = _solve_and_metrics(P, args.alpha)

        if P.get("is_shielded"):
            # Split combined psi_f into primary ([:n_fp]) and shield ([n_fp:])
            n_fp = P["n_primary"]
            psi_fp, psi_fs = psi_f[:n_fp], psi_f[n_fp:]
            psi_p = P["R_p"] @ psi_fp
            psi_s = P["R_s"] @ psi_fs
            gfu_p = GridFunction(P["fes"])
            gfu_s = GridFunction(P["shield_fes"])
            gfu_p.vec.FV().NumPy()[:] = psi_p
            gfu_s.vec.FV().NumPy()[:] = psi_s
            peak_J_p = _peak_current_density(P["fes"], P["coil"], gfu_p)
            peak_J_s = _peak_current_density(P["shield_fes"],
                                             P["shield_coil"], gfu_s)
            peak_J = max(peak_J_p, peak_J_s)
            # DSV-only fit (top Ac_f rows, not stacked with external)
            dsv_rms = float(np.linalg.norm(P["Ac_f"] @ psi_f - P["Bc_dsv"])
                            / (np.linalg.norm(P["Bc_dsv"]) + 1e-30))
            B_dsv_norm = np.linalg.norm(P["Bc_dsv"]) + 1e-30
            # HONEST stray = combined field at INDEPENDENT external measure
            # points (e_mpts, NOT the constraint vertices) -- relative to the
            # DSV target field.  This is the generalising metric (the
            # constraint-point residual is ~machine-zero and circular).
            stray_rms = float(np.linalg.norm(P["Aem_f"] @ psi_f) / B_dsv_norm)
            # the (circular) fit residual at the null-constraint points, kept
            # for transparency only -- NOT the headline shielding number
            stray_fit_rms = float(np.linalg.norm(P["Ae_f"] @ psi_f) / B_dsv_norm)
            gfu = gfu_p          # primary psi for harmonic decomp / msh
        else:
            psi = P["R"] @ psi_f
            gfu = GridFunction(P["fes"])
            gfu.vec.FV().NumPy()[:] = psi
            peak_J = _peak_current_density(P["fes"], P["coil"], gfu)
            dsv_rms = fit_rms
            stray_rms = None
            peak_J_p = peak_J_s = None

        t2 = time.perf_counter()
        result = {
            "method": "design",
            "target_cf": args.target_cf,
            "target_kind": "vector_B" if P["vector_b"] else "Bz",
            "coil_vol": os.path.basename(args.coil_vol),
            "eval_vol": os.path.basename(args.eval_vol),
            "order": args.order, "regularize": args.regularize,
            "alpha": args.alpha, "confine": P["confine"],
            "ndof": int(P["fes"].ndof), "ndof_free": int(P["n_free"]),
            "n_constraint": int(P["n_constraint"]),
            "n_measure": int(P["n_measure"]),
            "n_constraints": int(P["Ac_f"].shape[0]
                                 if P.get("is_shielded") else P["Ac"].shape[0]),
            "fit_residual_rms": fit_rms,
            "homogeneity_rms": homo,
            "rms": homo,                       # primary metric = true homogeneity
            "peak_J": peak_J,
            "t_mesh_s": round(t1 - t0, 3),
            "t_solve_s": round(t2 - t1, 3),
            "t_total_s": round(time.perf_counter() - t0, 3),
        }
        # Material-aware (iron) metadata
        if getattr(args, "iron_vol", None):
            result.update({
                "material_aware": True,
                "iron_vol": os.path.basename(args.iron_vol),
                "mu_r": float(args.mu_r),
                "iron_mat": getattr(args, "iron_mat", "iron"),
                "iron_exact_source": bool(getattr(args, "iron_exact_source",
                                                  False)),
            })
        # Active-shielding extra metrics
        if P.get("is_shielded"):
            result.update({
                "shielded": True,
                "shield_vol": os.path.basename(args.shield_vol),
                "shield_eval_vol": os.path.basename(
                    args.shield_eval_vol),
                "shield_weight": P["shield_weight"],
                "ndof_primary": int(P["fes"].ndof),
                "ndof_shield": int(P["shield_fes"].ndof),
                "ndof_free_primary": int(P["n_primary"]),
                "ndof_free_shield": int(P["n_shield"]),
                "n_external": int(P["n_external"]),
                "n_external_measure": int(P["n_external_measure"]),
                "dsv_rms": dsv_rms,
                # HONEST stray = field at INDEPENDENT external measure points,
                # relative to the DSV target field (a real generalisation
                # metric, NOT the circular constraint-point residual below).
                "stray_rms": stray_rms,
                "stray_fit_rms": stray_fit_rms,   # circular fit residual (info)
                "peak_J_primary": peak_J_p,
                "peak_J_shield": peak_J_s,
            })
        # spherical-harmonic decomposition of the ACHIEVED Bz over the DSV --
        # uses Am @ psi_f (Am is the DSV measure matrix for both single/shielded)
        if args.harmonic_lmax > 0 and not P["vector_b"]:
            bz = P["Am"] @ psi_f
            result["harmonics"] = _harmonic_decompose(
                P["mpts"], bz, args.harmonic_lmax,
                target_terms=getattr(args, "_harmonic_terms", None) or None)
        if args.regularize == "inductance":
            Sm = P["S"]
            result["stored_energy_J"] = 0.5 * float(psi_f @ (Sm @ psi_f))
        if args.msh_output:
            try:
                from radia.gmsh_post_export import GmshPostExport
                post = GmshPostExport(P["coil"], boundary=True)
                # primary psi at vertex nodes
                psi_v = (psi_p if P.get("is_shielded")
                         else psi)[:P["coil"].nv]
                post.add_field("psi", psi_v, ncomp=1)
                post.write(args.msh_output)
                result["msh"] = args.msh_output
            except Exception as e:                 # noqa: BLE001
                result["msh_error"] = str(e)
    return result


def _pareto_alpha(P, args):
    """alpha lever: Tikhonov L-curve (re-uses one factorisation per point)."""
    from ngsolve import GridFunction
    lams = np.concatenate([[0.0], np.logspace(
        np.log10(args.alpha_min), np.log10(args.alpha_max),
        max(1, args.n_alpha - 1))])
    front = []
    for lam in lams:
        psi_f, fit_rms, homo = _solve_and_metrics(P, lam)
        psi = P["R"] @ psi_f
        gfu = GridFunction(P["fes"]); gfu.vec.FV().NumPy()[:] = psi
        peak = _peak_current_density(P["fes"], P["coil"], gfu)
        front.append({"alpha": float(lam), "homogeneity_rms": homo,
                      "peak_J": peak})
    return front


def _pareto_geometry(args):
    """geometry lever: sweep the eval-region (DSV) scale about its centroid --
    the dual of the former-size lever; each scale -> one (homogeneity, peak)
    point at the design alpha.  Re-builds A per scale (the coil moves relative
    to the target)."""
    from ngsolve import GridFunction
    scales = np.linspace(args.geom_scale_min, args.geom_scale_max,
                         max(2, args.n_alpha))
    front = []
    for s in scales:
        P = _build_problem(args, eval_scale=float(s))
        psi_f, fit_rms, homo = _solve_and_metrics(P, args.alpha)
        psi = P["R"] @ psi_f
        gfu = GridFunction(P["fes"]); gfu.vec.FV().NumPy()[:] = psi
        peak = _peak_current_density(P["fes"], P["coil"], gfu)
        front.append({"eval_scale": float(s), "homogeneity_rms": homo,
                      "peak_J": peak})
    return front


def run_pareto(args):
    """(homogeneity, peak current density) Pareto front via a chosen lever:
    alpha (Tikhonov L-curve), linf (minimax IRLS trajectory) or geometry
    (eval-region scale sweep)."""
    from ngsolve import TaskManager
    t0 = time.perf_counter()
    with TaskManager():
        if args.pareto_lever == "geometry":
            P = _build_problem(args)        # for metadata (ndof etc.)
            front = _pareto_geometry(args)
        else:
            P = _build_problem(args)
            if args.pareto_lever == "alpha":
                front = _pareto_alpha(P, args)
            elif args.pareto_lever == "linf":
                front = _linf_irls_front(P, args.alpha, args.linf_iter)
            else:                            # pragma: no cover - argparse guards
                raise ValueError(f"unknown pareto-lever {args.pareto_lever}")
    return {
        "method": "pareto", "pareto_lever": args.pareto_lever,
        "target_cf": args.target_cf,
        "target_kind": "vector_B" if P["vector_b"] else "Bz",
        "coil_vol": os.path.basename(args.coil_vol),
        "eval_vol": os.path.basename(args.eval_vol),
        "order": args.order, "regularize": args.regularize,
        "alpha": args.alpha,
        "ndof": int(P["fes"].ndof), "n_measure": int(P["n_measure"]),
        "n_alpha": len(front), "n_points": len(front), "front": front,
        "t_total_s": round(time.perf_counter() - t0, 3),
    }


def run_manufacture(args):
    """psi iso-contours -> single-stroke wire -> (sheet-metal distort) ->
    CAD(STEP, --step-output) -> PEEC inductance (--peec).  Reports the
    manufactured-wire field homogeneity at the eval region (closes the design
    loop: does the discretised wire still hit the target?)."""
    from ngsolve import GridFunction, TaskManager
    # No-Fallback: a manufacture knob must FAIL LOUD on a value that would
    # silently degrade or CRASH, not produce a quiet wrong wire.  chain-ncut<0 ->
    # step=1 max resolution; chain-passes<0 / distort-iter<1 -> the optimiser
    # range() is EMPTY so it is SKIPPED and the wire silently un-optimised;
    # distort-grid<2 -> a degenerate control grid that crashes with a cryptic
    # numpy IndexError.
    _checks = [("chain_ncut", 0, "0 = auto-scale with loop count, positive = "
                "explicit cut candidates per loop"),
               ("chain_passes", 0, "0 = auto-scale, positive = explicit descent "
                "rounds")]
    if getattr(args, "distort", False):
        _checks += [("distort_grid", 2, "control-grid points per axis"),
                    ("distort_iter", 1, "Gauss-Newton iterations")]
    if getattr(args, "scan_map", False):
        _checks += [("scan_n", 2, "raster points per axis (n x n grid)")]
    for _nm, _lo, _hint in _checks:
        _v = int(getattr(args, _nm, _lo))
        if _v < _lo:
            return {"method": "manufacture",
                    "error": f"--{_nm.replace('_', '-')} must be >= {_lo} "
                             f"({_hint}); got {_v}"}
    t0 = time.perf_counter()
    with TaskManager():
        P = _build_problem(args)
        t1 = time.perf_counter()
        psi_f, fit_rms, homo = _solve_and_metrics(P, args.alpha)
        psi = P["R"] @ psi_f
        gfu = GridFunction(P["fes"]); gfu.vec.FV().NumPy()[:] = psi
        peak_J = _peak_current_density(P["fes"], P["coil"], gfu)
        t_psi = time.perf_counter()

        coil = P["coil"]
        verts = np.array([list(v.point) for v in coil.vertices])
        psi_v = psi[:coil.nv]
        tris = _surface_triangles(coil)
        clen = _char_len(verts, tris)
        tol = max(1e-12, 1e-6 * clen)

        # IH-resonance design: pick the turn count (nlevels) whose single-stroke
        # coil inductance resonates -- L_target = 1/(omega^2 C).  Unlike the
        # gradient-coil min-inductance objective, IH wants a SPECIFIC L (so the
        # work coil + tank cap resonate at the inverter frequency).  L_coil is
        # set by the turns (~N^2), so we search nlevels; the design (field) is
        # nlevels-independent, so this re-uses the one psi solve.
        resonance = None
        tgt_L = getattr(args, "target_inductance", None)
        cap = getattr(args, "resonance_cap", None)
        greedy_mode = bool(getattr(args, "greedy_turns", None))
        if tgt_L or cap:
            if cap:                                  # L from the tank cap + freq
                w = 2.0 * np.pi * args.peec_freq
                tgt_L = 1.0 / (w * w * cap)
            if greedy_mode:
                # GREEDY owns the turn count (the low-turn coil), so we CANNOT
                # also search nlevels to hit L_target -- the two levers both set
                # the turns and would fight (L_coil ~ N^2 pins N).  No-Fallback:
                # instead of silently dropping the resonance spec, report what
                # the delivered few-turn coil ACHIEVES + the tank capacitor it
                # needs to resonate at the operating frequency (filled from the
                # final PEEC L below).  Turns stay = greedy_turns; nlevels is
                # left untouched.
                resonance = {"mode": "from_greedy_coil",
                             "L_target_H": float(tgt_L),
                             "operating_freq_Hz": float(args.peec_freq),
                             "turns_source": "greedy"}
                if cap is not None:
                    resonance["resonance_cap_F"] = float(cap)
                args.peec = True                     # need the accurate final L
            else:
                resonance = _search_nlevels_for_L(
                    P, psi_v, verts, tris, tol, args, tgt_L,
                    nl_max=int(getattr(args, "nlevels_max", 60)))
                resonance["mode"] = "search_nlevels"
                args.nlevels = resonance["nlevels"]  # design AT the resonant N
                args.peec = True                     # need the accurate final L

        greedy_trace = None
        greedy_dict_used = None
        pin_tiling_used = False
        if getattr(args, "greedy_turns", None):
            # GREEDY CONSTRUCTIVE coil (low-turn strategy): lay turns one at a
            # time where each most reduces the residual (monotone by build) --
            # the antidote to the non-monotonic contour-quantise turn count.
            greedy_dict_used = getattr(args, "greedy_dict", "contour")
            if greedy_dict_used == "pin":
                cand_loops, cand_fields = _pin_loop_candidates(
                    verts, tris, P["mpts"], P["vector_b"])
            elif greedy_dict_used == "bubble":
                cand_loops, cand_fields = _bubble_pin_loops(
                    psi_v, verts, tris, P["mpts"], P["vector_b"],
                    n_target=int(getattr(args, "pin_tiling_pins", 80)),
                    frac=float(getattr(args, "pin_tiling_frac", 0.7)))
            else:
                cand_loops, cand_fields = _contour_loop_candidates(
                    psi_v, verts, tris, P["mpts"], P["vector_b"], tol)
            cand_cents = (np.array([p.mean(axis=0) for p in cand_loops])
                          if cand_loops else None)
            coil_diag = float(np.linalg.norm(verts.max(0) - verts.min(0)))
            sel, greedy_trace = _greedy_coil_turns(
                cand_fields, P["Bm"], int(args.greedy_turns),
                getattr(args, "greedy_target", None),
                cand_centroids=cand_cents,
                connector_weight=float(getattr(args, "greedy_connector_weight",
                                               0.0)),
                coil_diag=coil_diag)
            loops = [cand_loops[k] if s >= 0.0 else cand_loops[k][::-1]
                     for (k, s) in sel]
            levels = []; sub = 1                     # greedy: no contour levels
            lvl_rms0 = lvl_rms = None
            if not loops:
                return {"method": "manufacture",
                        "error": "greedy coil selected no turns (empty candidate "
                                 "dictionary, or no loop improves the target)"}
        elif getattr(args, "pin_tiling", False):
            # DENSE BUBBLE-TILING pin coil (many-turn, manufacturable): equal-
            # current loops packed with density ~ |K| = |grad psi| and sized to
            # their bubble -- the discrete-loop alternative to long iso-contours.
            pin_tiling_used = True
            loops, _ = _bubble_pin_loops(
                psi_v, verts, tris, P["mpts"], P["vector_b"],
                n_target=int(getattr(args, "pin_tiling_pins", 80)),
                frac=float(getattr(args, "pin_tiling_frac", 0.7)))
            levels = []; sub = 1
            lvl_rms0 = lvl_rms = None
            if not loops:
                return {"method": "manufacture",
                        "error": "pin-tiling produced no loops (psi is constant, "
                                 "or the surface mesh is too coarse for the k-ring "
                                 "bubble loops)"}
        else:
            if getattr(args, "optimize_levels", False):
                levels, lvl_rms0, lvl_rms = _optimize_contour_levels(
                    psi_v, verts, tris, args.nlevels, P["mpts"], P["Bm"],
                    P["vector_b"], tol)
            else:
                levels = _contour_levels(psi_v, args.nlevels)
                lvl_rms0 = lvl_rms = None
            sub = max(1, int(getattr(args, "contour_sub", 1)))
            if sub > 1:                              # order-p curved contour
                hp = _contour_segments_hp(coil, gfu, levels, sub)
                loops = []
                for lev in levels:
                    loops.extend(_segs_to_polylines(hp[lev], tol))
            else:                                    # vertex-linear (fast default)
                loops = []
                for lev in levels:
                    loops.extend(_segs_to_polylines(
                        _contour_segments(verts, psi_v, tris, lev), tol))
            if not loops:
                return {"method": "manufacture",
                        "error": "psi has no iso-contours at the requested levels "
                                 "(check the target / coil surface)"}
        wire_len = sum(float(np.sum(np.linalg.norm(np.diff(p, axis=0), axis=1)))
                       for p in loops)
        vector_b = P["vector_b"]
        target_flat = P["Bm"]
        den_t = float(np.linalg.norm(target_flat)) + 1e-30

        # consistent contour orientation by the surface-current winding
        # K = n_hat x grad psi (NOT per-loop field alignment) -- correct for a
        # multi-region (saddle, l>=2) psi where sign(f.B) flattens the reversals
        # (reduces to sign(f.B) for a monotone l=1 psi).
        kdt, tri_K = _build_gradpsi_kdt(verts, psi_v, tris)
        fields, signs = _gradpsi_signs(loops, kdt, tri_K, P["mpts"], vector_b,
                                       target_flat)
        loops_signed = [p if s >= 0.0 else p[::-1]
                        for p, s in zip(loops, signs)]
        # separate-turns best (per-loop currents) -- the discretisation's best
        Fm = np.column_stack([s * f for f, s in zip(fields, signs)]) \
            if fields else np.zeros((len(target_flat), 0))
        if Fm.shape[1] > 0:
            I_loops, *_ = np.linalg.lstsq(Fm, target_flat, rcond=None)
            loops_homo = float(np.linalg.norm(Fm @ I_loops - target_flat)
                               / den_t)
        else:
            loops_homo = float("nan")

        # open contours (hit the coil boundary) -> closing them with a rim
        # chord injects a spurious edge current; a high count means the former
        # is too small for this target (extend it).  This is the dominant
        # single-current error when it happens -- NOT the connectors.
        n_open = sum(1 for p in loops
                     if np.linalg.norm(p[0] - p[-1]) > 1e-9)
        if args.chain == "nn":
            chain, n_conn, conn_len = _single_stroke_chain(loops_signed)
        else:                                  # field_aware (default)
            chain, n_conn, conn_len = _field_aware_chain(
                loops, signs, P["mpts"], target_flat, vector_b,
                n_cut=(int(args.chain_ncut) or None),
                passes=(int(args.chain_passes) or None))
        chain = _dedupe(chain)
        t2 = time.perf_counter()

        # single-stroke wire (ONE series current) field homogeneity
        f_unit = _wire_field_flat(chain, P["mpts"], vector_b, 1.0)
        I_w = _best_fit_current(f_unit, target_flat)
        wire_homo = float(np.linalg.norm(I_w * f_unit - target_flat) / den_t)

        out_chain = chain
        result = {
            "method": "manufacture",
            "target_cf": args.target_cf,
            "target_kind": "vector_B" if vector_b else "Bz",
            "coil_vol": os.path.basename(args.coil_vol),
            "eval_vol": os.path.basename(args.eval_vol),
            "order": args.order, "regularize": args.regularize,
            "alpha": args.alpha, "nlevels": args.nlevels,
            "confine": getattr(args, "confine", "off"),
            "chain_method": args.chain,
            "ndof": int(P["fes"].ndof), "n_measure": int(P["n_measure"]),
            "homogeneity_rms": homo,             # continuous psi (design ref)
            "peak_J": peak_J,
            "n_levels": len(levels), "n_loops": int(len(loops)),
            "contour_sub": sub,
            "n_open_contours": int(n_open),
            "n_wire_points": int(len(chain)),
            "wire_length_m": wire_len,
            "n_connectors": int(n_conn),
            "connector_length_m": conn_len,
            "best_fit_current_A": I_w,
            "loops_homogeneity_rms": loops_homo,  # separate turns (best disc.)
            "wire_homogeneity_rms": wire_homo,    # single-stroke (delivered)
            "rms": wire_homo,                     # primary = delivered wire
            "note": "homogeneity_rms = continuous psi design; "
                    "loops_homogeneity_rms = N separate contour turns (best "
                    "discretisation); wire_homogeneity_rms = single-stroke "
                    "one-current wire. With closed contours (n_open_contours=0) "
                    "the field-aware chain reaches the separate-turns floor "
                    "without --distort; a high n_open_contours means the former "
                    "is too small (extend it) -- that, not the connectors, is "
                    "the dominant single-current error.",
        }
        if greedy_trace is not None:
            result.update({
                "greedy_turns": True,
                "greedy_dict": greedy_dict_used,
                "greedy_connector_weight": float(
                    getattr(args, "greedy_connector_weight", 0.0)),
                "n_turns_used": len(greedy_trace),
                # rms vs n_turns -- MONOTONE by construction; read off the fewest
                # turns meeting any field spec (the correct "min turns" answer).
                "greedy_trace": greedy_trace,
                "greedy_final_rms": (greedy_trace[-1]["rms"] if greedy_trace
                                     else None),
            })
            if getattr(args, "greedy_plot", "") and greedy_trace:
                try:
                    _greedy_trace_plot(greedy_trace,
                                       getattr(args, "greedy_target", None),
                                       args.greedy_plot)
                    result["greedy_plot"] = args.greedy_plot
                except Exception as e:                 # noqa: BLE001
                    result["greedy_plot_error"] = str(e)
        elif pin_tiling_used:
            result.update({
                "pin_tiling": True,
                "pin_tiling_pins": int(getattr(args, "pin_tiling_pins", 80)),
                "pin_tiling_frac": float(getattr(args, "pin_tiling_frac", 0.7)),
                "n_pin_loops": int(len(loops)),
                # HONEST regime note: bubble pins reproduce the target under
                # INDEPENDENT currents (loops_homogeneity_rms) -- the natural
                # metric for a driven pin / shim ARRAY -- but NOT under one series
                # current (wire_homogeneity_rms), because isolated equal-current
                # loops do not tile-cancel into the surface current K the way psi
                # iso-contours do.  Use --pin-tiling for a driven-array layout;
                # use iso-contours (default) / --greedy-turns for a single wire.
                "pin_tiling_note": "loops_homogeneity_rms = driven pin/shim ARRAY "
                                   "(independent currents, the pin-coil regime); "
                                   "wire_homogeneity_rms = single series wire "
                                   "(use iso-contours for that).",
            })
        elif getattr(args, "optimize_levels", False):
            result.update({
                "levels_optimized": True,
                # equal-current (single series current) discrete-loop field RMS:
                # uniform equal-deltaI levels vs the optimised levels
                "equal_current_rms_uniform": lvl_rms0,
                "equal_current_rms_optimized": lvl_rms,
            })

        if args.distort:
            cw, Iw2, rms0, rmsb, disp = _distort_wire(
                chain, P["mpts"], target_flat, vector_b,
                n_grid=args.distort_grid, n_iter=args.distort_iter)
            out_chain = _dedupe(cw)
            result.update({
                "distort": True,
                "distort_grid": args.distort_grid,
                "distort_iter": args.distort_iter,
                "wire_homogeneity_rms_before": rms0,
                "wire_homogeneity_rms": rmsb,
                "rms": rmsb,
                "distort_max_disp_m": disp,
                "best_fit_current_A": Iw2,
            })
        t3 = time.perf_counter()

        if getattr(args, "scan_map", False):
            # 3-component (Bx,By,Bz) planar raster scan of the DELIVERED wire --
            # the simulation twin of a 3-axis magnetometer (TMR/fluxgate) scan.
            # Matches a fixed sensor ARRAY or a single sensor on a moving stage
            # (both sample the same grid); compare against a measured scan taken
            # with current reversal (+-I averaged -> coil-only, Earth + sensor
            # offset cancelled).
            axn = {"x": 0, "y": 1, "z": 2}[args.scan_plane]
            ia, ib = [c for c in (0, 1, 2) if c != axn]
            coil_c = verts.mean(axis=0)
            ev = np.asarray(P["mpts"], float)
            ev_c = ev.mean(axis=0) if len(ev) else coil_c
            # standoff along the normal from the coil centroid (default = the
            # eval-region / DSV height: scan where the design cares)
            if args.scan_standoff is not None:
                standoff = float(args.scan_standoff)
            else:
                standoff = float(ev_c[axn] - coil_c[axn])
            # in-plane half-width (default = eval-region in-plane half-extent)
            if args.scan_half is not None:
                half = float(args.scan_half)
            elif len(ev):
                half = 0.5 * float(max(ev[:, ia].max() - ev[:, ia].min(),
                                       ev[:, ib].max() - ev[:, ib].min()))
            else:
                half = 0.0
            if not (half > 0.0):                       # eval degenerate -> coil
                half = 0.5 * float(max(verts[:, ia].max() - verts[:, ia].min(),
                                       verts[:, ib].max() - verts[:, ib].min()))
            if not (half > 0.0):
                return {"method": "manufacture",
                        "error": "--scan-map: could not derive a scan half-width "
                                 "from the eval region or coil; pass --scan-half "
                                 "(metres)"}
            plane_point = ev_c.astype(float).copy()
            plane_point[axn] = coil_c[axn] + standoff
            scan_current = (float(args.scan_current)
                            if args.scan_current is not None
                            else float(result["best_fit_current_A"]))
            scan = _scan_map_3component(out_chain, scan_current, axn,
                                        plane_point, half, int(args.scan_n))
            # Data Persistence: ALWAYS save the full grid + components.  Prefer
            # next to --output (the panel's result dir); else next to the coil.
            scan_path = (args.scan_output
                         or (os.path.splitext(args.output)[0] + "_scanmap.json"
                             if args.output
                             else os.path.splitext(args.coil_vol)[0]
                             + "_scanmap.json"))
            _write_scan_map_json(scan_path, args, scan)
            _big = ("Bx_T", "By_T", "Bz_T", "grid_u_m", "grid_v_m")
            result["scan_map"] = {k: v for k, v in scan.items()
                                  if k not in _big}
            result["scan_map"]["scan_output"] = scan_path
            result["scan_map"]["standoff_m"] = standoff
            result["scan_map"]["note"] = (
                "3-component planar scan of the delivered single-stroke wire at "
                "current_A -- simulation twin of a 3-axis TMR/fluxgate scan. "
                "Bnormal = along plane_normal (main uniform field); Bt = in-plane "
                "transverse leakage a 3-axis sensor sees (single-stroke "
                "connectors) that a 1-axis probe misses. Full grid saved to "
                "scan_output; compare vs a measured +-I current-reversal scan.")
            result["scan_map_output"] = scan_path
        t_scan = time.perf_counter()

        if args.step_output:
            step_poly = [out_chain] if args.distort else loops
            _write_step_polylines(step_poly, args.step_output)
            result["step"] = args.step_output

        if getattr(args, "former_output", ""):
            # printable FORMER: plate - wire channel (C2 deliverable) -- the
            # water-soluble-support 3D print of the out-of-surface groove the
            # single-stroke wire threads through.  Uses the DELIVERED wire
            # (distorted if --distort).  No-Fallback: negative knobs fail loud.
            for _nm in ("former_clearance", "former_margin", "former_wall",
                        "former_decimate"):
                if float(getattr(args, _nm)) < 0.0:
                    return {"method": "manufacture",
                            "error": f"--{_nm.replace('_', '-')} must be >= 0"}
            t_f0 = time.perf_counter()
            result["former"] = _write_former_step(
                out_chain, args.wire_diam, args.former_clearance,
                args.former_margin, args.former_wall, args.former_decimate,
                args.former_output)
            result["t_former_s"] = round(time.perf_counter() - t_f0, 3)

        if args.peec:
            result["peec"] = _peec_inductance(
                out_chain, args.wire_diam, 5.8e7, args.peec_freq)

        if resonance is not None:
            # accurate achieved L from the FINAL chain (search used a fast nn L)
            L_ach = float(result.get("peec", {}).get("L_H", float("nan")))
            resonance["achieved_inductance_H"] = L_ach
            w_op = 2.0 * np.pi * float(args.peec_freq)
            if resonance.get("mode") == "from_greedy_coil":
                # turns are fixed by --greedy-turns (the low-turn coil); the
                # resonance spec did NOT search turns.  Report the tank capacitor
                # this coil NEEDS to resonate at the operating frequency, and (if
                # a cap was given) the resonance frequency this coil + that cap
                # actually reach -- the "few-turn uniform coil -> required C" answer.
                if L_ach > 0:
                    resonance["required_cap_F"] = float(1.0 / (w_op * w_op * L_ach))
                if cap is not None and L_ach > 0:
                    resonance["resonance_freq_Hz"] = float(
                        1.0 / (2.0 * np.pi * (L_ach * cap) ** 0.5))
                    resonance["resonance_freq_rel_error"] = float(
                        abs(resonance["resonance_freq_Hz"] - args.peec_freq)
                        / (abs(args.peec_freq) + 1e-30))
                resonance["note"] = (
                    "turn count fixed by --greedy-turns (the low-turn coil); the "
                    "resonance spec was NOT used to search turns (L_coil ~ N^2 "
                    "pins the turns).  required_cap_F is the tank capacitor to "
                    "resonate THIS coil at operating_freq_Hz -- size the capacitor "
                    "to it; or drop --greedy-turns to let nlevels be searched for "
                    "L_target.")
            else:
                resonance["achieved_rel_error"] = float(
                    abs(L_ach - resonance["L_target_H"])
                    / (abs(resonance["L_target_H"]) + 1e-30))
                if cap is not None and L_ach > 0:
                    resonance["resonance_cap_F"] = float(cap)
                    resonance["resonance_freq_Hz"] = float(
                        1.0 / (2.0 * np.pi * (L_ach * cap) ** 0.5))
                    resonance["operating_freq_Hz"] = float(args.peec_freq)
                if not resonance["in_range"]:
                    resonance["note"] = (
                        "L_target is OUTSIDE the achievable range "
                        f"{resonance['L_range_H']} H over nlevels "
                        f"{resonance['nlevels_range']}; the closest coil was kept. "
                        "Change the former geometry or the tank capacitor C.")
            result["resonance"] = resonance

        if args.steps_plot:
            try:
                result.update(_steps_plot(
                    loops_signed, chain,
                    out_chain if args.distort else None,
                    args.wire_diam, args.steps_plot, args.target_cf))
            except Exception as e:                 # noqa: BLE001
                result["steps_plot_error"] = str(e)

        if args.flux_plot:
            try:
                result.update(_flux_line_plot(
                    out_chain, result["best_fit_current_A"], args.flux_plot,
                    plane=args.flux_plane))
            except Exception as e:                 # noqa: BLE001
                result["flux_plot_error"] = str(e)

        if args.msh_output:
            try:
                from radia.gmsh_post_export import GmshPostExport
                post = GmshPostExport(coil, boundary=True)
                post.add_field("psi", psi_v, ncomp=1)
                post.write(args.msh_output)
                result["msh"] = args.msh_output
            except Exception as e:                 # noqa: BLE001
                result["msh_error"] = str(e)

        result["t_mesh_s"] = round(t1 - t0, 3)
        result["t_solve_s"] = round(t_psi - t1, 3)
        result["t_chain_s"] = round(t2 - t_psi, 3)
        result["t_distort_s"] = round(t3 - t2, 3)
        result["t_scan_s"] = round(t_scan - t3, 3)
        result["t_total_s"] = round(time.perf_counter() - t0, 3)
    return result


def run(args):
    dispatch = {"design": run_design, "pareto": run_pareto,
                "manufacture": run_manufacture}
    return dispatch[args.method](args)


def build_argparser():
    ap = argparse.ArgumentParser(
        description="Stream-function coil design (design / pareto / manufacture).")
    # ---- shared inputs ----
    ap.add_argument("--coil-vol", required=True,
                    help="coil surface mesh (.vol, 2D-in-3D)")
    ap.add_argument("--eval-vol", required=True,
                    help="evaluation region mesh (.vol; surface or volume)")
    ap.add_argument("--target-cf", default="",
                    help="target field CoefficientFunction expr of x,y,z "
                         "(scalar -> Bz; 3-vector -> B). Give this OR "
                         "--target-harmonic.")
    ap.add_argument("--target-harmonic", default="",
                    help="target as a spherical harmonic (MRI gradient/shim "
                         "basis, l<=4): a name (Z0,X,Y,Z,Z2,ZX,ZY,C2,S2,Z3,"
                         "Z2X,..,Z4,Z3X,Z2C2,ZC3,C4,..) or 'l=L,m=M', "
                         "optionally weighted/summed, e.g. 'Z2', 'X' (Gx), "
                         "'Z2:1.0,Z:0.1'. Generates the solid-harmonic "
                         "--target-cf polynomial.")
    ap.add_argument("--harmonic-lmax", type=int, default=3,
                    help="max degree l for the achieved-Bz spherical-harmonic "
                         "decomposition (purity/contamination; l<=4). 0 = off.")
    ap.add_argument("--method", choices=["design", "pareto", "manufacture"],
                    default="design")
    # ---- solver (design + pareto) ----
    ap.add_argument("--order", type=int, default=3, help="psi FE order")
    ap.add_argument("--regularize", choices=["l2", "h1", "inductance"],
                    default="h1",
                    help="seminorm: l2 (min |psi|), h1 (min |grad psi|, a "
                         "smoothness proxy), or inductance (min 1/2 psi^T L "
                         "psi -- the physical MIN-STORED-ENERGY gradient-coil "
                         "objective via ngsolve.bem; dense, moderate N)")
    ap.add_argument("--confine", choices=["off", "on", "abe"], default="off",
                    help="current-confinement BC. off=none; on=psi 0 on the "
                         "former edge (gradient/shim; breaks solenoids); "
                         "abe=Abe edge-equipotential (one free constant per "
                         "boundary + ground) -- closes the contours AND works "
                         "for gradient and solenoid (recommended).")
    ap.add_argument("--alpha", type=float, default=0.0,
                    help="Tikhonov weight (design; 0 = exact-fit min-seminorm)")
    ap.add_argument("--eval-max", type=int, default=400,
                    help="cap on constraint / measure points (subsample)")
    # ---- pareto mode ----
    ap.add_argument("--pareto-lever",
                    choices=["alpha", "linf", "geometry"],
                    default="alpha",
                    help="front-pushing lever: alpha (Tikhonov L-curve), "
                         "linf (minimax IRLS peak-current trajectory), "
                         "geometry (eval-region scale sweep)")
    ap.add_argument("--alpha-min", type=float, default=1e-4,
                    help="pareto alpha-sweep lower bound (relative)")
    ap.add_argument("--alpha-max", type=float, default=3e1,
                    help="pareto alpha-sweep upper bound (relative)")
    ap.add_argument("--n-alpha", type=int, default=12,
                    help="pareto sweep number of points (alpha/geometry)")
    ap.add_argument("--linf-iter", type=int, default=8,
                    help="linf lever IRLS iterations (pareto)")
    ap.add_argument("--geom-scale-min", type=float, default=0.6,
                    help="geometry lever: eval-region scale lower bound")
    ap.add_argument("--geom-scale-max", type=float, default=1.6,
                    help="geometry lever: eval-region scale upper bound")
    # ---- manufacture mode ----
    ap.add_argument("--nlevels", type=int, default=12,
                    help="number of iso-contours = wire turns (manufacture)")
    ap.add_argument("--optimize-levels", action="store_true",
                    help="manufacture: optimise the N contour LEVELS (vs the "
                         "uniform equal-deltaI default) to minimise the equal-"
                         "current discrete-loop field error -- the few-turn "
                         "quantisation fix (keeps clean contours; complements "
                         "--distort's wire bend).")
    ap.add_argument("--greedy-turns", type=int, default=None,
                    help="LOW-TURN STRATEGY (manufacture): build the coil "
                         "GREEDILY -- lay up to this many turns one at a time, "
                         "each placed where it most reduces the field residual "
                         "(matching pursuit, monotone by construction; the "
                         "antidote to non-monotonic contour quantisation).  "
                         "Emits the rms-vs-turns trace (read off the fewest "
                         "turns meeting any spec).")
    ap.add_argument("--greedy-dict", choices=["contour", "pin", "bubble"],
                    default="contour",
                    help="greedy candidate loops: 'contour' = dense psi iso-"
                         "contour loops (A); 'pin' = on-surface k-ring loops at "
                         "surface pins (B, 'wind a loop on a pin'); 'bubble' = "
                         "bubble-packed pins, density ~ |grad psi| (C).")
    ap.add_argument("--greedy-connector-weight", type=float, default=0.0,
                    help="connector-aware greedy: bias selection toward loops "
                         "NEAR the ones already laid (score = field_residual + "
                         "weight * gap/coil_diag) so the single-stroke connector "
                         "('rung') path stays short -> manufacturable wire.  0 = "
                         "pure-field (default); ~0.3-1.0 trades a little field for "
                         "short rungs (fixes the isolated-pin chaining penalty).")
    ap.add_argument("--greedy-target", type=float, default=None,
                    help="greedy stop threshold: stop adding turns once the "
                         "field residual rms <= this (else use all --greedy-"
                         "turns).")
    ap.add_argument("--greedy-plot", default="",
                    help="turn-budget UI: write the greedy rms-vs-turns curve "
                         "(with the --greedy-target spec line + min-turns mark) "
                         "to this PNG.")
    ap.add_argument("--pin-tiling", action="store_true",
                    help="DENSE BUBBLE-TILING pin coil (manufacture, many-turn): "
                         "equal-current loops packed over the surface with "
                         "density ~ |K| = |grad psi| and sized to their bubble -- "
                         "the discrete-loop, manufacturable alternative to long "
                         "iso-contours (dense where the current is strong).")
    ap.add_argument("--pin-tiling-pins", type=int, default=80,
                    help="pin-tiling / bubble greedy dict: target number of "
                         "bubble loops (k is bisected to land ~this many).")
    ap.add_argument("--pin-tiling-frac", type=float, default=0.7,
                    help="pin-tiling / bubble greedy dict: loop radius as a "
                         "fraction of the bubble clearance radius (0.7 -> loops "
                         "tile without heavy overlap).")
    ap.add_argument("--target-inductance", type=float, default=None,
                    help="IH-RESONANCE design: search the turn count (nlevels) "
                         "so the single-stroke coil inductance equals this "
                         "value in Henries (the field design is unchanged; "
                         "L_coil ~ N^2 sets the resonance).  Forces --peec.  "
                         "WITH --greedy-turns: turns are fixed by greedy (the "
                         "low-turn coil), so nlevels is NOT searched; instead "
                         "the result reports required_cap_F (the tank capacitor "
                         "to resonate the delivered coil at --peec-freq).")
    ap.add_argument("--resonance-cap", type=float, default=None,
                    help="IH-RESONANCE design via the tank capacitor in Farads: "
                         "target L = 1/((2 pi f)^2 C) with f = --peec-freq "
                         "(the inverter frequency).  Alternative to "
                         "--target-inductance.  WITH --greedy-turns: turns are "
                         "fixed by greedy, so the result reports the resonance "
                         "frequency THIS coil + this cap reach (resonance_freq_Hz)"
                         " plus required_cap_F to hit --peec-freq exactly.")
    ap.add_argument("--nlevels-max", type=int, default=60,
                    help="upper bound for the IH-resonance nlevels search")
    ap.add_argument("--contour-sub", type=int, default=1,
                    help="order-p contour: subdivide each surface triangle "
                         "sub x sub and evaluate the full-order psi (sub=1 = "
                         "vertex-linear; sub=3 follows the curved order-2/3 "
                         "psi within each element) (manufacture)")
    ap.add_argument("--chain", choices=["field_aware", "nn"],
                    default="field_aware",
                    help="single-stroke chaining: field_aware (min net "
                         "connector field; reaches the separate-turns floor "
                         "for closed contours, no distort) or nn (nearest "
                         "neighbour)")
    ap.add_argument("--chain-ncut", type=int, default=0,
                    help="field_aware cut candidates per loop (0 = auto-scale "
                         "with loop count; raise for clustered/many-loop coils)")
    ap.add_argument("--chain-passes", type=int, default=0,
                    help="field_aware cut coordinate-descent rounds (0 = auto-"
                         "scale with loop count)")
    ap.add_argument("--distort", action="store_true",
                    help="single-current sheet-metal wire distortion (manufacture)")
    ap.add_argument("--distort-grid", type=int, default=3,
                    help="sheet-metal 3D control-grid size per axis (manufacture)")
    ap.add_argument("--distort-iter", type=int, default=5,
                    help="sheet-metal Gauss-Newton iterations (manufacture)")
    ap.add_argument("--step-output", default="",
                    help="STEP CAD output path for the wire (manufacture)")
    # ---- printable former (3D-print deliverable) ----
    ap.add_argument("--former-output", default="",
                    help="printable FORMER STEP (manufacture): a base plate with "
                         "the DELIVERED wire CHANNEL cut out (plate - swept "
                         "tube), for water-soluble-support 3D printing -- "
                         "dissolve the support to leave the out-of-surface "
                         "groove the single-stroke wire threads through. "
                         "Distinct from --step-output (the wire centreline). "
                         "STEP coords in metres; scale in the slicer.")
    ap.add_argument("--former-clearance", type=float, default=3e-4,
                    help="former channel radial clearance in m (channel radius = "
                         "wire-diam/2 + clearance; default 0.3 mm)")
    ap.add_argument("--former-margin", type=float, default=5e-3,
                    help="former plate XY margin beyond the wire footprint in m "
                         "(default 5 mm)")
    ap.add_argument("--former-wall", type=float, default=3e-3,
                    help="former minimum wall thickness around the channel in m "
                         "(sets the plate Z half-thickness beyond the wire "
                         "z-extent; default 3 mm)")
    ap.add_argument("--former-decimate", type=float, default=2.5e-3,
                    help="former channel spine resample step in m (coarsens the "
                         "wire polyline to bound the CAD primitive count / build "
                         "time; 0 = no decimation; default 2.5 mm)")
    ap.add_argument("--peec", action="store_true",
                    help="compute PEEC coil inductance L, R (manufacture)")
    ap.add_argument("--wire-diam", type=float, default=1e-3,
                    help="PEEC wire cross-section size in m (manufacture)")
    ap.add_argument("--peec-freq", type=float, default=1e5,
                    help="PEEC evaluation frequency in Hz (manufacture)")
    ap.add_argument("--flux-plot", default="",
                    help="bubble-system flux-line PNG of the coil field on a "
                         "cut-plane (manufacture)")
    ap.add_argument("--flux-plane", choices=["x", "y", "z"], default="y",
                    help="cut-plane normal for --flux-plot (manufacture)")
    ap.add_argument("--steps-plot", default="",
                    help="per-step manufacturing PNG (contours -> single-stroke "
                         "-> sheet-metal -> with thickness) (manufacture)")
    # ---- 3-component scan map (3-axis magnetometer twin) ----
    ap.add_argument("--scan-map", action="store_true",
                    help="MEASUREMENT TWIN (manufacture): 3-component (Bx,By,Bz) "
                         "planar raster scan of the delivered single-stroke wire "
                         "-- the simulation counterpart of a 3-axis magnetometer "
                         "(TMR/fluxgate) scan. Emits Bnormal (main uniform field), "
                         "Bt (in-plane transverse leakage from single-stroke "
                         "connectors -- what a 1-axis probe misses), homogeneity, "
                         "and saves the full grid to a JSON to compare against a "
                         "measured +-I current-reversal scan.")
    ap.add_argument("--scan-plane", choices=["x", "y", "z"], default="z",
                    help="scan-plane normal axis (default z: an x-y plane above a "
                         "planar coil). Bnormal is the component along this axis.")
    ap.add_argument("--scan-standoff", type=float, default=None,
                    help="scan-plane offset (m) along the normal from the coil "
                         "centroid = the sensor standoff height. Default: the "
                         "eval-region (DSV) height, so the scan sits where the "
                         "design targets uniformity.")
    ap.add_argument("--scan-half", type=float, default=None,
                    help="scan half-width (m) in the two in-plane axes (square "
                         "grid). Default: the eval-region in-plane half-extent "
                         "(the DSV footprint).")
    ap.add_argument("--scan-n", type=int, default=21,
                    help="scan raster points per axis (n x n grid; default 21)")
    ap.add_argument("--scan-current", type=float, default=None,
                    help="physical drive current (A) the map is scaled to (field "
                         "is linear in current). Default: the best-fit design "
                         "current. Set the bench current to match the sensor "
                         "reading directly.")
    ap.add_argument("--scan-output", default="",
                    help="JSON path for the full 3-component scan grid. Default: "
                         "next to --output (<base>_scanmap.json).")
    # ---- active shielding ----
    ap.add_argument("--iron-vol", default=None,
                    help="iron/material .vol (iron region + 'kelvin' ball + "
                         "periodic kelvin_int/kelvin_ext + 'GND' vertex -- the "
                         "calc_fem_kelvin pipeline). When given, the design is "
                         "MATERIAL-AWARE: the iron reaction (reduced scalar "
                         "potential, Kelvin open boundary) is ADDED to the free-"
                         "space kernel so a coil designed near a yoke/shield/"
                         "core hits the target IN the iron system. Free-space "
                         "by default (omit this flag).")
    ap.add_argument("--mu-r", type=float, default=1000.0,
                    help="relative permeability of the --iron-vol iron region "
                         "(linear; default 1000).")
    ap.add_argument("--iron-mat", default="iron",
                    help="material name of the iron region in --iron-vol "
                         "(default 'iron').")
    ap.add_argument("--iron-exact-source", action="store_true",
                    help="material-aware kernel: evaluate the coil source field "
                         "at the iron QUADRATURE points (analytic Biot-Savart, "
                         "no P1-vertex interpolation -> removes the source "
                         "plateau, ~14%% on a coarse shell).  Much heavier (one "
                         "assembly per quad point); opt-in.  Default is the fast "
                         "P1-vertex source.")
    ap.add_argument("--iron-quad-order", type=int, default=None,
                    help="quadrature order for --iron-exact-source (default "
                         "2*fem_order).")
    ap.add_argument("--shield-vol", default=None,
                    help="outer shield coil surface mesh (.vol). When given, "
                         "designs a PRIMARY + SHIELD coil jointly (Turner "
                         "1986 method): the shield cancels the stray field "
                         "outside --shield-eval-vol while primary + shield "
                         "together produce the target field inside --eval-vol "
                         "(DSV). Only 'design' mode; requires "
                         "--shield-eval-vol.")
    ap.add_argument("--shield-eval-vol", default=None,
                    help="external evaluation region (.vol) WHERE THE STRAY "
                         "FIELD MUST VANISH (used only with --shield-vol). "
                         "Typically a thin shell just outside the shield "
                         "cylinder.")
    ap.add_argument("--shield-weight", type=float, default=1.0,
                    help="weight on the stray-field null-constraint rows "
                         "relative to the DSV target rows (default 1.0). "
                         "Increase to push stray suppression; decrease if "
                         "DSV homogeneity degrades too much.")
    # ---- output ----
    ap.add_argument("--msh-output", default="",
                    help="optional GMSH .msh of psi on the coil surface")
    ap.add_argument("--output", default="", help="JSON output file")
    return ap


def main():
    return calc_main(run, build_argparser())


if __name__ == "__main__":
    main()
