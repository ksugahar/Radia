"""
Heat-only IH runtime solver driven by a precomputed 2D EM table.

Replaces the inner Karl-FEM EM iteration in ``calc_fem_kelvin.py``
(``--impedance esim``) for the engineering use case where the EM
problem is quasi-static at the heat timescale (heat tau ~ seconds vs
EM 1/omega ~ 100 us at 10 kHz -> ratio ~10^4).  Per timestep the
solver:

  1. Computes the current coil amplitude I(t) from the trajectory.
  2. Scales the reference surface tangential field |H_t_ref(r)| by
     I(t)/I_ref (linear in the weak-saturation limit -- the table
     itself captures any |H_t|-dependent Z_s nonlinearity).
  3. For each surface DOF, looks up q_surf_i = interp_qsurf(table,
     |H_t(r_i, t)|, T(r_i, t)) in the 2D table from calc_em_table.py.
  4. Steps the heat PDE forward (backward Euler) with the new q_surf
     as Neumann RHS.

Two sources for the reference |H_t_ref(r)| field on the workpiece
surface:

  --ht-source kelvin --ht-sol <Jsurf.sol> --em-vol <fem.vol>
        Use the ``<stem>_Jsurf.sol`` file produced by
        calc_fem_kelvin.py (which is exactly |J_s| = |H_t| on the
        SIBC face, H1 scalar).  Captures workpiece backreaction.
        Recommended.  Requires one upstream calc_fem_kelvin run at
        the reference coil current.

  --ht-source biot  --coil-step <step> [--biot-image-factor 2.0]
        Pure incident Biot-Savart |H_t| from the coil STEP, with NO
        workpiece backreaction.  The coil centerline is walked from
        the STEP solid (same extractor as the PEEC path) into one
        filament; the incident field is evaluated at each workpiece
        surface DOF and projected onto the local tangent plane.  Needs
        no upstream calc_fem_kelvin run, but UNDER-predicts the surface
        field by ~2x for a good flat conductor (the eddy-current image
        that doubles H_t is absent).  --biot-image-factor lets an
        informed user opt into that doubling explicitly; it is never
        applied silently.  Prefer --ht-source kelvin when the
        backreaction matters.

Coil current trajectory:
  --coil-current <scalar>            Constant amplitude (A_peak)
  --coil-current-csv <2col.csv>      ``t_s, I_A_peak``
                                     (linear interp, clamp ends)

This module is a pragmatic *engineering* replacement for the
sigma(T) coupling research-track that was dropped 2026-05-24.  Its
single dominant simplification is that the SPATIAL distribution of
|H_t(r)| is frozen at the reference solve; only the AMPLITUDE
follows I(t).  This is exact in the linear regime and a good
approximation when sigma(T) does not move enough to substantially
redistribute the surface currents.  When that assumption breaks
(strong Curie-point band moving across the workpiece), re-run
calc_fem_kelvin at a representative T to refresh ``--ht-sol``.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time

import numpy as np

_panels_dir = os.path.dirname(os.path.abspath(__file__))
if _panels_dir not in sys.path:
    sys.path.insert(0, _panels_dir)

from calc_common import setup_paths, progress, calc_main  # noqa: E402
from em_table import load_em_table, interp_qsurf  # noqa: E402
from calc_heat import THERMAL_PRESETS, _resolve_material, SIGMA_SB  # noqa: E402


def _log(msg):
    progress("HEAT_EMTAB", msg)


def _load_current_trajectory(scalar_val, csv_path):
    """Return a callable ``I(t)`` [A_peak] from either a scalar or CSV.

    CSV format: 2 cols ``t_s, I_A``.  Clamps outside the range
    (matches the sigma-T-curve convention).
    """
    if csv_path:
        data = np.loadtxt(csv_path, delimiter=",", comments="#")
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError(
                f"--coil-current-csv must be 2-col (t_s, I_A); "
                f"got shape {data.shape}")
        t_arr = data[:, 0].astype(float)
        I_arr = data[:, 1].astype(float)
        order = np.argsort(t_arr)
        t_arr = t_arr[order]
        I_arr = I_arr[order]

        def I_of_t(t):
            return float(np.interp(t, t_arr, I_arr,
                                   left=I_arr[0], right=I_arr[-1]))

        I_of_t.t_min = float(t_arr[0])
        I_of_t.t_max = float(t_arr[-1])
        I_of_t.kind = "csv"
        return I_of_t

    if scalar_val is None:
        raise ValueError(
            "Either --coil-current SCALAR or --coil-current-csv PATH "
            "is required.")
    val = float(scalar_val)

    def I_const(_t):
        return val

    I_const.t_min = 0.0
    I_const.t_max = float("inf")
    I_const.kind = "constant"
    return I_const


def _collect_surface_dofs(wp_mesh, surface_label_eff, fes_T):
    """Return ``(dof_indices, dof_xyz)`` for the H1 DOFs lying on the
    heating face.  H1 DOFs are vertices, so this is just the BND
    vertex set under the requested label.
    """
    from ngsolve import BND

    vnrs = set()
    for el in wp_mesh.Elements(BND):
        if surface_label_eff != ".*" and el.mat != surface_label_eff:
            continue
        for v in el.vertices:
            vnrs.add(v.nr)
    vnrs_sorted = sorted(vnrs)
    xyz = np.array([
        [float(c) for c in wp_mesh.vertices[vnr].point]
        for vnr in vnrs_sorted
    ], dtype=float) if vnrs_sorted else np.zeros((0, 3))
    return np.asarray(vnrs_sorted, dtype=np.int64), xyz


def _project_Ht_ref(wp_mesh, dof_xyz, em_vol, ht_sol, ht_order):
    """Sample |H_t_ref(r)| at each wp-surface DOF by pointwise
    evaluation of the EM-mesh scalar H1 GridFunction.

    Returns an (n_dof,) float array.  Vertices that fall outside the
    EM mesh are tagged with NaN; the caller masks them.
    """
    from ngsolve import Mesh, H1, GridFunction

    em_mesh = Mesh(em_vol)
    fes_J = H1(em_mesh, order=ht_order)
    gf_J = GridFunction(fes_J)
    gf_J.Load(ht_sol)

    n = dof_xyz.shape[0]
    Ht_ref = np.full(n, np.nan, dtype=float)
    for i in range(n):
        x, y, z = dof_xyz[i]
        try:
            mip = em_mesh(float(x), float(y), float(z))
            val = gf_J(mip)
            Ht_ref[i] = float(getattr(val, "real", val))
        except Exception:
            pass
    n_ok = int(np.isfinite(Ht_ref).sum())
    return Ht_ref, n_ok


def _polyline_to_filament_segments(poly, closed):
    """(M,3) centerline polyline -> (N_seg, 2, 3) filament endpoint pairs.

    Consecutive points become straight filaments; when ``closed`` the
    wrap-around segment ``(poly[-1], poly[0])`` is appended so the loop
    carries a divergence-free current.  Pure NumPy -- the unit-testable
    core of the STEP centerline -> Biot-Savart bridge.
    """
    poly = np.asarray(poly, dtype=float)
    if poly.ndim != 2 or poly.shape[0] < 2 or poly.shape[1] != 3:
        raise ValueError(
            f"coil centerline has shape {poly.shape}; need >=2 points "
            f"of (x,y,z).")
    segs = np.stack([poly[:-1], poly[1:]], axis=1)        # (M-1, 2, 3)
    if bool(closed):
        wrap = np.array([[poly[-1], poly[0]]], dtype=float)
        segs = np.concatenate([segs, wrap], axis=0)
    return segs


def _coil_segments_from_step(coil_step):
    """Return ``(segments, closed)`` for the coil centerline.

    ``segments`` is an (N_seg, 2, 3) float array of filament endpoints
    suitable for ``radia.biot_savart.h_segments_batch``.  This reuses
    the same ``extract_centerline`` walker the PEEC path uses, so the
    coil is described by ONE filament tracing the conductor centerline
    (a multi-turn solid traces all turns, giving the correct
    amp-turns).
    """
    from radia.coil_from_step import extract_centerline

    res = extract_centerline(coil_step)
    closed = bool(res.closed)
    return _polyline_to_filament_segments(res.polyline, closed), closed


def _surface_vertex_normals(wp_mesh, surface_label_eff, dof_vnrs):
    """Area-weighted unit surface normals aligned to ``dof_vnrs``.

    Returns an (n, 3) array.  Each boundary polygon's Newell normal
    (magnitude ~ 2*area) is accumulated to its vertices, then
    normalized -- so this is the area-weighted average face normal at
    each surface DOF.  Works for tri and quad boundary elements.

    Only the tangent PLANE matters for |H_t| = |H - (H.n)n|, so the
    normal SIGN (inward vs outward) is irrelevant.  A DOF with no
    incident boundary polygon keeps a zero normal, which makes
    |H_t| = |H| there (no tangent plane is known -- the full incident
    magnitude is the safe choice).
    """
    from ngsolve import BND

    nv = wp_mesh.nv
    acc = np.zeros((nv, 3), dtype=float)
    for el in wp_mesh.Elements(BND):
        if surface_label_eff != ".*" and el.mat != surface_label_eff:
            continue
        vs = [v.nr for v in el.vertices]
        m = len(vs)
        if m < 3:
            continue
        pts = np.array([list(wp_mesh.vertices[v].point) for v in vs],
                       dtype=float)
        nrm = np.zeros(3)
        for i in range(m):                       # Newell's method
            a = pts[i]
            b = pts[(i + 1) % m]
            nrm[0] += (a[1] - b[1]) * (a[2] + b[2])
            nrm[1] += (a[2] - b[2]) * (a[0] + b[0])
            nrm[2] += (a[0] - b[0]) * (a[1] + b[1])
        for v in vs:
            acc[v] += nrm
    out = np.zeros((len(dof_vnrs), 3), dtype=float)
    for k, vnr in enumerate(dof_vnrs):
        n = acc[int(vnr)]
        ln = float(np.linalg.norm(n))
        if ln > 1e-30:
            out[k] = n / ln
    return out


def _biot_Ht_on_surface(segments, obs_xyz, obs_normals, current,
                        image_factor=1.0):
    """|H_t| [A/m] at each obs point from coil ``segments`` (carrying
    ``current`` [A]), projected onto the tangent plane of ``obs_normals``.

    This is the pure, mesh-free core of the biot ht-source: it is the
    INCIDENT tangential field (no workpiece backreaction).  For a good
    conductor the eddy currents cancel the normal component and roughly
    double the tangential one (the perfect-conductor image), so the raw
    incident |H_t| under-predicts the true surface field by ~2x for a
    flat surface.  ``image_factor`` (default 1.0 = raw incident) lets an
    informed caller opt into that doubling explicitly; it is NEVER
    applied silently.

    Parameters
    ----------
    segments : (N_seg, 2, 3) array of filament endpoints [m]
    obs_xyz  : (n, 3) observation points [m]
    obs_normals : (n, 3) unit surface normals (zero -> no projection)
    current  : coil current [A]
    image_factor : multiply the result by this (default 1.0)

    Returns
    -------
    (n,) float |H_t| [A/m]
    """
    from radia.biot_savart import h_segments_batch

    H = h_segments_batch(segments, np.asarray(obs_xyz, dtype=float),
                         current=float(current))            # (n, 3)
    n_hat = np.asarray(obs_normals, dtype=float)
    Hn = np.sum(H * n_hat, axis=1)                          # (n,)
    H_t = H - Hn[:, None] * n_hat                           # tangential
    return float(image_factor) * np.linalg.norm(H_t, axis=1)


def solve_heat_em_table(wp_vol, em_table_path,
                         ht_source, ht_sol, em_vol, ht_order,
                         coil_current, coil_current_csv,
                         I_ref,
                         coil_step="", biot_image_factor=1.0,
                         material="steel", rho=None, cp=None, k=None,
                         h_conv=10.0, t_ext=20.0, t_initial=20.0, emissivity=0.0,
                         surface_label="",
                         dt=0.5, t_end=60.0,
                         fes_order=1,
                         linear_solver="sparsecholesky",
                         probe_point=None,
                         csv_output=""):
    """Backward-Euler heat solve with per-step 2D-table q_surf lookup."""
    setup_paths()
    t0 = time.perf_counter()

    from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                         Integrate, CF, ds, dx, BND, TaskManager,
                         InnerProduct, grad)

    if not os.path.isfile(wp_vol):
        return {"error": f"--wp-vol not found: {wp_vol}"}
    if not os.path.isfile(em_table_path):
        return {"error": f"--em-table not found: {em_table_path}"}

    tab = load_em_table(em_table_path)
    _log(f"TABLE:{os.path.basename(em_table_path)} "
         f"material={tab.material} f={tab.frequency:g}Hz "
         f"H_range=[{tab.H_grid[0]:.2e},{tab.H_grid[-1]:.2e}] A/m "
         f"T_range=[{tab.T_grid[0]:.1f},{tab.T_grid[-1]:.1f}] C")

    wp_mesh = Mesh(wp_vol).Curve(int(fes_order))
    _log(f"MESH:loaded {os.path.basename(wp_vol)} "
         f"materials={list(wp_mesh.GetMaterials())} "
         f"boundaries={list(wp_mesh.GetBoundaries())}")

    if surface_label:
        if surface_label not in wp_mesh.GetBoundaries():
            return {"error":
                    f"--surface-label {surface_label!r} not found in "
                    f"{wp_vol} boundaries={list(wp_mesh.GetBoundaries())}"}
        surface_label_eff = surface_label
    else:
        surface_label_eff = ".*"
    _log(f"BND:filter={surface_label_eff!r}")

    rho_v, cp_v, k_v = _resolve_material(material, rho, cp, k)
    _log(f"MATERIAL:{material} rho={rho_v} cp={cp_v} k={k_v}")

    fes_T = H1(wp_mesh, order=int(fes_order))
    u, v = fes_T.TnT()
    gfT = GridFunction(fes_T)
    gfT.vec[:] = float(t_initial)
    _log(f"FES:H1 order={fes_order} ndof={fes_T.ndof}")

    # -------------------- |H_t_ref(r)| --------------------
    dof_vnrs, dof_xyz = _collect_surface_dofs(
        wp_mesh, surface_label_eff, fes_T)
    n_surf = len(dof_vnrs)
    _log(f"BND_DOF:{n_surf} surface DOFs on filter {surface_label_eff!r}")

    if ht_source == "kelvin":
        if not ht_sol or not em_vol:
            return {"error":
                    "--ht-source kelvin requires both --ht-sol "
                    "(<stem>_Jsurf.sol from calc_fem_kelvin) and "
                    "--em-vol (<stem>_fem.vol the .sol was saved on)."}
        if not os.path.isfile(ht_sol):
            return {"error": f"--ht-sol not found: {ht_sol}"}
        if not os.path.isfile(em_vol):
            return {"error": f"--em-vol not found: {em_vol}"}
        Ht_ref_arr, n_ok = _project_Ht_ref(
            wp_mesh, dof_xyz, em_vol, ht_sol, ht_order)
        _log(f"HT_REF:projected {n_ok}/{n_surf} surface DOFs from "
             f"{os.path.basename(ht_sol)} (NaN tagged outside)")
        if n_ok == 0:
            return {"error":
                    "--ht-sol projection produced zero valid DOFs.  "
                    "Check that --em-vol matches the .sol that was saved."}
        # Replace NaN with 0 so the lookup returns 0 (clamped to H_min
        # which gives ~0 q_surf for steel BH below the first tabulated H).
        Ht_ref_arr = np.where(np.isfinite(Ht_ref_arr), Ht_ref_arr, 0.0)
    elif ht_source == "biot":
        if not coil_step:
            return {"error":
                    "--ht-source biot requires --coil-step <coil.step> "
                    "(the coil centerline is walked from the STEP solid)."}
        if not os.path.isfile(coil_step):
            return {"error": f"--coil-step not found: {coil_step}"}
        try:
            segs, closed = _coil_segments_from_step(coil_step)
        except Exception as e:
            return {"error":
                    f"coil centerline extraction failed for "
                    f"{os.path.basename(coil_step)}: {type(e).__name__}: "
                    f"{e}.  Use --ht-source kelvin, or regenerate the "
                    f"coil STEP with a cleaner swept cross-section."}
        _log(f"BIOT:centerline {segs.shape[0]} segments closed={closed} "
             f"from {os.path.basename(coil_step)}")
        normals = _surface_vertex_normals(
            wp_mesh, surface_label_eff, dof_vnrs)
        n_proj = int((np.linalg.norm(normals, axis=1) > 1e-9).sum())
        Ht_ref_arr = _biot_Ht_on_surface(
            segs, dof_xyz, normals, current=I_ref,
            image_factor=biot_image_factor)
        _log(f"HT_REF:biot incident |H_t| at {n_surf} surface DOFs "
             f"({n_proj} tangent-projected) I_ref={I_ref:.3f} A "
             f"image_factor={biot_image_factor:g}; |H_t| range "
             f"[{float(np.min(Ht_ref_arr)):.3e},"
             f"{float(np.max(Ht_ref_arr)):.3e}] A/m")
        # Defensive: zero out any non-finite (degenerate geometry).
        Ht_ref_arr = np.where(np.isfinite(Ht_ref_arr), Ht_ref_arr, 0.0)
    else:
        return {"error": f"Unknown --ht-source {ht_source!r}"}

    # Coil current trajectory.
    I_of_t = _load_current_trajectory(coil_current, coil_current_csv)
    _log(f"COIL:I_ref={I_ref:.3f} A trajectory={I_of_t.kind}")

    # -------------------- Forms --------------------
    K_cf = CF(float(k_v))
    rho_cp = CF(float(rho_v) * float(cp_v))

    a_form = BilinearForm(fes_T, symmetric=True)
    a_form += K_cf * InnerProduct(grad(u), grad(v)) * dx
    a_form += float(h_conv) * v * u * ds(surface_label_eff)
    a_form.Assemble()

    m_form = BilinearForm(fes_T, symmetric=True)
    m_form += rho_cp * u * v * dx
    m_form.Assemble()

    # Backward Euler: mstar = M + dt * K.
    mstar = m_form.mat.CreateMatrix()
    mstar.AsVector().data = (
        m_form.mat.AsVector() + float(dt) * a_form.mat.AsVector())
    inv = mstar.Inverse(freedofs=fes_T.FreeDofs(),
                        inverse=linear_solver)
    res_vec = gfT.vec.CreateVector()
    _log(f"SOLVER:{linear_solver} dt={dt} t_end={t_end}")

    # gf_q is the running q_surf GridFunction backing the LinearForm.
    fes_q = H1(wp_mesh, order=int(fes_order))
    gf_q = GridFunction(fes_q)
    gf_q.vec[:] = 0
    qv_fv = gf_q.vec.FV()
    Tv_fv = gfT.vec.FV()
    surface_region = wp_mesh.Boundaries(surface_label_eff)
    A_surf = float(Integrate(CF(1), wp_mesh, BND,
                              definedon=surface_region).real)

    # -------------------- Time loop --------------------
    n_steps = int(math.ceil(float(t_end) / float(dt)))
    t_arr = [0.0]
    T_probe_hist = []
    P_total_hist = []
    T_avg_hist = []
    T_max_hist = []

    if probe_point is not None:
        try:
            mip = wp_mesh(*[float(c) for c in probe_point])
            T_probe_hist.append(float(getattr(gfT(mip), "real", gfT(mip))))
        except Exception:
            T_probe_hist.append(float("nan"))

    Q_input_J = 0.0
    for step in range(1, n_steps + 1):
        t = step * float(dt)

        # 1. Refresh q_surf at the surface DOFs.
        I_now = I_of_t(t)
        amp_scale = I_now / max(I_ref, 1e-30)
        Ht_now = Ht_ref_arr * amp_scale  # (n_surf,)
        # T at each surface DOF: H1-on-volume DOFs include all vertices,
        # and dof_vnrs ARE H1 DOFs for fes_T (order=fes_order; ndof=
        # ngVertices + edge/face DOFs, but surface VERTEX DOFs come first).
        T_dof = np.asarray([float(Tv_fv[int(idx)]) for idx in dof_vnrs])
        q_dof = interp_qsurf(tab, Ht_now, T_dof)

        qv_fv[:] = 0.0
        for k_dof, idx in enumerate(dof_vnrs):
            qv_fv[int(idx)] = float(q_dof[k_dof])

        f_form = LinearForm(fes_T)
        f_form += gf_q * v * ds(surface_label_eff)
        f_form += float(h_conv) * float(t_ext) * v * ds(surface_label_eff)
        if float(emissivity) > 0.0:        # radiation (explicit, prev-step T, in K)
            _TK = gfT + 273.15
            f_form += -float(emissivity) * SIGMA_SB \
                * (_TK**4 - (float(t_ext) + 273.15)**4) * v * ds(surface_label_eff)
        f_form.Assemble()
        with TaskManager():
            res_vec.data = f_form.vec - a_form.mat * gfT.vec
            gfT.vec.data += float(dt) * (inv * res_vec)

        q_int = float(Integrate(gf_q, wp_mesh, BND,
                                 definedon=surface_region).real)
        Q_input_J += q_int * float(dt)

        T_vol = np.asarray(gfT.vec.FV().NumPy())
        T_max_now = float(np.max(T_vol))
        T_min_now = float(np.min(T_vol))
        # T_avg is the volume mean of the FE DOFs -- with H1 order=1 on
        # tets this is close to the mass-weighted spatial mean (off by
        # the regular-tet correction factor 1.0 / (1 - 1/(d+1))).  Good
        # enough as a tracking diagnostic; full mass-weighted average
        # if the user enables --integrate-T-avg in v2.
        T_avg_now = float(np.mean(T_vol))

        t_arr.append(t)
        T_avg_hist.append(T_avg_now)
        T_max_hist.append(T_max_now)
        P_total_hist.append(q_int)

        if probe_point is not None:
            try:
                mip = wp_mesh(*[float(c) for c in probe_point])
                val = gfT(mip)
                T_probe_hist.append(float(getattr(val, "real", val)))
            except Exception:
                T_probe_hist.append(float("nan"))
        _log(f"STEP:{step}/{n_steps} t={t:.3f}s I={I_now:.2f}A "
             f"P_in={q_int:.3e}W T_avg={T_avg_now:.2f}C "
             f"T_max={T_max_now:.2f}C")

    T_arr_final = np.asarray(gfT.vec.FV().NumPy())
    t_total = time.perf_counter() - t0
    _log(f"DONE:T_max={float(np.max(T_arr_final)):.2f}C "
         f"T_avg={float(np.mean(T_arr_final)):.2f}C "
         f"Q_in={Q_input_J:.4e}J t={t_total:.1f}s")

    # Optional CSV history.
    if csv_output:
        try:
            import csv
            with open(csv_output, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                w.writerow(["t_s", "I_A", "P_W", "T_avg_C", "T_max_C"])
                w.writerow([0.0, float(I_of_t(0.0)), 0.0,
                            float(t_initial), float(t_initial)])
                for ti, Pi, Tai, Tmi in zip(t_arr[1:], P_total_hist,
                                             T_avg_hist, T_max_hist):
                    w.writerow([f"{ti:.6f}", f"{I_of_t(ti):.6f}",
                                f"{Pi:.6f}", f"{Tai:.6f}", f"{Tmi:.6f}"])
            _log(f"CSV:wrote {os.path.basename(csv_output)}")
        except Exception as e:
            _log(f"CSV_ERROR:{type(e).__name__}: {e}")

    return {
        "T_max_C": float(np.max(T_arr_final)),
        "T_min_C": float(np.min(T_arr_final)),
        "T_avg_C": float(np.mean(T_arr_final)),
        "T_initial_C": float(t_initial),
        "T_probe_history_C": T_probe_hist if probe_point is not None else None,
        "t_history_s": t_arr,
        "P_total_history_W": [0.0] + P_total_hist,
        "T_avg_history_C": [float(t_initial)] + T_avg_hist,
        "T_max_history_C": [float(t_initial)] + T_max_hist,
        "Q_input_J": Q_input_J,
        "surface_area_m2": A_surf,
        "n_surface_dofs": int(n_surf),
        "n_steps": int(n_steps),
        "dt_s": float(dt),
        "t_end_s": float(t_end),
        "I_ref_A": float(I_ref),
        "ht_source": ht_source,
        "coil_step": (os.path.abspath(coil_step)
                      if (ht_source == "biot" and coil_step) else None),
        "biot_image_factor": (float(biot_image_factor)
                              if ht_source == "biot" else None),
        "em_table": os.path.abspath(em_table_path),
        "wp_vol": os.path.abspath(wp_vol),
        "material": material,
        "ndof": int(fes_T.ndof),
        "fes_order": int(fes_order),
        "t_total_s": round(t_total, 2),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Heat-only IH solver driven by a 2D EM table "
                    "(|H_t|, T) -> Z_s (no per-step Karl FEM iteration).")
    parser.add_argument("--wp-vol", required=True,
                        help="Workpiece volume mesh (.vol).")
    parser.add_argument("--em-table", required=True,
                        help="2D EM table .npz from calc_em_table.py.")

    parser.add_argument("--ht-source", default="kelvin",
                        choices=["kelvin", "biot"],
                        help="Source for |H_t_ref(r)| on the workpiece "
                             "surface.  kelvin: load <stem>_Jsurf.sol "
                             "from a calc_fem_kelvin reference run "
                             "(recommended -- captures backreaction).  "
                             "biot: pure incident Biot-Savart |H_t| from "
                             "--coil-step (no backreaction; ~2x low for a "
                             "good flat conductor unless --biot-image-factor "
                             "is set).")
    parser.add_argument("--ht-sol", default="",
                        help="kelvin path: <stem>_Jsurf.sol from "
                             "calc_fem_kelvin.py (= |H_t| on the SIBC "
                             "face).")
    parser.add_argument("--em-vol", default="",
                        help="kelvin path: <stem>_fem.vol the --ht-sol "
                             "was saved on.")
    parser.add_argument("--ht-order", type=int, default=1,
                        help="H1 order used when calc_fem_kelvin saved "
                             "<stem>_Jsurf.sol (must match fes-order of "
                             "that run).")
    parser.add_argument("--coil-step", default="",
                        help="biot path: coil STEP solid.  Its centerline "
                             "is walked (same extractor as the PEEC path) "
                             "into one filament; the incident Biot-Savart "
                             "|H_t| on the workpiece surface is the "
                             "reference field (NO workpiece backreaction).")
    parser.add_argument("--biot-image-factor", type=float, default=1.0,
                        help="biot path: multiply the incident |H_t| by "
                             "this factor (default 1.0 = raw incident).  "
                             "A good conductor's eddy currents roughly "
                             "DOUBLE the tangential field (perfect-conductor "
                             "image), so set 2.0 for a flat-surface "
                             "good-conductor estimate.  Applied explicitly, "
                             "never silently -- prefer --ht-source kelvin "
                             "when backreaction matters.")

    parser.add_argument("--coil-current", type=float, default=None,
                        help="Constant coil current amplitude [A_peak].")
    parser.add_argument("--coil-current-csv", default="",
                        help="2-col CSV (t_s, I_A_peak) for time-varying "
                             "current.  Overrides --coil-current.")
    parser.add_argument("--I-ref", type=float, required=True,
                        help="Coil current [A_peak] at which --ht-sol "
                             "was solved.  Used to scale the spatial "
                             "|H_t_ref(r)| by I(t)/I_ref each step.")

    # Thermal material.
    parser.add_argument("--material", default="steel",
                        choices=list(THERMAL_PRESETS) + ["custom"],
                        help="Thermal material preset.")
    parser.add_argument("--rho", type=float, default=None)
    parser.add_argument("--cp", type=float, default=None)
    parser.add_argument("--k", type=float, default=None)

    # BCs / IC.
    parser.add_argument("--h-conv", type=float, default=10.0)
    parser.add_argument("--t-ext", type=float, default=20.0)
    parser.add_argument("--emissivity", type=float, default=0.0,
                        help="Surface emissivity for radiation "
                             "eps*sigma*(T^4-T_ext^4) [0..1]; 0 = off "
                             "(radiation ambient = --t-ext).")
    parser.add_argument("--t-initial", type=float, default=20.0,
                        help="Initial workpiece temperature [degC].  "
                             "Same flag name as the --initial-T in the "
                             "spec; CLI uses the calc_heat-compatible "
                             "form.")

    # Time integration.
    parser.add_argument("--dt", type=float, default=0.5)
    parser.add_argument("--t-end", type=float, default=60.0,
                        help="End time [s] (matches the spec's "
                             "--time-span).")
    parser.add_argument("--fes-order", type=int, default=1)
    parser.add_argument("--linear-solver", default="sparsecholesky",
                        choices=["sparsecholesky", "umfpack", "pardiso"])

    parser.add_argument("--surface-label", default="",
                        help="Optional BND filter (default: all).")
    parser.add_argument("--probe-point", default="",
                        help="Optional 'x,y,z' [m] probe point.")
    parser.add_argument("--csv-output", default="",
                        help="Optional CSV of (t, I, P, T_avg, T_max).")

    def run(args):
        probe_point = None
        if args.probe_point:
            probe_point = [float(s) for s in args.probe_point.split(",")]
            if len(probe_point) != 3:
                return {"error": "--probe-point must be 'x,y,z'."}

        if args.I_ref <= 0:
            return {"error": "--I-ref must be positive."}

        return solve_heat_em_table(
            wp_vol=args.wp_vol,
            em_table_path=args.em_table,
            ht_source=args.ht_source,
            ht_sol=args.ht_sol,
            em_vol=args.em_vol,
            ht_order=int(args.ht_order),
            coil_current=args.coil_current,
            coil_current_csv=args.coil_current_csv,
            I_ref=float(args.I_ref),
            coil_step=args.coil_step,
            biot_image_factor=float(args.biot_image_factor),
            material=args.material,
            rho=args.rho, cp=args.cp, k=args.k,
            h_conv=args.h_conv, t_ext=args.t_ext,
            t_initial=args.t_initial, emissivity=args.emissivity,
            surface_label=args.surface_label,
            dt=args.dt, t_end=args.t_end,
            fes_order=args.fes_order,
            linear_solver=args.linear_solver,
            probe_point=probe_point,
            csv_output=args.csv_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
