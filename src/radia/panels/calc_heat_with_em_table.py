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

Production route
----------------
``calc_heat.py --em-table`` (and ``calc_heat_axisym.py --em-table``) is the
production temperature-dependent source: it keeps the EM solution exactly
at the reference temperature, checks the table against the EM run, uses
the enthalpy Newton integrator with temperature-dependent materials, and
reports thermal exposure.  This module remains for the incident
Biot-Savart source (``--ht-source biot``) and for its existing contracts.

Validity of the frozen |H_t| assumption
---------------------------------------
For a current-driven coil around a ferromagnetic workpiece the spatial
|H_t| is NOT frozen: it falls as sigma(T) falls and rises once the surface
passes the Curie band.  The true power stays nearly flat below the Curie
point and then rises, while a frozen-|H_t| table lookup predicts a rising
power below it and a false self-limit above it (validation_test/
induction_heating/results/coupled_curie_cylinder_frozen_ht.json).  Treat
results of this module as valid for non-magnetic workpieces only.

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
from calc_heat import (  # noqa: E402
    SIGMA_SB,
    THERMAL_PRESETS,
    _boundary_role_audit,
    _input_mesh_geometry_audit,
    _locate_probe,
    _probe_value,
    _resolve_boundary_role, _resolve_convection_roles,
    _resolve_material,
    _temperature_extrema,
)


def _log(msg):
    progress("HEAT_EMTAB", msg)


def _load_current_trajectory(scalar_val, csv_path):
    """Return a callable ``I(t)`` [A_peak] from either a scalar or CSV.

    CSV format: 2 cols ``t_s, I_A``.  The run must lie inside the CSV's
    time range; a time outside it is an error, not a held end value.
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
            if not (t_arr[0] - 1e-12 <= t <= t_arr[-1] + 1e-12):
                raise ValueError(
                    f"t={t:.6g} s is outside --coil-current-csv "
                    f"({t_arr[0]:.6g}..{t_arr[-1]:.6g} s)")
            return float(np.interp(t, t_arr, I_arr))

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


def _collect_surface_dofs(wp_mesh, heat_flux_boundary_names, fes_T):
    """Return ``(dof_indices, dof_xyz)`` for the H1 DOFs lying on the
    heating face.  H1 DOFs are vertices, so this is just the BND
    vertex set under the requested label.
    """
    from ngsolve import BND

    vnrs = set()
    for el in wp_mesh.Elements(BND):
        if el.mat not in heat_flux_boundary_names:
            continue
        for v in el.vertices:
            vnrs.add(v.nr)
    vnrs_sorted = sorted(vnrs)
    xyz = np.array([
        [float(c) for c in wp_mesh.vertices[vnr].point]
        for vnr in vnrs_sorted
    ], dtype=float) if vnrs_sorted else np.zeros((0, 3))
    return np.asarray(vnrs_sorted, dtype=np.int64), xyz


def _project_Ht_ref(wp_mesh, dof_xyz, em_vol, ht_sol):
    """|H_t_ref| at each heat-flux DOF from the EM run's verified _Ht.sol.

    The field is read through its sidecar and projected onto the EM heated
    surface; a DOF off that surface is an error.  Returns (values, audit).
    """
    from ngsolve import Mesh, H1, GridFunction
    from radia import ih_thermal

    em_mesh = Mesh(em_vol)
    pair = ih_thermal.verify_field_pair(ht_sol, em_vol, em_mesh, 1,
                                        quantity=ih_thermal.HT_QUANTITY)
    record = ih_thermal.read_field_sidecar(ht_sol)
    gf_J = GridFunction(H1(em_mesh, order=1))
    gf_J.Load(ht_sol)
    field = ih_thermal.SurfaceP1Field.from_gridfunction(
        em_mesh, gf_J, record["boundaries"])
    values, dist = field.evaluate(np.asarray(dof_xyz, float),
                                  max_distance=0.25 * field.h_median,
                                  what="heat-flux DOFs")
    return values, {**pair, "max_transfer_distance_m": float(dist.max()),
                    "frequency_Hz": record.get("frequency_Hz")}


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


def _surface_vertex_normals(wp_mesh, heat_flux_boundary_names, dof_vnrs):
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
        if el.mat not in heat_flux_boundary_names:
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
                         ht_source, ht_sol, em_vol,
                         coil_current, coil_current_csv,
                         I_ref,
                         coil_step="", biot_image_factor=1.0,
                         material="steel", rho=None, cp=None, k=None,
                         h_conv=10.0, t_ext=20.0, t_initial=20.0, emissivity=0.0,
                         heat_flux_boundaries="",
                         convection_boundaries="",
                         radiation_boundaries="",
                         dt=0.5, t_end=60.0,
                         fes_order=1,
                         linear_solver="sparsecholesky",
                         probe_point=None,
                         csv_output="",
                         allow_table_extrapolation=False,
                         allow_frozen_ht=False, convection_map=None):
    """Backward-Euler heat solve with per-step 2D-table q_surf lookup."""
    setup_paths()
    t0 = time.perf_counter()

    from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                         Integrate, CF, ds, dx, BND, TaskManager,
                         InnerProduct, grad)

    from radia import ih_thermal
    if allow_frozen_ht is not True:
        return {"error": ih_thermal.FROZEN_HT_REFUSAL}
    if not os.path.isfile(wp_vol):
        return {"error": f"--wp-vol not found: {wp_vol}"}
    if not os.path.isfile(em_table_path):
        return {"error": f"--em-table not found: {em_table_path}"}

    tab = load_em_table(em_table_path)
    _log(f"TABLE:{os.path.basename(em_table_path)} "
         f"material={tab.material} f={tab.frequency:g}Hz "
         f"H_range=[{tab.H_grid[0]:.2e},{tab.H_grid[-1]:.2e}] A/m "
         f"T_range=[{tab.T_grid[0]:.1f},{tab.T_grid[-1]:.1f}] C")

    wp_mesh = Mesh(wp_vol)
    if wp_mesh.dim != 3 or wp_mesh.ne == 0:
        return {"error": f"--wp-vol {os.path.basename(wp_vol)} must be a 3D "
                         "volume mesh of the workpiece solid"}
    if len(set(wp_mesh.GetMaterials())) > 1:
        return {"error": "thermal analysis targets the WORKPIECE ONLY, but "
                         f"--wp-vol has materials "
                         f"{sorted(set(wp_mesh.GetMaterials()))}"}
    try:
        mesh_geometry = _input_mesh_geometry_audit(wp_mesh, fes_order)
    except ValueError as exc:
        return {"error": str(exc)}
    _log(f"MESH:loaded {os.path.basename(wp_vol)} "
         f"materials={list(wp_mesh.GetMaterials())} "
         f"boundaries={list(wp_mesh.GetBoundaries())} "
         f"geometry_order={mesh_geometry['input_curve_order']} "
         f"field_order={fes_order} (input geometry preserved)")

    try:
        heat_flux_selector, heat_flux_names = _resolve_boundary_role(
            wp_mesh, heat_flux_boundaries, "--heat-flux-boundaries",
            required=True)
        convection_selector, convection_names, convection_terms, resolved_convection_map = _resolve_convection_roles(
            wp_mesh, convection_map, convection_boundaries, h_conv, t_ext)
        radiation_selector, radiation_names = _resolve_boundary_role(
            wp_mesh, radiation_boundaries, "--radiation-boundaries",
            required=float(emissivity) != 0.0)
    except ValueError as exc:
        return {"error": str(exc)}
    _log(f"BND:heat_flux={heat_flux_names} "
         f"convection={convection_names} radiation={radiation_names}")

    rho_v, cp_v, k_v = _resolve_material(material, rho, cp, k)
    _log(f"MATERIAL:{material} rho={rho_v} cp={cp_v} k={k_v}")

    fes_T = H1(wp_mesh, order=int(fes_order))
    u, v = fes_T.TnT()
    gfT = GridFunction(fes_T)
    # Set() interpolates the constant; a constant coefficient vector is not
    # a constant field for hierarchical H1 order >= 2.
    gfT.Set(CF(float(t_initial)))
    _log(f"FES:H1 order={fes_order} ndof={fes_T.ndof}")

    # -------------------- |H_t_ref(r)| --------------------
    dof_vnrs, dof_xyz = _collect_surface_dofs(
        wp_mesh, set(heat_flux_names), fes_T)
    n_surf = len(dof_vnrs)
    _log(f"BND_DOF:{n_surf} heat-flux surface DOFs on "
         f"{heat_flux_names}")

    ht_audit = None
    if ht_source == "kelvin":
        if not ht_sol or not em_vol:
            return {"error":
                    "--ht-source kelvin requires both --ht-sol "
                    "(<stem>_Ht.sol from calc_fem_kelvin, with its sidecar) "
                    "and --em-vol (<stem>_fem.vol the .sol was saved on)."}
        if not os.path.isfile(ht_sol):
            return {"error": f"--ht-sol not found: {ht_sol}"}
        if not os.path.isfile(em_vol):
            return {"error": f"--em-vol not found: {em_vol}"}
        try:
            Ht_ref_arr, ht_audit = _project_Ht_ref(wp_mesh, dof_xyz, em_vol,
                                                   ht_sol)
        except ValueError as exc:
            return {"error": str(exc)}
        f_ht = ht_audit.get("frequency_Hz")
        if f_ht is None or abs(f_ht - tab.frequency) > 1e-6 * f_ht:
            return {"error": f"the EM table is for {tab.frequency:g} Hz but "
                             f"the |H_t| field records {f_ht!r} Hz"}
        _log(f"HT_REF:projected {n_surf} surface DOFs from "
             f"{os.path.basename(ht_sol)} (max distance "
             f"{ht_audit['max_transfer_distance_m']:.3e} m)")
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
            wp_mesh, set(heat_flux_names), dof_vnrs)
        n_proj = int((np.linalg.norm(normals, axis=1) > 1e-9).sum())
        Ht_ref_arr = _biot_Ht_on_surface(
            segs, dof_xyz, normals, current=I_ref,
            image_factor=biot_image_factor)
        _log(f"HT_REF:biot incident |H_t| at {n_surf} surface DOFs "
             f"({n_proj} tangent-projected) I_ref={I_ref:.3f} A "
             f"image_factor={biot_image_factor:g}; |H_t| range "
             f"[{float(np.min(Ht_ref_arr)):.3e},"
             f"{float(np.max(Ht_ref_arr)):.3e}] A/m")
        if not np.all(np.isfinite(Ht_ref_arr)):
            return {"error": "the incident |H_t| is not finite at "
                             f"{int((~np.isfinite(Ht_ref_arr)).sum())} DOFs "
                             "(degenerate coil or surface geometry)"}
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
    for term in convection_terms:
        a_form += term["h_W_m2K"] * v * u * ds(term["selector"])
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
    heat_flux_region = wp_mesh.Boundaries(heat_flux_selector)
    A_surf = float(_boundary_role_audit(
        wp_mesh, heat_flux_names)["area_m2"])

    # -------------------- Time loop --------------------
    n_steps = int(math.ceil(float(t_end) / float(dt)))
    t_arr = [0.0]
    T_probe_hist = []
    P_total_hist = []
    T_avg_hist = []
    T_max_hist = []

    probe_mip = None
    if probe_point is not None:
        try:
            probe_mip = _locate_probe(wp_mesh, probe_point)
        except ValueError as exc:
            return {"error": str(exc)}
        T_probe_hist.append(_probe_value(gfT, probe_mip))
    from ngsolve import NodeId, VERTEX
    vertex_dofs_T = np.asarray([fes_T.GetDofNrs(NodeId(VERTEX, vx.nr))[0]
                                for vx in wp_mesh.vertices])
    volume = float(Integrate(CF(1.0), wp_mesh).real)
    lo_T, hi_T = float(tab.T_grid[0]), float(tab.T_grid[-1])
    max_T_excursion = 0.0
    max_H_excursion = 1.0

    Q_input_J = 0.0
    linear_true_residual_max = 0.0
    from radia.ih_heat_transient import checked_heat_increment
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
        over = max(0.0, float(T_dof.max()) - hi_T, lo_T - float(T_dof.min()))
        if over > 0.0:
            if not allow_table_extrapolation:
                return {"error":
                        f"t={t:.4g} s: surface temperature "
                        f"{T_dof.min():.1f}..{T_dof.max():.1f} C leaves the "
                        f"EM table {lo_T:.1f}..{hi_T:.1f} C; extend the table "
                        "or pass --allow-em-table-extrapolation"}
            max_T_excursion = max(max_T_excursion, over)
        live = Ht_now > 0
        H_lo = float(np.min(Ht_now[live])) if np.any(live) else 0.0
        H_hi = float(np.max(Ht_now))
        if (H_hi > float(tab.H_grid[-1]) or
                (np.any(live) and H_lo < float(tab.H_grid[0]))):
            if not allow_table_extrapolation:
                return {"error": f"t={t:.4g} s: |H_t| {H_lo:.3e}..{H_hi:.3e} "
                                 "A/m leaves the EM table "
                                 f"{tab.H_grid[0]:.3e}..{tab.H_grid[-1]:.3e}"}
            max_H_excursion = max(max_H_excursion,
                                  H_hi / float(tab.H_grid[-1]),
                                  float(tab.H_grid[0]) / max(H_lo, 1e-300))
        q_dof = interp_qsurf(tab, Ht_now, T_dof)

        qv_fv[:] = 0.0
        for k_dof, idx in enumerate(dof_vnrs):
            qv_fv[int(idx)] = float(q_dof[k_dof])

        f_form = LinearForm(fes_T)
        f_form += gf_q * v * ds(heat_flux_selector)
        for term in convection_terms:
            f_form += term["h_W_m2K"] * term["ambient_C"] * v \
                * ds(term["selector"])
        if float(emissivity) > 0.0:        # radiation (explicit, prev-step T, in K)
            _TK = gfT + 273.15
            f_form += -float(emissivity) * SIGMA_SB \
                * (_TK**4 - (float(t_ext) + 273.15)**4) * v \
                * ds(radiation_selector)
        f_form.Assemble()
        with TaskManager():
            res_vec.data = f_form.vec - a_form.mat * gfT.vec
            delta, relative = checked_heat_increment(
                mstar, inv, res_vec, fes_T.FreeDofs())
            linear_true_residual_max = max(linear_true_residual_max, relative)
            gfT.vec.data += float(dt) * delta

        q_int = float(Integrate(gf_q, wp_mesh, BND,
                                 definedon=heat_flux_region).real)
        Q_input_J += q_int * float(dt)

        T_vert = np.asarray(gfT.vec.FV().NumPy())[vertex_dofs_T]
        T_max_now = float(np.max(T_vert))
        T_avg_now = float(Integrate(gfT, wp_mesh).real) / volume

        t_arr.append(t)
        T_avg_hist.append(T_avg_now)
        T_max_hist.append(T_max_now)
        P_total_hist.append(q_int)

        if probe_mip is not None:
            T_probe_hist.append(_probe_value(gfT, probe_mip))
        _log(f"STEP:{step}/{n_steps} t={t:.3f}s I={I_now:.2f}A "
             f"P_in={q_int:.3e}W T_avg={T_avg_now:.2f}C "
             f"T_max={T_max_now:.2f}C")

    T_min_final, T_max_final, T_extrema = _temperature_extrema(
        gfT, wp_mesh, fes_order)
    T_avg_final = float(Integrate(gfT, wp_mesh).real) / volume
    heat_flux_audit = _boundary_role_audit(
        wp_mesh, heat_flux_names, q_cf=gf_q)
    convection_audit = _boundary_role_audit(wp_mesh, convection_names)
    radiation_audit = _boundary_role_audit(wp_mesh, radiation_names)
    t_total = time.perf_counter() - t0
    _log(f"DONE:T_max={T_max_final:.2f}C "
         f"T_avg={T_avg_final:.2f}C "
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
            raise RuntimeError(f"could not write {csv_output}: "
                               f"{type(e).__name__}: {e}") from e

    return {
        "T_max_C": T_max_final,
        "T_min_C": T_min_final,
        "T_extrema": T_extrema,
        "T_avg_C": T_avg_final,
        "table_T_excursion_C": max_T_excursion,
        "table_H_excursion_ratio": max_H_excursion,
        "ht_projection": ht_audit,
        "T_initial_C": float(t_initial),
        "T_probe_history_C": T_probe_hist if probe_point is not None else None,
        "t_history_s": t_arr,
        "P_total_history_W": [0.0] + P_total_hist,
        "T_avg_history_C": [float(t_initial)] + T_avg_hist,
        "T_max_history_C": [float(t_initial)] + T_max_hist,
        "Q_input_J": Q_input_J,
        "linear_true_residual_max": linear_true_residual_max,
        "surface_area_m2": A_surf,
        "heat_flux_boundaries": heat_flux_selector,
        "convection_boundaries": convection_selector,
        "convection_by_boundary": convection_terms,
        "radiation_boundaries": radiation_selector,
        "boundary_audit": {
            "heat_flux": heat_flux_audit,
            "convection": convection_audit,
            "radiation": radiation_audit,
        },
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
        "mesh_geometry": mesh_geometry,
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
                        help="kelvin path: <stem>_Ht.sol from "
                             "calc_fem_kelvin.py (P1 |H_t| with its "
                             "sidecar).")
    parser.add_argument("--em-vol", default="",
                        help="kelvin path: <stem>_fem.vol the --ht-sol "
                             "was saved on.")
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
                        choices=["sparsecholesky"])

    parser.add_argument("--heat-flux-boundaries", default="",
                        help="Required boundary expression receiving q_surf.")
    parser.add_argument("--convection-map", default="",
                        help="JSON map of exact boundary labels to h_W_m2K and ambient_C; "
                             "unlisted faces have no convection.")
    parser.add_argument("--convection-boundaries", default="",
                        help="Boundary expression receiving convection; "
                             "required when --h-conv is nonzero.")
    parser.add_argument("--radiation-boundaries", default="",
                        help="Boundary expression receiving radiation; "
                             "required when --emissivity is nonzero.")
    parser.add_argument("--surface-label", default=None,
                        help=argparse.SUPPRESS)
    parser.add_argument("--probe-point", default="",
                        help="Optional 'x,y,z' [m] probe point.")
    parser.add_argument("--allow-em-table-extrapolation", action="store_true",
                        help="Clamp surface temperatures / |H_t| outside "
                             "the EM table instead of failing; excursions "
                             "are reported.")
    parser.add_argument("--allow-frozen-ht", action="store_true",
                        help="Acknowledge that this solver freezes the EM "
                             "run's |H_t| (invalid for a ferromagnetic part "
                             "crossing its Curie band).")
    parser.add_argument("--csv-output", default="",
                        help="Optional CSV of (t, I, P, T_avg, T_max).")

    def run(args):
        if args.surface_label is not None:
            return {"error":
                    "--surface-label was removed because it coupled heat "
                    "input, convection, and radiation. Use the separate "
                    "--heat-flux-boundaries, --convection-boundaries, and "
                    "--radiation-boundaries options."}
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
            coil_current=args.coil_current,
            coil_current_csv=args.coil_current_csv,
            I_ref=float(args.I_ref),
            coil_step=args.coil_step,
            biot_image_factor=float(args.biot_image_factor),
            material=args.material,
            rho=args.rho, cp=args.cp, k=args.k,
            h_conv=args.h_conv, t_ext=args.t_ext,
            t_initial=args.t_initial, emissivity=args.emissivity,
            heat_flux_boundaries=args.heat_flux_boundaries,
            convection_boundaries=args.convection_boundaries,
            convection_map=args.convection_map,
            radiation_boundaries=args.radiation_boundaries,
            dt=args.dt, t_end=args.t_end,
            fes_order=args.fes_order,
            linear_solver=args.linear_solver,
            probe_point=probe_point,
            csv_output=args.csv_output,
            allow_table_extrapolation=args.allow_em_table_extrapolation,
            allow_frozen_ht=args.allow_frozen_ht,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
