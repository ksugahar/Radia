"""
Transient heat-transfer solver for IH workpieces (Phase B of the
Radia-NGSolve thermal pipeline).

Inputs
------
``--wp-vol``          Workpiece volume mesh (.vol).  MUST be a
                      WORKPIECE-ONLY 3D volume mesh: a single solid with
                      exactly one volume material region and a
                      explicitly named heat-flux, convection, and
                      radiation boundary roles.
                      The coil / air / Kelvin regions belong to the EM
                      mesh, NOT this thermal mesh.  In Kubota's flow this
                      mesh is SEPARATE from the EM mesh: the EM .vol
                      carries the workpiece as a hole with a SIBC face
                      (WP-HOLE policy), while this thermal mesh carries
                      the workpiece as a real solid with the same outer
                      surface.  A multi-material (coil+wp) or surface-only
                      mesh is rejected with a clear error (the thermal
                      step targets the workpiece solid only).

``--qsurf-sol``       q_surf [W/m^2] saved by ``calc_fem_kelvin.py``
                      (Phase A) on the EM mesh.  Spatial distribution
                      is preserved by per-vertex sampling onto this
                      thermal mesh's surface.

``--em-vol``          The EM .vol the qsurf-sol was computed on (required).
                      When the .sol has a ``.sol.json`` sidecar the mesh
                      digest, order and EM power are verified against it.

``--q-uniform``       Alternative to qsurf-sol: a UNIFORM scalar
                      heat flux [W/m^2].  Useful for testing and
                      rough optimization ladders ("does this much
                      power even get to soak temperature?").

Physics
-------
Heat equation with independently selected surface heat-flux and Newton
convection boundary roles:

    rho * cp * dT/dt = div(k grad T)                    (volume)
     k dT/dn = q_surf(x,y,z)                              (heating face)
     k dT/dn = -h_conv (T - T_ext)                        (cooling face)

The roles may intentionally overlap, but they are never coupled implicitly.

Discretization: backward Euler in time, H1 in space.

    (M + dt*K) T_{n+1} = M T_n + dt*(q_surf*v + h*T_ext*v)_ds

Output
------
JSON summary on stdout (last line) with peak temperature, probe
history and timings. The final spatial temperature and heat-flux fields are
written to the requested GMSH ``.msh v4.1`` artifact.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import time

import numpy as np

# Shared utilities
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import setup_paths, progress, calc_main  # noqa: E402


QSURF_HANDOFF_ORDER = 1


def _log(msg):
    progress("HEAT", msg)


def _input_mesh_geometry_audit(mesh, fes_order):
    """Record the geometry that the thermal solve actually consumes.

    The polynomial order of ``H1`` is independent of the geometry order stored
    in a Netgen ``.vol``.  In particular, calling ``mesh.Curve(fes_order)``
    after loading a ``.vol`` is not a valid way to raise the field order: the
    generating CAD association may be absent or incomplete, and NGSolve can
    silently replace valid imported geometry with a corrupt mapping.  Curving
    must therefore happen in the mesher, before the ``.vol`` is written.

    This helper is deliberately read-only.  The returned measures make that
    contract visible in ``result.json`` and give downstream audits a stable
    record of the input geometry used by the solve.
    """
    from ngsolve import BND, CF, Integrate

    dimension = int(mesh.dim)
    domain_measure = float(Integrate(CF(1), mesh).real)
    boundary_measure = float(Integrate(CF(1), mesh, BND).real)
    if not math.isfinite(domain_measure) or domain_measure <= 0.0:
        raise ValueError(
            "thermal input mesh has a non-positive or non-finite domain "
            f"measure ({domain_measure!r})"
        )
    if not math.isfinite(boundary_measure) or boundary_measure <= 0.0:
        raise ValueError(
            "thermal input mesh has a non-positive or non-finite boundary "
            f"measure ({boundary_measure!r})"
        )

    return {
        "policy": "preserve-input-vol-geometry",
        "post_load_curve_applied": False,
        "field_order": int(fes_order),
        "input_curve_order": int(mesh.GetCurveOrder()),
        "dimension": dimension,
        "domain_measure": domain_measure,
        "domain_measure_unit": "m^3" if dimension == 3 else "m^2",
        "boundary_measure": boundary_measure,
        "boundary_measure_unit": "m^2" if dimension == 3 else "m",
    }


def _temperature_extrema(gf_temperature, mesh, fes_order):
    """Return deterministic physical-field extrema and sampling metadata.

    H1 coefficients above order one are hierarchical coefficients, not
    point values.  Vertex values are sufficient for an order-one H1 field,
    while higher-order fields are also evaluated at volume and boundary
    integration points so edge, face, and cell modes contribute to the
    reported range.
    """
    from ngsolve import BND, IntegrationRule, NodeId, VERTEX, VOL

    vertex_dofs = [
        dof
        for vertex in mesh.vertices
        for dof in gf_temperature.space.GetDofNrs(NodeId(VERTEX, vertex.nr))
        if dof >= 0
    ]
    if not vertex_dofs:
        raise RuntimeError("the thermal H1 space has no vertex DOFs")

    coefficient_values = np.asarray(gf_temperature.vec.FV().NumPy())
    samples = [np.asarray(coefficient_values[vertex_dofs], dtype=float)]
    integration_order = 0

    if int(fes_order) >= 2:
        integration_order = max(6, 2 * int(fes_order) + 2)
        for vb in (VOL, BND):
            for element in mesh.Elements(vb):
                rule = IntegrationRule(element.type, integration_order)
                mapped_rule = mesh.GetTrafo(element)(rule)
                values = np.asarray(gf_temperature(mapped_rule)).reshape(-1)
                if np.iscomplexobj(values):
                    scale = max(1.0, float(np.max(np.abs(values))))
                    if float(np.max(np.abs(values.imag))) > 1.0e-12 * scale:
                        raise RuntimeError(
                            "the thermal GridFunction has non-real sample values"
                        )
                    values = values.real
                samples.append(np.asarray(values, dtype=float))

    values = np.concatenate(samples)
    if not np.all(np.isfinite(values)):
        raise RuntimeError("the thermal GridFunction contains non-finite values")
    return (
        float(np.min(values)),
        float(np.max(values)),
        {
            "method": (
                "vertices"
                if integration_order == 0
                else "vertices-and-volume-boundary-integration-points"
            ),
            "integration_order": int(integration_order),
            "sample_count": int(values.size),
        },
    )


# -----------------------------------------------------------------
# Material presets (workpiece thermal properties)
# -----------------------------------------------------------------
# rho [kg/m^3], cp [J/(kg.K)], k [W/(m.K)] at room temperature.
THERMAL_PRESETS = {
    "steel":    {"rho": 7800,  "cp": 467,  "k": 46.6},
    "aluminum": {"rho": 2700,  "cp": 900,  "k": 237.0},
    "copper":   {"rho": 8960,  "cp": 385,  "k": 401.0},
    "stainless":{"rho": 8000,  "cp": 500,  "k": 16.0},
    "brass":    {"rho": 8530,  "cp": 380,  "k": 109.0},
}


def _resolve_material(name, rho, cp, k):
    """Return (rho, cp, k) for the named preset, with CLI overrides."""
    if name in THERMAL_PRESETS:
        p = THERMAL_PRESETS[name]
        return (rho if rho is not None else p["rho"],
                cp  if cp  is not None else p["cp"],
                k   if k   is not None else p["k"])
    if name == "custom":
        if None in (rho, cp, k):
            raise ValueError(
                "--material custom requires explicit "
                "--rho, --cp and --k.")
        return rho, cp, k
    raise ValueError(
        f"Unknown thermal material preset {name!r}. "
        f"Choose from {list(THERMAL_PRESETS) + ['custom']}.")


def _resolve_boundary_role(mesh, selector, option_name, *, required):
    """Validate one thermal boundary-role selector against ``mesh``.

    NGSolve boundary selectors are regular expressions.  Return the selector
    unchanged for ``ds``/``Boundaries`` plus the concrete boundary names it
    matches, so result artifacts can state exactly where each physical term
    was applied.  Empty active roles and selectors matching no boundary fail
    before assembly.
    """
    text = str(selector or "").strip()
    if not text:
        if required:
            raise ValueError(
                f"{option_name} is required for this active thermal "
                "boundary condition. Pass an explicit mesh boundary name "
                "or NGSolve boundary expression such as 'side|top'.")
        return "", []

    available = sorted(set(mesh.GetBoundaries()))
    try:
        pattern = re.compile(text)
    except re.error as exc:
        raise ValueError(
            f"{option_name}={text!r} is not a valid boundary expression: "
            f"{exc}") from exc
    matched = [name for name in available if pattern.fullmatch(name)]
    if not matched:
        raise ValueError(
            f"{option_name}={text!r} matches no mesh boundary. "
            f"Available boundaries: {available}")
    return text, matched


def _boundary_role_audit(mesh, boundary_names, *, q_cf=None, weight=None):
    """Return deterministic per-boundary area and optional heat input."""
    from ngsolve import BND, CF, Integrate

    measure = CF(1.0) if weight is None else weight
    rows = []
    for name in boundary_names:
        region = mesh.Boundaries(name)
        area = float(Integrate(
            measure, mesh, BND, definedon=region).real)
        row = {"boundary": name, "area_m2": area}
        if q_cf is not None:
            power = float(Integrate(
                q_cf * measure, mesh, BND, definedon=region).real)
            row["heat_input_W"] = power
            row["q_mean_W_m2"] = power / area if area > 0.0 else 0.0
        rows.append(row)

    result = {
        "matched_boundaries": list(boundary_names),
        "area_m2": sum(row["area_m2"] for row in rows),
        "by_boundary": rows,
    }
    if q_cf is not None:
        result["heat_input_W"] = sum(
            row["heat_input_W"] for row in rows)
    return result


def _validate_qsurf_transfer_order(order):
    """Return the supported cross-mesh q_surf order or fail loudly."""
    value = int(order)
    if value != QSURF_HANDOFF_ORDER:
        raise ValueError(
            "Cross-mesh q_surf transfer currently supports "
            "--qsurf-order 1 only. Higher-order H1 coefficients are "
            "hierarchical and cannot be reconstructed by assigning surface "
            "vertex samples; refusing a silently distorted heat source.")
    return value


# -----------------------------------------------------------------
# q_surf source: spatial (.sol) or uniform scalar
# -----------------------------------------------------------------


def qsurf_args(**given):
    """The complete argument record of :func:`_build_qsurf_source`.

    Every field has the command-line default; unknown names are errors, so a
    caller cannot silently pass an option the builder does not read.
    """
    from types import SimpleNamespace
    fields = dict(
        q_uniform=None, qsurf_sol="", em_vol="", qsurf_order=1,
        heat_flux_boundary_names=(), rotation_axis="z", q_phi_average=False,
        q_phi_bin=None, em_heat_boundaries="", q_scale=1.0,
        transfer_tolerance=None, power_tolerance=None, em_table="",
        ht_sol="", ht_scale=1.0, em_reference_temperature=None,
        em_table_azimuths=64, allow_em_table_extrapolation=False,
        em_table_tolerance=0.05, allow_frozen_ht=False, rotor_states="")
    unknown = set(given) - set(fields)
    if unknown:
        raise TypeError(f"unknown q_surf source options {sorted(unknown)}")
    fields.update(given)
    return SimpleNamespace(**fields)


_SPATIAL_OPTIONS = (
    ("qsurf_sol", "--qsurf-sol"), ("em_vol", "--em-vol"),
    ("em_heat_boundaries", "--em-heat-boundaries"),
    ("transfer_tolerance", "--transfer-tolerance"),
    ("power_tolerance", "--power-tolerance"),
    ("q_phi_average", "--q-phi-average"), ("q_phi_bin", "--q-phi-bin"),
    ("em_table", "--em-table"), ("ht_sol", "--ht-sol"),
    ("allow_frozen_ht", "--allow-frozen-ht"),
    ("em_reference_temperature", "--em-reference-temperature"),
    ("rotor_states", "--rotor-states"))


def _reject_options_with_uniform(args):
    """--q-uniform has no EM source; options that shape one are errors."""
    given = [flag for attr, flag in _SPATIAL_OPTIONS
             if getattr(args, attr, None) not in (None, "", False)]
    for attr, flag in (("q_scale", "--q-scale"), ("ht_scale", "--ht-scale")):
        if float(getattr(args, attr, 1.0)) != 1.0:
            given.append(flag)
    if given:
        raise ValueError(
            f"--q-uniform is a constant flux; {', '.join(given)} only apply "
            "to a spatial --qsurf-sol source")


def _build_qsurf_source(wp_mesh, args):
    """Return ``(q_cf, resample_fn, audit)`` for the heating face.

    Modes:

    1. ``--q-uniform``: a constant flux; ``resample_fn`` is ``None``.  Any
       option that shapes a spatial source is an error.
    2. ``--qsurf-sol`` + ``--em-vol`` with ``--q-phi-average``: the exact
       circumferential average of the EM source; static.
    3. ``--qsurf-sol`` + ``--em-vol``: the EM source evaluated at the thermal
       heat-flux vertices by projection onto the EM heated surface.
    4. ``--rotor-states``: one EM solution per rotor angle of a part that is
       not a body of revolution; the source is written at the first state.

    ``resample_fn`` is always ``None``; a rotating run builds an
    :class:`ih_thermal.RotatingSurfaceSource` from ``audit["_source"]`` or
    ``audit["_rotor_states"]``.

    Every transfer goes through :meth:`ih_thermal.EMHeatSource.transfer`,
    which refuses vertices off the EM surface and a heat input that differs
    from the EM power by more than ``--power-tolerance``.
    """
    from ngsolve import H1, GridFunction, CoefficientFunction
    from radia import ih_thermal

    if args.q_uniform is not None:
        _reject_options_with_uniform(args)
        _log(f"Q_SURF:uniform {args.q_uniform:.4e} W/m^2")
        return (CoefficientFunction(float(args.q_uniform)), None,
                {"mode": "uniform"})

    if args.rotor_states:
        for attr, flag in (("qsurf_sol", "--qsurf-sol"), ("em_vol", "--em-vol"),
                           ("q_phi_average", "--q-phi-average"),
                           ("em_table", "--em-table"), ("ht_sol", "--ht-sol")):
            if getattr(args, attr):
                raise ValueError(f"--rotor-states lists its own EM solutions; "
                                 f"{flag} does not apply with it")
        _validate_qsurf_transfer_order(args.qsurf_order)
        states, meta = ih_thermal.load_rotor_states(
            args.rotor_states, em_boundaries=args.em_heat_boundaries or None,
            q_scale=float(args.q_scale),
            transfer_tolerance=args.transfer_tolerance)
        if meta["axis"] != str(args.rotation_axis).lower().strip():
            raise ValueError(f"--rotor-states turns about {meta['axis']}, "
                             f"--rotation-axis is {args.rotation_axis}")
        power_tol = (ih_thermal.DEFAULT_POWER_TOLERANCE
                     if args.power_tolerance is None else args.power_tolerance)
        gf_wp_q = GridFunction(H1(wp_mesh, order=1))
        a0, s0 = states[0]
        first = s0.transfer(
            wp_mesh, list(args.heat_flux_boundary_names), gf_wp_q,
            power_tolerance=power_tol,
            theta=a0 if meta["frame"] == "world" else 0.0)
        first.pop("vertex_numbers")
        audit = {"mode": "rotor-states", **meta,
                 "power_tolerance": float(power_tol), "initial": first,
                 "states": [{"angle_rad": a, **src.audit()}
                            for a, src in states],
                 "_rotor_states": (states, meta)}
        _log(f"Q_SURF:rotor states {meta['n_states']} from "
             f"{os.path.basename(args.rotor_states)} (frame {meta['frame']}, "
             f"period {meta['period_rad']:.6f} rad)")
        return gf_wp_q, None, audit
    if not args.qsurf_sol:
        raise ValueError(
            "Either --q-uniform, --qsurf-sol or --rotor-states is required "
            "(no heat input was supplied).")
    if not args.em_vol:
        raise ValueError(
            "--em-vol is required when --qsurf-sol is supplied.  "
            "NGSolve .sol files do not contain mesh information, "
            "so the EM .vol that the .sol was saved against must "
            "be passed explicitly.")
    _validate_qsurf_transfer_order(args.qsurf_order)

    axis_str = str(args.rotation_axis).lower().strip()
    if axis_str not in ("x", "y", "z"):
        raise ValueError(
            f"--rotation-axis must be one of x / y / z (got {axis_str!r}).")

    heat_names = list(args.heat_flux_boundary_names)
    phi_average = bool(args.q_phi_average)
    source = ih_thermal.EMHeatSource(
        args.qsurf_sol, args.em_vol,
        em_boundaries=args.em_heat_boundaries or None,
        mode="phi-average" if phi_average else "direct",
        axis=axis_str, q_scale=float(args.q_scale),
        transfer_tolerance=args.transfer_tolerance,
        bin_width=args.q_phi_bin)
    power_tol = (ih_thermal.DEFAULT_POWER_TOLERANCE
                 if args.power_tolerance is None else args.power_tolerance)
    _log(f"Q_SURF:{os.path.basename(source.qsurf_sol)} on "
         f"{os.path.basename(source.em_vol)}, EM boundaries "
         f"{source.em_boundaries} ({source.boundary_rule}), "
         f"P_EM={source.power:.6e} W")

    gf_wp_q = GridFunction(H1(wp_mesh, order=1))
    first = source.transfer(wp_mesh, heat_names, gf_wp_q,
                            power_tolerance=power_tol)
    vnrs = first.pop("vertex_numbers")
    audit = source.audit()
    audit.update({"power_tolerance": float(power_tol), "initial": first})
    _log(f"Q_SURF:{source.mode} transfer to {first['target_vertices']} "
         f"vertices, max distance {first['max_transfer_distance_m']:.3e} m, "
         f"P_thermal/P_EM - 1 = "
         f"{first['power_balance']['relative_error']:+.3e}")

    if args.em_table:
        if args.em_reference_temperature is None:
            raise ValueError("--em-reference-temperature (the workpiece "
                             "temperature of the EM run) is required with "
                             "--em-table")
        q_ref = np.asarray(gf_wp_q.vec.FV().NumPy())[vnrs].copy()
        xyz = ih_thermal.mesh_vertices(wp_mesh)[vnrs]
        coupled, coupled_audit = ih_thermal.build_temperature_dependent_source(
            source, table_path=args.em_table, target_xyz=xyz,
            q_ref=q_ref, gf_target=gf_wp_q, dofs=vnrs,
            T_ref=float(args.em_reference_temperature),
            H_scale=float(args.ht_scale), ht_sol=args.ht_sol,
            azimuths=int(args.em_table_azimuths) if phi_average else 0,
            allow_extrapolation=bool(args.allow_em_table_extrapolation),
            consistency_tolerance=args.em_table_tolerance,
            acknowledge_frozen_ht=bool(args.allow_frozen_ht))
        audit["temperature_dependent_source"] = coupled_audit
        audit["_coupled"] = coupled
        _log(f"Q_SURF:temperature-dependent source from "
             f"{os.path.basename(args.em_table)} (T_ref="
             f"{float(args.em_reference_temperature):g} C, |H_t| "
             f"{coupled_audit['H_t_range_A_m'][0]:.3e}.."
             f"{coupled_audit['H_t_range_A_m'][1]:.3e} A/m)")
        return gf_wp_q, None, audit
    for attr, flag in (("ht_sol", "--ht-sol"),
                       ("allow_frozen_ht", "--allow-frozen-ht")):
        if getattr(args, attr):
            raise ValueError(f"{flag} only applies with --em-table")
    if not phi_average:
        audit["_source"] = source
    return gf_wp_q, None, audit


# -----------------------------------------------------------------
# Time-domain heat solve
# -----------------------------------------------------------------

SIGMA_SB = 5.670374419e-8     # Stefan-Boltzmann constant [W/m^2/K^4]


def solve_heat(wp_vol,
               material="steel", rho=None, cp=None, k=None,
               h_conv=10.0, t_ext=20.0, t_initial=20.0, emissivity=0.0,
               heat_flux_boundaries="",
               convection_boundaries="",
               radiation_boundaries="",
               q_uniform=None, qsurf_sol="", em_vol="",
               qsurf_order=1,
               q_phi_average=False, q_phi_bin=None,
               em_heat_boundaries="", q_scale=1.0,
               transfer_tolerance=None, power_tolerance=None,
               exposure_thresholds=(), temperature_limit=None,
               material_table="", latent_heat=0.0, latent_range=None,
               em_table="", ht_sol="", ht_scale=1.0,
               em_reference_temperature=None, em_table_azimuths=64,
               allow_em_table_extrapolation=False,
               em_table_tolerance=0.05,
               allow_frozen_ht=False,
               allow_table_extrapolation=False,
               newton_tol=1.0e-3, newton_max_iter=25,
               max_halvings=0,
               dt=0.5, t_end=5.0,
               time_scheme="backward-euler",
               linear_solver="sparsecholesky",
               fes_order=1,
               rotation_rpm=0.0,
               rotation_axis="z",
               rotor_states="", angle_step_tolerance=None,
               probe_point=None,
               msh_output="",
               csv_output="",
               _wp_mesh=None,
               _write_solution=True):
    """Run the transient heat solve.  See module docstring for inputs."""
    setup_paths()
    t0 = time.perf_counter()
    from radia import ih_thermal

    # Lazy imports so --help is fast.
    from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                          Integrate, CF, CoefficientFunction, ds, dx, BND,
                          TaskManager)

    if _wp_mesh is None and not os.path.isfile(wp_vol):
        return {"error": f"--wp-vol not found: {wp_vol}"}

    wp_mesh = Mesh(wp_vol) if _wp_mesh is None else _wp_mesh

    # --- Workpiece-only volume mesh contract (radia-ih thermal) -------
    # The thermal step targets the WORKPIECE SOLID only.  The coil /
    # air / Kelvin regions live on the EM mesh (WP-HOLE policy); this
    # thermal mesh must be a single workpiece solid.  Fail loud rather
    # than silently heating the coil (keiko 2026-05-31: a coil+wp .vol
    # was passed to the thermal step and the coil diffused heat as
    # 'steel').  Strict-in-what-we-accept per CLAUDE.md No-Fallback.
    if wp_mesh.dim != 3:
        return {"error":
                f"--wp-vol {os.path.basename(wp_vol)} is {wp_mesh.dim}D; "
                f"calc_heat.py needs a 3D volume mesh of the workpiece "
                f"solid.  If this is a 2D (r,z) axisymmetric workpiece, "
                f"use calc_heat_axisym.py.  If it is a surface/SIBC mesh "
                f"(no solid), that belongs to the EM step -- export the "
                f"workpiece as a real 3D solid for the thermal step."}
    if wp_mesh.ne == 0:
        return {"error":
                f"--wp-vol {os.path.basename(wp_vol)} has 0 volume "
                f"elements (it looks like a surface-only / SIBC mesh).  "
                f"The thermal step needs a VOLUME mesh of the workpiece "
                f"solid; the SIBC-faced hole surface belongs to the EM "
                f"step (calc_fem_kelvin)."}
    _wp_mats = sorted(set(wp_mesh.GetMaterials()))
    if len(_wp_mats) > 1:
        return {"error":
                f"thermal analysis targets the WORKPIECE ONLY, but "
                f"--wp-vol {os.path.basename(wp_vol)} has {len(_wp_mats)} "
                f"volume materials {_wp_mats}.  Export a workpiece-only "
                f"volume mesh (a single solid) for the thermal step -- "
                f"the coil / air / Kelvin regions belong to the EM mesh, "
                f"not the thermal mesh."}

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
        convection_selector, convection_names = _resolve_boundary_role(
            wp_mesh, convection_boundaries, "--convection-boundaries",
            required=float(h_conv) != 0.0)
        radiation_selector, radiation_names = _resolve_boundary_role(
            wp_mesh, radiation_boundaries, "--radiation-boundaries",
            required=float(emissivity) != 0.0)
    except ValueError as exc:
        return {"error": str(exc)}
    _log(f"BND:heat_flux={heat_flux_names} "
         f"convection={convection_names} radiation={radiation_names}")

    rho_v, cp_v, k_v = _resolve_material(material, rho, cp, k)
    _log(f"MATERIAL:{material} rho={rho_v} cp={cp_v} k={k_v}")

    # H1 FES on the workpiece volume.  No Dirichlet BC -- the
    # surface flux + Robin convection give a well-posed problem.
    fes_T = H1(wp_mesh, order=int(fes_order))
    u, v = fes_T.TnT()
    gfT = GridFunction(fes_T)
    with TaskManager():
        # Uniform initial state.  ``gfT.vec[:] = T0`` is WRONG for
        # order >= 2 (hierarchical edge/face coefficients are not
        # nodal temperatures); Set() interpolates the constant exactly.
        gfT.Set(CF(float(t_initial)))
    _log(f"FES:H1 order={fes_order} ndof={fes_T.ndof}")

    # q_surf source -- needs a Namespace-like object to thread the
    # CLI args through helper.
    a_local = qsurf_args(
        q_uniform=q_uniform, qsurf_sol=qsurf_sol, em_vol=em_vol,
        qsurf_order=qsurf_order, heat_flux_boundary_names=heat_flux_names,
        rotation_axis=rotation_axis, q_phi_average=q_phi_average,
        q_phi_bin=q_phi_bin, em_heat_boundaries=em_heat_boundaries,
        q_scale=q_scale, transfer_tolerance=transfer_tolerance,
        power_tolerance=power_tolerance, em_table=em_table, ht_sol=ht_sol,
        ht_scale=ht_scale, em_reference_temperature=em_reference_temperature,
        em_table_azimuths=em_table_azimuths,
        allow_em_table_extrapolation=allow_em_table_extrapolation,
        em_table_tolerance=em_table_tolerance,
        allow_frozen_ht=allow_frozen_ht, rotor_states=rotor_states)
    heat_flux_region = wp_mesh.Boundaries(heat_flux_selector)
    try:
        q_cf, q_resample, qsurf_projection = _build_qsurf_source(
            wp_mesh, a_local)
    except (ValueError, FileNotFoundError) as exc:
        return {"error": str(exc)}
    coupled = qsurf_projection.pop("_coupled", None)
    single = qsurf_projection.pop("_source", None)
    rotor = qsurf_projection.pop("_rotor_states", None)
    if coupled is not None and float(rotation_rpm) > 0.0 and not q_phi_average:
        return {"error": "a temperature-dependent source on a rotating part "
                         "needs --q-phi-average (the fast-rotation limit)"}

    # Rotation: the source is tabulated in the body angle (one EM solution
    # turned with a body of revolution, or one EM solution per rotor angle)
    # and each step applies its exact average over the angles the step
    # sweeps.  Mesh / FES / stiffness / mass are held fixed; only the RHS
    # depends on q_cf.
    omega_mech = (2.0 * math.pi / 60.0) * float(rotation_rpm)
    if not (math.isfinite(omega_mech) and omega_mech >= 0.0):
        return {"error": f"--rotation-rpm must be >= 0, got {rotation_rpm!r}"}
    if rotor is not None and omega_mech == 0.0:
        return {"error": "--rotor-states describes a rotating part; give "
                         "--rotation-rpm > 0"}
    rotating = None
    if omega_mech > 0.0 and (single is not None or rotor is not None):
        try:
            if rotor is not None:
                states, meta = rotor
                rotating = ih_thermal.RotatingSurfaceSource(
                    wp_mesh, heat_flux_names, q_cf, states=states,
                    power_tolerance=qsurf_projection["power_tolerance"],
                    axis=meta["axis"], period=meta["period_rad"],
                    frame=meta["frame"],
                    angle_step_tolerance=angle_step_tolerance)
            else:
                rotating = ih_thermal.RotatingSurfaceSource(
                    wp_mesh, heat_flux_names, q_cf, states=[(0.0, single)],
                    power_tolerance=qsurf_projection["power_tolerance"],
                    axis=rotation_axis,
                    angle_step_tolerance=angle_step_tolerance)
        except ValueError as exc:
            return {"error": str(exc)}
        qsurf_projection["rotation"] = {
            **rotating.audit(), "omega_rad_s": omega_mech,
            "swept_per_step_rad": omega_mech * float(dt),
            "step_source": "exact average over the swept angles"}
        _log(f"ROTATION:rpm={float(rotation_rpm):g} {rotating.kind}, "
             f"{len(rotating.angles)} knots, {omega_mech * float(dt):.3f} rad "
             "per step (step-averaged)")
    elif omega_mech > 0.0:
        _log("ROTATION:rpm>0 has no effect here -- q_surf is azimuthally "
             "uniform (--q-uniform constant, or --q-phi-average already "
             "gives the rotation-averaged axisymmetric q).")
    rotation_active = rotating is not None

    # ------------------- Material model -------------------
    from radia import ih_heat_transient
    from radia import ih_thermal_material
    try:
        if material_table:
            thermal_material = ih_thermal_material.ThermalMaterial.from_csv(
                material_table, rho=rho_v, latent_heat=latent_heat,
                latent_range=latent_range,
                allow_extrapolation=allow_table_extrapolation)
        else:
            if latent_heat:
                raise ValueError("--latent-heat needs --material-table")
            thermal_material = ih_thermal_material.ThermalMaterial.constant(
                rho_v, cp_v, k_v)
    except (ValueError, OSError) as exc:
        return {"error": f"material table: {exc}"}
    nonlinear = (not thermal_material.is_constant) or coupled is not None
    if nonlinear and time_scheme != "backward-euler":
        return {"error": "temperature-dependent materials are integrated "
                         "with the enthalpy backward-Euler scheme only; "
                         "drop --time-scheme crank-nicolson"}
    if nonlinear:
        _log(f"MATERIAL:table {thermal_material.source} "
             f"k {thermal_material.k.min():g}..{thermal_material.k.max():g} "
             f"cp {thermal_material.cp.min():g}..{thermal_material.cp.max():g}"
             f" latent {thermal_material.latent_heat:g} J/kg")

    if not nonlinear:
        # ------------------- Bilinear forms -------------------
        K_cf = CF(float(k_v))
        rho_cp = CF(float(rho_v) * float(cp_v))

        a_form = BilinearForm(fes_T, symmetric=True)
        a_form += K_cf * grad_dot(u, v) * dx
        if float(h_conv) != 0.0:
            a_form += float(h_conv) * v * u * ds(convection_selector)
        a_form.Assemble()

        m_form = BilinearForm(fes_T, symmetric=True)
        m_form += rho_cp * u * v * dx
        m_form.Assemble()

        if time_scheme not in ("backward-euler", "crank-nicolson"):
            raise ValueError(
                f"Unsupported --time-scheme {time_scheme!r} "
                "(expected backward-euler or crank-nicolson).")

        # mstar = M + theta * dt * K
        theta = 1.0 if time_scheme == "backward-euler" else 0.5
        mstar = m_form.mat.CreateMatrix()
        mstar.AsVector().data = (
            m_form.mat.AsVector() + (theta * float(dt)) * a_form.mat.AsVector())
        inv = mstar.Inverse(freedofs=fes_T.FreeDofs(),
                            inverse=linear_solver)
        res_vec = gfT.vec.CreateVector()
        _log(f"SOLVER:{linear_solver} ({time_scheme}, dt={dt}, t_end={t_end})")

    else:
        stepper = ih_heat_transient.NonlinearHeatStepper(
            gfT, thermal_material,
            ih_heat_transient.HeatBoundaryTerms(
                heat_flux=heat_flux_selector,
                convection=convection_selector if float(h_conv) else "",
                h_conv=float(h_conv), t_ext=float(t_ext),
                radiation=radiation_selector if float(emissivity) else "",
                emissivity=float(emissivity)),
            q_source=q_cf,
            q_update=(None if coupled is None else
                      (lambda g: coupled.update(g, coupled_T_dofs))),
            q_damping=None if coupled is None else coupled.gf_damping,
            linear_solver=linear_solver,
            newton_tol_K=float(newton_tol),
            max_newton=int(newton_max_iter),
            max_halvings=int(max_halvings))
        if coupled is not None:
            from ngsolve import NodeId as _NodeId, VERTEX as _VERTEX
            coupled_T_dofs = np.asarray([
                fes_T.GetDofNrs(_NodeId(_VERTEX, int(vn)))[0]
                for vn in coupled.dofs])
            try:
                coupled.update(gfT, coupled_T_dofs)
            except ValueError as exc:
                return {"error": str(exc)}
        _log(f"SOLVER:{linear_solver} (enthalpy backward-euler + Newton, "
             f"dt={dt}, t_end={t_end})")

    # ------------------- Time loop -------------------
    t_arr = [0.0]
    T_probe = []
    probe_mip = None
    if probe_point is not None:
        try:
            probe_mip = _locate_probe(wp_mesh, probe_point)
        except ValueError as exc:
            return {"error": str(exc)}
        T_probe.append(_probe_value(gfT, probe_mip))

    from ngsolve import NodeId, VERTEX
    vertex_dofs_T = np.asarray([
        fes_T.GetDofNrs(NodeId(VERTEX, vtx.nr))[0]
        for vtx in wp_mesh.vertices])
    T_max_history = [float(t_initial)]
    n_steps = int(math.ceil(t_end / dt))
    Q_input_J = 0.0
    heat_flux_audit = _boundary_role_audit(
        wp_mesh, heat_flux_names, q_cf=q_cf)
    convection_audit = _boundary_role_audit(wp_mesh, convection_names)
    radiation_audit = _boundary_role_audit(wp_mesh, radiation_names)
    A_surf = float(heat_flux_audit["area_m2"])
    q_int = float(heat_flux_audit["heat_input_W"])
    _log(f"Q_SURF:int q_surf dA = {q_int:.4e} W "
         f"(area {A_surf:.4e} m^2)")

    for step in range(1, n_steps + 1):
        t = step * float(dt)
        if rotation_active:
            # the source averaged over the angles swept during this step
            rotating.average(omega_mech * (t - float(dt)), omega_mech * t)
        if nonlinear:
            try:
                with TaskManager():
                    stepper.advance(float(dt))
            except (ValueError, RuntimeError) as exc:
                return {"error": f"t={t:.4g} s: {exc}"}
        else:
            f_form = LinearForm(fes_T)
            f_form += q_cf * v * ds(heat_flux_selector)
            if float(h_conv) != 0.0:
                f_form += float(h_conv) * float(t_ext) * v \
                    * ds(convection_selector)
            if float(emissivity) > 0.0:    # radiation (explicit, prev-step T, in K)
                _TK = gfT + 273.15
                f_form += -float(emissivity) * SIGMA_SB \
                    * (_TK**4 - (float(t_ext) + 273.15)**4) * v \
                    * ds(radiation_selector)
            f_form.Assemble()
            with TaskManager():
                res_vec.data = f_form.vec - a_form.mat * gfT.vec
                gfT.vec.data += float(dt) * (inv * res_vec)
        if rotation_active:
            # the step-averaged source changes with the angle
            heat_flux_audit = _boundary_role_audit(
                wp_mesh, heat_flux_names, q_cf=q_cf)
            q_int = float(heat_flux_audit["heat_input_W"])
        Q_input_J += q_int * float(dt)
        t_arr.append(t)
        T_max_history.append(
            float(np.max(gfT.vec.FV().NumPy()[vertex_dofs_T])))
        if probe_mip is not None:
            T_probe.append(_probe_value(gfT, probe_mip))
        _log(f"STEP:{step}/{n_steps} t={t:.3f}s "
             f"T_probe={T_probe[-1] if probe_point is not None else 'n/a'}")

    if nonlinear:
        try:
            stepper.check_energy()
        except RuntimeError as exc:
            return {"error": str(exc)}
        Q_input_J = stepper.audit.energy_in_J
        q_int = stepper.last_heat_input_W
    if coupled is not None:
        qsurf_projection["temperature_dependent_source"].update(
            coupled.audit())

    # Final stats use physical field samples.  Raw order>=2 H1
    # coefficients are not temperatures, and vertices alone can miss an
    # edge-, face-, or cell-mode extremum.
    T_min, T_max, T_extrema = _temperature_extrema(
        gfT, wp_mesh, fes_order
    )
    from radia import ih_thermal_post
    thresholds = sorted(set(float(t) for t in (exposure_thresholds or ())))
    if temperature_limit is not None:
        thresholds = sorted(set(thresholds) | {float(temperature_limit)})
    thermal_exposure = ih_thermal_post.thermal_exposure(
        wp_mesh, gfT, thresholds, axis=rotation_axis)
    limit_record = (ih_thermal_post.limit_check(thermal_exposure,
                                                temperature_limit)
                    if temperature_limit is not None else None)
    for entry in thermal_exposure["thresholds"]:
        _log(f"EXPOSURE:T>{entry['T_C']:g} C volume "
             f"{entry['volume_m3'] * 1e9:.3f} mm^3 "
             f"({entry['volume_fraction']:.3%})")
    if "hottest_ring" in thermal_exposure:
        ring = thermal_exposure["hottest_ring"]
        if "spread_C" in ring:
            _log(f"EXPOSURE:hottest ring r={ring['r_m']:.4g} m spread "
                 f"{ring['spread_C']:.1f} C around the {rotation_axis} axis")
    # Volume-averaged mean temperature -- the integral quantity
    # (int T dV / int dV), the physically meaningful "average" rather
    # than a nodal mean (kubota 2026-05-29: report mean/max/min).
    _vol = float(Integrate(CF(1), wp_mesh))
    T_mean = float(Integrate(gfT, wp_mesh)) / _vol if _vol > 0 else T_max

    # GMSH .msh export of the final temperature field (Open GMSH).
    # Bundled fields:
    #   T_C    -- per-vertex temperature in degC (volume scalar)
    #   q_surf -- the input flux density on the heating face (surface
    #             scalar).  Saved alongside T so the user can confirm
    #             that the projected EM source landed where expected
    #             (q_surf is non-zero only on the SIBC vertices) and
    #             cross-check the integral against P_total visually.
    gmsh_file = ""
    T_sol_file = ""
    heat_vol_file = ""
    if msh_output:
        try:
            from gmsh_post_export import save_vol_sol_pair, vol2msh
            base_dir = os.path.dirname(os.path.abspath(msh_output))
            stem = os.path.splitext(os.path.basename(msh_output))[0]
            sol_T = os.path.join(base_dir, f"{stem}_T.sol").replace("\\", "/")
            vol_T = os.path.join(base_dir, f"{stem}_heat.vol").replace("\\", "/")
            save_vol_sol_pair(vol_T, sol_T, wp_mesh.ngmesh, gfT)
            T_sol_file = sol_T
            ih_thermal.write_field_sidecar(
                sol_T, mesh_path=vol_T, mesh=wp_mesh,
                fes_order=int(fes_order),
                quantity=ih_thermal.TEMPERATURE_QUANTITY,
                unit=ih_thermal.TEMPERATURE_UNIT,
                extra={"producer": "calc_heat", "geometry": "3d",
                       "t_end_s": float(t_end)})
            heat_vol_file = vol_T
            sol_entries = [
                {"sol": sol_T, "fes": "H1",
                 "fes_order": int(fes_order),
                 "fes_dim": 1,
                 "name": "T_C", "ncomp": 1},
            ]
            # Project q_cf onto a workpiece-mesh H1 GridFunction so
            # GMSH renders it.  The CF itself may be either a uniform
            # scalar (--q-uniform mode) or the cross-mesh projection
            # already living in gf_q from _build_qsurf_source; in either
            # case Set on the surface region is the right thing.
            try:
                fes_qg = H1(wp_mesh, order=int(fes_order))
                gf_qg = GridFunction(fes_qg)
                gf_qg.vec[:] = 0
                gf_qg.Set(q_cf,
                           definedon=heat_flux_region)
                sol_q = os.path.join(base_dir,
                                     f"{stem}_qsurf.sol").replace("\\", "/")
                gf_qg.Save(sol_q)
                sol_entries.append(
                    {"sol": sol_q, "fes": "H1",
                     "fes_order": int(fes_order),
                     "fes_dim": 1,
                     "name": "q_surf", "ncomp": 1})
            except Exception as e:
                raise RuntimeError(
                    f"could not write the q_surf overlay: "
                    f"{type(e).__name__}: {e}") from e
            vol2msh(msh_output, vol_T, sol_entries)
            gmsh_file = msh_output
            _log(f"GMSH:wrote {os.path.basename(msh_output)} "
                 f"({len(sol_entries)} fields)")
        except Exception as e:
            raise RuntimeError(
                f"thermal field export to {msh_output} failed: "
                f"{type(e).__name__}: {e}") from e
    elif _write_solution:
        # No --msh-output requested.  Still save the T GridFunction
        # next to the wp .vol so a later evaluation pass (e.g.
        # reload + post-process T at arbitrary points, or feed T
        # back into a second EM solve as a temperature-dependent
        # sigma) has access to it.  Mirrors the qsurf.sol contract
        # on the EM side: the .sol file lives ALONGSIDE the .vol it
        # was solved on, with a fixed naming convention.
        try:
            base_dir = os.path.dirname(os.path.abspath(wp_vol))
            stem = os.path.splitext(os.path.basename(wp_vol))[0]
            sol_T = os.path.join(
                base_dir, f"{stem}_heat_T.sol").replace("\\", "/")
            gfT.Save(sol_T)
            T_sol_file = sol_T
            ih_thermal.write_field_sidecar(
                sol_T, mesh_path=wp_vol, mesh=wp_mesh,
                fes_order=int(fes_order),
                quantity=ih_thermal.TEMPERATURE_QUANTITY,
                unit=ih_thermal.TEMPERATURE_UNIT,
                extra={"producer": "calc_heat", "geometry": "3d",
                       "t_end_s": float(t_end)})
            # heat_vol_file stays "" -- in this branch the wp_vol
            # itself IS the companion mesh (no separate _heat.vol
            # is written because there is no GMSH bundle to anchor).
            _log(f"T_SOL:wrote {os.path.basename(sol_T)} "
                 f"(no GMSH bundle requested; load with the same "
                 f"wp_vol + H1 order={fes_order})")
        except Exception as e:
            raise RuntimeError(
                f"could not save the temperature field: "
                f"{type(e).__name__}: {e}") from e

    # CSV export of probe history.
    if csv_output and probe_point is not None:
        try:
            import csv
            with open(csv_output, "w", newline="", encoding="utf-8-sig") as f:
                w = csv.writer(f)
                w.writerow(["t_s", "T_C"])
                for ti, Ti in zip(t_arr, T_probe):
                    w.writerow([f"{ti:.6f}", f"{Ti:.6f}"])
            _log(f"CSV:wrote {os.path.basename(csv_output)}")
        except Exception as e:
            raise RuntimeError(
                f"could not write {csv_output}: "
                f"{type(e).__name__}: {e}") from e

    t_total = time.perf_counter() - t0
    _log(f"DONE:T_max={T_max:.2f} C  Q_input={Q_input_J:.4e} J "
         f"t={t_total:.1f}s")

    return {
        "T_max_C": T_max,
        "T_min_C": T_min,
        "T_extrema": T_extrema,
        "T_mean_C": T_mean,
        "T_initial_C": float(t_initial),
        "T_probe_history_C": T_probe if probe_point is not None else None,
        "t_history_s": t_arr,
        "Q_input_J": Q_input_J,
        "q_surf_int_W": q_int,
        "surface_area_m2": A_surf,
        "n_steps": n_steps,
        "dt_s": float(dt),
        "t_end_s": float(t_end),
        "time_scheme": time_scheme,
        "linear_solver": linear_solver,
        "ndof": int(fes_T.ndof),
        "ne": int(wp_mesh.ne),
        "fes_order": int(fes_order),
        "mesh_geometry": mesh_geometry,
        "material": material,
        "thermal_material": thermal_material.audit(),
        "nonlinear_transient": stepper.audit.as_dict() if nonlinear else None,
        "rho_kg_m3": float(rho_v),
        "cp_J_kgK": float(cp_v),
        "k_W_mK": float(k_v),
        "rotation_rpm": float(rotation_rpm),
        "h_conv_W_m2K": float(h_conv),
        "t_ext_C": float(t_ext),
        "emissivity": float(emissivity),
        "heat_flux_boundaries": heat_flux_selector,
        "convection_boundaries": convection_selector,
        "radiation_boundaries": radiation_selector,
        "boundary_audit": {
            "heat_flux": heat_flux_audit,
            "convection": convection_audit,
            "radiation": radiation_audit,
        },
        "q_source": ("uniform" if q_uniform is not None
                     else ("qsurf_sol_phi_average" if q_phi_average
                           else "qsurf_sol")),
        "q_phi_average": bool(q_phi_average),
        "qsurf_projection": qsurf_projection,
        "thermal_exposure": thermal_exposure,
        "temperature_limit": limit_record,
        "T_max_history_C": T_max_history,
        "qsurf_sol": qsurf_sol if q_uniform is None else "",
        "em_vol": em_vol if q_uniform is None else "",
        "T_sol_file": T_sol_file,
        "heat_vol_file": heat_vol_file,
        "msh_file": gmsh_file,
        "csv_file": csv_output if (csv_output and probe_point is not None)
                     else "",
        "t_total_s": round(t_total, 2),
    }


def _locate_probe(mesh, point):
    """Return the mesh integration point for ``point`` or raise."""
    try:
        coords = [float(c) for c in point]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"--probe-point must be x,y,z numbers: {point!r}") \
            from exc
    try:
        mip = mesh(*coords)
    except Exception as exc:                       # NGSolve raises NgException
        raise ValueError(f"--probe-point {coords} cannot be located in the "
                         f"thermal mesh: {exc}") from exc
    if mip.nr < 0:
        raise ValueError(f"--probe-point {coords} lies outside the thermal "
                         "mesh")
    return mip


def _probe_value(gf, mip):
    val = gf(mip)
    return float(getattr(val, "real", val))


def grad_dot(u, v):
    """Return ``grad(u) . grad(v)`` as a NGSolve form expression.
    Helper to keep the BilinearForm body readable.
    """
    from ngsolve import grad, InnerProduct
    return InnerProduct(grad(u), grad(v))


# -----------------------------------------------------------------
# CLI
# -----------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Transient heat-transfer solver for IH workpieces.")
    parser.add_argument("--wp-vol", required=True,
                        help="Workpiece volume mesh (.vol).")
    parser.add_argument("--heat-flux-boundaries", default="",
                        help="Required boundary name or NGSolve boundary "
                             "expression receiving q_surf, for example "
                             "'heated_outer|heated_end'.")
    parser.add_argument("--convection-boundaries", default="",
                        help="Boundary expression receiving Newton "
                             "convection. Required when --h-conv is nonzero.")
    parser.add_argument("--radiation-boundaries", default="",
                        help="Boundary expression receiving radiation. "
                             "Required when --emissivity is nonzero.")
    parser.add_argument("--surface-label", default=None,
                        help=argparse.SUPPRESS)
    # Material thermal properties.
    parser.add_argument("--material", default="steel",
                        choices=list(THERMAL_PRESETS) + ["custom"],
                        help="Thermal material preset.")
    parser.add_argument("--rho", type=float, default=None,
                        help="Density [kg/m^3] (overrides preset).")
    parser.add_argument("--cp", type=float, default=None,
                        help="Specific heat [J/(kg.K)] (overrides preset).")
    parser.add_argument("--k", type=float, default=None,
                        help="Conductivity [W/(m.K)] (overrides preset).")
    parser.add_argument("--em-table", default="",
                        help="(|H_t|, T) surface-impedance table from "
                             "calc_em_table.py.  Makes q_surf follow the "
                             "local surface temperature: q = q_EM * "
                             "q_tab(s|H_t|, T) / q_tab(|H_t|, T_ref).")
    parser.add_argument("--ht-sol", default="",
                        help="|H_t| field (_Ht.sol) of the EM run on "
                             "--em-vol.  Without it |H_t| is inferred from "
                             "q_surf through the table at T_ref, and the "
                             "table-vs-EM consistency check is skipped.")
    parser.add_argument("--ht-scale", type=float, default=1.0,
                        help="Scale of |H_t| (coil current ratio I/I_ref) "
                             "for --em-table.")
    parser.add_argument("--em-reference-temperature", type=float,
                        default=None,
                        help="Workpiece temperature [degC] of the EM run "
                             "(default: --t-initial).")
    parser.add_argument("--em-table-azimuths", type=int, default=64,
                        help="Azimuths of |H_t| per ring for --em-table with "
                             "--q-phi-average.")
    parser.add_argument("--allow-em-table-extrapolation",
                        action="store_true",
                        help="Clamp temperatures / |H_t| outside --em-table "
                             "instead of failing; excursions are reported.")
    parser.add_argument("--allow-frozen-ht", action="store_true",
                        help="Acknowledge that --em-table freezes the EM "
                             "run's |H_t| (invalid for a ferromagnetic part "
                             "crossing its Curie band; see the validation "
                             "record coupled_curie_cylinder_frozen_ht.json).")
    parser.add_argument("--em-table-tolerance", type=float, default=0.05,
                        help="Allowed relative power difference between the "
                             "table at T_ref and the EM run (with --ht-sol).")
    parser.add_argument("--material-table", default="",
                        help="CSV with header T_C,k_W_mK,cp_J_kgK giving "
                             "temperature-dependent conductivity and "
                             "specific heat (density stays --rho/preset).  "
                             "Switches to the enthalpy Newton integrator.")
    parser.add_argument("--latent-heat", type=float, default=0.0,
                        help="Latent heat [J/kg] released uniformly over "
                             "--latent-range (needs --material-table).")
    parser.add_argument("--latent-range", default="",
                        help="T_start,T_end [degC] of the latent release.")
    parser.add_argument("--allow-table-extrapolation", action="store_true",
                        help="Hold the end values of --material-table "
                             "outside its range instead of failing; the "
                             "largest excursion is reported.")
    parser.add_argument("--newton-tol", type=float, default=1.0e-3,
                        help="Newton step tolerance [K] of the nonlinear "
                             "integrator.")
    parser.add_argument("--newton-max-iter", type=int, default=25,
                        help="Newton iterations per step; a step that does "
                             "not converge fails the run.")
    parser.add_argument("--max-halvings", type=int, default=0,
                        help="Allow a non-converged step to be split in "
                             "halves up to this many times (default 0: "
                             "fail).  Every split is recorded.")

    # Boundary conditions.
    parser.add_argument("--h-conv", type=float, default=10.0,
                        help="Newton convection coefficient [W/(m^2.K)].")
    parser.add_argument("--t-ext", type=float, default=20.0,
                        help="External temperature for convection [degC].")
    parser.add_argument("--emissivity", type=float, default=0.0,
                        help="Surface emissivity for radiation "
                             "eps*sigma*(T^4-T_ext^4) [0..1]; 0 = off "
                             "(radiation ambient = --t-ext).")
    parser.add_argument("--t-initial", type=float, default=20.0,
                        help="Initial workpiece temperature [degC].")

    # Heat source: pick exactly one mode.
    parser.add_argument("--q-uniform", type=float, default=None,
                        help="Uniform surface heat flux [W/m^2] "
                             "(testing / first-cut mode).")
    parser.add_argument("--qsurf-sol", default="",
                        help="q_surf .sol from calc_fem_kelvin.py "
                             "(spatial distribution).")
    parser.add_argument("--em-vol", default="",
                        help="EM .vol the qsurf-sol corresponds to.  "
                             "REQUIRED when --qsurf-sol is supplied; "
                             "NGSolve .sol is a coefficient vector only "
                             "(no embedded mesh), so the EM .vol must "
                             "be passed explicitly.  Auto-detection from "
                             "the .sol stem was removed 2026-05-20.")
    parser.add_argument("--qsurf-order", type=int, default=1,
                        help="Fixed H1 order for the cross-mesh qsurf.sol "
                             "handoff (only 1 is supported).")
    parser.add_argument("--q-phi-average", action="store_true",
                        help="Circumferentially (phi) average the spatial "
                             "--qsurf-sol into an AXISYMMETRIC q_surf "
                             "(uniform in phi) on the 3D mesh, and solve "
                             "without rotation time-stepping.  This is the "
                             "steady limit a fast-spinning workpiece "
                             "converges to.  Requires --qsurf-sol (not "
                             "--q-uniform).  The complementary 'no-rotation' "
                             "mode is simply --qsurf-sol with "
                             "--rotation-rpm 0 (the spatial q applied as-is, "
                             "non-axisymmetric).")
    parser.add_argument("--exposure-thresholds", default="",
                        help="Comma-separated temperatures [degC]; the "
                             "result reports the volume above each, where "
                             "it is and how it is distributed.")
    parser.add_argument("--temperature-limit", type=float, default=None,
                        help="Temperature limit [degC] (for example the "
                             "solidus of the workpiece alloy).  The result "
                             "reports whether and by how much it is "
                             "exceeded, as a constraint for optimisers.")
    parser.add_argument("--q-phi-average-n", default=None,
                        help=argparse.SUPPRESS)
    parser.add_argument("--q-phi-bin", type=float, default=None,
                        help="Meridian bin width [m] of the exact "
                             "--q-phi-average (default: a quarter of the "
                             "median EM surface edge).")
    parser.add_argument("--em-heat-boundaries", default="",
                        help="EM boundary expression carrying q_surf.  "
                             "Default: the .sol sidecar, else the thermal "
                             "--heat-flux-boundaries names when they exist "
                             "on the EM mesh.")
    parser.add_argument("--q-scale", type=float, default=1.0,
                        help="Multiply q_surf by this factor (for a current "
                             "sweep, (I/I_ref)^2).  Recorded in the result.")
    parser.add_argument("--transfer-tolerance", type=float, default=None,
                        help="Largest allowed distance [m] from a thermal "
                             "heat-flux vertex to the EM heated surface "
                             "(default: a quarter of the median EM edge).")
    parser.add_argument("--power-tolerance", type=float, default=None,
                        help="Largest allowed relative difference between "
                             "the thermal heat input and the EM source "
                             "power (default 0.02).")

    # Time integration.
    parser.add_argument("--dt", type=float, default=0.5,
                        help="Time step [s].")
    parser.add_argument("--t-end", type=float, default=5.0,
                        help="End time [s].")
    parser.add_argument("--time-scheme", default="backward-euler",
                        choices=["backward-euler", "crank-nicolson"],
                        help="Time integration scheme.")
    parser.add_argument("--linear-solver", default="sparsecholesky",
                        choices=["sparsecholesky"],
                        help="NGSolve direct solver for the inverse "
                             "of M + theta*dt*K.")

    # FES.
    parser.add_argument("--fes-order", type=int, default=1,
                        help="H1 polynomial order (default 1).")

    # Workpiece rotation: with a spatial --qsurf-sol the source is
    # re-evaluated on the rotated body every time step.
    parser.add_argument("--rotation-rpm", type=float, default=0.0,
                        help="Workpiece rotation speed [rpm] "
                             "(default 0 = stationary).  With a spatial "
                             "--qsurf-sol the source is re-evaluated on "
                             "the rotated body each time step; "
                             "--q-phi-average gives its fast-rotation "
                             "limit.")
    parser.add_argument("--rotor-states", default="",
                        help="JSON manifest (radia.ih-rotor-states/1) of one "
                             "EM solution per rotor angle, for a rotating "
                             "part that is not a body of revolution (cross "
                             "holes, flats, keys).  Replaces --qsurf-sol / "
                             "--em-vol; needs --rotation-rpm > 0.")
    parser.add_argument("--angle-step-tolerance", type=float, default=None,
                        help="Largest relative change of the source between "
                             "neighbouring rotor angles (default 0.25); a "
                             "larger change means the angles are too coarse.")
    parser.add_argument("--rotation-axis", default="z",
                        choices=["x", "y", "z"],
                        help="Workpiece rotation axis (default z).  "
                             "Positive --rotation-rpm gives CCW "
                             "rotation viewed from the positive end of "
                             "this axis (right-hand rule).  Pre-v4.78.0 "
                             "this was hardcoded to z; horizontal-axis "
                             "workpieces (billet along x) silently got "
                             "the wrong physics.")

    # Observation / output.
    parser.add_argument("--probe-point", default="",
                        help="Probe point 'x,y,z' [m] for the T(t) "
                             "history (default: none).")
    parser.add_argument("--msh-output", default="",
                        help="Optional .msh path for Open GMSH "
                             "(final-step T field).")
    parser.add_argument("--csv-output", default="",
                        help="Optional CSV path for the probe T(t) "
                             "history.")

    def run(args):
        if args.surface_label is not None:
            return {"error":
                    "--surface-label was removed because it coupled heat "
                    "input, convection, and radiation to one ambiguous "
                    "boundary set. Use --heat-flux-boundaries and "
                    "--convection-boundaries; also pass "
                    "--radiation-boundaries when emissivity is nonzero."}
        if (args.q_uniform is None) and (not args.qsurf_sol):
            return {"error":
                    "Either --q-uniform or --qsurf-sol is required."}
        try:
            thresholds = [float(t) for t in
                          str(args.exposure_thresholds).split(",")
                          if t.strip()]
        except ValueError:
            return {"error": "--exposure-thresholds must be comma-separated "
                             f"numbers, got {args.exposure_thresholds!r}"}
        latent_range = None
        if args.latent_range:
            try:
                latent_range = tuple(float(t) for t in
                                     args.latent_range.split(","))
                if len(latent_range) != 2:
                    raise ValueError
            except ValueError:
                return {"error": "--latent-range must be T_start,T_end"}
        if args.q_phi_average_n is not None:
            return {"error":
                    "--q-phi-average-n was removed: --q-phi-average now "
                    "integrates the EM source exactly over the azimuth "
                    "instead of sampling a fixed number of angles.  Use "
                    "--q-phi-bin to set its meridian resolution."}
        probe_point = None
        if args.probe_point:
            try:
                probe_point = [float(s) for s in args.probe_point.split(",")]
                if len(probe_point) != 3:
                    raise ValueError("probe_point must be x,y,z")
            except Exception as e:
                return {"error":
                        f"--probe-point parse error: {e}"}
        return solve_heat(
            wp_vol=args.wp_vol,
            material=args.material, rho=args.rho, cp=args.cp, k=args.k,
            h_conv=args.h_conv, t_ext=args.t_ext, t_initial=args.t_initial,
            emissivity=args.emissivity,
            heat_flux_boundaries=args.heat_flux_boundaries,
            convection_boundaries=args.convection_boundaries,
            radiation_boundaries=args.radiation_boundaries,
            q_uniform=args.q_uniform,
            qsurf_sol=args.qsurf_sol,
            em_vol=args.em_vol,
            qsurf_order=args.qsurf_order,
            dt=args.dt, t_end=args.t_end,
            time_scheme=args.time_scheme,
            linear_solver=args.linear_solver,
            fes_order=args.fes_order,
            rotation_rpm=args.rotation_rpm,
            rotation_axis=args.rotation_axis,
            rotor_states=args.rotor_states,
            angle_step_tolerance=args.angle_step_tolerance,
            q_phi_average=args.q_phi_average,
            q_phi_bin=args.q_phi_bin,
            em_heat_boundaries=args.em_heat_boundaries,
            q_scale=args.q_scale,
            transfer_tolerance=args.transfer_tolerance,
            power_tolerance=args.power_tolerance,
            exposure_thresholds=thresholds,
            temperature_limit=args.temperature_limit,
            material_table=args.material_table,
            em_table=args.em_table,
            ht_sol=args.ht_sol,
            ht_scale=args.ht_scale,
            em_reference_temperature=args.em_reference_temperature,
            em_table_azimuths=args.em_table_azimuths,
            allow_em_table_extrapolation=args.allow_em_table_extrapolation,
            em_table_tolerance=args.em_table_tolerance,
            allow_frozen_ht=args.allow_frozen_ht,
            latent_heat=args.latent_heat,
            latent_range=latent_range,
            allow_table_extrapolation=args.allow_table_extrapolation,
            newton_tol=args.newton_tol,
            newton_max_iter=args.newton_max_iter,
            max_halvings=args.max_halvings,
            probe_point=probe_point,
            msh_output=args.msh_output,
            csv_output=args.csv_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
