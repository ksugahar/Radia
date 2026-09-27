"""
2D axisymmetric transient heat solver for IH workpieces (Phase B v1.5).

When the workpiece is rotationally symmetric about the Z axis (a
cylinder, a stepped shaft, a disk, ...) the 3D thermal problem can
be reduced to a 2D problem in the (r, z) plane with all integrals
weighted by 2*pi*r.  This is 50-100x faster than the 3D solve and
matches Kubota's cylinder workflow exactly.

Geometry convention
-------------------
  - The workpiece thermal mesh is a 2D Netgen ``.vol`` whose
    coordinates are (r, z, 0) with r >= 0.
  - Heating, convection, and radiation boundary roles are selected
    independently and must be explicit when active.
  - The Z axis (r = 0) is the natural axisymmetric BC.  No DOFs
    are constrained there; the (2*pi*r) weight collapses to 0
    automatically.

q_surf cross-mesh transfer
--------------------------
The EM .vol is 3D.  Its heated surface is averaged exactly over the azimuth
by :class:`ih_thermal.AxisymmetricSurfaceProfile`: each EM triangle is cut
at meridian arc-length bin edges and integrated, so the average conserves
the EM power and has no azimuthal sampling error.  Each heat-flux vertex of
the (r, z) mesh takes the profile value at its meridian position; a vertex
farther than the transfer tolerance from the EM meridian is an error, and
so is a revolved heat input that differs from the EM power by more than
``--power-tolerance``.  ``--n-phi-samples`` now drives an independent
pointwise ring check that is reported alongside.

Output schema mirrors calc_heat.py so the Heat panel does not need
to know which solver produced the JSON.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

# Shared utilities
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import setup_paths, progress, calc_main  # noqa: E402

# Axisym shares its thermal preset table with the 3D solver.
from calc_heat import (  # noqa: E402
    SIGMA_SB,
    THERMAL_PRESETS,
    _boundary_role_audit,
    _resolve_boundary_role,
    _resolve_material,
    _input_mesh_geometry_audit,
    _locate_probe,
    _probe_value,
    _temperature_extrema,
    _validate_qsurf_transfer_order,
)


def _log(msg):
    progress("HEAT_AXI", msg)


def _reject_axis_boundary_roles(mesh, role_specs):
    """Reject active physical surface roles that include the ``r = 0`` axis.

    The axis has zero revolved area because every boundary integral is weighted
    by ``2*pi*r``.  Selecting it is therefore never a meaningful heat-flux,
    convection, or radiation configuration.  Inspect boundary *elements*, not
    just names, so a broad label shared by the axis and a physical surface also
    fails loudly.
    """
    from ngsolve import BND

    radial_scale = max(
        (abs(float(vertex.point[0])) for vertex in mesh.vertices),
        default=0.0,
    )
    axis_tol = max(1.0e-14, radial_scale * 1.0e-12)

    for option_name, boundary_names in role_specs:
        selected = set(boundary_names)
        axis_labels = set()
        for element in mesh.Elements(BND):
            if element.mat not in selected:
                continue
            radii = [
                abs(float(mesh.vertices[vertex.nr].point[0]))
                for vertex in element.vertices
            ]
            if radii and max(radii) <= axis_tol:
                axis_labels.add(element.mat)
        if axis_labels:
            labels = sorted(axis_labels)
            raise ValueError(
                f"{option_name} includes r=0 axis boundary elements through "
                f"{labels}. In an axisymmetric (r, z) heat model the center "
                "axis is the natural symmetry boundary and has zero-revolved-"
                "area (2*pi*r*ds = 0). Remove the axis from this selector and "
                "give the coil-facing or exposed physical surface a separate "
                "boundary label."
            )


# -----------------------------------------------------------------
# q_surf source for the axisym mesh
# -----------------------------------------------------------------

def _build_axisym_qsurf_gf(wp_mesh, heat_flux_boundary_names, args):
    """Return ``(gf, q_cf, audit)`` for the axisymmetric heating boundary.

    The 3D EM source is averaged exactly over the azimuth
    (:class:`ih_thermal.AxisymmetricSurfaceProfile`) and evaluated at the
    meridian position of each heat-flux vertex.  Every vertex must lie on
    the EM meridian within the transfer tolerance and the revolved heat
    input must match the EM power within ``power_tolerance``.

    ``args.n_phi_samples`` sets an independent check: the source is also
    sampled pointwise on ``n_phi_samples`` azimuths through each vertex and
    the ring mean is compared with the profile (reported, not enforced,
    because pointwise sampling aliases narrow azimuthal features).

    A scalar ``--q-uniform`` short-circuits the projection.
    """
    from ngsolve import (H1, GridFunction, CoefficientFunction, BND,
                         Integrate, x as r_coord)
    import ih_thermal

    if args.q_uniform is not None:
        _log(f"Q_SURF:uniform {args.q_uniform:.4e} W/m^2")
        return (None, CoefficientFunction(float(args.q_uniform)),
                {"mode": "uniform"})

    if not args.qsurf_sol:
        raise ValueError(
            "Either --q-uniform or --qsurf-sol is required.")
    # --em-vol must be explicit (NGSolve .sol is a coefficient vector
    # only -- no embedded mesh).
    if not args.em_vol:
        raise ValueError(
            "--em-vol is required when --qsurf-sol is supplied.  "
            "NGSolve .sol files do not contain mesh information, "
            "so the EM .vol that the .sol was saved against must "
            "be passed explicitly.")
    _validate_qsurf_transfer_order(args.qsurf_order)

    if isinstance(args.n_phi_samples, bool):
        raise ValueError("--n-phi-samples must be a positive integer")
    try:
        n_phi = int(args.n_phi_samples)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "--n-phi-samples must be a positive integer"
        ) from exc
    if n_phi < 1 or str(n_phi) != str(args.n_phi_samples).strip():
        raise ValueError("--n-phi-samples must be a positive integer")

    names = sorted(heat_flux_boundary_names)
    source = ih_thermal.EMHeatSource(
        args.qsurf_sol, args.em_vol, thermal_names=names,
        em_boundaries=getattr(args, "em_heat_boundaries", None) or None,
        mode="phi-average", axis="z",
        q_scale=float(getattr(args, "q_scale", 1.0) or 1.0),
        transfer_tolerance=getattr(args, "transfer_tolerance", None),
        bin_width=getattr(args, "q_phi_bin", None))
    power_tol = float(getattr(args, "power_tolerance", None)
                      or ih_thermal.DEFAULT_POWER_TOLERANCE)
    _log(f"Q_SURF:{os.path.basename(source.qsurf_sol)} on "
         f"{os.path.basename(source.em_vol)} "
         f"({source.pair_audit['provenance']}), EM boundaries "
         f"{source.em_boundaries} ({source.boundary_rule}), "
         f"P_EM={source.power:.6e} W")

    vnrs = ih_thermal.boundary_vertex_numbers(wp_mesh, names)
    rz = ih_thermal.mesh_vertices(wp_mesh)[vnrs][:, :2]
    values, dist = source.values_at_meridian(rz)
    gf_wp_q = GridFunction(H1(wp_mesh, order=1))
    vec = gf_wp_q.vec.FV().NumPy()
    vec[:] = 0.0
    vec[vnrs] = values
    region = wp_mesh.Boundaries("|".join(names))
    power = float(Integrate(gf_wp_q * 2.0 * math.pi * r_coord, wp_mesh, BND,
                            definedon=region).real)
    balance = ih_thermal.power_balance(
        power, source.power, power_tol,
        what="revolved thermal heat input against the EM source power",
        scale=source.surface.magnitude_power,
        hint=ih_thermal.TRANSFER_POWER_HINT)

    # Independent pointwise ring check.
    phis = 2.0 * math.pi * np.arange(n_phi) / n_phi
    ring = np.c_[
        (rz[:, :1] * np.cos(phis)[None, :]).ravel(),
        (rz[:, :1] * np.sin(phis)[None, :]).ravel(),
        np.repeat(rz[:, 1], n_phi)]
    tri, lam, d_ring = source.surface.locate(ring)
    ring_vals = np.einsum("ij,ij->i", lam,
                          source.surface.values[source.surface.tris[tri]])
    inside = (d_ring <= source.transfer_tolerance).reshape(-1, n_phi)
    ring_vals = ring_vals.reshape(-1, n_phi)
    full = inside.all(axis=1)
    scale = max(float(np.max(np.abs(values))), 1e-300)
    ring_dev = (float(np.max(np.abs(ring_vals[full].mean(axis=1)
                                    - values[full])) / scale)
                if np.any(full) else None)

    audit = source.audit()
    audit.update({
        "evaluation_region": "BND-meridian",
        "target_surface_vertices": int(len(vnrs)),
        "max_transfer_distance_m": float(dist.max()),
        "power_tolerance": power_tol,
        "power_balance": balance,
        "pointwise_ring_check": {
            "n_phi_samples": n_phi,
            "vertices_fully_on_source": int(full.sum()),
            "max_relative_deviation": ring_dev,
        },
    })
    _log(f"Q_SURF:phi-average to {len(vnrs)} vertices, max distance "
         f"{dist.max():.3e} m, revolved P/P_EM - 1 = "
         f"{balance['relative_error']:+.3e}, pointwise {n_phi}-azimuth ring "
         f"check {ring_dev if ring_dev is None else f'{ring_dev:.3e}'}")
    if getattr(args, "em_table", ""):
        q_ref = np.asarray(gf_wp_q.vec.FV().NumPy())[vnrs].copy()
        T_ref = getattr(args, "em_reference_temperature", None)
        if T_ref is None:
            raise ValueError("--em-reference-temperature is required with "
                             "--em-table")
        coupled, coupled_audit = ih_thermal.build_temperature_dependent_source(
            source, table_path=args.em_table,
            target_xyz=np.c_[rz[:, 0], np.zeros(len(rz)), rz[:, 1]],
            q_ref=q_ref, gf_target=gf_wp_q, dofs=vnrs, T_ref=float(T_ref),
            H_scale=float(getattr(args, "ht_scale", 1.0) or 1.0),
            ht_sol=getattr(args, "ht_sol", "") or "",
            azimuths=int(getattr(args, "em_table_azimuths", 64) or 64),
            allow_extrapolation=bool(getattr(
                args, "allow_em_table_extrapolation", False)),
            consistency_tolerance=float(getattr(
                args, "em_table_tolerance", None) or 0.05),
            acknowledge_frozen_ht=bool(getattr(args, "allow_frozen_ht",
                                               False)))
        audit["temperature_dependent_source"] = coupled_audit
        audit["_coupled"] = coupled
        _log(f"Q_SURF:temperature-dependent source from "
             f"{os.path.basename(args.em_table)} (T_ref={float(T_ref):g} C, "
             f"|H_t| {coupled_audit['H_t_range_A_m'][0]:.3e}.."
             f"{coupled_audit['H_t_range_A_m'][1]:.3e} A/m)")
    return gf_wp_q, gf_wp_q, audit


# -----------------------------------------------------------------
# Axisymmetric heat solve
# -----------------------------------------------------------------

def solve_heat_axisym(wp_vol,
                      material="steel", rho=None, cp=None, k=None,
                      h_conv=10.0, t_ext=20.0, t_initial=20.0, emissivity=0.0,
                      heat_flux_boundaries="",
                      convection_boundaries="",
                      radiation_boundaries="",
                      q_uniform=None, qsurf_sol="", em_vol="",
                      qsurf_order=1, n_phi_samples=128,
                      em_heat_boundaries="", q_scale=1.0,
                      transfer_tolerance=None, power_tolerance=None,
                      q_phi_bin=None,
                      exposure_thresholds=(), temperature_limit=None,
                      material_table="", latent_heat=0.0,
                      em_table="", ht_sol="", ht_scale=1.0,
                      em_reference_temperature=None,
                      em_table_azimuths=64,
                      allow_em_table_extrapolation=False,
                      em_table_tolerance=0.05,
                      allow_frozen_ht=False,
                      latent_range=None,
                      allow_table_extrapolation=False,
                      newton_tol=1.0e-3, newton_max_iter=25,
                      dt=0.5, t_end=5.0,
                      time_scheme="backward-euler",
                      linear_solver="sparsecholesky",
                      fes_order=2,
                      rotation_rpm=0.0,
                      probe_point=None,
                      msh_output="",
                      csv_output="",
                      _wp_mesh=None,
                      _write_solution=True):
    setup_paths()
    t0 = time.perf_counter()

    from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                          CF, ds, dx, x as r_coord,
                          TaskManager, InnerProduct, Integrate, grad)

    if _wp_mesh is None and not os.path.isfile(wp_vol):
        return {"error": f"--wp-vol not found: {wp_vol}"}

    wp_mesh = Mesh(wp_vol) if _wp_mesh is None else _wp_mesh
    if wp_mesh.dim != 2:
        return {"error":
                f"--wp-vol is {wp_mesh.dim}D; axisym needs a 2D mesh "
                f"in the (r, z) plane.  Use calc_heat.py for 3D."}
    # --- Workpiece-only mesh contract (radia-ih thermal, axisym) ------
    # The axisym thermal step targets the WORKPIECE (r, z) cross-section
    # ONLY -- a single region.  Reject an empty mesh or a multi-region
    # (coil+wp) mesh loudly rather than heating the coil.  See
    # calc_heat.py for the 3D twin of this guard.
    if wp_mesh.ne == 0:
        return {"error":
                f"--wp-vol {os.path.basename(wp_vol)} has 0 elements; "
                f"the axisym thermal step needs a 2D (r, z) workpiece "
                f"mesh."}
    _wp_mats = sorted(set(wp_mesh.GetMaterials()))
    if len(_wp_mats) > 1:
        return {"error":
                f"axisymmetric thermal analysis targets the WORKPIECE "
                f"ONLY, but --wp-vol {os.path.basename(wp_vol)} has "
                f"{len(_wp_mats)} material regions {_wp_mats}.  Use a "
                f"workpiece-only (r, z) mesh -- a single region."}
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
    try:
        active_roles = [
            ("--heat-flux-boundaries", heat_flux_names),
        ]
        if float(h_conv) != 0.0:
            active_roles.append(
                ("--convection-boundaries", convection_names)
            )
        if float(emissivity) != 0.0:
            active_roles.append(
                ("--radiation-boundaries", radiation_names)
            )
        _reject_axis_boundary_roles(wp_mesh, active_roles)
    except ValueError as exc:
        return {"error": str(exc)}
    _log(f"BND:heat_flux={heat_flux_names} "
         f"convection={convection_names} radiation={radiation_names}")

    rho_v, cp_v, k_v = _resolve_material(material, rho, cp, k)
    _log(f"MATERIAL:{material} rho={rho_v} cp={cp_v} k={k_v}")

    # In axisym mode the workpiece is rotation-symmetric by
    # construction, so a non-zero rotation_rpm justifies the
    # phi-averaging of the cross-mesh q_surf transfer.  Log the
    # value for the JSON output; the solve itself is unchanged
    # (rotation is implicit in the axisymmetric assumption).
    if float(rotation_rpm) > 0.0:
        _log(f"ROTATION:rpm={float(rotation_rpm):g} axisym -- the exact "
             "circumferential average models the long-time-average heat "
             "input from a spinning workpiece.")

    # Standard NGSolve H1 + 2 pi r weighting (FEMM-canonical; matches
    # the heat solver in FEMM's hsolv/prob1big.cpp -- standard P1
    # triangle on physical (r, z) with the 2 pi r Jacobian evaluated
    # at the element centroid).
    #
    # Per CLAUDE.md "Axisymmetric FE: Henrotte Basis Only" policy
    # (refined 2026-05-10 after surveying FEMM 4.2 source):
    #   * Henrotte basis is REQUIRED for axisymmetric MAGNETIC solves
    #     (curl operator brings 1/r axis singularity that standard FE
    #     cannot integrate accurately near the axis).
    #   * Standard H1 + 2 pi r is FINE for axisymmetric SCALAR solves
    #     (heat, electric potential, diffusion) because the weak form
    #     contains 2 pi r as a smooth Jacobian, NOT a 1/r integrand.
    #     Meeker uses standard P1 triangle for FEMM's heat solver and
    #     ships production accuracy.
    #
    # The `radia.axifem.AxiHenrotteHeat{Stiffness,Mass}BFI`
    # classes (added in radia 4.31.0) remain available as optional
    # parity-conscious infrastructure; they are not used here because
    # the FEMM convention says we don't need them for scalar T.
    #
    # Default order is 2 (2026-09-03 near-axis study): P1/Q1 cannot
    # represent the even-parity dT/dr = 0 at the r = 0 axis -- the
    # near-axis profile shows a cusp (apparent slope 16/3x the exact
    # secant, mesh-size-independent shape) and T(axis) is off by
    # O(h^2).  Order 2 contains r^2 and removes the cusp.
    fes_T = H1(wp_mesh, order=int(fes_order))
    u, v = fes_T.TnT()
    gfT = GridFunction(fes_T)
    _log(f"FES:H1 order={fes_order} ndof={fes_T.ndof}")
    if int(fes_order) == 1:
        # Order 1 stays available (an explicit user choice, e.g. to
        # reproduce an older result), but it must not look correct by
        # default: standard P1/Q1 cannot represent dT/dr = 0 at the
        # axis, so say so at run time rather than only in --help.
        _log("FES:WARNING order=1 cannot represent dT/dr = 0 at the "
             "r = 0 axis: the near-axis profile shows a cusp "
             "(apparent slope 16/3x the exact secant; the shape does "
             "NOT refine away) and T(axis) is off by O(h^2).  Use "
             "--fes-order 2 (the default) unless you are "
             "deliberately reproducing an order-1 result.")

    class _Args:
        pass
    a_local = _Args()
    a_local.q_uniform = q_uniform
    a_local.qsurf_sol = qsurf_sol
    a_local.em_vol = em_vol
    a_local.qsurf_order = qsurf_order
    a_local.n_phi_samples = n_phi_samples
    a_local.em_heat_boundaries = em_heat_boundaries
    a_local.q_scale = q_scale
    a_local.transfer_tolerance = transfer_tolerance
    a_local.power_tolerance = power_tolerance
    a_local.q_phi_bin = q_phi_bin
    a_local.em_table = em_table
    a_local.ht_sol = ht_sol
    a_local.ht_scale = ht_scale
    a_local.em_reference_temperature = (
        float(t_initial) if em_reference_temperature is None
        else float(em_reference_temperature))
    a_local.em_table_azimuths = em_table_azimuths
    a_local.allow_em_table_extrapolation = allow_em_table_extrapolation
    a_local.em_table_tolerance = em_table_tolerance
    a_local.allow_frozen_ht = allow_frozen_ht
    try:
        gf_q, q_cf, qsurf_projection = _build_axisym_qsurf_gf(
            wp_mesh, set(heat_flux_names), a_local)
    except (ValueError, FileNotFoundError) as exc:
        return {"error": str(exc)}
    coupled = qsurf_projection.pop("_coupled", None)

    # Uniform initial state.  ``gfT.vec[:] = T0`` is WRONG for order >= 2:
    # the hierarchical H1 edge/face coefficients are not nodal
    # temperatures, so a constant coefficient vector is not a constant
    # field.  Set() interpolates the constant exactly at every order.
    with TaskManager():
        gfT.Set(CF(float(t_initial)))
    weight = 2 * math.pi * r_coord

    import ih_heat_transient
    import ih_thermal_material
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

    if not nonlinear:
        # ------- Bilinear forms (axisym weight = 2*pi*r) -------
        # ``r_coord`` is NGSolve's x global coordinate (= radial coord).
        weight = 2 * math.pi * r_coord
        K_cf = CF(float(k_v))
        rho_cp = CF(float(rho_v) * float(cp_v))

        a_form = BilinearForm(fes_T, symmetric=True)
        a_form += K_cf * InnerProduct(grad(u), grad(v)) * weight * dx
        if float(h_conv) != 0.0:
            a_form += float(h_conv) * v * u * weight * ds(convection_selector)

        m_form = BilinearForm(fes_T, symmetric=True)
        m_form += rho_cp * u * v * weight * dx
        with TaskManager():
            a_form.Assemble()
            m_form.Assemble()

        if time_scheme not in ("backward-euler", "crank-nicolson"):
            raise ValueError(
                f"Unsupported --time-scheme {time_scheme!r}.")
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
            weight=weight, q_source=q_cf,
            q_update=(None if coupled is None else
                      (lambda g: coupled.update(g, coupled_T_dofs))),
            q_damping=None if coupled is None else coupled.gf_damping,
            linear_solver=linear_solver,
            newton_tol_K=float(newton_tol),
            max_newton=int(newton_max_iter))
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

    # ------- Time loop -------
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

    heat_flux_region = wp_mesh.Boundaries(heat_flux_selector)
    heat_flux_audit = _boundary_role_audit(
        wp_mesh, heat_flux_names, q_cf=q_cf, weight=weight)
    convection_audit = _boundary_role_audit(
        wp_mesh, convection_names, weight=weight)
    radiation_audit = _boundary_role_audit(
        wp_mesh, radiation_names, weight=weight)
    A_surf_axisym = float(heat_flux_audit["area_m2"])
    q_int = float(heat_flux_audit["heat_input_W"])
    _log(f"Q_SURF:int q dA = {q_int:.4e} W (axisym area "
         f"{A_surf_axisym:.4e} m^2)")
    Q_input_J = 0.0

    for step in range(1, n_steps + 1):
        t = step * float(dt)
        if nonlinear:
            try:
                with TaskManager():
                    stepper.advance(float(dt))
            except (ValueError, RuntimeError) as exc:
                return {"error": f"t={t:.4g} s: {exc}"}
        else:
            f_form = LinearForm(fes_T)
            f_form += q_cf * v * weight * ds(heat_flux_selector)
            if float(h_conv) != 0.0:
                f_form += float(h_conv) * float(t_ext) * v * weight \
                    * ds(convection_selector)
            if float(emissivity) > 0.0:    # radiation (explicit, prev-step T, in K)
                _TK = gfT + 273.15
                f_form += -float(emissivity) * SIGMA_SB \
                    * (_TK**4 - (float(t_ext) + 273.15)**4) * v * weight \
                    * ds(radiation_selector)
            with TaskManager():
                f_form.Assemble()
                res_vec.data = f_form.vec - a_form.mat * gfT.vec
                gfT.vec.data += float(dt) * (inv * res_vec)
        Q_input_J += q_int * float(dt)
        t_arr.append(t)
        T_max_history.append(
            float(np.max(gfT.vec.FV().NumPy()[vertex_dofs_T])))
        if probe_mip is not None:
            T_probe.append(_probe_value(gfT, probe_mip))
        _log(f"STEP:{step}/{n_steps} t={t:.3f}s "
             f"T_probe={T_probe[-1] if probe_point is not None else 'n/a'}")

    if nonlinear:
        Q_input_J = stepper.audit.energy_in_J
        q_int = stepper.last_heat_input_W
    if coupled is not None:
        qsurf_projection["temperature_dependent_source"].update(
            coupled.audit())

    # Raw order>=2 H1 coefficients are not temperatures, and vertices
    # alone can miss a higher-order field extremum.
    T_min, T_max, T_extrema = _temperature_extrema(
        gfT, wp_mesh, fes_order
    )
    import ih_thermal_post
    thresholds = sorted(set(float(t) for t in (exposure_thresholds or ())))
    if temperature_limit is not None:
        thresholds = sorted(set(thresholds) | {float(temperature_limit)})
    thermal_exposure = ih_thermal_post.thermal_exposure(
        wp_mesh, gfT, thresholds, axisymmetric=True)
    limit_record = (ih_thermal_post.limit_check(thermal_exposure,
                                                temperature_limit)
                    if temperature_limit is not None else None)
    for entry in thermal_exposure["thresholds"]:
        _log(f"EXPOSURE:T>{entry['T_C']:g} C volume "
             f"{entry['volume_m3'] * 1e9:.3f} mm^3 "
             f"({entry['volume_fraction']:.3%})")
    # Physical volume mean.  The computational mesh is a meridian section,
    # so both numerator and denominator require the same 2*pi*r Jacobian.
    revolved_volume = float(Integrate(weight, wp_mesh))
    T_mean = (
        float(Integrate(gfT * weight, wp_mesh)) / revolved_volume
        if revolved_volume > 0.0 else T_max
    )

    # Final-state field export.  GmshPostExport handles 2D meshes
    # natively (z is padded to 0 in the .msh nodes table) so the
    # axisym panel uses the same vol2msh path as the 3D thermal
    # panel.  q_surf overlay is bundled alongside T so the user
    # sees both the input flux and the resulting volume temperature
    # in one GMSH view.
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
            import ih_thermal
            ih_thermal.write_field_sidecar(
                sol_T, mesh_path=vol_T, mesh=wp_mesh,
                fes_order=int(fes_order),
                quantity=ih_thermal.TEMPERATURE_QUANTITY,
                unit=ih_thermal.TEMPERATURE_UNIT,
                extra={"producer": "calc_heat_axisym",
                       "t_end_s": float(t_end)})
            heat_vol_file = vol_T
            sol_entries = [
                {"sol": sol_T, "fes": "H1",
                 "fes_order": int(fes_order),
                 "fes_dim": 1,
                 "name": "T_C", "ncomp": 1},
            ]
            try:
                fes_qg = H1(wp_mesh, order=int(fes_order))
                gf_qg = GridFunction(fes_qg)
                gf_qg.vec[:] = 0
                gf_qg.Set(q_cf, definedon=heat_flux_region)
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
                 f"({len(sol_entries)} fields, 2D axisym)")
        except Exception as e:
            raise RuntimeError(
                f"thermal field export to {msh_output} failed: "
                f"{type(e).__name__}: {e}") from e
    elif _write_solution:
        # No --msh-output: still save T.sol alongside wp_vol so a
        # later evaluation pass can reload it (mirrors qsurf.sol
        # contract on the EM side).
        try:
            base_dir = os.path.dirname(os.path.abspath(wp_vol))
            stem = os.path.splitext(os.path.basename(wp_vol))[0]
            sol_T = os.path.join(
                base_dir, f"{stem}_heat_T.sol").replace("\\", "/")
            gfT.Save(sol_T)
            T_sol_file = sol_T
            import ih_thermal
            ih_thermal.write_field_sidecar(
                sol_T, mesh_path=wp_vol, mesh=wp_mesh,
                fes_order=int(fes_order),
                quantity=ih_thermal.TEMPERATURE_QUANTITY,
                unit=ih_thermal.TEMPERATURE_UNIT,
                extra={"producer": "calc_heat_axisym",
                       "t_end_s": float(t_end)})
            _log(f"T_SOL:wrote {os.path.basename(sol_T)} "
                 f"(no GMSH bundle requested; load with the same "
                 f"wp_vol + H1 order={fes_order}, axisym 2D)")
        except Exception as e:
            raise RuntimeError(
                f"could not save the temperature field: "
                f"{type(e).__name__}: {e}") from e

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
        "surface_area_m2": A_surf_axisym,
        "n_steps": n_steps,
        "dt_s": float(dt),
        "t_end_s": float(t_end),
        "time_scheme": time_scheme,
        "linear_solver": linear_solver,
        "ndof": int(fes_T.ndof),
        "ne": int(wp_mesh.ne),
        "fes_order": int(fes_order),
        "mesh_geometry": mesh_geometry,
        "revolved_volume_m3": revolved_volume,
        "n_phi_samples": int(n_phi_samples),
        "qsurf_projection": qsurf_projection,
        "thermal_exposure": thermal_exposure,
        "temperature_limit": limit_record,
        "T_max_history_C": T_max_history,
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
                     else "qsurf_sol"),
        "qsurf_sol": qsurf_sol if not q_uniform else "",
        "em_vol": em_vol if not q_uniform else "",
        "T_sol_file": T_sol_file,
        "heat_vol_file": heat_vol_file,
        "msh_file": gmsh_file,
        "csv_file": csv_output if (csv_output and probe_point is not None)
                     else "",
        "t_total_s": round(t_total, 2),
        "mesh_type": "axisymmetric",
    }


def main():
    parser = argparse.ArgumentParser(
        description="2D axisymmetric transient heat solver "
                    "for IH workpieces (Phase B v1.5).")
    parser.add_argument("--wp-vol", required=True,
                        help="2D axisymmetric workpiece mesh (.vol) in "
                             "the (r, z) plane.  r >= 0 required.")
    parser.add_argument("--heat-flux-boundaries", default="",
                        help="Required boundary name or NGSolve boundary "
                             "expression receiving q_surf.")
    parser.add_argument("--convection-boundaries", default="",
                        help="Boundary expression receiving Newton "
                             "convection. Required when --h-conv is nonzero.")
    parser.add_argument("--radiation-boundaries", default="",
                        help="Boundary expression receiving radiation. "
                             "Required when --emissivity is nonzero.")
    parser.add_argument("--surface-label", default=None,
                        help=argparse.SUPPRESS)
    parser.add_argument("--material", default="steel",
                        choices=list(THERMAL_PRESETS) + ["custom"],
                        help="Thermal material preset.")
    parser.add_argument("--rho", type=float, default=None,
                        help="Density [kg/m^3] (overrides preset).")
    parser.add_argument("--cp", type=float, default=None,
                        help="Specific heat [J/(kg.K)] (overrides preset).")
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
                        help="Newton iterations per step before the step is "
                             "halved.")
    parser.add_argument("--k", type=float, default=None,
                        help="Conductivity [W/(m.K)] (overrides preset).")
    parser.add_argument("--h-conv", type=float, default=10.0)
    parser.add_argument("--t-ext", type=float, default=20.0)
    parser.add_argument("--emissivity", type=float, default=0.0,
                        help="Surface emissivity for radiation "
                             "eps*sigma*(T^4-T_ext^4) [0..1]; 0 = off "
                             "(radiation ambient = --t-ext).")
    parser.add_argument("--t-initial", type=float, default=20.0)
    parser.add_argument("--q-uniform", type=float, default=None,
                        help="Uniform surface heat flux [W/m^2].")
    parser.add_argument("--qsurf-sol", default="",
                        help="q_surf .sol from calc_fem_kelvin.py.")
    parser.add_argument("--em-vol", default="",
                        help="EM .vol the qsurf-sol corresponds to.  "
                             "REQUIRED when --qsurf-sol is supplied; "
                             "auto-detection from the .sol stem was "
                             "removed 2026-05-20.")
    parser.add_argument("--qsurf-order", type=int, default=1)
    parser.add_argument("--n-phi-samples", type=int, default=128,
                        help="Azimuths of the independent pointwise ring "
                             "check reported next to the exact "
                             "circumferential average (default 128).")
    parser.add_argument("--exposure-thresholds", default="",
                        help="Comma-separated temperatures [degC]; the "
                             "result reports the volume above each, where "
                             "it is and how it is distributed.")
    parser.add_argument("--temperature-limit", type=float, default=None,
                        help="Temperature limit [degC] (for example the "
                             "solidus of the workpiece alloy).  The result "
                             "reports whether and by how much it is "
                             "exceeded, as a constraint for optimisers.")
    parser.add_argument("--q-phi-bin", type=float, default=None,
                        help="Meridian bin width [m] of the circumferential "
                             "average (default: a quarter of the median EM "
                             "surface edge).")
    parser.add_argument("--em-heat-boundaries", default="",
                        help="EM boundary expression carrying q_surf.  "
                             "Default: the .sol sidecar, else the thermal "
                             "heat-flux names when they exist on the EM "
                             "mesh.")
    parser.add_argument("--q-scale", type=float, default=1.0,
                        help="Multiply q_surf by this factor (for a current "
                             "sweep, (I/I_ref)^2).  Recorded in the result.")
    parser.add_argument("--transfer-tolerance", type=float, default=None,
                        help="Largest allowed distance [m] from a thermal "
                             "heat-flux vertex to the EM meridian "
                             "(default: a quarter of the median EM edge).")
    parser.add_argument("--power-tolerance", type=float, default=None,
                        help="Largest allowed relative difference between "
                             "the revolved heat input and the EM power "
                             "(default 0.02).")
    parser.add_argument("--dt", type=float, default=0.5)
    parser.add_argument("--t-end", type=float, default=5.0)
    parser.add_argument("--time-scheme", default="backward-euler",
                        choices=["backward-euler", "crank-nicolson"])
    parser.add_argument("--linear-solver", default="sparsecholesky",
                        choices=["sparsecholesky", "umfpack", "pardiso"])
    parser.add_argument("--fes-order", type=int, default=2,
                        help="H1 polynomial order (default 2).  Order 1 "
                             "cannot represent dT/dr = 0 at the r = 0 "
                             "axis (near-axis cusp; T(axis) error "
                             "O(h^2)); order 2 removes it.")
    parser.add_argument("--rotation-rpm", type=float, default=0.0,
                        help="Workpiece rotation [rpm] (default 0). "
                             "Recorded for metadata + justifies the "
                             "phi-averaging of cross-mesh q_surf "
                             "transfer when > 0.")
    parser.add_argument("--probe-point", default="",
                        help="Probe point 'r,z' [m] for the T(t) "
                             "history (axisym is 2D, so 2 coords).")
    parser.add_argument("--msh-output", default="")
    parser.add_argument("--csv-output", default="")

    def run(args):
        latent_range = None
        if args.latent_range:
            try:
                latent_range = tuple(float(t) for t in
                                     args.latent_range.split(","))
                if len(latent_range) != 2:
                    raise ValueError
            except ValueError:
                return {"error": "--latent-range must be T_start,T_end"}
        try:
            thresholds = [float(t) for t in
                          str(args.exposure_thresholds).split(",")
                          if t.strip()]
        except ValueError:
            return {"error": "--exposure-thresholds must be comma-separated "
                             f"numbers, got {args.exposure_thresholds!r}"}
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
        probe_point = None
        if args.probe_point:
            try:
                parts = [float(s) for s in args.probe_point.split(",")]
                # Accept either "r,z" (2D) or "r,z,0" (3D-style with
                # the third coord ignored).  The mesh is 2D so we
                # pass exactly 2 coords to wp_mesh().
                if len(parts) == 2:
                    probe_point = parts
                elif len(parts) == 3:
                    probe_point = parts[:2]
                else:
                    raise ValueError("probe_point must be r,z or r,z,0")
            except Exception as e:
                return {"error": f"--probe-point parse error: {e}"}
        return solve_heat_axisym(
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
            n_phi_samples=args.n_phi_samples,
            em_heat_boundaries=args.em_heat_boundaries,
            q_scale=args.q_scale,
            transfer_tolerance=args.transfer_tolerance,
            power_tolerance=args.power_tolerance,
            q_phi_bin=args.q_phi_bin,
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
            dt=args.dt, t_end=args.t_end,
            time_scheme=args.time_scheme,
            linear_solver=args.linear_solver,
            fes_order=args.fes_order,
            rotation_rpm=args.rotation_rpm,
            probe_point=probe_point,
            msh_output=args.msh_output,
            csv_output=args.csv_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
