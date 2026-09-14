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
The EM .vol from calc_fem_kelvin is 3D.  For each surface vertex
(r, z) on the axisym mesh we sample the 3D GridFunction on its boundary at
multiple azimuth angles (default 128) around the circle (r*cos(phi),
r*sin(phi), z) and average.  When the EM model is also rotationally
symmetric (the typical IH solenoid) the values agree across phi
within numerical noise; when it is not (gapped torus) the average
is the right physical quantity for the axisym thermal model.

Boundary evaluation is essential here.  Independently generated 2D and 3D
meshes generally describe the same CAD surface with slightly different
facets.  Volume-point lookup can therefore reject a valid surface point and
used to leave the corresponding thermal flux node silently at zero.  The
transfer now projects through NGSolve's boundary point locator, records its
coverage, and fails if the two surfaces are not compatible.

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
    """Return an H1 GridFunction on the axisym mesh whose values on
    ``heat_flux_boundary_names`` are the phi-averaged 3D q_surf.

    A scalar ``--q-uniform`` short-circuits the projection.
    """
    from ngsolve import (Mesh, H1, GridFunction, CoefficientFunction, BND)

    if args.q_uniform is not None:
        _log(f"Q_SURF:uniform {args.q_uniform:.4e} W/m^2")
        return (
            None,
            CoefficientFunction(float(args.q_uniform)),
            {
                "mode": "uniform",
                "evaluation_region": "not-applicable",
                "n_phi_samples": 0,
                "target_surface_vertices": 0,
                "requested_samples": 0,
                "accepted_samples": 0,
                "coverage_fraction": 1.0,
                "fully_sampled_vertices": 0,
                "partially_sampled_vertices": 0,
                "failed_vertices": 0,
                "minimum_samples_per_vertex": 0,
            },
        )

    if not args.qsurf_sol:
        raise ValueError(
            "Either --q-uniform or --qsurf-sol is required.")

    qsurf_sol = os.path.abspath(args.qsurf_sol)
    if not os.path.isfile(qsurf_sol):
        raise FileNotFoundError(f"--qsurf-sol not found: {qsurf_sol}")

    # --em-vol must be explicit (NGSolve .sol is a coefficient vector
    # only -- no embedded mesh); see calc_heat.py for the rationale and
    # the 2026-05-20 contract tightening.
    if not args.em_vol:
        raise ValueError(
            "--em-vol is required when --qsurf-sol is supplied.  "
            "NGSolve .sol files do not contain mesh information, "
            "so the EM .vol that the .sol was saved against must "
            "be passed explicitly.")
    em_vol = os.path.abspath(args.em_vol)
    if not os.path.isfile(em_vol):
        raise FileNotFoundError(f"--em-vol not found: {em_vol}")

    _log(f"Q_SURF:loading {os.path.basename(qsurf_sol)} on "
         f"{os.path.basename(em_vol)}")

    em_mesh = Mesh(em_vol)
    qsurf_order = _validate_qsurf_transfer_order(args.qsurf_order)
    fes_q_em = H1(em_mesh, order=qsurf_order)
    gf_q_em = GridFunction(fes_q_em)
    gf_q_em.Load(qsurf_sol)

    # Build axisym GridFunction.  Order matches qsurf order so the
    # round trip preserves the spatial detail of the EM solve.
    fes_wp_q = H1(wp_mesh, order=qsurf_order)
    gf_wp_q = GridFunction(fes_wp_q)
    gf_wp_q.vec[:] = 0

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
    phis = [2 * math.pi * k / n_phi for k in range(n_phi)]

    surf_vertex_nrs = set()
    for el in wp_mesh.Elements(BND):
        if el.mat not in heat_flux_boundary_names:
            continue
        for v in el.vertices:
            surf_vertex_nrs.add(v.nr)

    n_full = 0
    n_partial = 0
    n_fail = 0
    accepted_samples = 0
    minimum_samples = n_phi if surf_vertex_nrs else 0
    failed_coordinates = []
    for vnr in surf_vertex_nrs:
        v = wp_mesh.vertices[vnr]
        # 2D mesh stores axisym coords as (r, z, 0).
        r, z = float(v.point[0]), float(v.point[1])
        accum = 0.0
        n_local = 0
        for phi in phis:
            x3, y3, z3 = r * math.cos(phi), r * math.sin(phi), z
            try:
                # q_surf is a boundary field.  Independent 2D and 3D
                # meshes usually have slightly different facets even when
                # they came from the same CAD body, so volume lookup is not
                # a valid surface-to-surface transfer.
                em_mip = em_mesh(x3, y3, z3, BND)
                val = gf_q_em(em_mip)
                accum += float(getattr(val, "real", val))
                n_local += 1
            except Exception:
                # Individual misses are counted below.  A missing target
                # vertex or poor global coverage is an error, never a
                # zero-flux fallback.
                continue
        accepted_samples += n_local
        minimum_samples = min(minimum_samples, n_local)
        if n_local > 0:
            gf_wp_q.vec.FV()[vnr] = accum / n_local
            if n_local == n_phi:
                n_full += 1
            else:
                n_partial += 1
        else:
            n_fail += 1
            if len(failed_coordinates) < 5:
                failed_coordinates.append([r, z])

    target_vertices = len(surf_vertex_nrs)
    requested_samples = target_vertices * n_phi
    coverage = (
        accepted_samples / requested_samples
        if requested_samples else 0.0
    )
    audit = {
        "mode": "phi-averaged-3d-boundary",
        "evaluation_region": "BND",
        "n_phi_samples": n_phi,
        "target_surface_vertices": target_vertices,
        "requested_samples": requested_samples,
        "accepted_samples": accepted_samples,
        "coverage_fraction": coverage,
        "fully_sampled_vertices": n_full,
        "partially_sampled_vertices": n_partial,
        "failed_vertices": n_fail,
        "minimum_samples_per_vertex": minimum_samples,
    }
    _log(
        f"Q_SURF:phi-averaged boundary transfer: vertices="
        f"{target_vertices - n_fail}/{target_vertices}, samples="
        f"{accepted_samples}/{requested_samples} ({coverage:.1%}), "
        f"n_phi={n_phi}"
    )
    if target_vertices == 0:
        raise ValueError(
            "The selected heat-flux boundaries contain no boundary vertices."
        )
    if n_fail:
        raise ValueError(
            "Axisymmetric q_surf transfer could not locate the 3D EM "
            f"boundary for {n_fail}/{target_vertices} thermal boundary "
            f"vertices (first (r,z) coordinates: {failed_coordinates}). "
            "The 2D thermal meridian must describe the same physical "
            "workpiece surface as --em-vol; no zero-flux fallback was "
            "applied."
        )
    if coverage < 0.80:
        raise ValueError(
            "Axisymmetric q_surf boundary-transfer coverage is only "
            f"{coverage:.1%} ({accepted_samples}/{requested_samples} "
            "azimuth samples). The 2D thermal meridian does not match the "
            "3D EM workpiece closely enough; no partial-transfer fallback "
            "was applied."
        )

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
        _log(f"ROTATION:rpm={float(rotation_rpm):g} axisym -- "
             f"phi-averaging over {int(n_phi_samples)} samples "
             "models the long-time-average heat input from a "
             "spinning workpiece.")

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

    class _Args:
        pass
    a_local = _Args()
    a_local.q_uniform = q_uniform
    a_local.qsurf_sol = qsurf_sol
    a_local.em_vol = em_vol
    a_local.qsurf_order = qsurf_order
    a_local.n_phi_samples = n_phi_samples
    gf_q, q_cf, qsurf_projection = _build_axisym_qsurf_gf(
        wp_mesh, set(heat_flux_names), a_local)

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
        # Uniform initial state.  ``gfT.vec[:] = T0`` is WRONG for
        # order >= 2: the hierarchical H1 edge/face coefficients are
        # not nodal temperatures, so a constant coefficient vector is
        # not a constant field.  Set() interpolates the constant
        # exactly at every order.
        gfT.Set(CF(float(t_initial)))
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

    # ------- Time loop -------
    t_arr = [0.0]
    T_probe = []
    if probe_point is not None:
        try:
            mip = wp_mesh(*[float(c) for c in probe_point])
            T_probe.append(float(getattr(gfT(mip), "real", gfT(mip))))
        except Exception:
            T_probe.append(float("nan"))

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
        f_form = LinearForm(fes_T)
        f_form += q_cf * v * weight * ds(heat_flux_selector)
        if float(h_conv) != 0.0:
            f_form += float(h_conv) * float(t_ext) * v * weight \
                * ds(convection_selector)
        if float(emissivity) > 0.0:        # radiation (explicit, prev-step T, in K)
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
        if probe_point is not None:
            try:
                mip = wp_mesh(*[float(c) for c in probe_point])
                val = gfT(mip)
                T_probe.append(float(getattr(val, "real", val)))
            except Exception:
                T_probe.append(float("nan"))
        _log(f"STEP:{step}/{n_steps} t={t:.3f}s "
             f"T_probe={T_probe[-1] if probe_point is not None else 'n/a'}")

    # Raw order>=2 H1 coefficients are not temperatures, and vertices
    # alone can miss a higher-order field extremum.
    T_min, T_max, T_extrema = _temperature_extrema(
        gfT, wp_mesh, fes_order
    )
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
                _log(f"GMSH_qsurf overlay skipped: "
                     f"{type(e).__name__}: {e}")
            vol2msh(msh_output, vol_T, sol_entries)
            gmsh_file = msh_output
            _log(f"GMSH:wrote {os.path.basename(msh_output)} "
                 f"({len(sol_entries)} fields, 2D axisym)")
        except Exception as e:
            _log(f"GMSH_ERROR:{type(e).__name__}: {e}")
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
            _log(f"T_SOL:wrote {os.path.basename(sol_T)} "
                 f"(no GMSH bundle requested; load with the same "
                 f"wp_vol + H1 order={fes_order}, axisym 2D)")
        except Exception as e:
            _log(f"T_SOL_ERROR:{type(e).__name__}: {e}")

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
            _log(f"CSV_ERROR:{type(e).__name__}: {e}")

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
        "material": material,
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
                        help="Azimuth samples for phi-averaging the 3D "
                             "qsurf onto the axisym mesh (default 128).  "
                             "Use 1 if you know the EM problem is "
                             "exactly rotationally symmetric.")
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
