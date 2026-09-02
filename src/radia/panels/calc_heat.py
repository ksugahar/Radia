"""
Transient heat-transfer solver for IH workpieces (Phase B of the
Radia-NGSolve thermal pipeline).

Inputs
------
``--wp-vol``          Workpiece volume mesh (.vol).  MUST be a
                      WORKPIECE-ONLY 3D volume mesh: a single solid with
                      exactly one volume material region and a
                      heating-side surface label (``--surface-label``).
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

``--em-vol``          The EM .vol the qsurf-sol was computed on.
                      Needed to reconstruct the GridFunction on the
                      EM mesh before sampling.  Falls back to
                      ``<qsurf-sol stem>_fem.vol`` (the file
                      ``calc_fem_kelvin.py`` writes alongside) when
                      omitted.

``--q-uniform``       Alternative to qsurf-sol: a UNIFORM scalar
                      heat flux [W/m^2].  Useful for testing and
                      rough optimization ladders ("does this much
                      power even get to soak temperature?").

Physics
-------
Heat equation with surface heat flux Neumann BC + Newton convection
on the same surface (the SIBC twin):

    rho * cp * dT/dt = div(k grad T)                    (volume)
    -k dT/dn = q_surf(x,y,z) - h_conv (T - T_ext)        (heating face)

Discretization: backward Euler in time, H1 in space.

    (M + dt*K) T_{n+1} = M T_n + dt*(q_surf*v + h*T_ext*v)_ds

Output
------
JSON summary on stdout (last line) with peak temperature, probe
history and timings.  Per-step ``<msh-output stem>_T_NNN.vtu`` files
when ``--msh-output`` is supplied.
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


def _log(msg):
    progress("HEAT", msg)


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


# -----------------------------------------------------------------
# q_surf source: spatial (.sol) or uniform scalar
# -----------------------------------------------------------------

def _build_qsurf_cf(wp_mesh, args):
    """Return ``(q_cf, resample_fn)`` describing q_surf on the heating
    face of the workpiece thermal mesh.

    ``resample_fn(theta_rad)`` is a closure that, when called with a
    rotation angle ``theta_rad`` (radians), updates the underlying
    GridFunction backing ``q_cf`` so that q_surf is re-sampled with
    the workpiece body rotated by ``+theta_rad`` around the z axis
    relative to the EM frame.  Used by the time loop when
    ``--rotation-rpm > 0`` on the spatial-qsurf path.

    ``resample_fn`` is ``None`` for the uniform path (rotation
    has no effect on a constant source).

    Three modes, in priority order:
      1. ``--q-uniform`` : constant scalar everywhere on the face.
      2. ``--qsurf-sol`` + ``--em-vol`` : load the EM-mesh GridFunction,
         project onto wp_mesh's surface vertices.  Re-projectable.
      3. Both omitted : raise (no heat input).
    """
    from ngsolve import (Mesh, H1, GridFunction, CoefficientFunction, BND)

    if getattr(args, "q_phi_average", False) and args.q_uniform is not None:
        raise ValueError(
            "--q-phi-average averages the SPATIAL --qsurf-sol into an "
            "axisymmetric q; it cannot be combined with --q-uniform "
            "(a constant q is already azimuthally uniform).")

    if args.q_uniform is not None:
        _log(f"Q_SURF:uniform {args.q_uniform:.4e} W/m^2")
        return CoefficientFunction(float(args.q_uniform)), None

    if not args.qsurf_sol:
        raise ValueError(
            "Either --q-uniform or --qsurf-sol is required (no heat "
            "input was supplied).")

    qsurf_sol = os.path.abspath(args.qsurf_sol)
    if not os.path.isfile(qsurf_sol):
        raise FileNotFoundError(f"--qsurf-sol not found: {qsurf_sol}")

    # The EM .vol that the .sol was saved against MUST be supplied
    # explicitly.  Pre-2026-05-20 we auto-located a sibling
    # ``<stem>_fem.vol`` when --em-vol was omitted; that silent fallback
    # picked the wrong file when the user renamed the .sol or copied
    # it to a directory without the sibling, and violated CLAUDE.md
    # "No Fallbacks" / fail-fast policy.  NGSolve's .sol format is a
    # raw coefficient vector (no embedded mesh / no embedded fes
    # order), so the only safe contract is: both files explicitly.
    if not args.em_vol:
        raise ValueError(
            "--em-vol is required when --qsurf-sol is supplied.  "
            "NGSolve .sol files do not contain mesh information, "
            "so the EM .vol that the .sol was saved against must "
            "be passed explicitly.  Typically this is the "
            "``<stem>_fem.vol`` file that calc_fem_kelvin.py writes "
            "next to the ``<stem>_qsurf.sol`` -- pass that path.")
    em_vol = os.path.abspath(args.em_vol)
    if not os.path.isfile(em_vol):
        raise FileNotFoundError(f"--em-vol not found: {em_vol}")

    _log(f"Q_SURF:loading {os.path.basename(qsurf_sol)} on "
         f"{os.path.basename(em_vol)}")

    em_mesh = Mesh(em_vol)
    fes_q_em = H1(em_mesh, order=int(args.qsurf_order))
    gf_q_em = GridFunction(fes_q_em)
    gf_q_em.Load(qsurf_sol)

    # Project onto the wp_mesh surface vertices by point evaluation.
    # The wp thermal mesh and the EM mesh both describe the same
    # physical workpiece outer surface, possibly with different
    # element sizes -- pointwise sampling handles that.
    fes_wp_q = H1(wp_mesh, order=int(args.qsurf_order))
    gf_wp_q = GridFunction(fes_wp_q)
    gf_wp_q.vec[:] = 0

    # Enumerate surface vertices by walking boundary elements (NGSolve
    # does not expose "vertices on boundary X" directly).  Filter by the
    # requested boundary label (``el.mat`` is the BND name); an empty
    # ``surface_label`` means every boundary.
    surf_vertex_nrs = set()
    for el in wp_mesh.Elements(BND):
        if args.surface_label and el.mat != args.surface_label:
            continue
        for v in el.vertices:
            surf_vertex_nrs.add(v.nr)

    # Pre-extract surface-vertex (x,y,z) coordinates once so the
    # resample closure can apply rotation cheaply (no per-step
    # ``wp_mesh.vertices[vnr].point`` Python overhead).
    surf_vnrs_list = sorted(surf_vertex_nrs)
    surf_xyz = [
        tuple(float(c) for c in wp_mesh.vertices[vnr].point)
        for vnr in surf_vnrs_list
    ]

    # Rotation axis selection (v4.78.0+): allow the workpiece to spin
    # around x / y / z instead of the hardcoded +z that pre-2026-05-25
    # silently assumed.  Workpieces exported with a horizontal axis
    # (e.g. billet along x) would silently get the wrong physics under
    # the previous "spin around z" assumption.
    axis_str = getattr(args, "rotation_axis", "z") or "z"
    axis_str = str(axis_str).lower().strip()
    if axis_str not in ("x", "y", "z"):
        raise ValueError(
            f"--rotation-axis must be one of x / y / z (got "
            f"{axis_str!r}).  Workpiece spin around an arbitrary axis "
            f"is not yet supported.")

    def _project(theta_rad: float) -> tuple[int, int]:
        """In-place re-sampling at body-rotation angle ``theta_rad``.

        Returns (n_ok, n_fail) for the most recent projection.  When
        theta_rad == 0 this reproduces the original (single-shot)
        projection used pre-2026-05-20.

        Rotation axis is selected by ``args.rotation_axis`` (x / y / z,
        default z) -- positive ``theta_rad`` is CCW viewed from the
        positive end of the axis (right-hand rule).
        """
        c, s = math.cos(theta_rad), math.sin(theta_rad)
        gf_wp_q.vec[:] = 0
        n_ok = n_fail = 0
        fv = gf_wp_q.vec.FV()
        for vnr, (xb, yb, zb) in zip(surf_vnrs_list, surf_xyz):
            # Workpiece body rotates +theta_rad around the chosen axis;
            # compute the world coordinate of body point (xb, yb, zb).
            if axis_str == "z":
                xw = xb * c - yb * s
                yw = xb * s + yb * c
                zw = zb
            elif axis_str == "y":
                # Rotation around y: x'=c*x+s*z, y'=y, z'=-s*x+c*z
                xw = xb * c + zb * s
                yw = yb
                zw = -xb * s + zb * c
            else:  # axis_str == "x"
                # Rotation around x: x'=x, y'=c*y-s*z, z'=s*y+c*z
                xw = xb
                yw = yb * c - zb * s
                zw = yb * s + zb * c
            try:
                em_mip = em_mesh(xw, yw, zw)
                val = gf_q_em(em_mip)
                fv[vnr] = float(getattr(val, "real", val))
                n_ok += 1
            except Exception:
                n_fail += 1
        return n_ok, n_fail

    def _phi_average(n_angles: int) -> tuple[int, int]:
        """Circumferential (phi) average of the spatial q_surf.

        For each surface vertex, sample the EM q_surf at ``n_angles``
        body-rotation angles evenly spaced over 2*pi (around
        ``args.rotation_axis``) and average over the in-EM-mesh samples,
        producing an AXISYMMETRIC (phi-independent) q_surf written into
        gf_wp_q.  This is the steady limit a fast-spinning workpiece
        converges to -- computed once, without rotation time-stepping.
        """
        m = max(1, int(n_angles))
        accum = [0.0] * len(surf_vnrs_list)
        cnt = [0] * len(surf_vnrs_list)
        for kk in range(m):
            th = 2.0 * math.pi * kk / m
            c, s = math.cos(th), math.sin(th)
            for i, (xb, yb, zb) in enumerate(surf_xyz):
                if axis_str == "z":
                    xw = xb * c - yb * s; yw = xb * s + yb * c; zw = zb
                elif axis_str == "y":
                    xw = xb * c + zb * s; yw = yb; zw = -xb * s + zb * c
                else:  # x
                    xw = xb; yw = yb * c - zb * s; zw = yb * s + zb * c
                try:
                    val = gf_q_em(em_mesh(xw, yw, zw))
                    accum[i] += float(getattr(val, "real", val))
                    cnt[i] += 1
                except Exception:
                    pass
        gf_wp_q.vec[:] = 0
        fv = gf_wp_q.vec.FV()
        n_ok = 0
        for i, vnr in enumerate(surf_vnrs_list):
            if cnt[i] > 0:
                fv[vnr] = accum[i] / cnt[i]
                n_ok += 1
        return n_ok, len(surf_vnrs_list) - n_ok

    # uniform / phi-average path: write an axisymmetric q once and return
    # resample_fn=None (the time loop treats a None resampler as a static
    # source, exactly like --q-uniform).
    if getattr(args, "q_phi_average", False):
        n_avg = int(getattr(args, "q_phi_average_n", 48) or 48)
        n_ok, n_fail = _phi_average(n_avg)
        _log(f"Q_SURF:phi-averaged (uniform/axisymmetric) {n_ok}/"
             f"{n_ok + n_fail} surface vertices over n={n_avg} angles "
             f"around {axis_str} -- no rotation time-stepping.")
        if n_fail > n_ok:
            _log("Q_SURF:WARNING majority of wp surface vertices fell "
                 "outside the EM mesh -- check the two .vol files.")
        return gf_wp_q, None

    # Initial projection at theta=0 (original behaviour) so callers
    # that ignore the resampler get the same q_cf they had pre-2026-05-20.
    n_ok, n_fail = _project(0.0)
    _log(f"Q_SURF:projected {n_ok}/{n_ok + n_fail} surface "
         f"vertices ({n_fail} outside EM mesh, set to 0)")
    if n_fail > n_ok:
        _log("Q_SURF:WARNING majority of wp surface vertices fell "
             "outside the EM mesh -- check that the two .vol files "
             "describe the same physical workpiece geometry.")

    def _resample(theta_rad: float) -> None:
        """Public resampler — updates gf_wp_q.vec in place."""
        _project(theta_rad)

    return gf_wp_q, _resample


# -----------------------------------------------------------------
# Time-domain heat solve
# -----------------------------------------------------------------

SIGMA_SB = 5.670374419e-8     # Stefan-Boltzmann constant [W/m^2/K^4]


def solve_heat(wp_vol,
               material="steel", rho=None, cp=None, k=None,
               h_conv=10.0, t_ext=20.0, t_initial=20.0, emissivity=0.0,
               surface_label="",
               q_uniform=None, qsurf_sol="", em_vol="",
               qsurf_order=1,
               q_phi_average=False, q_phi_average_n=48,
               dt=0.5, t_end=5.0,
               time_scheme="backward-euler",
               linear_solver="sparsecholesky",
               fes_order=1,
               rotation_rpm=0.0,
               rotation_axis="z",
               probe_point=None,
               msh_output="",
               vtu_prefix="",
               csv_output=""):
    """Run the transient heat solve.  See module docstring for inputs."""
    setup_paths()
    t0 = time.perf_counter()

    # Lazy imports so --help is fast.
    from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                          Integrate, CF, CoefficientFunction, ds, dx, BND,
                          VTKOutput, TaskManager)

    if not os.path.isfile(wp_vol):
        return {"error": f"--wp-vol not found: {wp_vol}"}

    wp_mesh = Mesh(wp_vol)

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

    wp_mesh.Curve(int(fes_order))
    _log(f"MESH:loaded {os.path.basename(wp_vol)} "
         f"materials={list(wp_mesh.GetMaterials())} "
         f"boundaries={list(wp_mesh.GetBoundaries())}")

    # Resolve the BND filter.  Empty surface_label means "apply qsurf
    # + convection to ALL boundary elements" -- the common case for a
    # single-workpiece .vol where there is exactly one BND label (e.g.
    # 'sibc' for IH workpieces) and asking the user to retype it just
    # to match the panel default 'outer' is friction.  Pass ".*" to
    # NGSolve's regex-based Boundaries() / ds() to match every BND.
    if surface_label:
        if surface_label not in wp_mesh.GetBoundaries():
            return {"error":
                    f"--surface-label {surface_label!r} not found in "
                    f"{wp_vol} boundaries={list(wp_mesh.GetBoundaries())}"}
        surface_label_eff = surface_label
        _log(f"BND:filter={surface_label!r}")
    else:
        surface_label_eff = ".*"
        _log(f"BND:filter=ALL (surface_label empty -- all "
             f"{len(set(wp_mesh.GetBoundaries()))} BND labels: "
             f"{sorted(set(wp_mesh.GetBoundaries()))})")

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
    class _Args:
        pass
    a_local = _Args()
    a_local.q_uniform = q_uniform
    a_local.qsurf_sol = qsurf_sol
    a_local.em_vol = em_vol
    a_local.qsurf_order = qsurf_order
    a_local.surface_label = surface_label_eff
    a_local.rotation_axis = rotation_axis
    a_local.q_phi_average = q_phi_average
    a_local.q_phi_average_n = q_phi_average_n
    surface_region = wp_mesh.Boundaries(surface_label_eff)
    q_cf, q_resample = _build_qsurf_cf(wp_mesh, a_local)

    # Rotation control: when --rotation-rpm > 0 AND a resampler is
    # available (spatial qsurf only -- uniform is rotation-invariant),
    # re-project q_cf at the workpiece body's instantaneous angle each
    # timestep.  Mesh / FES / stiffness / mass are held fixed; only the
    # LinearForm RHS depends on q_cf and is re-Assembled per step.
    omega_mech = (2.0 * math.pi / 60.0) * float(rotation_rpm)
    rotation_active = (omega_mech > 0.0) and (q_resample is not None)
    if float(rotation_rpm) > 0.0 and q_resample is None:
        _log("ROTATION:rpm>0 has no effect here -- q_surf is azimuthally "
             "uniform (--q-uniform constant, or --q-phi-average already "
             "gives the rotation-averaged axisymmetric q).")
    elif rotation_active:
        _log(f"ROTATION:rpm={float(rotation_rpm):g} "
             f"omega={omega_mech:.4f} rad/s -- "
             f"resampling qsurf on the body frame each step.")

    # ------------------- Bilinear forms -------------------
    K_cf = CF(float(k_v))
    rho_cp = CF(float(rho_v) * float(cp_v))

    a_form = BilinearForm(fes_T, symmetric=True)
    a_form += K_cf * grad_dot(u, v) * dx
    a_form += float(h_conv) * v * u * ds(surface_label_eff)
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

    # ------------------- Time loop -------------------
    t_arr = [0.0]
    T_probe = []
    if probe_point is not None:
        try:
            mip = wp_mesh(*[float(c) for c in probe_point])
            T_probe.append(float(gfT(mip).real if hasattr(gfT(mip), "real")
                                  else gfT(mip)))
        except Exception:
            T_probe.append(float("nan"))

    n_steps = int(math.ceil(t_end / dt))
    Q_input_J = 0.0
    A_surf = float(Integrate(CF(1), wp_mesh, BND, definedon=surface_region).real)
    q_int = float(Integrate(q_cf, wp_mesh, BND,
                             definedon=surface_region).real)
    _log(f"Q_SURF:int q_surf dA = {q_int:.4e} W "
         f"(area {A_surf:.4e} m^2)")

    vtu_files = []
    if vtu_prefix:
        vtk = VTKOutput(wp_mesh, coefs=[gfT], names=["T"],
                        filename=f"{vtu_prefix}_000", subdivision=0,
                        legacy=False)
        try:
            vtk.Do()
            vtu_files.append(f"{vtu_prefix}_000.vtu")
        finally:
            del vtk

    for step in range(1, n_steps + 1):
        t = step * float(dt)
        if rotation_active:
            # Body has rotated by omega_mech * t around z at the start
            # of this step.  Re-sample qsurf on the wp surface at that
            # body orientation.  q_cf (a GridFunction) is updated in
            # place; q_int (the integrated heat input) tracks below.
            q_resample(omega_mech * t)
        f_form = LinearForm(fes_T)
        f_form += q_cf * v * ds(surface_label_eff)
        f_form += float(h_conv) * float(t_ext) * v * ds(surface_label_eff)
        if float(emissivity) > 0.0:        # radiation (explicit, prev-step T, in K)
            _TK = gfT + 273.15
            f_form += -float(emissivity) * SIGMA_SB \
                * (_TK**4 - (float(t_ext) + 273.15)**4) * v * ds(surface_label_eff)
        f_form.Assemble()
        with TaskManager():
            res_vec.data = f_form.vec - a_form.mat * gfT.vec
            gfT.vec.data += float(dt) * (inv * res_vec)
        if rotation_active:
            # q_int can drift with rotation if the EM-frame hotspot
            # only partially overlaps the wp surface at some angles.
            # Re-integrate to keep Q_input_J honest.
            q_int = float(Integrate(q_cf, wp_mesh, BND,
                                     definedon=surface_region).real)
        Q_input_J += q_int * float(dt)
        t_arr.append(t)
        if probe_point is not None:
            try:
                mip = wp_mesh(*[float(c) for c in probe_point])
                val = gfT(mip)
                T_probe.append(float(getattr(val, "real", val)))
            except Exception:
                T_probe.append(float("nan"))
        if vtu_prefix:
            vtk = VTKOutput(wp_mesh, coefs=[gfT], names=["T"],
                            filename=f"{vtu_prefix}_{step:03d}",
                            subdivision=0, legacy=False)
            try:
                vtk.Do()
                vtu_files.append(f"{vtu_prefix}_{step:03d}.vtu")
            finally:
                del vtk
        _log(f"STEP:{step}/{n_steps} t={t:.3f}s "
             f"T_probe={T_probe[-1] if probe_point is not None else 'n/a'}")

    # Final stats use physical field samples.  Raw order>=2 H1
    # coefficients are not temperatures, and vertices alone can miss an
    # edge-, face-, or cell-mode extremum.
    T_min, T_max, T_extrema = _temperature_extrema(
        gfT, wp_mesh, fes_order
    )
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
            # already living in gf_q from _build_qsurf_cf; in either
            # case Set on the surface region is the right thing.
            try:
                fes_qg = H1(wp_mesh, order=int(fes_order))
                gf_qg = GridFunction(fes_qg)
                gf_qg.vec[:] = 0
                gf_qg.Set(q_cf,
                           definedon=wp_mesh.Boundaries(surface_label_eff))
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
                 f"({len(sol_entries)} fields)")
        except Exception as e:
            _log(f"GMSH_ERROR:{type(e).__name__}: {e}")
    else:
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
            # heat_vol_file stays "" -- in this branch the wp_vol
            # itself IS the companion mesh (no separate _heat.vol
            # is written because there is no GMSH bundle to anchor).
            _log(f"T_SOL:wrote {os.path.basename(sol_T)} "
                 f"(no GMSH bundle requested; load with the same "
                 f"wp_vol + H1 order={fes_order})")
        except Exception as e:
            _log(f"T_SOL_ERROR:{type(e).__name__}: {e}")

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
        "surface_area_m2": A_surf,
        "n_steps": n_steps,
        "dt_s": float(dt),
        "t_end_s": float(t_end),
        "time_scheme": time_scheme,
        "linear_solver": linear_solver,
        "ndof": int(fes_T.ndof),
        "ne": int(wp_mesh.ne),
        "fes_order": int(fes_order),
        "material": material,
        "rho_kg_m3": float(rho_v),
        "cp_J_kgK": float(cp_v),
        "k_W_mK": float(k_v),
        "rotation_rpm": float(rotation_rpm),
        "h_conv_W_m2K": float(h_conv),
        "t_ext_C": float(t_ext),
        "emissivity": float(emissivity),
        "surface_label": surface_label,
        "q_source": ("uniform" if q_uniform is not None
                     else ("qsurf_sol_phi_average" if q_phi_average
                           else "qsurf_sol")),
        "q_phi_average": bool(q_phi_average),
        "q_phi_average_n": int(q_phi_average_n) if q_phi_average else 0,
        "qsurf_sol": qsurf_sol if not q_uniform else "",
        "em_vol": em_vol if not q_uniform else "",
        "T_sol_file": T_sol_file,
        "heat_vol_file": heat_vol_file,
        "msh_file": gmsh_file,
        "vtu_files": vtu_files,
        "csv_file": csv_output if (csv_output and probe_point is not None)
                     else "",
        "t_total_s": round(t_total, 2),
    }


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
    parser.add_argument("--surface-label", default="",
                        help="Boundary label where q_surf and Newton "
                             "convection are applied.  Leave empty "
                             "(default) to apply to ALL BND -- the "
                             "single-workpiece common case where "
                             "explicitly naming the sole BND label is "
                             "pure friction.  Pass a specific name "
                             "(e.g. 'top' / 'sides') only when the "
                             "workpiece has multiple BND sidesets and "
                             "heating + convection should be restricted "
                             "to a subset.")
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
                        help="H1 order used when calc_fem_kelvin.py "
                             "saved qsurf.sol (must match).")
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
    parser.add_argument("--q-phi-average-n", type=int, default=48,
                        help="Number of azimuthal samples for "
                             "--q-phi-average (default 48).")

    # Time integration.
    parser.add_argument("--dt", type=float, default=0.5,
                        help="Time step [s].")
    parser.add_argument("--t-end", type=float, default=5.0,
                        help="End time [s].")
    parser.add_argument("--time-scheme", default="backward-euler",
                        choices=["backward-euler", "crank-nicolson"],
                        help="Time integration scheme.")
    parser.add_argument("--linear-solver", default="sparsecholesky",
                        choices=["sparsecholesky", "umfpack", "pardiso"],
                        help="NGSolve direct solver for the inverse "
                             "of M + theta*dt*K.")

    # FES.
    parser.add_argument("--fes-order", type=int, default=1,
                        help="H1 polynomial order (default 1).")

    # Workpiece rotation (metadata for the 3D solver; the result is
    # "frozen at one azimuthal configuration" -- q_surf is not
    # rotated per timestep).  For a physically spinning workpiece in
    # the time average use the axisym solver with phi-averaging.
    parser.add_argument("--rotation-rpm", type=float, default=0.0,
                        help="Workpiece rotation speed [rpm] "
                             "(default 0 = stationary).  3D solver "
                             "treats this as metadata.")
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
    parser.add_argument("--vtu-prefix", default="",
                        help="Optional prefix for per-step VTU output "
                             "(produces <prefix>_NNN.vtu).")
    parser.add_argument("--csv-output", default="",
                        help="Optional CSV path for the probe T(t) "
                             "history.")

    def run(args):
        if (args.q_uniform is None) and (not args.qsurf_sol):
            return {"error":
                    "Either --q-uniform or --qsurf-sol is required."}
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
            surface_label=args.surface_label,
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
            q_phi_average=args.q_phi_average,
            q_phi_average_n=args.q_phi_average_n,
            probe_point=probe_point,
            msh_output=args.msh_output,
            vtu_prefix=args.vtu_prefix,
            csv_output=args.csv_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
