"""
Transient heat-transfer solver for IH workpieces (Phase B of the
Radia-NGSolve thermal pipeline).

Inputs
------
``--wp-vol``          Workpiece volume mesh (.vol).  Has at least one
                      volume material region (``--material-label``) and a
                      heating-side surface label (``--surface-label``).
                      In Kubota's flow this mesh is SEPARATE from the EM
                      mesh: the EM .vol carries the workpiece as a hole
                      with a SIBC face, while this thermal mesh carries
                      the workpiece as a real solid with the same outer
                      surface.

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

def _build_qsurf_cf(wp_mesh, surface_region, args):
    """Return a CoefficientFunction representing q_surf on the heating
    face of the workpiece thermal mesh.

    Three modes, in priority order:
      1. ``--q-uniform`` : constant scalar everywhere on the face.
      2. ``--qsurf-sol`` + ``--em-vol`` : load the EM-mesh GridFunction,
         project onto wp_mesh's surface vertices.
      3. Both omitted : raise (no heat input).
    """
    from ngsolve import (Mesh, H1, GridFunction, CoefficientFunction)

    if args.q_uniform is not None:
        _log(f"Q_SURF:uniform {args.q_uniform:.4e} W/m^2")
        return CoefficientFunction(float(args.q_uniform))

    if not args.qsurf_sol:
        raise ValueError(
            "Either --q-uniform or --qsurf-sol is required (no heat "
            "input was supplied).")

    qsurf_sol = os.path.abspath(args.qsurf_sol)
    if not os.path.isfile(qsurf_sol):
        raise FileNotFoundError(f"--qsurf-sol not found: {qsurf_sol}")

    # Locate the EM .vol that the .sol was saved against.
    em_vol = args.em_vol
    if not em_vol:
        # calc_fem_kelvin.py writes the .vol next to the .sol with
        # the suffix ``_fem.vol`` and the q_surf .sol with
        # ``_qsurf.sol`` -- swap them.
        stem = qsurf_sol[:-len("_qsurf.sol")] \
            if qsurf_sol.endswith("_qsurf.sol") else \
            os.path.splitext(qsurf_sol)[0]
        em_vol = stem + "_fem.vol"
    em_vol = os.path.abspath(em_vol)
    if not os.path.isfile(em_vol):
        raise FileNotFoundError(
            f"--em-vol could not be auto-located ({em_vol}). "
            f"Pass it explicitly.")

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

    surf_vertex_nrs = set()
    for el in wp_mesh.Elements(BND_REGION := __import__("ngsolve").BND):
        if surface_region is None or el.mat == "":
            pass  # skip filter check; iterate boundaries directly below
    # The cleanest enumeration is via mesh.vertices + boundary lookup;
    # NGSolve does not expose "vertices on boundary X" directly, so
    # iterate boundary elements and collect their vertices.
    for el in wp_mesh.Elements(__import__("ngsolve").BND):
        # Filter by boundary region (CF-region match).  surface_region is
        # the Region returned by mesh.Boundaries(label); el.mat is the
        # boundary's name string.
        if args.surface_label and el.mat != args.surface_label:
            continue
        for v in el.vertices:
            surf_vertex_nrs.add(v.nr)

    # Sample em_mesh's GridFunction at each surface vertex.
    n_ok = 0
    n_fail = 0
    for vnr in surf_vertex_nrs:
        v = wp_mesh.vertices[vnr]
        p = v.point
        try:
            em_mip = em_mesh(float(p[0]), float(p[1]), float(p[2]))
            val = gf_q_em(em_mip)
            gf_wp_q.vec.FV()[vnr] = float(getattr(val, "real", val))
            n_ok += 1
        except Exception:
            n_fail += 1

    _log(f"Q_SURF:projected {n_ok}/{n_ok + n_fail} surface "
         f"vertices ({n_fail} outside EM mesh, set to 0)")
    if n_fail > n_ok:
        _log("Q_SURF:WARNING majority of wp surface vertices fell "
             "outside the EM mesh -- check that the two .vol files "
             "describe the same physical workpiece geometry.")

    return gf_wp_q


# -----------------------------------------------------------------
# Time-domain heat solve
# -----------------------------------------------------------------

def solve_heat(wp_vol,
               material="steel", rho=None, cp=None, k=None,
               h_conv=10.0, t_ext=20.0, t_initial=20.0,
               surface_label="outer", material_label=".*",
               q_uniform=None, qsurf_sol="", em_vol="",
               qsurf_order=1,
               dt=0.5, t_end=5.0,
               time_scheme="backward-euler",
               linear_solver="sparsecholesky",
               fes_order=1,
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

    wp_mesh = Mesh(wp_vol).Curve(int(fes_order))
    _log(f"MESH:loaded {os.path.basename(wp_vol)} "
         f"materials={list(wp_mesh.GetMaterials())} "
         f"boundaries={list(wp_mesh.GetBoundaries())}")

    if surface_label not in wp_mesh.GetBoundaries():
        return {"error":
                f"--surface-label {surface_label!r} not found in "
                f"{wp_vol} boundaries={list(wp_mesh.GetBoundaries())}"}

    rho_v, cp_v, k_v = _resolve_material(material, rho, cp, k)
    _log(f"MATERIAL:{material} rho={rho_v} cp={cp_v} k={k_v}")

    # H1 FES on the workpiece volume.  No Dirichlet BC -- the
    # surface flux + Robin convection give a well-posed problem.
    fes_T = H1(wp_mesh, order=int(fes_order))
    u, v = fes_T.TnT()
    gfT = GridFunction(fes_T)
    gfT.vec[:] = float(t_initial)
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
    a_local.surface_label = surface_label
    surface_region = wp_mesh.Boundaries(surface_label)
    q_cf = _build_qsurf_cf(wp_mesh, surface_region, a_local)

    # ------------------- Bilinear forms -------------------
    K_cf = CF(float(k_v))
    rho_cp = CF(float(rho_v) * float(cp_v))

    a_form = BilinearForm(fes_T, symmetric=True)
    a_form += K_cf * grad_dot(u, v) * dx
    a_form += float(h_conv) * v * u * ds(surface_label)
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
        f_form = LinearForm(fes_T)
        f_form += q_cf * v * ds(surface_label)
        f_form += float(h_conv) * float(t_ext) * v * ds(surface_label)
        f_form.Assemble()
        with TaskManager():
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

    # Final stats: peak T over the volume.
    T_arr = np.asarray(gfT.vec.FV().NumPy())
    T_max = float(np.max(T_arr))
    T_min = float(np.min(T_arr))

    # GMSH .msh export of the final temperature field (Open GMSH).
    gmsh_file = ""
    if msh_output:
        try:
            from gmsh_post_export import save_vol_sol_pair, vol2msh
            base_dir = os.path.dirname(os.path.abspath(msh_output))
            stem = os.path.splitext(os.path.basename(msh_output))[0]
            sol_T = os.path.join(base_dir, f"{stem}_T.sol").replace("\\", "/")
            vol_T = os.path.join(base_dir, f"{stem}_heat.vol").replace("\\", "/")
            save_vol_sol_pair(vol_T, sol_T, wp_mesh.ngmesh, gfT)
            sol_entries = [
                {"sol": sol_T, "fes": "H1",
                 "fes_order": int(fes_order),
                 "fes_dim": 1,
                 "name": "T_C", "ncomp": 1},
            ]
            vol2msh(msh_output, vol_T, sol_entries)
            gmsh_file = msh_output
            _log(f"GMSH:wrote {os.path.basename(msh_output)}")
        except Exception as e:
            _log(f"GMSH_ERROR:{type(e).__name__}: {e}")

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
        "h_conv_W_m2K": float(h_conv),
        "t_ext_C": float(t_ext),
        "surface_label": surface_label,
        "q_source": ("uniform" if q_uniform is not None
                     else "qsurf_sol"),
        "qsurf_sol": qsurf_sol if not q_uniform else "",
        "em_vol": em_vol if not q_uniform else "",
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
    parser.add_argument("--surface-label", default="outer",
                        help="Boundary label where q_surf and Newton "
                             "convection are applied (default: outer).")
    parser.add_argument("--material-label", default=".*",
                        help="Volume material regex for the workpiece "
                             "(default: .* -- all materials).")

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
                        help="EM .vol the qsurf-sol corresponds to "
                             "(auto-detected from qsurf-sol stem if blank).")
    parser.add_argument("--qsurf-order", type=int, default=1,
                        help="H1 order used when calc_fem_kelvin.py "
                             "saved qsurf.sol (must match).")

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
            surface_label=args.surface_label,
            material_label=args.material_label,
            q_uniform=args.q_uniform,
            qsurf_sol=args.qsurf_sol,
            em_vol=args.em_vol,
            qsurf_order=args.qsurf_order,
            dt=args.dt, t_end=args.t_end,
            time_scheme=args.time_scheme,
            linear_solver=args.linear_solver,
            fes_order=args.fes_order,
            probe_point=probe_point,
            msh_output=args.msh_output,
            vtu_prefix=args.vtu_prefix,
            csv_output=args.csv_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
