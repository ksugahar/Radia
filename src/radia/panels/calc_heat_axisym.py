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
  - The "outer" boundary is the heating face -- usually the curve
    at r = R_workpiece.  Top / bottom curves participate when they
    are also heated.
  - The Z axis (r = 0) is the natural axisymmetric BC.  No DOFs
    are constrained there; the (2*pi*r) weight collapses to 0
    automatically.

q_surf cross-mesh transfer
--------------------------
The EM .vol from calc_fem_kelvin is 3D.  For each surface vertex
(r, z) on the axisym mesh we sample the 3D GridFunction at multiple
azimuth angles (default 8) around the circle (r*cos(phi),
r*sin(phi), z) and average.  When the EM model is also rotationally
symmetric (the typical IH solenoid) the values agree across phi
within numerical noise; when it is not (gapped torus) the average
is the right physical quantity for the axisym thermal model.

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
from calc_heat import THERMAL_PRESETS, _resolve_material  # noqa: E402


def _log(msg):
    progress("HEAT_AXI", msg)


# -----------------------------------------------------------------
# q_surf source for the axisym mesh
# -----------------------------------------------------------------

def _build_axisym_qsurf_gf(wp_mesh, surface_label, args):
    """Return an H1 GridFunction on the axisym mesh whose values on
    ``surface_label`` are the phi-averaged 3D q_surf.

    A scalar ``--q-uniform`` short-circuits the projection.
    """
    from ngsolve import (Mesh, H1, GridFunction, CoefficientFunction, BND)

    if args.q_uniform is not None:
        _log(f"Q_SURF:uniform {args.q_uniform:.4e} W/m^2")
        return None, CoefficientFunction(float(args.q_uniform))

    if not args.qsurf_sol:
        raise ValueError(
            "Either --q-uniform or --qsurf-sol is required.")

    qsurf_sol = os.path.abspath(args.qsurf_sol)
    if not os.path.isfile(qsurf_sol):
        raise FileNotFoundError(f"--qsurf-sol not found: {qsurf_sol}")

    em_vol = args.em_vol
    if not em_vol:
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

    # Build axisym GridFunction.  Order matches qsurf order so the
    # round trip preserves the spatial detail of the EM solve.
    fes_wp_q = H1(wp_mesh, order=int(args.qsurf_order))
    gf_wp_q = GridFunction(fes_wp_q)
    gf_wp_q.vec[:] = 0

    n_phi = max(1, int(args.n_phi_samples))
    phis = [2 * math.pi * k / n_phi for k in range(n_phi)]

    surf_vertex_nrs = set()
    for el in wp_mesh.Elements(BND):
        if surface_label and el.mat != surface_label:
            continue
        for v in el.vertices:
            surf_vertex_nrs.add(v.nr)

    n_ok = 0
    n_fail = 0
    for vnr in surf_vertex_nrs:
        v = wp_mesh.vertices[vnr]
        # 2D mesh stores axisym coords as (r, z, 0).
        r, z = float(v.point[0]), float(v.point[1])
        accum = 0.0
        n_local = 0
        for phi in phis:
            x3, y3, z3 = r * math.cos(phi), r * math.sin(phi), z
            try:
                em_mip = em_mesh(x3, y3, z3)
                val = gf_q_em(em_mip)
                accum += float(getattr(val, "real", val))
                n_local += 1
            except Exception:
                pass
        if n_local > 0:
            gf_wp_q.vec.FV()[vnr] = accum / n_local
            n_ok += 1
        else:
            n_fail += 1

    _log(f"Q_SURF:phi-averaged {n_ok}/{n_ok + n_fail} surface "
         f"vertices using n_phi={n_phi}")
    if n_fail > n_ok:
        _log("Q_SURF:WARNING majority of axisym surface vertices "
             "fell outside the EM mesh -- check that the axisym "
             "geometry mirrors the 3D EM workpiece.")

    return gf_wp_q, gf_wp_q


# -----------------------------------------------------------------
# Axisymmetric heat solve
# -----------------------------------------------------------------

def solve_heat_axisym(wp_vol,
                      material="steel", rho=None, cp=None, k=None,
                      h_conv=10.0, t_ext=20.0, t_initial=20.0,
                      surface_label="outer",
                      q_uniform=None, qsurf_sol="", em_vol="",
                      qsurf_order=1, n_phi_samples=8,
                      dt=0.5, t_end=5.0,
                      time_scheme="backward-euler",
                      linear_solver="sparsecholesky",
                      fes_order=1,
                      probe_point=None,
                      msh_output="",
                      vtu_prefix="",
                      csv_output=""):
    setup_paths()
    t0 = time.perf_counter()

    from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                          Integrate, CF, ds, dx, BND, x as r_coord,
                          VTKOutput, TaskManager, InnerProduct, grad)

    if not os.path.isfile(wp_vol):
        return {"error": f"--wp-vol not found: {wp_vol}"}

    wp_mesh = Mesh(wp_vol)
    if wp_mesh.dim != 2:
        return {"error":
                f"--wp-vol is {wp_mesh.dim}D; axisym needs a 2D mesh "
                f"in the (r, z) plane.  Use calc_heat.py for 3D."}
    wp_mesh.Curve(int(fes_order))
    _log(f"MESH:loaded {os.path.basename(wp_vol)} "
         f"materials={list(wp_mesh.GetMaterials())} "
         f"boundaries={list(wp_mesh.GetBoundaries())}")

    if surface_label not in wp_mesh.GetBoundaries():
        return {"error":
                f"--surface-label {surface_label!r} not in "
                f"{list(wp_mesh.GetBoundaries())}"}

    rho_v, cp_v, k_v = _resolve_material(material, rho, cp, k)
    _log(f"MATERIAL:{material} rho={rho_v} cp={cp_v} k={k_v}")

    fes_T = H1(wp_mesh, order=int(fes_order))
    u, v = fes_T.TnT()
    gfT = GridFunction(fes_T)
    gfT.vec[:] = float(t_initial)
    _log(f"FES:H1 order={fes_order} ndof={fes_T.ndof}")

    class _Args:
        pass
    a_local = _Args()
    a_local.q_uniform = q_uniform
    a_local.qsurf_sol = qsurf_sol
    a_local.em_vol = em_vol
    a_local.qsurf_order = qsurf_order
    a_local.n_phi_samples = n_phi_samples
    gf_q, q_cf = _build_axisym_qsurf_gf(wp_mesh, surface_label, a_local)

    # ------- Bilinear forms (axisym weight = 2*pi*r) -------
    # ``r_coord`` is NGSolve's x global coordinate, which is the
    # radial coordinate in our (r, z) convention.
    weight = 2 * math.pi * r_coord
    K_cf = CF(float(k_v))
    rho_cp = CF(float(rho_v) * float(cp_v))

    a_form = BilinearForm(fes_T, symmetric=True)
    a_form += K_cf * InnerProduct(grad(u), grad(v)) * weight * dx
    a_form += float(h_conv) * v * u * weight * ds(surface_label)
    a_form.Assemble()

    m_form = BilinearForm(fes_T, symmetric=True)
    m_form += rho_cp * u * v * weight * dx
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

    surface_region = wp_mesh.Boundaries(surface_label)
    A_surf_axisym = float(
        Integrate(weight, wp_mesh, BND, definedon=surface_region).real)
    q_int = float(
        Integrate(q_cf * weight, wp_mesh, BND,
                  definedon=surface_region).real)
    _log(f"Q_SURF:int q dA = {q_int:.4e} W (axisym area "
         f"{A_surf_axisym:.4e} m^2)")
    Q_input_J = 0.0

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
        f_form += q_cf * v * weight * ds(surface_label)
        f_form += float(h_conv) * float(t_ext) * v * weight \
            * ds(surface_label)
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

    T_arr = np.asarray(gfT.vec.FV().NumPy())
    T_max = float(np.max(T_arr))
    T_min = float(np.min(T_arr))

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
        "surface_area_m2": A_surf_axisym,
        "n_steps": n_steps,
        "dt_s": float(dt),
        "t_end_s": float(t_end),
        "time_scheme": time_scheme,
        "linear_solver": linear_solver,
        "ndof": int(fes_T.ndof),
        "ne": int(wp_mesh.ne),
        "fes_order": int(fes_order),
        "n_phi_samples": int(n_phi_samples),
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
        "mesh_type": "axisymmetric",
    }


def main():
    parser = argparse.ArgumentParser(
        description="2D axisymmetric transient heat solver "
                    "for IH workpieces (Phase B v1.5).")
    parser.add_argument("--wp-vol", required=True,
                        help="2D axisymmetric workpiece mesh (.vol) in "
                             "the (r, z) plane.  r >= 0 required.")
    parser.add_argument("--surface-label", default="outer",
                        help="Boundary curve where q_surf and Newton "
                             "convection are applied (default: outer).")
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
    parser.add_argument("--t-initial", type=float, default=20.0)
    parser.add_argument("--q-uniform", type=float, default=None,
                        help="Uniform surface heat flux [W/m^2].")
    parser.add_argument("--qsurf-sol", default="",
                        help="q_surf .sol from calc_fem_kelvin.py.")
    parser.add_argument("--em-vol", default="",
                        help="EM .vol the qsurf-sol corresponds to.")
    parser.add_argument("--qsurf-order", type=int, default=1)
    parser.add_argument("--n-phi-samples", type=int, default=8,
                        help="Azimuth samples for phi-averaging the 3D "
                             "qsurf onto the axisym mesh (default 8).  "
                             "Use 1 if you know the EM problem is "
                             "exactly rotationally symmetric.")
    parser.add_argument("--dt", type=float, default=0.5)
    parser.add_argument("--t-end", type=float, default=5.0)
    parser.add_argument("--time-scheme", default="backward-euler",
                        choices=["backward-euler", "crank-nicolson"])
    parser.add_argument("--linear-solver", default="sparsecholesky",
                        choices=["sparsecholesky", "umfpack", "pardiso"])
    parser.add_argument("--fes-order", type=int, default=1)
    parser.add_argument("--probe-point", default="",
                        help="Probe point 'r,z' [m] for the T(t) "
                             "history (axisym is 2D, so 2 coords).")
    parser.add_argument("--msh-output", default="")
    parser.add_argument("--vtu-prefix", default="")
    parser.add_argument("--csv-output", default="")

    def run(args):
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
            surface_label=args.surface_label,
            q_uniform=args.q_uniform,
            qsurf_sol=args.qsurf_sol,
            em_vol=args.em_vol,
            qsurf_order=args.qsurf_order,
            n_phi_samples=args.n_phi_samples,
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
