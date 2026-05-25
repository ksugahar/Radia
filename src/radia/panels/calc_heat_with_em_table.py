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

  --ht-source biot  --coil-step <step> ...
        Pure Biot-Savart from the coil STEP (no workpiece
        backreaction).  Not implemented in v1; raises with a clear
        pointer to use --ht-source kelvin until the upstream PEEC
        Biot-Savart-on-arbitrary-surface helper is wired through.

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
from calc_heat import THERMAL_PRESETS, _resolve_material  # noqa: E402


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


def solve_heat_em_table(wp_vol, em_table_path,
                         ht_source, ht_sol, em_vol, ht_order,
                         coil_current, coil_current_csv,
                         I_ref,
                         material="steel", rho=None, cp=None, k=None,
                         h_conv=10.0, t_ext=20.0, t_initial=20.0,
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
        return {"error":
                "--ht-source biot (pure Biot-Savart from --coil-step) "
                "is not implemented in v1.  Use --ht-source kelvin "
                "with a calc_fem_kelvin reference run.  See module "
                "docstring; v2 will wire the upstream PEEC "
                "Biot-on-surface helper."}
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
                             "biot: not yet implemented (raises).")
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
                        help="biot path: coil STEP (v1 not implemented).")

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
            material=args.material,
            rho=args.rho, cp=args.cp, k=args.k,
            h_conv=args.h_conv, t_ext=args.t_ext,
            t_initial=args.t_initial,
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
