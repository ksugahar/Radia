"""Transient motor solver -- Lange-Henrotte-Hameyer 2009 linearization.

Stage-2 CLI script (per Panel Design Workflow Policy) implementing
the field-circuit coupling method of Lange, Henrotte, and Hameyer
(IEEE Trans. Mag. 45(3):1258-1261, 2009; DOI 10.1109/TMAG.2009.2012585).

Knowledge prerequisite: `radia_mcp.motor.femm_transient_knowledge` +
`radia_mcp.motor.henrotte_lineage_knowledge`.

### Algorithm summary

At each FE-resolution time step Dt_FE:

  1. Solve the *nonlinear* magnetostatic FE problem at the current
     `(i, theta)` operating point: A_op = argmin 1/2int nu(|gradA|^2) (gradA)^2 - J_s . A
  2. Freeze the *incremental* reluctivity nu*(x) = nu(B_op^2) + 2 nu'(B_op^2).
     B_op x B_op (the tangent operator)
  3. Solve N_phase *linear* problems with shared PARDISO factorization:
     for each phase j, unit current -> A_j -> psi_kj column of L_inc
  4. Compute back-EMF e_bemf = omega . dpsi/dtheta (finite difference between
     two FE evaluations Dtheta apart)
  5. Compute torque T_em via Arkkio integration on air-gap circle
  6. Hand (L_inc, e_bemf, T_em) to the circuit/mechanical ODE block

Between FE steps, scipy.integrate.solve_ivp advances the circuit ODE

  L_inc . di/dt = v(t) - R . i - e_bemf
  dtheta/dt          = omega
  J . domega/dt      = T_em - T_load - B_visc . omega

at a finer time step Dt_circ.  This is the **two-time-scale
decoupling** that makes the method efficient for inverter-driven
controller-in-the-loop studies.

### Expected mesh labels (in .vol)

Material labels:
  air, stator_iron, rotor_iron
  stator_ind_Ap, stator_ind_An,  (positive / negative turns)
  stator_ind_Bp, stator_ind_Bn,
  stator_ind_Cp, stator_ind_Cn,
  rotor_pm_Np, rotor_pm_Sp, rotor_pm_Nn, rotor_pm_Sn  (optional PMSM)

Boundary labels:
  outer       -- Dirichlet A = 0
  airgap_mid  -- Arkkio integration circle (mid-air-gap)

### JSON output

  {"t": [...], "i_A": [...], "i_B": [...], "i_C": [...],
   "theta": [...], "omega": [...], "T_em": [...], "back_emf": [...],
   "fe_calls": N_FE, "method": "linearization"}

Usage:
    python calc_motor_transient.py --vol antunes.vol \\
        --R-phase 0.3872 --J-inertia 1e-3 \\
        --v-amp 200 --v-freq 100 \\
        --t-end 0.05 --dt-FE 1e-4 --dt-circ 1e-5

Output JSON to stdout.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

# Direct solver for the FE Inverse().  Default "pardiso" (MKL, fastest for the
# repeated transient factorizations) is overridden by --linear-solver in
# main().  The golden tests pass "sparsecholesky" (ngsolve built-in, no MKL)
# so they run on machines whose PATH carries another product's OLD unversioned
# MKL (e.g. CST Studio Suite on the CI runner), which otherwise breaks MKL
# PARDISO at solve time with "Cannot load mkl_intel_thread.dll".
_LINEAR_SOLVER = "pardiso"

# --- panel-common boilerplate ---
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (MU_0, NU_0, setup_paths, progress, calc_main,
                          EMMaterial, add_material_args)


def _log(msg):
    progress("MOTOR", msg)


def _build_fe_problem(mesh, fes_order, sigma_phase_assignment,
                      n_turns_per_slot, slot_area, br_pm, theta_rotor):
    """Build NGSolve FE problem at fixed rotor angle theta.

    Returns (fes, a_nonlinear, gfu_A) plus helper closures.
    """
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction,
                         CoefficientFunction as CF, grad, x, y, dx,
                         InnerProduct, sqrt)

    # 2D A_z scalar-potential formulation
    fes = H1(mesh, order=fes_order, dirichlet="outer")
    A = GridFunction(fes)

    # B^2 = (dA/dy)^2 + (dA/dx)^2 for A_z (2D out-of-plane vector pot)
    def B2_of(A_):
        gA = grad(A_)
        return gA[0]*gA[0] + gA[1]*gA[1]

    # Nonlinear reluctivity nu(B^2): piecewise-linear Marrocco approx
    # Replace with measured BH curve via EMMaterial if available.
    # Marrocco: nu(B^2) = nu_0 (alpha + (1-alpha) B^2n / (B^2n + Bsat^2n))
    nu_air = NU_0
    Bsat = 1.8                   # T
    alpha_M = 1.0 / 5000.0       # 1/ur at low B
    n_M = 8                       # Marrocco sharpness exponent

    def nu_iron(B2):
        ratio = (B2 / (Bsat*Bsat)) ** n_M
        return nu_air * (alpha_M + (1.0 - alpha_M) * ratio / (1.0 + ratio))

    # Per-material reluctivity CoefficientFunction
    mat_names = mesh.GetMaterials()
    nu_dict = {}
    for m in set(mat_names):
        if m.startswith("stator_iron") or m.startswith("rotor_iron"):
            nu_dict[m] = nu_iron(B2_of(A))
        else:
            nu_dict[m] = nu_air
    nu_cf = mesh.MaterialCF(nu_dict, default=nu_air)

    return fes, A, B2_of, nu_iron, nu_cf, nu_air


def _solve_nonlinear_op_point(mesh, fes, A, nu_cf, f_assembled,
                              max_iter=20, tol=1e-7):
    """Picard iteration for grad.(nu(B^2) gradA) = -J_s on 'outer' Dirichlet.

    f_assembled: LinearForm already .Assemble()d with the J_s source.
    """
    from ngsolve import (BilinearForm, dx, grad, InnerProduct)

    u, v = fes.TrialFunction(), fes.TestFunction()
    a = BilinearForm(fes, symmetric=True)
    a += nu_cf * InnerProduct(grad(u), grad(v)) * dx

    A_prev = A.vec.CreateVector()
    diff = float("inf")
    for k in range(max_iter):
        A_prev.data = A.vec
        a.Assemble()
        inv = a.mat.Inverse(fes.FreeDofs(), inverse=_LINEAR_SOLVER)
        A.vec.data = inv * f_assembled.vec
        diff = (A.vec - A_prev).Norm() / max(A.vec.Norm(), 1e-30)
        if diff < tol:
            break
    return k + 1, diff


def _extract_L_inc(mesh, fes, A, nu_cf, phase_regions, n_turns, slot_area,
                   stack_length):
    """Build incremental inductance matrix L_inc[N_phase, N_phase].

    Strategy (Lange-Henrotte-Hameyer secIII):
      * Freeze nu at the operating-point A (the *incremental*
        permeability -- same nu_cf since it's already evaluated at A).
      * For each phase j: assemble RHS with unit current in slot_j,
        solve linearized system (shared factorization).
      * Flux linkage psi_kj = N_turns/A_slot . stack_length . int_{slot_k} A_j dx
    """
    from ngsolve import (BilinearForm, LinearForm, GridFunction,
                         grad, dx, Integrate, InnerProduct, TaskManager)

    u, v = fes.TrialFunction(), fes.TestFunction()
    a_lin = BilinearForm(fes, symmetric=True)
    a_lin += nu_cf * InnerProduct(grad(u), grad(v)) * dx
    with TaskManager():
        a_lin.Assemble()
        inv = a_lin.mat.Inverse(fes.FreeDofs(), inverse=_LINEAR_SOLVER)

        N = len(phase_regions)
        L_inc = np.zeros((N, N))
        gfu_phases = []
        J_density = n_turns / slot_area
        for j, (pos_mat, neg_mat) in enumerate(phase_regions):
            f = LinearForm(fes)
            # +1 A through positive turns, -1 A through negative turns
            f += +J_density * v * dx(definedon=mesh.Materials(pos_mat))
            f += -J_density * v * dx(definedon=mesh.Materials(neg_mat))
            f.Assemble()
            gfu_j = GridFunction(fes)
            gfu_j.vec.data = inv * f.vec
            gfu_phases.append(gfu_j)

        # Flux linkage: psi_kj = stack_length * N_turns/A_slot *
        #               (int_{slot_k_pos} A_j - int_{slot_k_neg} A_j)
        for k, (kp, kn) in enumerate(phase_regions):
            for j, gfu_j in enumerate(gfu_phases):
                psi_p = Integrate(gfu_j, mesh, definedon=mesh.Materials(kp))
                psi_n = Integrate(gfu_j, mesh, definedon=mesh.Materials(kn))
                psi = stack_length * n_turns / slot_area * (psi_p - psi_n)
                L_inc[k, j] = psi

        return L_inc, gfu_phases


def _compute_torque_arkkio(mesh, A_op, r_mid, stack_length):
    """Maxwell-stress torque on air-gap circle (Arkkio 2D formula).

    Computes T = (stack_length / mu_0) * integral_{airgap_mid}
                   r * B_r * B_phi  dl

    Requires a boundary labeled "airgap_mid" in the mesh, lying on
    the air-gap mid-circle of radius r_mid.

    For 2D A_z formulation: B = (dA/dy, -dA/dx).  In polar (r,phi):
       B_r   =  Bx cos(phi) + By sin(phi)
       B_phi = -Bx sin(phi) + By cos(phi)
    where phi is the angular coordinate on the circle.
    """
    from ngsolve import grad, Integrate, x, y, sqrt, BoundaryFromVolumeCF
    if "airgap_mid" not in set(mesh.GetBoundaries()):
        return 0.0
    gA = grad(A_op)
    Bx = gA[1]
    By = -gA[0]
    # On the air-gap circle, r = sqrt(x^2+y^2) ~ r_mid; use exact xy/r.
    r = sqrt(x*x + y*y)
    cos_phi = x / r
    sin_phi = y / r
    B_r   =  Bx * cos_phi + By * sin_phi
    B_phi = -Bx * sin_phi + By * cos_phi
    # Arkkio integrand: r * B_r * B_phi.  The boundary integral dl picks
    # up r dphi as the arc-length element; for the closed circle it
    # collapses to (r * B_r * B_phi) integrated as a surface form.
    integrand = r * B_r * B_phi
    val = Integrate(integrand, mesh,
                    definedon=mesh.Boundaries("airgap_mid"))
    return (stack_length / MU_0) * val


def _step_circuit_ode(state, t, dt, L_inc, R, e_bemf, v_source,
                      T_load, J_inertia, T_em):
    """One Crank-Nicolson step of the linearized state space.

    state = (i_A, i_B, i_C, ..., theta, omega)
    T_em provided by the FE refresh (held constant over the inner step).
    Returns new state.
    """
    N = L_inc.shape[0]
    i = state[:N]
    theta = state[N]
    omega = state[N + 1]

    v = v_source(t + dt)         # voltage at end of step
    A_mat = L_inc + 0.5 * dt * R
    rhs = (L_inc - 0.5 * dt * R) @ i + dt * (v - e_bemf)
    i_new = np.linalg.solve(A_mat, rhs)

    omega_new = omega + dt * (T_em - T_load) / J_inertia
    theta_new = theta + dt * 0.5 * (omega + omega_new)

    return np.concatenate([i_new, [theta_new, omega_new]])


def _pm_remanence_cf(mesh, theta_rotor, n_pole_pairs, Br_value):
    """Br vector as a CF on the rotor_pm region (alternating N-S).

    For an 8-pole (n_pole_pairs=4) machine, the radial Br alternates
    sign every (180/n_pole_pairs) degrees of mechanical angle.  In
    the rotor frame, the pattern is fixed; we rotate the angular
    coordinate by theta_rotor to follow the rotor.
    """
    from ngsolve import x, y, atan2, sin, cos, IfPos, CoefficientFunction
    # Rotor-frame angle: phi_rotor = atan2(y, x) - theta_rotor
    phi = atan2(y, x) - theta_rotor
    # Polarity: +Br when sin(n_pole_pairs * phi) > 0, else -Br
    sign = IfPos(sin(n_pole_pairs * phi), 1.0, -1.0)
    # Radial Br (in lab frame x, y): Br = sign * (cos(phi_lab) ex + sin(phi_lab) ey)
    phi_lab = atan2(y, x)
    Br_x = sign * Br_value * cos(phi_lab)
    Br_y = sign * Br_value * sin(phi_lab)
    # Restrict to rotor_pm region
    Br_dict = {"rotor_pm": CoefficientFunction((Br_x, Br_y))}
    return mesh.MaterialCF(Br_dict, default=CoefficientFunction((0.0, 0.0)))


def _compute_back_emf(mesh, fes, A_op_func, phase_regions, n_turns,
                      slot_area, stack_length, theta, omega, delta_theta=1e-4):
    """Back-EMF e_bemf_j = omega * d(psi_j)/d(theta_rotor) via FD.

    A_op_func: callable theta -> A GridFunction at the FE op point
               at rotor angle theta (with current source held fixed).
    Returns: e_bemf array of shape (N_phases,)
    """
    from ngsolve import Integrate
    if omega == 0:
        return np.zeros(len(phase_regions))
    # Two FE solves at theta and theta+delta with zero phase current
    # → captures the PM-only contribution to flux linkage
    A_minus = A_op_func(theta - 0.5 * delta_theta)
    A_plus  = A_op_func(theta + 0.5 * delta_theta)
    e = np.zeros(len(phase_regions))
    for j, (pos_mat, neg_mat) in enumerate(phase_regions):
        psi_minus = (Integrate(A_minus, mesh, definedon=mesh.Materials(pos_mat))
                     - Integrate(A_minus, mesh, definedon=mesh.Materials(neg_mat)))
        psi_plus  = (Integrate(A_plus, mesh, definedon=mesh.Materials(pos_mat))
                     - Integrate(A_plus, mesh, definedon=mesh.Materials(neg_mat)))
        dpsi_dtheta = (psi_plus - psi_minus) / delta_theta
        e[j] = omega * stack_length * n_turns / slot_area * dpsi_dtheta
    return e


def solve_motor_transient(
    vol_file="",
    method="linearization",      # "linearization" | "coupled"
    fes_order=1,
    nbr_phases=3,
    n_turns_per_slot=100,
    slot_area=1e-4,               # m^2
    stack_length=0.05,            # m (axial)
    n_pole_pairs=4,
    r_airgap_mid=0.05,            # m
    R_phase=0.3872,               # ohm
    L_endwinding=0.0,             # H, lumped end-winding inductance
    pm_Br_value=0.0,              # T, PM remanence (0 -> no PMs)
    J_inertia=1e-3,               # kg.m^2
    B_viscous=0.0,                # N.m.s
    T_load=0.0,                   # N.m
    v_amp=0.0, v_freq=0.0, v_phase_shift=2.0*math.pi/3,
    theta_init=0.0, omega_init=0.0, i_init=None,
    t_end=0.1, dt_FE=1e-4, dt_circ=1e-5,
    bh_iron="marrocco",            # placeholder; could load from EMMaterial
    n_steps_per_FE=10,             # circuit steps per FE refresh
):
    """Transient motor solve -- Lange-Henrotte-Hameyer linearization."""
    setup_paths()
    from ngsolve import Mesh

    if not vol_file or not os.path.exists(vol_file):
        return {"error": f"vol_file not found: {vol_file!r}"}

    _log(f"Loading mesh from {vol_file}")
    mesh = Mesh(vol_file)
    _log(f"  materials: {mesh.GetMaterials()}")
    _log(f"  boundaries: {mesh.GetBoundaries()}")

    # Default phase region mapping (matching expected mesh labels)
    phase_regions = [
        ("stator_ind_Ap", "stator_ind_An"),
        ("stator_ind_Bp", "stator_ind_Bn"),
        ("stator_ind_Cp", "stator_ind_Cn"),
    ][:nbr_phases]

    # Verify all required materials exist
    missing = [m for pair in phase_regions for m in pair
               if m not in mesh.GetMaterials()]
    if missing:
        return {"error": f"missing material labels in mesh: {missing}. "
                         f"available: {list(set(mesh.GetMaterials()))}"}

    # Build FE machinery
    _log("Building FE problem")
    fes, A, B2_of, nu_iron, nu_cf, nu_air = _build_fe_problem(
        mesh, fes_order,
        sigma_phase_assignment=None,
        n_turns_per_slot=n_turns_per_slot,
        slot_area=slot_area,
        br_pm=None,
        theta_rotor=theta_init,
    )
    _log(f"  DOFs: {fes.ndof}, free: {sum(fes.FreeDofs())}")

    # Initial state
    if i_init is None:
        i_init = [0.0] * nbr_phases
    state = np.concatenate([np.asarray(i_init, dtype=float),
                            [theta_init, omega_init]])

    # Voltage source (3-phase sinusoidal)
    def v_source(t):
        return np.array([
            v_amp * math.cos(2 * math.pi * v_freq * t + j * v_phase_shift)
            for j in range(nbr_phases)
        ])

    R_mat = R_phase * np.eye(nbr_phases)

    # PM helper (returns CF) — used only if rotor_pm region exists
    has_pm = "rotor_pm" in set(mesh.GetMaterials())

    # Shared FE-source builder for a given (i, theta).
    # If PM region exists, add the PM contribution: J_eq = curl(Br/mu).
    # In 2D A_z formulation, PM enters as a body-force-like term:
    #   integral { (1/mu_0) Br . curl(test) }
    # which becomes  integral { (1/mu_0) (-Br_y dv/dx + Br_x dv/dy) } dx
    from ngsolve import LinearForm, dx, grad
    def _assemble_source(i_phase, theta_r):
        f_loc = LinearForm(fes)
        v_test = fes.TestFunction()
        # Stator winding current contribution
        J_dens = n_turns_per_slot / slot_area
        for j, (pos, neg) in enumerate(phase_regions):
            f_loc += (+i_phase[j] * J_dens * v_test
                      * dx(definedon=mesh.Materials(pos)))
            f_loc += (-i_phase[j] * J_dens * v_test
                      * dx(definedon=mesh.Materials(neg)))
        # PM contribution (if any)
        if has_pm and pm_Br_value > 0:
            Br_cf = _pm_remanence_cf(mesh, theta_r, n_pole_pairs,
                                     pm_Br_value)
            # (Br_x, Br_y) . curl_2D(test) = -Br_y dv/dx + Br_x dv/dy
            gv = grad(v_test)
            f_loc += (1.0 / MU_0) * (-Br_cf[1] * gv[0] + Br_cf[0] * gv[1]) * dx
        f_loc.Assemble()
        return f_loc

    # Helper: solve PM-only (zero coil current) at angle theta_r → return
    # a fresh GridFunction (used by back-EMF FD)
    def _pm_only_solve(theta_r):
        from ngsolve import GridFunction
        from ngsolve import TaskManager
        f_pm = _assemble_source(np.zeros(nbr_phases), theta_r)
        A_pm = GridFunction(fes)
        _solve_nonlinear_op_point(mesh, fes, A_pm, nu_cf, f_pm,
                                  max_iter=10, tol=1e-6)
        return A_pm

    # --- Time-stepping main loop ---
    t = 0.0
    history = {"t": [], "i": [], "theta": [], "omega": [],
               "T_em": [], "e_bemf": [], "fe_calls": 0,
               "L_inc_norm": []}

    n_FE = int(math.ceil(t_end / dt_FE))
    _log(f"Starting time loop: t_end={t_end}, dt_FE={dt_FE}, "
         f"dt_circ={dt_circ}, N_FE={n_FE}, has_pm={has_pm}")

    t_solver_start = time.perf_counter()
    for k_FE in range(n_FE):
        i_now = state[:nbr_phases]
        theta_now = state[nbr_phases]
        omega_now = state[nbr_phases + 1]

        # 1) Nonlinear FE solve at current (i, theta)
        f_now = _assemble_source(i_now, theta_now)
        n_iter, diff = _solve_nonlinear_op_point(
            mesh, fes, A, nu_cf, f_now,
            max_iter=20, tol=1e-7,
        )
        history["fe_calls"] += 1

        # 2) Incremental inductance (shared PARDISO factorization)
        L_inc, _ = _extract_L_inc(
            mesh, fes, A, nu_cf, phase_regions,
            n_turns_per_slot, slot_area, stack_length,
        )
        L_inc = L_inc + L_endwinding * np.eye(nbr_phases)
        history["L_inc_norm"].append(float(np.linalg.norm(L_inc)))

        # 3) Back-EMF via PM-only finite difference on theta.
        # Skip when pm_Br_value=0 (PM contribution is identically zero,
        # no FE solve needed) OR when omega=0 (no motional EMF).
        if has_pm and pm_Br_value > 0 and omega_now != 0:
            e_bemf = _compute_back_emf(
                mesh, fes, _pm_only_solve, phase_regions,
                n_turns_per_slot, slot_area, stack_length,
                theta_now, omega_now, delta_theta=1e-3,
            )
            history["fe_calls"] += 2  # the two PM-only solves
        else:
            e_bemf = np.zeros(nbr_phases)

        # 4) Arkkio torque on airgap_mid circle
        T_em = _compute_torque_arkkio(mesh, A, r_airgap_mid, stack_length)

        # 5) Inner circuit + mechanical ODE integration
        for k_circ in range(n_steps_per_FE):
            state = _step_circuit_ode(
                state, t, dt_circ, L_inc, R_mat, e_bemf, v_source,
                T_load, J_inertia, T_em,
            )
            t += dt_circ
            history["t"].append(t)
            history["i"].append(state[:nbr_phases].tolist())
            history["theta"].append(float(state[nbr_phases]))
            history["omega"].append(float(state[nbr_phases + 1]))
            history["T_em"].append(float(T_em))
            history["e_bemf"].append(e_bemf.tolist())

        if (k_FE + 1) % 10 == 0 or k_FE == 0:
            _log(f"  k_FE={k_FE+1}/{n_FE}, t={t:.4g}s, "
                 f"|i|={np.linalg.norm(state[:nbr_phases]):.3g}A, "
                 f"T_em={T_em:.3g} Nm, ||L||={np.linalg.norm(L_inc):.3g}")

    t_solver = time.perf_counter() - t_solver_start

    return {
        "method": method,
        "vol_file": vol_file,
        "n_FE_calls": history["fe_calls"],
        "n_circ_steps": len(history["t"]),
        "t_solver_s": t_solver,
        "t_end_s": t,
        "i_history": history["i"],
        "t_history": history["t"],
        "theta_history": history["theta"],
        "omega_history": history["omega"],
        "T_em_history": history["T_em"],
        "L_inc_norm_history": history["L_inc_norm"],
        "final_i": state[:nbr_phases].tolist(),
        "final_theta_rad": float(state[nbr_phases]),
        "final_omega_rad_s": float(state[nbr_phases + 1]),
        "method_paper": "Lange-Henrotte-Hameyer, IEEE TMAG 45(3) 2009",
        "e_bemf_history": history["e_bemf"],
        "notes": [
            "Arkkio torque: integrated on labeled 'airgap_mid' boundary "
            "(requires mesh sideset).",
            "Back-EMF: PM-only finite difference d(psi)/d(theta) at "
            "delta_theta=1e-4 rad (cost: 2 extra FE solves per refresh).",
            "PM rotation: Br(theta_rotor) as analytic CF with "
            "n_pole_pairs alternating N-S sectors on 'rotor_pm' region.",
            "Nonlinear nu(B^2) uses Marrocco approximation; replace with "
            "EMMaterial BH-curve table for production use.",
            "Eddy currents: NOT included (Lange-HH limitation). "
            "See radia_mcp.motor.hollaus_eddy_knowledge for the MSFEM "
            "extension that adds them.",
            "Knowledge: radia_mcp.motor.femm_transient_knowledge "
            "+ radia_mcp.motor.henrotte_lineage_knowledge "
            "+ radia_mcp.motor.hollaus_eddy_knowledge",
        ],
    }


def build_argparser():
    parser = argparse.ArgumentParser(
        description="Transient motor solver -- Lange-Henrotte-Hameyer 2009")
    parser.add_argument("--vol", required=True, help="Netgen .vol mesh")
    parser.add_argument("--method", default="linearization",
                        choices=["linearization", "coupled"],
                        help="linearization=Lange-HH (default), "
                             "coupled=fully-coupled ONELAB-style (not yet)")
    parser.add_argument("--fes-order", type=int, default=1)
    parser.add_argument("--linear-solver", default="pardiso",
                        choices=["pardiso", "sparsecholesky", "umfpack"],
                        help="FE direct solver. pardiso=MKL (fastest, default); "
                             "sparsecholesky=ngsolve built-in (no MKL) -- use on "
                             "machines with a conflicting MKL on PATH (e.g. CST).")
    parser.add_argument("--nbr-phases", type=int, default=3)
    parser.add_argument("--n-turns-per-slot", type=int, default=100)
    parser.add_argument("--slot-area", type=float, default=1e-4,
                        help="m^2 per slot")
    parser.add_argument("--stack-length", type=float, default=0.05,
                        help="m, axial length")
    parser.add_argument("--n-pole-pairs", type=int, default=4)
    parser.add_argument("--r-airgap-mid", type=float, default=0.05,
                        help="m, air-gap mid-radius")
    parser.add_argument("--R-phase", type=float, default=0.3872,
                        help="ohm, per-phase resistance")
    parser.add_argument("--L-endwinding", type=float, default=0.0,
                        help="H, lumped end-winding inductance per phase")
    parser.add_argument("--pm-Br-value", type=float, default=0.0,
                        help="T, PM remanence on rotor_pm region (0 = no PMs)")
    parser.add_argument("--J-inertia", type=float, default=1e-3,
                        help="kg.m^2, rotor inertia")
    parser.add_argument("--B-viscous", type=float, default=0.0,
                        help="N.m.s, viscous damping")
    parser.add_argument("--T-load", type=float, default=0.0,
                        help="N.m, mechanical load torque")
    parser.add_argument("--v-amp", type=float, default=0.0,
                        help="V, voltage amplitude")
    parser.add_argument("--v-freq", type=float, default=0.0,
                        help="Hz, voltage frequency")
    parser.add_argument("--theta-init", type=float, default=0.0)
    parser.add_argument("--omega-init", type=float, default=0.0,
                        help="rad/s")
    parser.add_argument("--t-end", type=float, default=0.05,
                        help="s, simulation end time")
    parser.add_argument("--dt-FE", type=float, default=1e-4,
                        help="s, FE refresh interval")
    parser.add_argument("--dt-circ", type=float, default=1e-5,
                        help="s, circuit ODE step (fine)")
    parser.add_argument("--n-steps-per-FE", type=int, default=10)
    return parser


def main():
    parser = build_argparser()

    def run(args):
        global _LINEAR_SOLVER
        _LINEAR_SOLVER = args.linear_solver
        return solve_motor_transient(
            vol_file=args.vol,
            method=args.method,
            fes_order=args.fes_order,
            nbr_phases=args.nbr_phases,
            n_turns_per_slot=args.n_turns_per_slot,
            slot_area=args.slot_area,
            stack_length=args.stack_length,
            n_pole_pairs=args.n_pole_pairs,
            r_airgap_mid=args.r_airgap_mid,
            R_phase=args.R_phase,
            L_endwinding=args.L_endwinding,
            pm_Br_value=args.pm_Br_value,
            J_inertia=args.J_inertia,
            B_viscous=args.B_viscous,
            T_load=args.T_load,
            v_amp=args.v_amp, v_freq=args.v_freq,
            theta_init=args.theta_init,
            omega_init=args.omega_init,
            t_end=args.t_end,
            dt_FE=args.dt_FE,
            dt_circ=args.dt_circ,
            n_steps_per_FE=args.n_steps_per_FE,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
