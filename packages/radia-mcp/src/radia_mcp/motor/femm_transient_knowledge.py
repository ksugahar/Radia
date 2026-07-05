"""FEMM transient motor solver — Lange-Henrotte-Hameyer 2009 method.

Source documents:
- https://www.femm.info/doku/doku.php?id=newbuild — David Meeker's
  "new build" page describing the transient FE-circuit coupling shipped
  in FEMM 4.2 (build dates 2018-2023).
- E. Lange, F. Henrotte, K. Hameyer, "An Efficient Field-Circuit Coupling
  Based on a Temporary Linearization of FE Electrical Machine Models",
  IEEE Trans. Mag. 45(3):1258-1261, 2009.
  DOI: 10.1109/TMAG.2009.2012585.
  (Institute of Electrical Machines, RWTH Aachen University.)
- Dave Meeker, private note (forwarded by Kengo Sugahara,
  2026-05-21).

The middle author **François Henrotte** is the same Henrotte cited in
`radia_mcp.differential_forms.forces_knowledge` for the 2004 EM Force
Density and the Bossavit-Henrotte 2004 Maxwell-stress Handbook —
this 2009 paper is the same author's machine-analysis contribution.

This is the **3rd transient-machine-FE approach** captured in radia_mcp.
The other two are:

| Approach | Subpackage / module | Strength |
|----------|---------------------|----------|
| ONELAB / GetDP magnetodynamic A (moving band) | `onelab_knowledge` | Industrial open-source reference |
| Kaimori-Wakao Darwin TD (Coulomb gauge) | `darwin_model_knowledge` | Capacitive coupling, PWM-EMI |
| **Lange-Henrotte-Hameyer 2009 (FEMM)** ⭐ | this module | Hybrid FE + circuit, fast per-step — **Sugahara-lab recommended path for transient motor analysis** |

The Lange-Henrotte-Hameyer approach is distinguished by **factoring
the nonlinear FE solve from the time integration**: solve the
nonlinear FE problem **once** at a chosen operating point, extract a
**locally linearized** state-space `(L_inc(i, θ), e_bemf(i, θ, ω),
T_em)`, then let the circuit simulator integrate the linearized ODE
with a *much smaller* electrical/PWM time step **without re-running
the FE** until the linearization becomes stale.

**The "no rotor rotation" insight** (Henrotte's core contribution):
between two FE evaluations, the rotor angle θ(t) advances by the
mechanical ODE `J dω/dt = T_em - T_load`, integrated **analytically**
in the lumped state space — there is NO mesh motion, NO remesh, and
no rotor coordinate change inside the FE.  The FE solver only sees
the rotor at the **discrete sequence of points** where the
linearization is refreshed, which is *one or two orders of magnitude*
fewer time steps than the circuit-level resolution.
"""

from __future__ import annotations


LANGE_2009 = """\
## Lange-Henrotte-Hameyer 2009 — Temporary linearization of FE machine models

Reference: E. Lange, **F. Henrotte**, K. Hameyer, "An Efficient
Field-Circuit Coupling Based on a Temporary Linearization of FE
Electrical Machine Models", IEEE Trans. Mag. 45(3):1258-1261, 2009.
DOI: 10.1109/TMAG.2009.2012585.

Institute of Electrical Machines, RWTH Aachen University.

### The setup

In a rotating machine with N_phase phases at instantaneous current
vector `i(t)` and rotor angle `θ(t) = ∫ω dt`, the flux linkage is

  ψ(i, θ)  =  ∫_Ω (∇ × A) · n_coil  dΩ

with A the magnetic vector potential and `n_coil` the per-phase turn
density.  Both `A` and `ψ` are nonlinear functions of `i` (through
B-H saturation) and `θ` (through geometry rotation).

The circuit equation in voltage-driven form:

  v(t)  =  R · i  +  dψ/dt
        =  R · i  +  ∂ψ/∂i · di/dt  +  ∂ψ/∂θ · ω
        =  R · i  +  L_inc(i, θ) · di/dt  +  e_bemf(i, θ, ω)

where:
- L_inc = **incremental inductance matrix** (N_phase × N_phase),
  `L_inc[j,k] = ∂ψ_j/∂i_k` at the current operating point.
- e_bemf = **back-EMF vector**, `e_bemf[j] = ∂ψ_j/∂θ · ω`.

Both depend on (i, θ), but **locally** (for a small δt around the
present operating point) they can be treated as constants.

### The trick

At each time step:

1. **Operating-point solve**: full nonlinear FE solve at the present
   `(i, θ)` to obtain B(x), H(x), reluctivity ν(B) per element.
2. **Incremental permeability**: freeze the reluctivity ν*(x) =
   dB/dH (the *incremental* slope of the BH curve at the current B).
   This is **linear** in the new variable δA — see Kaimori's
   "incremental permeability" derivation, identical to what
   FEMM and JMAG call "frozen permeability".
3. **N_phase linear solves**: for each phase j, set i_j = 1, all
   others = 0; solve the linearized FE problem with ν*(x).  Extract
   ψ_k(j) for k = 1..N_phase.  This row of ψ becomes column j of
   L_inc.
4. **Back-EMF**: Compute virtual-work torque at the operating point;
   then `∂ψ/∂θ · ω` is recovered from the torque using `T_em = i^T · ∂ψ/∂θ` (Park identity).
5. **Step circuit**: integrate

     L_inc · di/dt + R · i = v - e_bemf

   from t to t+Δt by a quadrature scheme (Backward Euler, Crank–
   Nicolson, ...) to get `i(t+Δt)`.
6. **Update rotor**: `θ(t+Δt) = θ(t) + ω · Δt`; update mechanical
   ODE `J dω/dt = T_em - T_load - B_visc · ω`.

### Cost per step

- 1 nonlinear FE solve + N_phase linear FE solves with the same
  Jacobian matrix → factorize ONCE, back-substitute N_phase times.
- Compared to a fully-coupled nonlinear FE+circuit step (ONELAB
  style), this is roughly (1 + 0.1·N_phase) FE-solve-equivalents
  per step vs ~3-10 for fully coupled (Newton iterations).

### Accuracy

The linearization is **exact** in the limit Δt → 0; the truncation
error scales as O(Δt²) for Crank-Nicolson, the same as fully-coupled
Newton.  In practice Lange showed step sizes up to ~T_pwm/20 give
indistinguishable results from fully-coupled.

### What it does NOT capture

Lange-2009 explicitly **drops the eddy-current term** in the FE
domain (σ ∂A/∂t = 0).  If the rotor cage / lamination eddy currents
matter, you must add them back in via **homogenized lamination
models** (Bertotti or Pry-Bean) at the post-processing stage; or you
must combine Lange with the standard ONELAB σ-coupled A-formulation
(more expensive).
"""

FEMM_NEWBUILD = """\
## FEMM "newbuild" — Dave Meeker's implementation

URL: https://www.femm.info/doku/doku.php?id=newbuild

FEMM (Finite Element Method Magnetics, David Meeker) ships the
Lange-2009 algorithm as the **transient motor solver** in its 2018+
builds.  Latest build: 2023-10-22 with VS2022 compatibility and
radiation/convection BC updates.

### Reference benchmark — `antunes.fem`

A 2-hp, 2000-rpm, 8-pole BLAC machine (Antunes-Bastos-Sadowski) ships
with FEMM as the canonical transient demo:

| Parameter | Value |
|-----------|-------|
| Rated power | 2 hp (1.5 kW) |
| Rated speed | 2000 rpm |
| Pole pairs  | 4 |
| Phase resistance | 0.3872 Ω |
| Phase inductance | 4.4 mH |
| PM flux linkage  | 0.1157 Wb |

Sliding-band BC at the air-gap mid-circle (NO remeshing as rotor
rotates) — Meeker's distinguishing feature vs ONELAB / JMAG which
both use a moving-band layer.

### Simulink S-function coupling

The shipped pattern is:

```
Simulink ODE solver (ode23t, ode45)
   ↓ v(t), i(t-) requests
S-function `femm_step.m`
   ↓ FEMM-Lua via OctaveFEMM or pyFEMM
FEMM transient analysis
   ↓ returns di/dt at present (i, θ, ω)
back to Simulink → integrates one step
```

The `mdlOutputs` callback in the S-function performs **exactly one**
nonlinear solve + N_phase incremental solves per requested di/dt
evaluation.  This is repeatable across machine geometries by editing
only the `.fem` file (FEMM problem) — the S-function never changes.

### Sliding-band BC

Implementation: introduce a **sliding-mesh boundary condition** at
the air-gap centerline.  The stator and rotor have separate meshes
that meet on a circle.  The matching is done by:

1. Project the rotor-side trace of A onto a Fourier basis around the
   air-gap circle, ψ_rotor(θ_local).
2. Project the stator-side trace, ψ_stator(θ_stator).
3. Equate them with a rotation: `ψ_rotor(θ_local) = ψ_stator(θ_stator
   - θ_rotor(t))`.
4. The Fourier truncation is treated as a (small) modeling error.

**Why this matters**: in moving-band approaches (ONELAB, JMAG) the
mesh changes every time step (the thin annulus is remeshed), which
introduces small **discretization noise** in the torque waveform —
visible as high-frequency ripple synchronized with the moving-band
remesh step.  Sliding-band has no such noise; the mesh is identical
across all time steps.

The cost is that the air-gap force/torque integration must use the
**stator-side trace of B** plus the projected rotor-side **B from
Fourier reconstruction** — slightly more complex than a single
volumetric integral.  Meeker chose Arkkio-formula style integration
on the stator side of the band to avoid the projection completely
on the force side.

### What this means for radia_mcp.motor

The `femm_transient_knowledge` topic gives the **fastest** transient
approach available for nonlinear FE machine simulation.  Apply when:

- You have a 2D machine geometry (FEMM is 2D only).
- You need transient torque ripple, harmonics, faulted operation.
- You DON'T need detailed cage eddy currents (use ONELAB instead).
- You want Simulink integration (controller-in-the-loop).

To replicate in radia_mcp.radia_ngsolve: implement an
`IncrementalPermeabilityFE` driver that

1. Holds an NGSolve `BilinearForm` for the nonlinear ν(B) problem.
2. After each nonlinear solve, freezes ν*(x) = dB/dH(B_op) and
   builds the linearized BilinearForm.
3. Solves N_phase linear problems (PARDISO with **shared factorization**
   — NGSolve's `Inverse(symmetric=True, freedofs=...)` caches the
   factor; just call it N_phase times with different RHS).
4. Computes ψ_jk = ∫ (∇×A_j) · n_coil_k dΩ → L_inc.
5. Computes T_em via Arkkio on the air-gap circle → e_bemf.
6. Hands (L_inc, R, e_bemf) to a SciPy/scikit-ODE integrator.
"""

SLIDING_BAND_VS_MOVING_BAND = """\
## Sliding-band vs Moving-band — practical comparison

Both techniques solve the same problem: how to couple a rotating
rotor mesh to a stationary stator mesh without remeshing the whole
machine at every time step.

| Aspect | Moving band (ONELAB) | Sliding band (FEMM) |
|--------|----------------------|---------------------|
| Stator mesh | Fixed | Fixed |
| Rotor mesh | Fixed | Fixed |
| Air-gap annulus | Thin layer remeshed each step | Single circle, no mesh |
| Air-gap DOFs | All elements in the annulus | Fourier modes on the circle |
| Implementation | Mesh editing primitive | BC primitive |
| Numerical noise | High-frequency ripple from remesh | None |
| Force integration | Volumetric in annulus | Stator-side surface |
| Torque ripple accuracy | Limited by annulus mesh density | Limited by Fourier truncation N_f |

For an 8-pole machine, N_f ≈ 100 Fourier modes is plenty (covers
harmonics up to ~12·8 = 96th); the projection cost is O(N_arc ·
N_f) per step, negligible vs the FE solve.

### Why moving-band is still the default in ONELAB

- 3D extension: moving-band generalizes to 3D (cylindrical-shell
  remesh of an annular ring); sliding-band 3D requires spherical-
  harmonic-style projection that is more complex.
- Multi-physics: if the air gap also carries heat or stress that
  matters (cooling, fan), the annulus needs a mesh for those, so
  you may as well use moving-band for E&M too.
- Slotted air gap: if the rotor has surface features that protrude
  into the air gap, sliding-band breaks (the air-gap circle is not
  well-defined); moving-band handles it.

For accelerator magnets (e.g. wiggler, undulator) with **rotational
periodicity but no actual rotation**, neither technique is needed
— just use periodic boundary conditions on the stator/rotor sector.
"""

NGSOLVE_RECIPE = """\
## NGSolve recipe for Lange-2009 transient

Skeleton (sketch — production code would live in
`radia_mcp.radia_ngsolve.motor_transient`):

```python
import numpy as np
from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                     CoefficientFunction, grad, x, y, dx, ds, Integrate,
                     TaskManager, PARDISO)
from radia_mcp.motor.onelab_knowledge import get_onelab_knowledge

# 1. Load 2D machine mesh (Cubit / Gmsh / Netgen-OCC)
mesh = Mesh("antunes.vol")
NbrPhases = 3

# 2. Define nonlinear ν(B) coefficient
def nu_iron(B2):
    # Spline of vacoflux50 or M19_29G — replace with material table
    Bsat = 1.8
    return 1.0 / 4e-7 / 3.14159 / (1 + 999/(1 + B2/Bsat**2))

# 3. Function space (2D scalar A_z formulation)
fes = H1(mesh, order=2, dirichlet="surface_inf")
A = GridFunction(fes)   # operating-point A

# 4. Nonlinear operating-point solve
def solve_nonlinear(i_phase, theta_rotor, A_init):
    a = BilinearForm(fes, symmetric=True)
    # ... assemble with ν(|grad A|²) using i_phase as current source
    # Newton with line search, return converged A
    return A_op

# 5. Incremental permeability matrix builder
def build_L_inc(A_op, theta_rotor):
    # Freeze ν*(x) = dB/dH at A_op
    B2_op = (grad(A_op)[0]**2 + grad(A_op)[1]**2)
    nu_star = CoefficientFunction(...)  # numerical d(nu)/dB^2

    a_lin = BilinearForm(fes, symmetric=True)
    a_lin += nu_star * grad(u) * grad(v) * dx
    a_lin.Assemble()

    inv = a_lin.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
    # ↑ factorize ONCE; back-substitute N_phase times.

    L_inc = np.zeros((NbrPhases, NbrPhases))
    e_bemf = np.zeros(NbrPhases)

    for j in range(NbrPhases):
        f = LinearForm(fes)
        f += 1.0 * v * dx("phase_%d" % j)   # unit current in phase j
        f.Assemble()

        gfu_j = GridFunction(fes)
        gfu_j.vec.data = inv * f.vec

        for k in range(NbrPhases):
            psi_jk = Integrate(
                grad(gfu_j)[1] * 1.0,   # 2D: A is A_z, ψ = ∫ ∂A/∂y dΩ
                mesh, definedon=mesh.Materials("phase_%d" % k))
            L_inc[k, j] = psi_jk

        # Back-EMF: dψ/dθ from finite difference (or virtual work)
        # e_bemf[j] = ω · dψ_j/dθ

    return L_inc, e_bemf

# 6. Time integrator (Crank-Nicolson)
def step_circuit(i_n, L_inc, R, v_n, v_np1, e_n, e_np1, dt):
    # (L_inc + dt/2 R) i_{n+1} = (L_inc - dt/2 R) i_n
    #                            + dt/2 (v_n + v_{n+1} - e_n - e_{n+1})
    A_mat = L_inc + 0.5 * dt * R
    rhs = (L_inc - 0.5 * dt * R) @ i_n + \\
          0.5 * dt * (v_n + v_np1 - e_n - e_np1)
    return np.linalg.solve(A_mat, rhs)

# 7. Main time loop
t = 0.0; dt = 1e-5
i = np.zeros(NbrPhases); theta = 0.0; omega = 2 * np.pi * 33.3
while t < 0.1:
    A_op = solve_nonlinear(i, theta, A_init=A.vec.FV().NumPy())
    A.vec.FV().NumPy()[:] = A_op.vec.FV().NumPy()
    L_inc, e_bemf = build_L_inc(A_op, theta)
    R = 0.3872 * np.eye(NbrPhases)
    v = v_source(t)
    i = step_circuit(i, L_inc, R, v, v, e_bemf, e_bemf, dt)
    theta += omega * dt
    t += dt
```

Open issues to resolve in a production implementation:
- ψ definition in 2D is per-meter; multiply by stack length L_z to
  get actual flux linkage.
- Per-phase coil regions need explicit + and - turn sub-regions
  (Stator_Ind_Ap / Stator_Ind_An) to get sign-correct ψ integration.
- PM contribution to ψ enters as a constant offset; compute it once
  at theta = 0 and rotate.
- Sliding-band BC: implement as a Lagrange-multiplier constraint
  in a second FESpace and assemble the coupled saddle-point system.
  Alternatively, use moving-band remesh in the air-gap annulus and
  accept the ripple.
"""

SIMULINK_COUPLING = """\
## Simulink coupling (controller-in-the-loop)

For controller verification (PI/MPC current control, vector control,
sensorless observers), the Lange-2009 method shines because it
exposes a clean small-signal model `L_inc · di/dt = v - R·i - e_bemf`
that Simulink can integrate.

### S-function pattern (MATLAB)

```matlab
function [sys, x0, str, ts] = femm_motor_sfun(t, x, u, flag)
% u = [v_a; v_b; v_c]   inputs from inverter
% x = [i_a; i_b; i_c; theta; omega]   states
switch flag
case 0  % mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 5;
    sizes.NumOutputs     = 5;
    sizes.NumInputs      = 3;
    sys = simsizes(sizes);
    x0 = [0;0;0; 0;0];
    str = []; ts = [0 0];
case 1  % mdlDerivatives
    i_phase = x(1:3);
    theta   = x(4);
    omega   = x(5);
    v       = u;
    [L_inc, e_bemf, T_em] = femm_op_point(i_phase, theta);
    R = 0.3872 * eye(3);
    di_dt = L_inc \\ (v - R*i_phase - e_bemf);
    dtheta_dt = omega;
    domega_dt = (T_em - T_load) / J;
    sys = [di_dt; dtheta_dt; domega_dt];
case 3  % mdlOutputs
    sys = x;
end
```

with `femm_op_point` calling FEMM (or NGSolve via Python) through
the OctaveFEMM wrapper.

For Python-only stacks, the equivalent uses `scipy.integrate.solve_ivp`
with `rhs=lambda t, x: motor_ode(t, x, ngsolve_op_point)`.

### Real-time considerations

A single op-point evaluation in 2D FEMM takes ~50-200 ms (depends on
mesh).  This is **not real-time** for inverter control loops (kHz
rate), but is well-suited for:

- Offline controller tuning (10x to 100x real-time slower OK).
- HIL setup where the motor model runs on a separate compute card.
- System-level verification before deployment.

For real-time, build a polynomial-fit or LUT model of `L_inc(i, θ)`,
`e_bemf(i, θ, ω)`, `T_em(i, θ)` from FE sweeps, then evaluate the
LUT at runtime.  This is the standard "FEA → mapping → real-time"
workflow (e.g. JMAG-RT, ANSYS Motor-CAD).
"""

HENROTTE_NO_ROTATION = """\
## Henrotte's "no rotor rotation" insight — why this method is clever

In Lange-Henrotte-Hameyer 2009 (and David Meeker's FEMM
implementation), the FE simulation **never sees the rotor in motion**
— it only sees the rotor at a sparse sequence of *frozen* angles.
This is structurally different from ONELAB / JMAG / standard FE-with-
moving-band, where the rotor mesh is updated every electrical-
resolution time step.

### Two-time-scale decoupling

There are two natural time scales in an inverter-driven machine:

| Scale | What changes here | Typical Δt |
|-------|-------------------|-----------|
| Mechanical (FE) | rotor angle θ, BH operating point, saturation | 100 μs – 1 ms |
| Electrical / PWM (circuit) | phase currents, switch events, di/dt | 1 – 10 μs |

ONELAB-style fully coupled FE+circuit must use the smaller Δt
*everywhere* — including in the FE solve where it is wasteful, because
saturation barely changes between PWM ticks.  Lange-Henrotte-Hameyer
decouple these:

- **FE time step Δt_FE**: chosen by the mechanical / saturation
  time scale.  At each Δt_FE: 1 nonlinear solve + N_phase incremental
  solves → refresh `(L_inc, e_bemf, T_em)`.
- **Circuit time step Δt_circ**: chosen freely by the circuit
  simulator (Simulink ode23t, ODE45, ...).  Between FE refreshes the
  rotor angle advances analytically as `θ(t+τ) = θ(t) + ω·τ + (τ²/2)·
  (T_em - T_load)/J`.

The "FE never rotates" property is consistent with the paper's claim:
*"This method allows decoupling the time constants of the field
problem from that of the circuit problem, which is typically one or
two orders of magnitude smaller."*  (Lange-Henrotte-Hameyer 2009,
introduction.)

### Why the rotor doesn't need to physically move in the FE

Even when Δt_FE is reached and a new FE solve is triggered, the rotor
is positioned at the **current** angle θ(t) by **rotating the source
distribution** (PM remanence direction, stator winding phase angles
relative to rotor d-axis), NOT by rotating the mesh.  In a 2D
axisymmetric mesh:

- Stator iron / coils: fixed
- Rotor iron / PM regions: fixed mesh, fixed material distribution
- PM remanence vector `B_r(θ_rotor)` rotated analytically in space
- Stator current pattern phased to `(θ_rotor + δ)` where δ is the
  torque angle

The result: **the same mesh is used at every Δt_FE** — only the
material function `B_r(x; θ)` and the source `J_s(x; θ)` change.  No
remeshing, no moving-band layer, no sliding-band BC (though FEMM does
use sliding-band for the trace-of-B coupling on the air-gap circle —
see `sliding_vs_moving_band` topic).

### Why this is structurally cleaner than ONELAB/JMAG

| Property | ONELAB (moving band) | FEMM-Lange (no rotation) |
|----------|----------------------|--------------------------|
| Mesh per step | Re-meshed in annulus | Identical mesh always |
| Discretization error | New per step | Time-constant per refresh interval |
| Time-step coupling | All physics share Δt | Decoupled by 1-2 orders of magnitude |
| Implementation | Mesh primitive (NGSolve `Identify`) | BC primitive (sliding band) + analytical PM rotation |
| Best for | Detailed eddy currents, 3D | Inverter-driven control / HIL |

### Energy-based linearization (Henrotte hallmark)

Section IV of the paper develops the linearization from an **energy
flow** argument: write the magnetic energy as a function of currents
`W(i, θ)`, expand to second order around (i₀, θ₀), identify L_inc
as the Hessian `∂²W/∂i² |_(i₀, θ₀)` and e_bemf from `∂W/∂θ`.  This
gives a *consistent* mapping between FE vector-potential variables
and the lumped circuit (i, ψ) variables — energy-conservative, exactly
the property Henrotte made foundational in his Bossavit-Henrotte 2004
Handbook.

The consistency means that **torque from the field calculation =
mechanical power supplied to the shaft = electrical power supplied
to the winding − resistive losses**, exactly, in the lumped model.
No spurious torque ripple from inconsistent torque extraction (a
common artifact in naive coupling schemes).
"""

LAB_RECOMMENDATION = """\
## Sugahara-lab recommendation for transient motor analysis

For new transient motor work in the Sugahara lab (radia + NGSolve
stack), **Lange-Henrotte-Hameyer 2009 is the recommended path**.
Rationale:

1. **Decouples FE from circuit time scales** — the actual physical
   structure of the problem.  Fully-coupled ONELAB-style integration
   wastes effort using PWM-resolution time steps in the FE solver.
2. **No moving mesh needed** — the FE mesh is identical at every
   refresh; only source distributions (PM `B_r(θ)`, stator current
   phases) rotate.  No mesh-noise contamination of torque waveforms.
3. **Same FE infrastructure as static / cogging studies** — the
   nonlinear FE solver and the incremental-permeability linearizer
   are exactly what you would write for cogging-torque maps or
   d/q inductance maps.  No motor-specific FE primitives.
4. **Clean Simulink / scipy-ODE hand-off** — the linearized
   `L_inc, e_bemf, T_em` triplet is a state-space block any
   controller-design framework can drive.  Liu Xinyao thesis-style
   SynRM control studies become tractable.
5. **Energy-consistent** (Henrotte hallmark) — torque from the field
   computation equals shaft power, by construction.  No spurious
   ripple from inconsistent extraction.

### Hameyer / Henrotte canonical references (Sugahara lab library)

All available at public-safe curated corpus
or public-safe curated corpus

| Year | Authors | Title | Relevance |
|------|---------|-------|-----------|
| 1993 | Henrotte-Hedia-Bamps-Genon-Nicolet-Legros | A New Method for Axisymmetrical Linear and Nonlinear Problems | Foundation of Henrotte basis (cf. `radia.axifem`) |
| 1997 | Dular-Henrotte-Robert-Genon-Legros | A Generalized Source Magnetic Field Calculation Method for Inductors of any Shape | Stage-3 source field construction (Whitney) |
| 2004 | Bossavit-Henrotte | EM Force Density in Ferromagnetic Material | `differential_forms.forces_knowledge` |
| 2006 | Henrotte | Energy-based vector hysteresis model | Energy-method core paper |
| 2009 | Lange-Henrotte-Hameyer | An Efficient Field-Circuit Coupling … | ⭐ this module |
| 2013 | Francois-Lavet-Henrotte-Stainier-Noels-Geuzaine | Energy-based variational model of ferromagnetic hysteresis for FE | Variational extension |
| ?    | Hameyer et al. | Eddy Currents in Windings of Switched Reluctance Machines (187p) | SRM cage eddy currents (monograph-style) |
| ?    | Henrotte (RWTH lecture) | Energy-Based Magnetic Hysteresis Models (254p) | Comprehensive monograph |

The arc is consistent: **energy-based formulations of E&M from FE
discretization through hysteresis through field-circuit coupling**.
Lange-Henrotte-Hameyer 2009 is the field-circuit-coupling capstone of
this research program — built on top of the energy-based hysteresis
work (2006/2013) and the Whitney/source-field machinery (1993/1997).

### Implementation roadmap (Stage-2 → Stage-3 per Panel Design Policy)

**Stage 1 — Enumerate** (one paragraph):
- knobs: mesh file, BH-curve material, N_phase, R, J_inertia,
  v_source(t), T_load(t), tspan, dt_FE, dt_circ
- solver switch: `--method lange-henrotte | onelab-coupled | darwin-td`

**Stage 2 — CLI script** (`calc_motor_transient.py` under
`src/radia/panels/`):
- argparse with above knobs
- Read `.vol` mesh
- Build nonlinear `BilinearForm` with ν(B)
- Outer loop: FE refresh at dt_FE intervals
- Inner loop: scipy.integrate.solve_ivp on linearized state space
- JSON stdout: i(t), θ(t), ω(t), T(t), losses
- Golden test: antunes-equivalent 8-pole BLAC

**Stage 3 — Panel widget**: `radia_motor.py` wrapping Stage-2 script.

This module is the **knowledge prerequisite** for Stage 2.
"""

SECTIONS = {
    "lab_recommendation": LAB_RECOMMENDATION,
    "lange_2009": LANGE_2009,
    "henrotte_no_rotation": HENROTTE_NO_ROTATION,
    "femm_newbuild": FEMM_NEWBUILD,
    "sliding_vs_moving_band": SLIDING_BAND_VS_MOVING_BAND,
    "ngsolve_recipe": NGSOLVE_RECIPE,
    "simulink_coupling": SIMULINK_COUPLING,
}


def get_femm_transient_knowledge(topic: str = "lange_2009") -> str:
    """Return FEMM transient-motor / Lange-2009 knowledge.

    Args:
        topic: section name; "all" returns everything.
    """
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
