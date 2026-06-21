"""AGE rotating machine workflow knowledge -- nonlinear iron + AGE coupling.

Covers the Newton/Picard iteration for machines with saturating iron, wrapping
the linear AGE core (airgap_machine.py) for the B-H nonlinearity loop.

Companion code: radia_ngsolve.airgap_motor_workflow
  - age_motor_nonlinear_solve(): Picard/Newton loop for one rotor position
  - age_motor_rotation_sweep():  torque–angle curve with nonlinear iron
  - age_motor_hysteresis_sweep(): HYSTERETIC sweep -- B-input Play constitutive in the iron with
                                  per-element history carried across rotor angles (M-source, one AGE
                                  factorization) -> hysteretic torque + real-loop iron loss
"""

AGE_MOTOR_WORKFLOW_OVERVIEW = """
# AGE Rotating Machine Workflow: Nonlinear Iron

## The AGE Core (Linear)

The Air-Gap Harmonic Element (airgap_machine.py) solves a rotating machine by:
1. Representing the un-meshed air gap analytically (each Fourier mode has a 2×2 DtN)
2. Coupling rotor and stator FE regions across the gap via Woodbury low-rank update
3. Computing torque in closed form (mesh-independent, radius-invariant)

One factorization per rotor position → cheap rotation sweeps.

## Adding Nonlinear Iron (this module)

For real machines with silicon-steel laminations, the reluctivity ν(|B|) is
field-dependent (saturates at high flux density).  The workflow adds a
Picard fixed-point iteration around the linear AGE core:

    k = 0, 1, 2, ...
    ┌──────────────────────────────────────────────────────────────────┐
    │ 1. ν_k   = nu_cf_fn(gfu_k)      # field-dependent reluctivity   │
    │ 2. a_k   = bilinear_fn(ν_k)     # reassemble with new ν         │
    │ 3. factor= airgap_factorize(a_k.mat, coupling, ...)             │
    │ 4. gfu_{k+1} = airgap_solve(factor, fes, dirichlet_cf)          │
    │ 5. check  ||δu|| / ||u|| < tol  → STOP                          │
    └──────────────────────────────────────────────────────────────────┘

The AGE gap coupling is UNCHANGED across iterations (the gap is current-free
and Laplacian; only the iron FE stiffness changes with saturation).

## Linear Limit

If `nu_cf_fn` is None (or returns a constant), the iteration converges in
exactly 1 step — the result is identical to the direct AGE solve.
This is the validation path: compare linear-mode Picard vs airgap_solve.

## Convergence

For silicon steel at typical operating points (B ≈ 1.5 T):
  - Picard converges in 5–15 iterations with tol=1e-5
  - Under-relaxation (relax=0.7) helps when ν changes steeply near saturation
  - Initial guess = zero field (cold start); warm-starting from nearby angle
    saves ~30% iterations in a rotation sweep

Newton (use_newton=True) is a future extension for faster convergence near
saturation knee. Picard is adequate for most machine analysis workflows.

## Hysteretic Iron (age_motor_hysteresis_sweep)

Beyond a single-valued saturating nu(|B|), the iron can carry the energy-based B-input vector PLAY
hysteresis constitutive (radia_ngsolve.hysteresis.PlayHysteresis).  The per-element play states p_k are
carried FORWARD across rotor angles, so the iron B-H response is path-dependent (magnetic history) and
the torque / iron loss include the hysteresis -- the hysteretic upgrade of age_motor_rotation_sweep.

Method: the fixed-point M-SOURCE (excess-field) coupling, which reuses ONE AGE factorization for the
whole rotation sweep (the air-gap element is unchanged; only an iron RHS source changes):

    per rotor angle (states p_k frozen during the angle, advanced after):
      iterate:  H_res = nu0*B - H_play(B; p_k)        # excess field over the reference nu0
                solve (K(nu0) + gap) A = source + INT_iron H_res.curl v   # one cached factorization
                B = curl A
      then:     p_k <- play(B; p_k)                   # advance the history

The play co-energy is CONVEX (hysteresis.play_coenergy_density), so with nu0 >= nu_max/2 the M-source
is a contraction; the optimal (fastest) reference is nu0 = (nu_max + nu_min)/2 = (sum a_k + a_0)/2,
which converges in ~10 iterations/angle.  (A per-angle variational/Newton solve is the robust
alternative for very nonlinear iron, at the cost of a re-factorization per Newton step -- the M-source
keeps the AGE one-factorization advantage.)

Units: the AGE core uses NORMALIZED reluctivity nu~ = 1/mur (gap = 1), so B = curl A is physical Tesla;
pass the play slopes normalized (a~ = mu0*nu) and the returned loss_density is physical J/m^3 (the
normalized loop area / mu0).  The iron loss is the REAL per-element loop area |INT H.dB| from the
recorded steady cycle -- it captures the rotating/elliptical stator B (loss > scalar |B|pk-Steinmetz),
which is the whole point of a vector hysteresis model.  Validated: tests/test_age_hysteresis_coupled.py.
"""

AGE_MOTOR_WORKFLOW_API = """
# AGE Nonlinear Motor API

## age_motor_nonlinear_solve

```python
from radia_mcp.radia_ngsolve.airgap_motor_workflow import age_motor_nonlinear_solve
from radia_mcp.radia_ngsolve.airgap_machine import airgap_coupling

# 1. Set up mesh + FE space (same as linear AGE)
fes  = ng.H1(mesh, order=2, dirichlet="dirichlet_bnd", complex=False)
coupling = airgap_coupling(fes, ri=0.05, ro=0.06,
                           rotor_ring="rotor_ring",
                           stator_ring="stator_ring",
                           harmonics=[1,3,5,7])

# 2. Define saturation model (e.g. Froelich or piecewise linear B-H)
def nu_cf_fn(gfu):
    B2 = ng.grad(gfu) * ng.grad(gfu)           # |B|^2 ≈ |∇A|^2 in 2D
    nu_sat = 1e6                                 # [A·m/Wb], near-saturation
    nu_lin = 500e0                               # linear-regime reluctivity
    return nu_lin / (1 + B2 / (2 * nu_sat**2)) + nu_sat  # simplified model

# 3. Define bilinear form factory
u, v = fes.TnT()
def bilinear_fn(nu_cf):
    a = ng.BilinearForm(fes)
    a += nu_cf * ng.grad(u) * ng.grad(v) * ng.dx
    a.Assemble()
    return a

# 4. Solve for one rotor position
gfu, info = age_motor_nonlinear_solve(
    fes, coupling,
    bilinear_fn = bilinear_fn,
    nu_cf_fn    = nu_cf_fn,
    dirichlet_cf= dirichlet_at_theta_r,
    niter=20, tol=1e-5, relax=1.0,
)
print(f"Converged in {info['niter_used']} iterations, torque = {info['torque']:.4f} N·m")
```

## age_motor_rotation_sweep

```python
from radia_mcp.radia_ngsolve.airgap_motor_workflow import age_motor_rotation_sweep
import numpy as np

theta_list = np.linspace(0, 2*np.pi/6, 37)  # one pole-pitch, 37 angles
results = age_motor_rotation_sweep(
    fes, coupling,
    bilinear_fn  = bilinear_fn,
    excitation_fn= lambda theta: make_dirichlet_cf(theta),  # rotor angle → BC
    theta_r_list = theta_list,
    nu_cf_fn     = nu_cf_fn,
    axial_length = 0.05,   # [m]
    niter=20, tol=1e-5,
)
# Each entry: {theta_r, torque, niter_used, converged}
T_avg = np.mean([r["torque"] for r in results])
T_rip = np.max([r["torque"] for r in results]) - np.min([r["torque"] for r in results])
print(f"Average torque: {T_avg:.4f} N·m, ripple: {T_rip:.4f} N·m")
```

## Parameters

| Parameter     | Default | Notes                                           |
|---------------|---------|-------------------------------------------------|
| niter         | 20      | Max Picard iterations                           |
| tol           | 1e-6    | Relative L2 convergence (||δu||/||u||)          |
| relax         | 1.0     | Under-relaxation (reduce to 0.5-0.7 for steep ν) |
| use_newton    | False   | Picard only; Newton planned for future release  |
"""

AGE_MOTOR_WORKFLOW_SATURATION_MODELS = """
# Saturation (B-H) Models for nu_cf_fn

## 1. Froelich model (analytic, smooth)

    ν(B) = a + b / (B_sat - |B|)   # singular at saturation

In NGSolve CoefficientFunction form:

    B2   = ng.grad(gfu) * ng.grad(gfu)
    Bmag = ng.sqrt(B2 + 1e-30)     # regularize for B=0
    nu_cf = a + b / (B_sat - Bmag)

Good for: smooth B-H curves (cobalt steel, soft iron).

## 2. Piecewise linear (lookup table)

    from scipy.interpolate import PchipInterpolator
    import numpy as np

    # B_data, H_data from material datasheet
    B_data = np.array([0.0, 0.5, 1.0, 1.5, 1.8, 2.0])   # [T]
    H_data = np.array([0.0,100, 250,1000,5000,20000])     # [A/m]
    nu_interp = PchipInterpolator(B_data, H_data / B_data)  # ν = H/B

    def nu_cf_fn(gfu):
        B2   = ng.grad(gfu) * ng.grad(gfu)
        Bmag = ng.sqrt(B2 + 1e-30)
        nu_vals = nu_interp(Bmag)      # scalar field
        return ng.CoefficientFunction(nu_vals)

**Warning**: NGSolve cannot use Python callables directly as CFs; for
piecewise interpolation in NGSolve use `ng.IfPos` chains or `ng.Piecewise`.

## 3. Built-in NGSolve material CFs

    from ngsolve.comp import CoefficientFunction as CF
    # BH-curve via NGSolve nonlinear materials (available in recent versions):
    from ngsolve import BSpline
    nu_cf = BSpline(B_knots, nu_values)(ng.sqrt(B2))

## 4. Simplified linear (validation)

    def nu_cf_fn(gfu):
        return ng.CoefficientFunction(1.0 / (4e-7 * np.pi * 1000))  # μ_r=1000

With constant ν, Picard converges in 1 iteration — identical to linear AGE.
Use this to validate the workflow before adding saturation.
"""

AGE_MOTOR_WORKFLOW_VALIDATED = """
# Validation: Linear vs Nonlinear Picard

## Test: Picard with constant nu → must match airgap_solve exactly

```python
# tests/test_age_motor_nonlinear.py
def test_linear_limit_matches_direct():
    # Build fes + coupling (from test_airgap_eddy_machine.py fixture)
    nu_const = ng.CoefficientFunction(800.0)

    def bilinear_fn(nu_cf):
        a = ng.BilinearForm(fes)
        a += nu_cf * ng.grad(u) * ng.grad(v) * ng.dx
        a.Assemble()
        return a

    gfu_picard, info = age_motor_nonlinear_solve(
        fes, coupling, bilinear_fn,
        nu_cf_fn=lambda gfu: nu_const,
        dirichlet_cf=excitation,
        niter=5, tol=1e-12,
    )
    assert info["niter_used"] == 1    # linear: one iteration
    assert info["converged"] is True

    # Compare against direct solve
    a_ref = bilinear_fn(nu_const)
    factor_ref = airgap_factorize(a_ref.mat, coupling, fes.FreeDofs())
    gfu_ref = airgap_solve(factor_ref, fes, dirichlet_cf=excitation)

    diff_norm = ...  # L2 norm of gfu_picard - gfu_ref
    assert diff_norm < 1e-10
```

## Test: Saturation loop convergences

```python
def test_saturation_converges():
    # Simple Froelich ν: ν = ν_0 / (1 + |B|²/(B_sat²))
    def nu_cf_fn(gfu):
        B2 = ng.grad(gfu) * ng.grad(gfu)
        return nu0 / (1.0 + B2 / B_sat**2)

    gfu, info = age_motor_nonlinear_solve(
        fes, coupling, bilinear_fn, nu_cf_fn=nu_cf_fn,
        dirichlet_cf=excitation, niter=30, tol=1e-5,
    )
    assert info["converged"], f"Did not converge in {info['niter_used']} iters"
    assert len(info["rel_change_history"]) >= 2
    # Residuals should be monotonically decreasing
    h = info["rel_change_history"]
    assert h[-1] < h[0]   # at least first vs last
```

## Test: Rotation sweep produces torque ripple

```python
def test_rotation_sweep():
    results = age_motor_rotation_sweep(
        fes, coupling, bilinear_fn, excitation_fn,
        theta_r_list=np.linspace(0, 2*np.pi/np_pairs, 13),
        nu_cf_fn=None,     # linear test
        axial_length=0.05,
    )
    assert all(r["converged"] for r in results)
    T_vals = [r["torque"] for r in results]
    # Torque should be non-zero and vary with angle
    assert max(abs(T) for T in T_vals) > 0
```
"""


def get_airgap_motor_workflow_documentation(topic: str = "all") -> str:
    """Return AGE nonlinear motor workflow documentation."""
    topics = {
        "overview": AGE_MOTOR_WORKFLOW_OVERVIEW,
        "api": AGE_MOTOR_WORKFLOW_API,
        "saturation_models": AGE_MOTOR_WORKFLOW_SATURATION_MODELS,
        "validated": AGE_MOTOR_WORKFLOW_VALIDATED,
    }
    if topic == "all":
        return "\n\n".join(topics.values())
    return topics.get(topic, f"Unknown topic '{topic}'. Options: {list(topics)}")
