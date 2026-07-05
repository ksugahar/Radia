"""Data assimilation methods for EM state estimation + sensor fusion."""


# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview": "(see knowledge file for details)",
    "kalman_family": "(see knowledge file for details)",
    "four_dvar": "(see knowledge file for details)",
    "all": "return \"\n\n\".join([OVERVIEW, KALMAN_FAMILY, FOUR_DVAR])",
}

OVERVIEW = r"""
# Data Assimilation Overview

Combine a **physical model** (forward simulation) with **noisy
sensor measurements** to estimate the system state better than
either alone.

The Bayesian view:

    p(state | observations) ∝ p(observations | state) * p(state | model)

DA methods differ in:
- How they represent the posterior (Gaussian, ensemble, particle)
- Whether they handle nonlinearity (EKF, UKF, EnKF, PF)
- Whether they go forward in time only (filter) or smooth over a
  window (smoother, 4D-Var)

## Decision table for EM problems

| Use case | Method |
|----------|--------|
| Linear FEM forward + Gaussian noise | **Kalman filter** |
| Nonlinear forward, small system (< 100 state) | **EKF** or **UKF** |
| Large state (FEM with > 1000 DOF) | **Ensemble Kalman Filter (EnKF)** |
| Strongly non-Gaussian / multi-modal | **Particle filter** |
| Smooth over a fixed time window | **4D-Var** |

## DA method family

```
Data Assimilation
├── Filtering (sequential, online)
│   ├── Kalman Filter (linear Gaussian)
│   ├── Extended KF (linearize around current estimate)
│   ├── Unscented KF (deterministic sigma points)
│   ├── Ensemble KF (random ensemble, Monte Carlo)
│   └── Particle Filter (sequential Monte Carlo)
├── Smoothing (window-based, can be offline)
│   ├── RTS smoother (Kalman backward pass)
│   ├── 4D-Var (variational over time window)
│   └── Ensemble smoother (EnKS)
└── Parameter estimation
    ├── Joint state-param filter (augment state)
    ├── Iterative-EnKF
```

## EM applications

- **Magnetic source tracking**: estimate position/strength of moving
  permanent magnet from external B sensors (Kalman or particle)
- **Eddy current state in inductor**: estimate distributed current
  in workpiece from coil voltage + position measurements (EnKF)
- **Hysteresis state tracking**: Play model state has K hysteron
  variables, noisy B sensor → EKF over Play state
- **BH curve online identification**: data assimilation + parameter
  estimation hybrid
- **Plasma equilibrium reconstruction** (fusion): magnetic
  diagnostics + Grad-Shafranov forward → EnKF (cross-link:
  `radia_mcp.fusion_reactor`)
"""


KALMAN_FAMILY = r"""
# Kalman Filter Family (for EM state estimation)

## Linear Kalman Filter (KF)

For a system:
    x_t = F * x_{t-1} + B * u_t + w_t,   w ~ N(0, Q)  (process noise)
    y_t = H * x_t + v_t,                  v ~ N(0, R)  (obs noise)

Two-step recursion:

    Predict:
        x_pred = F * x_prev + B * u
        P_pred = F * P_prev * F^T + Q

    Update (when y observed):
        K = P_pred * H^T * (H * P_pred * H^T + R)^-1   (Kalman gain)
        x = x_pred + K * (y - H * x_pred)
        P = (I - K * H) * P_pred

Optimal for linear Gaussian; closed-form, no iteration.

## EM example: 1D magnet position tracking

State: x = (position, velocity) of a moving magnet.
Observation: B field at fixed sensor (linear in x for small motion).

```python
import numpy as np

F = np.array([[1, dt], [0, 1]])           # constant-velocity model
B = np.array([[0.5*dt**2], [dt]])         # acceleration input
H = np.array([[dB_dx, 0]])                # linearized B(x)
Q = sigma_a**2 * np.array([[dt**4/4, dt**3/2], [dt**3/2, dt**2]])
R = np.array([[sensor_var]])

x_est = np.zeros(2)
P_est = np.eye(2) * 1e6

for t, y_obs in zip(times, B_measurements):
    # Predict
    x_pred = F @ x_est + B.flatten() * a_known
    P_pred = F @ P_est @ F.T + Q
    # Update
    S = H @ P_pred @ H.T + R
    K = P_pred @ H.T @ np.linalg.inv(S)
    x_est = x_pred + (K @ (y_obs - H @ x_pred)).flatten()
    P_est = (np.eye(2) - K @ H) @ P_pred
```

## Extended KF (EKF) for nonlinear EM

When the forward `y = h(x)` is nonlinear (full Biot-Savart, FEM):

    K = P_pred * H^T * (H * P_pred * H^T + R)^-1
        with H = ∂h/∂x evaluated at x_pred (Jacobian)

Caveats:
- Linearization can diverge if x_pred is far from truth
- Jacobian via autograd (PyTorch / JAX) or finite difference
- For full FEM forward, computing H is expensive → consider UKF instead

## Unscented KF (UKF)

Sample 2N+1 "sigma points" around x_est, propagate each through the
nonlinear h, fit a Gaussian to the propagated cloud. Avoids Jacobian.

Better than EKF for moderately nonlinear; same cost as ~5x KF.

```python
from filterpy.kalman import UnscentedKalmanFilter, MerweScaledSigmaPoints

points = MerweScaledSigmaPoints(n=state_dim, alpha=.1, beta=2., kappa=-1)
ukf = UnscentedKalmanFilter(dim_x=state_dim, dim_z=obs_dim, dt=dt,
                              fx=process_model, hx=observation_model,
                              points=points)
ukf.x = initial_state
ukf.P = initial_covariance
for y_t in observations:
    ukf.predict()
    ukf.update(y_t)
```

## Ensemble KF (EnKF) for large state

For high-dim state (FEM with 10^4 DOF), full KF is impossible
(NxN covariance). EnKF uses N ensemble members instead:

    {x_i^pred}_{i=1..N} ← propagate ensemble through model
    P_pred ≈ sample covariance of ensemble
    For each member, perturb obs by sensor noise; apply Kalman update

Cost: O(N * forward_solve) per step. Use N ~ 50-200.

```python
import numpy as np

def enkf_step(ensemble, y_obs, H, R, forward_model):
    N = ensemble.shape[1]
    # Predict each member
    ensemble_pred = np.column_stack([forward_model(ensemble[:, i])
                                       for i in range(N)])
    # Sample mean + covariance
    mean = ensemble_pred.mean(axis=1)
    A = ensemble_pred - mean[:, None]
    P = A @ A.T / (N - 1)
    # Kalman gain
    K = P @ H.T @ np.linalg.inv(H @ P @ H.T + R)
    # Update each member with perturbed obs
    ensemble_new = np.zeros_like(ensemble_pred)
    for i in range(N):
        y_perturbed = y_obs + np.random.multivariate_normal(np.zeros_like(y_obs), R)
        ensemble_new[:, i] = ensemble_pred[:, i] + K @ (y_perturbed - H @ ensemble_pred[:, i])
    return ensemble_new
```

## Libraries

| Library | Class | Notes |
|---------|-------|-------|
| **filterpy** | All KF variants | Pythonic, well-documented |
| **pykalman** | KF / smoother | Old but stable |
| **GeoStat-Framework/pyKrige** | Has Kalman | Geo applications |
| **PyDA** | EnKF + 4D-Var | Research-grade |
| **DAPPER** | DA benchmarks | Reproducible numerical experiments |

## Cross-references

- `radia_mcp.mor` — model order reduction for filter state
- `radia_mcp.fusion_reactor` — plasma equilibrium reconstruction
"""


FOUR_DVAR = r"""
# 4D-Var (Four-Dimensional Variational Data Assimilation)

`Recent Progress of Data Assimilation Methods in Meteorology.pdf`
(public-safe curated corpus) covers DA in geophysics where 4D-Var
is the gold standard. The principle transfers directly to EM.

## Idea

Minimize a cost function over a time window [0, T]:

    J(x_0) = 1/2 (x_0 - x_b)^T B^-1 (x_0 - x_b)        (background term)
           + 1/2 Σ_t (y_t - h(M_{0->t}(x_0)))^T R^-1 ...   (observation term)

where:
- x_0 = initial state (the unknown to optimize)
- x_b = background / prior estimate
- B = background covariance
- y_t = observations at time t
- M_{0->t} = forward model from time 0 to t
- h = observation operator
- R = observation noise covariance

Minimizing J(x_0) gives the initial condition that best matches
ALL observations over the window, while staying close to the prior.

## Tangent linear + adjoint

Gradient of J requires propagating sensitivity backward through M:

    ∂J/∂x_0 = B^-1 (x_0 - x_b) - M_{0->T}^T h^T R^-1 (y - h(M_{0->T}(x_0)))

The `M^T` is the **adjoint** of the forward model. For EM:
- Linear forward (Maxwell in linear material): adjoint = transpose
  of the FEM matrix
- Nonlinear forward: adjoint = tangent linear (TLM) linearized at the
  current trajectory; needs differentiation through nonlinearity

JAX / PyTorch make TLM + adjoint automatic.

## When 4D-Var beats EnKF

- Want smooth state trajectory (no jumps at update times)
- Have a known background / prior (B matrix)
- Window length ~ 6-24 hours (or in EM, ~ several time steps)
- Differentiable forward model

## When EnKF wins

- Non-differentiable forward (commercial FEM)
- Strongly non-Gaussian noise
- Online / no fixed window
- Easy parallelization (ensemble)

## EM example: workpiece eddy current trajectory

Given coil voltage + power measurements over T=1 second of IH heating,
estimate the time-evolving J(x,y,z,t) in the workpiece.

State x(t) = workpiece eddy current distribution (FEM HCurl DOFs).
Forward M: time-stepping FEM with coil source.
Observations: terminal V, I, P measurements.

```python
import jax.numpy as jnp
from jax import grad, jit

@jit
def cost(x_0, observations, background, R_inv, B_inv):
    # Forward simulate
    states = [x_0]
    for u in coil_drive_sequence:
        states.append(time_step(states[-1], u))
    # Predict observations
    obs_pred = jnp.array([h(s) for s in states])
    # Cost
    J_obs = jnp.sum((obs_pred - observations) * R_inv * (obs_pred - observations))
    J_bg = (x_0 - background) @ B_inv @ (x_0 - background)
    return 0.5 * (J_obs + J_bg)

# Minimize w.r.t. x_0 via L-BFGS
from jax.scipy.optimize import minimize
result = minimize(cost, x_0_initial, method='l-bfgs-b',
                   args=(observations, background, R_inv, B_inv))
```

## Cross-links

- `radia_mcp.pinn` — neural surrogate for TLM/adjoint
- `radia_mcp.mor` — reduce state dim before 4D-Var
"""


def get_data_assimilation_documentation(topic: str = "all") -> str:
    """Dispatch DA topics.

    Topics:
        overview         - DA family + when which method for EM
        kalman_family    - KF / EKF / UKF / EnKF for sensor fusion
        four_dvar        - 4D-Var with tangent linear + adjoint
        all              - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro"):
        return OVERVIEW
    if topic in ("kalman_family", "kalman", "kf", "ekf", "ukf", "enkf"):
        return KALMAN_FAMILY
    if topic in ("four_dvar", "4d_var", "4dvar", "variational"):
        return FOUR_DVAR
    if topic == "all":
        return "\n\n".join([OVERVIEW, KALMAN_FAMILY, FOUR_DVAR])
    return (f"Unknown topic '{topic}'. Available: overview, "
            "kalman_family, four_dvar, all.")
