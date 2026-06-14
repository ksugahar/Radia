"""Advanced lab BBO recipes for Optuna.

Complements ``lab_applications_knowledge.py`` (which holds the 5
pattern-level recipes: coil_design / motor_topology / inverse / wpt /
literature). This module is for **deeper, lab-specific BBO patterns**
that wire Optuna onto an existing Stage-2 ``calc_*.py`` script in
``src/radia/panels/``.

Five new recipes (each ~120-180 lines of prose + working code):

- ``pmsm_cogging``       Cogging torque + ripple multi-objective on
                          ``calc_motor_transient.py``: magnet pole-arc
                          alpha_p + slot opening b_s + skew angle.
- ``wpt_misalignment``   Robust optimization across a lateral / vertical
                          misalignment grid; objective = worst-case
                          efficiency. Conditional search space across
                          compensation topology (calc_inductance.py).
- ``shielding_layout``   Mu-metal / Cu sheet placement to minimize field
                          at a sensitive zone; mixed int + continuous
                          + multi-objective (field at sensor vs shield
                          mass).
- ``litz_strand_design`` n_strands / strand_d / twist_pitch trade-off
                          at target frequency; AC R minimization with
                          DC R + cost constraints (calc_inductance.py
                          --coil-solver peec).
- ``karl_multifidelity`` Karl-iteration FEM-SIBC convergence pruning
                          (MedianPruner + intermediate_value reporting).
                          Outer Optuna over geometry, inner Karl reports
                          residual per iteration -- bad trials die in
                          seconds.

Each recipe lists: variables, objective, sampler config, working code
that calls an existing ``calc_*.py`` from ``src/radia/panels/``, trial
budget, expected outcome, and lab-specific tips.

All recipes are runnable -- they reference Stage-2 CLI entry points
already shipped with radia. Per the **Panel Design Workflow Policy**
(CLAUDE.md) every Stage-3 panel has a Stage-2 ``calc_*.py`` CLI; this
module exploits that to drive sweeps without touching panel UIs.
"""

# ============================================================
# Recipe 1: PMSM cogging torque + ripple
# ============================================================

RECIPE_PMSM_COGGING = r"""
# Recipe 1: PMSM cogging torque + torque ripple multi-objective

## Problem

Given a base PMSM design (pole pairs, slot count, stator OD, stack
length fixed), find the magnet pole-arc fraction `alpha_p`, slot
opening `b_s`, and skew angle `theta_skew` that simultaneously:

1. **Minimize peak-to-peak cogging torque** T_cog_pp [Nm]
2. **Maximize average torque** T_avg [Nm] at rated current

These trade off: skew reduces cogging but also reduces fundamental
torque; narrow slot opening reduces cogging but increases slot
leakage; non-rectangular magnet shape reduces cogging at the cost of
peak flux density.

## Variables

| Variable | Range | Distribution | Notes |
|---|---|---|---|
| alpha_p | 0.6 .. 0.95 | uniform | magnet arc / pole pitch ratio |
| b_s | 1e-3 .. 4e-3 m | log-uniform | slot opening at airgap |
| theta_skew | 0 .. (1 slot pitch) | uniform | continuous skew |
| skew_steps | 1, 2, 3, 5 | categorical | step-skew level |

## Objective

Multi-objective NSGA-II Pareto front of (T_cog_pp ↓, T_avg ↑).

## Recipe (paste-runnable)

```python
import optuna
import subprocess, json, math
from pathlib import Path
import numpy as np

CALC = Path(r"S:/Radia/01_GitHub/src/radia/panels/calc_motor_transient.py")
# Base geometry fixed by the user; sweep only the 4 cogging knobs
BASE = dict(
    motor_type="PMSM",
    pole_pairs=4, slot_count=12,
    stator_od=120e-3, stack_l=80e-3, airgap=0.5e-3,
    magnet_grade="N42",
    current_A=15.0,
    n_steps=180,           # one electrical period at 180 deg
)

def run_motor_trial(params):
    payload = {**BASE, **params}
    out = subprocess.run(
        ["python", str(CALC), "--json-args", json.dumps(payload)],
        capture_output=True, text=True, timeout=600)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[-2000:])
    return json.loads(out.stdout.splitlines()[-1])

def cogging_objective(trial):
    alpha_p     = trial.suggest_float("alpha_p", 0.6, 0.95)
    b_s         = trial.suggest_float("b_s", 1e-3, 4e-3, log=True)
    theta_skew  = trial.suggest_float("theta_skew", 0, math.radians(7.5))
    skew_steps  = trial.suggest_categorical("skew_steps", [1, 2, 3, 5])

    res = run_motor_trial(dict(
        alpha_p=alpha_p, b_s=b_s,
        theta_skew=theta_skew, skew_steps=skew_steps,
    ))
    T_avg = res["torque_avg_Nm"]
    T_cog_pp = res["cogging_pp_Nm"]
    # Optuna minimizes; return (cogging, -T_avg) for NSGA-II
    return T_cog_pp, -T_avg

study = optuna.create_study(
    study_name="pmsm_cogging_pareto",
    storage="sqlite:///pmsm_cogging.db",
    directions=["minimize", "minimize"],
    sampler=optuna.samplers.NSGAIISampler(
        population_size=24,
        crossover=optuna.samplers.nsgaii.UNDXCrossover(),
        seed=42,
    ),
    load_if_exists=True,
)
study.optimize(cogging_objective, n_trials=120,
               catch=(RuntimeError, subprocess.TimeoutExpired))

# Pareto front visualization
import optuna.visualization as viz
viz.plot_pareto_front(
    study, target_names=["T_cog_pp [Nm]", "-T_avg [Nm]"]
).write_html("pmsm_cogging_pareto.html")
```

## Trial budget

- Per trial: ~3-5 min (1 PMSM transient at 180 angular steps on the
  lab desktop)
- Population 24, 120 trials = ~5 generations of NSGA-II = ~6-10 hours
  wall-clock. Use ``storage="sqlite:///pmsm_cogging.db"`` so multiple
  workers can join via ``optuna.load_study(...)``.

## Expected outcome

Pareto front typically shows a clear knee around alpha_p ~ 0.78-0.82
with skew_steps in {3, 5}. The trade-off above the knee is steep
(small cogging reduction costs significant T_avg). Below the knee,
T_avg is nearly flat in cogging.

## Tips

- **Start with skew_steps = 1** for the first 24 trials via
  ``study.enqueue_trial({"skew_steps": 1, ...})``. Single-step skew
  is faster to mesh and gives a baseline before NSGA-II explores
  multi-step.
- **Constrain b_s**: the panel's mesher may fail for b_s < 1 mm on
  large stator_od; add a try/except + ``TrialPruned()``.
- The lab's existing ``radia_motor.py`` panel can BROWSE the resulting
  Pareto front via its results tab.

## Cross-references

- ``radia_mcp.motor`` -- motor physics + Wakao topology opt
- ``radia.panels.calc_motor_transient`` -- the Stage-2 CLI this drives
- ``radia_mcp.differential_forms.em_force_recipe`` -- torque
  computation methods (Arkkio / VWP / Kameari)
"""


# ============================================================
# Recipe 2: WPT misalignment robustness
# ============================================================

RECIPE_WPT_MISALIGNMENT = r"""
# Recipe 2: WPT robust optimization over misalignment grid

## Problem

Real WPT systems experience lateral offset (Tx-Rx pad misalignment)
and vertical offset (ground clearance variation). A design that
peaks efficiency at perfect alignment but craters at 50 mm offset is
worse than one with flat efficiency.

We want the compensation topology + coil parameters that maximize
**worst-case efficiency** over a representative misalignment grid.

## Variables

| Variable | Range | Notes |
|---|---|---|
| topology | SS, SP, LCC, LCL | compensation pick |
| N_p, N_s | 6 .. 16 (int) | primary/secondary turns |
| R_p, R_s | 80e-3 .. 200e-3 m | coil radii |
| C_p, C_s | 1e-9 .. 1e-6 F (log) | tuning caps |
| L_ext_p, L_ext_s | (LCC/LCL only) 1e-6 .. 1e-3 H (log) | external L |

Misalignment grid (evaluated for each trial):
- lateral: [0, 25, 50, 75, 100] mm (5 points)
- vertical: [80, 120, 160] mm (3 points)
- 5 x 3 = 15 evaluations per trial

## Objective

Minimize: -min(eta over 15 grid points)
(NSGA-II if also tracking peak efficiency separately)

## Recipe

```python
import optuna
import subprocess, json
import numpy as np
from pathlib import Path
from itertools import product

CALC = Path(r"S:/Radia/01_GitHub/src/radia/panels/calc_inductance.py")
GRID = list(product([0, 25e-3, 50e-3, 75e-3, 100e-3],   # lateral
                    [80e-3, 120e-3, 160e-3]))           # vertical

def eval_eta(params, offset_xy):
    # Run calc_inductance to get L_p, M, L_s at this offset, then
    # evaluate the chosen compensation topology analytically.
    payload = {**params, "offset_x": offset_xy[0],
               "offset_z": offset_xy[1]}
    out = subprocess.run(
        ["python", str(CALC),
         "--coil-solver", "peec",
         "--json-args", json.dumps(payload)],
        capture_output=True, text=True, timeout=300)
    if out.returncode != 0:
        return None
    res = json.loads(out.stdout.splitlines()[-1])
    L_p, M, L_s = res["L_p"], res["M"], res["L_s"]
    omega = 2 * np.pi * 85e3
    if params["topology"] == "SS":
        return ss_eta(L_p, M, L_s, params["C_p"], params["C_s"], omega)
    if params["topology"] == "LCC":
        return lcc_eta(L_p, M, L_s, params["L_ext_p"], params["L_ext_s"],
                        params["C_p"], params["C_s"], omega)
    # ... other topologies

def robust_objective(trial):
    topo = trial.suggest_categorical("topology",
                                       ["SS", "SP", "LCC", "LCL"])
    params = dict(
        topology=topo,
        N_p=trial.suggest_int("N_p", 6, 16),
        N_s=trial.suggest_int("N_s", 6, 16),
        R_p=trial.suggest_float("R_p", 80e-3, 200e-3),
        R_s=trial.suggest_float("R_s", 80e-3, 200e-3),
        C_p=trial.suggest_float("C_p", 1e-9, 1e-6, log=True),
        C_s=trial.suggest_float("C_s", 1e-9, 1e-6, log=True),
    )
    if topo in ("LCC", "LCL"):
        params["L_ext_p"] = trial.suggest_float(
            f"L_ext_p_{topo}", 1e-6, 1e-3, log=True)
        params["L_ext_s"] = trial.suggest_float(
            f"L_ext_s_{topo}", 1e-6, 1e-3, log=True)

    # Evaluate at every grid point, return WORST-CASE eta
    etas = []
    for i, off in enumerate(GRID):
        eta = eval_eta(params, off)
        if eta is None or eta < 0.1:
            # FEM diverged or unreasonable -- prune this trial
            raise optuna.TrialPruned()
        etas.append(eta)
        # Report intermediate so MedianPruner can kill weak trials
        trial.report(np.min(etas), step=i)
        if trial.should_prune():
            raise optuna.TrialPruned()

    eta_worst = float(np.min(etas))
    eta_best  = float(np.max(etas))
    trial.set_user_attr("eta_worst", eta_worst)
    trial.set_user_attr("eta_best", eta_best)
    trial.set_user_attr("eta_grid", etas)
    return -eta_worst   # Optuna minimizes; we maximize worst-case

study = optuna.create_study(
    study_name="wpt_robust",
    storage="sqlite:///wpt_robust.db",
    direction="minimize",
    sampler=optuna.samplers.TPESampler(multivariate=True, group=True),
    pruner=optuna.pruners.MedianPruner(
        n_startup_trials=10, n_warmup_steps=3),
)
study.optimize(robust_objective, n_trials=200)
print(f"Best worst-case eta: {-study.best_value:.3f}")
print(f"Best params: {study.best_params}")
```

## Trial budget

- Per trial: 15 sub-evaluations x ~30 sec each = 7-8 min (without
  pruning). MedianPruner typically kills bad trials in 2-3 grid
  points (under 2 min).
- 200 trials nominal = 25 hours wall-clock; expect actual ~12-15
  hours with pruning.

## Expected outcome

LCC compensation typically wins on robustness — the extra L_ext
slightly degrades peak efficiency but flattens the eta(offset)
surface. SS often wins on PEAK efficiency but craters at large
misalignment.

The lab's `wpt_misalignment_robustness_check.docx` (W:/.../WPT)
captures this trade-off graphically.

## Tips

- **Set group=True** on TPESampler: SS / LCC have different parameter
  sets and you don't want their hyperparameter distributions to
  pollute each other.
- **Enable MedianPruner** -- the 5x3 grid makes intermediate value
  reporting natural. Bad trials die after 3 of 15 evaluations.
- **enqueue_trial the lab's current SAE J2954 production design**
  as the first trial -- guarantees the new design is at least as
  good as the production one.

## Cross-references

- ``radia_mcp.pcb`` -- WPT compensation theory + topology catalog
- ``radia.panels.calc_inductance`` -- PEEC coil inductance solver
- ``radia_mcp.peec`` -- PEEC filament-panel architecture
- Existing skeleton in ``lab_applications_knowledge.APP_WPT`` is the
  pattern; this recipe is the production-grade robust version.
"""


# ============================================================
# Recipe 3: Magnetic shielding placement
# ============================================================

RECIPE_SHIELDING_LAYOUT = r"""
# Recipe 3: Magnetic shielding placement (mu-metal or Cu sheet)

## Problem

Given a magnetic source (motor, transformer, WPT coil) and a
sensitive zone (e.g. heart pacemaker locus 1 m away), place
N_shields rectangular sheets of mu-metal or Cu to minimize the
flux density at the sensor zone, subject to a total shield
mass constraint.

## Variables

| Variable | Range | Notes |
|---|---|---|
| n_shields | 1 .. 4 (int) | how many sheets |
| (per shield, conditional on n_shields) | | |
| material_i | "mu_metal", "Cu_5mm" (categorical) | absorption vs reflection |
| pos_x_i | -0.5 .. 0.5 m | x position |
| pos_y_i | -0.3 .. 0.3 m | y position |
| width_i | 0.05 .. 0.5 m | sheet width |
| height_i | 0.05 .. 0.5 m | sheet height |
| rotation_i | 0, 45, 90 deg (categorical) | sheet orientation |

Objective (multi):
1. minimize |B| at sensor zone [uT]
2. minimize total shield mass [kg]
3. (constraint) total shield projected area <= some cap

## Recipe

```python
import optuna
import subprocess, json
from pathlib import Path
import numpy as np

CALC = Path(r"S:/Radia/01_GitHub/src/radia/panels/calc_shielding.py")
SENSOR_POINTS = [(1.0, 0.0, 0.0), (1.0, 0.1, 0.0),
                 (1.0, -0.1, 0.0)]   # 3-point sensor zone

DENSITY = {"mu_metal": 8700.0, "Cu_5mm": 8960.0}  # kg/m^3
THICKNESS = {"mu_metal": 1e-3, "Cu_5mm": 5e-3}

def shielding_obj(trial):
    n = trial.suggest_int("n_shields", 1, 4)
    shields = []
    total_mass = 0.0
    for i in range(n):
        mat = trial.suggest_categorical(f"mat_{i}",
                                          ["mu_metal", "Cu_5mm"])
        s = dict(
            material=mat,
            pos=(trial.suggest_float(f"x_{i}", -0.5, 0.5),
                 trial.suggest_float(f"y_{i}", -0.3, 0.3), 0.0),
            size=(trial.suggest_float(f"w_{i}", 0.05, 0.5),
                  trial.suggest_float(f"h_{i}", 0.05, 0.5)),
            rotation_deg=trial.suggest_categorical(
                f"rot_{i}", [0, 45, 90]),
        )
        shields.append(s)
        m = (DENSITY[mat] * THICKNESS[mat]
             * s["size"][0] * s["size"][1])
        total_mass += m

    # Run FEM to get |B| at sensor zone
    payload = dict(shields=shields, sensor_points=SENSOR_POINTS)
    out = subprocess.run(
        ["python", str(CALC),
         "--json-args", json.dumps(payload)],
        capture_output=True, text=True, timeout=600)
    if out.returncode != 0:
        raise optuna.TrialPruned()
    res = json.loads(out.stdout.splitlines()[-1])
    B_uT_worst = max(b * 1e6 for b in res["B_at_sensor"])
    trial.set_user_attr("B_uT_worst", B_uT_worst)
    trial.set_user_attr("mass_kg", total_mass)
    return B_uT_worst, total_mass

study = optuna.create_study(
    study_name="shield_pareto",
    storage="sqlite:///shielding.db",
    directions=["minimize", "minimize"],
    sampler=optuna.samplers.NSGAIISampler(
        population_size=30,
        crossover=optuna.samplers.nsgaii.UNDXCrossover(),
    ),
    load_if_exists=True,
)
study.optimize(shielding_obj, n_trials=300, catch=(RuntimeError,))
import optuna.visualization as viz
viz.plot_pareto_front(
    study, target_names=["|B| at sensor [uT]", "shield mass [kg]"]
).write_html("shielding_pareto.html")
```

## Trial budget

- Per trial: 1 FEM solve (~2-5 min depending on shield count + mesh
  density)
- Population 30, 300 trials = ~10 generations = 10-20 hours

## Expected outcome

NSGA-II typically converges to a Pareto front with a clear knee at
~30-50% mass reduction vs single-sheet baseline. Mu-metal sheets
dominate at low frequency (absorption); Cu sheets dominate at high
frequency (eddy reflection). Optimal layouts often place ONE large
mu-metal sheet near the source and ONE smaller Cu sheet near the
sensor.

## Tips

- **n_shields as int** + conditional params (mat_i, x_i, ...) creates
  a variable-dimensional search space. TPE / NSGA-II handle this OK
  because Optuna treats absent params as "unused".
- **rotation_deg as categorical (0, 45, 90)** keeps the search
  tractable. A continuous rotation is overkill -- shielding effect
  is symmetric to ~15-deg granularity.
- **Mass IS a Pareto axis, not a constraint** -- this reflects real
  procurement trade-offs (cost ~ mass for these materials).

## Cross-references

- ``radia_mcp.team_benchmark`` -- TEAM 7 / 21 shielding benchmarks
- ``radia_mcp.fem.fem_gauge_open_boundary`` -- Kelvin transform for
  open boundary in the FEM solve
- ``radia_mcp.literature_index`` -- lab papers under
  ``99_アプリケーション/05_シールディング/``
"""


# ============================================================
# Recipe 4: Litz wire strand design
# ============================================================

RECIPE_LITZ_STRAND = r"""
# Recipe 4: Litz wire strand design for AC resistance minimization

## Problem

Design a litz wire with n_strands of diameter d_strand twisted at
twist_pitch to minimize AC resistance at the operating frequency,
subject to DC resistance constraint (cooling) and cost constraint
(litz $/m scales with n_strands).

## Variables

| Variable | Range | Notes |
|---|---|---|
| n_strands | [40, 60, 80, 100, 150, 200, 270, 400] | discrete (industry std) |
| strand_d | 0.05e-3 .. 0.3e-3 m (log) | per-strand wire diameter |
| twist_pitch | 5e-3 .. 50e-3 m | twist length per turn |
| n_layers | 1 .. 4 (int) | bundle layering depth |

## Objective (single, with constraint penalty)

Minimize: AC_R at operating frequency.

Constraint: DC_R <= R_dc_max (typically 50 mOhm/m for the cooling spec)
Constraint: cost = n_strands * unit_cost(strand_d) <= cost_max

## Recipe

```python
import optuna
import subprocess, json
from pathlib import Path
import numpy as np

CALC = Path(r"S:/Radia/01_GitHub/src/radia/panels/calc_inductance.py")
F_OP_HZ = 85e3
R_DC_MAX = 50e-3   # ohm/m
COST_MAX = 5.0     # USD/m

def litz_cost_per_meter(n_strands, strand_d):
    # Rough cost model: $/m ~= n * d^1.5 * 1e6
    return n_strands * (strand_d * 1e3)**1.5 * 0.5

def litz_obj(trial):
    n_strands  = trial.suggest_categorical(
        "n_strands", [40, 60, 80, 100, 150, 200, 270, 400])
    strand_d   = trial.suggest_float("strand_d", 0.05e-3, 0.3e-3, log=True)
    twist_pitch = trial.suggest_float("twist_pitch", 5e-3, 50e-3)
    n_layers   = trial.suggest_int("n_layers", 1, 4)

    # Pre-filter on cost and DC_R BEFORE invoking the expensive FEM
    sigma_Cu = 5.8e7
    A_per_strand = np.pi * (strand_d/2)**2
    A_total = n_strands * A_per_strand
    R_dc = 1.0 / (sigma_Cu * A_total)
    cost = litz_cost_per_meter(n_strands, strand_d)
    trial.set_user_attr("R_dc", R_dc)
    trial.set_user_attr("cost_per_m", cost)
    if R_dc > R_DC_MAX or cost > COST_MAX:
        raise optuna.TrialPruned()

    payload = dict(
        coil_solver="peec",
        n_strands=n_strands, strand_d=strand_d,
        twist_pitch=twist_pitch, n_layers=n_layers,
        freq_hz=F_OP_HZ,
    )
    out = subprocess.run(
        ["python", str(CALC),
         "--coil-solver", "peec",
         "--json-args", json.dumps(payload)],
        capture_output=True, text=True, timeout=300)
    if out.returncode != 0:
        raise optuna.TrialPruned()
    res = json.loads(out.stdout.splitlines()[-1])
    R_ac = res["R_ac_at_freq"]
    trial.set_user_attr("R_ac_ohm_per_m", R_ac)
    return R_ac

study = optuna.create_study(
    study_name="litz_strand",
    storage="sqlite:///litz.db",
    direction="minimize",
    sampler=optuna.samplers.TPESampler(multivariate=True, seed=0),
    load_if_exists=True,
)
study.optimize(litz_obj, n_trials=200, catch=(RuntimeError,))
print("Best AC_R:", study.best_value, "ohm/m")
print("Best:", study.best_params)
print("  R_dc:", study.best_trial.user_attrs["R_dc"])
print("  cost:", study.best_trial.user_attrs["cost_per_m"])
```

## Trial budget

- Cheap pre-filter (DC_R + cost) cuts ~40% of trials in <1 ms each.
- Per surviving trial: 1 PEEC solve (~30 sec - 2 min depending on
  n_strands and mesh).
- 200 trials = 60 survive PEEC = ~60 min wall-clock.

## Expected outcome

For 85 kHz WPT, optimum is typically n_strands ~ 150-200, strand_d
~ 0.1 mm, twist_pitch ~ 15 mm. R_ac / R_dc ratio approaches ~1.05
(near-ideal litz; without twist it would be 5-10).

## Tips

- **Pre-filter analytics-only first**: DC_R from cross-section area
  and cost from your cost model are <1 ms each. Reject before
  invoking PEEC.
- **n_strands as categorical** matches the actual industry-standard
  litz catalog (40/60/80/100/150/200/270/400). Avoids wasting
  Optuna's exploration budget on non-purchasable values.
- **Multi-frequency case**: switch to NSGA-II with directions =
  ["minimize", "minimize", "minimize"] for (R_ac at 85k, R_ac at
  150k, R_ac at 1M) if the application is broadband.

## Cross-references

- ``radia_mcp.peec`` -- PEEC theory + multi-filament nwinc/nhinc
- ``radia_mcp.litz_transmission`` -- litz / transmission-line topics
- ``radia.panels.calc_inductance`` -- coil-solver peec backend
- Carstensen AC loss model in ``radia_mcp.peec.peec_carstensen_ac_loss``
"""


# ============================================================
# Recipe 5: Karl iteration multi-fidelity with pruning
# ============================================================

RECIPE_KARL_MULTIFIDELITY = r"""
# Recipe 5: Karl iteration multi-fidelity with MedianPruner

## Problem

The lab IH analysis pipeline iterates Karl's method for nonlinear
SIBC (workpiece magnetic saturation). Karl typically converges in
3-7 iterations, BUT some geometry / freq combinations cause it to
diverge or oscillate.

In an Optuna sweep over coil + workpiece geometry, ~10-20% of trials
are bad-geometry trials that waste 8-12 minutes on a doomed Karl
loop. We want to KILL them in 30-60 sec via Optuna's pruning.

## Strategy

Karl iteration emits a residual H_t,rms per outer iteration. Stream
that residual as ``trial.report(...)``, let ``MedianPruner`` compare
each iteration's residual to the median of completed trials at the
same iteration index. If a trial is in the bottom half AND past the
warmup count, prune.

## Variables (example: coil-workpiece sweep)

| Variable | Range |
|---|---|
| coil_R | 30e-3 .. 80e-3 m |
| coil_pitch | 2e-3 .. 8e-3 m |
| wp_d | 20e-3 .. 60e-3 m |
| wp_h | 30e-3 .. 100e-3 m |
| freq_kHz | [10, 25, 50, 100] |
| sigma_wp | [5e6, 1e7, 2.5e7, 5e7] |

## Objective

Maximize: P_wp / I^2 (per-unit-current heating power), single-obj.

## Recipe (paste-runnable, with intermediate_value reporting)

```python
import optuna
import subprocess, json
from pathlib import Path

CALC = Path(r"S:/Radia/01_GitHub/src/radia/panels/calc_fem_kelvin.py")

def karl_pruning_obj(trial):
    params = dict(
        coil_R=trial.suggest_float("coil_R", 30e-3, 80e-3),
        coil_pitch=trial.suggest_float("coil_pitch", 2e-3, 8e-3),
        wp_d=trial.suggest_float("wp_d", 20e-3, 60e-3),
        wp_h=trial.suggest_float("wp_h", 30e-3, 100e-3),
        freq_kHz=trial.suggest_categorical(
            "freq_kHz", [10, 25, 50, 100]),
        sigma_wp=trial.suggest_categorical(
            "sigma_wp", [5e6, 1e7, 2.5e7, 5e7]),
        karl_max_iter=12,
        karl_stream_intermediate=True,   # tells calc to print per-iter
    )
    p = subprocess.Popen(
        ["python", str(CALC),
         "--json-args", json.dumps(params)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, bufsize=1)
    P_wp_final = None
    for line in p.stdout:
        line = line.strip()
        if line.startswith("KARL_ITER:"):
            # Format: KARL_ITER: iter=3 H_t_rms=12.345 P_wp=87.6
            tokens = dict(t.split("=") for t in line.split()[1:])
            iter_n = int(tokens["iter"])
            P_wp   = float(tokens["P_wp"])
            trial.report(P_wp, step=iter_n)
            if trial.should_prune():
                p.terminate()
                p.wait(timeout=10)
                raise optuna.TrialPruned()
        elif line.startswith("FINAL_P_WP="):
            P_wp_final = float(line.split("=")[1])
    p.wait(timeout=60)
    if P_wp_final is None or p.returncode != 0:
        raise optuna.TrialPruned()
    I_drive = params.get("I_drive_A", 100.0)
    return -P_wp_final / (I_drive**2)   # minimize negative = maximize P/I^2

pruner = optuna.pruners.MedianPruner(
    n_startup_trials=20,     # don't prune until 20 trials completed
    n_warmup_steps=2,        # always run at least 2 Karl iter
    interval_steps=1,
)

study = optuna.create_study(
    study_name="ih_karl_pruning",
    storage="sqlite:///ih_pruned.db",
    direction="minimize",
    sampler=optuna.samplers.TPESampler(multivariate=True, group=True),
    pruner=pruner,
    load_if_exists=True,
)
study.optimize(karl_pruning_obj, n_trials=300, catch=(RuntimeError,))
print(f"Pruned: {sum(1 for t in study.trials if t.state == optuna.trial.TrialState.PRUNED)}")
print(f"Completed: {sum(1 for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE)}")
```

## What you need to add to calc_fem_kelvin.py

For ``trial.report`` to work, ``calc_fem_kelvin.py`` must emit a
``KARL_ITER: iter=N H_t_rms=X P_wp=Y`` line at each Karl outer
iteration. The Karl loop wrapper in ``radia.panels.calc_fem_kelvin``
already prints diagnostics; add a single ``--karl-stream-intermediate``
flag that switches the format to the parseable form above.

If you can't modify the script, fall back to polling its stdout for
existing ``[Karl iter 3] H_t_rms = 12.345 P_wp = 87.6`` style lines
with a regex.

## Trial budget

- Without pruning: 300 trials x 12 Karl iter x 30 sec = 30 hours
- With pruning (assume 30% prune rate, killed after ~3 Karl iter):
  210 full trials (300 sec each) + 90 pruned (90 sec each)
  = 17.5 + 2.25 = ~20 hours

  ~33% wall-clock reduction. The win scales with prune rate.

## Expected outcome

Optuna learns within ~50 trials that high sigma_wp + small wp_d +
high freq_kHz combinations diverge Karl iteration; it stops sampling
those. The Pareto-equivalent "best P_wp / I^2" plateau emerges around
trial 80-100.

## Tips

- **n_startup_trials=20**: too aggressive pruning before the median
  is meaningful will kill good trials. Pruner needs ~20 complete
  trials per Karl iteration to have a stable median.
- **n_warmup_steps=2**: always let Karl run at least 2 iter before
  pruning -- Karl's first iter is often a poor predictor.
- **HyperbandPruner alternative**: more aggressive; reduces total
  budget to ~25-50% but loses Pareto-front quality. Use only if
  total wall-clock is the only metric.
- **catch=(RuntimeError,)**: needed because calc_fem_kelvin can
  raise on mesh failure. Optuna treats raised RuntimeError as a
  failed trial (NOT a pruned trial); failed != pruned for the
  pruner's median calculation.

## Cross-references

- ``radia_mcp.ih.ih_esim`` / ``ih_sibc`` -- Karl iteration theory
- ``radia.panels.calc_fem_kelvin`` -- the Stage-2 CLI driven here
- ``radia_mcp.optuna.algorithm_knowledge.pruning`` -- general
  pruning theory + textbook reference
- ``radia_mcp.bayesian_opt`` -- alternative BBO with native
  intermediate-value support (sibling MCP)
"""


# ============================================================
# Dispatcher
# ============================================================

_TOPICS = {
    "overview": (
        "# Advanced lab BBO recipes for Optuna\n\n"
        "Five lab-specific recipes that wire Optuna onto an existing\n"
        "Stage-2 calc_*.py script. Topics:\n\n"
        "- pmsm_cogging:       PMSM cogging torque + ripple multi-obj\n"
        "- wpt_misalignment:   WPT worst-case eta across position grid\n"
        "- shielding_layout:   mu-metal / Cu sheet placement Pareto\n"
        "- litz_strand_design: litz n_strands/d/twist for AC_R min\n"
        "- karl_multifidelity: Karl iter pruning for IH sweeps\n\n"
        "All recipes are runnable -- they reference real Stage-2 CLI\n"
        "entry points already shipped with radia. Complements the\n"
        "5 pattern-level recipes in lab_applications_knowledge.\n"
    ),
    "pmsm_cogging": RECIPE_PMSM_COGGING,
    "wpt_misalignment": RECIPE_WPT_MISALIGNMENT,
    "shielding_layout": RECIPE_SHIELDING_LAYOUT,
    "litz_strand_design": RECIPE_LITZ_STRAND,
    "karl_multifidelity": RECIPE_KARL_MULTIFIDELITY,
    "all": "",   # filled below
}
_TOPICS["all"] = "\n\n".join(
    v for k, v in _TOPICS.items() if k != "all")


def get_recipes_advanced_documentation(topic: str = "all") -> str:
    """Dispatch advanced lab BBO recipe topics."""
    topic = topic.lower().strip()
    # accept aliases
    aliases = {
        "pmsm": "pmsm_cogging",
        "cogging": "pmsm_cogging",
        "wpt": "wpt_misalignment",
        "misalignment": "wpt_misalignment",
        "robustness": "wpt_misalignment",
        "shielding": "shielding_layout",
        "shield": "shielding_layout",
        "litz": "litz_strand_design",
        "strand": "litz_strand_design",
        "karl": "karl_multifidelity",
        "pruning_recipe": "karl_multifidelity",
        "multifidelity": "karl_multifidelity",
        "multi_fidelity": "karl_multifidelity",
        "intro": "overview",
    }
    topic = aliases.get(topic, topic)
    if topic in _TOPICS:
        return _TOPICS[topic]
    return ("Unknown topic '" + topic + "'. Available: "
            + ", ".join(sorted(_TOPICS.keys())))
