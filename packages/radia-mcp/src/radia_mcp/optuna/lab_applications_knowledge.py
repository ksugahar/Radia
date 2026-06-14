"""Lab-specific applications of Optuna for EM / Radia / NGSolve work.

This module documents how Optuna plugs into the existing radia-mcp
ecosystem -- which design problems are good BBO fits, and which are
better served by gradient-based / analytical methods.
"""

LAB_OVERVIEW = r"""
# Optuna in Sugahara-lab workflows

The lab's existing optimization machinery is gradient-based:
- `radia_mcp.topology_optimization` (SIMP / level-set / MMA)
- `radia_mcp.motor` Wakao autoencoder + LS SynRM (decoder gradient)
- `radia_mcp.accelerator` coil-end shape sensitivity (Gangl 2015)

Optuna fills the gaps where gradients are unavailable or expensive:

| Lab problem | Why BBO | Optuna feature |
|-------------|---------|----------------|
| Coil winding TOPOLOGY choice + per-topology params | Categorical + conditional | TPE + conditional search space |
| WPT compensation circuit (SS / SP / LCC / LCL) | 4 discrete topologies + L/C per | `suggest_categorical` + branch |
| Material grade (B-H curve selection) | Discrete lookup table | `suggest_categorical` |
| PEEC nwinc / nhinc tuning per litz strand | Discrete grid + skin/proximity | `suggest_int` + heavy compute |
| Coil + workpiece geometry simultaneously | Mixed continuous + discrete | TPE + storage |
| Bayesian inverse problem (BH curve fit) | Smooth but multi-modal | GP / TPE |
| Multi-objective: torque ↑ + loss ↓ | Pareto front needed | NSGA-II |
| Topology + sizing co-optimization | Outer BBO + inner gradient | Nested loop |

## When to use Optuna over the lab's gradient-based tools

Use **gradient-based** (`radia_mcp.topology_optimization`) when:
- The design variable is a continuous field over a fixed mesh
  (density, level set, shape velocity)
- The gradient is computable cheaply (adjoint method)
- A single objective with sensible scalar weighting exists

Use **Optuna** when:
- The design includes **categorical / discrete** choices
- The objective involves a **conditional** workflow (different
  physics models for different topologies)
- You need a **Pareto front** with multi-objective trade-offs
- Each evaluation is **expensive** (~minutes) but not so much that
  GP can be used (`< 100 trials`)
- You want to **explore** before committing to a gradient method

The two compose naturally:
- Optuna outer loop picks topology + initial sizing
- SIMP inner loop refines density field for that topology

## Cross-MCP map

| Sub-package | Tool name | When |
|-------------|-----------|------|
| `radia_mcp.optuna` | this server | All BBO design exploration |
| `radia_mcp.topology_optimization` | gradient-based density / shape | Fixed topology + continuous variables |
| `radia_mcp.motor` | Wakao decoder + LS | Motor cross-section topology |
| `radia_mcp.accelerator` | CoilBuilder + sensitivity | Coil end-region shape |
| `radia_mcp.pcb` | WPT design knowledge | Compensation topology pre-screening |
| `radia_mcp.pinn` | surrogate model | Replace FEM in inner Optuna loop |
| `radia_mcp.literature_index` | corpus search | Find similar BBO benchmarks |
| `radia_mcp.streamfunction` | `streamfunction(regularized)` | SF coil **(homogeneity, peak current density)** Pareto: NSGA-II over (geometry, Tikhonov alpha) + CMA-ES sheet-metal (板金) surface-forming shape optimisation (inner solve = cached folded-Tikhonov RegularizedTSVD) |
"""


APP_COIL_DESIGN = r"""
# Application 1: IH coil design (Optuna + Radia + PEEC + FEM)

## Problem

Design an IH coil + workpiece pair to maximize heating rate of a
given workpiece while minimizing coil power loss.

Variables:
- `n_turns` (int, 3..15) — coil turns
- `R_coil` (float, 20e-3 .. 80e-3 m) — coil mean radius
- `pitch` (float, 2e-3 .. 8e-3 m) — turn-to-turn pitch
- `wire_d` (float, 1e-3 .. 5e-3 m, log) — wire diameter
- `gap` (float, 1e-3 .. 10e-3 m, log) — coil-to-workpiece air gap
- `coil_material` (categorical, ["Cu", "litz_strand_100", "Al"])
- `freq_kHz` (categorical, [10, 25, 50, 100, 200])

Objectives (multi):
1. minimize coil power loss `P_coil` [W]
2. maximize heating rate `dT/dt` ∝ `P_wp` [W/kg]

Constraints:
- coil current density `J_coil <= J_max` (litz vs solid)
- coil voltage `V_coil <= V_drive` (capacitor rating)

## Recipe

```python
import optuna
from radia.panels.calc_inductance import calc_inductance_main
from radia.panels.calc_fem_kelvin import calc_fem_kelvin_main

def ih_objective(trial):
    params = {
        "n_turns": trial.suggest_int("n_turns", 3, 15),
        "R_coil": trial.suggest_float("R_coil", 20e-3, 80e-3),
        "pitch":  trial.suggest_float("pitch", 2e-3, 8e-3),
        "wire_d": trial.suggest_float("wire_d", 1e-3, 5e-3, log=True),
        "gap":    trial.suggest_float("gap", 1e-3, 10e-3, log=True),
        "coil_material": trial.suggest_categorical(
            "coil_material", ["Cu", "litz_strand_100", "Al"]),
        "freq_kHz": trial.suggest_categorical(
            "freq_kHz", [10, 25, 50, 100, 200]),
    }

    # 1) Build coil STEP from params; PEEC for L + R_coil
    peec_result = calc_inductance_main(coil_solver="peec", **params)
    R_coil_ohm = peec_result["R_dc"]
    L_air = peec_result["L_air"]

    # 2) FEM-SIBC for P_wp (workpiece dissipation)
    fem_result = calc_fem_kelvin_main(coil_params=params, **wp_params)
    P_wp = fem_result["P_wp"]

    # Inject derived constraint (current density)
    J_coil = compute_current_density(params, peec_result)
    trial.set_user_attr("constraints", [J_coil - J_max])

    # Compute P_coil from R and target I
    I_drive = solve_for_target_voltage(L_air, R_coil_ohm,
                                        params["freq_kHz"],
                                        V_target=V_drive)
    P_coil = R_coil_ohm * I_drive**2

    return P_coil, -P_wp   # minimize P_coil, maximize P_wp


sampler = optuna.samplers.NSGAIISampler(
    population_size=30,
    crossover=optuna.samplers.nsgaii.UNDXCrossover(),
    constraints_func=lambda t: t.user_attrs["constraints"],
)

study = optuna.create_study(
    study_name="ih_coil_pareto",
    storage="sqlite:///ih_coil.db",
    directions=["minimize", "minimize"],
    sampler=sampler,
    load_if_exists=True,
)

study.optimize(ih_objective, n_trials=500,
               catch=(RuntimeError, ValueError))

# Pareto front
import optuna.visualization as vis
vis.plot_pareto_front(study).write_html("ih_pareto.html")
```

## Tips for this lab application

- **Use Karl iteration pruning** (see `pruning` topic): if Karl
  diverges or Z_s blows up after 3 iter, `TrialPruned()`. Saves
  5-10 minutes per bad trial.
- **Heartbeat + retry** essential — Cubit licenses occasionally hang.
- **`load_if_exists=True`** so you can resume after sweeping all
  night.
- **Warm start with the production coil** via `enqueue_trial(...)`
  — guarantees the new design is at least as good as the old one.
"""


APP_MOTOR_TOPOLOGY = r"""
# Application 2: Motor topology + sizing (Optuna outer + gradient inner)

## Problem

Select motor topology (PMSM / IM / SRM / SynRM) and key dimensions;
let SIMP / shape opt refine the rotor cross-section.

This is a **2-level** problem:
- Outer (Optuna BBO): topology + global dimensions (slot count, pole
  pairs, stator OD, stack length, airgap)
- Inner (gradient): rotor barrier shape (SIMP density field) for the
  selected topology

## Recipe

```python
import optuna
from radia.panels.calc_motor_transient import run_topology_inner
from radia.topology_optimization import simp_inner_solve

def motor_outer(trial):
    motor_type = trial.suggest_categorical(
        "motor", ["PMSM", "SynRM", "IM", "SRM"])

    if motor_type == "PMSM":
        pp     = trial.suggest_int("pp", 3, 8)
        slot_p = trial.suggest_categorical("slot_p", ["9/8", "12/10", "12/14"])
        magnet_grade = trial.suggest_categorical("mag", ["N42", "N48", "N52"])
    elif motor_type == "SynRM":
        pp     = trial.suggest_int("pp", 2, 4)
        n_barriers = trial.suggest_int("n_barriers", 2, 5)
        magnet_grade = None
    elif motor_type == "IM":
        pp     = trial.suggest_int("pp", 1, 4)
        slip   = trial.suggest_float("slip", 0.005, 0.05)
        magnet_grade = None
    else:   # SRM
        pp     = trial.suggest_categorical("pp_srm", [(6,4),(8,6),(12,8)])
        magnet_grade = None

    stator_od = trial.suggest_float("stator_od", 80e-3, 200e-3)
    stack_l   = trial.suggest_float("stack_l", 50e-3, 200e-3)
    airgap    = trial.suggest_float("airgap", 0.3e-3, 1.0e-3, log=True)

    common = dict(stator_od=stator_od, stack_l=stack_l, airgap=airgap)

    # Inner gradient-based refinement (SIMP on rotor for SynRM,
    # magnet shape for PMSM, ...). Returns BEST after inner opt.
    if motor_type == "SynRM":
        inner = simp_inner_solve(motor=motor_type, pp=pp,
                                  n_barriers=n_barriers,
                                  **common, n_iter=50)
    else:
        inner = simp_inner_solve(motor=motor_type, pp=pp,
                                  magnet=magnet_grade,
                                  **common, n_iter=50)

    return -inner["torque_Nm"]   # maximize torque

sampler = optuna.samplers.TPESampler(
    multivariate=True, group=True, seed=42)

study = optuna.create_study(
    study_name="motor_topology_2level",
    storage="sqlite:///motor.db",
    direction="minimize",
    sampler=sampler,
)
study.optimize(motor_outer, n_trials=200, catch=(RuntimeError,))
```

## Cross-references

- `radia_mcp.motor` -- motor physics (Hollaus lamination, Carstensen
  losses, calc_motor_transient.py panel)
- `radia_mcp.topology_optimization` -- gradient inner loop (SIMP/MMA)
- Wakao 2025 autoencoder topology opt: outer Optuna picks the
  autoencoder latent, inner gradient picks density on the latent's
  decoded grid.

## Why this works

- Optuna excels at the discrete topology decision (4 motor types).
- SIMP excels at the continuous density refinement WITHIN a topology.
- Trying to use SIMP alone forces all topologies to share one mesh
  (impossible). Trying to use Optuna alone for the density would
  require ~thousands of trials per topology.
"""


APP_INVERSE = r"""
# Application 3: Inverse identification (Optuna + radia.MatPlayHysteresis)

## Problem

Given measured B(H) loop data, fit Play model parameters (K, eta_k,
shape function f_k(r)) so the model reproduces the loops.

This is multi-modal — many local minima — so gradient descent often
gets stuck. Optuna + TPE / CMA-ES handles it.

## Recipe

```python
import optuna
import numpy as np
import radia as rad

# Measured loop data
H_meas, B_meas = load_measurement("steel_loop_50Hz.csv")

def play_objective(trial):
    K = trial.suggest_int("K", 5, 30)

    # eta_k: monotonically increasing thresholds in T
    # Use log-uniform on increments to keep monotonic
    eta = [trial.suggest_float(f"eta_{k}", 1e-4, 2.0, log=True)
           for k in range(K)]
    eta = np.sort(eta).tolist()

    # Shape function f_k(r) -- piecewise linear with N_r breakpoints
    N_r = trial.suggest_int("N_r", 5, 30)
    f_k_tables = []
    for k in range(K):
        r = np.linspace(0, eta[k], N_r)
        slopes = [trial.suggest_float(f"slope_{k}_{j}", 0, 1e6)
                  for j in range(N_r-1)]
        f = np.concatenate([[0.0], np.cumsum(slopes) * np.diff(r)])
        f_k_tables.append((r.tolist(), f.tolist()))

    mat = rad.MatPlayHysteresis(K, eta, f_k_tables)

    # Apply H_meas sequence, get B_model, compare to B_meas
    B_model = np.array([rad.MatBvsH_PlayCommit(mat, [H, 0, 0])[0]
                        for H in H_meas])
    loss = np.mean((B_model - B_meas)**2)
    return loss

study = optuna.create_study(
    direction="minimize",
    sampler=optuna.samplers.TPESampler(multivariate=True, seed=42),
)
study.optimize(play_objective, n_trials=2000)

print("Best K:", study.best_params["K"])
print("Best fit RMSE:", np.sqrt(study.best_value))
```

## Cross-references

- `radia_mcp.magnetic_materials` -- Play / Energy model theory + lab core
- `radia_mcp.mathematica` -- symbolic check of Play model derivative
- Hysteresis fitting via gradient: only works if K is fixed.
  Optuna's strength is letting K itself be a hyperparameter.
"""


APP_WPT = r"""
# Application 4: WPT compensation topology selection

## Problem

Choose between SS / SP / LCC / LCL compensation topologies and tune
their L/C values for a given coil pair + load.

Variables:
- `topology` (categorical, ["SS", "SP", "LCC", "LCL"])
- For SS: `C_p`, `C_s` (2 caps)
- For SP: `C_p`, `C_s` (different position)
- For LCC: `L_ext_p`, `C_ext_p`, `C_p`, `C_s`, `L_ext_s`, `C_ext_s`
- For LCL: similar 6 params

Objective: maximize `eta = P_load / P_source` at the operating point.

## Recipe

```python
import optuna
import numpy as np

L_p, M, L_s = load_coil_params()   # from PEEC sweep
omega_op = 2*np.pi*85e3            # SAE J2954 operating

def wpt_obj(trial):
    topo = trial.suggest_categorical("topo", ["SS","SP","LCC","LCL"])

    if topo == "SS":
        C_p = trial.suggest_float("C_p", 1e-9, 1e-6, log=True)
        C_s = trial.suggest_float("C_s", 1e-9, 1e-6, log=True)
        eta = ss_efficiency(L_p, M, L_s, C_p, C_s, omega_op)
    elif topo == "SP":
        C_p = trial.suggest_float("Cp_sp", 1e-9, 1e-6, log=True)
        C_s = trial.suggest_float("Cs_sp", 1e-9, 1e-6, log=True)
        eta = sp_efficiency(L_p, M, L_s, C_p, C_s, omega_op)
    elif topo == "LCC":
        L_e = trial.suggest_float("L_ext", 1e-6, 1e-3, log=True)
        C_e = trial.suggest_float("C_ext", 1e-9, 1e-6, log=True)
        C_p = trial.suggest_float("Cp_lcc", 1e-9, 1e-6, log=True)
        C_s = trial.suggest_float("Cs_lcc", 1e-9, 1e-6, log=True)
        L_es = trial.suggest_float("L_ext_s", 1e-6, 1e-3, log=True)
        C_es = trial.suggest_float("C_ext_s", 1e-9, 1e-6, log=True)
        eta = lcc_efficiency(L_p, M, L_s, L_e, C_e, C_p,
                              C_s, L_es, C_es, omega_op)
    else:   # LCL
        ...

    return -eta   # maximize eta

study = optuna.create_study(direction="minimize",
                            sampler=optuna.samplers.TPESampler(
                                multivariate=True, group=True, seed=42))
study.optimize(wpt_obj, n_trials=500)
```

The conditional search space + categorical topology choice is exactly
the conditional pattern from textbook list 2.16.

## Cross-references

- `radia_mcp.pcb` -- WPT physics + compensation topology theory
- After Optuna picks the topology, run `radia_mcp.peec` to refine
  the actual coil for the chosen compensation.
"""


APP_LITERATURE = r"""
# Application 5: Find similar lab BBO benchmarks via literature_index

The lab corpus has ~3,889 PDF / JSON files. Use the literature_index
MCP to find similar BBO use cases before designing a new sweep:

```
mcp-server-literature-index.literature_search(
    "optimization motor topology")
mcp-server-literature-index.literature_search(
    "design optimization Pareto inductor")
mcp-server-literature-index.literature_search(
    "Bayesian optimization electromagnetic")
```

Likely useful lab folders (call `literature_folder_tree` to see):
- `04_機械学習と最適化/` -- this textbook + related
- `99_アプリケーション/01_誘導加熱/` -- IH coil design papers
- `99_アプリケーション/03_モータ/` -- motor topology opt
- `13_形状最適化/` -- shape opt (gradient-based comparison)

After finding 3-5 similar papers, use:

```
mcp-server-document.paper_writing_fetch_and_cite(
    "Bayesian optimization motor design 2023")
```

to pull recent literature for the introduction of any new
Optuna-based publication.

## Cross-references

- `radia_mcp.literature_index` -- corpus search
- `mcp-server-document.paper_writing_*` -- citation tools
"""


def get_lab_applications_documentation(topic: str = "all") -> str:
    """Dispatch lab application topics.

    Topics:
        overview         - Where Optuna fits in lab workflows
        coil_design      - IH coil + workpiece multi-objective design
        motor_topology   - 2-level: Optuna outer + SIMP inner
        inverse          - BH curve / Play model inverse identification
        wpt              - Compensation topology selection
        literature       - Finding similar BBO benchmarks in lab corpus
        all              - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "intro", "lab_overview"):
        return LAB_OVERVIEW
    if topic in ("coil_design", "coil", "ih", "induction_heating"):
        return APP_COIL_DESIGN
    if topic in ("motor_topology", "motor", "topology"):
        return APP_MOTOR_TOPOLOGY
    if topic in ("inverse", "inverse_problem", "hysteresis_fit",
                 "play_fit", "bh_fit"):
        return APP_INVERSE
    if topic in ("wpt", "wireless_power", "compensation"):
        return APP_WPT
    if topic in ("literature", "lit", "corpus"):
        return APP_LITERATURE
    if topic == "all":
        return "\n\n".join([LAB_OVERVIEW, APP_COIL_DESIGN,
                            APP_MOTOR_TOPOLOGY, APP_INVERSE,
                            APP_WPT, APP_LITERATURE])
    return (f"Unknown topic '{topic}'. Available: overview, "
            "coil_design, motor_topology, inverse, wpt, literature, all.")
