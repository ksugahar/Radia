"""Radia magnetic material implementation status: what is in MatXxx today,
how to use it, what is NOT yet implemented (TODO).

Cross-reference for users wanting to translate hysteresis_models /
iron_loss / silicon_steel / permanent_magnet knowledge into running code.
"""

OVERVIEW = r"""
# Radia magnetic material implementation status

★ **Lab core method = MatEnergyHysteresis (B-input Stop + Energy)**.
  For any new hysteresis-aware FE work, start with this.  See
  `magnetic_materials_hysteresis('lab_core')` for the formulation.

| Material class | Status | Use case | Code path |
|----------------|--------|----------|-----------|
| `rad.MatLin(mu_r)` | ◎ production | Linear isotropic | C++ |
| `rad.MatLin([mu_par,mu_perp], axis)` | ◎ production | Linear anisotropic (GO steel) | C++ |
| `rad.MatSatIsoTab(BH_data)` | ◎ production | Nonlinear DC BH curve | C++ |
| `rad.MatPlayHysteresis(K, eta, f_k)` | ◎ production | B-input Play hysteresis (alternative) | C++ |
| **`rad.MatEnergyHysteresis(K, eta, f_k, eps)`** | **★ production CORE** | **B-input Stop + Energy (lab primary)** | C++ |
| `rad.MatPM(Br, Hc, axis)` / direct `ObjHexahedron(verts, M)` | ◎ production | Permanent magnet (fixed M) | C++ |

## State management (Play / Energy)

```python
rad.MatApl(iron_obj, hysteresis_mat)
# ... solve quasi-static step ...
state = rad.MatHysSaveState(hysteresis_mat)   # save (n_elem × K × 9 array)
rad.MatHysRestoreState(hysteresis_mat, state)  # restore
rad.MatHysCommitState(hysteresis_mat)          # commit converged step
```

## .hys file format (Radia native)

```
# header line: K eta_1 eta_2 ... eta_K (in Tesla)
20 0.001 0.002 0.005 0.010 0.020 0.050 ... 1.500
# data rows: r f_1(r) f_2(r) ... f_K(r)
0.000  0.0  0.0  0.0  ...
0.001  10   8    6    ...
0.002  18   14   10   ...
...
2.000  1.7e5 1.5e5 ...
```

Load via:
```python
from radia.hysteresis_io import load_hys_file
K, eta, f_k_tables = load_hys_file('material.hys')
```

## Permanent magnet (no Solve needed for fixed M)

For NdFeB N42 (Br = 1.28 T → M = Br/mu_0 ≈ 1.02 MA/m):
```python
import radia as rad
# Rectangular PM, magnetized in +x
verts = [(0,0,0), (0.01,0,0), (0.01,0.01,0), (0,0.01,0),
         (0,0,0.005), (0.01,0,0.005), (0.01,0.01,0.005), (0,0.01,0.005)]
M_x = 1.28 / (4 * 3.14159 * 1e-7)  # ~1018592 A/m
pm = rad.ObjHexahedron(verts, [M_x, 0, 0])
B = rad.Fld(pm, 'b', [0.005, 0.005, 0.010])  # exterior B at offset
```

## PM + soft iron interaction (use Solve)

```python
pm = rad.ObjHexahedron(pm_verts, [0, 0, 954930])
iron = rad.ObjHexahedron(iron_verts, [0, 0, 0])
rad.MatApl(iron, rad.MatLin(1000))   # linear soft iron
assembly = rad.ObjCnt([pm, iron])
result = rad.Solve(assembly, 1e-4, 1000, 0)  # LU solver
B = rad.Fld(assembly, 'b', [0, 0, 0.05])
```
"""


HOW_TO_CHOOSE = r"""
# Decision tree: pick a Mat class for your scenario

```
Scenario
│
├── Linear isotropic (μ_r constant)
│   └── rad.MatLin(mu_r)
│
├── Linear anisotropic (GO steel, oriented permeability)
│   └── rad.MatLin([mu_par, mu_perp], axis_unit_vec)
│
├── Nonlinear BH curve (DC analysis only)
│   └── rad.MatSatIsoTab(BH_data)
│        where BH_data = [[H1,B1], [H2,B2], ...] in (A/m, T)
│
├── AC analysis with hysteresis loss
│   ├── B-input (A-formulation, FE-friendly)
│   │   └── rad.MatPlayHysteresis(K, eta, f_k_tables)
│   ├── Vector hysteresis / rotational losses
│   │   └── rad.MatEnergyHysteresis(K, eta, f_k_tables)
│   │        (energy-consistent; better for 2D / 3D)
│   └── H-input (Omega-formulation)
│       └── Stop model (TODO; not yet in Radia)
│
├── Permanent magnet (fixed M)
│   └── rad.ObjHexahedron(verts, [Mx, My, Mz])  -- no Mat needed
│        M in A/m, convert from Br via M = Br/mu_0
│
├── Permanent magnet + demagnetization knee
│   ├── Approximation: rad.MatLin(mu_recoil) + offset M
│   └── Full: TODO (planned via MatPM; the MatMagCurve skeleton was removed 2026-06-26)
│
├── Effective material for laminated stack (μ_eff complex per freq)
│   └── radia.calc_motor_lamination (Hollaus MSFEM, see motor_hollaus_eddy)
│        Not via Mat class — via panel script
│
└── Time-domain transient with hysteresis (Picard / Newton)
    └── radia.calc_motor_transient + rad.MatPlayHysteresis state mgmt
        (commit_state after each converged time step)
```

## Common pitfalls

| Pitfall | Symptom | Fix |
|---------|---------|-----|
| MatLin called with anisotropic [μ_par, μ_perp] but axis missing | crash | Pass axis as 2nd arg |
| MatSatIsoTab missing H=0 point | wrong incremental μ at small H | Include [0,0] row |
| Play hysteresis without committed state | state ≠ closed loop | rad.MatHysCommitState after each step |
| Permanent magnet specified as [Br_x, Br_y, Br_z] in Tesla | wrong M (off by μ_0) | Use M in A/m, NOT Br in T |
| Wrong material assigned to body (rad.MatApl missed) | linear default mu_r=1 | Verify rad.MatGetMat() returns expected handle |
"""


TODO_LIST = r"""
# Radia magnetic material TODO list (research-grade additions)

## High priority
1. **Jiles-Atherton Mat class** (`MatJilesAtherton(Ms, k, c, alpha, a)`)
   — many users have JA-fitted data from JMAG / Femtet
2. **Stop model** (`MatStopHysteresis`) — H-input for Omega-formulation
3. **Vector Energy-Based** (2D / 3D vector mode for MatEnergyHysteresis)
   — currently only scalar superposition

## Medium priority
4. **Preisach with Everett interpolation** — for PM FORC analysis
5. **Vector Preisach (Mayergoyz)** — full 2D rotational hysteresis
6. **PM full demag (via MatPM)** — knee detection + irreversible demag tracking (MatMagCurve skeleton removed 2026-06-26)
7. **Thermal extension wrapper** — wrap any Mat with T-dependent params

## Low priority
8. **Bouc-Wen for piezo-magnetic coupling**
9. **LLG micromagnetics** (very small scale, research only)
10. **E&S 2D model** — RSST-calibrated rotational loss
11. **Cellular Automaton hysteresis** — Mayergoyz coarse-grained alternative

## Iron loss post-process TODO
12. **Bertotti 3-term post** (P_iron = k_h f B^n + k_c (Bf)² + k_a (Bf)^1.5)
13. **MSE for non-sinusoidal** (Reinert 2001)
14. **iGSE for asymmetric PWM** (Venkatachalam 2002)

## Adjacent tools
15. **BH curve fitting helper** (fit Frohlich / Brauer / piecewise-cubic)
16. **.hys file generation from major loop measurement** (Play decomposition)
17. **FORC analysis pipeline** (raw measurement → Everett surface)
18. **Steel grade database lookup tool** (JIS code → constants)

Items 1-3 are the highest-value next additions.  Item 17 (FORC) would
enable proper PM characterisation (currently the lab does this in
external Mathematica notebooks).

## Cross-references

- `motor_henrotte_lineage` MCP — Energy-based theoretical foundations
- `motor_hollaus_eddy` MCP — Effective material for laminated stacks
- `peec_carstensen_ac_loss` MCP — AC copper + per-layer iron loss
- `magnetic_materials_*` MCP (this subpackage) — what to read for each topic
"""


def get_radia_status_knowledge(topic: str = "overview") -> str:
    """Dispatch by topic.

    Topics:
        overview         - Mat class table + state mgmt + .hys format
        how_to_choose    - Decision tree for picking a Mat class
        todo_list        - Research-grade additions still missing
        all              - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "status", "mat_class"):
        return OVERVIEW
    if topic in ("how_to_choose", "choose", "decision", "decision_tree"):
        return HOW_TO_CHOOSE
    if topic in ("todo_list", "todo", "future", "missing"):
        return TODO_LIST
    if topic == "all":
        return "\n\n".join([OVERVIEW, HOW_TO_CHOOSE, TODO_LIST])
    return (f"Unknown topic '{topic}'. Available: overview, how_to_choose, "
            "todo_list, all.")
