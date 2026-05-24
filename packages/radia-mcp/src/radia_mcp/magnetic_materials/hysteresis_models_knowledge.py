"""Hysteresis model catalog & decision tree.

Covers 13+ models documented in the lab library
(W:/03_文献・論文/00_電磁界解析/磁気特性/ヒステリシス/ + subfolders):
  Jiles-Atherton, Play, Stop, Energy-Based (Henrotte/Egger), Preisach,
  Chua, Chan, Bouc-Wen, ES (E&S), Lee (李), Potter-Schmulian, Zirka,
  LLG (micromagnetics), Inverse-Distribution Function method,
  Cellular Automaton, Thermal extension, Bouc-Wen friction.

★ **LAB CORE METHOD**: B-input Stop model based energy formulation.
  See topic 'lab_core' for the recommended starting point.
"""


LAB_CORE_METHOD = r"""
# ★ Sugahara Lab CORE method: B-input Stop model based energy formulation

This is the **PRIMARY** hysteresis method used in current lab work
(2026-05).  Start here if you need hysteresis-aware FE in any lab
project (motor, IH, accelerator magnet, transformer).

## What it is (in one sentence)

A vector hysteresis model where:
  - **B (magnetic flux density) is the INPUT** (driven by FE A-formulation)
  - **Stop hysterons** (rate-independent operators with B as the
    "displacement") give H_k(B) per operator
  - Each operator wrapped in a **convex energy function** U_k → the
    constitutive law H = ∂U_total/∂B is thermodynamically consistent

## Why this specific combination?

| Feature | B-input | Stop hysterons | Energy wrapper |
|---------|---------|----------------|----------------|
| FE-natural input | ✓ A-formulation gives B | ✓ stop sees B directly | ✓ U(B) is convex |
| Iron loss correct | ✓ B-driven excitation matches motor reality | ✓ rate-independent → P_hyst per cycle | ✓ 2nd-law dissipation |
| FE convergence | ✓ symmetric system matrix | ✓ monotone operator | ✓ convex → global Newton |
| Vector hysteresis | ✓ B is the natural vector | ✓ stop is dual of Play with same vector extension | ✓ each U_k is 2D/3D natural |
| Calibration | ✓ measured B(H) → invert to H(B) | ✓ direct fit from major loop | ✓ monotone f_k constraint = physical |

**Alternatives and why they were NOT chosen**:
- Jiles-Atherton: 5-parameter, poor minor loops, parameter
  initialisation sensitive, NO vector extension
- Preisach: O(K²) cost, Everett surface calibration heavy, vector
  Preisach (Mayergoyz) only approximate
- Chua: limited minor loop accuracy
- LLG (micromagnetics): too expensive for device scale

**Play vs Stop**: identical mathematical structure (Bobbio-Bertotti
1997 duality), but the lab chose Stop because it pairs naturally with
B-input + energy:
  - Play hysteron: y = max(min(y_prev, x + η), x − η) — H-input mental
    model, M-output
  - Stop hysteron: y = max(−η, min(η, x − x_prev + y_prev)) — B-input
    mental model, H-output
  - In the energy wrapper, the convex U_k(r) with monotone
    f_k(r) = U_k'(r) makes the stop hysteron naturally satisfy:
    H_k(B) = f_k(B − x_k) where x_k is the state.

## Formulation

```
U_total(B) = sum_{k=1}^{K} U_k(|B − x_k|)         (convex energy)
H(B) = grad_B U_total = sum_k f_k(|B − x_k|) · (B − x_k)/|B − x_k|
state evolution: x_k = arg min over y in body of  (U_k(|B − y|) − η_k * y)
                       (Egger Schur complement Newton solve)
```

where:
- K = number of stop operators (typically 10-50, lab uses 20-30)
- η_k = stop thresholds in Tesla (log-spaced, e.g., 0.001 to 1.5)
- U_k(r) = ∫_0^r f_k(s) ds   (convex iff f_k monotone non-negative)
- f_k(r) = scalar shape function (1D in operator-coord r, used 3D-wise)

## Radia implementation

```python
import radia as rad
from radia.hysteresis_io import load_hys_file

K, eta, f_k_tables = load_hys_file('material.hys')

# The "B-input Stop model based energy model" = MatEnergyHysteresis
# with the conventional `f_k` shape functions
hys_mat = rad.MatEnergyHysteresis(K, eta, f_k_tables, eps=1e-6)

# Apply to a body
iron_obj = rad.ObjHexahedron(verts, [0, 0, 0])
rad.MatApl(iron_obj, hys_mat)

# Solve quasi-static step, commit state, repeat
result = rad.Solve(iron_obj, 1e-4, 1000, 0)
rad.MatHysCommitState(hys_mat)   # commit current step's state
# ... next time step ...
```

**Note**: Radia's `MatEnergyHysteresis` class implements the
*energy-based hysteresis* formulation.  Whether the internal hysterons
are "play-like" or "stop-like" is a matter of variable choice — the
SHAPE FUNCTIONS f_k in the .hys file encode the lab's specific
Stop-operator convention.

## Calibration workflow (.hys file)

The lab's standard procedure to generate a .hys file from measured
B(H):
1. Measure B(H) major loop on Epstein frame or ring core (DC bias OK)
2. Choose K (typically 20-30) and log-space η_k from 0.001 to ~B_sat
3. Run LSQ fit (lab tool: `decompose_stop.m` or Python equivalent) to
   determine f_k(r) per operator from the measured major loop
4. Verify by re-evaluating H(B) on the same major loop — should
   reproduce within < 1% error
5. Save as 2-column-per-operator .hys file

**Calibration constraint**: For energy wrapping, the f_k must be:
  - **non-negative** (for convexity of U_k)
  - **monotone non-decreasing in r** (for U_k convexity)
  - **continuous** (no jumps; smoothing typically required at small η_k)

If the raw fit gives negative or non-monotone f_k, the energy wrapping
fails and the FE solver becomes unstable.  The lab's calibration
script enforces these constraints with a soft penalty.

## Key references (lab library)

In `W:/磁気特性/ヒステリシス/02_Play_model/` (the historical Play
naming; Stop operators sit in the same family):
- **北尾純士 PhD thesis 「ヒステリシス特性を考慮した有限要素磁界解析
   の実用化に関する研究」** (Osaka U, 102 MB) — full development of
  the lab's B-input + energy formulation
- **161224 博士論文公聴会 v4** — PhD defense slides (297 MB, public
  reference for the method)
- 11_ヒステリシス解析2_40.pptx — internal training slides
- Models_of_magnetic_hysteresis_based_on_play_and_stop_hysterons.pdf
- Application_of_stop_and_play_models_to_the_representation_of_
  magnetic_characteristics_of_silicon_steel_sheet.pdf
- A General Approach to Hysteresis.pdf (Bobbio 1997 — the seminal
  paper on play / stop unification)
- Equivalent_Circuit_Modeling_of_DC_and_AC_Ferrite_Magnetic_
  Properties_Using_H-Input_and_B-Input_Play_Models.pdf

In `W:/磁気特性/ヒステリシス/12_Energy_Based/`:
- Energy-Based Hysteresis Model for Ferromagnetic Materials and
  Efficient Finite Element Formulation (Henrotte 2013)
- Energy-Based Magnetic Hysteresis Models (Jacques 2018, 254p
  monograph)

In `W:/磁気特性/ヒステリシス/03_Stop_model/`:
- ストップモデルとプレイモデルによるヒステリシス特性表現に関する検討

## Cross-reference

- MCP `motor_henrotte_lineage('energy_hys_2006')` — Henrotte energy
  hysteresis lineage (Bergqvist 1997 → Henrotte 2006 → Jacques 2018)
- MCP `magnetic_materials_hysteresis('play')` — Play model (B-input,
  same FE-friendly direction, alternative formulation)
- MCP `magnetic_materials_hysteresis('energy_based')` — Energy-based
  formulation theory
- MCP `magnetic_materials_radia_status('overview')` — MatEnergyHysteresis
  Python API
- Radia source: `src/radia/hysteresis_io.py` (`.hys` file format),
  `src/core/rad_mat_energy_hysteresis.cpp` (Egger Schur solver)

## When to deviate from the lab core method

| Scenario | Use this instead |
|----------|------------------|
| Permanent magnet FORC analysis | Preisach via external tool |
| Magnetic recording head, thin film | Potter-Schmulian |
| Domain wall dynamics, nm scale | LLG micromagnetics |
| SPICE circuit-coupled simulation | Chan (Nonlinear Transformer model) |
| Piezo-magnetic actuator | Bouc-Wen multiphysics |
| Linear DC analysis only (no hysteresis) | MatSatIsoTab |

For most engineering FE work on motors, transformers, IH workpieces,
electromagnets, **start with the lab core method**.
"""


CATALOG_OVERVIEW = r"""
# Hysteresis model catalog (W:/磁気特性/ヒステリシス/)

★ **Lab CORE method**: B-input Stop model based energy formulation.
  Use `magnetic_materials_hysteresis('lab_core')` for the production
  recipe.  The catalog below is the REFERENCE LANDSCAPE; for actual
  lab work, start with lab_core.

## Why so many models?

Magnetic hysteresis is empirically rich and not derivable from first
principles (except micromagnetics LLG, which is too expensive for FE
device analysis).  Each model captures DIFFERENT subset of features
with DIFFERENT computational cost.  Choosing the right one is a
trade-off between accuracy, calibration cost, and FE compatibility.

## 13 models in the lab library

| # | Model | Input | Output | Major loop | Minor loop | Vector | FE-friendly | Lab folder |
|---|-------|-------|--------|------------|------------|--------|-------------|------------|
| 1 | **Jiles-Atherton (JA)** | H | M | ◎ | ○ (modified) | ✗ | ○ | 01_Jiles_Atherton |
| 2 | **Preisach** | H | B | ◎ | ◎ | ✗ scalar | △ slow | Preisach_モデル |
| 3 | **Play (B-input)** | B | H | ◎ | ◎ | ○ | ◎ | 02_Play_model |
| 4 | **Stop (H-input dual)** | H | B | ◎ | ◎ | ○ | ◎ | 03_Stop_model |
| 5 | **Energy-Based (Henrotte/Bergqvist)** | B | H | ◎ | ◎ | ◎ | ◎ | 12_Energy_Based |
| 6 | **Chua-type** | dM/dt | M | ○ | △ (limited) | ✗ | △ | Chua型 |
| 7 | **Chan (Nonlinear Transformer)** | H | M | ◎ | ○ | ✗ | ○ SPICE | Chanモデル |
| 8 | **Bouc-Wen** | (force) | (disp) | ◎ piezo | ○ | ✗ | △ | Bouc-Wenモデル |
| 9 | **E&S (Enokizono-Soda)** | H_x, H_y | M_x, M_y | ◎ | ◎ | ◎ 2D | △ | ESモデル |
| 10 | **Lee (李) model** | H | M | ○ | △ | ✗ | △ | 李モデル |
| 11 | **Potter-Schmulian** | H | M | ◎ recording | △ | ✗ | × | Potter_Schmulianモデル |
| 12 | **Zirka extended JA** | H | M | ◎ | ◎ improved | ✗ | ○ | Zirka論文 |
| 13 | **LLG (micromagnetics)** | H_eff | M | ◎ | ◎ | ◎ | × too small | LLG |

Adjacent:
- Cellular Automaton (`セルラー・オートマトン`): coarse-grained alternative for Preisach
- Inverse-distribution-function method (`逆分布関数法`): numerical alternative for Preisach inversion
- Thermal_model: thermal-extension wrapper for any base model
- Friction model (`ヒステリシスの摩擦モデル`): unified play/stop interpretation

## Reference papers (lab)

- `Mathematical_Models_of_Hysteresis_invited.pdf` (Mayergoyz invited review)
- `Review_of_Hysteresis_Models_for_Magnetic_Materials.pdf` (top-down comparison)
- `Model_and_Simulations_of_Hysteresis_in_Magnetic_Cores.pdf`
- `Formulation_of_the_Everett_function_using_least_square_method.pdf` (Preisach calibration)
- `Iron_Loss_Estimation_Method_for_a_General_Hysteresis_Loop_With_Minor_Loops.pdf`

## Textbooks (W:/磁気特性/00_教科書/)

- **Bertotti, "Hysteresis in Magnetism"** (579 MB scan; the canonical 1998 monograph)
- **Mathematical Models of Hysteresis and Their Applications** (Mayergoyz; 260 MB)
- **Magnetic Hysteresis** (Della Torre + appendices on Play/Stop)
"""


DECISION_TREE = r"""
# Decision tree: which hysteresis model to use?

★ **DEFAULT for lab work**: B-input Stop + Energy formulation
  (`lab_core` topic).  Below tree only if you have a SPECIFIC reason
  to deviate.

```
Need: magnetic hysteresis in FE solver
│
├── DEFAULT (motor / IH / transformer / electromagnet):
│   └── Lab core method: B-input Stop + Energy → MatEnergyHysteresis
│        (see magnetic_materials_hysteresis('lab_core'))
│
├── Speed CRITICAL, scalar approximation OK, no rotational losses?
│   └── Use Play model (B-input, O(K) forward eval; in Radia as
│       MatPlayHysteresis). Vector mode also supported by directional
│       superposition.
│
├── Need TRUE vector hysteresis (rotational losses, motor 2D)?
│   ├── 2D only:  E&S model (Enokizono-Soda)
│   └── 2D / 3D:  Energy-Based Henrotte/Bergqvist
│       In Radia: MatEnergyHysteresis (B-input, Egger Schur Newton)
│
├── Need MOST ACCURATE minor loops (e.g., DC-bias on AC excitation)?
│   ├── Preisach (gold standard for scalar)
│   │     - Cost: O(K²) Everett function evaluation
│   │     - Calibration: needs full FORC measurement
│   │     - Inversion (H→B): non-trivial, use 逆分布関数法
│   └── Energy-Based: cleaner FE integration, similar minor-loop fidelity
│
├── Need GROUP-DELAY (frequency-domain hysteretic eddy effects)?
│   └── Complex permeability mu = mu' - j mu''
│       (NOT a true hysteresis model — see iron_loss_models topic)
│
├── Permanent magnet demagnetization (PM rotor at high T)?
│   ├── Use MatMagCurve (skeleton in Radia; full demag = TODO)
│   ├── Linear approximation: MatLin(mu_recoil) + Br offset
│   └── Reference: 99 永久磁石/FORC解析と永久磁石材料への適用.pdf
│
├── Electromagnetic launcher / SPICE-coupled circuit?
│   └── Chan model (designed for SPICE)
│
├── Piezoelectric / mechanical-electric coupling?
│   └── Bouc-Wen (multiphysics, both EM and mechanical hysteresis)
│
├── Magnetic recording head / thin film?
│   └── Potter-Schmulian (designed for thin-film media)
│
├── Need FULL micromagnetics (domain wall dynamics)?
│   └── LLG (Landau-Lifshitz-Gilbert)
│       - Too expensive for >μm scale
│       - See radia.lab.LLG examples (research only)
│
└── DC only, no hysteresis needed?
    └── MatSatIsoTab (B-H curve, no minor loops)
```

## Calibration data needed

| Model | Min measurements | Typical |
|-------|------------------|---------|
| MatSatIsoTab (DC) | DC initial mag curve | 5-10 (B, H) points |
| Jiles-Atherton | 5 params (Ms, k, c, alpha, a) | 1 major loop fit |
| Play / Stop | K play operators with f_k shape functions | 1 major loop + decomposition |
| Preisach | Everett surface | FORC measurement (~50 reversals) |
| Energy-Based | K play operators + non-negative shape | 1 major loop + monotone fit |
| E&S 2D | Vector measurements | RSST or 2D-SST device |
| LLG | M_s, exchange, anisotropy | Material constants |

## Standard reference: Mayergoyz "Mathematical Models of Hysteresis"

The 1991/2003 monograph (W:/磁気特性/00_教科書/Mathematical_Models...pdf,
260 MB) is the definitive overview.  Key chapters:
- Ch 1-2: Preisach formalism
- Ch 3-4: Vector Preisach
- Ch 5: Eddy-current hysteresis coupling
- Ch 7-8: JA, Stoner-Wohlfarth, magnetic recording

For energy-based hysteresis specifically, the Jacques 2018 PhD
monograph (254p, in S:/lab library, also covered in MCP
`motor_henrotte_lineage('jacques_2018')`) is the canonical reference.
"""


JA_MODEL = r"""
# Jiles-Atherton (JA) model

## History
- D. C. Jiles, D. L. Atherton, J. Magn. Magn. Mater. 61:48-60, 1986.
- Generalised to vectors: Bergqvist 1996, Sablik-Jiles thermal 1993.
- Modified JA for minor loops: Carpenter 1991, Calvert 2010.

## Formulation

5 parameters (Ms, k, c, alpha, a) define the anhysteretic
magnetisation M_an plus the reversible/irreversible split:
```
M_an(H_e) = M_s * [coth(H_e/a) - a/H_e]            (Langevin)
H_e = H + alpha * M                                  (effective field)
dM_irr/dH = (M_an - M_irr) / (delta * k - alpha * (M_an - M_irr))
M_rev = c * (M_an - M_irr)
M = M_irr + M_rev
```

where `delta = sign(dH/dt)` selects ascending/descending branch.

## Parameter ranges (typical Fe-Si)
- Ms ≈ 1.7 MA/m
- a ≈ 200-500 A/m
- alpha ≈ 1e-4 - 1e-3
- k ≈ 50-200 A/m
- c ≈ 0.05 - 0.3

## Pros
- Few parameters, physically interpretable
- Fast forward evaluation
- SPICE-compatible (Chan / Chua variants)

## Cons
- Poor minor-loop accuracy (without modification)
- NO vector formulation in original; Bergqvist 1996 vector extension
  is non-trivial
- Calibration sensitive to initialisation in fitting

## Modified JA for minor loops

`ヒステリシス/01_Jiles_Atherton/Minor loops modelling with a modified
Jiles-Atherton model and comparison with the Preisach model.pdf` —
applies an artificial "k(B)" or "delta" modulator to fix minor-loop
closure.  Common production fix in commercial codes (JMAG, FEMM).

## Radia status
**NOT implemented**.  Recommended path:
1. Use `MatSatIsoTab` for DC analysis
2. Use `MatPlayHysteresis` for AC analysis (better FE behaviour than JA)
3. Use Energy-Based via `MatEnergyHysteresis` for vector AC

JA was the lab-historical first choice but Play / Energy-Based have
won on FE convergence and calibration robustness.
"""


PLAY_MODEL = r"""
# Play model (B-input)

## Status: PRODUCTION in Radia

`rad.MatPlayHysteresis(K, eta, f_k_tables)` — B-input play model.

## Formulation (Bobbio et al. 1997 / Bertotti 1998 / Henrotte 2014)

Decompose H(B) into K play operators in parallel:
```
H(B, history) = sum_{k=1}^{K} f_k(B - state_k)
where state_k = B - eta_k * sign(B - state_k^prev),  if |B - state_k^prev| > eta_k
              = state_k^prev,                          otherwise
```
- K: number of play operators (typically 10-50)
- eta_k: play thresholds in Tesla (typically log-spaced 0.001 to 1.5)
- f_k(r): shape functions (per-operator force-displacement curves)

## Calibration (.hys file format)

Radia's hysteresis_io.py:
```
# K eta_1 eta_2 ... eta_K
# r_1  f_1(r_1)  f_2(r_1)  ...  f_K(r_1)
# r_2  f_1(r_2)  ...
# ...
```

Load via `K, eta, f_k = load_hys_file(path)`.

Procedure to calibrate:
1. Measure major loop B(H) up/down branches.
2. Decompose using LSQ fit to play operators.
3. K=20-30 usually sufficient for engineering accuracy.

## Pros
- **O(K) forward evaluation** (much faster than Preisach O(K²))
- **Vector extension natural** (project on K orthogonal axes)
- Direct B-input matches FE A-formulation (where B is the natural output)
- Shape functions can be NEGATIVE (no constraint), so easy to fit
- Convergent FE coupling (Henrotte's "energy-consistent" property)

## Cons
- Less physically interpretable than JA (it's a decomposition, not a model)
- Calibration data file (.hys) can grow with K

## Lab references

Major paper trail in `W:/磁気特性/ヒステリシス/02_Play_model/`:
- Bobbio 1997 (original)
- Henrotte 2014 (energy-consistent FE)
- 北尾純士 PhD thesis (Osaka U, 102 MB) — production FE implementation
- Application of Internal Reaction Field paper
- A General Approach to Hysteresis (review)

## Stop model (H-input dual)

The Stop model is the Legendre dual of Play with H input:
```
B(H) = sum_{k=1}^{K} g_k(H - state_k)   with stop thresholds in A/m
```
Use Stop when H is the natural input (e.g., Omega-reduced scalar
potential formulation).  Use Play when B is the natural input
(A-formulation).
"""


ENERGY_BASED = r"""
# Energy-Based hysteresis (Henrotte / Bergqvist)

## Status: PRODUCTION in Radia

`rad.MatEnergyHysteresis(K, eta, f_k_tables, eps=1e-6)` — energy-based
B-input, Egger Schur complement Newton solver.

## History
- Bergqvist 1997 IEEE TMAG 33(5):4109 (original energy formulation)
- Henrotte 2006 J. Appl. Phys. 99:08D505 (vector extension)
- Henrotte 2013 IEEE TMAG 49(5):2353 (variational FE formulation)
- Egger thesis (Linz) — Schur complement solver
- Jacques 2018 ULiège PhD (254p, S: library) — canonical reference

## Formulation

Decompose H(B) via K convex energy functions U_k:
```
U_total(B) = sum_k U_k(|B - x_k|)             (energy)
H(B, x_k) = grad_B U_total                    (constitutive)
x_k = arg min over x in body of (U_k(|B - x|) - eta_k * x)   (state evolution)
```

The energy functions U_k must be **convex** (non-negative, monotone f_k
derivatives).  This guarantees:
- Existence and uniqueness of state x_k at each time step
- Energy dissipation per cycle (2nd law of thermodynamics satisfied)
- FE convergence (constitutive operator is monotone)

## Pros
- **Thermodynamically consistent** (Bergqvist 1997)
- Vector hysteresis NATURAL (each U_k can be 2D / 3D)
- Better FE coupling than Play (the constitutive operator is monotone)
- Convex U_k allows GLOBAL Newton convergence

## Cons
- Per-eval cost: K independent 1D Newton solves (slower than Play's O(K))
- Calibration constrained (U_k must be convex → f_k monotone)
- Less data-fit flexibility than Play

## Comparison: Play vs Energy-Based

| Feature | Play | Energy-Based |
|---------|------|--------------|
| Forward (B→H) | O(K) direct | K Newton solves |
| Inverse (H→B) | Newton + analytical Jac | K independent Newton |
| Vector | Project on axes (approx) | Native (full vector U_k) |
| Calibration | Any shape (incl. negative) | Convex U_k required |
| Speed | 4-9 μs / eval | 100-500 μs / eval |
| FE convergence | Empirical (works) | Theoretical (monotone) |

For motor analysis: **Play is faster, Energy-Based is more reliable**
for nonlinear iteration on coarse meshes.

## Lab references

`W:/磁気特性/ヒステリシス/12_Energy_Based/`:
- Energy-Based Hysteresis Model for Ferromagnetic Materials and
  Efficient Finite Element Formulation (Henrotte 2013)
- Energy-Based Magnetic Hysteresis Models (Jacques 2018 monograph)
- Masters_Thesis_Alexander_Sauseng

Cross-ref MCP: `motor_henrotte_lineage('jacques_2018')` and
`motor_henrotte_lineage('energy_hys_2006')`.
"""


PREISACH = r"""
# Preisach model

## History
- F. Preisach, Z. Physik 94:277-302, 1935 ("Uber die magnetische
  nachwirkung" — IN W:/磁気特性/ヒステリシス/Preisach_モデル/)
- I. D. Mayergoyz formalized as Everett surface representation (1986)
- Vector extension: classical Preisach is SCALAR; Mayergoyz vector
  Preisach 1986; Della Torre moving Preisach 1990.

## Formulation (scalar, Mayergoyz form)

B(t) is the integral of "hysterons" (rectangular elementary loops) over
the Preisach plane (α, β):
```
B(t) = integral integral over (alpha, beta) plane of
       gamma_{alpha,beta}(H(t)) * mu(alpha, beta) dalpha dbeta
```
where:
- `gamma_{alpha,beta}` = elementary hysteron (B=±1 step)
- `mu(alpha, beta)` = Preisach distribution density
- The plane is restricted to alpha >= beta (geometric constraint)

The Everett function E(alpha, beta) = ∫∫ mu over a triangle gives
direct evaluation:
```
B(t) = -E(alpha_max, beta_min) + 2 * sum_k E(alpha_k, beta_k)
```
where (alpha_k, beta_k) is the staircase of return points.

## Calibration
- Measure FORC (First-Order Reversal Curves) at ~50 reversal points
- Compute Everett surface by second mixed partial derivative
- Lab paper: `Formulation_of_the_Everett_function_using_least_square_
  method.pdf`

## Pros
- **Best minor-loop accuracy** of all engineering hysteresis models
- Thermodynamically equivalent to Energy-Based (under certain conditions)
- Independence from initial state (closure property)

## Cons
- O(K²) cost vs O(K) for Play
- Calibration data heavy (FORC, not just major loop)
- Vector extension non-trivial (Mayergoyz vector Preisach is approximate)
- Inverse (H → B) requires solving for staircase numerically — see
  `逆分布関数法/Problems_in_practical_finite_element_analysis_using_
   Preisach_hysteresis_model.pdf`

## Radia status
**NOT implemented as Mat class**.  Workaround: use Energy-Based
(`MatEnergyHysteresis`) which gives similar minor-loop fidelity with
better FE convergence.

If true Preisach needed for permanent magnet FORC analysis:
- Use external Preisach library (e.g., Python `pyrose`)
- Compute B(H) trajectories offline
- Tabulate as MatSatIsoTab approximation
"""


VECTOR_HYSTERESIS = r"""
# Vector hysteresis (E&S, vector Preisach, vector Energy-Based)

## When needed
- Rotating electric machines: B rotates in stator/rotor laminations
  → rotational hysteresis loss (often 30-50% larger than alternating)
- 2D round-rotor analysis: B trajectory traces ellipses, not lines
- 3D end-region of transformers, motor end-windings

## E&S (Enokizono-Soda) model

`ヒステリシス/ESモデル/E＆Sモデルによる2次元磁気特性のヒステリシス
モデリング.pdf`

Splits 2D magnetization into:
```
M_x(t) = f_x(H_x(t)) + cross-coupling term
M_y(t) = f_y(H_y(t)) + cross-coupling term
```
Each f_x, f_y is a scalar hysteresis model (often Play / Preisach).
Cross-coupling captures rotational losses.

## Vector Energy-Based

Each play operator becomes 2D / 3D:
```
U_k(B_vec - x_k_vec) = U_k(|B - x_k|)
x_k_vec evolves in 2D / 3D plane
```
For motor stator iron, this naturally captures the elliptical B
trajectory.  No cross-coupling parameters needed (rotational loss
emerges from vector geometry).

## Rotational SST measurement

The lab has rotational single sheet tester (RSST) papers in
02_Play_model:
- `RSST_Field homogeneity in a two-phase round rotational single
  sheet tester...`
- `RSST_誘導電動機固定子を用いた電磁鋼板2次元磁気特性計測...`
- `RSST_圧電アクチュエータにより印加された応力下における電磁鋼板の
  ベクトル磁気ヒステリシス特性の測定`

These provide calibration data for vector models (alternating + circular
+ elliptical B trajectories).

## Decision

For motor analysis with rotational losses:
1. **Default**: scalar Play with rotational-loss correction factor
   (typical 1.3-1.7× alternating loss)
2. **Better**: Vector Energy-Based (`MatEnergyHysteresis` 2D mode —
   actually not yet implemented in Radia, TODO)
3. **Research-grade**: Full E&S model with RSST calibration

For transformer / inductor (uni-directional B), scalar Play is
sufficient.
"""


def get_hysteresis_models_knowledge(topic: str = "lab_core") -> str:
    """Dispatch by topic.

    Topics:
        lab_core       - ★ Sugahara Lab CORE: B-input Stop+Energy
                          (start here for any hysteresis FE work)
        catalog        - 13-model overview table + paper references
        decision_tree  - When to use which model (decision tree)
        jiles_atherton - JA model formulation + parameters
        play           - Play model (PRODUCTION in Radia)
        energy_based   - Energy-Based formulation theory
        preisach       - Preisach formalism (Everett function, FORC)
        vector         - Vector hysteresis (E&S, vector Preisach, RSST)
        all            - Everything
    """
    topic = topic.lower().strip()
    if topic in ("lab_core", "core", "lab_core_method", "stop_energy",
                  "b_input_stop"):
        return LAB_CORE_METHOD
    if topic in ("catalog", "overview", "13_models"):
        return CATALOG_OVERVIEW
    if topic in ("decision_tree", "decision", "which", "tree"):
        return DECISION_TREE
    if topic in ("jiles_atherton", "ja", "jiles"):
        return JA_MODEL
    if topic in ("play", "play_model"):
        return PLAY_MODEL
    if topic in ("energy_based", "energy", "henrotte", "bergqvist"):
        return ENERGY_BASED
    if topic in ("preisach", "preisach_model", "everett"):
        return PREISACH
    if topic in ("vector", "vector_hysteresis", "es", "rotational"):
        return VECTOR_HYSTERESIS
    if topic == "all":
        return "\n\n".join([
            LAB_CORE_METHOD, CATALOG_OVERVIEW, DECISION_TREE, JA_MODEL,
            PLAY_MODEL, ENERGY_BASED, PREISACH, VECTOR_HYSTERESIS,
        ])
    return (f"Unknown topic '{topic}'. Available: lab_core, catalog, "
            "decision_tree, jiles_atherton, play, energy_based, preisach, "
            "vector, all.")
