"""Hysteresis model catalog & decision tree.

Covers 13+ models documented in the lab library at
``W:/03_文献・論文/00_電磁界解析/30_磁気特性/02_ヒステリシスモデル/``:
  Jiles-Atherton (01_), Preisach (02_), Play (03_), Stop (04_),
  Energy-Based (05_; Henrotte / Bergqvist / Francois-Lavet / Egger),
  Chua (06_), Chan (07_), Bouc-Wen (08_), Potter-Schmulian (09_),
  ES = E&S Enokizono-Soda (10_), Lee 李 (11_), Zirka (12_),
  LLG micromagnetics (13_), Inverse-Distribution-Function (80_),
  Thermal extension (81_), research-group surveys (90_).
  Reviews + cellular automaton + friction model in 00_review/.
  Textbooks in 30_磁気特性/00_教科書/ (Bertotti, Mayergoyz, Della
  Torre, Visintin "Differential Models of Hysteresis" / Applied
  Mathematical Sciences, Hashimoto motor-drive magnetic materials).

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

In `W:/03_文献・論文/00_電磁界解析/30_磁気特性/02_ヒステリシスモデル/
03_Play_model/` (the historical Play naming; Stop operators sit
in the same family; this folder has the largest paper collection,
1000+ files):
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

In `.../02_ヒステリシスモデル/05_Energy_Based/`:
- 2006_Henrotte_EnergyBased.pdf — Henrotte/Nicolet/Hameyer COMPEL
  25:71-80, foundational vector model with spring+friction analogy
- 2013_Energy-based_VariationalModel_FrancoisLavetHenrotte.pdf
  — Francois-Lavet/Henrotte/Stainier/Noels/Geuzaine J Comp Appl Math
  246:243-250, fully-implicit time-incremental version
- Energy-Based Hysteresis Model for Ferromagnetic Materials and
  Efficient Finite Element Formulation (Henrotte 2013 IEEE TMAG)
- Energy-BasedMagnetic Hysteresis Models.pdf (Jacques 2018 ULiège
  PhD, 254p monograph — canonical reference)
- Efficient evaluation of forward and inverse energy-based magnetic
  hysteresis operators.pdf (Egger/Engertsberger/Schafelner 2025
  arXiv:2507.15289 — **regularization + Schur complement → O(K)
  inverse operator**, see `efficient_eval` topic)
- On_Forward_and_Inverse_Energy-Based_Magnetic_Vector_Hysteresis_
  Operators.pdf (Egger et al., IEEE Trans Magn 61:7300207, 2025)
- Masters_Thesis_Alexander_Sauseng.pdf (Linz, rotational losses)
- COMPUMAG2023_poster_Domenig.pdf — 3-D scalar potential FE
  formulation incorporating energy-based vector hysteresis

In `.../02_ヒステリシスモデル/04_Stop_model/`:
- ストップモデルとプレイモデルによるヒステリシス特性表現に関する検討

In `.../02_ヒステリシスモデル/00_review/`:
- Mathematical Models of Hysteresis_invited.pdf (Mayergoyz)
- Model and Simulations of Hysteresis in Magnetic Cores.pdf
- Review_of_Hysteresis_Models_for_Magnetic_Materials.pdf
- Iron_Loss_Estimation_Method_for_a_General_Hysteresis_Loop_With_
  Minor_Loops.pdf (Taitoda/Takahashi/Fujiwara 2015 IEEE TMAG
  51:8112304 — DC-biased coefficients + apparent frequency, see
  also `iron_loss_knowledge` topic)
- Formulation_of_the_Everett_function_using_least_square_method.pdf
- ヒステリシスの摩擦モデルとパラメータ.pdf (friction interpretation)
- セルラー・オートマトンによる磁化特性の表現.pdf

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
# Hysteresis model catalog
# (W:/03_文献・論文/00_電磁界解析/30_磁気特性/02_ヒステリシスモデル/)

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
| 2 | **Preisach** | H | B | ◎ | ◎ | ✗ scalar | △ slow | 02_Preisach |
| 3 | **Play (B-input)** | B | H | ◎ | ◎ | ○ | ◎ | 03_Play_model |
| 4 | **Stop (H-input dual)** | H | B | ◎ | ◎ | ○ | ◎ | 04_Stop_model |
| 5 | **Energy-Based (Bergqvist/Henrotte/Francois-Lavet/Egger)** | B | H | ◎ | ◎ | ◎ | ◎ | 05_Energy_Based |
| 6 | **Chua-type** | dM/dt | M | ○ | △ (limited) | ✗ | △ | 06_Chua |
| 7 | **Chan (Nonlinear Transformer)** | H | M | ◎ | ○ | ✗ | ○ SPICE | 07_Chan |
| 8 | **Bouc-Wen** | (force) | (disp) | ◎ piezo | ○ | ✗ | △ | 08_Bouc_Wen |
| 9 | **Potter-Schmulian** | H | M | ◎ recording | △ | ✗ | × | 09_Potter_Schmulian |
| 10 | **E&S (Enokizono-Soda)** | H_x, H_y | M_x, M_y | ◎ | ◎ | ◎ 2D | △ | 10_ES |
| 11 | **Lee (李) model** | H | M | ○ | △ | ✗ | △ | 11_Lee |
| 12 | **Zirka extended JA** | H | M | ◎ | ◎ improved | ✗ | ○ | 12_Zirka |
| 13 | **LLG (micromagnetics)** | H_eff | M | ◎ | ◎ | ◎ | × too small | 13_LLG |

Adjacent:
- Inverse-distribution-function method (80_Inverse_Distribution):
  numerical alternative for Preisach inversion
- Thermal extension (81_Thermal_extension): temperature-dependent wrapper
- Cellular Automaton (in 00_review/): coarse-grained alternative for Preisach
- Friction model (in 00_review/ヒステリシスの摩擦モデルとパラメータ.pdf):
  unified play/stop interpretation
- Research-group surveys (90_research_groups, 21 files)

## Reference papers (00_review/)

- `Mathematical Models of Hysteresis_invited.pdf` (Mayergoyz invited)
- `Review_of_Hysteresis_Models_for_Magnetic_Materials.pdf` (top-down)
- `Model and Simulations of Hysteresis in Magnetic Cores.pdf`
- `Formulation_of_the_Everett_function_using_least_square_method.pdf`
  (Preisach calibration)
- `Iron_Loss_Estimation_Method_for_a_General_Hysteresis_Loop_With_
   Minor_Loops.pdf` (Taitoda 2015 — DC-biased + apparent frequency)
- `ヒステリシスの摩擦モデルとパラメータ.pdf`
- `セルラー・オートマトンによる磁化特性の表現.pdf`
- `高周波大振幅励磁におけるMn-Znフェライトコア内部の磁化特性解析.pdf`

## Textbooks (00_教科書/, in W:/03_文献・論文/00_電磁界解析/30_磁気特性/)

| File | Description |
|------|-------------|
| Hysteresis_in_magnetism.pdf (579 MB) | **Bertotti 1998** canonical monograph (full + chunked) |
| Magnetic_Hysteresis.pdf (152 MB) | **Della Torre** full text |
| Magnetic_Hysteresis_*.pdf (8 chunks) | Della Torre per-chapter (Preisach Model, Preisach Applications, Moving / Product Model, Irreversible / Reversible, Physics of Magnetism, Frontmatter, Appendix A: Play and Stop Models, Appendix B: LogNormal) |
| Mathematical Models of Hysteresis and Their Applications.pdf (266 MB) | **Mayergoyz** + applications |
| Mathematical Models of Hysteresis.pdf (750 KB) | Mayergoyz abstract / TOC excerpt |
| Applied Mathematical Sciences.pdf (34 MB) | **Visintin "Differential Models of Hysteresis"** (Springer AMS series 111, 1994) — Preisach mathematical theory |
| Hysteresis_in_magnetism.pdf | Bertotti chapter-by-chapter scan |
| General properties of power losses in soft ferromagnetic material.pdf | Bertotti 1988 power-loss decomposition foundation |
| magnetic-material-for-motor-drive-systems-fusion-technology-of-electromagnetic-fields...pdf (24 MB) | **Hashimoto et al. 2019** Springer (978-981-32-9905-4) — fusion technology for motor magnetic materials |
| 電磁鋼板の磁気特性と取扱方法.pdf (97 MB) | Silicon steel sheet handling — Japanese practical reference |
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
│   ├── Full PM demag = TODO (planned via MatPM; MatMagCurve skeleton removed 2026-06-26)
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

## Modified JA for minor loops (Benabou 2008, JMMM 320:e1034)

The classical JA model **cannot close minor loops** — the
magnetisation variation rate is too high on the return branch
(point 2 → 3 in their Fig 1), so the loop does not close.

**Benabou 2008 fix** (simple, no prior knowledge of H evolution
needed — unlike other approaches):
- Minor loops are mainly the result of **irreversible** changes
- → INCREASE the irreversible parameter `k` and DECREASE the
   reversible parameter `c` during the minor-loop return branch

```
k = k_ini + Δk / (1 + (H/H_0)^γ)        (sigmoid, smooth — eq 7)
c = c_ini − Δc                          (step — eq 8)
```

- Sigmoid for k (abrupt change would give non-physical
   discontinuity in M); step for c (no noticeable impact if abrupt)
- `H_0` = inflexion point = (H_1 + H_2)/2 (the two turning points)
- `γ` = slope, set to impose k continuity at loop closure
- Applied ONLY on the return branch (point 2 → 3); behaviour from
   1 → 2 is already correct
- Minor-loop turning points detected via the **directional
   parameter δ** sign change (δ = +1 ascending, −1 descending)

**Concrete parameters** (Fe-Si 3% non-oriented 0.5 mm sheet,
1 Hz with 3rd harmonic, Benabou 2008 Table 1):
- M_s = 1.48e6 A/m, k_ini = 65.5 A/m, c_ini = 0.316,
   a = 71.33 A/m, α = 1.47e-4, Δk = 20.0 A/m, Δc = 0.218

**Accuracy vs Preisach** (Benabou 2008 Table 2, loss error):
| Case | Modified JA | Preisach |
|------|-------------|----------|
| Without minor loops | < 1% | 5% |
| With minor loops | 2.4% | 2.9% |

Modified JA matches Preisach accuracy on losses while keeping JA's
simplicity (only 5+2 params vs full Everett surface) and adding NO
computational overhead.  Recommended over full Preisach when only
**loss estimation** (not detailed minor-loop shape) is needed.

Inverse (B-input) JA for A-formulation FE (Sadowski 2002,
Benabou eq 5-6): dM/dB form with effective B_e = μ_0 H_e.

**Lab note**: JA + Benabou minor-loop fix is a viable LIGHTWEIGHT
alternative to the lab core Energy-based method when (a) only iron
loss matters, (b) vector hysteresis is NOT needed, (c) calibration
data is limited to a few symmetric loops.  For vector / motor work,
still prefer the Energy-based method (`lab_core_production`).

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

Major paper trail in `W:/03_文献・論文/00_電磁界解析/30_磁気特性/
02_ヒステリシスモデル/03_Play_model/` (1000+ files, the lab's
densest collection):
- Bobbio 1997 (original)
- Henrotte 2014 (energy-consistent FE)
- 北尾純士 PhD thesis (Osaka U, 102 MB) — production FE implementation
- Application of Internal Reaction Field paper
- A General Approach to Hysteresis (review)
- Plus RSST measurement papers, motor application papers, ESIM
  coupling papers, etc.

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
# Energy-Based hysteresis (Bergqvist / Henrotte / Francois-Lavet / Egger)

## Status: PRODUCTION in Radia

`rad.MatEnergyHysteresis(K, eta, f_k_tables, eps=1e-6)` — energy-based
B-input, Egger Schur complement Newton solver.

## Lineage (foundational papers in 05_Energy_Based/)

| Year | Authors | Contribution |
|------|---------|--------------|
| 1997 | Bergqvist | Original energy formulation, dry-friction-like pinning, Physica B 233:342-347 |
| 2006 | **Henrotte/Nicolet/Hameyer** | Vector model with h = h_r + h_i (spring + friction analogy), COMPEL 25:71-80 — see "Henrotte 2006 mechanical analogy" below |
| 2013 | **Francois-Lavet/Henrotte/Stainier/Noels/Geuzaine** | Fully-implicit time-incremental variational form, J Comp Appl Math 246:243-250 — this is the version solved in Radia |
| 2015 | Jacques/Sabariego/Geuzaine/Gyselinck | Direct + inverse energy-consistent hysteresis in dual FE formulations, IEEE TMAG 52:7300304 |
| 2016 | Prigozhin/Sokolovsky/Barrett/Zirka | Duality approach for unregularized vector model, IEEE TMAG 52:1-11 |
| 2018 | **Jacques** | ULiège PhD monograph (254p) — canonical reference |
| 2022 | Sauseng/Domenig/Roppert/Kaltenbacher | Rotational losses corrections, CEFC 2022 |
| 2024 | Domenig/Roppert/Kaltenbacher | 3-D scalar potential FE incorporation, IEEE TMAG |
| 2025 | **Egger/Engertsberger/Schafelner** | **Regularization + Schur complement → O(K) inverse operator** (arXiv:2507.15289, see `efficient_eval` topic) |
| 2025 | Egger/Engertsberger/Domenig/Roppert/Kaltenbacher | Semi-smooth Newton method, arXiv:2506.01499 |

## Henrotte 2006 mechanical analogy (key intuition)

Decompose the applied field into reversible + irreversible parts:
```
h = h_r + h_i      with   |h_i| <= κ   (subgradient G, "grey circle")
```
- **h_r** = reversible spring force (nonlinear spring), drives the
  anhysteretic curve M = M(h_r) via Langevin / atan / saturating fn
- **h_i** = irreversible friction force, bounded by pinning strength κ
  (= Bloch wall pinning by impurities/inclusions)
- Memory effect originates from the NON-DIFFERENTIABILITY of the
  dissipation functional ρ̇^Q = κ|Ṁ|: the subgradient is a whole set,
  not a single gradient → h_i is multi-valued inside the circle.

Update rule (simplified, efficient):
```
if |h^{n+1} - h_r^n| > κ:
    h_r^{n+1} = h^{n+1} - κ * (h^{n+1} - h_r^n) / |h^{n+1} - h_r^n|
```
This naturally extends to **VECTOR hysteresis** without further
assumption (`h, h_r, h_i, M` all 2-D / 3-D vectors).

Combined model with n+1 fractions for better fit:
```
h = Σ_k ω^k h_r^k + Σ_k ω^k h_i^k       (k = 0, ..., n)
where κ^0 = 0 (reversible Bloch wall bending fraction)
```
Identification via measured coercivity staircase:
```
h_coer(h_peak) = Σ_{k=0}^{m(h_peak)} ω^k κ^k  /  Σ_{k=0}^{m(h_peak)} ω^k
```
where m(h_peak) is the highest k with κ^k < h_peak.  Lab uses 5+
fractions for engineering accuracy.

## Modern variational form (Francois-Lavet 2013, fully-implicit)

Internal polarizations J = Σ_k J_k defined by **convex non-smooth**
minimization at each time step:
```
J_k = argmin_J  U_k(J) − ⟨H, J⟩ + χ_k |J − J_{k,p}|       (forward)
J_k = argmin    (ν_0/2)|B − ΣJ_k|² + Σ U_k(J_k) + χ_k |J_k - J_{k,p}|
                                                            (inverse)
B = μ_0 H + Σ_k J_k       (forward)
H = ν_0 (B − Σ_k J_k)     (inverse)
```
Internal energy density (Jacques 2018 form):
```
U_k(J) = −(2 A_{s,k} J_{s,k} / π) · log( cos( (π/2) · J / J_{s,k} ) )
```
The pinning strengths χ_k ≥ 0 play the role of κ^k above.

## Pros
- **Thermodynamically consistent** (Bergqvist 1997)
- Vector hysteresis NATURAL (each U_k can be 2D / 3D)
- Better FE coupling than Play (the constitutive operator is monotone)
- Convex U_k allows GLOBAL Newton convergence
- Provably convergent FE solvers (Egger 2025 semi-smooth Newton)

## Cons
- Naive per-eval cost: K independent 1D Newton solves (forward)
  or full K×K dense system (inverse, naive)
- Inverse non-trivial because polarizations are coupled via ν_0|B−ΣJ_k|²
- Calibration constrained (U_k must be convex → f_k monotone)
- Less data-fit flexibility than Play

**These cons are largely solved by Egger 2025** — see `efficient_eval`.

## Comparison: Play vs Energy-Based (updated 2025)

| Feature | Play | Energy-Based |
|---------|------|--------------|
| Forward (B→H) | O(K) direct | O(K) Newton (Egger 2025) |
| Inverse (H→B) | Newton + analytical Jac | **O(K) Schur complement (Egger 2025)** |
| Vector | Project on axes (approx) | Native (full vector U_k) |
| Calibration | Any shape (incl. negative) | Convex U_k required |
| Speed (K=20) | 4-9 μs / eval | ~10 ms / 500-step trajectory (Egger 2025) |
| FE convergence | Empirical (works) | Theoretical (monotone, provably convergent) |

For motor analysis: with Egger 2025 efficient implementations,
**Energy-Based is now competitive in speed** with Play AND retains
its theoretical advantages.

## Lab references

`W:/03_文献・論文/00_電磁界解析/30_磁気特性/02_ヒステリシスモデル/
05_Energy_Based/`:
- 2006_Henrotte_EnergyBased.pdf
- 2013_Energy-based_VariationalModel_FrancoisLavetHenrotte.pdf
- Energy-Based Hysteresis Model for Ferromagnetic Materials and
  Efficient Finite Element Formulation.pdf
- Energy-BasedMagnetic Hysteresis Models.pdf (Jacques 2018 monograph)
- Efficient evaluation of forward and inverse energy-based magnetic
  hysteresis operators.pdf (Egger/Engertsberger/Schafelner 2025)
- On_Forward_and_Inverse_Energy-Based_Magnetic_Vector_Hysteresis_
  Operators.pdf (Egger et al. IEEE TMAG 61:7300207, 2025)
- Masters_Thesis_Alexander_Sauseng.pdf (rotational losses, Linz)
- COMPUMAG2023_poster_Domenig.pdf (3-D scalar potential)
- Contribution_202_a.pdf

Cross-ref MCP: `motor_henrotte_lineage('jacques_2018')` and
`motor_henrotte_lineage('energy_hys_2006')` and
`magnetic_materials_hysteresis('efficient_eval')`.
"""


EFFICIENT_EVAL = r"""
# Efficient evaluation of forward and inverse Energy-Based operators
# (Egger / Engertsberger / Schafelner 2025, arXiv:2507.15289)

## Why this matters

The Francois-Lavet 2013 energy-based variational model gives implicit
hysteresis through convex but **non-smooth** minimization problems
at every material point.  Two pain points historically:

1. The non-smooth norm |J − J_p| forbids standard Newton solvers
   (subgradient methods required, e.g. Prigozhin 2016 duality).
2. The inverse operator H = H(B; J_p) was naively O(K³) per eval
   because the K partial polarizations are COUPLED through the
   shared term (ν_0/2)|B − ΣJ_k|².

Egger 2025 fixes BOTH issues with a clean recipe.

## The fix: regularization + Schur complement

### Step 1 — Regularize the norm

```
|x|_ε := sqrt(|x|² + ε),   ε > 0
```
- For ε = 0: recovers the original norm
- For ε > 0: infinitely differentiable everywhere (no kink at x=0)
- Error bound: |J − J^ε| = O(√ε), in practice **observed O(ε)**
- ε = 1e−6 already gives engineering accuracy

### Step 2 — Damped Newton on regularized problem

Forward operator (decoupled):
```
J_k^ε = argmin_J  U_k(J) − ⟨H, J⟩ + χ_k |J − J_{k,p}|_ε
```
- F_ε(J) is strongly convex with parameter γ > 0 (γ depends on
  internal energy parameters A_{s,k}, J_{s,k})
- Standard damped Newton + Armijo backtracking (σ=0.1, ρ=0.5)
  → global + locally quadratic convergence (O(1) iterations)
- Cost: **O(K)** for K partial polarizations
- Tolerance tol = ε (matched to regularization)

Inverse operator (coupled — needs more care):
```
{J_k^ε} = argmin_{J_1, ..., J_K}  (ν_0/2) |B − ΣJ_k|² + Σ_k [U_k(J_k) + χ_k |J_k − J_{k,p}|_ε]
H = ν_0 (B − Σ_k J_k^ε)
```
The shared (ν_0/2)|B − ΣJ_k|² couples all K polarizations.

### Step 3 — Schur complement trick

The Newton system has block structure:
```
[ ν_0 H_K     E_K^T  ]  [ δJ̃  ]   [ −∇G_ε(J̃) ]
[ E_K        −I      ]  [ δJ   ] = [   0        ]
```
- H_K = block-diagonal Hessian (K × K blocks of size d×d, d ≤ 3)
- E_K = K-fold column-wise concatenation of d×d identity I
- δJ̃ = Newton update for all polarizations (length Kd)
- δJ = Σ_ℓ δJ_ℓ = total update (length d, auxiliary)

Eliminate δJ̃ via block diagonality of H_K:
```
δJ̃ = H_K^(−1) (b_K − E_K^T δJ)
```
Substitute → small d×d Schur complement system:
```
[ E_K H_K^(−1) E_K^T + I ] δJ = E_K H_K^(−1) b_K
```
- Setup of Schur system: **O(K)** (H_K is block-diagonal)
- Solve d×d Schur system: **O(1)** (d ≤ 3)
- Back-substitute for δJ̃: **O(K)**

→ **Inverse operator now O(K) per evaluation**, same complexity as
   forward.  Newton iteration count is O(1).

## Benchmark numbers (from Egger 2025 Table I)

500 time-steps, AMD Ryzen 5 PRO 5675U @ 4.3 GHz, uniform excitation:

| K | Prigozhin duality | Forward (regularized Newton) | Inverse, standard | Inverse, Schur |
|---|-------------------|------------------------------|-------------------|----------------|
| 5 | 0 ms | 2 ms | 6 ms | **2 ms** |
| 10 | 1 ms | 5 ms | 14 ms | **5 ms** |
| 20 | 2 ms | 10 ms | 40 ms | **10 ms** |
| 50 | 5 ms | 27 ms | 289 ms | **27 ms** |
| 100 | 11 ms | 54 ms | 1071 ms | **59 ms** |

→ Schur trick recovers forward-operator speed for the inverse.
   Standard inversion is ~20× slower at K=100.

## Implications for Radia

`rad.MatEnergyHysteresis(K, eta, f_k_tables, eps=1e-6)`:
- The `eps=1e-6` default is exactly the regularization parameter ε.
- The C++ implementation (`rad_mat_energy_hysteresis.cpp`) already
  uses Egger Schur complement for the inverse (this is the "Egger
  Schur complement Newton solver" mentioned in `energy_based` topic).
- Choice of K: 10-30 typical for engineering; cost scales linearly.

## Jacobian access (for outer Newton in FE solver)

Egger 2025 §V points out: regularized smoothness ⇒ standard implicit
function theorem applies ⇒ ∇_H B(H; J_p) and ∇_B H(B; J_p) can be
computed by formal differentiation of the optimality KKT system.
This unblocks **outer Newton on the magnetic FE problem** (without
nested Picard) — see arXiv:2506.01499 and arXiv:2506.01499 follow-up.

## Reference

- H. Egger, F. Engertsberger, A. Schafelner, "Efficient evaluation
  of forward and inverse energy-based magnetic hysteresis operators",
  arXiv:2507.15289v1 [math.NA] 21 Jul 2025.
- File: `W:/03_文献・論文/00_電磁界解析/30_磁気特性/02_ヒステリシスモデル/
   05_Energy_Based/Efficient evaluation of forward and inverse
   energy-based magnetic hysteresis operators.pdf`

Cross-ref: `energy_based`, `lab_core`, `motor_henrotte_lineage`.
"""


PREISACH = r"""
# Preisach model

## History
- F. Preisach, Z. Physik 94:277-302, 1935 ("Uber die magnetische
  nachwirkung" — IN W:/03_文献・論文/00_電磁界解析/30_磁気特性/
  02_ヒステリシスモデル/02_Preisach/)
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

## Mathematical foundations textbook

Visintin's **"Differential Models of Hysteresis"** (Applied
Mathematical Sciences vol 111, Springer 1994) provides the
rigorous mathematical theory for Preisach as a hysteresis
*operator* (variational inequality / monotone operator
framework).  File: `W:/03_文献・論文/00_電磁界解析/30_磁気特性/
00_教科書/Applied Mathematical Sciences.pdf`.  Use when
proving well-posedness of Preisach + PDE coupled systems.

## Preisach extensions — the 2 limitations + their fixes
## (Mayergoyz "Mathematical Models of Hysteresis" 2003 + Della Torre
##  "Magnetic Hysteresis" 1999 — the two canonical Preisach monographs)

Classical Preisach has **two properties that real magnetic materials
VIOLATE** (Mayergoyz Ch 1.3 representation theorem proves these are
the *necessary + sufficient* conditions for classical Preisach):

1. **Congruency property**: all minor loops between the same two
   field extrema have IDENTICAL shape regardless of history.  Real
   steel minor loops depend on the dc-bias position (cf.
   `stop_vs_play_silicon_steel` "equal vertical chords" — same issue).
2. **Deletion (wiping-out) property**: closing a minor loop erases
   its memory completely in ONE traversal.  Real materials show
   **accommodation** — several traversals needed to reach the stable
   minor loop.

Each limitation has a named fix in the monographs:

### Congruency fixes (Mayergoyz Ch 2, Della Torre Ch 4)

| Extension | Idea | Source |
|-----------|------|--------|
| **Moving model** | effective field H_eff = H + αM (mean-field interaction shifts the Preisach plane) — breaks congruency | Mayergoyz §2.1, Della Torre Ch 4 |
| **Input-dependent measure** | Preisach density μ(α,β,H) depends on current input — Matsuo's input-dependent shape function (`input_dependent_shape`) is the play/stop analog | Mayergoyz §2.2 |
| **Product model** | state-dependent measure as a product of two factors | Della Torre Ch 4 |
| **DOK / CMH reversible** | add magnetization-dependent (DOK) / state-dependent (CMH) REVERSIBLE component for accurate susceptibility | Della Torre Ch 3 |

### Deletion fix (Mayergoyz §2.5, Della Torre Ch 5)

- **Accommodation model** (a.k.a. aftereffect / reptation): the
  minor loop drifts over repeated cycles before stabilizing.  Removes
  the deletion property.  Essential for accurate sub-loop loss when
  the excitation is not symmetric.

### Vector Preisach (Mayergoyz Ch 3)

The classical Preisach is SCALAR.  Mayergoyz's vector generalization
superposes scalar Preisach operators over all directions:
```
M(H) = ∫_{|e|=1} e · P̂[e·H] dΩ_e      (integral over the unit sphere)
```
- Isotropic identification (§3.4): from a single scalar measurement
- Anisotropic (§3.5): needs the group-of-rotations machinery
- This is the rigorous alternative to E&S (`vector` topic) and vector
  energy-based (`lab_core_production`) for rotational losses.

### Eddy-current + Preisach core losses (Mayergoyz Ch 6)

Couples the Preisach operator to Maxwell diffusion in a lamination.
§6.3 derives **excess hysteresis losses** from sharp magnetic
transitions — the Preisach-formalism counterpart of Bertotti's
statistical excess loss (`iron_loss('bertotti')`).  §6.5 handles
circularly-polarized fields (rotational eddy losses).

### Superconducting hysteresis (Mayergoyz Ch 5)

Striking result: the **Bean critical-state model** for type-II
superconductors is a SPECIAL CASE of the Preisach model.  Nonlinear
flux diffusion with linear/circular/elliptical polarization.

### Della Torre Appendix A — Play and Stop

Della Torre treats play/stop as the mechanical-backlash limit
(gear-train deadzone).  Note (Della Torre's own caution): the play
model output does NOT saturate naturally → must be fed through an
external saturating nonlinearity, which he argues limits its
usefulness vs Preisach.  (The lab's energy-based method sidesteps
this by wrapping the stop in a convex energy — `lab_core`.)

## Which Preisach variant for what

| Need | Variant |
|------|---------|
| Textbook scalar reference | Classical Preisach (Mayergoyz Ch 1) |
| dc-biased silicon steel minor loops | Moving model / input-dependent measure |
| Accommodation (drifting minor loops) | Accommodation model (Della Torre Ch 5) |
| Rotational motor losses | Vector Preisach (Mayergoyz Ch 3) — or vector energy-based |
| Superconductor magnetization | Bean = Preisach special case (Mayergoyz Ch 5) |
| FE production (lab) | Energy-based, NOT Preisach (better FE convergence + O(K) cost) |

For Radia's lab work the energy-based method (`lab_core_production`)
is preferred over any Preisach variant — Preisach's O(K²) cost,
FORC-heavy calibration, and approximate vector extension lose to the
energy-based model's O(K) cost, virgin-curve calibration, and native
vector form.  Preisach remains the GOLD STANDARD for scalar minor-loop
fidelity and the reference against which other models are validated.

## References (the two Preisach monographs)

- I. D. Mayergoyz, "Mathematical Models of Hysteresis and Their
  Applications", 2nd ed., Elsevier 2003 (492p).  6 chapters:
  classical / generalized scalar / vector / stochastic-viscosity /
  superconducting / eddy-current core losses.  File:
  `00_教科書/Mathematical Models of Hysteresis and Their Applications.pdf`.
- E. Della Torre, "Magnetic Hysteresis", IEEE Press 1999 (232p).
  Physics-first (not math-rigorous): classical Preisach + DOK/CMH
  reversible + moving/product (congruency) + accommodation (deletion)
  + vector + applications (eddy/magnetostriction/inverse) + Appendix A
  play/stop, Appendix B lognormal.  File: `00_教科書/Magnetic_Hysteresis.pdf`.

Cross-ref: `mathematical_foundations` (Visintin operator theory),
`input_dependent_shape` (play/stop analog of input-dependent measure),
`stop_vs_play_silicon_steel` (the congruency-violation problem),
`vector` (E&S alternative to vector Preisach), `iron_loss('bertotti')`
(excess-loss counterpart), `lab_core_production` (why the lab prefers
energy-based over Preisach).
"""


VECTOR_HYSTERESIS = r"""
# Vector hysteresis (E&S, vector Preisach, vector Energy-Based)

## When needed
- Rotating electric machines: B rotates in stator/rotor laminations
  → rotational hysteresis loss (often 30-50% larger than alternating)
- 2D round-rotor analysis: B trajectory traces ellipses, not lines
- 3D end-region of transformers, motor end-windings

## E&S (Enokizono-Soda) model — the reluctivity-tensor form

`ヒステリシス/10_ES/E＆Sモデルによる2次元磁気特性のヒステリシス
モデリング.pdf` (Soda & Enokizono, J. Magn. Soc. Japan 24:827, 2000,
Oita University).

E&S works in the **H = ν·B** direction (reluctivity), splitting H into
an **in-phase** reluctivity term + a **90°-lagging** hysteresis term:

```
H_x = ν_xr · B_x + ν_xi · dB_x/dτ
H_y = ν_yr · B_y + ν_yi · dB_y/dτ
```

- `ν_xr, ν_yr` = magnetic reluctivity coefficients (in-phase with B)
- `ν_xi, ν_yi` = magnetic hysteresis coefficients (in-phase with dB/dτ,
   i.e. 90° lead) — these encode the hysteresis loss
- All four ν are FUNCTIONS of (B_max, θ_B, axis-ratio α=B_min/B_max,
   phase τ=ωt) — determined from 2D vector measurements (RSST)

**Why the dB/dτ term**: a pure reluctivity tensor `H = [ν]B` (eq 1)
CANNOT represent hysteresis, because when B=0 it forces H=0 (no
remanence) and a point on increasing-B vs decreasing-B is
ambiguous with B alone.  Adding the dB/dτ derivative term lets the
model distinguish ascending/descending → captures both alternating
AND rotating-flux hysteresis.

Identification: apply sinusoidal (B_x, B_y) → measure (H_x, H_y),
extract ν_r from the in-phase Fourier component and ν_i from the
90°-lag component (up to 3rd harmonic typically).  Calibrated from
elliptical B-trajectory measurements at various axis ratios.

**E&SS variant**: replace dB/dτ with ∫B dτ (Matsuo 2014 slides).

This is the production 2D vector model for **motor rotational iron
loss** (Enokizono/Oita lineage).  For the lab, the energy-based
vector method (`lab_core_production`) is preferred (thermodynamically
consistent, no separate RSST calibration needed), but E&S remains the
measurement-driven reference for anisotropic 2D characterisation.

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


INPUT_DEPENDENT_SHAPE = r"""
# Input-dependent shape function for stop and play models
# (Matsuo & Shimasaki, Kyoto Univ; 04_Stop_model/ + 80_Inverse_Distribution/)

## Why input-dependent shape function?

A stop or play model with **input-INDEPENDENT** shape function
has the **property of equal vertical chords regardless of dc bias**:
for back-and-forth input variation of amplitude `a` centered on
arbitrary bias `B_0`, the vertical chord Δh(a, B_0) only depends
on `a`, NOT on `B_0`.

**Silicon steel sheet does NOT satisfy this property** -- vertical
chords are LARGE near `B = B_s` and SMALL near `B = 0`, so the
basic stop / play model gives ~5% representation error for
asymmetric (dc-biased) loops even after correct fit of symmetric
loops.

Ni-Zn ferrite ring DOES satisfy equal vertical chords → basic
stop / play works fine.

## The fix (Matsuo 2004, 2005): product form `g(η, s, B) = w(B) · g_0(η, s)`

Use a B-dependent **weighting function** `w(B)` multiplying the
input-independent shape `g_0`:

```
H = S(B) = ∫_0^{B_s} g(η, s_η(B), B) dη         (general)
        = w(B) · S_0(B)                          (product form)
```

where `S_0(B) = ∫_0^{B_s} g_0(η, s_η(B)) dη` is the input-independent
stop model.

This is **equivalent to Mayergoyz's nonlinear Preisach model**
(Matsuo & Shimasaki 2005, Representation Theorems): same necessary
and sufficient conditions (wiping-out + equal vertical chords for
the rescaled function H/w).

## Weighting function choices

Three practical forms, in increasing sophistication:

| `w(B)` | When useful | Reference |
|--------|-------------|-----------|
| `w(B) = 1` | Ni-Zn ferrite (equal vertical chords already) | baseline |
| `w(B) = Δh(B_s, B) = h_M^+(B) - h_M^-(B)` | Major-loop width; simple + effective | Matsuo 2004 |
| `w(B) = Δh(B_s, 0) · [Δh(B_s, B)/Δh(B_s, 0)]^c`, c=1.5 | Best fit for silicon steel | Matsuo 2005 |

The `c=1.5` form often outperforms `c=1.0` for grain-oriented
silicon steel along both rolling and transverse directions.

## Five identification methods (Matsuo 2004 stop model comparison)

| Method | Inputs used | Notes |
|--------|------------|-------|
| (a) Loop widths of symmetric loops only | symmetric | Simplest; least accurate |
| (b) Loop widths of asymmetric loops + symmetrization | asymmetric | Needs asymmetric measurements |
| (c) Descending curves of symmetric loops + symmetrization | symmetric | **Best if only symmetric data available** |
| (d) Least-squares over symmetric + asymmetric | both | **Most accurate overall** |
| (e) Least-squares over symmetric only | symmetric | Poor on asymmetric (overfits) |

Result for grain-oriented silicon steel (JIS 35p105, rolling
direction): method (d) achieves ~1% total error vs ~4% baseline.

## Play model variant (Matsuo & Shimasaki 2005, IEEE TMAG 41:3112)

For the play model, the weighting function is determined by a
**congruency property** constraint: the Everett function from
first-order reversal curves must match the Everett function from
the symmetric major-loop branches.  Least-squares minimization
over `w(B)` enforces this.

**Data correction trick** (Matsuo 2005, Section III.B):
Before identification, replace boundary values to ensure the
distribution function μ_{2m, M-m+1} stays NON-positive
(corresponds to physical irreversibility).  Specifically replace
h_S^+(a_{m-1}, B_{m-1}) by h_S^+(a_m, B_{m-1}) and the same on
the descending branch.  Without this correction, simulated loops
that were NOT used for identification can become wider than the
major loop -- a clear physical violation.

## Cross-reference with lab core

The lab's `MatPlayHysteresis` and `MatEnergyHysteresis` Radia
classes use the **input-INDEPENDENT** shape function by default
(simpler calibration).  For silicon steel where the dc-bias
asymmetry matters (motor stator with PWM-driven harmonics, DC-bias
choke), use the input-dependent extension by passing a
B-dependent `f_k_tables`.  Currently not exposed in the public
Radia API -- on the development roadmap.

## References (lab library)

`W:/03_文献・論文/00_電磁界解析/30_磁気特性/02_ヒステリシスモデル/`:
- `04_Stop_model/Stop Model With Input-Dependent Shape Function
   and Its Identification Methods.pdf` (Matsuo, Terada, Shimasaki,
   IEEE TMAG 40:1776-1783, 2004)
- `04_Stop_model/Representation Theorems for Stop and Play Models
   With Input-Dependent Shape Functions.pdf` (Matsuo & Shimasaki,
   IEEE TMAG 41:1548-1551, 2005)
- `80_Inverse_Distribution/An_identification_method_of_play_model_
   with_input-dependent_shape_function.pdf` (Matsuo & Shimasaki,
   IEEE TMAG 41:3112-3114, 2005)
- `04_Stop_model/Application_of_stop_and_play_models_to_the_
   representation_of_magnetic_characteristics_of_silicon_steel_
   sheet.pdf` (Matsuo, Shimode, Terada, Shimasaki, IEEE TMAG
   39:1361-1364, 2003)
"""


LAMINATION_HOMOGENIZATION = r"""
# Lamination homogenization: two-step pragmatic approach
# (Henrotte, Steentjes, Hameyer, Geuzaine, IET SMT 9:152-159, 2015)
# 81_Thermal_extension/Pragmatic two-step homogenisation technique
# for ferromagnetic laminated cores.pdf

## Problem

Electrical machines (motors, transformers) use **stacks of thin
ferromagnetic laminations** (~0.3 mm thickness) separated by
insulating coating to suppress eddy currents.  Direct 2D / 3D FE
simulation of each individual lamination is computationally
prohibitive -- you'd need to resolve the lamination thickness +
the skin depth (~50 μm at 50 Hz, much smaller at higher freq) in
a 3D macroscale device model with overall dimensions of ~100s of mm.

Need: a **HOMOGENIZED** material model that hides the laminated
structure from the macroscale solver but still captures eddy
currents + saturation + hysteresis correctly.

## Why HMM is too expensive (monolithic)

Heterogeneous Multiscale Method (HMM, E & Vanden-Eijnden 2007):
- Solve mesoscale 1D problem at EACH macroscale Gauss point
- Macroscale problem size × mesoscale problem size = huge

Niyonzima-Sabariego-Dular-Henrotte-Geuzaine 2013 implemented HMM
for ferromagnetic laminations but the macroscale model becomes
1-2 orders of magnitude larger than the non-homogenized one.
Unacceptable for desktop engineering use.

## The two-step pragmatic approach

**Step 1 -- Mesoscale 1D cross-lamination model**:
Solve 1D magnetodynamics with hysteresis across HALF the
lamination thickness (use symmetry: `curl h(0) × n = 0`).

```
∫_ω [ḃ(h, history) · h' + σ^{-1} curl h · curl h'] dω = 0,  ∀h'
```

with:
- `b(h, history)` = Bergqvist-Henrotte energy-based hysteresis
   (see `energy_based` topic)
- Typically 50 equidistant nodes, 360 time steps per period
- h-field formulation (h is the natural driving quantity)
- Power density: P = (e × h)(d) · n at the lamination surface

This 1D model is "exact" given the hysteresis material law.

**Step 2 -- Algebraic approximation by system identification**:

Build an algebraic `H(B, Ḃ, p_k)` that REPRODUCES the 1D
mesoscale response for a chosen excitation pattern, region by
region:

```
H(B, Ḃ, p_k) = (p_0 + p_1 · |B|^{2 p_2}) · B
               + (p_3 + p_4 / sqrt(p_5^2 + |Ḃ|^2)) · Ḃ
```

- `p_0, p_1, p_2`: polynomial anhysteretic saturation curve
   (reversible part h_r)
- `p_3`: eddy current term (~ σ d^2 from analytical thin-lamination
   solution)
- `p_4`: hysteresis amplitude
- `p_5`: regularization in dissipation term (prevents
   non-smoothness at Ḃ = 0; cf. Egger 2025 trick in
   `efficient_eval`)

Identification: `scipy.optimize.leastsq` on chosen excitation
signal (triangle, sinusoidal, or REAL waveform from the device
simulation at a chosen probe point).

## Region-by-region identification

In a motor (e.g. 4-pole PMSM):
- STATOR teeth: alternating B at line frequency + PWM harmonics
   → use real waveform at `P_1` (center of stator tooth)
- ROTOR yoke: B fluctuating around average value that depends on
   load → use real waveform at `P_2` (below PM)

Identified parameters (Henrotte 2015 Table 1, IEW JIS-50A270):

| Signal | p_0 | p_1 | p_2 | p_3 | p_4 | p_5 |
|--------|-----|-----|-----|-----|-----|-----|
| Triangular | 93.3 | 0.26 | 11.5 | 0.035 | 34.5 | 8.32 |
| Sinusoidal | 95.9 | 0.29 | 11.4 | 0.041 | 28.6 | 8.03 |
| Tooth waveform | 77.8 | 0.80 | 9.87 | 0.036 | 36.7 | 8.68 |
| Rotor waveform | 88.6 | 13.6 | 3.11 | 0.053 | 7.65 | 9.52 |

Stator and rotor parameter sets differ substantially -- region
distinction is ESSENTIAL.

## Comparison with HMM and other homogenization techniques

| Method | Model size | Accuracy | Implementation |
|--------|-----------|----------|----------------|
| Monolithic HMM (Niyonzima 2013) | 10-100× larger | High | Heavy code |
| Asymptotic homogenization (Krähenbühl 2004) | Same | Linear only | Needs periodicity |
| Algebraic (Giordano 2014) | Same | No hysteresis | Limited to energy-density potentials |
| **Two-step (Henrotte 2015)** | **Same** | **Good for engineering** | **Non-invasive (no FE coupling needed)** |

## Implementation in Radia

The 2-step approach can be implemented as:
1. Pre-process: run mesoscale 1D for each region (motor stator
   tooth, rotor back-iron) with representative waveform → fit
   `{p_0, ..., p_5}` per region
2. Macroscale FE: use the algebraic `H(B, Ḃ)` law as a Jacobian-
   compatible nonlinear constitutive relation

For Radia + NGSolve, this could be wrapped as a
`MatLaminationHomogenized(p_k_per_region)` material class.  Not
yet implemented -- on the roadmap for motor lamination support.

For now, the lab uses the `MatEnergyHysteresis` class directly on
the lamination (treating each lamination as a body) for small
test problems, but this is too expensive for full motor analysis.

## Cross-reference

- `energy_based` topic -- the underlying Bergqvist-Henrotte
   hysteresis model
- `efficient_eval` -- Egger 2025 regularization (same `sqrt(p_5^2
   + |Ḃ|^2)` smoothing pattern)
- MCP `motor_hollaus_eddy('overview')` -- Hollaus effective
   material approach (alternative homogenization)
- MCP `motor_henrotte_lineage('jacques_2018')` -- Jacques PhD
   monograph, full energy-based theory

## References

- F. Henrotte, S. Steentjes, K. Hameyer, C. Geuzaine,
   "Pragmatic two-step homogenisation technique for ferromagnetic
   laminated cores", IET Sci. Meas. Technol. 9(2):152-159, 2015
   (CEM 2014 special issue, DOI 10.1049/iet-smt.2014.0201)
- I. Niyonzima, R. V. Sabariego, P. Dular, F. Henrotte, C. Geuzaine,
   "Computational homogenization for laminated ferromagnetic cores
   in magnetodynamics", IEEE TMAG 49(5):2049-2052, 2013
- S. Steentjes, F. Henrotte, C. Geuzaine, K. Hameyer, "A dynamical
   energy-based hysteresis model for iron loss calculation in
   laminated cores", Int. J. Numer. Model. 27(3):433-443, 2013
"""


CELLULAR_AUTOMATON = r"""
# Cellular Automaton model of magnetic hysteresis
# (Miyasaka & Saito 兆古, Hosei University, SICE 2008)
# 02_ヒステリシスモデル/00_review/セルラー・オートマトンによる
#   磁化特性の表現.pdf

## Key insight

The classical Preisach model is mathematically EQUIVALENT to a
**2D cellular automaton** (CA) over the Preisach plane (α, β):
each Preisach hysteron `γ_{α,β}` is a cell that switches between
+1 (black) and -1 (gray) according to a local update rule driven
by the external H field.  The total magnetization is the SUM of
all cell states, weighted by the Preisach distribution.

This was implicit in Preisach's 1935 original paper but not
articulated as CA until Wolfram (1984) formalized cellular
automata and Saito's group (Hosei) connected the two.

## Three magnetization regions (classical domain theory)

The CA framework naturally decomposes the magnetization process
into 3 regions matching domain physics:

| Region | Mechanism | CA behavior |
|--------|-----------|-------------|
| 1. Reversible domain wall motion | Small H, walls oscillate | All cells reversible (no permanent state change) |
| 2. **Irreversible domain wall motion** | Medium H, walls snap past pinning | Snap transitions on rule threshold |
| 3. Rotation magnetization | High H (near saturation) | Final coherent flip of remaining cells |

The Miyasaka-Saito 2008 paper visualizes domains via **Bitter
method** (colloidal Fe3O4 suspension on the sample surface), then
binarizes the images to construct the CA initial state directly
from observed domain structure.

## The "rule" abstraction

A Wolfram-style local CA rule (e.g. Rule 30 for 1D, or a 2D
analog) maps the current cell state + neighborhood to the next
state.  For Preisach-type hysteresis CA, the rule is:

```
if H > α_cell and current_state == -1:  next_state = +1
if H < β_cell and current_state == +1:  next_state = -1
else:                                    next_state = current_state
```

where `(α_cell, β_cell)` are the cell's switching thresholds
(its position in the Preisach plane).

**Slow / fast regions** observed experimentally:
- Slow magnetization when |H| is small (few cells exceed
   threshold)
- Rapid magnetization in the "knee" region (many cells just below
   threshold flip together)
- Slow approach to saturation (few cells remain to flip)

CA naturally reproduces this S-curve magnetization characteristic
without needing a parameterized analytical model.

## Connection to other models

| Model | CA interpretation |
|-------|-------------------|
| Preisach | 2D CA on (α, β) plane, sum of cell states |
| Play hysterons | 1D chain of CAs (each cell = one play operator) |
| Energy-based | CA + convex energy regularization |
| LLG micromagnetics | CA on nm-scale mesh with continuous-time evolution |

## Why this matters for the lab

1. **Intuition for the Bergqvist-Henrotte friction analogy**:
   the "grey circle" (admissible h_i region) is the cell's
   threshold; switching ≡ snap past pinning ≡ CA transition.
2. **1/f fluctuation in magnetization** is a natural CA
   prediction (avalanche dynamics) -- explains Barkhausen noise.
3. **Sub-grid microstructure** can be captured by CA without
   explicit RVE meshing.
4. Saito 兆古 lab (Hosei) is the same group that developed the
   **Cauer Ladder Network (CLN)** reduction now in
   `radia_mcp.mor` -- cross-reference for nonlinear material
   reduction.

## CA-based Preisach implementation pseudocode

```python
import numpy as np

class PreisachCA:
    def __init__(self, N_alpha, N_beta, mu_density):
        self.cells = np.full((N_alpha, N_beta), -1)  # init -1
        self.alpha_grid = np.linspace(0, B_s, N_alpha)
        self.beta_grid  = np.linspace(-B_s, 0, N_beta)
        self.mu = mu_density  # Preisach distribution

    def step(self, H_external):
        for i, a in enumerate(self.alpha_grid):
            for j, b in enumerate(self.beta_grid):
                if H_external > a and self.cells[i,j] == -1:
                    self.cells[i,j] = +1
                elif H_external < b and self.cells[i,j] == +1:
                    self.cells[i,j] = -1

    def magnetization(self):
        return np.sum(self.cells * self.mu) * d_alpha * d_beta
```

Cost: O(N²) per step (same as standard Preisach).  Wavefront
acceleration techniques (Lazarus-Bardi 2012) reduce this to
O(N · staircase_length).

## References

- S. Miyasaka, Y. Saito, "Representation of Magnetization
   Characteristics by Cellular Automaton", 51st SICE Joint
   Conference, Yamagata University, Nov 2008, pp. 68-69
- T. Sunaga, Y. Saito, K. Horii, "Visualization of iron loss
   distribution in magnetic materials by Bitter method",
   Visualization Information Society Symposium 2007, C207
- S. Wolfram, "Cellular Automata", 1984 (Rule 30 + classification)
- Cross-ref MCP: `magnetic_materials_hysteresis('preisach')`
   for the analytical Preisach form;
   `mor_cln('overview')` for Saito lab's CLN reduction
"""


MATHEMATICAL_FOUNDATIONS = r"""
# Mathematical foundations of hysteresis operators
# (Krasnoselskii-Pokrovskii 1989 + Visintin 1994 — the two canonical
#  textbooks cited by all the engineering papers in the lab library)

## Why a math-foundations topic?

Engineering papers on Play, Stop, Energy-Based, and even Preisach models
all open with "let `P` be the play operator (cf. Krasnosel'skii-Pokrovskii
[155; Sect. 2])" or "the stop hysteron of Bobbio (Visintin Ch III.3)".
Engineers usually skip the proofs.  But three results from the
mathematical theory matter for **practical use**:

1. **Lipschitz continuity of Play/Stop in C^0([0,T])** — guarantees
   the FE Picard iteration converges (no anomalous "non-convergence
   due to play hysteron" failure modes).
2. **Hilpert's inequality** (Visintin Thm III.2.6) — guarantees
   uniqueness of solutions to PDEs with Play/Stop coupling (i.e.
   the lab core method does not have spurious multi-valued
   solutions).
3. **Equivalence theorem** (Krasnoselskii §35, Visintin Ch III.4) —
   the Prandtl-Ishlinskii model (parallel combination of plays) is
   equivalent to a continual system of hysterons (Ishlinskii's
   material) — explains why the Bergqvist-Henrotte K-fraction
   decomposition is the right discretization.

## Formal definition (Visintin Ch III.1)

A **hysteresis operator** F : Dom(F) ⊂ C^0([0,T]) × X → C^0([0,T])
satisfies:

| Property | Statement |
|----------|-----------|
| **Causality** (1.1) | F(u, ξ_0)(t) depends only on u\|_{[0,t]} and ξ_0 |
| **Rate independence** (1.2) | F(u ∘ φ, ξ_0) = F(u, ξ_0) ∘ φ for any continuous monotone time change φ : [0,T] → [0,T] |
| **Semigroup property** (1.4 or 1.5) | F can be restarted from any intermediate state — the state ξ(t) is a complete summary of the history |

The state ξ(t) generally takes values in a metric space X (not just R^2).
This generality matters for the Preisach model, where ξ = staircase of
reversal points.

## Generalized Play (Krasnoselskii §2, Visintin §III.2)

**Mechanical picture** (Krasnoselskii Fig 2.1):
piston-in-cylinder with gap h.  Input u = piston position,
output x = cylinder position.  State `(u, x)` lives in the strip
`u ≤ x ≤ u + h`.  Cylinder moves only when piston touches the wall.

**Generalized play** replaces the straight walls with two
nondecreasing continuous curves γ_r, γ_l with γ_r ≤ γ_l (Visintin
(2.7)).  Set J(σ) := [γ_r(σ), γ_l(σ)].  Then the play operator P
is defined by the **differential inclusion** (Visintin (2.8)):

```
ε̇ ∈ -∂I_J(σ)(ε)
```

where I_J(σ) is the indicator function of J(σ).  Equivalently, the
**variational inequality** (Visintin (2.9)):

```
ε ∈ J(σ),       ε̇·(ε - v) ≤ 0    for all v ∈ J(σ)
```

with initial condition `ε(0) = min{max{γ_r(σ(0)), ε^0}, γ_l(σ(0))}`
(projection onto the admissible set).

**Theorem (Lipschitz continuity, Visintin Thm 2.2 + 2.3)**:
If γ_r, γ_l are Lipschitz continuous, then P : C^0([0,T]) × R →
C^0([0,T]) is Lipschitz continuous.  P also maps W^{1,p}(0,T) →
W^{1,p}(0,T) for any p ∈ [1, ∞].

**Practical use**: this is the formal guarantee that Radia's
`MatPlayHysteresis` will not blow up in a transient solve as long
as the f_k shape functions in the .hys file are Lipschitz.

## Stop (Krasnoselskii §3, Visintin §III.3)

**Mechanical picture** (Visintin Ch III.3 Fig 4):
elastic-plastic spring E-P in series.  Input ε = total
deformation, output σ = stress.  σ saturates when |σ| reaches the
yield K.

**Stop operator** Q is the dual of play under the change of
variables `u := ε + σ, w := ε - σ`.  Variational inequality
(Visintin (3.11)):

```
σ ∈ K,       (σ̇ - ε̇)·(σ - v) ≤ 0    for all v ∈ K
```

with K = [A_l(ε), A_u(ε)] (Visintin (3.10)).

**Bobbio-Bertotti 1997 duality**: Play with input σ + output ε is
isomorphic to Stop with input ε + output σ under (u,w) = (ε+σ, ε-σ).
This is why the lab core method can use either — the choice depends
on whether the FE formulation gives H or B as natural input.

## Hilpert's inequality (Visintin Thm III.2.6 — KEY result)

For (σ_1, ε_1^0), (σ_2, ε_2^0) ∈ W^{1,1}(0,T) × R, set
ε_i := P(σ_i, ε_i^0) and δε := ε_1 - ε_2.  Let h(t) be a
measurable selection of the Heaviside graph H(σ_1 - σ_2).  Then:

```
d|δε|/dt ≤ (σ̇_1 - σ̇_2)·h         a.e. in ]0, T[
```

**Why this matters**: integrating in time over [0,t]:

```
|ε_1(t) - ε_2(t)| ≤ |ε_1^0 - ε_2^0| + ∫_0^t (σ̇_1 - σ̇_2)·h dτ
```

This is the **comparison principle** that gives:
- uniqueness of FE solutions for the coupled PDE+hysteresis problem
   (motor magnetic field with B-input Play model)
- continuous dependence on the input data
- convergence proofs for Picard / Newton outer iterations

Used in Visintin Ch VIII (semigroups) + Ch XI (PDE with discontinuous
hysteresis) to prove well-posedness of magnetic FE problems with
hysteresis.

## Prandtl-Ishlinskii Models (Visintin §III.4, Krasnoselskii §35)

Parallel combination of plays with different gaps:

```
σ = Σ_{k=1}^{K} f_k · P[γ_{l,k}, γ_{r,k}](σ, ε_k^0)
```

where each play has its own thresholds γ_{l,k}, γ_{r,k}.

**Limit K → ∞ = Ishlinskii's material** (Krasnoselskii §35):
continuum of hysterons indexed by yield threshold η ∈ [0, η_max],
with weight function w(η).  This is **mathematically equivalent
to the classical Preisach model** when w(η) is suitably chosen.

**Engineering implication**: the K=20-30 used in `MatPlayHysteresis`
is a finite-element-style discretization of an Ishlinskii continuum.
Doubling K halves the calibration error (similar to FE p-refinement).

## Vector extension (Visintin Thm 2.11, Krasnoselskii §16)

For K ⊂ R^N convex, the **vector play** is defined by the variational
inequality (Visintin (2.24)):

```
σ - ε ∈ K,    c(ε̇ - ε̇^free, ε - v - σ + ε^free) ≤ 0    ∀v ∈ K
```

**Theorem 2.11**: Vector play P(·, ε^0) : C^0([0,T]; R^N) → C^0([0,T]; R^N)
is well-defined and Lipschitz continuous.

This is the formal foundation for vector hysteresis in 2D / 3D motor
FE (E&S model, vector Energy-Based).  The constant c can be promoted
to a positive-definite tensor for anisotropic materials.

## Homogenization of the Preisach model (Visintin §XI.7)

Theorem: under a fast-oscillating Preisach input
`σ_ε(x, t) = σ̄(x, t) + σ_osc(x/ε, t/ε)`, the limit ε → 0 of the
output ε_ε(x, t) is governed by a **homogenized constitutive law**
of the form `ε̄(x, t) = P̄(σ̄)(x, t)` where P̄ is again a Preisach
operator but with a modified weight function depending on the
multi-scale geometry.

**Practical connection**: this rigorous result is the mathematical
underpinning for the **Henrotte-Steentjes 2015 two-step lamination
homogenization** approach (see `lamination_homogenization` topic).
The pragmatic 2-step is essentially a one-step Riemann-sum
approximation of Visintin's asymptotic homogenization.

## What this means for the lab

| Lab code/method | Mathematical foundation cite |
|-----------------|------------------------------|
| `MatPlayHysteresis` (Radia C++) | Krasnoselskii §2 + Visintin §III.2 |
| `MatEnergyHysteresis` (Radia C++) | Bergqvist 1997 + Visintin §III.4 (PI model) |
| Bergqvist friction analogy | Krasnoselskii §39 + Visintin §III.4 |
| Forward eval `H(B) = Σ f_k P_k(B)` | Prandtl-Ishlinskii (Visintin §III.4) |
| Inverse Egger 2025 Schur | Regularization in Visintin sense + monotone-operator theory |
| Motor FE convergence guarantee | Hilpert inequality (Visintin Thm III.2.6) |
| Lamination Henrotte 2015 | Visintin §XI.7 (Preisach homogenization) |
| Vector E&S / Vector Energy | Visintin Thm 2.11 (vector play) |

## References

- **A. Visintin**, "Differential Models of Hysteresis", Applied
   Mathematical Sciences vol 111, Springer-Verlag, 1994
   (420 pages, file `Applied Mathematical Sciences.pdf` in
   `W:/03_文献・論文/00_電磁界解析/30_磁気特性/00_教科書/`).
   Comprehensive treatment of play/stop/Preisach as variational
   inequalities + PDE coupling theory.
- **M. A. Krasnosel'skii & A. V. Pokrovskii**, "Systems with
   Hysteresis", Springer-Verlag, 1989 (English translation of 1983
   Russian original; 429 pages, file `System with Hysteresis.pdf` in
   the same folder).  Foundational work — defines play/stop as
   mechanical transducers, proves identification theorem,
   Ishlinskii's material.  Sections 2, 3, 35, 39 are the canonical
   references cited in every engineering paper.
- **S. Bobbio, G. Miano, C. Serpico, C. Visone**, "Models of
   magnetic hysteresis based on play and stop hysterons", IEEE TMAG
   33:4417-4426, 1997 — bridges Krasnoselskii-Visintin to the
   electrical-engineering audience.  Cited in `play` topic.
- **M. Brokate**, "Some mathematical properties of the Preisach
   model for hysteresis", IEEE TMAG 25:2922-2924, 1989 — proves
   equivalence of play model (with input-independent shape) to
   scalar static Preisach.

## Reading guide

For Radia / lab work, the bare minimum from each book:

**Visintin** (start with these sections):
- Introduction (10 pages) — overview
- I.1 "What is Hysteresis?" — concept
- III.1, III.2, III.3 — Play and Stop formal definitions
- III.4 — Prandtl-Ishlinskii (engineering name for what
   `MatEnergyHysteresis` uses)
- IV.1, IV.2 — Preisach operator + Preisach plane
- VIII.1-VIII.3 — Hysteresis and semigroups (only if proving
   well-posedness of your specific PDE)

**Krasnoselskii** (start with these sections):
- §1-§3 — Transducer + Generalized Play + Stop
- §39 — Rheological models (mechanical analogy — see also
   `rheological_models` topic)
- §35 — Ishlinskii's material (limit K → ∞)

Total essential reading: ~60 pages from Visintin + ~40 pages from
Krasnoselskii = readable in 1-2 sessions.  The rest of both books
is reference material for PDE proofs.

Cross-ref: `play`, `lab_core`, `preisach`, `rheological_models`,
`input_dependent_shape`, `lamination_homogenization`.
"""


RHEOLOGICAL_MODELS = r"""
# Rheological models of hysteresis
# (Krasnoselskii §39 + Visintin §III.4 — the mechanical analogy
#  that motivates Bergqvist-Henrotte energy-based hysteresis)

## Why rheological models matter

The Bergqvist 1997 / Henrotte 2006 energy-based hysteresis model
splits the magnetic field h into reversible h_r (spring) and
irreversible h_i (friction).  This is **NOT** a pure mathematical
formalism — it is a direct analog of a **mechanical rheological
model** with elastic spring elements (E) and plastic friction
elements (P).

If you grasp the rheological picture, the formulas in Henrotte
2006 + Henrotte 2015 + Jacques 2018 become natural rather than
ad-hoc.

## The elementary elements (Krasnoselskii §39.1)

**Elastic element E_i** (spring with modulus k_i):
```
u_i(t) = k_i · x_i(t)                                (39.2)
```
where u_i = stress on element, x_i = deformation.

**Plastic element E_j** (rigid-perfectly-plastic with yield T_j):
```
ẋ_j(t) ≥ 0   if u_j(t) = T_j       (yielding upward)
ẋ_j(t) = 0   if |u_j(t)| < T_j      (frozen, inside yield)
ẋ_j(t) ≤ 0   if u_j(t) = -T_j      (yielding downward)         (39.1)
```

The plastic element acts like dry friction: it locks until the
applied force exceeds the yield threshold T_j, then slips.

## Graph representation (Krasnoselskii §39.2)

The whole rheological model is encoded as an **oriented graph G**:
- m vertices (bars S_i)
- (n+1) edges (one element E_j per edge, plus a designated
   "input-output" edge e_0)
- n_c cycles
- Incidence matrix {α_iq} and cycle matrix {β_jq}

The graph structure determines all the algebra.  No physical
intuition is needed once the graph is drawn — the **incidence
matrix gives the equations of motion automatically**.

## The two transducers M and W (Krasnoselskii §39.3-39.5)

The same rheological model defines TWO dual transducers depending
on which variable is taken as the input:

| Transducer | Input | Output | Variational structure |
|------------|-------|--------|------------------------|
| **M** (force-driven) | force v(t) on input edge | deformation w(t) | multidim PLAY (Thm 39.2) |
| **W** (deformation-driven) | deformation w(t) | reaction force v(t) | multidim STOP (dual) |

**Theorem 39.2** (Krasnoselskii):
```
W[t_0, x_0] v(t) = b + L(Z) v(t)
```
where L(Z) is a multi-dimensional play with polyhedral
characteristic Z = K (the prism of admissible (v, x)).

This means: **any rheological model with n elastic + p plastic
elements is equivalent to a multidim play in R^p with appropriate
polyhedral admissible set**.  Reverse: any multidim play can be
realized as some rheological model.

This is the **equivalence theorem** that gives meaning to the
"K play operators in parallel" decomposition used in
`MatPlayHysteresis` and `MatEnergyHysteresis`.

## Bergqvist-Henrotte mechanical analogy (the EM translation)

Bergqvist 1997 maps the rheological elements onto magnetism:

| Rheological | Magnetism | Variable |
|-------------|-----------|----------|
| Spring (E_i) with stiffness k_i | Reversible domain wall bending | h_r,i |
| Plastic friction (P_j) with yield T_j | Bloch wall pinning (impurities) | h_i,j with κ_j = T_j |
| Series spring + friction (E-P cell) | Single hysteron with anhysteretic + pinning | h^k_r + h^k_i |
| Parallel of n+1 E-P cells | Sum of n+1 fractions with weights ω^k | Σ_k ω^k h^k |
| Total deformation | Magnetic polarization J | J = Σ_k J^k |
| Total stress | Applied field h | h = h_r + h_i |

Henrotte 2006 mechanical analogy (the "grey circle" picture):
- h_r is the spring force (1-1 with J via anhysteretic curve)
- h_i is the friction force, **bounded** |h_i| ≤ κ (subgradient
   of an indicator function of a sphere of radius κ)
- Inside |h_i| < κ: J is locked (frozen)
- At |h_i| = κ: J slips parallel to (h - h_r) until equilibrium

This is **isomorphic to Krasnoselskii §39 transducer W** with one
spring + one friction element + cycle of one (n_c = 1).

## Why this matters for FE convergence

The rheological viewpoint clarifies WHY energy-based hysteresis has
better FE properties than pure-mathematical Preisach:

1. **Convexity**: each spring's elastic energy U_k(J_k) is convex.
   Sum of convex = convex.  Hence the variational problem at each
   time step is a **convex minimization** → unique minimizer →
   Newton converges globally.
2. **Monotonicity**: the constitutive map h ↦ J is **maximal
   monotone** (Krasnoselskii §39.4 Theorem 39.3).  This is the
   right structure for FE + Newton.
3. **Vector extension automatic**: the rheological model has no
   preferred direction — it works equally well in R^1, R^2, R^3
   (Krasnoselskii §16-§22).  Hence Bergqvist-Henrotte vector
   hysteresis "for free".

## Identification = graph parameter estimation

Calibrating Bergqvist-Henrotte from measured B-H major loop:

**Krasnoselskii §39.3** gives the **inverse problem**:
given the input-output curve of transducer M (or W), recover the
graph topology G + spring constants {k_i} + yield thresholds {T_j}.

The answer is non-unique in general (different graphs can give the
same input-output curve — congruency classes).  The standard
practical choice is the **n+1 fraction model** (one spring k_0
representing the elastic component, plus n series E_k - P_k cells
representing pinning at thresholds κ_1 < κ_2 < ... < κ_n).

This is exactly what `MatEnergyHysteresis` does with the .hys
file's K = n+1 operators + log-spaced thresholds eta_k = κ_k.

## Reading Krasnoselskii §39

Sections that matter:
- §39.1 (p384-386): Elementary elements + graph definition
- §39.3 (p389-390): Transducer M = force-driven, output deformation
- §39.4 (p390-393): Properties of M, multidim play equivalence
- §39.5 (p393-): Transducer W = deformation-driven (dual)

Cross-refs:
- Visintin §III.4 (Prandtl-Ishlinskii models) — same algebra but in
   variational-inequality language
- Bobbio 1997 IEEE TMAG 33:4417 — first translation of this to
   magnetic hysteresis
- Henrotte 2006 IEEE TMAG 43:899 — Bergqvist + vector extension
- Jacques 2018 ULiège PhD monograph (254 pages) — full development
   for FE motor analysis
- `lab_core` topic — the lab's production implementation
- `energy_based` topic — Bergqvist/Henrotte/Francois-Lavet/Egger lineage
"""


LAB_CORE_PRODUCTION = r"""
# ★★ Lab CORE method — production deep dive
# (Henrotte 2006 COMPEL + Francois-Lavet 2013 J Comp Appl Math +
#  Matsuo 2003 IEEE TMAG, with concrete formulas, parameters, FE algorithm)
#
# Companion to `lab_core` (the high-level overview).  Use THIS topic
# when actually implementing or calibrating the B-input Stop+Energy
# model in a Radia / NGSolve FE solve.

## 1. The four-parameter elementary model (Henrotte 2006)

Inputs: applied magnetic field `h(t)` (vector).
Internal: reversible "spring" field `h_r`, polarisation `M`.
State: stored at each FE element/Gauss point.

**Decomposition** (Bergqvist 1997, fundamental relation):

```
h = h_r + h_i              (vector sum; h_r, h_i, M all 2D/3D)
h_r = ∂_M ρ^Ψ              (gradient of stored magnetic energy density)
h_i ∈ G                    (irreversible "friction force", subgradient)

G = {h_i : |h_i| ≤ κ if Ṁ = 0,
            h_i = κ e_Ṁ   if Ṁ ≠ 0}    (the "grey circle" of Fig. 2)
```

**Equilibrium** (eq 10 of Henrotte 2006):

```
if |h - h_r| < κ:  h_i = h - h_r,  Ṁ = 0       (frozen, inside circle)
if |h - h_r| = κ:  h_i = κ e_Ṁ,    Ṁ ≠ 0       (slipping on boundary)
```

The non-univocity of h_i when Ṁ = 0 is the **memory effect**.

**Saturation curves**: four candidates (Henrotte 2006 Fig 1),
all scaled so L'(0) = L(∞) = 1:

| Function | Notes |
|----------|-------|
| Langevin `coth(x) - 1/x` | Boltzmann-derived; not invertible analytically |
| `tanh(x)` | Simple, analytically invertible |
| `atan(πx/2)/(π/2)` | Heavier-tailed |
| `(-1 + √(1+4x²))/(2x)` | Bounded, smooth |
| **`atanh(x/J_S)`** | **Francois-Lavet 2013 choice** — analytically integrable AND invertible |

**Constitutive law**:

```
M(h_r) = M_s · L(h_r / h_o) · ê_{h_r}        (vector, parallel to h_r)
b      = M(h_r) + μ_0 (1 + χ_e) h            (induction)
```

Four parameters: `M_s, h_o, χ_e, κ`.

## 2. The N+1 fraction combined model (Henrotte 2006 §V)

Single κ gives rectangular loops (Fig 3) — insufficient for engineering.
Decompose M into n+1 fractions with different friction forces:

```
M = Σ_{k=0}^{n} ω^k · M^k          with  Σ ω^k = 1
M^k subjected to (10) independently with friction κ^k
κ^0 = 0 (reversible Bloch wall bending fraction)
0 < κ^1 < κ^2 < ... < κ^n          (log-spaced typical)
```

Per-fraction equilibrium:

```
h = h_r^k + h_i^k              for k = 0, ..., n
M^k = ω^k M(h_r^k)
```

Total magnetisation = sum: `M = Σ_k ω^k M(h_r^k)`.

**Parameter count**: 2n + 4 (M_s, h_o, χ_e, plus n+1 values of κ^k and
n values of ω^k since Σ ω^k = 1).  Typical lab choice: **n = 5 fractions**
(Henrotte 2006 Fig 6 shows excellent match to electrical steel).

## 3. Variational form (Francois-Lavet 2013) — production version

Reformulate as a **minimization problem** at each time step.
Internal energy density `u(J)` and dissipation `d = χ |J̇|`:

```
u(J)  =  ∫_0^J h_r(x) dx                          (stored)
d     =  χ |J̇| = h_i · J̇                         (friction power)
```

**Pseudo-potential** trick: since d depends on J̇ but we want a J-function,
define D(J, J_p) := χ |J - J_p| (Francois-Lavet eq 12).
Then ∂_J D ≈ χ J̇/|J̇| = h_i for small Δt.

**The functional to minimize at each time step** (eq 15):

```
Ω(h, J, J_p) = u(J) - h·J + χ |J - J_p|
```

**N-cell variant** (eq 23):

```
Ω = Σ_{k=1}^{N} Ω^k(h, J^k, J_p^k)
Ω^k = u^k(J^k) - h·J^k + χ^k |J^k - J_p^k|
```

**Crucial property**: each Ω^k is INDEPENDENT of the others (no coupling).
Minimise N separate convex functionals — embarrassingly parallel.

**Saturation curve choice** (Francois-Lavet eq 16-18):

```
h_r^k(x) = α · atanh(x / J_S^k)                  (Bergqvist scalar form)
u^k(J^k) = α J_S^k [(J^k/J_S^k) atanh(J^k/J_S^k)
            + (1/2) ln|(J^k/J_S^k)² - 1|]       (analytically integrable)
```

α is shared across fractions (empirical result, Francois-Lavet §3.1).
Per-fraction parameters: (J_S^k, χ^k).  Total = **2N + 1 parameters**.

## 4. FE algorithm (Francois-Lavet 2013 Algorithm 1) — production

Magnetic vector potential `a` formulation, b = curl a, weak Ampere:

```
∫_Ω (1/μ_0) curl a · curl a' dΩ
  - Σ_{k=1}^{N} ∫_Ω (1/μ_0) J^k · curl a' dΩ
  = ∫_Ω j_s · a' dΩ        ∀ test functions a'
```

The N internal variables J^k are the result of the IMPLICIT relation
J^k = Update(h, J_p^k) with h = (curl a - Σ J^k) / μ_0.

**Algorithm** (per time step):

```python
# At time t:
while ||Δa|| > tol_picard:                 # Picard outer iteration
    for elem in all_elements:
        h_elem = (curl(a)[elem] - sum_k J_k[elem]) / mu_0
        for k in range(N):
            # Minimize Ω^k for J^k:
            #   Ω^k = u^k(J^k) - h·J^k + χ^k |J^k - J_p^k|
            # First-order descent, init guess J = J_p, < 20 iter
            J_k[elem] = argmin_J( u_k(J) - h_elem·J
                                   + chi_k * np.linalg.norm(J - Jp_k[elem]) )
        assemble_into_linear_system()
    solve_linear_system_for_a()
# After convergence, commit state:
for k in range(N):
    Jp_k = J_k.copy()                       # snapshot for next time step
```

**Storage**: 2N additional vector unknowns per element (J^k current + J_p^k previous).
For 3D, that's 6N scalars per element.

**Convergence**: Picard "always converges in less than 20 iterations" per
Francois-Lavet's tests on a 3-phase transformer T-joint.  Inner J^k
minimization uses J = J_p as initial guess (avoids the angular point at
h_i ≠ 0).

**Why this is FE-stable**:
- Each Ω^k is **convex** (u^k convex from monotone f_k, χ^k|·| convex)
  → global Newton convergence guaranteed
- Cell-based decomposition → **embarrassingly parallel**
- Symmetric linear system (no Jacobian asymmetry)
- Storage independent of K (only J^k per elem, not Preisach staircase)

## 5. Identification: two routes

### Route A — Virgin curve (Francois-Lavet 2013 §3.1, RECOMMENDED)

Inputs: first magnetisation (virgin) curve B(H) measured on Epstein /
Single Sheet Tester.

Method: increment H from 0, identify (J_S^k, χ^k) pairs systematically.
- For H below first threshold (H < κ^1): only k=0 fraction active
  → fit α and ω^0 J_S^0 from initial slope
- Above first knee (κ^1 < H < κ^2): k=0 and k=1 active
  → fit κ^1 and J_S^1
- Continue for higher fractions...

**Concrete example: M250-50A non-oriented electrical steel**
(Francois-Lavet 2013 Table 1, α = 65 [A/m] shared):

| Fraction k | J_S^k [T] | χ^k [A/m] |
|------------|-----------|-----------|
| 0 (Bloch wall reversible) | 0.11 | 0 |
| 1 | 0.3 | 10 |
| 2 | 0.44 | 20 |
| 3 | 0.33 | 40 |
| 4 | 0.04 | 60 |

(Variant: 3-fraction subset for fast tests = (0.11, 0), (0.8, 16), (0.31, 47).)

### Route B — Coercivity staircase (Henrotte 2006 §VI)

Inputs: family of SYMMETRIC B-H loops at increasing amplitude.
Plot coercive field h_coer vs peak field h_peak (Henrotte Fig 5).

The model predicts:

```
h_coer(h_peak) = Σ_{k=0}^{m(h_peak)} ω^k κ^k / Σ_{k=0}^{m(h_peak)} ω^k
```

where m(h_peak) = highest k with κ^k < h_peak.
This is a **staircase function** in h_peak; fit by adjusting (ω^k, κ^k).
Henrotte 2006 Fig 5 shows 5-fraction fit to electrical steel.

**Compare**: Route A simpler (only virgin curve needed) but requires
careful α fit. Route B more data-hungry (multiple symmetric loops) but
more robust to virgin-curve noise.  Lab default: **Route A** for new
materials, **Route B** for re-validation.

## 6. Vector extension — automatic

Since the formulation is vector from the start (h, h_r, h_i, J all
2D / 3D), no extra calibration is needed for vector hysteresis.
The 1D-identified parameters {J_S^k, χ^k, α} carry over directly.

Francois-Lavet 2013 §3.2 demonstrates this on a 3-phase transformer
T-joint: the same 3-fraction set (Table 1a above) gives correct
elliptical b/h trajectories at all 11 sample points around the
T-joint, **identified from 1D Epstein measurements alone**.

This is THE killer feature vs Preisach (vector Preisach needs separate
calibration with RSST measurements).

## 7. Higher harmonics + minor loops — automatic

Same N parameters reproduce minor loops driven by superimposed higher
harmonics (Francois-Lavet 2013 Fig 4).  Going from 3 → 5 fractions
significantly improves minor loop fidelity.  Going from 5 → 8 gives
small improvement (Fig 5).  Engineering default: **5 fractions**.

## 8. Inverse evaluation: Egger 2025 Schur complement (production)

The **inverse operator** H ↦ J(H) is naively O(N³) because the N J^k
share the term (ν_0/2) |B - ΣJ^k|² in the alternative B-input
formulation.  Egger 2025 (arXiv:2507.15289, see `efficient_eval` topic)
solves this with a Schur complement reducing the Newton system to a
d × d block (d = 2 or 3 = dim of J) → **O(N) per evaluation**.

Combined with regularization |x|_ε := √(|x|² + ε), ε = 1e-6, the
non-smoothness at J = J_p disappears and standard damped Newton
converges globally + locally quadratically.

This is what `MatEnergyHysteresis` in Radia uses internally.

## 9. Connection to mathematical foundations

The Ω^k functional minimization IS the variational inequality
formulation of a play / stop hysteron (Visintin §III.2-III.3,
Krasnoselskii §2-§3).  Each cell is a "spring + slider" mechanical
analog (Krasnoselskii §39, Henrotte 2006 Fig 2).  The N-cell sum is
Ishlinskii's material (Krasnoselskii §35) / Prandtl-Ishlinskii model
(Visintin §III.4) — see `mathematical_foundations` and
`rheological_models` topics.

## 10. Implementation in Radia today

`rad.MatEnergyHysteresis(K, eta, f_k_tables, eps=1e-6)`:
- `K = N` = number of fractions (10-30 typical)
- `eta[k] = κ^k` = pinning thresholds in Tesla (log-spaced)
- `f_k_tables[k] = (r_array, f_array)` = shape functions f_k(r)
  derived from the (J_S^k, atanh) parametric form
- `eps` = Egger 2025 regularization

State management:
- `rad.MatHysSaveState(mat)` → ndarray of shape (K × 3, ) for vector
- `rad.MatHysRestoreState(mat, state)` for time-step restart
- `rad.MatHysCommitState(mat)` after Picard converges (snapshot J_p^k)

For .hys file format see `play` topic.  For radia hysteresis_io.py
see `radia_status('overview')`.

## References (chronological lineage)

| Year | Authors | Contribution | Lab library |
|------|---------|--------------|-------------|
| 1997 | Bergqvist | Dry-friction-like pinning, Physica B 233:342 | 05_Energy_Based |
| 2003 | Matsuo, Shimode, Terada, Shimasaki | Stop+Play combination for silicon steel, IEEE TMAG 39:1361 | 04_Stop_model |
| 2006 | **Henrotte, Nicolet, Hameyer** | h = h_r + h_i vector hysteresis, COMPEL 25:71 | 05_Energy_Based |
| 2013 | **Francois-Lavet, Henrotte, Stainier, Noels, Geuzaine** | **Variational fully-implicit, Ω = u - h·J + χ|J-J_p|, M250-50A params** | 05_Energy_Based |
| 2015 | Henrotte, Steentjes, Hameyer, Geuzaine | 2-step lamination homogenization, IET SMT 9:152 | 81_Thermal_extension |
| 2018 | **Jacques** | ULiège PhD monograph (254p) — canonical reference | 05_Energy_Based |
| 2025 | Egger, Engertsberger, Schafelner | Regularization + Schur → O(N) inverse, arXiv:2507.15289 | 05_Energy_Based |

## Cross-references

- `lab_core` — high-level overview (this topic = production deep dive)
- `efficient_eval` — Egger 2025 inverse operator algorithm
- `energy_based` — Bergqvist/Henrotte/Jacques lineage in catalog form
- `mathematical_foundations` — Visintin variational inequality view
- `rheological_models` — Krasnoselskii §39 spring+friction roots
- `stop_vs_play_silicon_steel` — when neither stop nor play alone suffices
- `input_dependent_shape` — Matsuo trilogy fix for silicon steel
- `lamination_homogenization` — Henrotte 2015 macroscale homogenization
"""


STOP_VS_PLAY_SILICON_STEEL = r"""
# Stop vs Play for silicon steel: when neither alone suffices
# (Matsuo, Shimode, Terada, Shimasaki, IEEE TMAG 39:1361-1364, 2003)
# 04_Stop_model/Application_of_stop_and_play_models_to_the_
#   representation_of_magnetic_characteristics_of_silicon_steel_sheet.pdf

## Setup

Both stop and play models can give H from B (B-input → H-output):

```
Stop: H(B) = Σ_{k=1}^M  g_k(s_k(B))          (eq 1)
        s_k(B) = max(min(B - B^0 + s_k^0, η_k), -η_k)        (stop hysteron op)
Play: H(B) = Σ_{k=1}^M  f_k(p_k(B))          (eq 4)
        p_k(B) = max(min(p_k^0, B + ζ_k), B - ζ_k)           (play hysteron op)
```

Both are equivalent to Preisach for SCALAR static hysteresis (Bobbio 1997).
But for actual silicon steel sheet (JIS 50A290 non-oriented, B_s = 1.46 T),
direct comparison reveals which works.

## The four test cases (Matsuo 2003 Figs 2-5)

20 measured loops on a single sheet tester (1 Hz):
- **Symmetric loops** (Fig 1a): amplitudes k B_s / N, k=1..20, centered at B=0
- **Asymmetric loops** (Fig 1b): first-order reversal curves from negative
   saturation — these probe the dc-bias dependence

### Stop model identified from symmetric loops (Fig 2)

- Symmetric: reproduced reasonably (small artifacts from piecewise-linear
   shape function in saturation knee)
- **Asymmetric: BADLY WRONG** — width on B-axis is much larger than
   measured

### Stop model identified from asymmetric loops + symmetrization (Fig 3)

- Symmetric: also wrong (the symmetrization correction does not
   compensate enough)
- Asymmetric: better than Fig 2 but still inaccurate

### Play model identified from symmetric loops (Fig 4)

- Symmetric: GOOD (better than stop because uses 2M hysterons vs M for
   stop — same number of measurement points yields more hysterons)
- Asymmetric: small discrepancy but USABLE

### Play model identified from asymmetric loops + symmetrization (Fig 5)

- Symmetric: asymmetric representation (symmetrization improves it but
   not perfectly)
- Asymmetric: best of all stop/play variants alone

## Root cause: silicon steel violates "equal vertical chords"

The property required for representation by stop with input-independent
shape function:

```
Δh(a, B_0) = h^+(a, u + B_0) - h^-(a, u + B_0)
            depends ONLY on amplitude a, NOT on dc bias B_0
```

Silicon steel asymmetric loops (Fig 1b) have **chord widths that DEPEND
ON dc bias** — at B near saturation the chord is wide, at B near zero
it's narrower.  Therefore pure stop model cannot represent silicon steel
correctly under arbitrary loading.

(Play with symmetric-loop identification works better because it uses
the **congruency property** instead of equal vertical chords; this is
a weaker constraint that silicon steel satisfies more closely.)

## The combination model C_i (Matsuo 2003 §IV) — the lab solution

Iteratively refine by alternating Play and Stop corrections:

```
C_0 = P(h)                                    (play, symmetric-id)
C_1 = C_0 + S(h - C_0)                        (stop residual correction)
C_{2i}   = C_{2i-1} + P(h - C_{2i-1})         (play correction at step 2i)
C_{2i+1} = C_{2i}   + S(h - C_{2i})           (stop correction at step 2i+1)
```

where:
- `P(f(B))` = play model identified from symmetric loops of f(B)
- `S(f(B))` = stop model identified from asymmetric loops of f(B) +
   symmetrization

The corrections converge: C_9 and C_10 give visually-perfect match to
all measured loops (Figs 8, 9).  Lab production: use **(C_9 + C_10)/2**
for symmetric/asymmetric averaging (Matsuo 2003 Fig 10).

In practice C_3 or C_5 is sufficient for engineering applications.

## When this matters in the lab

If you are simulating:
- **Pure AC (no dc bias)** → scalar Play is enough; ignore combination
- **DC + AC superposition** (PWM-driven inductor, DC-biased choke,
   permanent-magnet biased core) → combination model needed OR use
   `input_dependent_shape` (Matsuo 2004/2005 alternative)
- **Vector hysteresis** (motor stator rotational losses) → use Energy-Based
   `MatEnergyHysteresis` instead — naturally vector + handles minor loops

## Implementation note

The combination model is currently **NOT exposed as a public Radia API**.
For dc-biased silicon steel, options in priority order:
1. Use `input_dependent_shape` with weighting w(B) (Matsuo 2004/2005) —
   the only one of these that PROVABLY removes the equal-vertical-chords
   limitation (it is the necessary fix per the Matsuo-Shimasaki 2005
   representation theorem).
2. Implement Matsuo combination C_5 manually by alternating
   `MatPlayHysteresis` and a future `MatStopHysteresis`.
3. `MatEnergyHysteresis` ONLY if it uses a Henrotte-style nonlinear
   reversible spring with the friction bound in H (kappa in A/m). That
   structure CAN break equal vertical chords. A B-input stop whose
   states are clipped to a fixed box in B with an input-INDEPENDENT
   shape function does NOT -- see CAUTION below.

**CAUTION (corrected 2026-05-29).** Do NOT assume "energy-based => dc
bias handled." The energy/variational (Bergqvist) wrapper does NOT by
itself remove the equal-vertical-chords limitation. A B-input stop with
states clipped to a fixed box [-eta_k, eta_k] in B has the
equal-vertical-chords property for ANY shape function g_k (linear or
nonlinear) AND under the gamma-regularised variational update, because
the clip box is bias-independent (Matsuo-Shimasaki 2005 theorem). This
was verified numerically on the IGTE'26 Bergqvist energy-stop model:
the minor-loop vertical chord Delta_h(a, B0) is bias-independent to
0.0-0.1% over B0 in [0, 1.2] T
(W:/.../2026_09_IGTE_Symposium/.../python_stop/demo_equal_vertical_chords.py).
Equal vertical chords is removed ONLY by (a) an input-DEPENDENT shape
w(B)*g_0 -- which makes the per-hysteron coenergy B-dependent,
U_k(s_k,B)=w(B) integral g_0, breaking the *simple* fixed-domain
potential. RESOLVED (route i, 2026-05-29): a generalised B-dependent
energy IS thermodynamically consistent if H is taken as the free-energy
gradient H=d/dB[w(B)U0]=w*H0+w'*U0 -> D=w*sum_clipped g0(+-eta)*Bdot>=0
(the w' terms cancel; the bare Matsuo H=w*H0 drops w'U0 and is NOT a
gradient). See the `bergqvist_binput_stop` route-(i) section and docs
HYSTERESIS_THEORY.md §7.5.9. Or (b) a Henrotte
H-friction model with a nonlinear reversible spring (kappa in H, NOT a
B-clip). The earlier "spring+friction naturally handles asymmetric
loading without the equal-vertical-chord limitation" applied to (b) and
was WRONG when read as covering the B-clip energy stop.

## Cross-references

- `play` — Play model basics
- `lab_core_production` — Energy-based recipe (recommended path #1)
- `input_dependent_shape` — Matsuo trilogy w(B) extension (path #3)
- `preisach` — Preisach formalism (mathematical equivalent of pure stop/play)
- `mathematical_foundations` — Hilpert inequality for FE+hysteresis
   well-posedness

## Reference

T. Matsuo, D. Shimode, Y. Terada, M. Shimasaki, "Application of stop
and play models to the representation of magnetic characteristics of
silicon steel sheet", IEEE Trans. Magn. 39(3):1361-1364, May 2003.
DOI: 10.1109/TMAG.2003.810171.  File:
`W:/03_文献・論文/00_電磁界解析/30_磁気特性/02_ヒステリシスモデル/
 04_Stop_model/Application_of_stop_and_play_models_to_the_
 representation_of_magnetic_characteristics_of_silicon_steel_sheet.pdf`.
"""


FORWARD_INVERSE_VARIATIONAL = r"""
# Forward and Inverse energy-based vector hysteresis operators
# (Egger / Engertsberger / Domenig / Roppert / Kaltenbacher,
#  IEEE Trans. Magn. 61(4):7300207, 2025 — open access)

This paper provides the **rigorous mathematical foundation** for the
energy-based hysteresis model used by `MatEnergyHysteresis` in Radia.
It formalizes BOTH the forward operator (used in scalar-potential FE)
AND the inverse operator (used in vector-potential FE) as
**convex minimization problems**, and proves they describe the SAME
material law via convex duality.

## The four Assertions (theorems)

### Assertion 1 — Forward operator via co-energy density

The Bergqvist-Henrotte material law `J = J(H; J_p)` defined as

```
J(H; J_p) = arg min_J [ U(J) - ⟨H, J⟩ + χ |J - J_p| ]    (3)
B = μ_0 H + J(H; J_p)                                      (4)
```

is equivalent to

```
B = ∂_H w*(H; J_p)                                         (5)
```

where the **co-energy density** is

```
w*(H; J_p) = (μ_0/2)|H|² - min_J [U(J) - ⟨H,J⟩ + χ|J-J_p|]  (6)
```

**Proof sketch**: Apply Danskin's theorem to F(H) = min_J f(H,J)
where f(H,J) = U(J) + ⟨H,J⟩ + χ|J - J_p|.  Then
∂_H F(H) = ∂_H f(H, J(H;J_p)) = -J(H;J_p).  Substitute into (6) +
elementary computation.

**Use case**: scalar potential FE formulation
`H = H_s - ∇ψ` needs `B = B(H)` constitutive law → use (5).

### Assertion 2 — Inverse operator via energy density (single pinning)

By convex duality, the inverse `H = H(B; J_p)` has the compact form

```
H = ∂_B w(B; J_p)                                          (8)
w(B; J_p) = min_J [ (1/(2μ_0))|B - J|² + U(J) + χ|J - J_p| ]  (10)
```

**Proof**: Apply Danskin's theorem to f(B,J) = (1/2μ_0)|B-J|² +
U(J) + χ|J-J_p|.  Then ∂_B w = (1/μ_0)(B - J̃(B;J_p)) where J̃
is the minimizer.  Using B = μ_0 H + J, comparing with (4) shows
J̃(B;J_p) = J(H;J_p) — same internal state.

**Use case**: vector potential FE formulation `B = curl a`
needs `H = H(B)` constitutive law → use (8).

### Assertion 3 — Forward operator with N pinning forces

Decompose J = Σ_k J_k and J_p = Σ_k J_p,k.  Co-energy density:

```
w*(H; {J_p}) = (μ_0/2)|H|² - Σ_k min_{J_k} [U_k(J_k) - ⟨H,J_k⟩
                                              + χ_k |J_k - J_p,k|]   (13)
```

The N inner minimizations are **INDEPENDENT** — embarrassingly parallel.

Material law:

```
B = μ_0 H + Σ_k J_k(H; J_p,k)                              (14)
```

### Assertion 4 — Inverse operator with N pinning forces

Energy density (the inner minimization is now COUPLED across k):

```
w(B; {J_p}) = min_{J_k} [ (1/(2μ_0))|B - Σ_k J_k|²
                          + Σ_k (U_k(J_k) + χ_k |J_k - J_p,k|) ]    (15)

H = (1/μ_0) (B - Σ_k J̃_k(B; {J_p}))                       (16)
```

The coupling through `(1/2μ_0)|B - ΣJ_k|²` means inverse is more
expensive than forward (Table I below).

## Computational complexity (Table I, Egger 2025)

500 time steps, AMD Ryzen 5 PRO 5675U @ 4.3 GHz:

| N_χ | Forward (ms) | iter | Inverse (ms) | iter | Ratio |
|-----|--------------|------|--------------|------|-------|
| 2   | 0.253        | 6.57 | 0.724        | 6.05 | 2.99 |
| 5   | 0.386        | 7.99 | 1.74         | 7.56 | 4.66 |
| 10  | 0.420        | 8.29 | 4.10         | 8.31 | 9.84 |
| 15  | 0.483        | 8.33 | 7.94         | 8.80 | 16.3 |
| 20  | 0.523        | 8.42 | 12.4         | 8.88 | 25.1 |

→ **Forward**: O(N_χ).  **Inverse**: closer to O(N_χ²) or O(N_χ³).
→ Newton iteration count ~8 INDEPENDENT of N_χ (good).
→ Egger 2025 arXiv:2507.15289 reduces inverse to O(N_χ) via **Schur
   complement** trick — see `efficient_eval` topic.

## Standard internal energy density (Prigozhin 2016, used by Egger 2025)

```
U_k(J) = -(A · w_k · J_s / π) · log(cos( (π/2) · |J|_ε / (w_k J_s) ))   (17)
```

where:
- `J_s` = saturation polarization (e.g. 1.54 T for Prigozhin test,
   1.573 T for Egger TEAM 32)
- `A` = steepness of anhysteretic curve (e.g. 38 single-pinning,
   50 multi-pinning, 90.302 for TEAM 32)
- `w_k` = weight factor for pinning force χ_k (Σ w_k = 1)
- `|J|_ε := (|J|² + ε)^(1/2)` regularized norm with ε = 1e-8

## Regularization for non-smooth |J - J_p|

The functional Ω has an **angular point at J = J_p** (the |·| norm
is not differentiable there).  Standard Newton fails.  Two fixes:

1. **Regularization**: replace `|J - J_p|` by `|J - J_p|_ε` with
   ε ≈ 1e-8.  Adds smooth approximation; Newton with line search
   converges globally (Egger 2025 §V.A).
2. **Specialized solver**: Prigozhin 2016 / Domenig 2024 use
   subgradient methods on the non-smooth form directly.

Egger 2025 recommends route #1 for unified analysis + good numerics.

## TEAM 32 benchmark validation

Egger 2025 §VI validates BOTH formulations on the
**TEAM problem 32** (Bottauscio 2018):
- 2D iron yoke + 2 coils + air domain
- Driving current alternated periodically (sinusoidal + harmonics)
- Test points C_6, C_{1,2}, C_{3,4} for B_y(t) comparison
- 5 pinning forces, A = 90.302, J_s = 1.573 T

Result (Fig 5-6): **Scalar potential (forward) ≡ Vector potential
(inverse) within discretization error**.  Hysteresis loops at C_6
overlap perfectly.

## FE formulations

**Scalar potential** (H = H_s - ∇ψ, source field H_s = H_0·I_s(t)):

```
min_ψ ∫_Ω w*(H_s - ∇ψ; {J_p^i}) dΩ                          (23)
```

H^1-conforming elements for ψ.  Globally convergent generalized
Newton (Egger 2024 arXiv:2409.01015).

**Vector potential** (B = curl a, 2D scalar a out-of-plane):

```
min_a ∫_Ω [w(curl a; {J_p^i}) - j_s · a] dΩ                 (26)
```

H^1-conforming elements for a.  Same generalized Newton works.

## Implications for Radia

Radia's `MatEnergyHysteresis(K, eta, f_k_tables, eps=1e-6)` uses:
- The forward operator (Assertion 1 + 3) for scalar-potential
- The inverse operator (Assertion 2 + 4) for vector-potential (A-form)
- The Schur complement trick (`efficient_eval` topic) for O(N) inverse
- ε = 1e-6 regularization (matching Egger 2025 §V.A)

When choosing which formulation to use in your FE solve:

| FE formulation | Constitutive needed | Cost | Use case |
|----------------|---------------------|------|----------|
| Scalar potential (H = H_s - ∇ψ) | B(H; J_p) forward | Fast (O(N_χ)) | Magnetostatic, simple sources |
| Vector potential (B = curl A) | H(B; J_p) inverse | Slower (Egger Schur: O(N_χ)) | Eddy currents, source field from coils |

Both give **identical results** at convergence (Assertion 1 + 2 +
TEAM 32 validation).

## Theoretical properties (Egger 2025 §VII discussion)

Both operators are:
- **Strongly monotone** (uniqueness of solutions)
- **Lipschitz continuous** (continuous dependence on data)
- **NOT classically differentiable** (require generalized Newton)

→ Standard Newton fails; **generalized Newton** (Egger 2024
arXiv:2409.01015) converges globally for both.

## References

- H. Egger, F. Engertsberger, L. Domenig, K. Roppert, M. Kaltenbacher,
   "On Forward and Inverse Energy-Based Magnetic Vector Hysteresis
   Operators", IEEE Trans. Magn. 61(4):7300207, April 2025
   (open access, CC-BY 4.0).
- V. François-Lavet, F. Henrotte, L. Stainier, L. Noels, C. Geuzaine,
   "An energy-based variational model of ferromagnetic hysteresis
   for finite element computations", J. Comp. Appl. Math. 246:243-250,
   2013 — the foundational paper Egger 2025 extends.
- L. Prigozhin, V. Sokolovsky, J. W. Barrett, S. E. Zirka,
   "On the energy-based variational model for vector magnetic
   hysteresis", IEEE Trans. Magn. 52(12):1-11, 2016 — the internal
   energy density (17) is from this paper.
- L. D. Domenig, K. Roppert, M. Kaltenbacher, "Incorporation of a
   3-D energy-based vector hysteresis model into the finite element
   method using a reduced scalar potential formulation", IEEE Trans.
   Magn. 60(6):1-8, 2024 — 3D scalar potential implementation.
- H. Egger, F. Engertsberger, L. Domenig, K. Roppert, M. Kaltenbacher,
   "On nonlinear magnetic field solvers using local quasi-Newton
   updates", arXiv:2409.01015, 2024 — generalized Newton theory.

Cross-ref: `efficient_eval` (Schur trick for O(N) inverse),
`lab_core_production` (Bergqvist-Henrotte production formulas),
`energy_based` (lineage overview).
"""


TP_EEC_STEADY_STATE = r"""
# Time-Periodic Explicit Error Correction (TP-EEC) for Play model
# (Kitao / Hashimoto / Takahashi / Fujiwara / Ishihara / Ahagon / Matsuo
#  IEEE Trans. Magn. 48(11):3375-3378, November 2012)

## The problem

For motor / reactor analysis with **dc bias** superposed on ac, normal
time-stepping FE with the play model needs MANY periods to reach
steady state — the dc error component decays slowly because the
play model is **rate-independent** (no eddy-current-like damping of
the dc component).

For a sinusoidal voltage driving a ring core (Fig 1 of Kitao 2012)
with inrush current:
- Normal step-by-step calculation: **dc error remains in steady state**
  (Fig 5a, b)
- The hysteresis loop becomes asymmetric with permanent dc offset

Two cures: (1) gradually ramp the voltage over 3 periods to avoid
the inrush current (Kitao §III.B), or (2) **TP-EEC method** to force
convergence to symmetric steady state regardless of inrush.

## Conventional TP-EEC (Takahashi 2010, 2011)

For the system S(x) + C ∂x/∂t = f, divide one (or half) period into N
time steps.  Correction vector q from linearized auxiliary system:

```
[Σ_{i=1}^{N} S_i + (1 ∓ 1) C] q = -C (x_0 ± x_N)            (8)
```

When the C matrix dominates, the simplified form is:

```
q = (1/2) (-x_0 - x_N),    x_N ← x_N + q                    (9)
```

Apply (9) every N steps.  Drastically accelerates convergence
(Fig 5a — TP-EEC vs normal step-by-step).

## The pitfall with play model — hysteron states not corrected

Naive TP-EEC corrects the UNKNOWN VECTOR x (= flux density,
current, etc.) but **NOT the hysteron states p^ζn**.  So even after
x converges to symmetric, the play hysterons retain dc bias from
the inrush — and the next half-period diverges back to asymmetric.

Result (Fig 6a, b): TP-EEC alone gives **WRONG** dc-biased steady
state when inrush is large.

## Kitao 2012 fix — Novel TP-EEC for Play

Apply the SAME correction principle to the hysteron states:

```
r^ζn = (1/2)(-p_0^ζn - p_N^ζn)        for each hysteron n at each time step
p_N^ζn ← p_N^ζn + r^ζn                                       (10)
```

Procedure:
1. Run a transient analysis starting with x_0
2. Compute (8) or (9) for unknown vector correction q AND (10) for
   hysteron correction r^ζn
3. Update x_N ← x_N + q AND p_N^ζn ← p_N^ζn + r^ζn
4. If not converged, substitute x_N → x_0 and go to step 1

Result (Fig 6): converges to **correct symmetric steady state**,
inrush-induced dc error removed.  Confirmed on JIS 50A470 electrical
steel ring core.

## Ring core circuit equations (Kitao §III.A)

For a ring core with mean radius r_a, n turns, resistance R:

```
∮_C h · ds = 2π r_a · h = n · i                              (5)
v = n · dφ/dt + R · i = n · S · db/dt + R · i                (6)
```

Solve for unknowns (b, i) with backward Euler.

## Simple homogenized eddy current term

For ac analysis with eddy currents (sheet thickness d, conductivity σ):

```
H_ac(B, t) = H_dc(B) + H_E(B, t)                             (3)
H_E(B, t) = k_E · σ d² / 12 · dB/dt                          (4)
```

`k_E` is an empirical anomaly factor tuned so total loss matches
measurement.  Works well at 50 Hz, breaks down at 200 Hz (need
better homogenization — see `lamination_homogenization` topic).

## Identification numerology

Standard Matsuo recipe: **40 symmetric loops** at intervals 0.05 T
from 0.05 T to 2.0 T → identifies play model with M = 40 hysterons
adequate for engineering accuracy on:
- Dust core (no eddy current, k_E = 0): perfect match Fig 3
- JIS 50A470 electrical steel: 50 Hz perfect, 200 Hz needs better
   eddy current model

## Why this matters for Radia

If you're running a long transient with `MatPlayHysteresis` (or
the energy-based variant) and seeing slow convergence to steady
state — implement the Kitao 2012 TP-EEC scheme:

1. Add a snapshot mechanism that grabs (b, i, p_k^state) at start
   and end of period
2. Apply (9) to unknowns + (10) to each hysteron state
3. Iterate periods until ||q|| + ||r|| < tol

For a Radia + NGSolve transient solve:
- `MatHysSaveState(mat)` returns the K×3 state array
- After full period, snapshot state at t = 0 and t = T
- Apply correction (10) to each hysteron
- `MatHysRestoreState(mat, corrected_state)`
- Continue next period

This pattern is NOT yet implemented in Radia; flagged as a future
optimization for long IH transient runs and motor steady-state
analyses.

## References

- J. Kitao, K. Hashimoto, Y. Takahashi, K. Fujiwara, Y. Ishihara,
   A. Ahagon, T. Matsuo, "Magnetic field analysis of ring core taking
   account of hysteretic property using play model", IEEE Trans.
   Magn. 48(11):3375-3378, Nov. 2012.
- Y. Takahashi, T. Tokumasu, A. Kameari, H. Kaimori, M. Fujita,
   T. Iwashita, S. Wakao, "Convergence acceleration of time-periodic
   electromagnetic field analysis by the singularity decomposition-
   explicit error correction method", IEEE Trans. Magn. 46(8):
   2947-2950, Aug. 2010 — original TP-EEC for nonlinear (no hysteresis).
- Y. Takahashi et al., "Convergence acceleration in steady state
   analysis of synchronous machines using time-periodic explicit
   error correction method", IEEE Trans. Magn. 47(5):1422-1425,
   May 2011 — practical TP-EEC on PMSM.
- T. Matsuo, M. Shimasaki, "Time-periodic finite element method for
   hysteretic Eddy-current analysis", IEEE Trans. Magn. 38(2):
   549-552, Mar 2002 — TPFEM with hysteresis foundation.

Cross-ref: `play`, `lab_core_production`,
`lamination_homogenization` (for the eddy current homogenization
needed above 50 Hz).
"""


H_INPUT_B_INPUT_COMBINATION = r"""
# H-input AND B-input Play models combined for ferrite eq circuit
# (Ito / Mifune / Matsuo / Watanabe / Igarashi / Kawano / Iijima /
#  Suzuki / Uehara / Furuya, IEEE Trans. Magn. 49(5):1985-1988, 2013)

## The duality use case

NiZn ferrite for dc-dc converter / reactor:
- Used at high frequency (100 kHz - 2 MHz)
- High resistivity (10^7 Ω·m) → no eddy current → loss = hysteresis
   + residual (domain wall damping)
- Often dc-biased — minor loop accuracy matters
- B is the natural quantity (measured via secondary winding integral)
- H is the natural circuit quantity (proportional to current)

Need BOTH B-input AND H-input play models in the SAME workflow.

## The decomposition

**AC field** (sinusoidal B excitation):

```
H_AC(B) = h_fast(B) + R_1 · dB/dt                            (1)
```

- `h_fast(B)` = rate-independent hysteretic property (static loop)
- `R_1 · dB/dt` = dynamic residual loss (rate-dependent, linear)
- `R_1` determined from coercive-force vs ωB_a slope (Fig 2)

**DC field** (dc-biased measurement):

```
B_DC(H) = b_slow(H) + b_fast(H)                              (2)
```

- `b_fast(H)` = inverse of h_fast(B) → fast component, responds at
   high frequency
- `b_slow(H)` = slow component → does NOT respond above ~1 MHz cutoff

## Equivalent circuit (Fig 3)

Two parallel LCR-type cells:

```
       φ_0 = b_slow(I)         φ_1 = b_fast(I)
        ┌────[L]────[R_0]      ┌────[L]────[R_1]
   I ───┤                      ├──── I
        └────[C_0]──┘          └────[C_1]──┘
```

with I ↔ H and (φ_0 + φ_1) ↔ B.

Circuit equations:

```
I = I_0 + V_0/R_0 + I_C0 = I_1 + V_1/R_1 + I_C1
V_0 = dφ_0/dt,  V_1 = dφ_1/dt
I_C0 = C_0 · dV_0/dt,  I_C1 = C_1 · dV_1/dt
I_0 = h_slow(φ_0),  I_1 = h_fast(φ_1)                         (3)
```

The LEFT LCR represents slow component (cutoff f_0 < operating freq).
The RIGHT LCR represents fast component (domain-wall resonance > 10 MHz).

## The 4-step identification flow (Fig 4)

```
[A1] AC: H_AC measured with B-control  → -R_1·dB/dt
[A2] h_fast(B): B-input play model (handles ac data)
[A3] h_fast(H): inverse via Newton-Raphson with given H amplitudes
[a2] B-input play → inversion to H-input play

[B1] DC: H_DC measured with B-control
[b1] B_DC(H): invert by Newton method (also via h_fast)
[B2] B_DC(H) - b_fast(H) = b_slow(H)
[B3] b_slow(H): H-input play model
[b3] H-input play → inversion to B-input play
[B4] b_slow(B): B-input play model
```

**Key insight**: B-input play (from ac data) and H-input play (from
dc data) are inverted into each other via Newton-Raphson on
F(B) = H - h_fast(B) = 0.

## Newton-Raphson inversion (Appendix)

H-input play:

```
B(H) = Σ_{n=1}^{N_p} f_n(p_n(H))                              (8)
p_n(H) = max(min(p_n^0, H + ζ_n), H - ζ_n)                    (9)
```

Derivative for Newton:

```
dB/dH = Σ_n (df_n/dp_n) · (dp_n/dH)                          (10)
dp_n/dH = 0  if  p_n = p_n^0  (frozen — inside dead band)     (11a)
dp_n/dH = 1  if  p_n ≠ p_n^0  (slipping)                     (11b)
```

B-input play identical structure with (B, ζ, p, f) replaced by
(H, η, s, g).

## Validation (Ito 2013 results)

NiZn ferrite ring core measurements 100 kHz - 800 kHz + dc:

| Frequency | Core loss error | Notes |
|-----------|-----------------|-------|
| DC | 1.7% | h_slow + b_fast |
| 100 kHz | 1.7% | h_fast + R_1·dB/dt |
| 400 kHz | 0.5% | best agreement |
| 2 MHz | 4.6% (modified) | first-order lag element added (Fig 9) |
| DC bias 40 A/m | 5.1% | minor BH loop |
| DC bias 80 A/m | 10.6% | larger discrepancy |

Note (Ito conclusion): "The input-dependent shape function [7] is
required for the play model to achieve more accurate representation
of the minor BH loops" — see `input_dependent_shape` topic.

## Implications for Radia

This paper demonstrates that **both forward (B-input) AND inverse
(H-input) play models** are needed in a SINGLE production workflow:

- B-input play: when input data is B-controlled (most measurements,
   most A-formulation FE)
- H-input play: when input data is H-controlled (some measurements,
   Omega-reduced scalar potential FE)
- **Newton-Raphson inversion** between them takes ~5-10 iterations
   per FE element per time step

For `MatPlayHysteresis` in Radia:
- The forward direction (B → H, using stored f_k shape functions)
   is direct O(K)
- The inverse direction (H → B) requires Newton-Raphson with the
   df_k/dp_k Jacobian — algorithm given in Ito 2013 Appendix
- Both directions are O(K) per evaluation

For ferrite at MHz frequencies, the equivalent circuit decomposition
(slow + fast + residual loss) is the production recipe — implement
as a wrapper around `MatPlayHysteresis` for power electronics
simulations.

## Reference

S. Ito, T. Mifune, T. Matsuo, K. Watanabe, H. Igarashi, K. Kawano,
Y. Iijima, M. Suzuki, Y. Uehara, A. Furuya, "Equivalent Circuit
Modeling of DC and AC Ferrite Magnetic Properties Using H-Input
and B-Input Play Models", IEEE Trans. Magn. 49(5):1985-1988, May 2013.
DOI: 10.1109/TMAG.2013.2243717.

Cross-ref: `play`, `lab_core_production`, `input_dependent_shape`
(for dc-bias minor loop improvement), `mathematical_foundations`
(Krasnoselskii / Visintin formal definitions).
"""


BOBBIO_1997_UNIFICATION = r"""
# Bobbio 1997 — Play-Type and Stop-Type Models unification
# (S. Bobbio, G. Miano, C. Serpico, C. Visone,
#  IEEE Trans. Magn. 33(6):4417-4426, November 1997)

This is **THE foundational paper** cited by every subsequent play /
stop / energy-based hysteresis paper in the engineering literature.

## Key contributions

1. Formalizes **PTM** (Play-Type Model) and **STM** (Stop-Type Model)
   as parallel / series compositions of finite hysterons.
2. Proves **series-PTM ≡ parallel-PTM ≡ Bobbio-Marrucci 1993 model**.
3. Shows PTM is a discretization of the **Prandtl-Ishlinskii model**
   (PIM) and **equivalent to a specific Preisach distribution**.
4. Gives the **identification procedure** via pseudo-inversion of
   piecewise-linear shape function matrix.
5. Demonstrates that **STM is numerically the inverse of PTM**
   (within ~1% error on series connection).
6. Validates on **Vitrovac 7505Z** amorphous alloy toroidal core +
   ferroresonant circuit application.

## PTM and STM definitions

**Parallel PTM** (Play-Type):

```
B(t) = Σ_{k=1}^{N} f_k[h_k(t)]                                (1)
h_k(t) = P_σ_k[H(·)](t)                                       (2)
```

where P_σ = ordinary play hysteron with threshold σ ≥ 0,
σ_k ≥ σ_{k-1}, and f_k are memoryless nonlinear functions.

State vector: x = (X, x_1, ..., x_N)^T = (H, h_1, ..., h_N).
For physical realism: set σ_1 = 0 (first hysteron memoryless).

**Parallel STM** (Stop-Type):

```
H(t) = Σ_{k=1}^{N} g_k[b_k(t)]                                (3)
b_k(t) = S_η_k[B(·)](t)                                       (4)
```

where S_η = ordinary stop hysteron with stop value η, η_k ≤ η_{k-1}.
For physical realism: set η_1 → ∞ (first hysteron memoryless).

**Spins** (loop orientation):
- PTM: positive spin (B-vs-H loop counterclockwise from origin) →
  representing **B = B(H)**
- STM: negative spin (H-vs-B loop clockwise from origin) →
  representing **H = H(B)**

## Equivalence of series and parallel PTM (Bobbio 1997 §II)

**Theorem** (operator composition identity, eq 8):

```
P_ζ1 ∘ P_ζ2[H(·)](t) = P_{ζ1 + ζ2}[H(·)](t)
```

The composition of two play operators with thresholds ζ_1, ζ_2 equals
a single play with threshold ζ_1 + ζ_2.

Applied to series chain `h_0 = H → P_ζ1 → P_ζ2 → ... → P_ζN`:

```
h_k = P_σ_k[H(·)](t)    with    σ_k = Σ_{j=1}^{k} ζ_j         (10)
```

Comparing series eq (6) with parallel eq (1) shows the **series PTM
is equivalent to the parallel PTM** with σ_k = cumulative sum of ζ_j.

This is also equivalent to the **Bobbio-Marrucci 1993** model.

## Equivalence to Prandtl-Ishlinskii model

PIM (continuous superposition of plays) is given by:

```
B(t) = ∫_0^∞ f(ϑ, P_ϑ[H(·)](t)) dϑ                            (A.9)
```

For the choice f(ϑ, z) = Σ_k f_k(z) δ(ϑ - σ_k) (eq 12):

```
B(t) = Σ_{k=1}^{N} f_k(P_σ_k[H](t))                           (13)
```

So **PTM is a discretization of PIM** with delta-function weight.
Similarly STM is a discretization of the H-input PIM (A.10).

## Equivalence to Preisach (Bobbio 1997 §III)

Following [Krasnosel'skii], rotate the Preisach plane (α, β) by 45°:

```
ϑ = (α - β)/√2,   v = (α + β)/√2                              (16)
```

Then the Preisach operator can be rewritten:

```
B(t) = Γ[H(·)](t) = ∫_0^{+∞} f(ϑ, P_ϑ[H(·)](t)) dϑ            (18)
```

with f(ϑ, z) = ∫_{-∞}^z μ̃(ϑ, v) dv - ∫_z^{+∞} μ̃(ϑ, v) dv,
μ̃(ϑ, v) = μ((ϑ+v)/√2, (v-ϑ)/√2)/2.

So **Preisach ≡ PIM of play type** with distribution function
μ̃(ϑ, v) (eq 19).

Conclusion: PTM ⊂ PIM ⊂ Preisach (specific case).

## Identification procedure (Bobbio 1997 §IV)

PTM/STM identification reduces to a **linear least-squares problem**:

For piecewise-linear shape functions u_k, let s = (s_1, ..., s_P) be
the vector of nodal values.  Then the model output is **linear in s**:
Y = A · s where A is an M × P matrix.

Minimization of ||A·s - Q||² (where Q = measured output) reduces to
**pseudo-inversion of A** (M > P, overdetermined).

Tested on Vitrovac 7505Z (B_s ≈ 0.5 T, H_c ≈ 20 A/m, remanence
ratio 0.36):

| Model | Identification residual |
|-------|------------------------|
| PTM (B-input) | 1.5% |
| STM (H-input) | 4.2% |

PTM fits BETTER — Bobbio attributes to its closer match to physical
hysteresis spin direction (counterclockwise for B(H)).

## Numerical inversion: PTM and STM are mutually inverse

**Claim** (Bobbio 1997 §V): Given a PTM, the STM identified on the
SAME material data approximates its inverse operator H = H(B).

Test: connect PTM(H) → STM in series, starting from zero state.
Plot output vs input for several waveforms (Figs 12-14):

| Input | Slope of fitted line | R^2 |
|-------|---------------------|-----|
| i_1 (multi-extrema) | 0.972 | 0.9985 |
| i_2 (smooth oscillation) | 0.994 | 0.9993 |
| Preisach model + STM | 1.000 | very high |

So **STM ≈ PTM^(-1)** numerically (small error from finite
identification residuals).  This is the foundation for:
- Using PTM for B-input FE (A-formulation needs B → H mapping)
- Using STM for H-input FE (Omega scalar potential needs H → B)
- Or composing PTM with STM for complementary FE formulations

## Why the analytical inverse is hard

Bobbio 1997 §V proves: starting from PTM with distribution μ(ϑ,v) = f(ϑ),
the inverse Preisach is well-defined but has a NON-PHYSICAL virgin
curve `B(H) = φ(H)` with sign(φ''(H)) = sign(f(H)) ≥ 0 — constant
convexity, which real magnetic materials do not exhibit.

So the **simple analytical inverse via Preisach distribution
manipulation FAILS** for engineering materials.  The **numerical
PTM ↔ STM identification + Newton-Raphson inversion** (Ito 2013) is
the practical workaround.

## Ferroresonance application

Bobbio §IV closes with a non-trivial nonlinear circuit application:
**RLC ferroresonant circuit** with hysteretic inductor (Fig 9).
Measured + simulated ferroresonance frequency curves (Fig 10) show
the **PTM correctly predicts the non-uniqueness** of steady-state
solutions characteristic of nonlinear LC resonance.

This is the only paper in the lab library that applies play model to
**circuit simulation** (vs FE simulation) — relevant for Cubit
Plugin's future SPICE-export feature.

## Why Bobbio 1997 matters

This paper:
1. Provides the **rigorous foundation** for the lab core method
   (Stop+Energy formulation is Stop-Type Model + energy density).
2. Justifies the **discretization choice** of K = 20-30 hysterons
   (any finite K is a PIM discretization → arbitrarily accurate as K↑).
3. Explains **why PTM ≡ STM^(-1) numerically** (the basis for Ito 2013
   H-input + B-input combination workflow).
4. Establishes **Preisach ⊃ PIM ⊃ PTM/STM** hierarchy — PTM/STM are
   the simplest engineering subset that captures all features.

## References

- S. Bobbio, G. Miano, C. Serpico, C. Visone, "Models of magnetic
   hysteresis based on play and stop hysterons", IEEE Trans. Magn.
   33(6):4417-4426, November 1997.
- S. Bobbio, G. Marrucci, "A possible alternative to Preisach's
   model of static hysteresis", Il Nuovo Cimento 15-D(5):723-734, 1993
   — the originating "Bobbio-Marrucci" model.
- M. A. Krasnosel'skii, A. V. Pokrovskii, "Systems with Hysteresis",
   Springer-Verlag 1989 — where play / stop hysterons were first
   defined (see `mathematical_foundations` topic).
- M. Brokate, "Some mathematical properties of the Preisach model
   for hysteresis", IEEE Trans. Magn. 25(4):2922-2924, 1989 — proves
   Brokate equivalence (Preisach ≡ PIM of play type, the rigorous
   theorem behind §III above).

Cross-ref: `play`, `lab_core_production`, `mathematical_foundations`,
`preisach`, `h_input_b_input_combination` (Ito 2013 uses Bobbio's
numerical inversion claim).
"""


ASP_MODEL = r"""
# ASP model — Asymmetric transition Probability (function-based,
# NO hysterons)
# (Lee, Miyata, Kaido, Matsuo, IEEE Trans. Magn. 50(4):7300104, 2014)

This is the "ASPモデル" reviewed in Matsuo's 2014 consortium
presentation.  It is a **function-based** hysteresis model that
needs NO assembly of hysterons and NO integration — much cheaper
than Preisach or play for the same accuracy on minor loops.

## Core idea

Macroscopic magnetisation = accumulation of microscopic transitions,
described by a **transition probability** p±(H):

```
dB±(H)/dH = p±(H)               (ascending +, descending −)
p−(H) = p+(−H)                  (asymmetric about H = 0)
```

The pair of ASYMMETRIC probabilities p±(H) creates hysteresis (Fig 1).

## BH curve from first-order reversal curves (FORCs)

Every BH curve defined by the last two **branching points**:
initial Pi = (H_i, B_i), terminal Pf = (H_f, B_f).

```
B±(H : Pi, Pf) = B_i + (B_f − B_i) · φ±(H : H_i, B_i, H_f)     (4)
φ±(H) = [g±(H : B_i) − g±(H_i : B_i)]
        / [g±(H_f : B_i) − g±(H_i : B_i)]                      (3)
```

where `g±(H : B_0)` = ascending/descending FORC starting from B = B_0.

**Geometric meaning**: the BH curve is the FORC **compressed along
the B-direction** by ratio

```
c± = (B_f − B_i) / [g±(H_f : B_i) − g±(H_i : B_i)]             (5)
```

## Three rules (memory management)

1. Every BH curve defined by the last 2 branching points (Pi, Pf)
2. Branching points = extrema of H (local max/min create new Pi,
   former Pi becomes Pf)
3. When H passes through H_f, branching points Pi, Pf are DELETED
   and the previous two return (wiping-out / deletion property,
   same as Preisach)

History stored concisely as the **sequence of branching points** —
no hysteron state array needed.

## Normal magnetisation curve (derived, not input)

Via the demagnetising process, the normal curve b(H) satisfies an ODE
(Lee 2014 eq 13):

```
db(H)/dH = 2 b(H) g'_−(−H : b(H)) / [g_−(H : b(H)) − g_−(−H : b(H))]
```

Integrate to get b(H); symmetric loops obtained by setting branching
points ±(H_j, B_j) on b(H).  The match with measured normal curve
(Fig 3) validates the FORC-based construction.

## B-input version (for A-formulation FE)

Using the 1-1 B↔H correspondence (eq 4), redefine with input B:

```
H±(B : H_i, B_i, H_f) = g±^(-1)( B_i0 + (B − B_i)/c± : B_i )    (15)
```

Express FORCs as inverse functions g±^(-1) of B.

## NO congruency property (vs Preisach)

Because p± depends on B_i and B_f (eq 2), the ASP model does NOT
have the congruency property of the classical Preisach model.  This
is actually an ADVANTAGE — silicon steel does not obey congruency
either (cf. `stop_vs_play_silicon_steel`), so ASP fits real
materials better.

## Accuracy vs play model (Lee 2014, JIS 35A300 silicon steel)

Compared against B-input play with 60 hysterons identified from
FORCs (play-F) and from symmetric loops (play-S):

| Test | ASP | play-F | play-S |
|------|-----|--------|--------|
| FORCs (Fig 4) | accurate | disagrees (symmetrization error) | n/a |
| Symmetric loops (Fig 5) | accurate | poor | accurate |
| Asymmetric minor loops (Fig 6,7) | **best** | worse | worse |

**ASP wins on minor loops** because it fully exploits the FORC
information without symmetrization.  Play model needs 60 hysterons
+ symmetrization; ASP needs only the FORCs + a branching-point
sequence.

## When to use ASP

- **Minor-loop-critical applications**: MRI scanner iron cores (the
   paper's motivation — minor loops affect image resolution),
   variable-flux memory motors, DC-biased reactors
- **FORCs available**: ASP is FORC-native; if you only have symmetric
   loops, play / energy-based may be easier
- **Scalar only** (so far): a vector ASP was promised in the 2014
   conclusion but the scalar version is what's published

## ASP vs lab core (Energy-based)

| Feature | ASP | Energy-based (lab core) |
|---------|-----|--------------------------|
| Hysterons | None (function-based) | K cells |
| Calibration | FORCs | virgin curve / coercivity staircase |
| Minor loops | Excellent | Good |
| Vector | Not yet (2014) | Native |
| Thermodynamic consistency | No | Yes |
| FE convergence guarantee | Empirical | Monotone (provable) |

For motor / vector work → Energy-based.  For scalar minor-loop-heavy
(MRI, memory motor) with FORC data → ASP is the more accurate,
cheaper choice.

## Reference

C. Lee, K. Miyata, C. Kaido, T. Matsuo, "Phenomenological Hysteresis
Modeling Based on Asymmetric Transition Probability of
Magnetization", IEEE Trans. Magn. 50(4):7300104, April 2014.
DOI: 10.1109/TMAG.2013.2287030.  Authors: Hitachi Research Lab +
Kitakyushu National College + Kyoto University (Matsuo).

Cross-ref: `play`, `stop_vs_play_silicon_steel` (why congruency-free
helps), `preisach` (the congruent alternative ASP improves on),
`zirka_hdhm` (another FORC-direct non-Preisach model).
"""


ZIRKA_HDHM = r"""
# Zirka congruency-based History-Dependent Hysteresis Model (HDHM)
# (Zirka, Moroz, Marketos, Moses, IEEE Trans. Magn. 40(2):390-399, 2004)

A **non-Preisach** history-dependent model that uses experimental
first-order reversal curves (FORCs) DIRECTLY — no Preisach
distribution function, no differentiation of noisy data.  The
"transplantation method" builds higher-order reversal curves from
internal segments of FORCs.

## Madelung's 5 rules (the refined statement, Zirka 2004 §II)

1. Path of any FORC uniquely determined by its reversal point
   coordinates.
2. **Return-point memory**: curve 3-4-2 originating at point 3
   returns to the initial point 2.
3. **Wiping-out**: if curve 4-3 extends beyond point 3, it follows
   section 3-1 as if loop 3-4-3 never existed.
4. **B-congruency** of reversal curves: all transition curves from
   reversal points with same B-ordinate + sign of dB are congruent
   by parallel translation along H axis.
5. Symmetric reversal from the initial magnetisation curve (rules
   2-3 still hold).

## Why CPM-E fails (the key insight)

CPM-E = "Classical Preisach Model with Experimental FORCs"
(shifting FORC patterns instead of evaluating the Preisach
double integral).  Two variants both FAIL on real materials:
- CPM-E [Mayergoyz]: shift pattern so LEFT endpoint matches →
   right endpoint does NOT match → **discontinuity** (jump)
- CPM-E [Parker 1995]: shift so the matching endpoint coincides →
   **minor loops with crossed branches** / unlimited accommodation

Root cause: real materials do NOT obey the H-congruency (CPM),
B-congruency (Della Torre M model), skew congruency (moving model),
or nonlinear congruency (product model).  So uncritically extending
ANY congruency to higher-order curves produces nonphysical results.

**Most unwanted symptom**: a NEGATIVE SLOPE on the reversal curve →
instability of the magnetodynamic equations.

## The transplantation method (Zirka's fix)

To build a 2nd+ order reversal curve with known overall dimensions
(ΔH, ΔB) from a turning point:

1. Search the WHOLE M-H plane for a **transplant** = an internal
   segment of some FORC with the required (ΔH, ΔB)
2. Two transplants always exist:
   - **α**: found in the column of width ΔH (Fig 6)
   - **β**: found in the strip of height ΔB (Fig 7)
3. Weighted sum gives the actual reversal curve:

```
M(H) = (1 − w) M_α(H) + w M_β(H)            (1)
w = [(M_n − M_rev)/(M_n − M_{n-1})]
    · [1 − (M_n − M_{n-1})/(M_{n-2} − M_{n-1})]      (3, refined)
```

where M_n, M_{n-1} = last two reversal points, M_rev = current
reversal.  w gradually increases as the reversal point recedes
from the minor-loop vertex.

**Negative slopes are impossible by construction** (transplants are
internal FORC segments which always have positive slope) — this is
the model's main advantage over CPM-E.

## Memory organisation

Reversal curve sequence stored as code (Enderby-style):

```
[ R_1  R_3  ... R_{19} ]
[ -R_s R_2  ... R_{20} ]
```

Each R_i contains the (H_i, M_i) turning point coords + (optionally)
the cubic spline coefficients of the reversal curve.  Closed minor
loops are wiped out per rule 3.

- **N_r = 100 reversal curves** sufficient for the most complex
   transients
- **5-10 FORCs** desirable for engineering accuracy (can even be
   hand-drawn once per material)

## HDHM vs HIHM (history-independent)

HIHM = simplest version, trajectory determined only by current
reversal point → leads straight to major loop tip (no return-point
memory).

| Excitation | HDHM vs HIHM loss difference |
|-----------|------------------------------|
| Cosine (smooth) | ~4% (HIHM acceptable) |
| PWM (many minor loops) | ~18% (HIHM fails) |

HIHM is adequate for smooth excitations where minor loops are
comparable to the major loop, but cannot reproduce the small minor
loops in PWM excitation (where the dB/dt sign changes many times
per period).

## Application: PWM transient in conducting sheet (Zirka 2004 §V)

1D Maxwell across half the sheet thickness:

```
∂²H/∂x² = σ ∂B/∂t                                            (5)
```

Neumann BC at surface from winding voltage; fully finite-difference
scheme (eq 8), Gear method for ODE integration (stable down to
fractions of 1 Hz).  Each grid node (sheet layer) has its OWN
magnetisation history → the HDHM is recalled per node per time step.

For a 0.5 mm NO steel under PWM (f=50Hz, carrier 2400Hz, m=0.8):
- Static HDHM loss: 564.3 J/m³
- Dynamic HDHM (with excess loss): 781.1 J/m³
- Measured: 734.2 J/m³

The static model UNDERESTIMATES (need dynamic / excess loss term);
minor-loop areas decrease from sheet surface (node 1) to centre
(absent at node 7) — Fig 14.

## When to use Zirka HDHM

- **Transient simulation with complex history** (PWM, pulse
   transformer, surge) where return-point memory matters
- **FORC data available** (or willing to hand-draw 5-10 curves)
- **Negative-slope instability** observed with CPM-E or Preisach
   inversions
- NOT for vector hysteresis (scalar model)

A free demo + 5 test-material datasets at www.zirka.dp.ua.

## Zirka HDHM vs ASP vs lab core

| Model | FORC-direct | Memory | Negative-slope-free | Vector |
|-------|-------------|--------|---------------------|--------|
| Zirka HDHM | ✓ | return-point | ✓ (by construction) | ✗ |
| Lee ASP | ✓ | branching points | mostly | ✗ (2014) |
| Lab core (Energy) | ✗ (uses virgin curve) | K cell states | ✓ (monotone) | ✓ |

All three avoid the Preisach negative-slope pathology.  For scalar
transients with strong history dependence → Zirka HDHM.  For scalar
minor-loop accuracy → ASP.  For vector / motor → Energy-based.

## Reference

S. E. Zirka, Y. I. Moroz, P. Marketos, A. J. Moses, "Congruency-Based
Hysteresis Models for Transient Simulation", IEEE Trans. Magn.
40(2):390-399, March 2004.  DOI: 10.1109/TMAG.2004.824137.
Dnepropetrovsk National University (Ukraine) + Cardiff University
Wolfson Centre for Magnetics.

Key concept: Madelung 1905 "On the magnetization by fast currents",
Ann. Phys. 17(5):861-890 — the original 5 magnetisation rules.

Cross-ref: `asp_model` (Lee FORC-direct alternative), `preisach`
(the congruent model Zirka improves on), `play`, `lamination_homogenization`
(the conducting-sheet transient application).
"""


CHUA_MODEL = r"""
# Chua-type hysteresis model (dynamic, circuit-oriented)
# (Chua & Bass, IEEE Trans. Circuit Theory CT-19(1):36-48, 1972;
#  building on Chua & Stromsmoe 1970)

A fundamentally DIFFERENT animal from Preisach/play/stop/energy:
the Chua model is a **first-order ODE** designed for **circuit
simulation** (SPICE-style), and it is **rate-dependent (dynamic)**
from the start — the frequency behaviour is built in, not bolted on.

## The model (Chua-Bass 1972, eq 1)

```
dy/dt = w∘(dx/dt) · h∘(y(t)) · g∘(x(t) − f(y(t)))
```

where ∘ denotes functional composition, and:
- `x(t)` = input (e.g. H or current i), `y(t)` = output (e.g. B or flux φ)
- `f, g` : R→R **bijective, nonzero slope** (C¹) — shape functions
- `h(y)` : bounded 0 < α ≤ h(y) ≤ β < ∞ — output-dependent modulation
- `w(dx/dt) ≥ 0` : **the frequency-control weight** (the key innovation)

The `f` curve in the y-x plane divides it: right of f → dy/dt > 0,
left → dy/dt < 0.  The trajectory crosses f only at extrema of y(t).

## The w(dx/dt) weight — controls frequency behaviour (Table I)

| w shape | Frequency behaviour | Physical example |
|---------|---------------------|------------------|
| a: w = dx/dt | Loop WIDENS with freq, collapses at dc | Chua-Stromsmoe 1970 (no dc loop) |
| b: w = \|dx/dt\| slope 1 | Frequency-INDEPENDENT (same loop ω→0 and ω→∞) | when widening not wanted |
| c: w with [−K,K] band | dc loop + widening above threshold K/A | iron-core inductors |
| d: w steeper outside [−K,K] | Loop NARROWS with freq | **fluorescent lamp i-v curves** |
| e: multi-segment | dc loop + moderate-freq widening + HF-insensitive | real iron-core |

This is the unique Chua selling point: by choosing w you reproduce
loop WIDENING (most ferromagnets) OR loop NARROWING (fluorescent
lamps, some dielectrics) OR frequency-insensitivity — all in one
framework.

## Chua-Stromsmoe (C&S) 1970 — the predecessor

Setting w ≡ 1 (or w = dx/dt) recovers the C&S model
(Chua-Stromsmoe, Int. J. Eng. Sci. 9:435, 1971; IEEE TCT-17:564,
1970): `dy/dt = h(y) g(x − f(y))`.  Its FAILING: NO dc hysteresis
loop (loop collapses to a single curve as ω→0), so it can't model
ferrite memory cells under pulse input.  Chua-Bass 1972 fixes this
with the w weight.

## Key mathematical properties (Chua-Bass 1972)

- Prop 1: unique C¹ solution y(t; y_0, t_0) for bounded x ∈ C¹
- Prop 4: for periodic nonconstant x, a UNIQUE periodic y exists
   (independent of initial conditions) — important for steady-state
- Prop 5: y(t) has NO MORE extrema than x(t) per period
- Prop 8: as freq → ∞, peak-to-peak y(t) → 0 (skin-effect analogy)

## Arbitrary excitation (inductor modeling, eq 4)

For inductor / transformer use (port flux φ, current i, voltage v):

```
L di/dt = v − g∘(i − f(φ − Li)) · w[(i − i_2)R]
G dφ/dt = ... (parasitic)
```

Adding parasitic leakage inductance L + large resistance R = 1/G
makes the 2nd-order model 3rd-order.  Drawback: small time
constants → slow at low frequency.  N-winding transformer model
(Fig 10) handles multi-winding; material modeled ONCE, scaled to
core dims (A, l, n) via Table II.

## Applications validated (Chua-Bass 1972)

- 50-50 NiFe orthonol (square B-H loop) — w type c, widening > 3 Hz
- 4-79 Mb-Permalloy — w type c
- BaTiO3 + PbTiO3 + CaTiO3 ferroelectric (asymmetric P-E loop)
- Ferrite memory cell (pulse response, remanent states)

## Why Chua matters / when to use it

| Feature | Chua | Lab core (Stop+Energy) |
|---------|------|------------------------|
| Domain | Circuit (SPICE) | FE field |
| Rate | Dynamic (w weight) | Static (rate-independent) |
| dc loop | Yes (Chua-Bass) | Yes |
| Loop narrowing | Yes (type-d w) | No |
| Vector | No (scalar ODE) | Yes |
| Minor loops | Limited | Good |
| Identification | f, g, h, w nonlinearities | major loop / virgin curve |

**Use Chua** for: SPICE circuit-level transformer/inductor/ferrite
models, ferroresonance circuit analysis, fluorescent-lamp i-v
modeling (loop narrowing), pulse memory cells.  **Do NOT use Chua**
for: FE field analysis, vector hysteresis, accurate minor loops →
use the lab core energy-based method.

## Lab connection — Saito Hosei lineage + ferroresonance

The lab library 06_Chua folder (18 papers, the largest single-model
collection) is dominated by:
- **Saito 兆古 (Hosei University) Chua-type model identification**
   (斎藤流Chua-model同定.pdf, Hybrid model papers) — same group as
   the Cauer Ladder Network (CLN) work in `radia_mcp.mor` and the
   cellular automaton model (`cellular_automaton` topic)
- **Ferroresonance circuit analysis** (Parallel Ferroresonance
   Circuit Analysis by Chua-type Magnetization Model; 鉄共振現象)
   — nonlinear LC resonance + chaos (cf. Bobbio 1997 PTM
   ferroresonance, `bobbio_1997_unification`)
- Fourier-series representation of hysteresis (HBM — harmonic
   balance method)
- A REPRESENTATION OF MAGNETIC HYSTERESIS (Chua foundational)

This makes Chua the lab's go-to for **circuit-coupled / SPICE**
hysteresis (vs FE field, where Stop+Energy reigns).  Relevant to a
future Cubit plugin SPICE-export feature.

## References

- L. O. Chua, S. C. Bass, "A Generalized Hysteresis Model", IEEE
   Trans. Circuit Theory CT-19(1):36-48, January 1972.
- L. O. Chua, K. A. Stromsmoe, "Lumped-circuit models for nonlinear
   inductors exhibiting hysteresis loops", IEEE Trans. Circuit Theory
   CT-17:564-574, Nov. 1970.
- L. O. Chua, K. A. Stromsmoe, "Mathematical model for dynamic
   hysteresis loops", Int. J. Eng. Sci. 9:435-450, 1971.

Cross-ref: `bobbio_1997_unification` (ferroresonance with play
model), `cellular_automaton` (Saito Hosei), `mor_cln` MCP topic
(Cauer Ladder Network, same Saito lineage), `catalog` (model
landscape).
"""


OTHER_MODELS = r"""
# Niche / application-specific hysteresis models
# (Chan SPICE transformer, Bouc-Wen structural, Potter-Schmulian
#  recording media, LLG micromagnetics)

These four models are in the lab catalog but serve **specific
application domains** outside the lab's core soft-iron FE work.
Documented here for completeness + decision-tree depth.

## Chan model (SPICE transformer) — Chan 1991, IEEE TCAD 10:476

The **production SPICE / PSpice transformer-core hysteresis model**
(implemented in DSPICE = Daisy's SPICE; now standard in most SPICE
transformer libraries).

Major loop = **two hyperbolic branches**:

```
Upper (increasing H):  B(H) = B_s · (H + H_c)/(|H + H_c| + H_k) + μ_0 H
Lower (decreasing H):  B(H) = B_s · (H − H_c)/(|H − H_c| + H_k) + μ_0 H
```

- `B_s` saturation, `H_c` coercive field, `H_k` controls loop knee
- Inversion symmetry: B_+(H) = −B_−(−H)
- **3 parameters** (B_s, H_c, H_k from datasheet) — extremely cheap
- Minor loops via a **transition algorithm** between branches
- Air gap modeled by extending the effective magnetic length
- Frequency + temperature dependence, wire skin effect add-ons
- Parasitic C + leakage L for the full transformer model

**Use Chan** for: SPICE circuit-level transformer/inductor models in
power-electronics simulation.  Cheapest hysteresis model that exists;
3 parameters straight from a datasheet.  **Limitation**: hyperbolic
branches are an analytic fit, NOT history-accurate for complex minor
loops (use Preisach/play if minor-loop fidelity matters).
Contrast with Chua (`chua_model`): Chan = static hyperbolic branches,
Chua = dynamic ODE.  Both circuit-oriented.

## Bouc-Wen model (structural / piezo) — Ikhouane-Rodellar 2007 (Wiley)

Originally for **mechanical / civil-engineering** hysteresis
(seismic isolators, structural dampers), and **piezoelectric /
magnetostrictive actuators**.  Smooth ODE form:

```
ż = α ẋ − β |ẋ| |z|^(n-1) z − γ ẋ |z|^n          (hysteretic variable z)
Φ(x, ẋ) = a k x + (1 − a) k z                    (restoring force/output)
```

- `α, β, γ, n` shape the loop; `a` = post-yield/pre-yield stiffness ratio
- **Smooth (C¹)** everywhere → ideal for control-system design + ID
- BIBO-stability classification (Ikhouane Ch 2): not all parameter
   sets give bounded output — must check
- Passivity conditions for closed-loop stability

**Use Bouc-Wen** for: piezo/magnetostrictive actuator hysteresis
compensation (the lab's `08_Bouc_Wen/改良Bouc-Wen` piezo paper),
multiphysics with mechanical coupling, control-oriented models.
**NOT** for soft-iron FE field analysis (no thermodynamic energy
balance, no vector magnetic extension).  The lab's energy-based
method IS the magnetic analog of Bouc-Wen's spring+friction idea
but done thermodynamically (`rheological_models`).

## Potter-Schmulian (magnetic recording media) — Potter 1971

Self-consistent magnetization in **thin magnetic recording media**
(hard-disk write process):
- 1D grid medium under a moving write head (Karlqvist fringe field)
- Head motion via ~100 sequential locations + convergent iterations
- Demagnetizing-field reduction from the high-permeability head;
   image term enforces head-surface boundary condition
- M-H model with hysteresis **squareness** S = M_r/M_s
- NRZ digital recording, single/double/quadruple transition analysis

**Use Potter-Schmulian** for: magnetic recording head/media
co-simulation (write process, transition broadening).  Entirely
outside the lab's EM-device scope — listed for completeness.

## LLG micromagnetics — Landau-Lifshitz-Gilbert

The **first-principles** model: solve the magnetisation dynamics on
a nm-scale mesh.

```
dM/dt = −γ M × H_eff + (α/M_s) M × dM/dt          (Gilbert form)
H_eff = H_applied + H_exchange + H_anisotropy + H_demag + H_thermal
```

- γ gyromagnetic ratio, α Gilbert damping
- Resolves Weiss domains + Bloch walls explicitly
- Hysteresis EMERGES from the energy landscape (no phenomenological
   loop input)
- Lab reference: `13_LLG/まぐね連載講座 マイクロマグネティクス
   シミュレーション III` (Japanese tutorial series)

**Use LLG** for: domain-wall dynamics, nm-scale magnetics, spintronics,
understanding the MICROSCOPIC origin of the macroscopic models above.
**Never** for device-scale FE — the mesh (nm) vs device (mm-m) scale
gap makes it ~10^9× too expensive.  LLG is the "ground truth" that
phenomenological models (play/stop/energy/Preisach) approximate.

## Decision summary

| Model | Domain | Cost | Lab use |
|-------|--------|------|---------|
| Chan | SPICE transformer | Trivial (3 params) | power-electronics circuit sim |
| Bouc-Wen | structural / piezo actuator | Low | piezo/magnetostrictive control |
| Potter-Schmulian | recording media | High (head iterations) | none (out of scope) |
| LLG | nm micromagnetics | Extreme | research / ground-truth only |
| **Lab core (Stop+Energy)** | **soft-iron device FE** | **moderate** | **everything else** |

## References

- J. H. Chan, A. Vladimirescu, X.-C. Gao, P. Liebmann, J. Valainis,
   "Nonlinear Transformer Model for Circuit Simulation", IEEE Trans.
   CAD 10(4):476-482, April 1991.
- F. Ikhouane, J. Rodellar, "Systems with Hysteresis: Analysis,
   Identification and Control using the Bouc-Wen Model", Wiley, 2007.
- R. I. Potter, R. J. Schmulian, "Self-Consistently Computed
   Magnetization Patterns in Thin Magnetic Recording Media", IEEE
   Trans. Magn. MAG-7(4):873, Dec 1971.
- LLG: see `13_LLG` folder; standard micromagnetics references
   (Brown 1963; OOMMF / mumax3 solvers).

Cross-ref: `chua_model` (the other circuit-oriented model),
`catalog` (full landscape), `rheological_models` (Bouc-Wen's
spring+friction done thermodynamically = lab core).
"""


JACQUES_MONOGRAPH = r"""
# Jacques 2018 PhD monograph — the canonical Energy-Based reference
# (K. Jacques, "Energy-Based Magnetic Hysteresis Models — Theoretical
#  Development and Finite Element Formulations", PhD, ULiège, 2018, 254p)

The definitive single document for the lab core method.  Derives the
energy-based (EB) model rigorously from thermodynamics, then gives
THREE discrete FE implementations + inversion methods with all the
numerical pitfalls spelled out.  Use this topic when you need the
WHY behind `lab_core_production`'s formulas.

## 1. Thermodynamic derivation (Ch 2)

**First law** (energy conservation, local form):
```
u̇ = Ẇ − div q
Ẇ = h·J̇          (magnetic work; empty-space J_0 = μ_0 h removed as
                   reversible, b = J_0 + J)
```

**Second law**: ṡ = ṡ_ext + ṡ_int = −div(T⁻¹q) + ṡ_int, with ṡ_int ≥ 0
(internal entropy production from irreversible processes like
hysteresis).

**Helmholtz free energy** f := u − Ts.  Combining 1st + 2nd law gives
the **Clausius-Duhem inequality**:
```
T ṡ_int = h·J̇ − grad T·(T⁻¹q) − (ḟ + s Ṫ) ≥ 0
```
Neglecting thermal terms → dissipation `D = (h − h_rev)·J̇ ≥ 0`.

## 2. Reversible part — anhysteretic via free energy (Ch 2 §2)

```
h_rev = ∂f/∂J          (reversible "spring" field)
J = J_an(h) ↔ h = J_an⁻¹(J)        (anhysteretic 1-to-1 law)
Δf = ∫ h·dJ            (free energy = area under polarization curve)
```
Isotropic: `J = J_an(|h|)·h/|h|`.

## 3. Irreversible part — dry-friction pinning (Ch 2 §3)

```
D = |κ J̇| = h_irr·J̇       (κ = pinning field, A/m = friction magnitude)
h_irr = ∂_J̇ D             (subgradient of the non-smooth |κ J̇|)
```

**The origin of hysteresis** is the NON-UNIVOCITY of h_irr at J̇ = 0:
- At J̇ = 0: subdifferential = the whole ball |h_irr| ≤ κ (many possible
   gradients → memory)
- At J̇ ≠ 0: subdifferential = singleton {κ J̇/|J̇|} (classical gradient)

## 4. Single-cell equilibrium (Ch 2 §4)

```
h − h_rev − h_irr = 0
⇔  h − ∂f/∂J ∈ ∂D(J̇) = { |h_irr| ≤ κ        if J̇ = 0  (FROZEN)
                          { h_irr = κ J̇/|J̇|   otherwise (SLIPPING)
```

Mechanical analogy: spring (h_rev) ∥ friction slider (h_irr).  The
"grey sphere" of radius κ centered at h_rev: while tip of h stays
inside, |h_irr| < κ → J̇ = 0 → magnetization frozen.  Slider unlocks
when |h − h_rev| = κ.

The elementary J-h curve is **non-differentiable at reversal points**
(angular points where the regime switches frozen↔slipping) — this is
the source of all numerical difficulty.

## 5. Multi-cells model + homogenization (Ch 2 §5)

Single cell gets the MAJOR loop right but minor loops are unrealistic.
Fix: decompose into λ-regions weighted by probability ζ(λ), each with
its own pinning κ*(λ):
```
f*(J*(λ)) = ζ(λ) f(J*(λ)/ζ(λ))
```
Each λ-region obeys the same evolution law (2.29), independently.

**Two homogenization strategies**:
- (a) **Energy-balance homogenization** (PREFERRED, used by all EB
   literature): ⟨f⟩ = ∫ ζ(λ) ḟ(J*(λ)/ζ(λ)) dλ — net magnetization =
   infinite sum of cell contributions, exact thermodynamic balance
- (b) **Magnetic-field homogenization**: h_rev = ∫ ζ(λ) h*_rev(λ) dλ —
   loses the anhysteretic 1-to-1 at macroscale (artificial dissipation
   term), but useful for parameter identification (Ch 6)

Choosing one always involves SOME approximation; the lab uses (a).

## 6. THREE discrete implementations (Ch 3) — the production choice

| Approach | Idea | Cost | Robustness |
|----------|------|------|------------|
| **Vector Play** (simplified differential) | explicit update rule (Henrotte 2006 simplification: J̇ direction ∥ h − h_rev(p) at PREVIOUS step) | cheapest, O(N) explicit | approximate (explicit) |
| **Variational** (var, functional min) | minimize Ω^k = g_L^k + D̃^k per cell | medium | robust if line-search good |
| **Angle-searching** (diff, full differential) | solve for h_rev on sphere radius κ, find angle α s.t. h_irr ∥ J̇ | medium | exact (no extra approx beyond time disc.) |

### Variational approach (var) details

Landau free energy `g_L^k(h,J^k) = f^k(J^k) − h·J^k` (analog of
plasticity / Saint-Venant kinematic hardening; h↔stress, J↔strain).
Pseudo-potential `D̃^k(J^k, J^k_(p)) = κ^k |J^k − J^k_(p)|`.  Then:
```
J^k = argmin_{J^k} Ω^k(h, J^k, J^k_(p))    where  Ω^k = g_L^k + D̃^k
b = μ_0 h + Σ_k J^k(h, J^k_(p))
```
**Only minimize when |h − h^k_rev(p)| > κ^k** (else frozen, J^k = J^k_(p),
zero cost).  Minimization technique matters:
- Descent direction: steepest descent (sd) SLOW vs conjugate gradient
   (cg) FAST (cg recommended)
- Line search: naive < Wolfe < **Brent** (Brent finds the true minimum)
- **Numerical loss of significance** at ~1e-8 T scale in J → use a
   2nd-order Taylor approx near the minimum when step gets too small

### Angle-searching approach (diff) details

If |h − h^k_rev(p)| > κ^k, the solution h^k_rev lies on a sphere of
radius κ^k centered at h:
```
h^k_rev(α^k) = h − κ^k ê_h_irr(α^k),   ê_h_irr = (cos α^k, sin α^k)  (2D)
```
Find α^k so the irreversible-field direction ∥ J̇^k, with
```
e_J̇ ≈ [J_an(h^k_rev) − J_an(h^k_rev(p))] / |...|
```
(3D needs two angles.)  No approximation beyond time discretization →
the most accurate of the three, but requires a root-find on α.

## 7. Inversion (Ch 5) — Newton-Raphson + analytical permeability tensor

For the **inverse** operator H = H(B) (vector-potential FE), Jacques
2018 Ch 5 gives:
- Classical Newton-Raphson with the **analytical differential
   permeability tensor** ∂B/∂H (closed-form from the cell equations)
- Practical convergence + relaxation technique for the angular points
- Approximated / quasi-Newton variants (numerical Jacobian) for when
   the analytical tensor is too expensive

This is the predecessor of the Egger 2025 Schur-complement O(N)
inverse (`forward_inverse_variational`, `efficient_eval`), which
superseded the naive Newton inversion.

## Why this monograph matters

It is the SINGLE most complete derivation+implementation reference for
the lab core method.  The three implementations map to:
- **Vector Play** → `rad.MatPlayHysteresis` (the fast explicit path)
- **Variational** → Francois-Lavet 2013 + Egger 2025 (the FE-robust path)
- **Angle-searching** → the "Full differential approach" / Henrotte
   2006 exact update

When deciding which to implement in Radia: Variational (var) is the
production choice because its convex functional gives provable
convergence (cf. `forward_inverse_variational`); Vector Play is the
cheap explicit fallback; Angle-searching is most accurate per step
but needs a robust root-finder.

## Reference

K. Jacques, "Energy-Based Magnetic Hysteresis Models — Theoretical
Development and Finite Element Formulations", PhD thesis, Université
de Liège, 2018 (254 pages).  https://orbi.uliege.be/handle/2268/229596
File: `W:/03_文献・論文/00_電磁界解析/30_磁気特性/02_ヒステリシスモデル/
05_Energy_Based/Energy-BasedMagnetic Hysteresis Models.pdf`.

Cross-ref: `lab_core_production` (the formulas this monograph derives),
`energy_based` (lineage), `forward_inverse_variational` (Egger 2025
successor to Ch 5 inversion), `rheological_models` (the spring+friction
analogy of §4), `mathematical_foundations` (Visintin variational view).
"""


CONGRUENCY_BINPUT_SELECTION = r"""
# Why the lab core is B-input ENERGY Stop: the congruency selection argument
# (Sugahara Lab flagship logic; Hane-Sugahara-Morishita-Ahagon 2026)
#
# Accepted: IEEE Trans. Magn., "Experimental Verification of Congruency
# Property in Magnetic Material and Its Reproducibility Using B-input Play
# Model" (Manuscript TMAG-25-11-0894, COMPUMAG 2025 regular paper).
# File: W:/02_学会資料/2025年度/2025_11_21_IEEE_Magnetics_論文投稿/ (Final
# Submission); copy stored in 02_ヒステリシスモデル/03_Play_model/.

## TL;DR (the selection chain)

Requiring BOTH (1) reproduction of the experimentally-measured congruency
behavior AND (2) thermodynamic (energy) consistency UNIQUELY selects the
**B-input ENERGY-based STOP** model.  Neither requirement alone forces it;
together they leave exactly one survivor.

```
(1) match measured congruency  -->  B-input  (forced, NOT a free choice)
(2) thermodynamic consistency  -->  convex energy U_k
                                     -> monotone shape function f_k = U_k'
among B-input {play, stop}:
    B-input STOP : f_k monotone increasing   -> convex U_k OK  -> energy-based
    B-input PLAY : f_k NON-monotone (measured)-> no convex U_k -> NOT energy-based
=> B-input energy-based STOP is the unique model satisfying (1) AND (2).
```

## Step 1 - the experiment forces B-input (Hane 2026, the accepted paper)

Ring core of non-oriented silicon steel **35A360**, minor loops measured
under several excitation conditions (l=0.550 m, S=1.00e-3 m^2, N1=615
exciting turns, Ns=700 search turns; Hioki PW6001 + Denshijiki BH analyzer).

- **H-axis congruency IS observed**: minor loops swept between the same two
  H extrema have identical shape regardless of operating point (Figs 10-13).
- **B-axis congruency is NOT observed** ("B-axis incongruency", Fig 14):
  minor loops between the same two B extrema change shape with operating
  point.
- **Physical reason** (paper Sec III): H is the external field set by the
  winding current -- it knows nothing about the material -> H-axis
  congruency.  B is the average flux density INSIDE the material -- it is
  shaped by the magnetization process -> B-axis incongruency.
- **The two congruencies cannot hold simultaneously** (paper conclusion).

A congruency-based model (Preisach / play / stop) has congruency built in
along its INPUT axis.  To match the measurement (H-axis congruent, B-axis
NOT), the model input must be **B**:
  - B-input model -> H-axis congruency  (matches measurement)  OK
  - H-input model -> B-axis congruency  (contradicts measurement)  WRONG
The accepted paper proves this directly: B-input play reproduces both
H-axis congruency and B-axis incongruency; **H-input play fails** -- it
forces a B-axis congruency that does not exist, giving wrong minor-loop
shapes, an unphysical mid-curve overshoot (Fig 18b), and wrong loop loss
(Tables 2-6).  So **matching reality forces B-input** -- this is dictated by
the measured congruency, not an "input-argument preference".

## Step 2 - energy-consistency, among B-input, forces STOP (the lab argument)

Lab production wants a **thermodynamically consistent** model (energy known
at all times, convex inversion, Picard fast inverse -- see
`lab_core_production`, `energy_based`).  Energy consistency requires the
per-cell energy density U_k to be **convex**, i.e. its derivative (the
shape function) f_k = U_k' must be **monotone non-decreasing, non-negative**.

Inside the B-input class that Step 1 forces:

- **B-input STOP**: the natural energy form.  Shape functions for the stop
  operator are monotone increasing -> U_k convex -> a valid energy density
  exists -> the model IS energy-based.  (Exactly the
  Bergqvist/Henrotte/Francois-Lavet/Egger B-input energy/stop lineage --
  see `energy_based`, `forward_inverse_variational`.)
- **B-input PLAY** (the inverse model used in the accepted paper): the shape
  functions identified from measured loops via the Everett function are
  **NOT monotone** -- they go negative and intersect (accepted paper
  Fig 6(b) negative branch; Fig 7(b),(c) crossing curves; the text
  attributes the H-input overshoot to "intersection of the shape functions
  for different play hysterons").  A non-monotone f_k has no convex U_k, so
  **B-input play cannot be recast as an energy-based model.**  (Deeper
  structural reason: the B-input Play state domain [B-eta_k, B+eta_k]
  translates rigidly with B, so no fixed-domain coenergy U_k(p_k) exists at
  all -- Bergqvist's variational principle is structurally inapplicable to
  Play, before convexity is even at issue.  See `bergqvist_binput_stop`.)

Therefore, within the B-input class that Step 1 forces, **only STOP admits
an energy-based (thermodynamically consistent) realization.**  That is why
the lab core is **B-input energy-based STOP**, not B-input play.

## CRITICAL scoping (do NOT overstate at IGTE)

State it as a conjunction, or the energy-based community
(Henrotte / Egger / Jacques -- IGTE regulars) will object on sight:

- Do NOT say "play cannot be energy-based".  **Energy-based PLAY exists**
  (Bergqvist 1997 vector play, Henrotte 2006, Jacques 2018) -- but it is
  **H-input** (co-energy form).  By the convex duality
  (energy<->co-energy = stop<->play = B-input<->H-input; see
  `forward_inverse_variational`) the energy-consistent play is the H-input
  one.
- The accepted paper shows the H-input model (= that energy-consistent play)
  fails the congruency test (forces B-axis congruency).
- Correct statement (the conjunction):
  **"The B-input play (inverse) cannot be made energy-based (non-monotone
  shape function); the energy-based play that DOES exist is H-input, which
  fails the measured congruency.  Hence energy-consistency + measured
  congruency leaves only B-input STOP."**

## Positioning vs prior congruency literature (answers TMAG Reviewer 3)

Reviewer 3 of the accepted paper flagged heavy prior art on congruency-based
models and warned the novelty was "limited" ("identification not new,
measurements not new").  The same risk hits any follow-up.  Prior art:

- Davino-Giustiniani-Visone 2008 (Physica B 403:414) -- H vs B input to a
  Preisach solver, congruency conjectures.
- Zirka-Moroz-Marketos-Moses 2004 (IEEE TMAG 40:390) -- congruency-based
  history-dependent models (HDHM) for transient simulation.
- Zeinali-Krop-Lomonova 2020 (IEEE TMAG 56:6700904) / 2021
  (IEEE TMAG 57:7300304) -- (anisotropic) congruency-based vector hysteresis
  for non-oriented laminated steel.

What ALL of these share: they are **phenomenological congruency-based
models**.  None is **thermodynamically consistent** (none is built from a
convex energy/co-energy with a guaranteed energy balance).  The accepted
Hane 2026 paper itself repositioned to "experimental verification" to clear
this prior art.

**Lab differentiator for the NEXT contribution = thermodynamic consistency
+ uniqueness.**  Not "another congruency-based model for NOES" (Zeinali owns
that), but: "the measured-congruency requirement and the energy-consistency
requirement TOGETHER select a UNIQUE model -- B-input energy STOP -- which
reproduces the measured H-axis congruency / B-axis incongruency that
phenomenological congruency-based models reproduce, while ALSO carrying a
thermodynamically consistent energy balance that they do not."  Keywords to
hammer: **uniqueness (selection)** + **thermodynamic consistency**.

## Out of scope (lab decision 2026-05-28)

- **Open boundary / Kelvin transformation is NOT part of this story.**  The
  novelty is the selection argument + energy-consistency on the material
  LAW, demonstrated on the measured congruency dataset.  Do NOT bolt on a
  Kelvin open-boundary FE demo as the headline.

## Concrete deliverable (so it is more than a positioning essay)

Show that **B-input ENERGY STOP** (convex U_k, monotone f_k) reproduces the
SAME measured congruency battery (Hane 2026 Figs 10-19, 35A360 minor loops)
that B-input PLAY did -- i.e. energy-consistency costs nothing in congruency
fidelity, while H-input (the energy-consistent play) still fails.
The differentiator vs the Hane 2026 play paper is **thermodynamic
consistency** (energy balance, non-negative dissipation), NOT a dc-bias
accuracy win. Demonstrate the energy stop on the SAME congruency battery
(symmetric + H-axis-congruency minor loops) and show energy-consistency
costs nothing in congruency fidelity.
DO NOT claim a dc-bias advantage for the (input-independent) energy stop:
on asymmetric / dc-biased loops it has equal vertical chords -- the SAME
limitation as pure stop (verified, see `stop_vs_play_silicon_steel` and
demo_equal_vertical_chords.py), and worse than play (which has
congruency). Recovering dc-bias accuracy needs the input-dependent w(B)
extension, and whether w(B) preserves the fixed-domain energy is open.

## Cross-references

- `lab_core_production` -- the B-input energy Stop implementation recipe
- `energy_based` -- Bergqvist/Henrotte/Jacques energy lineage
- `forward_inverse_variational` -- Egger 2025 energy<->co-energy duality
  (= the B-input/H-input <-> stop/play map)
- `stop_vs_play_silicon_steel` -- Matsuo 2003 (equal-vertical-chords limit,
  dc bias): why input-independent stop (incl. the energy/variational
  B-clip form) cannot represent dc-biased silicon steel; needs w(B)
- `bobbio_1997_unification` -- PTM/STM; B-input play ~= stop
- `preisach` -- congruency as a Preisach property
- `zirka_hdhm` -- congruency-based HDHM (Zirka 2004, prior art)

## Reference

Y. Hane, K. Sugahara, H. Morishita, A. Ahagon, "Experimental Verification of
Congruency Property in Magnetic Material and Its Reproducibility Using
B-input Play Model", IEEE Trans. Magn. (accepted; COMPUMAG 2025), Manuscript
TMAG-25-11-0894.  Companion: H. Morishita, K. Sugahara, Y. Hane, A. Ahagon,
"Identification of the Shape Function of a Play Model Using Parametric
Polynomial Approximation", COMPUMAG 2025, MP2B-13.
"""


BERGQVIST_BINPUT_STOP = r"""
# Bergqvist energy-based B-input Stop model (Sugahara Lab, IGTE 2026 digest)
# Title: "Bergqvist Energy-Based B-Input Stop Model" (Sugahara, Hane)
# Target: 22nd IGTE Symposium, Graz, Sept 2026, "Magnetic Losses in
#   Laminated Cores" special session, 1-page A4 digest.
# This is the v5-Bergqvist framing (variational update + EM identification),
#   a rewrite of the v4 "adaptive saturation heights" digest.

## The structural thesis (why STOP, not Play -- the rigorous version)

This SHARPENS `congruency_selection` Step 2.  The reason B-input Play
cannot be energy-based is NOT merely "the identified shape function is
non-monotone" (that is a symptom); it is a DOMAIN argument:

- **B-input PLAY**: state p_k is confined to [B - eta_k, B + eta_k], which
  TRANSLATES RIGIDLY with B.  The admissible domain has no fixed reference,
  so NO B-independent per-hysteron coenergy U_k(p_k) = int f_k exists.
  Bergqvist's variational principle (which needs a fixed-domain convex
  potential per hysteron) is therefore STRUCTURALLY INAPPLICABLE to
  B-input Play -- before any question of convexity even arises.
- **B-input STOP**: state s_k = B - p_k is confined to [-eta_k, +eta_k],
  which is B-INDEPENDENT (fixed).  Hence U_k(s_k) = int_0^{s_k} g_k is a
  genuine fixed-domain potential and Bergqvist's variational principle
  applies DIRECTLY.

So requiring (i) B-input (from the measured H-axis congruency / B-axis
incongruency, Hane 2026) AND (ii) the Bergqvist energy formulation selects
the STOP model UNIQUELY, by the fixed-domain property -- not by a
shape-function fitting accident.  (This is the cleanest statement of the
"play cannot be energy-based" point: scope it to B-input play via the
moving domain; the H-input energy-based play of Bergqvist/Henrotte/Jacques
still exists -- see `congruency_selection` CRITICAL scoping.)

## The Bergqvist variational update (replaces the clip rule)

Standard Stop uses a hard clip s_k <- clip(s_k + dB, [-eta_k, +eta_k]).
Energy-based Stop replaces it with a per-hysteron variational minimization:

    s_k* = argmin_{s in [-eta_k, +eta_k]} [ U_k(s)
             + (s - s_p - dB)^2 / (2 gamma) ]

where s_p is the previous state and gamma -> 0 recovers the clip rule.
Bisection solves the stationarity condition.  The penalty term
(s* - c)^2 / (2 gamma) >= 0 guarantees NON-NEGATIVE DISSIPATION at every
step (thermodynamic consistency by construction).  This is Bergqvist 1997
dry-friction-like pinning written per Stop hysteron.

## Self-consistent identification by EM

Circularity: g_k both (a) produces the output H = sum_k g_k(s_k) AND (b)
IS the gradient dU_k/ds_k that drives the variational state update.  g_k is
the unknown yet it drives the states used to identify it.  The only
constraint is g_k monotone non-decreasing (= U_k convex).  Resolve by
Expectation-Maximisation:

- **E-step**: with current g_k, solve the variational update (gamma = 1e-5)
  over the training waveform to collect states {s_k(t)}.
- **M-step**: regress the target H on {s_k} by bounded-variable least
  squares (BVLS); g_k a degree-4 I-spline (n_basis = 13) with NON-NEGATIVE
  coefficients (enforces monotone non-decreasing g_k).
Iterate E/M to convergence.

## Validation (Potter-Schmulian synthetic benchmark)

Analytical Potter-Schmulian loop (S = 0.9, M_s = 1.91 T, H_c = 500 A/m) as
the identification target; N_a = 20 symmetric loops, dB = 0.095 T, K = 25
cosine-spaced thresholds eta_k in [0.045, 2.10] T.

- Loop **rms = 13.7 A/m** across the 20 amplitudes.
- All K = 25 identified shape functions g_k are monotone non-decreasing,
  confirming g_k = dU_k/ds_k exactly -> field output equals coenergy
  gradient -> non-negative dissipation every step -> a thermodynamically
  CLOSED B-input description with no post-hoc correction.

(Potter is a synthetic test bench, not a target material; measured
silicon-steel validation against Hane data is reserved for the long paper.)

## Route (i): input-dependent shape keeps energy-consistency AND breaks EVC (2026-05-29)

The v5 digest's stop uses an input-INDEPENDENT g_k. Congruency follows
the INPUT axis, so the B-input stop has H-axis congruency just like
B-input play (loops within the same B-range are congruent for both --
hysteron states follow the same history within the range). What
DISTINGUISHES the stop (Matsuo-Shimasaki 2005 N&S: input-indep stop <=>
wiping-out + equal vertical chords) is the ADDITIONAL equal-vertical-
chords constraint: same B-amplitude at any dc bias => same H-chord.
Silicon steel HAS H-axis congruency but LACKS EVC (Matsuo 2003:
asymmetric loops are wider in H). So the digest's stop, as-is, cannot
represent dc-biased
silicon steel -- and the symmetric Potter/loop validation CANNOT expose
this (all loops at the same bias). The energy/variational wrapper does
NOT remove EVC (the B-clip box is bias-independent): verified, chord
bias-independent to 0.0-0.1% (demo_equal_vertical_chords.py).

FIX (route i): scale the coenergy by a weighting on the FIXED box,
Psi(s,B)=w(B)U0(s), and take H as its B-gradient:
    H = dPsi/dB = w(B)H0(s) + w'(B)U0(s)         (NOT Matsuo's bare w*H0)
    D = H Bdot - Psidot = w(B) sum_{k clipped} g0(+-eta) Bdot >= 0
The w' terms cancel => D = w*D0 >= 0 (w>=0): thermodynamic consistency
preserved. The box [-eta_k,+eta_k] stays B-independent => the domain
argument (why STOP, not B-input PLAY, can be energy-based) SURVIVES; only
the energy magnitude scales with w(B). The chord becomes bias-dependent
Delta_h = w(B0)Delta_h0 + w'(B0)Delta_U0 => EVC broken => matches the
material. (The bare Matsuo H=w*H0 drops w'U0, is not a free-energy
gradient, and is not guaranteed D>=0.)

KAN usage: parametrise Psi (monotone g0 + positive w(B), e.g. RBF +
softplus) and derive H = dPsi/dB by autograd; do NOT free-fit H or D>=0
is lost. URN (KAN-inspired) is NOT applicable -- it is linear
frequency-domain relaxation (impedance), not rate-independent hysteresis;
only its KAN edge-function machinery is reused for w(B).

Validation (2026-05-29):
- synthetic strong-EVC recovery on real-data g0: w_KAN ~ w_true to ~1%,
  fit 100.95 -> 7.83 A/m (92%), chord 61->101 A/m vs flat 61 for w=1
  (demo_kan_w_recovery.py);
- real symmetric loops (35A360-class B-input set, 30 amplitudes
  0.05-1.5 T): energy-consistent monotone g0 fits 3.1% rel RMS
  (fit_binput_symmetric.py) -- energy-consistency costs ~nothing on
  symmetric data;
- SPCC real EVC violation: pure stop mid-B relative error sym 16% vs
  asym 42% (fit_spcc_stop.py).
OPEN: a scalar w(B) is insufficient on messy real data (trades symmetric
for asymmetric fidelity); the full input-dependent SHAPE g_eta(s,B) =
w(B) g0(s) / a 2-D KAN on (s,B) is the next step. Full write-up: docs
HYSTERESIS_THEORY.md §7.5.9.

## Two digest framings on record (do not confuse / do not blend)

| Version | Title | Mechanism | Identification |
|---|---|---|---|
| v4 (logged in CHANGES.md) | "Energy-Consistent B-Input Stop" | hard clip + per-hysteron coenergy psi_k = int g_k; ADAPTIVE saturation heights (half-grid offset + arc-length knee refinement, N=23); H-input/MMM use + Egger Schur | constrained lsqlin, descending branch |
| v5-Bergqvist (THIS topic; copy folder, NOT yet logged) | "Bergqvist Energy-Based B-Input Stop" | VARIATIONAL update (eq above) + Bergqvist dry-friction; the DOMAIN argument front-and-center | EM (E: variational, M: BVLS I-spline) |

v5 is the cleaner story (domain argument + a genuine Bergqvist variational
update + self-consistent EM), rms 13.7 A/m.  v4 leads with the H-input/MMM
use case and adaptive hysteron placement.  Pick ONE for the IGTE
submission; do NOT blend (1-page limit; the comparison does not fit --
CHANGES.md lesson 3).

## Cross-references

- `congruency_selection` -- why B-input (measured congruency) + the
  selection logic this topic sharpens with the domain argument
- `energy_based` -- Bergqvist/Henrotte/Jacques lineage
- `forward_inverse_variational` -- Egger 2025 energy<->co-energy duality
- `lab_core_production` -- the production B-input energy Stop recipe
- `stop_vs_play_silicon_steel` -- Matsuo 2003 equal-vertical-chords limit

## References

- A. Bergqvist, "Magnetic Vector Hysteresis Model with Dry Friction-Like
  Pinning", Physica B 233(4):342-347, 1997.
- T. Matsuo, D. Shimode, Y. Terada, M. Shimasaki, "Application of Stop and
  Play Models to the Representation of Magnetic Characteristics of Silicon
  Steel Sheet", IEEE Trans. Magn. 39(3):1361-1364, 2003.
- H. Egger, F. Engertsberger, A. Schafelner, "Forward and Inverse
  Energy-Based Hysteresis Operators", arXiv:2507.15289, 2025.  (NOTE: a
  THIRD Egger-2025 paper, distinct from arXiv:2507.14521 vector-potential
  and arXiv:2410.11705 energy-consistent vector operators.)
- R. I. Potter, "Self-Consistently Computed Magnetization Patterns in Thin
  Magnetic Recording Media", IEEE Trans. Magn. 7(4), 1971.
- Y. Hane, K. Sugahara, "Experimental Verification of Congruency Property
  ... B-Input Play Model", IEEE Trans. Magn. (accepted), 2026.
"""


PEELING_IDENTIFICATION = r"""
# Peeling identification of play / stop shape functions
# (Sugahara Lab / Ahagon -- analytic, NON-least-squares calibration)
# Code: W:/999_菅原賢悟/19_磁気ヒステリシス/2024_IGTE_共同研究/
#       VectorPlayModel/supplement/playmodel.py  (shapeFunction)
#       + ...VectorPlayModel/supplement/阿波根@形状関数同定.py

## What it is

"Peeling" (剥ぎ取り) identification recovers the PER-HYSTERON shape
functions of a B-input play (or stop) model DIRECTLY from measured
first-order reversal curves -- the descending branches of symmetric
B-input loops measured at an equidistant ladder of amplitudes
B_max = dB, 2dB, ..., N*dB.  It is a closed-form difference scheme:
NO least-squares, NO matrix inversion, NO iteration.

It is the play/stop analogue of the Preisach Everett-function
identification (mixed second differences of FORC data): the difference
array `mu` plays the role of the Preisach density.

## Why "peeling"

The nested descending curves form a triangular structure: the loop of
amplitude k*dB activates exactly k hysterons (radii dB ... k*dB).
Going from amplitude k*dB to (k+1)*dB adds exactly ONE more active
hysteron.  Differencing successive-amplitude curves therefore "peels
off" one hysteron's contribution at a time, isolating each shape
function.

## Algorithm (playmodel.py `shapeFunction(Hdata)`)

Input `Hdata` = descending-curve look-up table (HdataLUT), shape
(2N+1, N): column j is the descending branch of the symmetric loop of
amplitude B_max,j sampled at equidistant B from +B_max to -B_max
(step dB).  Lab column order: col 1 = 1.5 T ... col 30 = 0.05 T.

1. First / cross differences -> "density" mu:
     mu[i,0]   = Hdata[0,i] - Hdata[1,i]
     mu[i,j+1] = (Hdata[j+1,i]   - Hdata[j+2,i])
               - (Hdata[j,i+1]   - Hdata[j+1,i+1])      (recursive)
   plus a closing term so each hysteron row sums consistently.
2. Expand mu into the staircase (each entry halved), flip, and
   cumulatively sum: descending shape Hdown = -cumsum(mu2),
   ascending shape Hup = -flip(Hdown), assemble the full ODD shape
   function  Hfunc = [ Hup ; 0 ; Hdown ].
3. Hfunc[:,i] is the shape function of hysteron i, tabulated on the
   Bplay axis (consumed by playHysteron -- see `vector_play_model`).

## Why NOT least squares  (verified 2026-05)

Per-hysteron LEAST-SQUARES fitting of the play model on the play
variable p_k = B - s_k is ILL-CONDITIONED: every p_k shares the common
input B, so the design columns are strongly collinear.  Naive ridge
LSQ therefore mis-attributes the fit and can under-perform even the
stop model -- a pure CONDITIONING ARTEFACT, NOT a statement that "play
is worse than stop".  The peeling difference scheme sidesteps the
collinearity entirely (it reads each hysteron off the FORC ladder), so
it is the CORRECT way to calibrate a play model.  (Empirically: a crude
per-hysteron I-spline LSQ play landed ~30% while the well-conditioned
stop fit was ~21% on the same SPCC data -- the gap was the method, not
the operator.  Use peeling, not LSQ, for play.)

## Requirements / caveats

- EQUIDISTANT thresholds (dB ladder) are mandatory -- the difference
  scheme assumes uniform amplitude spacing.
- The descending curves must be clean (monotone, de-noised).  The lab
  pipeline smooths each measured loop and re-parametrises the shape as
  a constrained polynomial before tabulating (pipeline below).
- Reproduces the identification curves EXACTLY by construction; minor /
  biased loops are then PREDICTED via the play operator (congruency).

## Full lab calibration pipeline
## (W:/.../2024_10_16_B入力用ループ/prototype_step0-4.m)

step0: load measured symmetric loops data_XXXXmT.mat (B_max =
       0.05:0.05:1.5 T; raw 5 Hz columns: t, V_exc, I_exc, V_search,
       H, B); split ascending / descending, polyfit each, average ->
       smoothed symmetric loop Hm(Bm); store tip values (B_max, H_max).
step1: fit the loop-TIP locus H0(B0) by a rational (Pade) function
       B2H(b) = polyval(p,b)/polyval(q,b) (anhysteretic-like backbone
       through the loop tips).
step2: normalise each loop x = B/B_max, y = H/H_max and fit a degree-17
       polynomial shape p(x) with p(+-1)=+-1 and a MONOTONICITY
       constraint (no p' roots in [-1,1]); build p_interp = the loop
       shape INTERPOLATED across amplitude.  NB: the play model DOES
       interpolate its normalised loop shape across amplitude here --
       a nuance vs the loose claim "play cannot interpolate shape".
step3/4: tabulate descending curves h(b) = B2H(B_max)*polyval(p_interp,
       b/B_max) into HdataLUT (and the dual BdataLUT for H-input).
playmodel.py: shapeFunction(HdataLUT) -> Hfunc (peeling); then forward
       via playHysteron (see `vector_play_model`).

## Cross-reference

- `vector_play_model`  - the forward model that consumes Hfunc
- `play`               - general play-model theory
- `preisach`           - Everett-function (FORC 2nd-difference) analogue
- `input_dependent_shape`, `stop_vs_play_silicon_steel` - why silicon
  steel needs congruency (play) not equal-vertical-chords (stop)
- `bergqvist_binput_stop` - the lab energy STOP (EVC + thermo + Egger
  inverse); play and stop are complementary (p = B - s)
"""


VECTOR_PLAY_MODEL = r"""
# Canonical B-input vector play model (Sugahara Lab / Ahagon)
# Code: W:/999_菅原賢悟/19_磁気ヒステリシス/2024_IGTE_共同研究/
#       VectorPlayModel/  (+ supplement/playmodel.py)
#       and .../2024_10_16_B入力用ループ/  (data + calibration pipeline)

## Constitutive form (scalar, B-input)

  H(B, history) = sum_{i=1}^{Np} Hfunc_i( p_i )

with the per-hysteron PLAY (backlash) operator
  p_i = clip( p_i^prev , B - zeta_i , B + zeta_i )
      = max( min( p_i^prev, B + zeta_i ), B - zeta_i ).

- zeta_i = play radius of hysteron i, EQUIDISTANT: zeta = 0, dB, 2dB,
  ... up to ~B_max  (lab: np.arange(0, B_max, dB)).
- Hfunc_i = per-hysteron shape function, identified by PEELING from the
  descending-curve LUT (see `peeling_identification`); tabulated on a
  Bplay axis, evaluated by 1-D interpolation at p_i.
- Saturation: B is clipped to +-B_max (Blimit); for |B| > B_max a
  linear extrapolation term is added to hysteron 0.

Forward cost: O(Np) per step (one clip + one table interp per
hysteron).  Reference: VectorPlayModel/supplement/playmodel.py
(`playHysteron`).

## Key property: CONGRUENCY

The play operator yields the CONGRUENCY property: dc-biased minor loops
of a given B-excursion have the SAME shape regardless of bias.  This
MATCHES H-axis-congruent silicon steel (Hane 2026, NOES 35A360), so the
B-input play model reproduces dc-biased minor loops well -- the regime
where the input-independent STOP model (equal vertical chords) fails.

## Vector / rotational extension

The isotropic vector play model superposes scalar play operators over
directions (rotational hysteresis, elliptical / circular B traject-
ories).  Lab reference implementations: VectorPlayModel/prototype4*.m
(vector), VectorPlayModel/BQM1, BQM2 (Alex "BQM" circular / ellipsoid /
minor-loops-with-bias), Evaluate_Rotation_Hys_Loss.m; Jiles-Atherton
cross-check in JilesAtherton.m.

## Calibration (two-stage; see `peeling_identification`)

1. Build HdataLUT from measured symmetric loops (prototype_step0-4.m:
   smooth -> tip rational fit B2H -> normalised constrained-poly shape
   p_interp -> descending-curve table).
2. Peeling difference scheme -> Hfunc (shapeFunction in playmodel.py).
NO least-squares (ill-conditioned for play -- see
`peeling_identification` "Why NOT least squares").

## Data (W:/.../2024_10_16_B入力用ループ/)

- data_XXXXmT.mat : measured symmetric loops, B_max = 50 ... 1500 mT in
  50 mT steps (cols: t, V_exc, I_exc, V_search, H, B; raw 5 Hz).
- BmHm.mat, B0H0.mat : loop tips + tip-locus rational fit.
- p_interp.mat : normalised loop shape interpolated across amplitude.
- HdataLUT.mat / BdataLUT.mat : descending-curve LUTs (B-input/H-input).

## Relation to the lab energy STOP (IGTE digest model) + inverse trichotomy

Play and stop are complementary (p_i = B - s_i), so the play model is
EXACTLY an input-dependent stop:
   H = sum f_k(p_k) = sum f_k(B - s_k) =: sum g_k(s_k;B), g_k(s;B)=f_k(B-s).
But g_k is non-monotone in s (~83% of hysterons), so the coenergy
U_k = int g_k is NON-convex.  The clean Egger-Schur (vector-FE) inverse
needs CONVEX U_k (monotone g).  Hence the trichotomy (2026-05 study,
"does the STOP reproduce the canonical PLAY?", oracle = peeling play
validated to 0.01-0.60% on the symmetric loops):

  | construction                       | play-equiv. | inverse              |
  | input-dependent stop g=f(B-s)      | EXACT       | non-convex -> Newton |
  | monotone Bergqvist stop            | ~9.6%       | convex -> Egger Schur|
  | route-(i) input-dependent convex   | shrinks 9.6%| convex (open)        |

=> EXACT play-equivalence and the clean convex inverse CANNOT both hold;
   the ~9.6% is the *cost of the easy (convex) inverse*.  Concretely, a
   per-hysteron stop calibrated on SYMMETRIC play loops and tested on
   dc-biased minor loops reproduces the canonical play to: non-monotone
   stop 7.6%, monotone Bergqvist stop 9.6% (mean biased; denser-high-eta
   K=60 best config; hardest small loops amp=0.1T: 7.4% / 9.9%).  The
   monotone (digest) model thus buys thermo-consistency (D>=0) + Egger
   O(K) inverse for ~+2pt vs the free stop.  CAVEAT: the SCALAR B<->H
   inverse (1 point) is easy for ANY model with H monotone in B (total
   dH/dB>0 even if individual f_k'<0), so convexity is required only for
   the vector-FE convex-Schur inverse.  An earlier SYNTHETIC PI-oracle
   gave a spurious "50% wall" -- artefact of an adversarial anhysteretic-
   heavy mu(eta); real silicon steel is far milder (~10%).
   (Lab scripts: port_canonical_play.py, does_stop_reproduce_play.py,
   refine_and_plot.py, distinguish_mono_vs_nonmono.py.)

## Cross-reference

- `peeling_identification` - how Hfunc is calibrated (the companion topic)
- `play` - general play theory;  `vector` - E&S / vector Preisach
- `congruency_selection`, `bergqvist_binput_stop` - the lab play-vs-stop
  selection logic and the energy STOP digest
"""


def get_hysteresis_models_knowledge(topic: str = "lab_core") -> str:
    """Dispatch by topic.

    Topics:
        lab_core             - ★ Sugahara Lab CORE: B-input Stop+Energy
                                (start here for any hysteresis FE work)
        catalog              - 13-model overview table + paper references
        decision_tree        - When to use which model (decision tree)
        jiles_atherton       - JA model formulation + parameters
        play                 - Play model (PRODUCTION in Radia)
        energy_based         - Energy-Based formulation theory (Bergqvist
                                / Henrotte / Francois-Lavet / Jacques lineage)
        efficient_eval       - Egger 2025 regularization + Schur complement
                                → O(K) inverse operator (computational
                                recipe behind MatEnergyHysteresis)
        preisach             - Preisach formalism (Everett function, FORC)
        vector               - Vector hysteresis (E&S, vector Preisach, RSST)
        input_dependent_shape - Matsuo trilogy: stop/play with B-dependent
                                shape function w(B)·g_0; required for
                                silicon steel (does not satisfy equal
                                vertical chords); 5 identification methods
        lamination_homogenization - Henrotte-Steentjes-Hameyer-Geuzaine 2015
                                two-step approach: 1D mesoscale + algebraic
                                H(B, Ḃ, p_k) for macroscale motor FE
        cellular_automaton   - Saito 兆古 lab (Hosei) — Preisach ≡ 2D CA;
                                Bitter-method domain visualization;
                                connection to Barkhausen noise / 1/f
        mathematical_foundations - Krasnoselskii-Pokrovskii 1989 +
                                Visintin 1994 (Applied Mathematical
                                Sciences vol 111): formal play/stop
                                variational inequalities, Hilpert's
                                inequality (PDE well-posedness),
                                Prandtl-Ishlinskii equivalence,
                                Preisach homogenization theorem
        rheological_models   - Krasnoselskii §39 + Visintin §III.4:
                                spring + friction mechanical analogy,
                                graph representation, transducers M/W,
                                origin of Bergqvist friction analogy
        lab_core_production  - ★★ Production deep dive: Henrotte 2006
                                (h = h_r + h_i, 4 params elementary,
                                5-fraction combined, coercivity staircase
                                identification) + Francois-Lavet 2013
                                (variational Ω = u - h·J + χ|J - J_p|,
                                atanh saturation curve, M250-50A concrete
                                parameters, Picard FE algorithm) + Egger
                                2025 inverse Schur. Use when implementing
                                the model.
        stop_vs_play_silicon_steel - Matsuo 2003 IEEE TMAG 39:1361:
                                silicon steel violates "equal vertical
                                chords" → pure stop fails on asymmetric
                                loops. Play from symmetric loops better
                                but imperfect. Combination model C_i
                                alternating play+stop corrections
                                converges to visually-perfect fit.
                                Lab production path = Energy-based
                                (handles dc-bias naturally) > combination
                                > input-dependent shape.
        congruency_selection - ★★ Sugahara Lab flagship logic (Hane 2026
                                accepted IEEE TMAG): measured H-axis
                                congruency / B-axis incongruency in NOES
                                35A360 FORCES B-input; energy-consistency
                                then forces STOP (B-input play shape
                                function is non-monotone -> no convex
                                energy).  => B-input energy STOP is the
                                unique (congruency-correct AND
                                thermodynamically consistent) model.
                                Includes the IGTE positioning vs Zeinali /
                                Zirka / Davino + the "do not say play is
                                impossible" scoping caveat.
        bergqvist_binput_stop - IGTE 2026 digest (Bergqvist energy-based
                                B-input Stop): the DOMAIN argument (Play
                                state domain [B-eta,B+eta] translates with
                                B -> no fixed coenergy -> Bergqvist
                                inapplicable; Stop [-eta,+eta] fixed ->
                                applies), the variational update replacing
                                the clip rule, EM self-consistent g_k
                                identification (BVLS I-spline), Potter
                                benchmark rms 13.7 A/m.  Sharpens
                                congruency_selection Step 2; notes the v4
                                (adaptive-eta) vs v5 (Bergqvist) digest split.
        peeling_identification - Analytic (剥ぎ取り) shape-function
                                identification of a B-input play/stop
                                model from descending FORC ladder via a
                                difference scheme (Everett-density mu),
                                NO least-squares.  Why per-hysteron LSQ
                                on the play variable is ill-conditioned
                                (p_k share B -> collinear).  Full lab
                                pipeline prototype_step0-4 + playmodel.py
                                shapeFunction.  Sugahara/Ahagon code.
        vector_play_model    - Canonical Sugahara/Ahagon B-input vector
                                play model: H=sum Hfunc_i(p_i), play
                                operator p_i=clip(p_i,B-zeta_i,B+zeta_i),
                                equidistant zeta, Hfunc by peeling.
                                Congruency (matches H-axis-congruent
                                silicon steel); vector/rotational BQM
                                variants; complementary to the energy
                                STOP (p=B-s).  VectorPlayModel/ code.
        all                  - Everything
    """
    topic = topic.lower().strip()
    if topic in ("lab_core", "core", "lab_core_method", "stop_energy",
                  "b_input_stop"):
        return LAB_CORE_METHOD
    if topic in ("lab_core_production", "lab_core_deep", "production",
                  "francois_lavet_2013", "henrotte_2006",
                  "variational_form", "n_fraction", "n_sphere",
                  "picard_fe_algorithm", "m250_50a", "implementation"):
        return LAB_CORE_PRODUCTION
    if topic in ("stop_vs_play_silicon_steel", "stop_vs_play",
                  "matsuo_2003", "silicon_steel_combination",
                  "combination_model", "equal_vertical_chords"):
        return STOP_VS_PLAY_SILICON_STEEL
    if topic in ("congruency_selection", "congruency", "congruency_property",
                  "hane_2026", "hane", "binput_selection",
                  "b_input_selection", "why_binput", "why_b_input",
                  "why_stop", "h_axis_congruency", "b_axis_incongruency",
                  "selection_argument", "35a360", "noes_congruency"):
        return CONGRUENCY_BINPUT_SELECTION
    if topic in ("bergqvist_binput_stop", "bergqvist_stop", "bergqvist_binput",
                  "igte2026", "igte_2026", "igte_digest", "stop_digest",
                  "variational_stop", "em_identification",
                  "dry_friction_pinning", "fixed_domain",
                  "play_domain_translates", "domain_argument"):
        return BERGQVIST_BINPUT_STOP
    if topic in ("peeling_identification", "peeling", "shape_function_id",
                  "shape_identification", "shapefunction", "ahagon",
                  "forc_difference", "everett_difference",
                  "play_calibration", "hagitori"):
        return PEELING_IDENTIFICATION
    if topic in ("vector_play_model", "vector_play", "canonical_play",
                  "playmodel", "play_hysteron", "playhysteron",
                  "ahagon_play", "binput_play", "b_input_play"):
        return VECTOR_PLAY_MODEL
    if topic in ("forward_inverse_variational", "forward_inverse",
                  "egger_2025", "co_energy", "energy_density",
                  "convex_duality", "team_32", "scalar_potential_hyst",
                  "vector_potential_hyst", "danskin"):
        return FORWARD_INVERSE_VARIATIONAL
    if topic in ("tp_eec_steady_state", "tp_eec", "kitao_2012",
                  "time_periodic", "steady_state_acceleration",
                  "ring_core_play", "inrush_current"):
        return TP_EEC_STEADY_STATE
    if topic in ("h_input_b_input_combination", "h_input_b_input",
                  "ferrite_equivalent_circuit", "ito_2013",
                  "newton_raphson_inversion", "nizn_ferrite"):
        return H_INPUT_B_INPUT_COMBINATION
    if topic in ("bobbio_1997_unification", "bobbio_1997",
                  "ptm_stm", "play_type_model", "stop_type_model",
                  "ptm_stm_inversion", "vitrovac"):
        return BOBBIO_1997_UNIFICATION
    if topic in ("asp_model", "asp", "asymmetric_probability",
                  "lee_2014", "transition_probability", "forc_model"):
        return ASP_MODEL
    if topic in ("zirka_hdhm", "zirka", "hdhm", "transplantation",
                  "congruency_based", "madelung", "history_dependent"):
        return ZIRKA_HDHM
    if topic in ("chua_model", "chua", "chua_bass", "chua_stromsmoe",
                  "dynamic_circuit", "spice_hysteresis", "ferroresonance"):
        return CHUA_MODEL
    if topic in ("other_models", "chan", "chan_transformer", "bouc_wen",
                  "boucwen", "potter_schmulian", "potter", "llg",
                  "micromagnetics", "recording_media", "niche_models"):
        return OTHER_MODELS
    if topic in ("jacques_monograph", "jacques_2018", "jacques",
                  "clausius_duhem", "three_approaches", "angle_searching",
                  "vector_play_approach", "variational_minimization",
                  "thermodynamic_derivation"):
        return JACQUES_MONOGRAPH
    if topic in ("catalog", "overview", "13_models"):
        return CATALOG_OVERVIEW
    if topic in ("decision_tree", "decision", "which", "tree"):
        return DECISION_TREE
    if topic in ("jiles_atherton", "ja", "jiles"):
        return JA_MODEL
    if topic in ("play", "play_model"):
        return PLAY_MODEL
    if topic in ("energy_based", "energy", "henrotte", "bergqvist",
                  "francois_lavet", "jacques"):
        return ENERGY_BASED
    if topic in ("efficient_eval", "efficient", "egger", "schur",
                  "regularization", "inverse_operator"):
        return EFFICIENT_EVAL
    if topic in ("preisach", "preisach_model", "everett"):
        return PREISACH
    if topic in ("vector", "vector_hysteresis", "es", "rotational"):
        return VECTOR_HYSTERESIS
    if topic in ("input_dependent_shape", "input_dependent",
                  "weighting_function", "matsuo", "silicon_steel_shape"):
        return INPUT_DEPENDENT_SHAPE
    if topic in ("lamination_homogenization", "homogenization", "lamination",
                  "two_step", "henrotte_2015", "steentjes"):
        return LAMINATION_HOMOGENIZATION
    if topic in ("cellular_automaton", "ca", "saito", "bitter",
                  "barkhausen"):
        return CELLULAR_AUTOMATON
    if topic in ("mathematical_foundations", "math", "visintin",
                  "krasnoselskii", "krasnosel_skii", "hilpert",
                  "variational_inequality", "prandtl_ishlinskii"):
        return MATHEMATICAL_FOUNDATIONS
    if topic in ("rheological_models", "rheological", "rheology",
                  "mechanical_analogy", "spring_friction", "friction"):
        return RHEOLOGICAL_MODELS
    if topic == "all":
        return "\n\n".join([
            LAB_CORE_METHOD, LAB_CORE_PRODUCTION, CATALOG_OVERVIEW,
            DECISION_TREE, JA_MODEL, PLAY_MODEL, ENERGY_BASED,
            EFFICIENT_EVAL, PREISACH, VECTOR_HYSTERESIS,
            INPUT_DEPENDENT_SHAPE, LAMINATION_HOMOGENIZATION,
            CELLULAR_AUTOMATON, MATHEMATICAL_FOUNDATIONS,
            RHEOLOGICAL_MODELS, STOP_VS_PLAY_SILICON_STEEL,
            FORWARD_INVERSE_VARIATIONAL, TP_EEC_STEADY_STATE,
            H_INPUT_B_INPUT_COMBINATION, BOBBIO_1997_UNIFICATION,
            ASP_MODEL, ZIRKA_HDHM, CHUA_MODEL, OTHER_MODELS,
            JACQUES_MONOGRAPH, CONGRUENCY_BINPUT_SELECTION,
            BERGQVIST_BINPUT_STOP, PEELING_IDENTIFICATION,
            VECTOR_PLAY_MODEL,
        ])
    return (f"Unknown topic '{topic}'. Available: lab_core, "
            "lab_core_production, catalog, decision_tree, "
            "jiles_atherton, play, energy_based, efficient_eval, "
            "preisach, vector, input_dependent_shape, "
            "lamination_homogenization, cellular_automaton, "
            "mathematical_foundations, rheological_models, "
            "stop_vs_play_silicon_steel, forward_inverse_variational, "
            "tp_eec_steady_state, h_input_b_input_combination, "
            "bobbio_1997_unification, asp_model, zirka_hdhm, "
            "chua_model, other_models, jacques_monograph, "
            "congruency_selection, bergqvist_binput_stop, "
            "peeling_identification, vector_play_model, all.")
