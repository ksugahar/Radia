"""
COMSOL practical tips harvested from Sugahara Lab (Kindai University)
NAS -- 5+ years (2020-2025) of COMSOL Multiphysics use.

This module is the radia-mcp surface for the same lab-discovered
COMSOL tips that ship in the COMSOL fork's markdown:
  S:/COMSOL/mcp-server/docs/LAB_TIPS.md

The two surfaces are intentionally kept in sync:
  * COMSOL fork (LAB_TIPS.md) -- ships in ksugahar/COMSOL_Multiphysics_MCP
    chromadb when --extra-pdfs S:/COMSOL/mcp-server/docs is enabled
  * radia-mcp (this file)     -- ships on PyPI as radia-mcp; available
    in any Sugahara Lab MCP client without requiring COMSOL fork install

Why this lives in radia-mcp.interop:
  * The tips are COMSOL-specific but lab-discovered
  * Sugahara Lab routinely cross-validates Radia (PEEC + MMM) / NGSolve
    (FEM) results against COMSOL solves -- the interop layer is the
    natural home for cross-tool numerical recipes
  * Lab users authoring with build123d / Cubit / Radia / NGSolve
    sometimes need to switch to COMSOL for a particular study
    (e.g. acoustic-structure coupling, electrochem); having the lab's
    tribal knowledge in the MCP catalog avoids re-discovering it

Source files (all under S:/COMSOL/, Japanese filenames given in romaji
+ translation; see LAB_TIPS.md for the verbatim Japanese names):
  * 2021_12_01_COMSOL_Solver_no_know_how.txt -- 12-tip seed list
    (Solver no Know-how)
  * 99_Tips/2023_10_20_SimmetryPlane_notVer.5.5.txt -- vendor reply
  * 99_Tips/Yousosuu_no_kakunin.txt -- DOF check shortcut
    (Element count check, ASCII filename: yososuu_no_kakunin)
  * 99_Tips/Comsol_Tips.docx -- coil current/voltage scaling
  * 99_Tips/Kasokuki_you_jishaku_no_settei.pptx -- accelerator magnet
    setup (kasokuki = accelerator)
  * 99_Tips/2024_04_04_Comsol_toiawase/ -- support inquiry artifacts
    (toiawase = inquiry)
  * 2023_01_29_Netsudendou_houteishiki_ni_kansuru_gimon.docx --
    heat-transfer BC bug (Netsudendou houteishiki = heat conduction
    equation)
  * 2024_03_08_SquarCoil_wo_mochiita_mesh_type_no_kentou/ -- mesh
    study (kentou = investigation)

All quoted Japanese is preserved as ASCII-safe romaji + English
translation in the docstring strings below (Windows cp932 safe).
The original Japanese (correct kanji, full filename, vendor quote
verbatim) is in the COMSOL fork's LAB_TIPS.md.

Authoring guide for additions:
  * Every tip MUST cite a concrete source file under S:/COMSOL/
  * Every solver/UI claim MUST be reproducible on COMSOL 6.1 or 6.2
    (the lab's deployed versions); flag version-specific behavior
  * Vendor inquiry quotes MUST cite the KESCO Japan / COMSOL support
    person and date (e.g. "Barbara, 2023-10-20")
  * ASCII-only in this file (Windows cp932 safe); put verbatim
    Japanese in the COMSOL fork's LAB_TIPS.md companion
"""

from __future__ import annotations


# ---------------------------------------------------------------------------
# Overview
# ---------------------------------------------------------------------------

LAB_TIPS_OVERVIEW = """
# Sugahara Lab COMSOL tips -- overview

This catalog distils Sugahara Lab's COMSOL know-how into ~20 topics
that an AI agent or new lab member can query for fast onboarding.

The single most important meta-tip:

  Sugahara Lab's primary EM stack is **Radia (MMM/MSC/PEEC) + NGSolve
  (FEM)**. COMSOL is used for:
    1. Validation cross-checks (third independent solver)
    2. Multi-physics couplings we don't yet support natively
       (acoustic-structure with PML, electrochem, viscoelastic)
    3. Reproducing reviewer results from papers that only published
       a `.mph`

  For pure magnetostatics / eddy current / accelerator-magnet work,
  the lab pipeline is **Radia + NGSolve** (faster, scriptable,
  hex-mesh-friendly). Use COMSOL when the physics genuinely requires
  it.

## Topic index

| Topic name                | What it covers |
|---------------------------|----------------|
| solver_choice             | MUMPS up to 1M DOF; GMRES for AC/DC > 1M |
| iteration_limit           | Fully Coupled 25 -> 250 default bump |
| acdc_iterative_solvers    | GMRES OK, BiCGStab/CG erratic in AC/DC |
| convergence_curve         | Reading the residual plot; lab heuristics |
| symmetry_plane_bc         | Magnetic Insulation vs Perfect Magnetic Conductor |
| version_55_unsupported    | Why old .mph files need re-validation |
| dof_count_check           | numberofdofs Help shortcut; quick sanity check |
| tokyo_tech_hpc            | FNL license + Batch Run for cluster parameter sweeps |
| coil_current_scaling      | Symmetry-plane factor-of-2 rules (CRITICAL!) |
| homogenized_multiturn_coil| Turns x current setup ; the (3) x (4) trap |
| coil_geometry_analysis    | When the Coil Geometry Analysis study step pays off |
| infinite_element_domain   | Accelerator-magnet open-boundary recipe |
| element_order             | 1/2/3-ji yousou: cost vs accuracy, lab default p=1 |
| mesh_statistics           | Free Tetrahedral Information 1/2 readout |
| stripping_mph_for_support | Save Solutions=OFF before emailing KESCO |
| heat_transfer_bc_dimcheck | Heat-Transfer-in-Solids BC dimensional bug |
| square_coil_mesh_study    | 2024-03-08 ObjectLayer vs AirLayer vs hex study |
| cubit_to_comsol_handoff   | Nastran .bdf pitfalls (compress, hanging nodes) |
| lab_file_index            | Where every quoted source lives on S:/COMSOL/ |

For the long-form bilingual version with full Japanese quotes, see
the COMSOL fork's LAB_TIPS.md.
"""


# ---------------------------------------------------------------------------
# Section 1: Solver choice -- MUMPS up to 1M DOF
# ---------------------------------------------------------------------------

LAB_TIPS_SOLVER_CHOICE = """
# Solver selection -- MUMPS up to 1M DOF

Lab seed note (S:/COMSOL/2021_12_01_COMSOL_Solver_no_know_how.txt):

  "Up to one million DOFs, the direct solver MUMPS is the right
  choice."  (Original Japanese romaji: "Hyaku-man jiyuu-do made wa,
  Direct Solver de MUMPS ga yoi." -- see LAB_TIPS.md for the
  verbatim kanji.)

## Rule of thumb

| Problem size (DOF) | Recommended           | Notes                              |
|--------------------|-----------------------|------------------------------------|
| < 100k             | Any direct            | Fast, exact                        |
| 100k -- 1M         | **MUMPS direct**      | Best balance memory/speed          |
| 1M -- 10M          | Iterative (GMRES)     | MUMPS RAM blows up                 |
| > 10M              | DD + GMRES            | Cluster only; needs FNL (see HPC) |

Default lab choice up to 1M DOF: MUMPS.  PARDISO when memory is
tight (PARDISO is more memory-efficient on Intel; MUMPS is slightly
more numerically robust on AC/DC complex matrices).

## Quick check before solving

In COMSOL: hit F1, search "numberofdofs", read.  Then aim for:
  100k -- 1M DOF -> click Solve, use MUMPS
  > 1M           -> switch to GMRES, plan for 10-100x more iter
  > 10M          -> stop, plan for HPC + FNL license

## See also

  ih_iteration_limit       -- iteration cap for nonlinear/coupled
  ih_acdc_iterative_solvers -- which iterative in AC/DC works
"""


# ---------------------------------------------------------------------------
# Section 2: Iteration limit for Fully Coupled
# ---------------------------------------------------------------------------

LAB_TIPS_ITERATION_LIMIT = """
# Iteration limit for Fully Coupled

Lab seed note:

  "For Fully Coupled, increasing Maximum #iteration from 25 to 250
  is often necessary."  (Original Japanese romaji: "Fully Coupled
  no Maxium#interation wo 25 -> 250 nado fuyasu hou ga yoi." -- see
  LAB_TIPS.md for kanji.)

## Why the default fails

For nonlinear / multi-physics problems (BH curve, hysteresis,
electrothermal coupling), COMSOL's default Fully Coupled iteration
cap of 25 is too low.  Symptoms:

  * "Maximum number of iterations reached" with residual stuck at
    1e-3 (not 1e-6)
  * Result looks plausible but is the partial-converged solution at
    iter 25
  * Re-running with maxiter = 100 gives a different (correct) answer

## Lab default

Bump Fully Coupled maxiter to **100--250** as the FIRST tuning knob,
BEFORE changing solver or preconditioner.

## Path in COMSOL GUI

  Study -> Solver Configuration -> Stationary/Time Dependent Solver
       -> Fully Coupled -> Maximum number of iterations

For Java automation: see
  Java/MATLAB recipes in COMSOL_JAVA_AUTOMATION.md (COMSOL fork)

## When 250 is still not enough

Switch to Newton with line search damping, or split the problem into
quasi-stationary steps.  See COMSOL_SOLVER_TIPS.md section 4
(convergence curve) before tuning anything else.
"""


# ---------------------------------------------------------------------------
# Section 3: AC/DC iterative solver gotchas
# ---------------------------------------------------------------------------

LAB_TIPS_ACDC_ITERATIVE_SOLVERS = """
# AC/DC iterative solver gotchas

Lab seed note:

  "In AC/DC, GMRES works but BiCGStab and CG behave suspiciously."
  (Original Japanese romaji: "AC/DC de wa, GMRES wa yoi ga,
  BiCGStab ya CG wa dousa ga ayashii." -- see LAB_TIPS.md for kanji.)

| Solver    | AC/DC reliability     | When to use                        |
|-----------|----------------------|------------------------------------|
| GMRES     | OK reliable           | Default iterative for AC/DC        |
| BiCGStab  | BAD erratic            | Avoid for AC/DC; OK other modules  |
| CG        | BAD only SPD systems   | Only for SPD (rare in AC/DC)       |

## Why

AC/DC Magnetic Fields produces non-symmetric matrices when:
  * Coil Geometry Analysis is active (non-symmetric current
    constraint)
  * Periodic BC with phase shift
  * Lumped port with non-Hermitian impedance

BiCGStab is theoretically applicable but breaks down silently on
non-symmetric matrices with bad conditioning.  CG requires SPD,
which AC/DC almost never delivers.

GMRES is non-symmetric-friendly and the lab default for AC/DC.

## See also

  COMSOL Reference Manual, chapter "Linear System Solvers"
  Sugahara Lab cross-check: same problem solved in NGSolve with
    radia.sparsesolv_ngsolve.COCRSolver (Complex Symmetric COCR --
    different math, used as third-party validator).
"""


# ---------------------------------------------------------------------------
# Section 4: Convergence-curve reading
# ---------------------------------------------------------------------------

LAB_TIPS_CONVERGENCE_CURVE = """
# Convergence-curve reading

Lab seed notes:

  "Iterative methods: look at the convergence curve. If it doesn't
  drop steadily, give up and change parameters."
  (Original Japanese romaji: "Hanpukuhou wa, shuusoku-kyokusen wo
  mite, ston-ston to ochinai baai ni wa, akiramete parameter no
  henkou." -- see LAB_TIPS.md for kanji.)

  "Memorize what a working convergence curve looks like."
  (Original Japanese romaji: "Umaku ugoiteiru toki no shuusoku-curve
  wo oboe-oku." -- see LAB_TIPS.md for kanji.)

## How to enable

  Study -> Solver Configuration -> Plot While Solving = ON
  -> Results While Solving plots residual every iter

## Curve diagnostics

| Curve shape              | Meaning                          | Action                              |
|--------------------------|----------------------------------|-------------------------------------|
| Straight line down (log) | Healthy                          | Wait                                |
| Steps / plateaus         | Stiff problem, solver re-evaluating | Wait longer                      |
| Oscillating              | Numerical instability            | Reduce step size, switch solver     |
| Flat (no progress)       | Hopelessly ill-conditioned       | Change preconditioner / re-mesh     |

## Lab heuristic -- the "shape memory" rule

For each problem class you study regularly (e.g. eddy current
heating, magnetostatic in iron, accelerator magnet sweep), remember
what the WORKING convergence curve looks like:

  * Roughly how many iterations to drop 4 orders of magnitude?
  * Linear-on-log or staircase-on-log?
  * Does the first iter drop a lot then plateau?

Anomalies stand out instantly.  If today's solve hits iter 80 but
your stored memory says "should converge by iter 30", something
changed (geometry, BC, material).  Stop and diagnose BEFORE waiting
to see if it eventually converges -- it usually won't.
"""


# ---------------------------------------------------------------------------
# Section 5: Symmetry Plane BC
# ---------------------------------------------------------------------------

LAB_TIPS_SYMMETRY_PLANE_BC = """
# Symmetry Plane BC -- Magnetic Insulation vs Perfect Magnetic Conductor

Vendor reply (Barbara at KESCO / COMSOL Japan, 2023-10-20):

  "The Symmetry Plane boundary condition is a feature that imposes
  either the Symmetry condition for magnetic flux density --
  equivalent to the Magnetic Insulation BC -- OR the Antisymmetry
  condition for B -- equivalent to the Perfect Magnetic Conductor
  BC."

(Verbatim Japanese is in the COMSOL fork's LAB_TIPS.md section 5;
omitted here for ASCII safety.)

## Decision rule

| Field vs plane                                          | BC choice                              |
|---------------------------------------------------------|----------------------------------------|
| B **parallel** to plane (tangential B preserved, normal B = 0) | **Magnetic Insulation** (Symmetry)     |
| B **perpendicular** to plane (normal B preserved, tangential B = 0) | **Perfect Magnetic Conductor** (Antisymmetry) |

## Vendor-cited references

  https://www.comsol.com/blogs/exploiting-symmetry-simplify-magnetic-field-modeling/
  https://doc.comsol.com/6.1/docserver/#!/com.comsol.help.acdc/acdc_ug_magnetic_fields.08.040.html

## CRITICAL companion -- current/voltage scaling

When you cut a coil with a symmetry plane, the input/output
quantities scale.  See topic = "coil_current_scaling".  Getting the
BC right but the scaling wrong gives an answer that is "almost
right" -- the failure mode for accidentally-published wrong results.

## Radia counterpart

For Radia / Sugahara Lab MMM/MSC pipeline, IMA (Image Method of
Analysis) is the analog:
  rad.Solve(container, ..., image='+x-z')  # parallel/perp signs
See Radia CLAUDE.md "IMA Sign Selection Policy".
"""


# ---------------------------------------------------------------------------
# Section 6: Version 5.5 unsupported
# ---------------------------------------------------------------------------

LAB_TIPS_VERSION_55_UNSUPPORTED = """
# Version 5.5 unsupported

Vendor reply (Barbara, 2023-10-20):

  "Version 5.5 is no longer supported; we recommend the latest
  COMSOL Multiphysics."

(Verbatim Japanese is in the COMSOL fork's LAB_TIPS.md section 6.)

## What broke between 5.5 -> 6.x

  * Symmetry Plane BC didn't exist in 5.5; manual coil BC was used
    (from S:/COMSOL/99_Tips/Kasokuki_you_jishaku_no_settei.pptx
    slide 1: lab note "5.5 does not have this setting; learn the
    coil specification method [in the old way] instead.")
  * Some BH-curve defaults changed (spline extrapolation behavior)
  * AC/DC Magnetic Fields default discretization moved Linear ->
    Quadratic (silent default change)

## Action when finding a 5.5 .mph

  1. Open via `File -> Open` -- COMSOL prompts to migrate
  2. Re-validate the migrated model against a known reference
     (analytical, or a linear-mu cross-check, or a Radia/NGSolve
     independent solve)
  3. Do NOT trust the migrated answer until validated -- the silent
     default changes can shift results by 10s of percent
"""


# ---------------------------------------------------------------------------
# Section 7: DOF count check
# ---------------------------------------------------------------------------

LAB_TIPS_DOF_COUNT_CHECK = """
# Quick DOF count check

Lab note (S:/COMSOL/99_Tips/Yousosuu_no_kakunin.txt, full content):

  https://kesco.co.jp/wp-content/uploads/2021/09/IntroductionToCOMSOLMultiphysics_53a_jp.pdf
  enter word
  numberofdofs

## How to use

In COMSOL:
  1. Hit F1 (Help)
  2. In the search box type `numberofdofs`
  3. Help opens the Reference Manual page on counting DOFs per
     physics

## Why it matters

ALWAYS look at #DOFs BEFORE clicking Solve.  Sweet spots:

  100k DOF  -- fast iteration loop (sub-1 min solves)
  1M DOF    -- MUMPS direct sweet spot (5-30 min solves)
  10M DOF   -- iterative solver (GMRES) territory, hours
  > 10M DOF -- cluster + FNL license needed

Also useful: `Mesh -> Statistics` shows element count + skewness
histogram.  See topic = "mesh_statistics".

## Lab convention for paper / poster reporting

Always quote in figures' captions:
  * #elements
  * #DOFs (after BC application)
  * min element quality (skewness > 0.1 OK; < 0.05 = red flag)
  * solver type (MUMPS / GMRES + preconditioner)
  * wall-clock time + RAM
"""


# ---------------------------------------------------------------------------
# Section 8: Tokyo Tech HPC + Batch Run
# ---------------------------------------------------------------------------

LAB_TIPS_TOKYO_TECH_HPC = """
# Tokyo Tech HPC + Batch Run (FNL license required)

Lab seed note:

  "Some Tokyo Tech faculty run COMSOL on the supercomputer. To use
  that license, FNL is required (extra cost). With Batch Run, 100
  cases run in parallel (proportional to node count?)."

(Verbatim Japanese in the COMSOL fork's LAB_TIPS.md section 8;
romaji summary: Toukoudai no senseigata de COMSOL wo SuperComputer
de ugokashite-iru hito ga iru. Sono license wo mochiiru tame ni
wa, FNL ga hitsuyou (tsuika hiyou). Batch Run no kinou wo mochiite,
100-case no jikkou ga heiretsu (node-suu ni hirei?) ni kanou.)

## Key facts

  * FNL (Floating Network License) is the **mandatory add-on** for
    HPC / cluster use.  Single-user or Class Kit licenses CANNOT
    be used on a cluster.
  * Confirm with KESCO (COMSOL Japan distributor) before planning
    Tokyo Tech / RIKEN / other shared HPC usage.
  * Cost: discuss with KESCO; FNL is significantly more expensive
    than single-user.

## Batch Run feature

  Study -> Job Sequence -> Batch
  -> Parametric sweep with one .mph per parameter case
  -> Each case becomes its own MPI / slurm job

Speedup is roughly linear in node count for embarrassingly parallel
sweeps (each case independent).  For tightly-coupled problems,
MUMPS / PARDISO stay node-local; only the parameter sweep
parallelizes across nodes.

## Lab alternative (free / open-source)

For large parametric sweeps where COMSOL is overkill, the Sugahara
Lab uses **mdx (Open Cluster)** + **Radia/NGSolve** running native
Python parallel sweeps.  See:
  - tools/push_pyds_to_mdx.py (LAB-private)
  - examples/*/bench_*.py JSON output protocol

This avoids the FNL cost for the "parameter sweep" use case.  Use
COMSOL Batch Run when you specifically need COMSOL's physics
modules (acoustic-structure, electrochem, etc.).
"""


# ---------------------------------------------------------------------------
# Section 9: Coil current scaling under symmetry
# ---------------------------------------------------------------------------

LAB_TIPS_COIL_CURRENT_SCALING = """
# Coil current scaling under symmetry (CRITICAL!)

From S:/COMSOL/99_Tips/Comsol_Tips.docx -- 1/8 symmetry rectangular
coil model in COMSOL's AC/DC Module documentation.

  "The Magnetic Insulation (magenta) and Perfect Magnetic Conductor
  (cyan) boundary conditions are applied along the appropriate
  symmetry planes for this problem."

## Scaling rules -- direction matters

| BC plane                                                  | Excitation    | Where to scale            | Factor       |
|-----------------------------------------------------------|---------------|---------------------------|--------------|
| **Magnetic Insulation** (current flows **normal** to plane) | **Voltage in**  | INPUT voltage divide by 2 | input / 2    |
| **Magnetic Insulation** (current flows **normal** to plane) | **Current in**  | OUTPUT voltage multiply by 2 | output * 2 |
| **Perfect Magnetic Conductor** (current flows **tangential** to plane) | **Current in**  | INPUT current divide by 2 | input / 2    |
| **Perfect Magnetic Conductor** (current flows **tangential** to plane) | **Voltage in**  | OUTPUT current multiply by 2 | output * 2 |

Apply ONCE per such symmetry plane.  A 1/8 model with three
symmetry planes can need a 2 * 2 * 2 = 8-fold scaling on whichever
side is post-processed.

## Trip wire (why this is in the lab tip catalog)

Getting one factor wrong silently gives an answer that is "almost
right" -- easy to publish, hard to retract.  COMSOL does NOT warn.
The result looks plausible because the field SHAPE is correct; only
the magnitude is off by 2^N.

## Lab checklist

Before publishing any half/quarter/eighth model result:

  1. Identify each symmetry plane and its BC type (Mag Insul vs PMC)
  2. Identify the excitation type (voltage in OR current in)
  3. Walk through the table above for each plane
  4. Compute the total scaling factor for the post-processed quantity
  5. Cross-check by solving the FULL geometry once (coarse mesh OK)
     and comparing post-processed quantities

## See also

  symmetry_plane_bc -- which BC for which field direction
  homogenized_multiturn_coil -- companion N*I trap
"""


# ---------------------------------------------------------------------------
# Section 10: Homogenized Multiturn Coil
# ---------------------------------------------------------------------------

LAB_TIPS_HOMOGENIZED_MULTITURN_COIL = """
# Homogenized Multiturn Coil -- turns * current factor

From S:/COMSOL/99_Tips/Kasokuki_you_jishaku_no_settei.pptx, slide 2
(romaji slide title: "Homogenized Multiturn Coil no settei no
shikata (2) Taisetsu!" -- "How to set up the Homogenized Multiturn
Coil (2), IMPORTANT!")

## The 4 input fields

The Homogenized Multiturn Coil node takes:

  (1) Coil current I_c [A]
  (2) Number of turns N
  (3) (cross-section variable)
  (4) (current-density derivation parameter)

The effective excitation is the PRODUCT N * I_c.  Fields (3) and
(4) also multiply into the effective NI.

  "(3) and (4) are multiplied into the setting value.
   In this case, set such that 2000/2 = 1000."
  (romaji of original: "(3) to (4) wa kakezan de setteichi.
   Kono toki wa 2000/2 no 1000 ni naru you ni settei shita."
   -- "(3) and (4) are multiplied into the setting value;
   in this case set so it becomes 1000 = 2000/2.")

## Worked example

For a 1/2 symmetry model of a 2000-A coil (real coil amps,
full-coil specification):

  Type into COMSOL: 1000 A  (= 2000 / 2)
  Because one symmetry plane cuts the coil in half

This is the same factor-of-2 from topic = "coil_current_scaling".
The Homogenized Multiturn node makes the trap easier to step on
because the input box just says "Coil current" with no annotation
that the value should be the symmetry-reduced amount.

## Lab convention

Always document in the .mph file Comments field:
  "Geometry: 1/2 symmetry on YZ plane (Mag Insul)
   Full-coil amps: 2000 A
   Input box `Coil current`: 1000 A"

So future-you (or the reviewer reproducing) can verify in 30 seconds.

## See also

  coil_current_scaling -- the master rule for symmetry-plane scaling
  coil_geometry_analysis -- when Homogenized fails and you need the
                            Coil Geometry Analysis study step
"""


# ---------------------------------------------------------------------------
# Section 11: Coil Geometry Analysis
# ---------------------------------------------------------------------------

LAB_TIPS_COIL_GEOMETRY_ANALYSIS = """
# Coil Geometry Analysis -- when to use the Geometry Analysis sub-node

From S:/COMSOL/99_Tips/Kasokuki_you_jishaku_no_settei.pptx, slide 2 (English
snippet quoted verbatim from AC/DC Module User's Guide):

  "The Coil Geometry Analysis study and study step are available
  for 3D models using the Magnetic Fields interface and the Coil
  node (which requires the AC/DC Module). [...] The best results
  are obtained when the coil has a smoothly varying or constant
  cross section without sharp bends and bottlenecks. The local
  current directions are solved for by the specialised Coil
  Geometry Analysis study step.

  The Geometry Analysis subnode to the Coil feature automatically
  appears to set up the automatic computation of the current flow
  in the coil. The boundary conditions for the Geometry Analysis
  are specified using the Input and Output subnodes available with
  the node."

## Lab rule of thumb

| Coil shape                                | Recommended            |
|-------------------------------------------|------------------------|
| Simple prismatic / toroidal block          | Homogenized Multiturn  |
| Single-turn round wire                     | Single-Conductor Coil  |
| Bends, splits, varying cross-section       | **Coil Geometry Analysis** |
| Saddle / racetrack / accelerator end-winding | **Coil Geometry Analysis** |
| Complex motor end-winding                  | **Coil Geometry Analysis** |

## Cost

Coil Geometry Analysis adds a **prepass study step** that solves a
small magnetostatic problem to determine local current direction.
For a typical accelerator-magnet model:

  Prepass: ~30 s to ~5 min (depends on coil DOF count)
  Main solve: unchanged

If the prepass diverges, the coil geometry is too convoluted; split
it into multiple coil nodes connected end-to-end.

## See also

  homogenized_multiturn_coil -- simpler when the coil is prismatic
  AC/DC Module User's Guide, "Coil Geometry Analysis" chapter
"""


# ---------------------------------------------------------------------------
# Section 12: Infinite Element Domain -- accelerator-magnet setup
# ---------------------------------------------------------------------------

LAB_TIPS_INFINITE_ELEMENT_DOMAIN = """
# Infinite Element Domain -- accelerator-magnet setup

From S:/COMSOL/99_Tips/Kasokuki_you_jishaku_no_settei.pptx, slide 3
(romaji slide title: "Infinite Element Domain no settei ni tsuite"
-- "About the Infinite Element Domain setup")

## Lab convention for accelerator magnet open-boundary modelling

| Layer       | Geometry                              | Material                              |
|-------------|---------------------------------------|---------------------------------------|
| Innermost   | Iron yoke + coil regions              | Iron / Cu / air as appropriate        |
| Middle      | Air sphere or box around the magnet   | Air (mu_r = 1)                        |
| Outermost   | **Infinite Element Domain** shell on the air boundary | Air (mu_r = 1); COMSOL stretches internally |

Slide labels (translated):

  * kuuki (bishou mesh fukumu) -- air (with refined mesh near the magnet)
  * kyuu no gaimen            -- outer surface of the air sphere
  * tetshin                   -- iron core
  * denpeki (men de shitei)   -- Electric Wall (specify by face)
  * jiheki (men de shitei)    -- Magnetic Wall (specify by face)
  * coil ryuunyuushutsu danmen hyouji-men -- coil in/out cross-section face

Lab annotations from the slide:

  * (romaji) "BH ni kanshite wa, ELF to no hikaku no tame you-kakunin,
    to iu ka, senkei de hikaku"
    -- For BH validation, cross-check vs ELF; in practice compare
    in the linear regime first.
  * (romaji) "Denkai-vector wa denki-kabe ni suichoku de, jikai-
    vector wa denki-kabe ni heikou"
    -- E vector is perpendicular to the electric wall; H vector is
    parallel to it.

## When to prefer Kelvin over Infinite Element Domain

| Method                     | Pros                  | Cons                                  |
|----------------------------|-----------------------|---------------------------------------|
| Dirichlet on truncation sphere | Simple             | Big mesh; truncation error 1/R         |
| Infinite Element Domain    | Built-in COMSOL       | Stretching hurts conditioning at HF    |
| PML                        | Standard for waves    | Heavy conditioning; tuning required    |
| **Kelvin transform**       | Exact open BC         | Setup tricky; needs Identity Pair      |

For magnetostatic / low-frequency eddy current, Kelvin transform is
the lab's preferred method.  See the COMSOL fork's
KELVIN_TRANSFORM.md (291 lines, 5 years of practice).

The Radia / NGSolve counterpart: Periodic BC + Kelvin reluctivity
scaling, see Radia CLAUDE.md "IMA" section + sparsesolv_ngsolve
ComplexCompactAMSPreconditioner.

## See also

  Kelvin transform tutorial (COMSOL fork): docs/KELVIN_TRANSFORM.md
  symmetry_plane_bc -- complementary BC for non-open boundaries
"""


# ---------------------------------------------------------------------------
# Section 13: Element-order setting (1-ji / 2-ji / 3-ji yousou)
# ---------------------------------------------------------------------------

LAB_TIPS_ELEMENT_ORDER = """
# Element-order setting (1jisei / 2jisei / 3jisei -- 1st / 2nd / 3rd)

From S:/COMSOL/99_Tips/Kasokuki_you_jishaku_no_settei.pptx, slide 7 (visual
comparison) and slide 15 (how to change the setting).

## How to change

  Component -> Magnetic Fields -> Discretization
            -> Magnetic vector potential = Linear / Quadratic / Cubic

Shortcut (slide 15, romaji: "2-ji yousou tou settei no shikata"
-- "How to set 2nd-order elements"):

  Component -> Mag (icon) -> here change 1 -> 2

## Cost vs accuracy

| Order        | Cost per element | Typical accuracy boost | When to use                                     |
|--------------|------------------|------------------------|-------------------------------------------------|
| 1 (Linear)   | 1x               | baseline               | Default; AC/DC eddy current; coarse-mesh studies|
| 2 (Quadratic)| ~5x              | ~10x                   | Mesh refinement too expensive; smooth fields    |
| 3 (Cubic)    | ~15x             | ~100x                  | Very smooth (electrostatic in air); rarely needed in AC/DC |

## Lab default

  Accelerator magnets: order 1 with a fine mesh.
  Order 2 when mesh refinement hits RAM limits.
  Order 3 almost never (cost grows super-linearly in 3D).

## Reference models cited in the slide

  https://www.comsol.com/model/vector-hysteresis-modeling-20671
  https://www.comsol.jp/model/electrodynamics-of-a-power-switch-33511

## Radia / NGSolve counterpart

  NGSolve: fes = HCurl(mesh, order=p)
  Radia + NGSolve curving: order 1-5 via export netgen
  Radia accelerator panel default: p=1 (mesh refinement preferred)

Same lab rule applies: p=1 with fine mesh first, p=2 only when
necessary, p=3+ is research territory.
"""


# ---------------------------------------------------------------------------
# Section 14: Mesh statistics
# ---------------------------------------------------------------------------

LAB_TIPS_MESH_STATISTICS = """
# Mesh statistics -- reading `Free Tetrahedral Information`

From S:/COMSOL/99_Tips/Kasokuki_you_jishaku_no_settei.pptx, slide 9
("Mesh1 > Free Tetrahedral 1 > Information 1 / Information 2") and
slide 16 ("how to check element count").

## How to find it

  Mesh -> Statistics
  (or right-click `Free Tetrahedral` -> `Statistics`)

The two `Information` reports show:

  Information 1: per-domain element count + skewness histogram +
                 quality measures (min / avg / max)
  Information 2: per-physics DOF count (per node / edge / face /
                 volume) and total DOFs after BC application

## Lab convention for paper / poster reporting

ALWAYS quote in figures / captions:

  * #elements (total + per-material breakdown)
  * #DOFs (after BC application)
  * min element quality (skewness)
    - > 0.1 = OK
    - 0.05 -- 0.1 = caution
    - < 0.05 = red flag (re-mesh or accept higher solver iter cost)
  * solver type + preconditioner
  * wall-clock time + peak RAM

## Lab cross-check protocol

For any non-trivial AC/DC study, also report:

  * Effective preconditioner cost (memory + time per iter), not
    just MUMPS factor time
  * Convergence-curve shape (cite topic = "convergence_curve")
  * Independent solver cross-check (Radia / NGSolve / FEMM)

## See also

  dof_count_check -- numberofdofs Help shortcut
  cubit_to_comsol_handoff -- pre-mesh in Cubit to get clean stats
"""


# ---------------------------------------------------------------------------
# Section 15: Stripping .mph for support
# ---------------------------------------------------------------------------

LAB_TIPS_STRIPPING_MPH_FOR_SUPPORT = """
# Stripping results before sending .mph to KESCO support

From S:/COMSOL/99_Tips/Kasokuki_you_jishaku_no_settei.pptx, slide 15
(romaji slide title: "Comsol-sha ni shitsumon tou, kekka wo
torinozoite souzuku shitai toki" -- "When sending a question /
inquiry to COMSOL and you want to strip out the results first")

## Procedure when emailing a .mph to KESCO / COMSOL Forum

  1. File -> Save As with `Save Solutions = OFF`
  2. Confirm saved file is < 50 MB (typical attachment limit)
  3. Include a one-paragraph reproduction recipe in the email body:
     "Open file -> click Compute on Study 1 -> error X in step Y"

Without stripping the solution dataset, a typical AC/DC model is
500 MB -- 5 GB (mostly raw GridFunction storage) and bounces from
the support inbox.

## Actual lab support inquiry (2024-04-04)

  S:/COMSOL/99_Tips/2024_04_04_Comsol_toiawase/
    BH16.txt            (660 B, 50-point BH curve, hand-built)
    M3GB_to_COMSOL.csv  (50.5 MB, raw BH from MAT-3GB measurement)
    model_forQA.mph     (175 MB, stripped model)

The lab sent the 175 MB stripped model to KESCO; the unstripped
version would have been ~3 GB (typical for a meshed iron core with
hysteresis state stored).

## BH16.txt contents (representative steel BH curve)

50 points covering H = 0 to H = 6.2e6 A/m, plus 2 extrapolation
points at H = 6.4e5 A/m (B = 3.01 T) and H = 6.2e6 A/m (B = 10.57 T)
to keep COMSOL's spline extrapolator linear past saturation.

Example excerpt:
  0          0
  10         0.0156
  16         0.02925
  ...
  100000     2.246281
  640599     3.01226       <- linear extrapolation kicks in
  6207043    10.57212      <- mu_0 H + saturation

The lesson: NEVER let COMSOL extrapolate a BH curve in log-mode --
it diverges exponentially.  Always include 2 linear-regime points
past the highest expected H.

## See also

  Radia: rad.MatSatIsoTab(BH_DATA) accepts [[H, B], ...] in A/m, T
  Same lesson: include linear-regime extrapolation points
"""


# ---------------------------------------------------------------------------
# Section 16: Heat-transfer BC dimensional check
# ---------------------------------------------------------------------------

LAB_TIPS_HEAT_TRANSFER_BC_DIMCHECK = """
# Heat-transfer boundary equation -- dimensional consistency

From S:/COMSOL/2023_01_29_Netsudendou_houteishiki_ni_kansuru_gimon.docx (lab Q&A).

The lab found that COMSOL's `Heat Transfer in Solids (ht)` module
documentation displayed a boundary equation whose second term FAILS
dimensional analysis:

  "Applying dimensional analysis to equation (1), only the second
  term is wrong. If instead we use equation (4), dimensional
  analysis confirms consistency."

Original:
  (romaji of original Japanese: "Shiki (1) ni taishite jigen-kaiseki
  wo shitemiru to, dai-2-kou nomi okashii koto ga wakarimasu.
  Shiki (1) no kawari ni, shiki (4) to shite jigen-kaiseki wo
  okonau to, mumujun de aru koto ga kakunin dekimasu."
  -- "Performing dimensional analysis on equation (1), only the
  2nd term is wrong. Using equation (4) instead, dimensional
  analysis confirms consistency.")

## Lab rule when implementing custom heat-transfer BCs in COMSOL

  1. Always run dimensional analysis on the equation symbol-by-symbol
     before believing the COMSOL documentation
  2. The General Form PDE node lets you check units interactively
     (`Units` column turns red if inconsistent)
  3. If in doubt, compare against the Egger / Henrotte axisymmetric
     conventions in Radia (CLAUDE.md "Axisymmetric FE" section) --
     those are derived by hand and unit-checked

## Companion file

  S:/COMSOL/2023_01_29_Netsudendou_mondai/Netsudendou_hou_houteishiki.docx -- the lab's
  by-hand derivation: 2D + axisymmetric + variable-separation form
  of dT/dt = alpha * Laplacian(T) (radiation neglected).  This is
  the preferred reference for cross-checking COMSOL's `ht` output.

## Radia / NGSolve counterpart

  Axisymmetric heat: standard NGSolve H1 + 2*pi*r weighting
  (FEMM-canonical convention, see Radia CLAUDE.md "Axisymmetric FE:
  Henrotte for Magnetic, Standard H1 for Scalar").

  The convention split is documented in
  memory/reference_femm_source_axisym_conventions.md and verified
  against FEMM 4.2 source (S:/FEMM/02_source/femm42src_22Oct2023/
  hsolv/prob1big.cpp).
"""


# ---------------------------------------------------------------------------
# Section 17: Square-coil mesh-type investigation (2024-03-08)
# ---------------------------------------------------------------------------

LAB_TIPS_SQUARE_COIL_MESH_STUDY = """
# Square-coil mesh-type investigation (2024-03-08)

From S:/COMSOL/2024_03_08_SquareCoil_wo_mochiita_mesh_type_no_kentou/

Lab study comparing 8 mesh strategies for a square coil + air model.

## Variants tested

| Variant                              | Geometry interface                  | Mesh scheme | Density |
|--------------------------------------|-------------------------------------|-------------|---------|
| 01_AirLayer/coarse, fine             | Air block surrounds coil, **air-layer interface** | Free tet | coarse / fine |
| 02_ObjectLayer/coarse, fine          | Air + coil share **object-layer interface** (subtract A from B keep) | Free tet | coarse / fine |
| 03_tet                               | Full tet, single material            | tetmesh     | medium  |
| 04_hex/roundHex-corse, fine          | Hex sweep with rounded coil corners  | hex         | coarse / fine |
| 04_hex/straightHex                   | Hex sweep with straight (chamfered) corners | hex  | medium  |

## Key findings (from the diff PNGs, see jou folder)

| Scheme                              | Diff from all-tet ref | Recommendation                  |
|-------------------------------------|----------------------|----------------------------------|
| AirLayer + free tet                 | ~1% diff in Bz       | Good general default             |
| ObjectLayer + free tet              | Cleanly imprints; matches ref | **Best for sharp boundary fields** |
| Hex sweep, round corners            | Smooth Bz, matches ref at fine density | Best for cylinder-like coils |
| Hex sweep, straight corners         | ~3% bump at corner-pyramid transition | **Avoid for accurate B near corners** |

## Practical rule (lab take-away)

For COMSOL:

  * Prefer **ObjectLayer + free tet** unless the coil cross-section
    is round (then hex sweep with rounded corners wins).
  * Avoid hex sweep with **straight** (chamfered) corners -- the
    pyramid transition causes a visible Bz step.
  * The all-tet reference (`for_refarence_all-tet.jou`) is the
    fastest path to a correct baseline; only optimize mesh type
    once you have it.

## Reference Cubit recipe

S:/COMSOL/2024_03_08_SquareCoil_wo_mochiita_mesh_type_no_kentou/jou_file/for_refarence_all-tet.jou

```cubit
reset
brick x 1   y 1   z 0.2   # air block
brick x 0.6 y 0.6 z 0.2   # coil block
brick x 0.2 y 0.1 z 0.2   # coil inner
move Volume 3 x 0.4 y 0.05 include_merged
brick x 3   y 3   z 2.6   # outer air

subtract volume 2 from volume 1 keep_tool
subtract volume 3 from volume 1 keep_tool
subtract volume 1 2 3 from volume 4 keep_tool

imprint all
merge all
compress

volume all scheme tetmesh
volume 2 3 1 size auto factor 3
mesh   volume 2 3 1
volume 4   size auto factor 4
mesh   volume 4

block 1 add volume 1 ; block 1 name "V_coil1"
block 2 add volume 3 ; block 2 name "V_coil2"
block 3 add volume 2 4 ; block 3 name "V_air"

nodeset 1 add surface 15 16 17 18 19 20 ; nodeset 1 name "S_outer"
nodeset 2 add surface 13                ; nodeset 2 name "S_in"
nodeset 3 add surface 11                ; nodeset 3 name "S_out"

export nastran "O:/SquareCoil_tet.bdf" dimension 3 overwrite everything large
# Number of elements = 434229
```

## Investigation artifacts

  03_tet/SquareCoil_tet.mph     (241 MB, COMSOL all-tet baseline)
  03_tet/SquareCoil_tet.bdf     (32.5 MB, Nastran export)
  Comsol_kaiseki_hikaku_jissekihyou.xlsx    (Excel comparison ledger)
  difference_plot_comsol_*.py (x4)     (Python plot scripts -- generate PNGs)

## See also

  cubit_to_comsol_handoff -- Nastran .bdf gotchas from this study
  Radia CLAUDE.md "AI-Driven Cubit: Probe, Don't Guess" -- the
    after-subtract-keep duplicate-volume warning rediscovered here
"""


# ---------------------------------------------------------------------------
# Section 18: Cubit-to-COMSOL handoff
# ---------------------------------------------------------------------------

LAB_TIPS_CUBIT_TO_COMSOL_HANDOFF = """
# Cubit-to-COMSOL handoff (Nastran .bdf)

Full pipeline is documented in the COMSOL fork's CUBIT_TO_COMSOL.md.
This topic only collects lab-specific gotchas discovered while
debugging the section-17 square-coil investigation.

## Gotchas

| Pitfall                                          | How it appears in COMSOL                                | Workaround                                                              |
|--------------------------------------------------|---------------------------------------------------------|-------------------------------------------------------------------------|
| Sideset name lost on Nastran round-trip          | COMSOL shows sideset as `Boundary System Group 1`       | Keep a manifest: block 1 = V_coil, sideset 1 = S_in, etc.               |
| Block ID renumbered after `compress`             | COMSOL imports `Material 1` not equal to original block 1 | Run `compress` BEFORE assigning blocks, or skip `compress`              |
| Hanging nodes in hex+tet boundary                | COMSOL shows `Invalid mesh: node X not connected`        | Use ObjectLayer scheme (topic 17); or imprint the hex-tet interface     |
| Pyramid transition wedges silently rejected       | Mesh import OK, but Bz has 5% step at the wedge layer    | Pre-mesh the pyramid layer in Cubit with `volume X scheme tetmesh`     |
| `# Number of elements = N` Cubit `.jou` comment   | Cubit ignores `#` in Aprepro; harmless                   | Use as a soft contract: re-running should hit +/- 1% of recorded count  |
| `subtract A from B keep` leaves duplicate of B    | Both COMSOL-imported volumes have correct geometry; cube list has 1 extra entry | Delete duplicate explicitly via centroid + volume inspection (Radia "AI-Driven Cubit" rule) |

## Radia / NGSolve counterpart

For the Radia pipeline:

  Cubit -> export netgen "f.vol" order N -> NGSolve Mesh()

This bypasses Nastran .bdf entirely and:
  * Preserves block/sideset names as material/boundary labels
  * Supports curved high-order elements (order 1-5)
  * Includes companion .vol.json with CAD volumes/areas for
    consistency checking (`check-vol model.vol`)

When you need COMSOL specifically (e.g. acoustic-structure with
PML), the .bdf pipeline is the bridge; the gotchas above are the
hard-won knowledge from the lab using it in production.

## See also

  COMSOL fork docs/CUBIT_TO_COMSOL.md -- full pipeline (motor / WPT examples)
  Radia CLAUDE.md "Cubit Block/Sideset Label Convention"
"""


# ---------------------------------------------------------------------------
# Section 19: Lab file index
# ---------------------------------------------------------------------------

LAB_TIPS_LAB_FILE_INDEX = """
# Lab file index (sources of these tips on S:/COMSOL/)

This index lists every file under S:/COMSOL/ that a tip in this
catalog cites by name.  All paths are NAS-internal (Sugahara Lab
file server).  External users see the COMSOL fork's LAB_TIPS.md for
the same list with web-friendly framing.

## Plain text seed notes

  S:/COMSOL/2021_12_01_COMSOL_Solver_no_know_how.txt       (766 B)
    -- the original 12-line seed list (topics 1, 2, 3, 4, 8)

  S:/COMSOL/99_Tips/2023_10_20_SimmetryPlane_notVer.5.5.txt  (864 B)
    -- KESCO Japan support reply re: Symmetry Plane BC + 5.5
       sunset (topics 5, 6)

  S:/COMSOL/99_Tips/Yousosuu_no_kakunin.txt                      (116 B)
    -- the 3-line numberofdofs Help shortcut (topic 7)

  S:/COMSOL/99_Tips/2024_04_04_Comsol_toiawase/BH16.txt    (660 B)
    -- the 50-point representative steel BH curve used in the
       2024-04-04 KESCO inquiry (topic 15)

## Word / PowerPoint sources

  S:/COMSOL/99_Tips/Comsol_Tips.docx                     (14.8 KB)
    -- coil current/voltage scaling under symmetry (topic 9)

  S:/COMSOL/99_Tips/Kasokuki_you_jishaku_no_settei.pptx              (1.3 MB, 18 slides)
    -- accelerator magnet setup tips (topics 10, 11, 12, 13, 14, 15)

  S:/COMSOL/2023_01_29_Netsudendou_houteishiki_ni_kansuru_gimon.docx     (31.9 KB)
    -- heat-transfer BC dimensional bug write-up (topic 16)

  S:/COMSOL/2023_01_29_Netsudendou_mondai/Netsudendou_hou_houteishiki.docx    (33.8 KB)
    -- companion: lab's by-hand axisymmetric heat derivation (topic 16)

## Folder-level artifacts

  S:/COMSOL/99_Tips/2024_04_04_Comsol_toiawase/             (~225 MB)
    -- KESCO support inquiry: model_forQA.mph (175 MB),
       M3GB_to_COMSOL.csv (50 MB) -- topic 15

  S:/COMSOL/2024_03_08_SquareCoil_wo_mochiita_mesh_type_no_kentou/  (~290 MB)
    -- 8 mesh variants + 4 diff PNGs + 1 xlsx ledger -- topic 17
       jou_file/for_refarence_all-tet.jou (1.1 KB) -- reference
       03_tet/SquareCoil_tet.mph (241 MB) -- baseline COMSOL model
       03_tet/SquareCoil_tet.bdf (32.5 MB) -- Nastran export
       Comsol_kaiseki_hikaku_jissekihyou.xlsx (114 KB) -- comparison ledger

  S:/COMSOL/88_Kelvin_henkan/                              (large)
    -- 5-year Kelvin transform sub-project history
       (2020-06 to 2023-03; ~89 GB across 6 subfolders).
       Documented in the COMSOL fork's KELVIN_TRANSFORM.md.

  S:/COMSOL/2020_07_07_Clawpole/                          (medium)
    -- claw-pole motor mesh study; lab's reference for the
       Cubit-to-COMSOL Nastran .bdf pipeline.
       Documented in the COMSOL fork's CUBIT_TO_COMSOL.md.

## How to read this catalog

If you are an AI / new lab member trying to USE the tip:

  1. Call `interop_comsol_lab_tips(topic="overview")` to see the
     topic index
  2. Call `interop_comsol_lab_tips(topic="<topic_name>")` for any
     specific topic
  3. The topic body is self-contained -- no need to open the source
     file for normal use

If you are debugging or extending a tip:

  1. Find the source file in this index
  2. Open it on the Sugahara Lab NAS (S:\\COMSOL\\)
  3. If you find new info, add a topic here AND update the COMSOL
     fork's LAB_TIPS.md in lockstep
"""


# ---------------------------------------------------------------------------
# Topic registry
# ---------------------------------------------------------------------------

TOPICS = {
    "all":                          None,
    "overview":                     LAB_TIPS_OVERVIEW,
    "solver_choice":                LAB_TIPS_SOLVER_CHOICE,
    "iteration_limit":              LAB_TIPS_ITERATION_LIMIT,
    "acdc_iterative_solvers":       LAB_TIPS_ACDC_ITERATIVE_SOLVERS,
    "convergence_curve":            LAB_TIPS_CONVERGENCE_CURVE,
    "symmetry_plane_bc":            LAB_TIPS_SYMMETRY_PLANE_BC,
    "version_55_unsupported":       LAB_TIPS_VERSION_55_UNSUPPORTED,
    "dof_count_check":              LAB_TIPS_DOF_COUNT_CHECK,
    "tokyo_tech_hpc":               LAB_TIPS_TOKYO_TECH_HPC,
    "coil_current_scaling":         LAB_TIPS_COIL_CURRENT_SCALING,
    "homogenized_multiturn_coil":   LAB_TIPS_HOMOGENIZED_MULTITURN_COIL,
    "coil_geometry_analysis":       LAB_TIPS_COIL_GEOMETRY_ANALYSIS,
    "infinite_element_domain":      LAB_TIPS_INFINITE_ELEMENT_DOMAIN,
    "element_order":                LAB_TIPS_ELEMENT_ORDER,
    "mesh_statistics":              LAB_TIPS_MESH_STATISTICS,
    "stripping_mph_for_support":    LAB_TIPS_STRIPPING_MPH_FOR_SUPPORT,
    "heat_transfer_bc_dimcheck":    LAB_TIPS_HEAT_TRANSFER_BC_DIMCHECK,
    "square_coil_mesh_study":       LAB_TIPS_SQUARE_COIL_MESH_STUDY,
    "cubit_to_comsol_handoff":      LAB_TIPS_CUBIT_TO_COMSOL_HANDOFF,
    "lab_file_index":               LAB_TIPS_LAB_FILE_INDEX,
}


def get_comsol_lab_tips_documentation(topic: str = "all") -> str:
    """Return Sugahara Lab COMSOL practical tips for the requested topic.

    Args:
        topic: One of the keys in TOPICS, or "all" for everything.

    Returns:
        Markdown-formatted documentation string.  If topic is unknown,
        returns a help message listing available topics.
    """
    topic = (topic or "all").strip().lower()
    if topic == "all":
        return "\n\n".join(
            f"## {key}\n{value}"
            for key, value in TOPICS.items()
            if key != "all" and value
        )
    if topic in TOPICS:
        value = TOPICS[topic]
        if value is None:
            # Defensive: only "all" should map to None
            return "(topic body not loaded)"
        return value
    available = ", ".join(k for k in TOPICS if k != "all")
    return (
        f"Unknown topic: {topic!r}.\n\n"
        f"Available topics: {available}.\n"
        f"Pass topic='all' for the full Sugahara Lab COMSOL tips catalog,\n"
        f"or topic='overview' for the topic index with one-line summaries."
    )
