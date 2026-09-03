"""BEM/MoM overview, decision tree, and lab stack."""

LAB_STACK = r"""
# Sugahara Lab BEM/MoM stack (production, 2026-07)

The lab uses BEM/MoM techniques in MULTIPLE places, each with a different
formulation and tool:

| Use case | Formulation | Tool | Files (this folder) |
|----------|-------------|------|---------------------|
| **Permanent magnets / soft iron** ★ | HDiv-VIM / charge Gram | radia HDiv-VIM + NGSolve | HDiv-VIM / reduced FEM |
| **IH workpiece + Kelvin** | FEM-SIBC + Robin BC (NOT pure BEM, but BEM-like) | radia `calc_fem_kelvin.py` | 03_bem_eddy_current/, 05_low_freq_* |
| **PEEC filament+panel** | FastImp-style PEEC | radia `peec_matrices` | 80_applications/Fast_Impedance_* |
| **PEEC conductor analysis plus HDiv-VIM / reduced-FEM magnetic-material coupling** | PEEC + HDiv-VIM / reduced FEM | radia + NGSolve | application notebooks |
| **HCurl Galerkin BEM** | EFIE/MFIE/PMCHWT with RWG | `ngsolve.bem` | 04_efie_mfie_cfie/, 05_low_freq_* |
| **H-matrix acceleration** | ACA on charge-Gram / BEM blocks | radia `HACApK` | 06_h_matrix_aca/ |
| **FEM-BEM hybrid (transformer)** | FEM inside, BEM outside | (not lab production, ref only) | 08_fem_bem_hybrid/ |

★ = Radia's current magnetic-material route.  Follow NGSolve terminology and
prefer HDiv-VIM / reduced-FEM coupling over proprietary legacy vocabularies.

## What the lab does NOT use

- **FMM** — Removed from Radia 2026-03 (see CLAUDE.md).  HACApK with ACA
  handles the same large-N problems with better practical performance.
- **Calderón preconditioners** — Only relevant for ngsolve.bem
  high-frequency scattering (not a lab focus).
- **Null-field / TDS** — Reference only; lab does not use thin-dielectric
  approximations.
- **Wire-grid models** — NEC-era; superseded by RWG triangular surface
  meshes everywhere.

## Why HDiv-VIM is the Radia magnetic-material route

HDiv-VIM keeps the magnetic-material unknowns in the same finite-element
language used by NGSolve.  This is the important advantage over old
collocation-style integral routes:

| Property | HDiv-VIM / reduced FEM | Galerkin RWG-BEM (ngsolve.bem) |
|----------|-----------------|--------------------------------|
| Basis order | NGSolve HDiv order, curved elements | RWG (linear edge basis) |
| Element type | FEM volume mesh, reduced-FEM compatible | Triangles on surface only |
| Material | Nonlinear iron, PM, reduced-FEM coupling | Linear conductors |
| Open boundary | Kelvin / image / charge-Gram route | Natural (surface integral) |
| Kernel | Laplace 1/r only | Laplace OR Helmholtz |
| Acceleration | HACApK charge Gram | H-matrix or FMM |
| Lab production | YES ★ | Optional (for high-freq scattering) |

The choice is deliberate: Radia should align with NGSolve's abstractions,
especially where reduced FEM and higher-order curved meshes matter.  Galerkin
BEM remains appropriate for high-frequency scattering, where the lab uses
NGSolve's `bem` module instead.
"""


DECISION_TREE = r"""
# Decision tree: which BEM/MoM technique for which problem

```
1. What's the SOURCE?
   ├── Permanent magnet (fixed M)
   │   → HDiv-VIM / Radia field evaluator
   │     → radia_mcp.radia_ngsolve.hdiv_vim
   │
   ├── Soft iron (M unknown, BH curve)
   │   → HDiv-VIM + reduced-FEM-compatible nonlinear solve
   │     → radia_mcp.electromagnet / radia_mcp.radia_ngsolve
   │
   ├── Coil current (known I)
   │   → Biot-Savart analytical or PEEC filaments
   │     → CoilBuilder or PEECBuilder
   │
   └── Eddy current source (frequency-domain σ, ω)
       └── See "What's the FREQUENCY?"

2. What's the FREQUENCY?
   ├── DC / static
   │   → HDiv-VIM / reduced FEM
   │
   ├── Quasi-static (MQS, DC to 100 kHz)
   │   ├── Thin conductor with skin effect
   │   │   → PEEC + SIBC (radia.peec_matrices)
   │   │     → HDiv-VIM / reduced-FEM coupling
   │   ├── Bulk conductor (workpiece in IH)
   │   │   → FEM-SIBC + Kelvin in radia
   │   │     → calc_fem_kelvin.py --vol model.vol
   │   └── Need surface BEM only (no FEM volume)
   │       → ngsolve.bem Weggler single-trace (low-freq stable)
   │
   ├── 1 MHz - 100 MHz (transition, Darwin approximation)
   │   → PEEC with full radiation
   │     → or ngsolve.bem Calderón-preconditioned EFIE
   │
   └── > 100 MHz (full wave, Helmholtz)
       → ngsolve.bem (EFIE/MFIE/CFIE with FMM or H-matrix)
       → OUT OF Radia/lab scope per CLAUDE.md Laplace-kernel policy

3. Which OPERATOR owns the solve?
   ├── HDiv-VIM magnetic material
   │   → use vim.Solve / vim.HDivSolver.Solve named solver and Gram options
   │   → HACApK compresses charge-Gram operators where that route selects it
   ├── PEEC or ngsolve.bem
   │   → use that formulation's named solver/preconditioner API
   └── Legacy Radia C++ relaxation object
       → method=0 dense LU only; methods 1 and 2 are retired

4. What's the GEOMETRY topology?
   ├── Simply connected
   │   → Standard MoM/BEM (any formulation)
   ├── Multiply connected (loops in conductor)
   │   ├── A-V Bíró-Preis gauge (FEM) or PEEC loop currents
   │   └── Tree-cotree decomposition for ngsolve.bem
   └── Open boundary (unbounded domain)
       ├── HDiv-VIM: Kelvin / image / charge-Gram route
       ├── FEM: Kelvin transformation + Radia evaluator
       └── BEM (ngsolve.bem): natural via radiation BC
```

## Common pitfalls

| Pitfall | Symptom | Fix |
|---------|---------|-----|
| Using EFIE at low freq | nan, ill-conditioned | Switch to Loop-Star or Weggler ST |
| MFIE with RWG below 1 MHz | wrong solution | Use Vico-Greengard MFIE-LF stabilization |
| Non-uniform M expected | wrong low-order approximation | Use higher-order HDiv / curved elements |
| Treating HACApK as `rad.Solve(method=2)` | unsupported legacy selector | Use the owning HDiv/PEEC/BEM named API |
| Galerkin BEM on nonlinear iron | awkward nonlinear coupling | Use HDiv-VIM / reduced FEM |
| FMM on Radia | Removed 2026-03 | Use HACApK ACA instead |
"""


HISTORY = r"""
# Genealogy: BEM/MoM lineage 1968 -> 2024

```
1939  Stratton-Chu  — surface integral representation
1968  Harrington — "Field Computation by Moment Methods" (THE textbook)
1969  Müller — Foundations of EM wave theory (CFIE-like precursor)
1980s NEC code (Burke-Poggio) — wire-grid MoM
1982  Rao-Wilton-Glisson RWG basis (IEEE TAP 30:409) ★
      — triangular surface basis, foundation of all surface BEM
1987  Greengard-Rokhlin FMM (J Comp Phys 73:325) ★
      — N^2 -> N log N for N-body / BEM
1988  Bossavit edge elements (also key for FEM-BEM hybrid)
1989  Bíró-Preis A-V gauge (IEEE IntMag) — eddy current BEM
1993  Saad iterative methods (textbook)
1997  Bebendorf-Tyrtyshnikov ACA — early matrix compression
2000  Bebendorf ACA original (Num Math 86:565) ★
      — practical low-rank approximation of BEM matrices
2003  Bebendorf-Rjasanow adaptive low-rank
2008  Andriulli Calderón preconditioner EFIE (IEEE TAP)
2008  Hsiao-Wendland Boundary Integral Equations textbook
2011  Sauter-Schwab Boundary Element Methods textbook
2012  Vico et al. MFIE at very low freq (Vico-Greengard family)
2013  Vico-Greengard-Gimbutas-Cools MFIE-LF (IEEE TAP 61:1285)
2013  Weggler single-trace formulation low-freq
2021  Ostrowski-Hiptmair Two-Step Maxwell freq-stable (SIAM SISC 43)
2007-25 HACApK library (CREST Post-Peta) — distributed H-matrix
2024  Lattice H-matrix, BLR-on-tensor-cores (modern hardware)
```

## Sugahara lab specific lineage

```
- ELF_MAGIC (1980s-2000s, lab predecessor of Radia)
  └── 02_mmm_surface_charge/lab_predecessor_ELF_MAGIC/
- Ishibashi group (Tohoku) — eddy current BEM with SIBC
  └── 03_bem_eddy_current/ishibashi/
- Hitachi/Koizumi — VIE with edge elements for nonlinear+eddy
  └── 02_mmm_surface_charge/hitachi_koizumi/
- Radia (Chubar-Elleaume-Chavanne 1998) ★
  └── 02_mmm_surface_charge/radia/
- HACApK (Ida-Iwashita) — Post-Peta CREST project
  └── 06_h_matrix_aca/HACApK_*
```

## Cross-references

- `radia_mcp.radia_ngsolve` — radia code usage, ngsolve.bem integration
- `radia_mcp.peec` — PEEC filament/panel + Carstensen + HOIBC
- `radia_mcp.ih` — IH-specific SIBC + ESIM
- `radia_mcp.matrix_solvers` — solver layer (CG/COCR/HACApK)
- `radia_mcp.motor` — Hollaus effective material (uses BEM-like cell problem)
"""


def get_overview_knowledge(topic: str = "decision_tree") -> str:
    """Dispatch BEM/MoM overview topics.

    Topics:
        lab_stack       - Production BEM/MoM components + which to pick
        decision_tree   - Solver+formulation decision tree (DEFAULT)
        history         - Genealogy 1968 Harrington -> 2024
        all             - Everything
    """
    topic = topic.lower().strip()
    if topic in ("lab_stack", "stack", "production"):
        return LAB_STACK
    if topic in ("decision_tree", "decision", "choose", "tree"):
        return DECISION_TREE
    if topic in ("history", "genealogy", "lineage"):
        return HISTORY
    if topic == "all":
        return "\n\n".join([LAB_STACK, DECISION_TREE, HISTORY])
    return (f"Unknown topic '{topic}'. Available: lab_stack, decision_tree, "
            "history, all.")
