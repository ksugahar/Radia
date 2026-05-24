"""BEM/MoM overview, decision tree, and lab stack."""

LAB_STACK = r"""
# Sugahara Lab BEM/MoM stack (production, 2026-05)

The lab uses BEM/MoM techniques in MULTIPLE places, each with a different
formulation and tool:

| Use case | Formulation | Tool | Files (this folder) |
|----------|-------------|------|---------------------|
| **Permanent magnets / soft iron** ★ | MSC (surface charge) | radia C++ `rad_polyhedron.cpp` | 02_mmm_surface_charge/ |
| **Volume magnetization** ★ | MMM (volume integral, M as DOF) | radia C++ `rad_interaction.cpp` | 02_mmm_surface_charge/ |
| **IH workpiece + Kelvin** | FEM-SIBC + Robin BC (NOT pure BEM, but BEM-like) | radia `calc_fem_kelvin.py` | 03_bem_eddy_current/, 05_low_freq_* |
| **PEEC filament+panel** | FastImp-style PEEC | radia `peec_matrices` | 80_applications/Fast_Impedance_*, 02/rna_mmm/PEEC_MMM_* |
| **PEEC + BEM coupling** | weak coupling, ΔL_telegen | radia `calc_inductance.py` | 02/rna_mmm/PEEC_MMM_Coupling.pdf |
| **HCurl Galerkin BEM** | EFIE/MFIE/PMCHWT with RWG | `ngsolve.bem` | 04_efie_mfie_cfie/, 05_low_freq_* |
| **H-matrix acceleration** | ACA on MMM/MSC blocks | radia `HACApK` (method=2) | 06_h_matrix_aca/ |
| **FEM-BEM hybrid (transformer)** | FEM inside, BEM outside | (not lab production, ref only) | 08_fem_bem_hybrid/ |

★ = Radia's CORE methods (MMM + MSC).  See CLAUDE.md "Terminology" for
the distinction from Galerkin BEM (ngsolve.bem).

## What the lab does NOT use

- **FMM** — Removed from Radia 2026-03 (see CLAUDE.md).  HACApK with ACA
  handles the same large-N problems with better practical performance.
- **Calderón preconditioners** — Not needed because MMM/MSC are
  well-conditioned by construction.  Only relevant for ngsolve.bem
  high-frequency scattering (not a lab focus).
- **Null-field / TDS** — Reference only; lab does not use thin-dielectric
  approximations.
- **Wire-grid models** — NEC-era; superseded by RWG triangular surface
  meshes everywhere.

## Why MMM/MSC instead of Galerkin BEM

MMM/MSC use **constant-per-element** M or σ (analogous to RWG with
constant-tangential basis).  Trade-off vs Galerkin RWG-BEM:

| Property | MMM/MSC (Radia) | Galerkin RWG-BEM (ngsolve.bem) |
|----------|-----------------|--------------------------------|
| Basis order | 0 (constant per element) | RWG (linear edge basis) |
| Element type | Tet/hex/wedge volumes (MMM) or faces (MSC) | Triangles on surface only |
| Material | Nonlinear iron, PM | Linear conductors |
| Open boundary | Natural (volume integral) | Natural (surface integral) |
| Kernel | Laplace 1/r only | Laplace OR Helmholtz |
| Acceleration | HACApK (ACA) | H-matrix or FMM |
| Lab production | YES ★ | Optional (for high-freq scattering) |

The choice was made because Radia's heritage is accelerator magnets
(permanent magnets + iron yokes), where MMM is naturally fast and
nonlinear-material-friendly.  Galerkin BEM is appropriate for high-frequency
scattering, where the lab uses NGSolve's `bem` module instead.
"""


DECISION_TREE = r"""
# Decision tree: which BEM/MoM technique for which problem

```
1. What's the SOURCE?
   ├── Permanent magnet (fixed M)
   │   → MSC in Radia (no Solve needed)
   │     → rad.ObjHexahedron(verts, [Mx, My, Mz])
   │
   ├── Soft iron (M unknown, BH curve)
   │   → MMM in Radia + nonlinear Solve
   │     → rad.MatSatIsoTab(BH_DATA) + rad.MatApl
   │
   ├── Coil current (known I)
   │   → Biot-Savart analytical or PEEC filaments
   │     → CoilBuilder or PEECBuilder
   │
   └── Eddy current source (frequency-domain σ, ω)
       └── See "What's the FREQUENCY?"

2. What's the FREQUENCY?
   ├── DC / static
   │   → MMM/MSC (Radia)
   │
   ├── Quasi-static (MQS, DC to 100 kHz)
   │   ├── Thin conductor with skin effect
   │   │   → PEEC + SIBC (radia.peec_matrices)
   │   │     → CoupledPEECSolver
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

3. What's the SIZE?
   ├── N < 500 elements
   │   → LU direct (radia method=0)
   ├── 500 < N < 2000
   │   → BiCGSTAB (radia method=1)
   └── N > 2000
       → HACApK H-matrix (radia method=2) — lab core for big BEM
         → see matrix_solvers MCP for solver theory

4. What's the GEOMETRY topology?
   ├── Simply connected
   │   → Standard MoM/BEM (any formulation)
   ├── Multiply connected (loops in conductor)
   │   ├── A-V Bíró-Preis gauge (FEM) or PEEC loop currents
   │   └── Tree-cotree decomposition for ngsolve.bem
   └── Open boundary (unbounded domain)
       ├── Radia MMM/MSC: natural (no PML/Kelvin needed)
       ├── FEM: Kelvin transformation + Radia evaluator
       └── BEM (ngsolve.bem): natural via radiation BC
```

## Common pitfalls

| Pitfall | Symptom | Fix |
|---------|---------|-----|
| Using EFIE at low freq | nan, ill-conditioned | Switch to Loop-Star or Weggler ST |
| MFIE with RWG below 1 MHz | wrong solution | Use Vico-Greengard MFIE-LF stabilization |
| MMM with non-uniform M expected | wrong integral | Use higher-order MMM or subdivide elements |
| HACApK on small N (<500) | overhead dominates | Use LU (method=0) instead |
| Galerkin BEM on nonlinear iron | doesn't converge | Use MMM (Radia native nonlinear) |
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
