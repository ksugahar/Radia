# NGSolve 6.2.2604+ `ngsolve.bem` integration

NGSolve 6.2.2604 (released 2026-04-30) ships a substantially upgraded
`ngsolve.bem` module with both **FMM-style** kernel-specific
acceleration AND a **kernel-agnostic block-entry callback bridge**
(via Pierre Marchand's HTool work).  Our SF framework's callback
contract is structurally compatible with the latter.

## Status (LAB, 2026-05-30)

```
> pip show ngsolve
Name: ngsolve
Version: 6.2.2604
Location: C:\Program Files\Python312\Lib\site-packages
```

## `ngsolve.bem` public API

```python
import ngsolve.bem as bem
# 36 public symbols:
```

**Operators** (Galerkin BEM bilinear forms):

  - `LaplaceSL`, `LaplaceDL`
  - `HelmholtzSL`, `HelmholtzDL`, `HelmholtzCombinedFieldOperator`,
    `HelmholtzHypersingularOperator`
  - `MaxwellSL`, `MaxwellDL`, `MaxwellSingleLayerPotentialOperator`,
    `MaxwellSingleLayerPotentialOperatorCurl`,
    `MaxwellDoubleLayerPotentialOperator`
  - `LameSL`
  - `HypersingularOperator`
  - `SumOfPotentialOperators`, `SumOfPotentialOperatorsAndTest`
  - `IntegralOperator`, `PotentialOperator`, `BasePotentialOperatorAndTest`
  - `SingleLayerPotentialOperator`, `DoubleLayerPotentialOperator`

**FMM-style kernel CFs**:

  - `BiotSavartCF(order, kappa, center, rad)`
  - `BiotSavartRegularMLCF`, `BiotSavartSingularMLCF`
  - `RegularMLExpansion`, `SingularMLExpansion`, `SingularMLExpansion3`
  - `RegularExpansionCF`, `SingularExpansionCF`
  - `SphericalHarmonicsCF`, `Sphericalharmonics`
  - `HelmholtzCF`, `PotentialCF`

**H-matrix bridge** (added 2026-04-20 commit `d90c59e` by Lackner,
crediting Marchand):

  - `IntegralOperator.CalcSubMatrix(rowids, colids) -> Matrix`
  - `IntegralOperator.CalcSubMatrixCapsule() -> capsule`
  - `IntegralOperator.NearFieldMatrix() -> BaseMatrix`
  - `IntegralOperator.mat` (full matrix)
  - `IntegralOperator.GetPotential(gf, intorder, nearfield_experimental)`

**Helpers**:

  - `GetDofCoordinates`

## Canonical usage pattern

From `tests/pytest/test_bem.py` upstream:

```python
import ngsolve.bem as bem
from ngsolve import H1, ds

fes = H1(mesh, order=2, dirichlet="boundary")
u, v = fes.TnT()

# SURFACE-TO-SURFACE Galerkin matrix
op = bem.LaplaceSL(u * ds, use_fmm=False) * v * ds
mat = op.mat                                 # full assembly

# Block extraction for external H-matrix library
rows = np.asarray([0, 5, 10, 15], dtype=np.int32)
cols = np.asarray([0, 1, 2, 3], dtype=np.int32)
submat = np.asarray(op.CalcSubMatrix(rows, cols))
```

## Scope limitation

The `CalcSubMatrix` callback pattern is **surface-to-surface Galerkin**:
both `rows` and `cols` index the SAME surface FES.  For SF INVERSE
design with off-surface point targets (rows = M target points, cols =
source FES DOFs), the `CalcSubMatrix` interface does NOT apply
directly.  Off-surface evaluation goes through

```python
potop(gfu, target_boundary)
```

which uses the FMM backend, NOT the H-matrix callback.

This means our SF inverse design (off-surface targets) stays on the
**(A) path**: callback into HACApK from a per-`(target, basis)`
integration (our current `demo_planar_uniform_fem_psi_aca.py`).  The
callback contract is identical in shape — `entry(i, j) -> float` —
just hosted by `radia.stream_function.aca_tsvd` instead of
`ngsolve.bem.IntegralOperator`.

For SURFACE-TO-SURFACE applications (= inductance extraction, mutual
inductance, FEM-BEM coupling), our framework can DIRECTLY use
`ngsolve.bem` operators.  These are the natural use cases for the
2604 H-matrix bridge.

## FMM vs ACA+ — do NOT conflate

Both are H-matrix-style accelerators but fundamentally different math:

| Aspect       | FMM                                  | ACA+                                    |
|--------------|--------------------------------------|-----------------------------------------|
| Math         | Kernel-specific analytic multipole   | Algebraic / kernel-agnostic             |
| Kernel       | One implementation per kernel        | Same code for any low-rank smooth kernel |
| Good for     | Smooth surface, far-field-dominated  | Compact / near-field-heavy / material   |
| In 2604      | `BiotSavartCF`, `*MLCF`, `*MLExpansion`, `SphericalHarmonicsCF` | `CalcSubMatrix`-via-HTool                |
| In Radia     | Removed (CLAUDE.md 2026-03-06)       | HACApK (kept and developed)             |

CLAUDE.md "FMM Removed (2026-03-06)" applies to Radia's own
HDiv-VIM volume integral, NOT to `ngsolve.bem`'s surface BEM
acceleration.  Different layer, different geometry class.

## How upstream activity affects us

If/when Joachim's H-matrix support in `ngsolve.bem` is fleshed out:

  - For surface-to-surface (standard BEM): drop-in replacement of
    `radia.stream_function.aca_tsvd` with `ngsolve.bem` H-matrix.
    Saves the per-entry numerical integration cost (huge for
    material kernels).
  - For SF inverse design (off-surface targets): KEEP the (A) callback
    path.  Our entry-function abstraction is independent of any
    H-matrix lib choice.
  - Either way, the rest of the pipeline (Path-A, single-stroke
    spiral, deformation outer loop) is unchanged.

## How to swap matrix assembly to ngsolve.bem (when applicable)

For surface-to-surface problems, the swap is ~5 lines:

```python
# Interim: LinearForm per target (current implementation)
A_free = build_fem_matrix(...)[:, free_idx]
def entry(i, j):
    return float(A_free[i, j])

# After ngsbem H-matrix matures, for surface-to-surface ONLY:
op = bem.LaplaceSL(integrand_for_source).Operator("name")
def entry(i, j):
    M = op.CalcSubMatrix(np.array([i], dtype=np.int32),
                         np.array([j], dtype=np.int32))
    return float(np.asarray(M)[0, 0])

# Rest (radia.stream_function.aca_tsvd, Path-A, ...) unchanged
```

## Cross-reference

  - FMM vs ACA+ deep dive: MCP topic
    `streamfunction(topic=session_2026_05_30)` sections 6–8
  - Joachim's GitHub commit history: 2026-04-20 `d90c59e` (HTool bridge),
    April-May erdieee FMM work
  - Memory entry: `feedback_fmm_vs_aca_distinction`
  - CLAUDE.md scope clarification: "FMM Removed from Radia core
    (2026-03-06) / SCOPE CLARIFICATION (2026-05-30)"
