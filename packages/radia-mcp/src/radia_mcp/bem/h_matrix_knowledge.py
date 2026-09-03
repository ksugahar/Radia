"""H-matrix and ACA acceleration for BEM / HDiv-VIM."""

OVERVIEW = r"""
# H-matrix / ACA / BLR for BEM acceleration

Dense BEM matrices have N² memory and N²·M cost per MatVec, prohibitive
for N > 10,000.  Two main acceleration families:

| Family | Memory | MatVec | Lab use |
|--------|--------|--------|---------|
| **FMM** (Greengard-Rokhlin) | O(N) | O(N) | Removed (CLAUDE.md) |
| **H-matrix + ACA** | O(N log N) | O(N log² N) | ★ HACApK (lab core) |
| **BLR (Block Low-Rank)** | O(N^{3/2}) | O(N^{3/2}) | for direct solver |

## Why lab uses H-matrix over FMM

Per CLAUDE.md "FMM (Fast Multipole Method): Removed (2026-03-06)":
1. FMM dipole approximation is inaccurate for distributed face charges
2. Compact geometries (87% near-field pairs) eliminate FMM's advantage
3. HACApK provides same O(N log N) memory with better practical performance
4. FMM tree-build overhead dominates for typical Radia N < 10,000

So Radia ships HACApK only (Greengard-Rokhlin FMM doc kept here for
reference completeness).

## Reference (foundational)

- Hackbusch 2015 "Hierarchical Matrices: Algorithms and Analysis"
  [LOCAL] 06_h_matrix_aca/Hackbusch_Hierarchical_Matrices.pdf
- Bebendorf 2000 "Approximation of boundary element matrices"
  (the original ACA paper, [DOI] in 00_references.bib)
"""


ACA = r"""
# Adaptive Cross Approximation (ACA) — Bebendorf 2000

Reference: M. Bebendorf, "Approximation of boundary element matrices",
Numerische Mathematik 86(4):565-589, 2000.  DOI: 10.1007/PL00005410.

## The idea

For an off-diagonal block A_block ∈ R^{m×n} of a BEM matrix where the
source cluster and target cluster are well-separated (admissibility
criterion), the kernel G(r, r') is smooth.  So A_block is **low rank**:

```
A_block ≈ Σ_{k=1}^{r} u_k v_k^T,    r << min(m, n)
```

ACA finds u_k, v_k in O(r·(m+n)) by greedily selecting "pivot" rows
and columns without computing the full A_block.

## Algorithm sketch

```
A_residual = A_block
for k = 1, 2, ...:
    Choose pivot (i, j) (e.g. max residual)
    u_k = A_residual[:, j] / A_residual[i, j]
    v_k = A_residual[i, :]
    A_residual -= u_k v_k^T
    if ||u_k|| · ||v_k|| < ε · ||A_block||: break
```

Termination: rank r determined automatically by tolerance ε
("adaptive").

## Properties

| Property | Value |
|----------|-------|
| Storage | O(r·(m+n)) vs O(m·n) full |
| Setup cost | O(r²·(m+n)) — does NOT need full matrix |
| Algebraic (kernel-agnostic) | YES (works for any G(r,r')) |
| Convergence | Geometric in r for asymptotically smooth kernels |
| Lab uses | HACApK applies ACA on each admissible block |

## Failure modes

| Issue | Symptom | Fix |
|-------|---------|-----|
| Bad pivot choice | Slow convergence in r | Use ACA+ (better pivots) |
| Non-smooth kernel | Rank doesn't decay | Use HCA or other compression |
| Near-singularity (touching clusters) | Mediocre rank | Use FMM for near-field instead |

References (this folder):
- [LOCAL] 06_h_matrix_aca/Bebendorf_ACA_solver_v1.0.pdf
- [LOCAL] 06_h_matrix_aca/Bebendorf_Approximation_BE_matrices.pdf
- [LOCAL] 06_h_matrix_aca/ACA_for_ill_posed_problems.pdf
- [LOCAL] 06_h_matrix_aca/Improvement_H_matrix_ACA_LargeScale.pdf
"""


HACAPK = r"""
# HACApK library — Sugahara lab production H-matrix

★ THIS IS THE LAB CORE H-MATRIX LIBRARY.

Source: `src/ext/HACApK_LH-Cimplm/` in the Radia monorepo (MIT license).
Origin: Post-Peta CREST project (Ida-Iwashita-Tohoku U).

Reference (this folder):
[LOCAL] 06_h_matrix_aca/HACApK_manual_v1.0.0.pdf — official manual
[LOCAL] 06_h_matrix_aca/HACApK_PostPeta_CREST.pdf — CREST overview
[LOCAL] 06_h_matrix_aca/Lattice_H_Matrices_Distributed.pdf — distributed
[LOCAL] 06_h_matrix_aca/HACApK_階層型行列法とライブラリ.pdf — JP intro

## How Radia uses HACApK

```python
import ngsolve as ng
from radia import vim

with ng.TaskManager():
    result = vim.Solve(
        mesh,
        mu_r=1000.0,
        H_ext=applied_field,
        gram_eps=1e-4,
        leaf=10,
        eta=2.0,
        linear_solver="auto",
        preconditioner="auto",
    )
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `gram_eps` | route default | Charge-Gram truncation tolerance |
| `leaf` | 32 | Minimum cluster size |
| `eta` | 2.0 | Admissibility parameter |

## Selection boundary

HACApK remains a production H-matrix implementation for charge-Gram, PEEC,
and BEM operators. It is no longer selected through the legacy
`rad.Solve(method=2)` integer. The owning formulation chooses H-matrix
construction and its linear solver through named APIs and reports the selected
backend in its result diagnostics.

## Performance

Per CLAUDE.md "H-Matrix Acceleration (HACApK)":
- Memory: O(N log N) vs O(N²) full
- MatVec: O(N log² N) vs O(N²)
- Tested on Radia models up to N ≈ 100,000 elements
- ~10-100× faster than ExaFMM for compact geometries

## Sign convention (critical, see CLAUDE.md)

All kernel functions return **+N** (positive physical demag tensor).
Sign flip to **-N** for the system matrix happens ONLY in
`ComputeEntry()` in `rad_hacapk.cpp`.

## BLR (Block Low-Rank) for direct solver

Recent extension: use BLR (a simpler H-matrix variant) for direct
LU factorization on the H-matrix.

References (this folder):
- [LOCAL] 06_h_matrix_aca/BLR_MBGS_QR_Python.pdf
- [LOCAL] 06_h_matrix_aca/QR_BLR_Weak_Admissibility.pdf
- [LOCAL] 06_h_matrix_aca/TensorCore_BLR_TallSkinny_BEM.pdf
- [LOCAL] 06_h_matrix_aca/Distributed_memory_lattice_H_factorization.pdf

Not yet in lab production; use the validated owning formulation and its
reported backend rather than assuming one Krylov method.
"""


def get_h_matrix_knowledge(topic: str = "overview") -> str:
    """Dispatch H-matrix / ACA topics.

    Topics:
        overview    - H-matrix family overview vs FMM (DEFAULT)
        aca         - Adaptive Cross Approximation (Bebendorf)
        hacapk      - ★ HACApK library (lab production)
        all         - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "family", "h_matrix"):
        return OVERVIEW
    if topic in ("aca", "adaptive_cross", "bebendorf"):
        return ACA
    if topic in ("hacapk", "lab", "lab_core"):
        return HACAPK
    if topic == "all":
        return "\n\n".join([OVERVIEW, ACA, HACAPK])
    return f"Unknown topic '{topic}'. Available: overview, aca, hacapk, all."
