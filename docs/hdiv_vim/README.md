# HDiv-Type VIM

HDiv-VIM is Radia's production soft-iron demagnetization route.  It keeps the
magnetization, material state, charge map, and reduced-FEM handoff in NGSolve
mesh/function-space vocabulary, then uses Radia's C++ charge-Gram H-matrix for
the open-boundary integral operator.

## Live API

```python
import radia as rad
import radia.vim as vim

iron = vim.MeshSoftIron(mesh, mu_r=1000)
model = rad.ObjCnt([iron, source])
result = rad.Solve(model, 1e-3, 100, 2, demag_backend="hdiv")
field = rad.Fld(model, "b", [0, 0, 0.02])
```

For direct VIM use, call `radia.vim.Solve(mesh, mu_r=... | bh_table=...,
H_ext=..., image=...)`.

## Why This Is The Main Route

- `N = B^T G B` is symmetric and loop-free by construction: loops are
  `ker(B)`, so charge-free modes are field-null without hand-built loop bases.
- NGSolve `Mesh`, `GridFunction`, `CoefficientFunction`, `BilinearForm`, and
  `TaskManager` are shared with reduced FEM and motor workflows.
- Curved/high-order geometry uses the same finite-element geometry path instead
  of translating through a separate object discretization.
- TET/HEX/WEDGE and 2D planar validation live in `validation_test/feec/`.

## Validation Expectations

Fast tests belong in `tests/`; heavier numerical checks belong in
`validation_test/feec/` and mdx-labelled JSON artifacts.  Required HDiv gates:

- analytic demag factors where available;
- nonlinear BH convergence metadata;
- image symmetry checked against an explicitly mirrored full model on truly
  symmetric meshes;
- `rad.Fld` after `rad.Solve(..., demag_backend="hdiv")`;
- charge-Gram H-matrix build stats and memory/timing on mdx for large runs;
- 2D planar motor saliency checks for the motor lane.

## Public Docs

Current user-facing notebooks in this directory should be result-bearing and
paired with synchronized JSON sidecars.  Old migration notes and comparison
archives are not part of the public docs surface; recover them from git history
only if needed.
