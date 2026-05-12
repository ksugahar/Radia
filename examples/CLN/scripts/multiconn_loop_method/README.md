# Multiply-Connected T-Ω: Hiptmair-Ostrowski Loop Method

Reference implementation of the **Loop / cohomology-basis method** for
T-Ω eddy-current FE on multiply-connected conductors. Sourced from
`S:/NGSolve/EMPY/EMPY_Analysis/` (EMPY analysis framework, M. Tanimoto et al.)
on 2026-05-12.

## Why T-Ω needs special handling for multiply-connected geometry

The T-Ω formulation writes the current density as `J = curl T` and the
remaining magnetic field as `H = T - grad Ω`. For this to be globally
well-defined the conductor's interior must be **simply connected**.

For a multiply-connected conductor (a torus, a plate with N holes, a
geometry of higher genus g), the de Rham first cohomology group
`H¹(conductor; ℝ)` has dimension g ≠ 0. This means there are g independent
**non-exact closed 1-forms** that cannot be written as `grad Ω` of any
single-valued Ω. Each missing cohomology dimension corresponds to a "loop
current" that threads through a hole / handle of the conductor.

**Without g extra DOFs in the FE system, the T-Ω matrix is singular**
(or worse, silently gives the wrong answer when the loop current is
non-trivial physically).

The Hiptmair-Ostrowski Loop Method augments the FE system with `g`
extra scalar unknowns — the *amplitudes* of an explicit basis of the
cohomology classes — and couples them to the standard FE block.

## Files

| File | Purpose |
|---|---|
| `include/LoopField.py` | Genus computation + cohomology basis construction (`LoopFields()`) and FE coupling matrix builder (`loopFieldCouplings()`). |
| `include/MatrixSolver.py` | Bordered-system solver. `AddCoupling()` builds the augmented sparse matrix; `SolveCoupled2()` uses ICCG on it. |
| `include/HtoOmega.py` | Helper used in the notebook for the H-to-Ω reduction (reduced/total decomposition). |
| `T_Omega2_BathPlate_with_Holes_Bn-reg.ipynb` | Canonical user-facing example: bath plate with N circular holes. |

## Algorithm (LoopField.py:LoopFields)

1. **Genus** via Euler characteristic on the conductor boundary surface:
   ```
   χ = nv − ne + nf
   g = p − χ/2      (p = number of connected components of ∂conductor)
   ```
   `surface_genus(smesh, p)` returns g.

2. **For each k = 0..g−1, construct one cohomology basis field**:
   - Pick a random edge of the conductor surface and fix its HCurl DOF to 1
     (`fes.FreeDofs().__setitem__(edge_dofs, False)`).
   - Solve `∫_air curl u · curl v = 0` with that fixed DOF and outer-boundary
     Dirichlet. Gives a **closed** (curl-free in air) but non-exact field
     `gfu` — it has non-trivial circulation around hole `k`.
   - Helmholtz-decompose: solve `∫ ∇φ·∇ψ = ∫ ∇ψ·gfu` and let
     `gfw = gfu − ∇φ`. This removes the exact (gradient) part.
   - Gram-Schmidt orthogonalize `gft = gfw − Σ_{j<k} ⟨gfw, loops[j]⟩ loops[j]`.
   - Normalize `gft ← gft / ‖gft‖_{L²}`.
   - `loops.append(gft)`.

3. **Coupling** (`loopFieldCouplings`):
   For each loop k, build an RHS vector `fv[k]` (column of g extra rows
   coupling to the FE T-Ω unknowns) and a g × g Gram matrix `fafv` from
   the bilinear form
   ```
   1/(s σ) ∫_cond curl T · curl gfT[k]
     + μ ∫_cond T · gfT[k]
     + μ ∫_air loopField[n] · grad ψ
   ```
   where the linear test is `(W + grad ψ)` for `(T, Ω)` and the
   amplitudes `amp[0..g−1]` are coupled to the standard FE block.

4. **Solve** the bordered system `[A | B; Bᵀ | C][u; amp] = [f; 0]` via
   `MatrixSolver.SolveCoupled2`, which calls `AddCoupling` to assemble the
   augmented sparse matrix and then runs ICCG on it.

5. **Assemble the T field on the conductor boundary**:
   ```python
   loopsum = sum(amp[k] * loops[k] for k in range(g))
   gfT.Set(loopsum, BND, mesh.Boundaries(conductor_boundary))
   ```
   Then add the FE interior solution on top via `+=`.

## Dependencies

These scripts require the EMPY C++ bindings:
- `SparseSolvPy` (JP-MARs sparse ICCG solver)
- `EMPY_Solver` (EMPY internal solver)

Both live in `S:/NGSolve/EMPY/EMPY_Analysis/bin/Release/` (not bundled
here). For a portable re-implementation, the bordered system in
`AddCoupling` can be solved directly with `scipy.sparse.linalg.cg` or
`spsolve` — only the C++ solver layer is EMPY-specific.

The notebook also imports `BathPlateModel2` from EMPY's `..\model` —
again not bundled. Treat these scripts as algorithm reference rather
than directly executable here.

## When to use vs. when not to

**Need the Loop method when:**
- T-Ω formulation on a torus, plate with through-holes, gear teeth ring,
  any geometry of genus ≥ 1.
- Eddy currents that physically circulate around a hole (vs. just
  surface eddies).

**Don't need it when:**
- Conductor is simply connected (sphere, cylinder, cuboid, disk). For
  these the `tanimoto_canonical/CLN_T-Omega.ipynb` baseline is sufficient.
- The A-T or A-Φ formulation is acceptable — those formulations don't
  have the simply-connected restriction.
- The CLN extraction is only for the leading rung and the loop current
  is known to be physically zero (e.g., axisymmetric problem with
  azimuthally uniform driving).

## See also

- `../tanimoto_canonical/CLN_T-Omega.ipynb` — single-connected T-Ω
  reference (1 cm Cu cylinder).
- Memory: `reference_loop_method_multiconnected_TOmega.md` for the
  algorithm summary and cross-references.
- Hiptmair, R., "Discrete Hodge operators" (Numer. Math. 90, 2001).
- Bossavit, A., "Computational electromagnetism", Academic Press 1998 —
  classical reference on cohomology basis construction.
