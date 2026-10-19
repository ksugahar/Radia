# Nonlinear BH HDiv GPU acceptance contract

Status: implementation and numerical acceptance pending. Existing 10,344 DoF
seconds reported in this directory are linear, not nonlinear BH results.

## Preserve the CPU formulation

The independent reference is `_solve_nonlinear_energy_cpp` in
`src/radia/vim/_solve.py`, with `_EnergyMaterialQuadrature` and
`_bh_inverse_funcs` in `src/radia/vim/_nonlinear.py`.
For fixed geometry, retain the resident charge Gram action N = B^T G B.
With M the NGSolve mapped HDiv field, use

- R(m) = integral H(M).v dx + N m - rhs;
- J(m) = W_tangent(m) + N;
- W_tangent = integral [nu_sec I + (nu_diff-nu_sec) Mhat Mhat^T] u.v dx;
- energy = integral W_co(|M|) dx + 0.5 m^T N m - rhs^T m.

Use the same NGSolve integration rule for energy, residual and tangent:
IntegrationRuleSpace order max(6, 2*p+2). Do not infer nodal field values from
HDiv coefficients, replace the tensor tangent by a scalar susceptibility,
or replace near/singular Green integration with a point cloud.
NGSolve owns mapped shapes, orientations, weights and sparse assembly.
The demagnetizing quadrature must meet the CPU energy/Newton contract;
linear-only Gram evidence cannot certify nonlinear energy positivity.

Preserve the CPU monotone forward magnetization PCHIP, retained strictly
monotone inverse sample grid, inverse PCHIP antiderivative, zero-field tangent
and continuous hard-saturation barrier. Validate MATLAB/native law
construction against independent Python; Python is never a GPU user adapter.
A fixed sparse SPD polynomial may precondition inner Newton PCG only when
its assumptions hold. It does not replace the material tangent.

## Staged implementation

1. Add bounded native NGSolve material evaluation/assembly for an explicit
   coefficient vector and explicit inverse-BH polynomial contract. First
   verify energy/gradient/Hessian consistency and mapped FE parity. State the
   host/device placement and all per-Newton transfers; no implicit fallback.
2. Reuse one checked resident H-matrix across accepted Newton iterates.
   Solve (W_tangent+N) delta = -R with original-equation inner true residual.
   Keep candidate coefficients separate from accepted coefficients.
3. Preserve Armijo energy descent, saturation handling, failed-line-search
   rejection and assembled nonlinear residual relative to the fixed source.
   A small step or settled energy alone cannot certify convergence.
4. Test reset, failed material input, failed inner solve and failed line
   search without overwriting accepted state. Destroy only owned contexts.

## Public input and comparison

Start with the already public 13-tetrahedron P2 sphere, HDiv order2 (288 DoF).
Use an independently authored synthetic monotone SI curve, for example
M(H)=Ms*chi0*H/(Ms+chi0*H), Ms=1e6 A/m, chi0=1000, with declared positive
H knots [0,1,10,100,1000,10000,100000] A/m and B=mu0*(H+M) tesla.
Both routes consume the same declared table and CPU production interpolation
contract; do not use private engineering or commercial comparison data.
Cover weak, intermediate and saturated excitation separately, including
non-collinear vector excitation and tensor-tangent directional checks.

Require independent CPU coefficient/field agreement, final assembled
nonlinear residual, every inner true residual and energy/line-search history.
Use existing declared tolerances, not iteration count as acceptance. Begin
with CPU/RTX agreement, then same-artifact A100 runs under verified idle load.
After acceptance, scale through existing public meshes; precheck geometry,
charge counts, cache, sparse staging and solver memory before any cap change.

Benchmark cold initialization separately from warmed Newton workflows.
Whole time includes FE/material setup, initial H-matrix construction, every
material update, transfer, inner solve, line-search evaluation, result
recovery and destruction. Also report amortized reuse for fixed-geometry
load sequences with construction counted explicitly. Save source/input/law/
binary hashes, iteration counts, original residuals, exact owned buffers and
sampled memory limitations in JSON; payloads above2GiB require HDF5.
No claim that ten Newton steps cost ten times the linear whole-workflow time.
