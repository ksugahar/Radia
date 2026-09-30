# Open-boundary method map

Radia separates spatial truncation, evaluation of an exterior boundary operator,
and model-order reduction. Choose each from the physical problem. The core
magnetic interaction kernels remain Laplace/MQS kernels; the analytic wave DtN
utilities do not turn the core into a full-wave solver.

## Select the exterior operator

| Problem | Maintained route | Qualification |
|---|---|---|
| Static magnetic field in exterior air | Kelvin transformation or a Laplace DtN | No frequency-dependent memory is needed. |
| Spherical/circular, separable diffusion exterior | Analytic multipole DtN | Exact for each retained mode under the declared homogeneous-material assumptions. |
| Non-separable or material-containing exterior | Kelvin FEM DtN | Accuracy depends on mesh, geometry, material representation and retained boundary modes. |
| Separable radiating exterior | Analytic wave DtN utilities | A separate wave-boundary utility, not a Radia core wave formulation. |
| General radiating geometry | A separately validated wave solver and boundary treatment | Do not infer applicability from a spherical analytic symbol. |

A finite modal expansion still has truncation error. Agreement for individual
multipoles does not certify an arbitrary geometry or an unlimited frequency band.

## Exact separable DtN

`radia.open_boundary.eddy_dtn` evaluates the exact multipole symbol for a
homogeneous spherical diffusion exterior. `wave_dtn` evaluates the separate
wave symbol. Diffusion reverse-Bessel roots are roots in `q = sqrt(s)`,
not physical-time relaxation rates.

```python
import radia.open_boundary as ob
value = ob.eddy_dtn(2, 50j, R0=0.1, mu_sigma=1.0)
```

## Time-domain approximations and model reduction

`dtn_exact.sqrt_s_passive_poles` and `eval_sqrt_poles` provide a pole-residue
approximation over a declared band. Verify approximation error and stability in
that band. A finite fit of diffusion memory is not an exact all-frequency model.

`kelvin_dtn.kelvin_dtn_matrix` builds a material-aware boundary operator;
`steklov_spectrum` characterizes it and `band_rational_fit` fits a sampled scalar
response. A sampled rational fit does not by itself establish passivity.

For finite-dimensional electromagnetic systems, Radia uses PRIMA projection or
Foster modal realizations through their maintained APIs. These reductions are
validated against the full-order system over the operating band.

## Validation and interfaces

- [Continued-fraction validation](../../validation_test/open_boundary/test_dtn_continued_fraction.py)
- [Exact symbols and pole fits](../../validation_test/open_boundary/test_dtn_exact.py)
- [Kelvin DtN implementation](../../src/radia/open_boundary/kelvin_dtn.py)
- [Kelvin spectral guidance](../kelvin/DTN_SPECTRUM_COARSE_MESH.md)
- [Executed open-boundary notebook](open_boundary_demo.ipynb)
- MCP: `dtn_coarse_mesh` and `kelvin_transformation` expose the maintained
  boundary guidance. Consult their topic index for current names.

The former `dtn_cln` module, circuit-named helpers, and `mor_cln` tool are retired.
The analytic continued-fraction representation is restored: use
`continued_fraction_stages` and `eval_continued_fraction` for exact/high-order
non-reflecting boundaries. It terminates at n+1 partial quotients in sqrt(s). Historical
experiments remain in Git history; their measurements are not acceptance results
for a changed implementation. The numerical scope above is not a patent or
freedom-to-operate determination.
