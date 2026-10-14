# P1 panel ESIM in weak PEEC/BEM coupling

For a genus-0 workpiece, `calc_inductance.py --esim-per-panel` uses one
piecewise-constant impedance per BND triangle. It supports `intree-dense` and
`hacapk`, requires a B-H file, and keeps the coil current distribution fixed
(weak coupling). Surface P2, strong nonlinear coupling, and genus-1 nonlinear
ESIM remain unsupported. This does not restrict the thermal FEM order.

The scalar potential is P1; its tangential gradient is constant on each flat
triangle. The assembled source-side form is

```
K_Z[i,j] = integral Zs(x) grad_s(phi_i).grad_s(phi_j) dS
A = M/2 - DL + SL M^-1 K_Z / (i omega mu0)
P = real(phi^H Re(K_Z) phi) / 2
```

Thus Zs belongs inside the surface integral before the single-layer operator.
Multiplying BEM observation rows by local Zs is not this discretization.
VIM's impedance Gram follows the same weighted-integral principle, but its
current-mode matrices cannot be substituted for the BIE potential matrices.
NGSolve assembles the mapped P1 forms and local heat.

The API uses `radia.surface_impedance.PanelSurfaceImpedance(values)`, where
`values` follow exactly `mesh.Elements(BND)` order. Untagged variable nodal
arrays remain rejected. Output includes `esim_impedance_layout` and panel
centroids alongside panel impedance and field arrays; their length is the
number of BND triangles, not the number of vertices. The MATLAB Python BEM
adapter classifies this as an initialization/batch Python operation, not a
new native MEX kernel.

Each cell receives the solved panel peak |Ht|. The outer convergence check
uses the undamped constitutive mismatch; small relaxation alone cannot make
it pass. Unconverged cells or outer iterations fail before exporting heat.
Negative resistance and nonfinite impedances are rejected.

Local heat is `q = Re(Zs) |Ht|^2 / 2` for peak phasors. It is projected to a
positive lumped P1 field for the existing thermal handoff, preserving its
integrated power. Reaction includes the bilinear electric term
`phi_inc.T K_Z phi / (i omega I^2)` without conjugation. The independent
reaction/heat balance gate remains in force.

The saved `Ht.sol` contains sqrt(lumped mean |Ht|^2). Across a discontinuous
panel impedance, dividing projected heat by panel Zs cannot reconstruct it.
Likewise, applying a nonlinear law to this smoothed Ht is an approximation
to panel/quadrature evaluation.

## Accuracy and temporary comparisons

The previous copper, 1 kHz, mu_r=1 acceptance is a certificate for its declared
uniform SIBC cases. It does not establish an error bound for a different
material, geometry, frequency or nonlinear B-H law. A converged solution and
power balance are necessary checks but do not certify material accuracy.
The current variable-coefficient tests establish assembly, source placement,
constant reduction and conservative file handoff. A prescribed passive
nonlinear cell law tests orchestration; it is not a ferromagnetic calibration.

Reintegrating `Re(Zs(Ht)) |Ht|^2 / 2` on an existing uniform-field solution
is useful as a sensitivity estimate. Record it as frozen-field postprocessing:
it omits the field redistribution caused by local saturation and impedance,
has no guaranteed upper/lower bound, and cannot inherit a reaction-balance
or nonlinear-solution accuracy claim. Use peak/RMS conventions consistently
and area/quadrature weights, not unweighted averages. Prefer the unsmoothed
panel field when available.

An axisymmetric volumetric comparison should use identical geometry, current,
frequency, conductivity, B-H data and temperature, resolve the skin layer and
test mesh/exterior convergence. `radia.ih_axisym_coupled` currently updates
sigma(T) and linear mu_r(T), not B-H saturation. The separate
`radia.panels.calc_axisym_volumetric` implements elementwise B-H Picard updates
and is the existing starting point for that comparison. Its phasor model is
not a hysteretic or multi-harmonic reference. A nearly axisymmetric coil may
still have lead or pitch effects which require a separate 3D check.

Nonlinear first-harmonic SIBC models are approximations; magnetic saturation
can distort the waveform. See the primary
[finite BEM/SIBC formulation](https://www.compumag.org/Proceedings/2015_Montreal/papers/PA5-5.pdf)
and [nonlinear B-H impedance study](https://zaguan.unizar.es/record/147674?ln=en).
Validate local ESIM against a skin-resolved volumetric model over the intended
field range before claiming a magnetic-material power error.

## Focused validation

Run on an available computation/test host with job-local sources and inputs:

```powershell
python -m pytest tests/test_panel_surface_impedance.py tests/test_bem_sibc_acceptance.py tests/test_sibc_reaction_contract.py
python -m pytest validation_test/induction_heating/test_panel_esim_workflow.py
```

The second lane executes surface extraction, coil excitation, BEM solve,
reaction, positive heat projection and SOL save/read. The file-load boundary
and explicitly prescribed cell law are test substitutes. It does not certify
coil CAD, production VOL export, a real magnetic material, or tracked SLX/MEX
release acceptance.
