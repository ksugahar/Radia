# Coupled HDiv-MMM x HCurl-VIM physical validation lanes

End-to-end physics locks for the coupled magnetization + eddy-current hybrid
VIM (`radia.vim` `NgsolveBDMEddyBubbleVIM` -> `CoupledHDivHybridVIMSystem`),
at the physical harmonic coupling scales (magnetic row `-K/mu0`, eddy row
`s K^H`; `solve_frequency` defaults since commit 6ed8d39d1).

| Lane | Reference | What it locks |
|------|-----------|---------------|
| `test_sphere_alpha_lane.py` | exact complex polarizability `alpha(omega)` of a conducting permeable sphere (`radia.analytical_formulas.sphere_complex_polarizability`) | eddy-branch transition shape (mu_r=1), coupled static + transition (mu_r=100), vanishing DC current from a static magnetization |
| `test_torus_cohomology_tau.py` | short-circuited ring law `(R+sL) I = -s Phi_ext` with the closed-form ring R and L evaluated on the MEASURED mesh cross-section | the genus-1 (H1 cohomology) loop class: single-pole shape, the CALIBRATION-FREE time constant `tau = L/R`, flux-freeze plateau, and the cycle-class DC zero |
| `test_hdiv_hcurl_transient_validation.py` | eleven-step moving manufactured transient with a changing magnetic operator | production-form snapshot/artifact contract, Joule-loss observables, and the discrete energy-balance gate |
| `validate_hcurl_eddy_bubble_disk.py` | public axisymmetric BEM modal spectrum plus live Q1/Q2 `radia.axifem` checks | 3-D TET HCurl eddy-bubble h/p convergence, passive CLN extraction, and the response-basis completeness gate |
| `validate_magnetic_conductor_disk.py` | regenerated mapped-HEX mesh plus axisymmetric Q2 and full 3-D HCurl references | quick/full live replay of the magnetic-conductor adjudication without tracking `.vol` files |
| `validate_mapped_hex_bdm2_hodge_reference.py` | finite-domain H1 Omega discrete Hodge projection on the same two-cell mapped BDM2 body | separates a valid BDM2 approximation space from the out-of-range mapped volume/surface charge operator, then feeds the body-restricted BDM2 response into a shared-mesh HCurl eddy-bubble solve |
| `validate_h1_hodge_bdm2_disk.py` | independently regenerated static axisymmetric Q2 disk | heavy 3-D TET BDM1/BDM2 plus H1 h/p ladder; the final BDM2/H1-p4 observable reaches 0.98% error |
| `validate_nonlinear_iron_esim_coupling.py` | regenerated same-region magnetic/conductive TET cube plus local 1-D nonlinear B-H cells | full HDiv-MMM/HCurl solve inside the local-ESIM Karl loop, field-amplitude dependence, passivity, Joule loss, and fixed-Gram replay |
| `results_magnetic_conductor_disk_adjudication.json` | fine axisymmetric Q2 volume reference plus independently refined 3-D HCurl and mapped-HEX HDiv solves | thin magnetic-conductor formulation adjudication, BDM1 h-convergence, BDM2 negative control, and frequency-route rejection |

Protocol notes (both lanes are documented compromises, re-examine before
tightening):

* The sphere lane fits ONE complex amplitude calibration `c` at the smallest
  `a/delta` point (it absorbs the reduced-basis drive/extraction amplitude
  convention, |c| ~ 1.43); every later number is calibration-frozen.  The
  torus lane's primary gate `tau = L/R` needs no calibration at all (the
  single-pole model is invariant under joint (A,R,L) rescaling).
* The volumetric eddy basis saturates at strong skin (`a/delta >~ 3`); the
  affected points carry explicit error CEILINGS documenting the limitation
  (the SIBC surface branch is the designed strong-skin path, not tested here).
* The torus mesh is capped under the analytic HCurl interaction limit
  (512 tets), which under-fills the wire by ~25%: the ring reference uses the
  measured cross-section (inscribed-polygon correction), and the sampled
  bridge-bridge kernel epsilon is set EXPLICITLY to ~half the sample spacing
  (the default heuristic mis-sizes the loop inductance ~3x at this
  coarseness).
* The conducting-disk lane freezes a useful negative control: one uniform
  vector-potential port remains about 5% wrong when only Krylov depth is
  increased.  The accepted lane uses `A`, `r^2 A`, `r^4 A`, and `z^2 A` as
  training ports, then identifies the physical pole by port residue rather
  than selecting the numerically largest time constant.  Spatial h-refinement
  is not a substitute for response-basis enrichment.
* The magnetic-conductor adjudication is deliberately bounded.  The fine
  axisymmetric Q2 reference and full 3-D HCurl A-form agree within 2%, and the
  mapped-HEX BDM1 static ladder converges to 1.22%.  The recorded BDM2 HEX
  ladder is non-monotone and is not promoted.  Both the explicit sampled and
  production direct-Q2 coupled HDiv-MMM plus HCurl eddy-bubble lanes decrease
  from about 5.38% to 3.34% to 1.90% error over 32/96/384 HEX while preserving
  non-negative loss and stale-route rejection.  At the fine level, direct-Q2
  differs from sampled by 8.03e-7 in the observable and 2.93e-6 in Joule loss;
  the maximum direct residual is 3.90e-16.  The direct path prunes only
  projection-roundoff zero monomials and replaces the sampled diagonal through
  a restricted global H-matrix apply, avoiding a dense sample-pair matrix.  The
  ladder must use the 18432-element fine Q2 reference: the 6656-element quick
  reference moves the fine coupled error to 2.19% and is a frozen reject.
* The mapped-HEX BDM2 failure is localized to the charge operator, not to the
  approximation space.  On the same 207 active BDM2 DoFs, H1 Omega orders 2
  and 3 produce discrete Hodge spectra bounded by 0.992 and 0.988 with no
  out-of-range modes; the mapped charge diagnostic reaches 1.2175 with three
  out-of-range modes.  `H1HodgeDemagOperator` also drives a two-mode BDM2
  response reduction on the 54-HEX air/body mesh; its generic CG residual is
  below `3e-12`, and the coupled HCurl eddy-bubble solve has residual below
  `3e-16` with positive Joule loss.  The generic CG path honors the
  body-restricted HDiv `FreeDofs()`.  The finite H1 box and sampled mixed
  interaction are contraction/mechanics gates, not open-boundary accuracy
  oracles.
* The independent static-disk accuracy lane separates that mechanics result
  from observable accuracy.  Against the 18432-element axisymmetric Q2 value
  `1.0916977441`, the 3-D finite-domain H1-Hodge sequence reduces relative
  error from 4.65% (coarse BDM1/H1-p2) to 2.24% (coarse BDM2/H1-p3), 1.44%
  (fine BDM2/H1-p3), and 0.98% (fine BDM2/H1-p4).  This is a strict h/p ladder
  for this disk, not a universal open-boundary or solver-superiority claim.
* The nonlinear iron local-ESIM lane gives one regenerated TET body the same
  magnetic and conductive label.  At 50, 1000, and 5000 A/m, every Karl loop
  converges in 2--5 updates, all local surface Grams remain passive, all Joule
  losses are positive, and fixed-Gram replay agrees with the converged mixed
  solution to machine precision.  It establishes nonlinear skin coupling to
  a fixed low-field HDiv bulk operator.  It does not yet establish a
  simultaneous ordinary bulk nonlinear B-H update, hysteretic/rotational skin
  state, or multidimensional corner cell.

Runtime: ~1-2 min per lane (LAB-runnable; correctness, not timing).
