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

Runtime: ~1-2 min per lane (LAB-runnable; correctness, not timing).
