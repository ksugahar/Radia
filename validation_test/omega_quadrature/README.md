# Mixed Omega quadrature re-solves

Run this heavy validation on mdx or hibino, not as part of fast CI. This lane
does not grant #6 or three-engine field acceptance. No p=4 run is needed.
It uses the production mixed solver, requests the assembled system explicitly,
and solves again for every assembly bonus and response order. The source Hodge
projection is constructed once and held fixed across the sweep.

`bonus_intorder` controls volume stiffness, source loads AND interface terms.
This is not an isolated load-only intervention. The independent term audit
splits air, iron and Kelvin loads, reconstructs the default assembled primal J,
then contracts and integrates with an explicit fixed volume quadrature rule.
The latter checks W-J-C on the newly solved field; it is not a claim that the
fixed-rule system was solved. Repeat with increasing evaluation order before
claiming quadrature convergence. Observed changes are not error bounds.

Example (paths refer to the compute host):

```text
python run.py --mode research --factory ctype.py --mesh C:/temp/radia-ctype-family/coarse/kelvin_domain.vol --orders 1 2 3 --bonuses 4 8 12 --evaluation-order 16 --output C:/temp/omega-audit/result.json
```

The C-type factory consumes the sibling mesh_result.json. Use only a mesh that
has already passed check-vol and its label contract. A factory for another case
must return the same keys and expose every physical setting in `controls`.
Factories must build the physical source, never fit a previously computed field.

`--mode wheel` requires a non-editable distribution and matching import path.
It forbids `--research-solver`. Research mode may explicitly load a Python solver
from a clean checkout with `--research-solver <path>`; this does not change the
installation or copy native binaries. Such output is research evidence ONLY.
Each output records the native hash, every package Python hash, solver/factory/
runner hashes, physical controls and mesh hash. Source identity is checked again
at each checkpoint. An incomplete file has no `completed: true` flag.

Embedding diagnostics include volume value, gradient, interface trace, all
nonfree coefficients (inactive coefficients and the zero gauge), multiplier
residual, and self-embedding controls. The interface residual must not be read
as a non-nesting proof until all embedding/gauge checks are satisfactory. The
current factory assumes homogeneous gauge; nonzero prescribed potential needs
a separate factory and lifting audit. Embedding integration uses order 16.

Diagnostic gates retain relative residual 1e-8 and energy consistency 1e-10.
Missing residual denominators are undecidable, not automatic passes. Exit 2
means a diagnostic gate failed; exit 0 means these diagnostics passed, not that
field agreement, mesh convergence, or release acceptance was established.

Before #6 release acceptance: re-solve at converged assembly/evaluation rules,
check region loads and embedding, decide the multiplier-space contract, and
compare HDiv, reduced-A and mixed Omega field probes using the SAME formal
candidate wheel and physical model. C-type is a control, not an #6 substitute.

## Small research validation, 2026-09-11

Saved JSON files in this directory are research evidence, NOT formal-wheel
acceptance. mdx2 used Radia 4.95.81 staging with the clean-main b5a54e473
kelvin_solver.py explicitly overlaid (Python only), NGSolve 6.2.2606, 24 threads.
The C-type coarse mesh has 24,134 elements. Source order was 3, evaluation order
8; this low-cost check verifies the lane, not final quadrature convergence.
Each JSON includes the exact runner/native/Python hashes used. After these runs
the runner additionally gained a mesh/dependent-source hash recheck at saves;
those checks do not retroactively appear in these older result schemas.

* `omega_quadrature_small_20260911.json`: actual p=1 solves at bonuses 4 and 8.
  Solve times were 2.49 and 5.94 s, including neither projection nor subsequent
  audits. W-J-C was 5.73e-14 and 3.40e-14. The bonus-4 iron-block gate FAILS:
  RHS norm 1.02e-15, residual norm 2.37e-16, relative 0.233. This is recorded
  without relaxing the gate; a near-zero block RHS is a normalization issue
  to inspect, not grounds for quietly marking the run passed.
* `omega_quadrature_embedding_p12_20260911.json`: bonus 8, p=1 and 2 re-solved;
  all algebraic/energy gates passed. Cross-order value error <=8.36e-16,
  interface trace error <=1.19e-15, gradient error <=9.09e-13 (relative).
  Nonfree coefficients <=7.24e-15. Interface constraint relative violation
  remains 3.3562e-4; self controls are 1.85e-16 and 4.10e-16. This supports a
  non-nested assembled constraint set in this C-type control; it is not #6
  acceptance or a variational-bound proof.

The original lane had 8 passing gate tests. This describes the scope of the
2026-09-11 small runs only, not the status of subsequent remote campaigns.
The formal candidate remains the release coordinator's responsibility.

## Recovery review, 2026-09-12

`--algebraic-only` is an explicit short diagnostic: it still re-solves every
requested assembly rule, records the original RHS-relative residual and the
norm of each column-block action on each free row block, and samples the
ESRF shared eight-point magnetic-field observable. It omits energy and embedding
audits, records their absence, and exits 2 even when the solve succeeds. It is
not an alternative release gate. The additional residual scale is
`||b_i|| + sum_j ||A_ij x_j||`; this measures cancellation between block actions,
not componentwise backward error, discretization error, or a field-error bound.
No existing residual threshold is changed and historical JSON is unchanged.

The gate now rejects nonfinite energy operands and negative residual norms.
An infinite reference energy must not turn an invalid comparison into a pass.
Saved results retain their original hashes and gate values; changing the
validator does not retroactively certify them. The small-run JSON files above
remain HOLD evidence. No new solve or native validation was performed during
this recovery review.
