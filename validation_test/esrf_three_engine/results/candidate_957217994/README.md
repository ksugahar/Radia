# ESRF6 candidate input evidence

Status: **HOLD**. These files identify the candidate and selected meshes; they
do not certify field agreement, quadrature convergence, or performance.

## Candidate

- Native build source: `95721799420952a446d44793fc1baf83fa6c92f9`.
- CI run: `34703027075`, successful build/import/native smoke/wheel production.
- Wheel SHA256: `a4c64938ac814fc12b6631d7ae1e2bddf5da01efd25238ddda052ea080ccafe7`.
- Native SHA256: `f99ec26e15d880431ed560c3f5cac7a537056e66f1bea0872186cfe32c2158a0`.
- The isolated Hibino venv passed `pip check`. `preflight.json` records the
  installed Python/native files checked against the actual wheel bytes.
- No native binary was copied into an existing installation.

## Which mesh was used?

The old nonconforming iron mesh had SHA256
`0d4449962b60757cba78e771a136b2805a4323e26a1ee6d7ae482751bb7e5a5f`.
Commit `caa4294ab1fc47df82dc4899c2fdf8fd71ce9f09` fixed the unmerged Cubit
partitions and replaced it with the 2408-HEX mesh `fdd13872...` recorded in
`case6_mesh_contract.json`. The selected input matches that post-fix artifact,
not the old 4768-HEX mesh with duplicated faces and hanging contacts.

The FEM input is a separate conforming TET generation route, SHA256
`dfc12b84...`, with 134686 elements. Its September 7 generation record references
the corrected iron asset manifest `b2bbe5db...`; the mesh bytes also equal the
September 4 local-gap artifact. The generation record, byte identity, and
present topology checks are separate evidence. They do not prove when those
bytes were originally generated or that their resolution is sufficient.

`mesh_revision_audit.json` records a fresh read-only check on Hibino:

- Both selected meshes: zero hanging facets and duplicated boundary-face pairs.
- FEM: 1393 Kelvin node pairs; maximum translation error `9.60e-15 m`.
- Periodic H1, orders 1 and 2: 1393 and 5566 constrained free DOFs.
- The **constant** trace squared-integral ratios differ from one by at most
  `3.6e-15`. This sanity check is not a proof for all nonconstant trace modes.

`fem_check.json` and `iron_check.json` are strict structural/label checks.
Topology, volume and area agreement do not prove geometric fidelity or field
convergence. The C-type N=6/12/24 gap family is a different geometry and must
not be substituted for this ESRF6 quadrupole.

## Input selection gate

Run the standard-library-only checker before launching this fixed-input lane:

```text
python validation_test/esrf_three_engine/check_case6_inputs.py --iron <selected-iron.vol> --fem <selected-fem.vol> --output <new-run>/input_identity.json
```

It rejects unregistered bytes regardless of filename. There is no force option.
A new mesh family needs a reviewed contract and new numerical evidence.
`input_identity.json` here was also measured during the running candidate audit;
it is not falsely labelled a pre-launch check. The two existing runs had already
checked their input hashes through their runner/preflight provenance.

Old scratch files are not deleted blindly because other tasks may reference
them. Their mere presence is not evidence that the current runner opened them.

## Completed candidate runs: acceptance remains HOLD

Both foreground Hibino jobs completed with exit code zero. The original result
JSON files are retained without changing their reported flags:

- `three_engine_case6_bonus8_p2.json`, SHA256
  `40ab46e2d2aa2676633cc8731005b666b073df44ff8ff41347770cc809332e21`.
- `omega_bonus8_p1_q16_22.json`, SHA256
  `af29048e10e2cb19959b2ec6ac67b45ec6cdcdc2481b056101b70f9bf7ecbf65`.

The nonlinear run used BDM2 for HDiv, FEM order 2, mixed source order 3,
exact exterior source, mixed assembly bonus 8, and tolerance `2e-5`.
The core observations are the runner's 27 gap points, not the whole field.

| Pair | Core relative RMS | Maximum core vector difference (T) |
| --- | ---: | ---: |
| HDiv / reduced-A | 0.0036015580 | 0.0010756475 |
| HDiv / mixed total/reduced Omega | 0.0070167627 | 0.0017637700 |
| reduced-A / mixed total/reduced Omega | 0.0078074982 | 0.0020484297 |

### Release blocker: HDiv false nonlinear convergence

The raw runner reports `passed=true` and `nonlinear_converged=true`. These are
**not accepted**: HDiv reported one Newton iteration, 34 Armijo backtracks and
relative step `1.091285e-12`. In this candidate, exhausting the line search
below `1e-10` still applies the tiny unaccepted step, then declares convergence
from step size. Thus the runner's Boolean is not evidence of a nonlinear root.
The separate settled-iteration exit also lacks a nonlinear residual gate.

The repair must address fail-loud line-search exhaustion and independently
check residual convergence. Energy/residual/tangent consistency is an additional
investigation: the inverse BH energy uses a separate interpolated trapezoidal
integral, while the residual and tangent use other derivatives; the HDiv path
also evaluates coefficients from an elementwise projected magnetization norm.
These are source-level concerns, not proven sole causes of this run's failure.
A repaired candidate requires new provenance and a new numerical run.

Reduced-A reported 31 iterations and final relative change `1.859221e-5`.
Mixed Omega reported 14 iterations and relative B change `7.876665e-6`.
Neither observation repairs the HDiv acceptance defect or proves absolute
accuracy. The complete three-engine acceptance remains **HOLD**.

### Frozen-state quadrature audit

This is a separate **linear**, order-1, permeability-1000 run with assembly
bonus 8. It solved once and evaluated that same state at q=16 and q=22.
Source and mesh identity remained unchanged; every recorded algebraic and
same-rule identity gate passed.

| Evaluation q | W | J (same rule) | W - J - source offset |
| --- | ---: | ---: | ---: |
| 16 | 8.47725441135 | 8.80958313240 | 7.67e-14 |
| 22 | 8.47699289981 | 8.80932051468 | 6.17e-14 |

The observed W change is `-2.6151154e-4` (about `3.085e-5` relative to q=16).
This is an observed evaluation change, not an error bound. The assembled
functional is `8.80306880342`, still different from the q=22 functional by
`0.00625171126`. The air source-load contraction changes from
`-0.65772528574` (assembly) to `-0.66397947200` (q=22), whereas the quadratic
part changes by only about `1.16e-8`. Source-load assembly quadrature must
therefore be varied with **new solves** before declaring it converged.
This frozen-state audit neither establishes p-nesting nor nonlinear acceptance.

### Timing scope

The jobs overlapped on Hibino, so these are run diagnostics, not comparative
performance claims. HDiv took 6246.7 s, including about 6134.1 s for ChargeGram
construction and 96.8 s for its reported solve. Reduced-A took 2253.2 s and
mixed Omega 16194.3 s. The high-order HEX Gram build is the dominant cost in
this run; the falsely converged HDiv solve time is not a valid solver benchmark.
