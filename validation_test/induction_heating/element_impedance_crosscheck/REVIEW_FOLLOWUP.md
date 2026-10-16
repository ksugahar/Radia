# Temperature-BH review follow-up

The reviewed source is `266580d23`; publication remains blocked by B1–B3.
The historical 20261013 results are provenance, not qualification of a new
wheel. The preparation branch `66a8ef16e` is excluded from integration.

| Finding | Implementation owner | Required closure |
|---|---|---|
| B1 | IH | Independent port-power/face-loss identity and nonmagnetic reactive-sign gate; proposed criteria in reaction_gate_plan.json |
| B2 | IH | Fresh installed-wheel numerical evidence with source, wheel and loaded MEX hashes |
| B3 | GPU | Exact-source IH CI build, same-wheel native/Simulink/Toolkit, Native Release |
| N1 | IH | Explain that the block preconditioner is built once from the base operator, without ESIM terms; changing Zs does not refresh it; true-residual checks still apply |
| N2 | IH, adapter change pending ownership confirmation | Reject observe_element_impedance=true unless material_mode is element-esim-operators during setup |
| N3 | IH | Mark zero initial/reset Zs as nonphysical placeholders in block/model descriptions; omit sample zero when plotting physical impedance |
| N4 | IH | Correct stale material-connection claim and distinguish the older preview from the temperature-BH workflow |
| N5 | GPU | Match required native-source inventory to the actual build inventory |
| N6 | GPU | Replace unverified fixed receipt assertions with checked evidence and remove private host/compiler paths from published metadata |

All rows are open until the implementation and its validation evidence are
reviewed. N2 touches
`packages/radia-gpu/matlab/+radia/+simulink/ihResidentSFunction.m`; N3 may touch
`buildTemperatureBHResearchModel.m`. Coordinate those exact files before
editing the candidate. Documentation targets are
`packages/radia-gpu/docs/IH_TEMPERATURE_BH.md`,
`packages/radia-gpu/matlab/IH_RELEASE.md` and
`validation_test/ih_cuda/README.md`.

The reactive-sign case deliberately excludes magnetic permeability effects.
A passive magnetic sample is not guaranteed to have a negative net change in
reactance. Its mandatory identity is real port power versus dissipated heat.
The separate nonmagnetic case makes a negative eddy-current reactance a
well-defined expectation; independent FEM supplies the second analysis route.

Execution waits for the final artifact and coordinated compute admission.
Small source edits and validation preparation can proceed separately.

## Recorded port-work diagnostic

`reaction_gate.py` evaluates records using the thresholds in
`reaction_gate_plan.json`, which remains proposed pending formula review.
It reads no assembled operator, impedance map or reaction-row coefficients.
For a new evidence directory, run:

```text
python reaction_gate.py records.json --output reaction-result.json --wheel candidate.whl --loaded-mex /absolute/path/radia_cuda_ih.mexw64 --source-commit FULL_SHA
```

The input object declares `phasor_convention` exactly as in the plan,
`case_kind` (`temperature-bh-zero-resistance` or `nonmagnetic-reactive`),
`supplied_coil_impedance_ohm` as [real, imaginary], and a nonempty
`accepted_updates` list. The nonmagnetic case additionally declares
`material_contract="constant-mu0-positive-resistivity-no-remanence-no-bh"`;
the acquisition driver and independent reviewer must verify that declaration
against actual inputs, not infer it from the resulting sign.

Each accepted update contains positive increasing `update_index`, matching
`voltage_update_index`, `current_update_index` and `heat_update_index`,
`state="accepted"`, `current_A` and `port_voltage_V` as complex component
pairs, nonnegative `face_heat_W` (integrated watts per face), and
`linear_true_relative_residual`. Initialization/reset placeholders must not
be submitted. Voltage must come from the accepted native electromagnetic
port unknown, not a reconstruction using the same reaction assembly. Match
the prescribed current to that update, not the next Simulink input sample.

The input also requires `artifact_identity` with `source_commit` (full SHA),
`source_dirty=false`, `wheel_sha256`, `loaded_mex_sha256`, `loaded_mex_path`
(the resolved absolute path from the acquisition process), and `plan_sha256`.
The CLI hashes the supplied wheel and actual module, verifies that the wheel
contains exactly one IH MEX with matching bytes, and compares all identity
fields before evaluation. The acquisition process must obtain the loaded
path from MATLAB, not assume it from the path argument. Raw input records
containing local paths stay in private operational storage; the result omits
that path. Source-to-wheel build receipt verification and fresh source hashes
remain the installed-wheel runner's responsibility; this evaluator does not
turn a caller-supplied commit into proof of build provenance.

Outputs contain per-update work and sign verdicts, plus input, plan and
evaluator hashes. The output path must be new; existing evidence is never
overwritten. The evaluator always reports `b1_closed=false` and
`release_accepted=false`: wheel/MEX/source identity, native extraction,
material declarations, official Toolkit lifecycle and independent FEM with
refinement are separate obligations. No physical solver execution has been
performed for this preparation. The unit-test numbers are synthetic algebra
checks, not numerical qualification of a device or release.

Remaining implementation: an acquisition driver for the final installed
wheel must export coherent accepted voltage/current/heat records for both
dedicated cases, with artifact receipts. The nonmagnetic case also needs
independent FEM reaction results and both BEM refinements before review can
close B1. B2 waits for the final wheel and fresh same-artifact runs.

## Installed MATLAB acquisition call contract

`recordReactionGPU.m` now implements the two-case native acquisition. The
installed-wheel runner supplies the original public temperature-BH sample's
`cfg` and `assembly`, after verifying its clean exact-source build receipt:

```matlab
paths = recordReactionGPU(cfg, assembly, wheelPath, fullSourceSHA, newOutputDir);
```

It keeps the supplied geometry/material/thermal configuration for the B-H
case, but explicitly sets the single supplied coil impedance to zero before
context creation. For the second case it additionally removes the B-H table,
sets every initial permeability to mu0 and saturation magnetization to zero,
and requires positive resistivity. These are declared input experiments, not
corrections to resulting voltage or heat. The default excitation is the same
synthetic 10 kA peak as the public sample; five updates are recorded. The runner
may explicitly supply `CurrentA`, `AmbientK` and bounded `Updates` options.

After each native update, the driver reads `outputs`, `electromagnetic_state`
and `element_impedance` from its own checked handle, requires the native step
count to match, and reads voltage directly from the declared port unknown.
It preserves each completed update JSON and the actual case input MAT before
proceeding. It never changes MATLAB paths, unloads modules or closes another
session. A different resident IH MEX is rejected; its owned handle is cleaned
up on success or failure. The two returned records paths are CLI inputs for
`reaction_gate.py`. Raw acquisition includes local paths and stays private.

No acquisition has been executed. The driver does not replace independent
FEM or two-level BEM refinement: its geometry is whatever verified public
configuration the runner supplies. Missing final wheel, formula approval and
compute admission still block numerical execution. The previously listed
acquisition implementation gap is addressed in source, pending execution and
review; all physical and release conclusions remain open.
