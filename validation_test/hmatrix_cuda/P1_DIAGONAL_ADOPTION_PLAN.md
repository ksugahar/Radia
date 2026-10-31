# P1 self-diagonal validation and promotion

Status: source-reviewed implementation; **not adopted**. The criterion proposal
is [p1_charge_diagonal.criterion.proposed.json](p1_charge_diagonal.criterion.proposed.json).
Its numbers are proposals only; the driver never uses them to declare a pass.

## Re-executable production entries

`validate_p1_charge_diagonal.py` executes the production CPU
`_ChargeGramHMatrix.entry(i,i)` and `charge_sigma()` for the 13 P1 charge modes
on six geometries. It obtains B/rules through the existing NG-owned charge
staging under a caller-owned TaskManager. It does not compute expected values
through the new diagonal arithmetic. FE DOFs (12) are distinct from charge rows
(13). Physical target weights exclude the monomial; the GPU entry applies it.

The six cases are base, a nonsymmetric affine shape, translation, rotation,
scale 2 and a reference-vertex 1/2 permutation. A reference permutation also
permutes the two linear modes; it is not an unchanged basis. Transformed cases
are charge-constructor diagnostics, not independently checked FE meshes or
complete solves. The original one-TET fixture and label contract must be
supplied explicitly with an externally pinned portable mesh-check receipt.
Missing inputs fail; no optional skip or mesh repair is implemented. Verify
publishable fixture provenance before adding a fixture to the repository.

CPU reference command (paths/pins are supplied by the owned execution frame):

```text
<absolute-python> validation_test/hmatrix_cuda/validate_p1_charge_diagonal.py
  --backend cpu --git <absolute-git>
  --source-root <frozen-independent-cpu-source> --cpu-source-sha <full-sha>
  --cpu-native <actual-cpu-module> --cpu-native-sha256 <sha256>
  --mesh <one-tet.vol> --contract <labels.json>
  --receipt <portable-receipt> --receipt-sha256 <sha256>
  --receipt-provenance <producer-provenance.json>
  --prior-oracle <frozen-original.json> --prior-oracle-sha256 <sha256>
  --out <new-cpu-result.json>
```

The original frozen reference remains attached to its own source and hash. A new
driver execution gets a new result and driver SHA; its replay comparison does
not relabel old evidence. The result records actual CPU module/receipt hashes,
NG version/provider path, checked input hashes and all 78 values/scales.

GPU mode consumes an externally pinned CPU-result JSON. It requires a separately
approved, hash-pinned manager admission (source/native/oracle SHA, device UUID
and absolute deadline), plus a fresh external resource frame. It invokes the
production `radia.gpu._native('_hdiv_cuda')._P1ChargeDiagonal` entry with
`order=1, backend='gpu'`; the same-source CPU/GPU loader guard stays enabled.
It never imports a module by file path to bypass that guard. Source qualification
may provide a clean same-source root; installed exact-wheel qualification must
use the accepted installation and release artifact/RECORD checks separately.

```text
<absolute-python> validation_test/hmatrix_cuda/validate_p1_charge_diagonal.py
  --backend gpu --git <absolute-git>
  --oracle <new-frozen-cpu-result.json> --oracle-sha256 <sha256>
  --admission <manager-frame.json> --admission-sha256 <sha256>
  --nvidia-smi <absolute-device-probe>
  --mesh <same-one-tet.vol> --contract <same-labels.json>
  --receipt <same-portable-receipt> --receipt-sha256 <sha256>
  --receipt-provenance <producer-provenance.json> --out <new-gpu-result.json>
```

No GPU execution is authorized by these commands or this document. No fallback
to CPU is provided. CPU/GPU import/provider compatibility must be qualified on
each actual target. Prior constant-only or released solver evidence does not
qualify this 13-mode constructor.

## Required checks and current deferrals

| Obligation | Entry/evidence | Promotion condition |
|---|---|---|
| Source behavior | C++ input/math controls and copied mutants | Guards only; not numerical adoption or a GPU surrogate |
| Independent implementation parity | This driver's CPU and GPU modes, 13 modes x 6 cases | Actual same-source GPU artifacts and separately approved limits; no pending/missing case can pass |
| Formulation independence | Proposed singular self-integral/analytic reference | Not yet implemented/qualified; CPU and GPU share primitive potentials, so parity alone does not establish it |
| Convergence | Proposed fixed-P1 outer-rule study, orders 3/5/6 | Not executed. Resolve reference, theory/rate and finest threshold, then user approval; do not count rule aliases 3/4 as distinct levels |
| Whole-method h/p convergence | Existing/future full HDiv workflow validation | Separate from this self-entry component; no claim from diagonal parity or a fake P2 run |
| Runtime lifetime | Actual native construct/download/close and failure controls | Execute on real candidate; host instantiation is insufficient |
| Release acceptance | Owning distribution, exact wheel/RECORD and adopted-case gates | Required before publication; source precheck is not a replacement |

`execution_complete=true` and exit 0 mean only that the requested campaign ran.
`numerically_accepted` stays false and `criterion_status` stays
`pending-user-approval`. A failed execution records an error and exits nonzero.
No `adoption_result` or acceptance record is emitted. Raw JSON stays in owned
evidence storage; no result is committed as a golden baseline by this change.

The current constructor-only quantity has no linear solve residual. Existing
adoption schema requires residual metadata. Do not invent a zero residual to
fit it. Resolve applicability with the schema owner/user, or link adoption to a
future actual production solve and its true residual, before recording adoption.

## Promotion map, only after qualification and user approval

1. Review the complete criterion (quantity/units, independent reference and
   shared components, named convergence series, theory/rate/finest thresholds,
   failure/plateau policy and residual applicability). Record actual user
   approval of its hash; recomputing a hash is not approval.
2. Execute the declared campaign on the exact native/wheel source, actual CPU
   and required GPU targets. Preserve original result JSON, artifact hashes,
   runtime/provider identity and all failures. Source changes require new
   execution/provenance; do not attach old results to a new SHA.
3. Create the coarsest production-entry regression only then, proposed path
   `tests/test_hdiv_p1_diagonal_surrogate.py`: checked one-TET P1 base case,
   all 13 modes, the same actual `_P1ChargeDiagonal` entry and approved stored
   baseline. It must not reimplement the kernel, substitute the five constant
   modes, or run the host-only arithmetic as GPU proof. Record its runtime.
4. Demonstrate that a temporary production source-moment omission and a target
   monomial omission fail that surrogate on actual rebuilt candidate kernels,
   and that restoring production passes. Existing host mutant kills are useful
   guards but are not those GPU surrogate mutation results. Requires separate
   build/GPU grants; do not update expected values to make mutants pass.
5. Only after those gates, add a proposed record under
   `validation_records/hdiv_p1_diagonal/<decision-id>/` with approved criterion,
   rationale, user approval identifier/hash, exact execution source/version,
   result/script/native/wheel hashes, reference independence/shared components,
   conditions/units, covered-content hashes and `protected_by` surrogate path.
   Use the actual repository schema/checker. Preserve stale/failed history.

No surrogate baseline, adoption block or adoption record is created now. The
remaining independent formulation reference, convergence, residual applicability,
user approval, native/GPU and exact-wheel steps are explicit deferrals, not
waivers. The present claim is reproducible CPU reference production-entry
execution plus a source-reviewed GPU campaign path, with adoption still pending.
