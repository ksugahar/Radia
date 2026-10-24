# Numerical CI consumer plan

`python tools/ci_numerical_plan.py --base <commit> --head <commit> --out <json>`
reads immutable Git blobs and reports separate native, AMS, Python-consumer and
heavy-case selections. Base and head dependencies are combined, including package
initializers, conftest files and statically identified Python subprocess scripts.
Batch blob loading and a shared parsed graph avoid per-case subprocess/parsing
cost. Each case records dependency blob identities, components and unresolved inputs.

Unknown changed inputs select all cases; unresolved case closures select that
case. Native/build/control changes and non-PR events retain conservative execution.
A proven Gmsh-only dependency change selects post-processing consumers rather than
classifying Gmsh as a solver dependency. Mixed solver dependency changes remain
conservative. Numerical criteria and registered execution routes are unchanged.
Gmsh remains post-processing support; checked `.vol` files remain solver inputs.

`tools/ci_numerical_inputs.json` declares finite dynamic loaders, CLI inputs and
fixtures. A declaration resolves only the named issue in the named source and is
sealed to its Git-blob SHA-256. Runtime-generated coil code is additionally scoped
to the exact test that writes it. Stale declarations retain conservative execution;
they never silently extend to new callers. Explicit consumer edges retain the
viewer/native-field preflight checks without assigning the viewer to AMS/heavy.
Docstrings and generated template fragments are data, not executable Python.
Actual `exec`, file loaders and unbounded subprocesses remain unresolved unless
explicitly declared. New declarations require review of their finite inputs.

The plan remains a static graph, not execution evidence: `qualification=false`
and `execution_wired=false` describe the JSON alone. Workflow consumers verify
its complete recomputed selection, runner event/base, HEAD and current input blob
hashes before selecting executable pytest files. Git clean filters permit normal
Windows line endings without trusting stat-cache or skip-worktree flags. Empty
execution is an error, not a successful gate.

SparseSolv execution uses separate AMS, Python and heavy conditions. Selected
native consumers still build current native dependencies and verify provenance.
Gmsh installs alone from its dedicated official dev index with an installer report;
other dependencies use their existing separate install. A plan and results are
uploaded before runner scratch cleanup. Thread limits, including child
`NGS_NUM_THREADS`, remain explicit and unchanged.

Heavy source receipts fingerprint selected case dependencies, native build inputs,
profile contracts and relevant numerical/toolchain packages. Gmsh versions matter
for Gmsh consumers; unrelated venv packages do not invalidate a proved closure.
Unknown closures retain the full environment fallback. Source receipt reuse still
requires trusted successful main ancestry, exact hashes, non-skipped JUnit and the
existing age limit. It never certifies a newly built artifact. PRs no longer force
heavy validation unconditionally; explicit manual forcing remains available.

Conservative selections cannot use a scoped source fingerprint. Unclassified
inputs, global/control selection or unknown selected heavy closures fingerprint
all tracked files and the entire environment, including the workflow or runtime
contract that caused fallback. Proven closures retain their static third-party
imports; CMake and Ninja are bound to the configured executable hashes and version
outputs. Unresolved aliased subprocess/shell calls remain conservative. Joined
literal documentation paths read by a case are inputs, not docs-only exemptions;
an unresolvable documentation read marks the case unresolved. This is not a
general proof of arbitrary filesystem access.


SparseSolv does not reuse a PR's whole-tree job success on main: after per-case
selection that success no longer proves the complete lane ran. Its parent receipt
explicitly records `case-plan-not-full-lane`. The full main plan must execute;
only the separately verified, scoped heavy source receipt may avoid duplicate
heavy work. Other distributions' existing tree-reuse rules are unchanged.


Execution targets come from the current tier registry plus the explicit current
Gmsh-consumer supplement; retired paths are recorded
separately rather than passed to pytest after deletion. Both Git snapshots still
contribute dependencies of retained cases. This does not authorize retirement:
the existing route-retention and review gates still govern registry removals.
