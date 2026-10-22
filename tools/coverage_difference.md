# Manual coverage difference report

`coverage_difference.py` is a stdlib-only comparator over explicit JSON exports.
It does not capture coverage, run either command, check out source, modify a
baseline, or approve moving a validation. It is not a CI gate. Its synthetic
contracts test set arithmetic and failure handling, not physical correctness.

```text
python tools/coverage_difference.py --heavy heavy.json --light light.json --out new-report.json
```

The output must not already exist. `tests/test_coverage_difference.py` contains
a complete hand-authored minimal export pair. Each input has these fields:

| Field | Content |
| --- | --- |
| `schema`, `role` | `radia.coverage-evidence.v1`; `heavy` or `light` |
| `identity.source_sha` | Full generating source SHA; provenance, not ancestry |
| `identity.artifacts` | Artifact name to SHA256 map; include relevant Python/native identities |
| `identity.backend`, `route` | Explicit backend and identical entry point, formulation, element, order, nonlinear/boundary/constitutive branches |
| `identity.production_scope` | Identical normalized repository-relative file list |
| `identity.instrumentation` | Tool/version, branch tracing, exclusions and subprocess requirement/sentinel result |
| `run` | Recorded argv, model/mesh/hardware/thread/affinity conditions, runtime metadata, original numerical gate status and evidence SHA256 |
| `files` | Per-file content hash, executed positive lines, executed branch arcs (negative entry/exit lines permitted) |
| `native` | Named entry-point invocation count, declared artifact and configuration; empty for Python-only evidence |

The comparator requires matching source/artifact/route/scope identities and
instrumentation versions/exclusions. Different file hashes fail loudly. Missing
file capture, branch tracing, required subprocess sentinel or successful original
numerical gate remains UNKNOWN. Empty production capture is insufficient.
Line and arc intersections and both differences are reported separately: equal
line sets cannot hide different branch arcs. No coverage ratio or elapsed-time
threshold is used. The 30-second reduced-model budget belongs to model selection,
not this tool's pass/fail logic. Numerical gate statuses are copied unchanged.

Native entry hits never establish internal C++/CUDA/JIT branch coverage. This
first slice does not discover enclosing Python functions, attest artifact hashes,
inspect model semantics, or verify capture subprocess identities. Those are
producer/domain-review responsibilities; no missing attestation is invented.
The JSON inputs are declared evidence, not trusted execution receipts. Invalid
schema or incompatible identities raise an error. A well-formed report returns
normally even with UNKNOWNs: command exit status is not a qualification verdict.

The owner must review each heavy-only branch, retain the needed independent
validation, and attach evidence before recommending a reduced regression lock.
No report automatically downgrades a test or creates a lifecycle registry.
Capture integration and optional native counters are separate future work.
