# Numerical CI consumer plan (S1)

`python tools/ci_numerical_plan.py --base <commit> --head <commit> --out <json>`
reads immutable Git blobs and reports separate native, AMS, Python-consumer and
heavy-case selections. Base and head dependencies are combined, including package
initializers, conftest files and statically identified Python subprocess scripts.
Each case records dependency blob identities and unresolved inputs.

Unknown changed inputs select all cases; unresolved case closures select that
case. Native/build/control changes and non-PR events retain conservative execution.
A proven Gmsh-only dependency change selects post-processing consumers rather than
classifying Gmsh as a solver dependency. Mixed solver dependency changes remain
conservative. Numerical criteria and registered execution routes are unchanged.

This is a diagnostic static graph, not proof of execution or a complete Python
program analysis. Dynamic loaders, runtime commands and missing imports require
declarations or conservative fallback. The CLI result has `qualification=false`
and `execution_wired=false`: it cannot authorize skipping any existing gate.
S2 must separately review workflow wiring, per-case fingerprints and real-source
replays before removing existing forced execution. Execution must propagate the
admitted thread budget to child processes explicitly.
