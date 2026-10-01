# CLN implementation policy

Radia and radia-mcp must not contain implementations that infringe CLN patents.
This applies to Python, native code, MATLAB/Simulink, MCP tools, examples,
notebooks, and distributed artifacts. Renaming, wrapping, or moving an
implementation within these products does not satisfy this policy.

## Removal scope

The current removal program excludes CLN reduction implementations from these
products and their published Git history. Identify implementations from their
algorithms, including alternating field recurrences; do not rely on filenames
or API names alone. Preserve unrelated changes when rewriting history.

Paper citations and bibliographic records may remain. Independent Foster modal
models, PRIMA projection, snapshot POD, ordinary Lanczos/Arnoldi methods, and
classical scalar rational/continued-fraction algebra are not removed merely
because they previously shared a file or were described using CLN terminology.
DtN continued fractions represent exact/high-order non-reflecting boundaries
and remain supported. Retired reduction APIs are not production/comparison
routes. Select reduced orders by measured operating-band error.

## Review and evidence

Review implementation bodies and callers, preserve independent numerical
methods and citations, and verify affected behavior. Record the exact reviewed
revision and unresolved findings. Reciprocal review is required for the current
removal and history-rewrite campaign.

Algorithm names, literature citations, tests, or code removal do not establish
patent clearance. Do not claim that an implementation is non-infringing solely
from those facts; unresolved patent-scope questions require qualified review.
This policy records the project restriction, not a legal conclusion.
