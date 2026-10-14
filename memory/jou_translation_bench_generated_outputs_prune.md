# JOU Translation Bench Generated Outputs Prune

Date: 2026-06-27

The `examples/jou_translation_bench/*.jou` fixtures are protected source
artifacts: they encode Cubit CAD/mesh intent and are useful for future
translator evaluation. The generated `examples/jou_translation_bench/out/*.py`
files are not fixtures; they are one run's LLM output.

The generated Python outputs were pruned rather than promoted to user-facing
docs; recover the old generated files from git history if a forensic comparison
is ever needed.

Future translation runs should write fresh output under an ignored or temporary
location, then summarize the verdict in docs or memory. Do not keep generated
LLM translation attempts as canonical example scripts unless they become a
tested translator fixture.
