# Documentation maintenance evidence

This directory contains migration ledgers, cleanup records, release evidence,
and other internal documentation-maintenance artifacts. They are retained for
traceability and regression work, but they are not public capability pages.

The supported user entry point is the Radia MCP server. Its tool descriptions,
schemas, workflow guidance, and knowledge resources form the canonical manual.
GitHub-facing content under `docs/` is the discovery and evidence layer: it
should show the engineering problem, what Radia can do, the expected result or
artifact, and the owning MCP tool family.

MATLAB, MEX, and Simulink assets currently serve implementation, integration,
validation, and parity work while their product role is evaluated. An `.mlx`
file must not be the sole public explanation because GitHub does not render it
as a readable notebook page; use Markdown or a saved-output `.ipynb` for the
GitHub-facing account.
