# Documentation maintenance evidence

This directory contains migration ledgers, cleanup records, release evidence,
and other internal documentation-maintenance artifacts. They are retained for
traceability and regression work, but they are not public capability pages.

The supported user entry point is the Radia MCP server. Its tool descriptions,
schemas, workflow guidance, and knowledge resources form the canonical manual.
GitHub-facing content under `docs/` is the discovery and evidence layer: it
should show the engineering problem, what Radia can do, the expected result or
artifact, and the owning MCP tool family.

Radia is MCP-native. LLM agents drive the Python implementation through Radia
MCP, while the formal Simulink UI requires MathWorks' official MATLAB MCP
Server. Missing MCP connectivity is an explicit failure, not a reason to route
users to Python or private MATLAB automation.

Masked Simulink blocks are Radia's formal human-facing UI. MATLAB and MEX assets
support that UI, integration, validation, and parity work; a separate standalone
MATLAB edition remains undecided. MATLAB `.m` and `.mlx` files are implementation
assets under `matlab/`, not public `docs/` content. Use Markdown or an executed
Python `.ipynb` for the GitHub-facing capability account.
