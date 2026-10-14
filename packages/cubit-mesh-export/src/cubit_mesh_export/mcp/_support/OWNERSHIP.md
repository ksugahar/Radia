# Cubit-owned runtime support

These are intentionally independent implementations, not synchronized copies
or generated files. Cubit maintains and releases them under its own contracts.
Radia MCP is neither the source of truth nor a build/runtime dependency for them.
No cross-product byte-equality test or automatic synchronization is intended.

When fixing a defect in inherited logic, inspect the other product for the same
defect and notify its owner; each owner decides and tests its own correction.
Intentional differences in environment variables, provenance, state directories,
example providers and mesh contracts must remain owner-specific. Cubit's tests
cover these contracts independently, including operation without Radia installed.

Support derived from the former shared runtime retains the copyright and full
BSD-3-Clause terms in LICENSE-BSD-3-Clause.txt. The containing exporter uses MIT;
this does not relicense the inherited support or discard its notices.

Cubit controls use CUBIT_MCP_* environment variables. RADIA_MCP_* variables no
longer configure this product. Its state directory is independent; old Radia
failure logs, caches and learned example indexes are not copied or deleted.
Operators may explicitly copy needed records into CUBIT_MCP_STATE_DIR after
checking their format and ownership. There is no automatic history migration.
