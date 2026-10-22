# Frozen calculation regression records

These JSON bytes are copied from the Git blobs identified in `manifest.json`.
Their original execution revision was not recorded. Reproduction at the exact
`reference_source_sha` with the stated conditions matched all eighteen records;
that verified baseline is distinct from their unknown original execution.

The integration test compares calculation outputs only. Independent analytic
acceptance remains in `validation_test/radia_mcp/test_netgen_vol_analytic_examples.py`
with nightly classification. Run that module explicitly with the selected MCP
source on PYTHONPATH; this change does not claim an existing scheduled job runs it.
Historical documentation records stay at their original paths.

To propose a reference update, run
`python -m validation_test.radia_mcp.regenerate_netgen_vol_regression --out <new-directory> --reason <reason>`.
It creates per-case `baseline.json` (the original nested record schema) and
common-helper `proposal.json` files plus a top-level index. It records the clean
generating HEAD, Python version, runtime metadata, original fixture SHA/blob,
conditions and differences. Dirty source is rejected; it never overwrites these
records. Obtain user approval of changed references
and independent review before replacing approved fixture bytes or the manifest.

The examples require an installed `radia-mcp`, or its selected source on
PYTHONPATH. For direct CLI runs, always pass `--out <job-owned-output.json>`;
the default writes beside the example source and must not be used in automation.
For example: `python validation_test/radia_mcp/validation_p1_surface_triangle_element_summary.py --out <job-owned-output.json>`.
Unchanged JSON values retain the approved bytes in regeneration proposals,
including whitespace and signed-zero spelling; changed proposals include the
raw approved-to-proposed textual diff and separate semantic/byte flags.

Comparisons use the common helper with `atol=rtol=0`. Numerical leaves become
dotted names only in the in-memory adapter; structure and nonnumerical fields
are exact conditions. Historical runtime is unknown (the adapter's zero runtime
is explicitly an unmeasured placeholder), never compared. Fixture JSON values,
bytes and manifest hashes remain unchanged.
