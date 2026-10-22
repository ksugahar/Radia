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
`python validation_test/radia_mcp/regenerate_netgen_vol_regression.py --out <new-directory> --reason <reason>`.
It creates a separate proposal, producing source SHA, conditions and differences;
it never overwrites these records. Obtain user approval of changed references
and independent review before replacing approved fixture bytes or the manifest.
