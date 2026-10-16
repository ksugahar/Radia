# Eqnedit64 obsolete PyPI releases: owner procedure

Keep `3.1.4`. Deletion is performed by the project owner in the authenticated
PyPI web interface; this procedure and its checker perform no deletion or
credential handling. Do not delete the project or recreate existing artifacts.

## Prepare the exact list

Read `https://pypi.org/pypi/eqnedit64/json`, save the response privately with its
retrieval time and SHA-256, and enumerate actual `releases` keys. Do not fill
gaps in version numbering. Separate the retained release from obsolete releases
and record the number of distribution files for each. Review the concrete list
through both repository-management leads before the owner acts.

## Owner clicks

1. Sign in directly to PyPI with the owner's account and complete any requested
   authentication yourself.
2. Open <https://pypi.org/manage/project/eqnedit64/releases/>.
3. For each reviewed obsolete version, select its **Options**, then **Delete**
   for that release. Verify the version in the confirmation and follow the
   page's requested confirmation. Leave `3.1.4` and the project intact.
4. Return the completion result through repository management for verification.

These steps follow [PyPI's release management documentation](https://docs.pypi.org/project-management/storage-limits/#freeing-up-storage-on-an-existing-project).
Deletion is irreversible and filenames cannot be reused; see
[PyPI help](https://pypi.org/help/#file-name-reuse). Do not interpret a missing
Delete option or authentication failure as permission to bypass the interface.

## Read-only post-check

```powershell
python tools/eqnedit64/build/check_pypi_latest_only.py --keep 3.1.4 --report <private-result.json>
```

PASS requires both the production project JSON and uncached production
`pip index versions eqnedit64` to list only `3.1.4`, with retained files present.
FAIL means the observed state is not latest-only; UNVERIFIED means retrieval or
parsing failed. A pre-deletion FAIL is a baseline, not a failed deletion attempt.
This check does not certify GitHub Release deletion or remove local artifacts.
