# Grant-writing corpus regression

The fast suite in `packages/radia-mcp/tests/test_grant_writing.py` locks each
check against small synthetic inputs. This lane locks the whole detector set
against **real proposals**, because that is where the false positives were.

## Why the corpus is not in this repository

Real proposals belong to their authors. Several in the working corpus are
colleagues' work, and this repository is public. Neither the documents nor
their paths are committed: the lane reads a manifest named by the
`GRANT_WRITING_CORPUS` environment variable and skips entirely when it is
unset, which is what happens in CI and in a fresh clone.

## Setting up a corpus

Put the manifest, the text files and the baseline together outside the tree:

```
<somewhere private>/
  manifest.json
  baseline.json          # written by sweep.py --write-baseline
  texts/
    adopted-example.txt
    rejected-example.txt
```

`manifest.json`:

```json
{
  "documents": [
    {"label": "adopted-example", "path": "texts/adopted-example.txt",
     "outcome": "adopted", "program": "kaken_oss"}
  ]
}
```

Paths resolve against the manifest's directory. `program` selects the
program-specific checks (`generic`, `kaken_oss`, `kddi_digital`). `outcome` is
recorded but never scored — four measurements have found no relationship
between these checks and adoption, so a document's outcome is context for the
reader, not a target for the tool.

Documents extracted from Word or PDF should be converted to UTF-8 text once
and stored in `texts/`; the lane does not run Word or a PDF reader.

## Running

```powershell
$env:GRANT_WRITING_CORPUS = "<somewhere private>/manifest.json"
python validation_test/grant_writing/sweep.py
python -m pytest validation_test/grant_writing
```

`sweep.py` prints one row per document and the pattern table across the whole
corpus. That table is the working surface: read every pattern, look at the
excerpt behind it, and decide whether a reader would agree the tool found a
real defect. Running every detector over every document at once found eight
false-positive families in a single pass, after several sessions of finding
one at a time.

`sweep.py --write-baseline` records the current counts. Only do that once the
counts have been adjudicated — the baseline's value is that every number in it
was read and judged.

## What the tests assert

- every baselined document is still in the manifest
- each document's finding count matches the baseline
- no finding pattern appears that was absent when the corpus was adjudicated

A count that moves is either a fix worth re-baselining or a false positive
coming back. The test cannot tell which, and does not try: it forces the
question to be asked in the session that caused it.

Verified by deliberately disabling the form-instruction stripping, which
raised the adopted 科研費 proposal from 5 findings to 6 and introduced
`international_standing/no_named_counterpart` — the form's own 人権 boilerplate
mentions 国際共同研究, and without the stripping that reads as the applicant
claiming a foreign collaboration. Both tests named the document and the
pattern.
