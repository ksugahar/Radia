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
     "outcome": "adopted", "program": "kaken_oss",
     "source_paths": ["../submitted/application.pdf"],
     "outcome_basis": "award notice identifies the submitted project",
     "outcome_evidence": ["../results/award-notice.pdf"]},
    {"label": "own-draft",
     "paths": ["../draft/purpose.tex", "../draft/abilities.tex"],
     "pdf": "../draft/proposal.pdf", "program": "kaken_oss"}
  ]
}
```

Use `path` for one frozen text snapshot. Use `paths` for the ordered source
files of a live proposal; the lane joins them before running the checks. Name
exactly one of the two. Live sources avoid a stale extracted-text copy making
the regression suite pass after the proposal itself has changed.

An optional `pdf` names the compiled document. A page limit is a property of
the rendered page, and it is the only defect class that gets a proposal
returned before anyone reads it, so a document that has one is locked field by
field even while the check reports nothing.

Paths resolve against the manifest's directory. `program` selects the
program-specific checks (`generic`, `kaken_generic`, `kaken_oss`,
`kddi_digital`). Use `kaken_generic` for ordinary KAKENHI applications and
reserve `kaken_oss` for the current OSS-platform theme. `outcome` is
recorded but never scored — four measurements have found no relationship
between these checks and adoption, so a document's outcome is context for the
reader, not a target for the tool.

For documents labelled `adopted` or `rejected`, record provenance when it is
available. `source_paths` names the immutable submitted files from which the
text snapshot was extracted. `outcome_basis` explains the classification, and
`outcome_evidence` names award notices, review results, or other files that
support it. These fields stay outside the public repository, but the loader
validates every named file so a moved or guessed source cannot silently remain
in the corpus.

Documents extracted from Word or PDF should be converted to UTF-8 text once
and stored in `texts/`; the lane does not run Word or a PDF reader.

## Running

```powershell
$env:GRANT_WRITING_CORPUS = "<somewhere private>/manifest.json"
python validation_test/grant_writing/sweep.py
python validation_test/grant_writing/sweep.py --audit
python validation_test/grant_writing/sweep.py --compare-outcomes
python -m pytest validation_test/grant_writing
```

`sweep.py` prints one row per document and the pattern table across the whole
corpus. That table is the working surface: read every pattern, look at the
excerpt behind it, and decide whether a reader would agree the tool found a
real defect. Running every detector over every document at once found eight
false-positive families in a single pass, after several sessions of finding
one at a time.

`sweep.py --audit` reports, per check, how often it applied to a real proposal
and how often it said anything. A false positive announces itself; a check that
has quietly stopped working, or that was aimed at a genre this suite does not
serve, says nothing and looks exactly like a clean document. On the first run
eight checks were silent on all eight documents, and one of them
(`literature_gap_evidence_check`) turned out to be applicable to nothing at all.
Run it whenever the corpus grows.

The audit honors each manifest entry's `program`: KDDI-only checks are not
counted against a KAKENHI or generic document, and the KAKENHI OSS-platform
check is not counted against unrelated applications. `reported` is derived
only from explicit finding containers and known defect counts. Measurement
metadata such as `statement_count` is not itself a finding.

`sweep.py --compare-outcomes` reports detector prevalence separately for
recorded adopted and rejected documents, both across the corpus and for
ordinary KAKENHI entries. It is deliberately descriptive: programme, year,
panel, scientific maturity, and competition are uncontrolled, so the output
must not be used as an adoption score or a causal explanation.

`sweep.py --write-baseline` records the current counts. Only do that once the
counts have been adjudicated — the baseline's value is that every number in it
was read and judged.

## What the tests assert

- every baselined document is still in the manifest
- each source fingerprint matches the text whose findings were adjudicated
- each document's finding count matches the baseline
- no finding pattern appears that was absent when the corpus was adjudicated
- no field runs past its page allowance, and page usage per field is unchanged

A count that moves is either a fix worth re-baselining or a false positive
coming back. The test cannot tell which, and does not try: it forces the
question to be asked in the session that caused it.

Verified by deliberately disabling the form-instruction stripping, which
raised the adopted 科研費 proposal from 5 findings to 6 and introduced
`international_standing/no_named_counterpart` — the form's own 人権 boilerplate
mentions 国際共同研究, and without the stripping that reads as the applicant
claiming a foreign collaboration. Both tests named the document and the
pattern.
