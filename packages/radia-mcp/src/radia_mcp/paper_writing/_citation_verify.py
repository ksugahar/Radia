"""Citation-verification workflow for radia_mcp.paper_writing.

Added 2026-05-26 (radia-mcp v0.91.0).  Enforces the lab policy:

  **ALWAYS use the bundled canonical references.bib AND back-check every
  citation via Crossref / Semantic Scholar / arXiv before inserting it.**

Why this matters: the most common AI-assisted-writing failure mode
is *citation hallucination*  -- the AI invents a plausible-looking
DOI / author / year that does NOT correspond to any real paper.
Reviewers catch this 100% of the time (the DOI doesn't resolve, or
resolves to a different paper than the one cited), and one
hallucinated citation taints the entire submission.

This module composes the EXISTING radia-mcp tools into a hard
verify-first workflow:

  1. Read canonical references.bib (the SINGLE SOURCE OF TRUTH).
  2. Search Crossref / Semantic Scholar / arXiv for the claim.
  3. Match candidates against references.bib (already-cited? skip).
  4. Verify the top candidate exists via DOI resolve.
  5. Generate a ready-to-paste BibTeX entry only AFTER verification.
  6. Recommend whether to insert into references.bib (or report that
     the user is already citing it).

If steps 1-4 cannot complete cleanly, the tool returns an explicit error or
no-candidate verdict rather than a plausible-looking unverified citation.

Cross-references:
  * paper_writing_resolve_doi (Crossref backbone)
  * paper_writing_doi_to_bibtex (DOI -> BibTeX entry generator)
  * paper_writing_semantic_scholar_lookup (DOI/arXiv/etc resolver)
  * paper_writing_arxiv_search (free-text arXiv search)
  * paper_writing_lint_reference_format (post-insertion hygiene)
  * paper_writing_check_citation_usage (.tex cross-validation)
"""

from __future__ import annotations

import difflib
import os
import pathlib
import re
import unicodedata
from typing import Optional

from ..bibliography._bibparse import parse_bib


# ============================================================
# Behavioral knowledge: the workflow recipe
# ============================================================


_CITATION_VERIFICATION_RECIPE = r"""
# Lab POLICY: citation verification workflow (radia_mcp.paper_writing)

This recipe is the SINGLE BEHAVIORAL RULE for inserting any
citation into a lab paper:

  **NEVER invent a citation. ALWAYS use the bundled canonical
  references.bib and ALWAYS verify the source via search FIRST.**

## Why this is enforced

The most common AI-assisted-writing failure is HALLUCINATED
CITATIONS -- a plausible-looking DOI / author / year that does NOT
correspond to a real paper.  Reviewers catch every single one
(typical workflow: reviewer clicks the DOI, it doesn't resolve,
they conclude the rest of the paper is suspect, manuscript is
rejected).

One hallucinated citation taints the whole submission.  The cost
of preventing the failure is one extra search-and-verify per
citation (~5 seconds of human time, 0 seconds of AI time).

## The mandatory 6-step workflow

### Step 1: READ the canonical references.bib first

```python
from radia_mcp.bibliography.plans.T14_canonical import CANONICAL
from radia_mcp.paper_writing import paper_writing_lint_reference_format
report = paper_writing_lint_reference_format(str(CANONICAL))
# report contains: every existing entry's citation_key, DOI, title, author, year
```

The .bib is the **single source of truth** for what the paper
already cites.  Before suggesting ANY new citation, verify whether
the user is already citing the source.  If yes, REUSE the existing
citation_key; do NOT introduce a duplicate entry.

### Step 2: SEARCH for grounding via Crossref / Semantic Scholar /
arXiv

When the claim references a specific result / method / dataset,
search to find the canonical paper:

```python
from radia_mcp.paper_writing import (
    paper_writing_semantic_scholar_lookup,
    paper_writing_arxiv_search,
    paper_writing_resolve_doi,
)

# If you have a DOI guess:
meta = paper_writing_resolve_doi("10.1109/TAP.2006.880747")
# meta contains the Crossref-authoritative title, authors, year, venue

# If you have an arXiv ID:
lookup = paper_writing_semantic_scholar_lookup("arXiv:2603.17339")
# lookup['metadata'] has authoritative title, abstract, externalIds

# Free-text search by title / topic:
hits = paper_writing_arxiv_search(
    "Sommerfeld closed form impedance plane", max_results=5
)
# hits['papers'] is a ranked list with arxiv_id, title, authors, year
```

**Never make up a DOI** like "10.1109/TAP.20XX.XXXXXXX".  Either
resolve a candidate DOI via Crossref, OR find the paper via
arXiv search, OR mark the citation as "TODO: verify with user".

### Step 3: MATCH candidates against references.bib

```python
from radia_mcp.paper_writing import paper_writing_verify_citation

result = paper_writing_verify_citation(
    claim="Koh-Yook derived an exact closed form for impedance plane",
    candidate_doi="10.1109/TAP.2006.880747",
)
# Returns:
#   "already_in_bib": bool          -- if yes, reuse the existing key
#   "matching_key":    str|None     -- which key in bib (if matched)
#   "candidates":      list[dict]   -- search-verified candidates
#   "suggested_bibtex": str|None    -- ready-to-paste entry (if new)
#   "verification_method": str      -- "crossref" | "arxiv" | "s2"
#   "verdict":         "found_in_bib" | "ready_to_insert" |
#                       "needs_disambiguation" | "no_candidate_found"
#   "advice":          str          -- next step
```

### Step 4: VERIFY the DOI resolves (defensive double-check)

```python
# Even if the BibTeX was auto-generated, double-check the DOI:
v = paper_writing_resolve_doi(suggested_bibtex_doi)
assert v["ok"], f"Suggested DOI does not resolve: {suggested_bibtex_doi}"
```

If the DOI doesn't resolve, the BibTeX is unreliable -- discard
the candidate and search again.

### Step 5: INSERT into references.bib (only if verdict ==
"ready_to_insert")

Append the human-reviewed entry to the canonical `.bib`; never create a local
manuscript copy and never overwrite an existing entry.

### Step 6: LINT after insertion

```python
post_check = paper_writing_lint_reference_format(bib_path)
# Verify no duplicate keys, no missing required fields, no DOI typos
```

Then run `bibliography_make_bbl(tex_path)` for the manuscript. Submit the
generated `.bbl`, not a copied `.bib` database.

## When verification cannot complete

If after Steps 2-3 NO candidate has a resolvable DOI:

  * **DO NOT** generate a fake / placeholder DOI.
  * **DO NOT** invent author names or years.
  * **DO** mark the citation in the manuscript with a clearly
    visible TODO:
        ``\\cite{TODO: verify -- claim about Sommerfeld closed form}``
  * **DO** report this to the user explicitly:
        "I could not verify a citation for this claim.  Please
         provide the DOI or paper title."

## What this prevents

  * Hallucinated DOIs that don't resolve (catches 100% of these).
  * Wrong-author citations (the Crossref metadata is authoritative).
  * Mixed-up year/venue citations (e.g. confusing a 2006 paper
    with a 2007 conference version).
  * Duplicate citations in references.bib (different keys pointing
    to the same DOI).
  * Cite-but-not-in-bib (\\cite{xyz} where xyz isn't in the bib).
  * Bib-but-not-cited (entries in bib that never appear in \\cite).
    -- the existing paper_writing_check_citation_usage tool catches
    this on the .tex/.bib pair.

## Quick-reference: which tool for which job

  paper_writing_lint_reference_format(bib)
        Bib hygiene: missing fields, duplicate keys, DOI format.

  paper_writing_resolve_doi(doi)
        Crossref API lookup; authoritative metadata.

  paper_writing_doi_to_bibtex(doi)
        Generate BibTeX entry from DOI via Crossref.

  paper_writing_arxiv_search(query, categories="")
        Free-text search arXiv; EM-tuned default categories.

  paper_writing_semantic_scholar_lookup(paper_id)
        Resolve DOI / arXiv ID / S2 ID into metadata + abstract.

  paper_writing_semantic_scholar_references(paper_id)
        List of papers CITED BY the given paper -- find canonical
        predecessors.

  paper_writing_semantic_scholar_citations(paper_id)
        List of papers CITING the given paper -- find recent extensions.

  paper_writing_verify_citation(claim, bib_path, ...)
        The composite workflow: SEARCH + MATCH + VERIFY + SUGGEST.
        **Use this for any new citation.**

  paper_writing_check_citation_usage(tex, bib)
        Post-write check: every \\cite resolves; no orphan entries.

  bibliography_make_bbl(tex)
        Export only the cited canonical entries as the submission .bbl.

  paper_writing_fetch_and_cite(doi, dest_dir)
        Download PDF + generate BibTeX (IEEE / Crossref).
"""


def paper_writing_citation_workflow_recipe() -> str:
    """Return the lab's mandatory citation-verification workflow recipe.

    This is the BEHAVIORAL RULE for inserting any citation into a
    lab paper: NEVER invent, ALWAYS verify via Crossref / S2 /
    arXiv FIRST, ALWAYS check against canonical references.bib.

    Read this BEFORE generating any \\cite{} or BibTeX entry.
    """
    return _CITATION_VERIFICATION_RECIPE


# ============================================================
# Helpers
# ============================================================


def _canonical_bib() -> str:
    """The lab keeps one bibliography; papers do not copy it."""
    from ..bibliography.plans.T14_canonical import CANONICAL
    return str(CANONICAL)


def _parse_bib_lightweight(bib_path: str) -> dict:
    """Parse references.bib with the package's balanced-brace parser.

    Returns dict keyed by citation_key with sub-dict of fields.
    ``@string``, ``@preamble``, and ``@comment`` blocks are deliberately
    excluded because they are not citation targets.
    """
    bib_path = bib_path or _canonical_bib()
    text = pathlib.Path(bib_path).read_text(encoding="utf-8", errors="strict")
    return {
        entry.key: {"type": entry.kind, **entry.fields}
        for entry in parse_bib(text)
        if entry.key and not entry.kind.startswith("@")
    }


def _normalize_doi(value: str) -> str:
    """Return a resolver-independent DOI token for equality checks."""
    value = value.strip().lower()
    value = re.sub(r"^doi\s*:\s*", "", value)
    value = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", value)
    value = re.sub(r"^doi\.org/", "", value)
    return value.rstrip(".,; ")


def _doi_already_cited(bib_entries: dict, target_doi: str) -> Optional[str]:
    """Return the citation_key in bib_entries whose DOI matches
    target_doi (case-insensitive), or None if no match.
    """
    if not target_doi:
        return None
    t = _normalize_doi(target_doi)
    for key, fields in bib_entries.items():
        bib_doi = _normalize_doi(fields.get("doi", ""))
        if bib_doi == t:
            return key
    return None


def _title_already_cited(bib_entries: dict, target_title: str) -> Optional[str]:
    """Return the first citation_key whose title approximately matches
    target_title (case-insensitive, normalised whitespace).
    """
    if not target_title:
        return None
    def _norm(s):
        s = unicodedata.normalize("NFKC", s).casefold()
        s = "".join(
            char if char.isspace() or unicodedata.category(char)[0] in {"L", "N"}
            else " "
            for char in s
        )
        return re.sub(r"\s+", " ", s).strip()

    t_norm = _norm(target_title)
    if len(t_norm) < 10:
        return None  # too short, would match many things
    for key, fields in bib_entries.items():
        bib_title = _norm(fields.get("title", ""))
        if not bib_title:
            continue
        target_tokens = set(t_norm.split())
        bib_tokens = set(bib_title.split())
        token_overlap = len(target_tokens & bib_tokens) / max(len(target_tokens), 1)
        if bib_title == t_norm or (
            token_overlap >= 0.8
            and difflib.SequenceMatcher(None, bib_title, t_norm).ratio() >= 0.92
        ):
            return key
    return None


# ============================================================
# paper_writing_verify_citation -- the composite workflow tool
# ============================================================


def paper_writing_verify_citation(
    claim: str,
    bib_path: str = "",
    candidate_doi: str = "",
    candidate_arxiv_id: str = "",
    candidate_title: str = "",
    search_arxiv_if_no_doi: bool = True,
    arxiv_max_results: int = 5,
) -> dict:
    """Verify a citation BEFORE inserting it into the paper.

    Lab POLICY: this is the mandatory workflow for any new
    \\cite{} or BibTeX entry.  Composes the existing search +
    bib-lint tools into a single pass-fail check.

    Args:
        claim: free-text description of WHAT is being cited
            (used for search if no DOI / arXiv ID).  Example:
            "Koh-Yook derived an exact closed form for the
             impedance plane Sommerfeld integral."
        bib_path: optional external bibliography override. The bundled
            canonical references.bib is the default single source of truth.
        candidate_doi: known DOI candidate (skip search step).
        candidate_arxiv_id: known arXiv ID candidate.
        candidate_title: known paper title.  Used both for
            duplicate detection in .bib AND for arXiv search if
            no DOI / arXiv ID was provided.
        search_arxiv_if_no_doi: True (default) to fall back to
            arXiv keyword search when no DOI / arXiv ID is given.
        arxiv_max_results: cap on arXiv-search candidates.

    Returns:
        dict with:
          "verdict": one of:
              "found_in_bib"          -- already cited; reuse existing key
              "ready_to_insert"       -- verified, suggested_bibtex valid
              "needs_disambiguation"  -- multiple candidates, user pick
              "no_candidate_found"    -- search failed; mark as TODO
              "error"                 -- bib read failed, etc.
          "matching_key": str|None   -- existing bib key (if found_in_bib)
          "candidates": list[dict]    -- search-verified candidates
          "suggested_bibtex": str|None
          "verification_method": str  -- "doi-direct" | "arxiv-search"
                                          | "title-match" | "bib-only"
          "advice": str               -- actionable next step

    Lab usage:
        Before inserting ANY \\cite{xxx}:
        1. Call paper_writing_verify_citation(claim, bib_path, ...).
        2. If verdict == "found_in_bib": use matching_key, no edit.
        3. If verdict == "ready_to_insert": append suggested_bibtex
           to references.bib, then \\cite{citation_key}.
        4. If verdict == "needs_disambiguation": ask user to pick.
        5. If verdict == "no_candidate_found": mark TODO; ask user.
        6. If verdict == "error": investigate; do NOT proceed.
    """
    # --- Step 1: read the canonical bibliography (or explicit override) ---
    bib_path = bib_path or _canonical_bib()
    if not os.path.exists(bib_path):
        return {
            "verdict": "error",
            "advice": (f"bib_path not found: {bib_path}.  Ask the user "
                       f"to create references.bib or supply the correct path."),
            "candidates": [],
            "matching_key": None,
            "suggested_bibtex": None,
            "verification_method": None,
        }
    try:
        bib_entries = _parse_bib_lightweight(bib_path)
    except (OSError, UnicodeError, ValueError) as exc:
        return {
            "verdict": "error",
            "advice": f"failed to read references.bib: {exc}",
            "candidates": [],
            "matching_key": None,
            "suggested_bibtex": None,
            "verification_method": None,
        }

    # --- Step 2a: DOI fast-path (if user gave us one) ---
    if candidate_doi:
        matching = _doi_already_cited(bib_entries, candidate_doi)
        if matching:
            return {
                "verdict": "found_in_bib",
                "matching_key": matching,
                "candidates": [],
                "suggested_bibtex": None,
                "verification_method": "doi-direct",
                "advice": (f"Use existing \\cite{{{matching}}}.  Do NOT "
                            f"introduce a duplicate entry."),
            }
        # Verify DOI resolves, then generate BibTeX
        from .paper_download import (
            paper_writing_resolve_doi,
            paper_writing_doi_to_bibtex,
        )
        v = paper_writing_resolve_doi(candidate_doi)
        if not v.get("ok") and v.get("temporary_failure"):
            return {
                "verdict": "error",
                "matching_key": None,
                "candidates": [],
                "suggested_bibtex": None,
                "verification_method": "doi-direct",
                "advice": (
                    f"DOI verification is temporarily unavailable: {v.get('error')}. "
                    "Do not infer that the citation is fabricated; retry later."
                ),
            }
        if not v.get("ok"):
            return {
                "verdict": "no_candidate_found",
                "matching_key": None,
                "candidates": [],
                "suggested_bibtex": None,
                "verification_method": "doi-direct",
                "advice": (f"DOI {candidate_doi!r} was not found via "
                            f"Crossref.  It may be incorrect or absent from "
                            f"that registry.  Search for the paper via "
                            f"paper_writing_arxiv_search(claim) or ask "
                            f"the user for the correct DOI.  DO NOT insert "
                            f"this DOI into the bib."),
            }
        bib_result = paper_writing_doi_to_bibtex(candidate_doi)
        if not bib_result.get("ok"):
            return {
                "verdict": "error",
                "matching_key": None,
                "candidates": [{k: value for k, value in v.items() if k != "ok"}],
                "suggested_bibtex": None,
                "verification_method": "doi-direct",
                "advice": (f"DOI resolved but BibTeX generation failed: "
                            f"{bib_result.get('error')}"),
            }
        meta = {k: value for k, value in v.items() if k != "ok"}
        return {
            "verdict": "ready_to_insert",
            "matching_key": None,
            "candidates": [meta],
            "suggested_bibtex": bib_result["bibtex"],
            "verification_method": "doi-direct",
            "advice": (f"Verified.  Append the suggested_bibtex to "
                        f"{os.path.basename(bib_path)}, then use "
                        f"\\cite{{{bib_result['citation_key']}}}."),
        }

    # --- Step 2b: arXiv-ID fast-path ---
    if candidate_arxiv_id:
        from ._arxiv_source import (
            paper_writing_semantic_scholar_lookup,
            _normalize_arxiv_id,
        )
        aid = _normalize_arxiv_id(candidate_arxiv_id)
        # Try DOI lookup via S2 first (gives us a DOI to feed Crossref)
        lookup = paper_writing_semantic_scholar_lookup(f"arXiv:{aid}")
        if "error" not in lookup:
            ext = lookup.get("metadata", {}).get("externalIds", {})
            doi = ext.get("DOI") or ""
            if doi:
                return paper_writing_verify_citation(
                    claim=claim, bib_path=bib_path, candidate_doi=doi,
                )
            # No DOI but we have S2 metadata -- check title match
            title = lookup.get("metadata", {}).get("title", "")
            matching = _title_already_cited(bib_entries, title)
            if matching:
                return {
                    "verdict": "found_in_bib",
                    "matching_key": matching,
                    "candidates": [lookup.get("metadata", {})],
                    "suggested_bibtex": None,
                    "verification_method": "arxiv-title-match",
                    "advice": (f"arXiv:{aid} title matches existing bib "
                               f"entry \\cite{{{matching}}}.  Reuse."),
                }
            # No DOI, no title match -- emit a basic arXiv BibTeX
            authors = ", ".join(
                a.get("name", "?")
                for a in lookup.get("metadata", {}).get("authors", [])[:4]
            )
            year = lookup.get("metadata", {}).get("year", "????")
            key = (f"arxiv_{aid.replace('/', '_').replace('.', '_')}")
            bibtex = (
                f"@misc{{{key},\n"
                f"  author       = {{{authors}}},\n"
                f"  title        = {{{title}}},\n"
                f"  year         = {{{year}}},\n"
                f"  archivePrefix = {{arXiv}},\n"
                f"  eprint       = {{{aid}}},\n"
                f"}}"
            )
            return {
                "verdict": "ready_to_insert",
                "matching_key": None,
                "candidates": [lookup.get("metadata", {})],
                "suggested_bibtex": bibtex,
                "verification_method": "arxiv-direct",
                "advice": (f"Verified via arXiv (no DOI on file).  Append "
                            f"suggested_bibtex, then \\cite{{{key}}}."),
            }
        # S2 lookup failed
        return {
            "verdict": "no_candidate_found",
            "matching_key": None,
            "candidates": [],
            "suggested_bibtex": None,
            "verification_method": "arxiv-direct",
            "advice": (f"arXiv ID {candidate_arxiv_id!r} could not be "
                        f"resolved via Semantic Scholar.  Verify the ID "
                        f"is correct and try again, or supply a DOI."),
        }

    # --- Step 2c: title-only fast-path (check bib first) ---
    if candidate_title:
        matching = _title_already_cited(bib_entries, candidate_title)
        if matching:
            return {
                "verdict": "found_in_bib",
                "matching_key": matching,
                "candidates": [],
                "suggested_bibtex": None,
                "verification_method": "title-match",
                "advice": (f"Title approximately matches existing bib "
                            f"entry \\cite{{{matching}}}.  Reuse this key "
                            f"OR verify with user that they intend a "
                            f"different paper."),
            }

    # --- Step 3: fall back to arXiv keyword search ---
    if search_arxiv_if_no_doi and (candidate_title or claim):
        from ._arxiv_source import paper_writing_arxiv_search
        query = candidate_title or claim
        hits = paper_writing_arxiv_search(query, max_results=arxiv_max_results)
        if "error" in hits or hits.get("n_results", 0) == 0:
            return {
                "verdict": "no_candidate_found",
                "matching_key": None,
                "candidates": [],
                "suggested_bibtex": None,
                "verification_method": "arxiv-search",
                "advice": (f"No candidates found via arXiv search for "
                            f"{query!r}.  Mark this citation as a TODO "
                            f"in the manuscript: "
                            f"\\cite{{TODO: verify -- {claim[:60]}...}}.  "
                            f"Ask the user for the DOI or paper title."),
            }
        cands = hits["papers"]
        # Filter out already-cited
        for c in cands:
            t = c.get("title", "")
            matching = _title_already_cited(bib_entries, t)
            if matching:
                c["already_in_bib_key"] = matching
        return {
            "verdict": "needs_disambiguation",
            "matching_key": None,
            "candidates": cands,
            "suggested_bibtex": None,
            "verification_method": "arxiv-search",
            "advice": (f"Found {len(cands)} candidates via arXiv search "
                        f"for {query!r}.  Ask the user to pick one and "
                        f"re-call paper_writing_verify_citation with the "
                        f"chosen arxiv_id."),
        }

    # --- Step 4: insufficient info ---
    return {
        "verdict": "no_candidate_found",
        "matching_key": None,
        "candidates": [],
        "suggested_bibtex": None,
        "verification_method": None,
        "advice": ("Insufficient information to verify.  Provide at least "
                    "one of: candidate_doi, candidate_arxiv_id, "
                    "candidate_title.  Or pass a more detailed claim "
                    "(>20 chars) and set search_arxiv_if_no_doi=True."),
    }


# ============================================================
# paper_writing_check_citation_keys_exist  (2026-05-27)
# Static check: every \cite{key} resolves to a .bib entry.
# ============================================================

# All cite-family commands.  Captures the keys-list inside {...}.
_CITE_KEYS_RE = re.compile(
    r"\\(?:cite|citet|citep|citeyear|citealt|citealp|citeauthor|"
    r"nocite|fullcite|textcite|parencite|footcite|autocite|"
    r"citeNP|citeA|citeauthorNP)"
    r"(?:\[[^\]]*\])?"            # optional [page]
    r"(?:\[[^\]]*\])?"            # optional second [prenote]
    r"\{([^\}]+)\}"
)


def _collect_cited_keys(tex_source: str) -> dict[str, list[int]]:
    """Pull every key out of every \\cite{a,b,c} occurrence.

    Returns dict ``{key: [char_offset_first_seen, ...]}``.
    A single key may be cited many times; we record only the FIRST
    offset (used for reporter context).
    """
    out: dict[str, list[int]] = {}
    for m in _CITE_KEYS_RE.finditer(tex_source):
        keys_blob = m.group(1)
        for raw in keys_blob.split(","):
            k = raw.strip()
            if not k:
                continue
            out.setdefault(k, []).append(m.start())
    return out


def paper_writing_check_citation_keys_exist(
    tex_path: str,
    bib_path: str = "",
    auto_resolve_inputs: bool = True,
    encoding: str = "utf-8",
    max_reported: int = 50,
) -> dict:
    """Static check: every ``\\cite{key}`` in the .tex resolves to a
    .bib entry, and vice-versa report unused bib entries.

    This is the cheapest sanity check for the reviewer-pattern
    "reference [12] does not exist" — bibtex / biber will mark the
    citation as ``???`` in the rendered PDF, but only LaTeX-compile
    errors that the author may have ignored.  This MCP tool surfaces
    them up-front without compiling.

    Args:
        tex_path: main .tex file (paper entry point).  If the project
            uses ``\\input{}`` to split content across files, set
            ``auto_resolve_inputs=True`` (default) to merge the chain
            before scanning.
        bib_path: optional external bibliography override. The bundled
            canonical ``references.bib`` is used when omitted.
        auto_resolve_inputs: True (default) inlines ``\\input{}``
            recursively before scanning ``\\cite{}``.
        encoding: tex / bib file encoding.  ``"utf-8"`` default;
            ``"cp932"`` for legacy Japanese tex.
        max_reported: cap the per-list report length in the response
            (full lists still under ``all_missing_in_bib`` /
            ``all_unused_in_tex``).

    Returns:
        dict with::

            {
              "tex_path": echoed,
              "bib_path": echoed,
              "n_cited_keys": int,           # unique keys cited in tex
              "n_bib_entries": int,          # unique keys in bib
              "n_missing_in_bib": int,       # cited but bib doesn't have
              "n_unused_in_tex": int,        # bib has but never cited
              "missing_in_bib": [
                  {"key": "Smith2024",
                   "n_occurrences": 3,
                   "first_offset": 12345,
                   "first_context": "... as shown in \\cite{Smith2024} ..."},
                  ...
              ],
              "unused_in_tex": ["Foo2020", ...],  # alphabetical
              "all_missing_in_bib": [str, ...],   # alphabetical
              "all_unused_in_tex": [str, ...],
              "truncated_missing": bool,
              "truncated_unused": bool,
              "status": "clean" | "warning" | "fail",
              "advice": str,
            }

    Status logic:
      * ``fail``:    any key cited in .tex is missing from .bib
                     (compile-time undefined reference).
      * ``warning``: every cited key resolves, but the .bib has
                     unused entries (minor cleanup).
      * ``clean``:   one-to-one match.
    """
    bib_path = bib_path or _canonical_bib()
    if not os.path.exists(tex_path):
        return {"error": f"tex_path not found: {tex_path}"}
    if not os.path.exists(bib_path):
        return {"error": f"bib_path not found: {bib_path}"}

    # 1. Read tex (with optional \input chain inlining)
    try:
        raw = pathlib.Path(tex_path).read_bytes()
        try:
            src = raw.decode(encoding, errors="strict")
        except UnicodeDecodeError:
            if encoding.casefold().replace("-", "") != "utf8":
                raise
            src = raw.decode("cp932", errors="strict")
    except (OSError, UnicodeError) as e:
        return {"error": f"failed to read tex: {e}"}

    if auto_resolve_inputs:
        try:
            from ._tex_resolver import resolve_input_chain
            r = resolve_input_chain(tex_path, encoding=encoding)
            if not r.get("ok"):
                return {"error": f"failed to resolve TeX inputs: {r.get('error')}"}
            src = r["merged_tex"]
        except Exception as exc:  # noqa: BLE001
            return {"error": f"failed to resolve TeX inputs: {exc}"}

    cited = _collect_cited_keys(src)
    try:
        bib_entries = _parse_bib_lightweight(bib_path)
    except (OSError, UnicodeError, ValueError) as exc:
        return {"error": f"failed to read bib: {exc}"}
    bib_keys = set(bib_entries.keys())

    missing = sorted(set(cited.keys()) - bib_keys)
    unused = sorted(bib_keys - set(cited.keys()))

    # Build per-missing-key context snippet
    missing_records = []
    for k in missing[:max_reported]:
        first = cited[k][0]
        a = max(0, first - 50)
        b = min(len(src), first + 60)
        ctx = src[a:b].replace("\n", " ")
        missing_records.append({
            "key": k,
            "n_occurrences": len(cited[k]),
            "first_offset": first,
            "first_context": ctx,
        })

    if missing:
        status = "fail"
        advice = (
            f"{len(missing)} citation key(s) used in .tex but NOT in {os.path.basename(bib_path)}. "
            "These will render as '[?]' in the PDF and bibtex will emit "
            "'undefined references'.  Fix by either (a) adding a BibTeX entry "
            "for each missing key, or (b) correcting the typo in the cite key.  "
            "See `paper_writing_arxiv_fetch_latex_source` / `paper_writing_verify_citation` "
            "for help building new entries."
        )
    elif unused:
        status = "warning"
        advice = (
            f"All {len(cited)} cited key(s) resolve, but {len(unused)} bib entry(s) "
            "are never cited.  Optional cleanup: remove unused entries from "
            f"{os.path.basename(bib_path)} to keep the bibliography lean, OR add "
            "the missing \\cite{} where the prior work is mentioned."
        )
    else:
        status = "clean"
        advice = (
            f"All {len(cited)} cited key(s) resolve to .bib entries with no "
            "leftover unused entries — bibliography is in sync."
        )

    return {
        "tex_path": tex_path,
        "bib_path": bib_path,
        "n_cited_keys": len(cited),
        "n_bib_entries": len(bib_keys),
        "n_missing_in_bib": len(missing),
        "n_unused_in_tex": len(unused),
        "missing_in_bib": missing_records,
        "unused_in_tex": unused[:max_reported],
        "all_missing_in_bib": missing,
        "all_unused_in_tex": unused,
        "truncated_missing": len(missing) > max_reported,
        "truncated_unused": len(unused) > max_reported,
        "status": status,
        "advice": advice,
    }
