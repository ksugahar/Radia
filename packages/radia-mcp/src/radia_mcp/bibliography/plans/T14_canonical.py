"""The lab's single bibliography, and the two operations that make it usable.

Papers do not keep their own ``references.bib``. They cite the canonical file
here, and what travels with a submission is the generated ``.bbl``, which holds
only the works actually cited. That keeps one place to fix a DOI while leaving
each manuscript self-contained for a publisher.

Two operations follow from that:

``bibliography_make_bbl``
    Read the citation keys out of a ``.tex``, then let BibTeX format exactly
    those entries from the canonical file. A throwaway driver document is
    compiled rather than the manuscript itself, so a paper whose figures are
    missing still produces a correct bibliography.

``bibliography_refresh_unpublished``
    An entry written while a paper was in preparation carries no DOI, and stays
    that way long after the paper appears -- one entry here still said "to be
    published" for something published in 2020. This asks Crossref about every
    DOI-less entry and reports close title matches for review. It never writes.
"""
from __future__ import annotations

import difflib
import pathlib
import re
import shutil
import subprocess
import tempfile

from .._bibparse import read_bib_file

CANONICAL = pathlib.Path(__file__).resolve().parents[1] / "data" / "references.bib"

# A Japanese title must survive normalisation; stripping to [a-z0-9]
# erases it completely and two copies of the same paper stop matching.
CJK_RANGE = "\u3040-\u30ff\u3400-\u4dbf\u4e00-\u9fff\uff66-\uff9f"
_CITE = re.compile(r"\\(?:no)?cite[a-zA-Z]*\s*(?:\[[^\]]*\])*\s*\{([^}]*)\}")


def bibliography_canonical_path() -> str:
    """Return the path of the lab's canonical bibliography and its size."""
    if not CANONICAL.exists():
        return f"canonical bibliography not found at {CANONICAL}"
    n = sum(1 for e in read_bib_file(CANONICAL) if not e.kind.startswith("@"))
    return (f"bibliography_canonical_path\n  {CANONICAL}\n"
            f"  {n} entries, {CANONICAL.stat().st_size/1024:.0f} KB\n"
            "  Papers cite this file; they do not keep a copy. Ship the .bbl.")


def _keys_in_order(tex: str) -> list[str]:
    """Citation keys in first-appearance order; BibTeX styles depend on it."""
    seen, out = set(), []
    for group in _CITE.findall(tex):
        for k in (x.strip() for x in group.split(",")):
            if k and k not in seen:
                seen.add(k)
                out.append(k)
    return out


def bibliography_make_bbl(tex_path: str, style: str = "",
                          out_path: str | None = None) -> str:
    """Build a .bbl for one manuscript from the canonical bibliography.

    Once the .bbl exists beside the manuscript, LaTeX needs no .bib at all --
    \\bibliography{} pulls in \\jobname.bbl, and the .bib is read only when
    BibTeX itself runs. That is what lets a paper folder hold no bibliography
    while staying self-contained for a publisher.

    tex_path : the manuscript; only its \\cite keys are read.
    style    : BibTeX style; taken from the manuscript's
               \\bibliographystyle{} when omitted, IEEEtran if it has none.
    out_path : where to write the .bbl; defaults to beside the manuscript.
    """
    src = pathlib.Path(tex_path)
    if not src.exists():
        return f"no such file: {src}"
    if not style:
        m = re.search(r"\\bibliographystyle\s*\{([^}]*)\}",
                      src.read_text(encoding="utf-8", errors="replace"))
        style = (m.group(1).strip() if m else "") or "IEEEtran"
    if not CANONICAL.exists():
        return f"canonical bibliography missing: {CANONICAL}"
    if shutil.which("bibtex") is None:
        return "bibtex is not on PATH; a TeX installation is required"

    keys = _keys_in_order(src.read_text(encoding="utf-8", errors="replace"))
    if not keys:
        return f"no \\cite keys found in {src.name}"

    known = {e.key for e in read_bib_file(CANONICAL)}
    missing = [k for k in keys if k not in known]

    with tempfile.TemporaryDirectory() as tmp:
        d = pathlib.Path(tmp)
        shutil.copy(CANONICAL, d / "references.bib")
        # a journal's own .bst is not in the TeX tree; without it BibTeX writes
        # an empty bibliography and still exits cleanly
        for bst in list(src.parent.glob("*.bst")) + \
                list(src.parent.parent.glob("*.bst")):
            shutil.copy(bst, d / bst.name)
        # a driver, not the manuscript: missing figures must not break this
        (d / "drv.tex").write_text(
            "\\documentclass{article}\\begin{document}\n"
            + "".join(f"\\nocite{{{k}}}\n" for k in keys if k in known)
            + f"\\bibliographystyle{{{style}}}\n"
              "\\bibliography{references}\n\\end{document}\n",
            encoding="utf-8")
        # never decode with the console codepage: a log line outside cp932
        # raises from a reader thread once the work is already finished
        subprocess.run(["latex", "-interaction=nonstopmode", "drv.tex"],
                       cwd=d, capture_output=True, timeout=180)
        r = subprocess.run(["bibtex", "drv"], cwd=d, capture_output=True,
                           text=True, encoding="utf-8", errors="replace",
                           timeout=180)
        bbl = d / "drv.bbl"
        if not bbl.exists():
            return "bibtex produced no .bbl:\n" + (r.stdout or "")[-1200:]
        if "\\bibitem" not in bbl.read_text(encoding="utf-8", errors="replace"):
            return (f"bibtex produced an EMPTY bibliography for {src.name} "
                    f"(style {style}). The .bbl was NOT written.\n"
                    + (r.stdout or "")[-800:])
        dst = pathlib.Path(out_path) if out_path else src.with_suffix(".bbl")
        dst.write_text(bbl.read_text(encoding="utf-8", errors="replace"),
                       encoding="utf-8")
        n = bbl.read_text(encoding="utf-8", errors="replace").count("\\bibitem")

    lines = [f"bibliography_make_bbl: {src.name} -> {dst}",
             f"  cited {len(keys)} keys, wrote {n} bibitems, style {style}"]
    if missing:
        lines.append(f"  NOT in the canonical file ({len(missing)}): "
                     + ", ".join(missing[:12]))
        lines.append("  add them there; do not create a local references.bib")
    err = [l for l in (r.stdout or "").splitlines()
           if l.startswith(("I was", "You're", "Repeated"))]
    if err:
        lines.append("  bibtex errors: " + "; ".join(err[:3]))
    return "\n".join(lines)


def bibliography_verify_dois(limit: int = 0, delay: float = 0.3) -> str:
    """Check that every DOI in the canonical file resolves, and to the right work.

    A closed bibliography only stops invented citations if the bibliography
    itself is clean. Two failures are caught here:

      unresolvable -- Crossref has no such DOI. Either a typo or a fabrication;
        this repository carried ``10.1109/TMAG.2018.2848702`` on the Kameari
        2018 CLN paper, which does not exist. The real one is
        ``10.1109/TMAG.2017.2743224``.
      mismatched -- the DOI resolves, but to a different paper than the entry
        claims. The nastier case, because the reference list still compiles and
        looks plausible.
    """
    import time

    from .T1_doi_to_bibtex import _fetch_crossref

    def plain(s: str) -> str:
        s = re.sub(r"\\[a-zA-Z]+", "", s or "")
        return re.sub(r"[^a-z0-9]", "", s.lower())

    entries = [e for e in read_bib_file(CANONICAL) if e.fields.get("doi")]
    if limit:
        entries = entries[:limit]

    dead, wrong, ok, skipped = [], [], 0, 0
    for e in entries:
        doi = e.fields["doi"].strip().rstrip("}").replace("https://doi.org/", "")
        # Registries other than Crossref: arXiv mints through DataCite, and
        # JSIAM and its neighbours through JaLC. An absent record proves
        # nothing about these, so they are not counted as dead.
        if doi.lower().startswith(("10.48550/", "10.11540/", "10.14947/")):
            skipped += 1
            continue
        try:
            msg = _fetch_crossref(doi)
        except Exception:
            msg = None
        time.sleep(delay)
        if not msg:
            dead.append((e.key, doi))
            continue
        theirs = msg.get("title") or [""]
        theirs = plain(theirs[0] if isinstance(theirs, list) else theirs)
        ours = plain(e.fields.get("title", ""))
        # a book listed without its subtitle is the same book, so containment
        # counts as agreement before falling back to a similarity ratio
        # A Japanese title and its English translation share no characters, so
        # a ratio of zero here means the languages differ, not the works.
        ja = re.search(r"[\u3040-\u30ff\u4e00-\u9fff]",
                       e.fields.get("title", "") or "")
        same = (ours and theirs
                and (bool(ja) and theirs.isascii()
                     or ours.startswith(theirs) or theirs.startswith(ours)
                     or difflib.SequenceMatcher(
                         None, ours[:80], theirs[:80]).ratio() >= 0.70))
        if ours and theirs and not same:
            wrong.append((e.key, doi, e.fields.get("title", "")[:52],
                          (msg.get("title") or [""])[0][:52]))
        else:
            ok += 1

    out = [f"bibliography_verify_dois: {len(entries)} entries with a DOI",
           f"  resolved and matching : {ok}",
           f"  DOI does not resolve  : {len(dead)}",
           f"  resolves to another work: {len(wrong)}",
           f"  other registry (arXiv/JaLC, not in Crossref): {skipped}"]
    for k, d in dead[:20]:
        out.append(f"  DEAD  {k}: {d}")
    for k, d, o, t in wrong[:20]:
        out.append(f"  WRONG {k}: {d}")
        out.append(f"          ours: {o}")
        out.append(f"          xref: {t}")
    return "\n".join(out)


def bibliography_refresh_unpublished(min_year: str = "2024",
                                     threshold: float = 0.80) -> str:
    """Report canonical entries that Crossref now knows about.

    Reports only. Review each match before editing the entry: a short or
    generic title can match a different paper.
    """
    from .T3_search_crossref import bibliography_search_crossref

    def plain(s: str) -> str:
        s = re.sub(r"\\[a-zA-Z]+", "", s or "")
        return re.sub(r"[^a-z0-9 ]", "", s.lower()).strip()

    todo = [e for e in read_bib_file(CANONICAL)
            if not e.fields.get("doi") and e.fields.get("title")
            and e.fields.get("year", "0") >= min_year
            and not re.search(r"(?i)web ?site|github", e.fields["title"])]

    out = [f"bibliography_refresh_unpublished: {len(todo)} DOI-less entries "
           f"from {min_year}"]
    hits = 0
    for e in todo:
        title = plain(e.fields["title"])
        if not title:
            continue
        try:
            res = bibliography_search_crossref(title[:110], 3)
        except Exception as exc:                      # network / API trouble
            out.append(f"  ? {e.key}: lookup failed ({exc})")
            continue
        best, line = 0.0, ""
        for cand_line in res.splitlines():
            m = re.match(r"\s*\d+\.\s*\[(\d{4})\]\s*(.+?)\s+—\s+(.+)", cand_line)
            if not m:
                continue
            r = difflib.SequenceMatcher(None, title[:90],
                                        plain(m.group(3))[:90]).ratio()
            if r > best:
                best, line = r, cand_line.strip()
        if best >= threshold:
            hits += 1
            out.append(f"  * {best:.2f} {e.key}")
            out.append(f"      ours: {e.fields['title'][:70]}")
            out.append(f"      xref: {line[:100]}")
    out.append(f"  {hits} candidate(s); verify before editing the entry.")
    return "\n".join(out)


def bibliography_find_stray_bibs(root: str = r"W:\02_学会資料",
                                 max_files: int = 200) -> str:
    """List .bib files outside the canonical one, and whether they can go.

    The lab keeps one bibliography. A local copy in a paper folder is the thing
    that lets a DOI be right in one place and wrong in another, so this reports
    them and, for each, how many entries the canonical file already covers:

      absorbed   -- every entry is already central; the file can be deleted once
                    its .tex points at the canonical bibliography.
      to merge   -- entries not present centrally; merge those first, or they
                    are lost with the file.

    Reports only; it deletes nothing.
    """
    base = pathlib.Path(root)
    if not base.exists():
        return f"no such directory: {base}"

    def ident(e):
        t = re.sub(r"\\[a-zA-Z]+", "", (e.fields.get("title") or "").lower())
        t = re.sub("[^a-z0-9" + CJK_RANGE + "]", "", t)[:80]
        return t or (e.fields.get("doi") or e.key).lower()

    canon = {ident(e) for e in read_bib_file(CANONICAL)
             if not e.kind.startswith("@")}
    rows, unread, total_missing = [], [], 0
    for f in sorted(base.rglob("*.bib")):
        # _retired_bib holds files already taken out of paper folders;
        # counting them would make the migration look stalled
        # _retired_bib holds files already taken out of paper folders;
        # counting them would make the migration look stalled
        if (f.resolve() == CANONICAL.resolve() or ".git" in f.parts
                or "_retired_bib" in f.parts):
            continue
        try:
            entries = [e for e in read_bib_file(f) if not e.kind.startswith("@")]
        except Exception as exc:
            unread.append((f, exc))
            continue
        missing = [e for e in entries if ident(e) not in canon]
        total_missing += len(missing)
        rows.append((f, len(entries), missing))

    out = [f"bibliography_find_stray_bibs: {base}",
           f"  {len(rows)} 個のローカル .bib（正典以外）",
           f"  正典に無いエントリの総数: {total_missing}"]
    ready = [r for r in rows if not r[2]]
    out.append(f"  そのまま削除できる（全件が正典にある）: {len(ready)} ファイル")
    out.append("  ※ 削除前に、その .tex が正典を参照しているか確認すること。")

    for f, n, missing in rows[:max_files]:
        try:
            rel = f.relative_to(base)
        except ValueError:
            rel = f
        if missing:
            out.append(f"  [要マージ {len(missing):3d}/{n:3d}] {rel}")
            for e in missing[:3]:
                t = re.sub(r"[{}]", "", e.fields.get("title", ""))[:56]
                out.append(f"        - {e.key}: {t}")
            if len(missing) > 3:
                out.append(f"        ... 他 {len(missing)-3} 件")
        else:
            out.append(f"  [削除可   {n:3d}件] {rel}")
    if len(rows) > max_files:
        out.append(f"  ... 他 {len(rows)-max_files} ファイル")
    for f, exc in unread:
        out.append(f"  [読めない] {f.name}: {exc}")
    return "\n".join(out)


def bibliography_check_keys(bib_path: str = "") -> str:
    """Report cite keys that cannot survive a real BibTeX run.

    Two failures are silent, which is what makes them worth a check:

      non-ASCII -- classic BibTeX drops the entry with no error, so the
        reference is missing from the printed list while everything looks fine.
      LaTeX-special characters -- & % $ # ~ ^ backslash and braces are acted on
        while \\cite{} is scanned, so the citation never resolves.

    Also lists keys that carry no information about what they cite (POD, ref3),
    which are not broken but collide as soon as a second such paper arrives.
    """
    path = pathlib.Path(bib_path) if bib_path else CANONICAL
    entries = [e for e in read_bib_file(path) if not e.kind.startswith("@")]

    special = set("&%$#{}~^" + chr(92))
    nonascii = [e.key for e in entries if not e.key.isascii()]
    hostile = [(e.key, "".join(sorted(set(e.key) & special)))
               for e in entries if set(e.key) & special]
    vague = [e.key for e in entries
             if re.fullmatch(r"[A-Za-z]{2,6}\d*|ref\d+|bib\d+", e.key)]

    out = [f"bibliography_check_keys: {path.name}, {len(entries)} entries",
           f"  non-ASCII (BibTeX drops these silently): {len(nonascii)}",
           f"  LaTeX-special characters               : {len(hostile)}",
           f"  uninformative (collide on the next one): {len(vague)}"]
    for k in nonascii[:20]:
        out.append(f"  DROPPED  {k}")
    for k, ch in hostile[:20]:
        out.append(f"  UNCITABLE {k}  [{ch}]")
    for k in vague[:20]:
        out.append(f"  vague    {k}")
    if not (nonascii or hostile):
        out.append("  every key survives BibTeX.")
    return "\n".join(out)
