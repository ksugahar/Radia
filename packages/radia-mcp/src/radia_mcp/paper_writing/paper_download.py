"""IEEE Xplore + Emerald + Crossref paper download tools.

Provides programmatic PDF fetching for academic papers from
publishers that rely on cookie-seeded session anti-bot (IEEE Xplore,
Emerald Publishing).  Works from a host that already has institutional
access — no proxy / auth setup needed beyond what the publisher
sees from the calling host.

Tools:
  paper_writing_resolve_doi          - Crossref metadata lookup
  paper_writing_ieee_doi_to_arnumber - DOI redirect → IEEE arnumber
  paper_writing_ieee_download_pdf    - DOI or arnumber → PDF on disk
  paper_writing_emerald_download_pdf - Emerald article URL → PDF on disk

Origin: distilled from `feedback_ieee_emerald_pdf_download_works.md`
in the radia-mcp lab memory after the 2026-05-21 successful download
of 6 Hollaus papers from IEEE TMAG via curl+cookies.

CAUTION: only call these from a machine with legitimate institutional
access to the target publisher.  These tools do NOT bypass paywall
authentication — they accept whatever access the caller's IP has.
"""

from __future__ import annotations

import os
import re
import pathlib
import html
import urllib.parse
from typing import Optional

# requests is optional (lazy-imported via _require_requests below) --
# this module only runs when the user invokes a download tool. Keeping
# it lazy means the rest of radia_mcp.paper_writing imports cleanly on
# the lightweight CI matrix install (`pip install radia-mcp --no-deps`)
# without dragging in requests / urllib3 / charset-normalizer.
def _require_requests():
    """Import requests lazily with a clear install hint."""
    try:
        import requests  # type: ignore
        return requests
    except ImportError as e:
        raise ImportError(
            "paper_writing download tools require requests: "
            "pip install requests  (or pip install radia-mcp[rag])."
        ) from e


_CHROME_UA = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
    "AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/120.0.0.0 Safari/537.36"
)

_CROSSREF_UA = (
    "radia-mcp (mailto:ksugahar@ele.kindai.ac.jp) Crossref-Plus"
)

_HTML_HEADERS = {
    "User-Agent": _CHROME_UA,
    "Accept": (
        "text/html,application/xhtml+xml,application/xml;q=0.9,"
        "image/avif,image/webp,*/*;q=0.8"
    ),
    "Accept-Language": "en-US,en;q=0.9,ja;q=0.8",
    "Accept-Encoding": "gzip, deflate, br",
}

_PDF_HEADERS = {
    "User-Agent": _CHROME_UA,
    "Accept": "application/pdf,*/*",
    "Accept-Language": "en-US,en;q=0.9,ja;q=0.8",
}


def _verify_pdf(path: str) -> dict:
    """Return basic verification of a downloaded PDF.

    Uses PyMuPDF if available; otherwise just checks magic bytes
    + page count via a minimal parse.
    """
    if not os.path.exists(path):
        return {"ok": False, "error": f"file not found: {path}"}
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        head = f.read(8)
    if not head.startswith(b"%PDF"):
        return {
            "ok": False, "error": "not a PDF (magic bytes mismatch)",
            "first_bytes": head.hex(), "size_bytes": size,
        }
    try:
        import fitz
        with fitz.open(path) as doc:
            result = {
                "ok": True,
                "size_bytes": size,
                "page_count": doc.page_count,
                "title": doc.metadata.get("title", "") or "",
                "author": doc.metadata.get("author", "") or "",
            }
            # First-page text preview helps caller verify they got the
            # RIGHT paper (not a stale arnumber for a different article)
            if doc.page_count > 0:
                preview = doc.load_page(0).get_text()[:300]
                result["page1_preview"] = preview.replace("\n", " ")[:300]
            return result
    except ImportError:
        return {
            "ok": True, "size_bytes": size,
            "page_count": None,
            "note": "PyMuPDF not installed; size + magic OK",
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "ok": False,
            "error": f"invalid or unreadable PDF: {exc}",
            "size_bytes": size,
        }


def paper_writing_resolve_doi(doi: str) -> dict:
    """Look up a DOI's metadata via the Crossref public API.

    Useful for verifying that a DOI exists and getting the canonical
    title / authors / journal / year before attempting to download.
    Crossref is free, no auth required, no rate limit for polite use.

    Args:
        doi: DOI string like "10.1109/TMAG.2019.2943463" (with or
             without "https://doi.org/" prefix).

    Returns:
        {ok, doi, title, authors, year, journal, url} on success,
        {ok: False, error: ...} on failure.
    """
    doi = re.sub(r"^doi\s*:\s*", "", doi.strip(), flags=re.IGNORECASE)
    doi = re.sub(
        r"^(?:https?://(?:dx\.)?)?doi\.org/", "", doi, flags=re.IGNORECASE
    )
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/.")
    requests = _require_requests()
    try:
        r = requests.get(url, headers={"User-Agent": _CROSSREF_UA},
                         timeout=20)
    except Exception as e:
        return {
            "ok": False,
            "error": f"network error: {e}",
            "error_kind": "temporary",
            "temporary_failure": True,
        }
    if r.status_code != 200:
        # 404 is a stable "not found" answer. Rate limits, server errors,
        # access blocks, and other responses do not prove the DOI is absent.
        temporary = r.status_code != 404
        return {
            "ok": False,
            "error": f"crossref HTTP {r.status_code}",
            "error_kind": "temporary" if temporary else "not_found",
            "temporary_failure": temporary,
        }
    try:
        msg = r.json().get("message", {})
    except Exception as e:
        return {
            "ok": False,
            "error": f"json parse error: {e}",
            "error_kind": "temporary",
            "temporary_failure": True,
        }

    title = (msg.get("title") or [""])[0]
    authors = []
    for author in msg.get("author", []):
        personal = (
            author.get("family", "") + ", " + author.get("given", "")
        ).strip(", ")
        name = personal or author.get("name", "")
        if name:
            authors.append(name)
    container = (msg.get("container-title") or [""])[0]
    date_record = msg.get("published") or msg.get("issued") or msg.get("created") or {}
    parts = date_record.get("date-parts", [[None]])
    year = parts[0][0] if parts and parts[0] else None

    return {
        "ok": True,
        "doi": doi,
        "title": title,
        "authors": authors,
        "year": year,
        "journal": container,
        "url": msg.get("URL", f"https://doi.org/{doi}"),
        "type": msg.get("type", ""),
    }


def paper_writing_ieee_doi_to_arnumber(doi: str) -> dict:
    """Resolve an IEEE DOI to its IEEE Xplore arnumber.

    IEEE DOIs (10.1109/...) redirect via doi.org to a URL of the
    form `https://ieeexplore.ieee.org/document/<arnumber>/`.  This
    tool follows the redirect and extracts the arnumber.

    Args:
        doi: DOI like "10.1109/TMAG.2019.2943463" (with or without
             https://doi.org/ prefix).

    Returns:
        {ok, doi, arnumber, abstract_url} on success,
        {ok: False, error: ...} on failure.
    """
    doi = doi.strip().removeprefix("https://doi.org/").removeprefix("doi.org/")
    if not doi.lower().startswith("10.1109/"):
        return {"ok": False, "error": f"not an IEEE DOI: {doi}"}

    requests = _require_requests()
    session = requests.Session()
    try:
        r = session.get(f"https://doi.org/{doi}",
                        headers=_HTML_HEADERS,
                        allow_redirects=True, timeout=30)
    except Exception as e:
        return {"ok": False, "error": f"network error: {e}"}
    if r.status_code != 200:
        return {"ok": False, "error": f"DOI redirect HTTP {r.status_code}",
                "url_attempted": str(r.url)}

    # Final URL should be https://ieeexplore.ieee.org/document/<arnumber>/
    m = re.search(r"ieeexplore\.ieee\.org/document/(\d+)", r.url)
    if not m:
        return {
            "ok": False,
            "error": f"could not extract arnumber from URL: {r.url}",
        }
    arn = m.group(1)
    return {
        "ok": True,
        "doi": doi,
        "arnumber": arn,
        "abstract_url": f"https://ieeexplore.ieee.org/document/{arn}/",
    }


def paper_writing_ieee_download_pdf(
    doi_or_arnumber: str,
    dest_path: str,
    overwrite: bool = False,
) -> dict:
    """Download an IEEE Xplore PDF via cookie-seeded curl-like session.

    Requires the caller's IP to have institutional access to IEEE
    Xplore (i.e. a subscribing institution's network).  The cookie-seeded pattern:

    1. Visit the abstract page → publisher sets session cookies
    2. Request stampPDF/getPDF.jsp with Referer header → PDF returned

    Args:
        doi_or_arnumber: Either an IEEE DOI ("10.1109/TMAG...") or
                          a bare arnumber ("8954674").
        dest_path:        Absolute path where the PDF should be saved.
                          Directory must already exist.
        overwrite:        If False (default), refuses to overwrite an
                          existing file.

    Returns:
        {ok, dest_path, arnumber, size_bytes, page_count, title,
         page1_preview} on success.
        {ok: False, error: ...} on failure.

    Important: After download, the caller should verify
    `page1_preview` contains the expected title — wrong arnumbers
    give wrong PDFs (no error from the server).
    """
    inp = doi_or_arnumber.strip()
    # Determine arnumber
    if inp.isdigit():
        arn = inp
        doi = None
    else:
        # Treat as DOI
        res = paper_writing_ieee_doi_to_arnumber(inp)
        if not res["ok"]:
            return res
        arn = res["arnumber"]
        doi = res["doi"]

    # File-overwrite guard (safety per executing-actions-with-care policy)
    if os.path.exists(dest_path) and not overwrite:
        return {
            "ok": False,
            "error": f"dest_path already exists: {dest_path}. "
                     f"Pass overwrite=True to replace.",
        }
    parent = os.path.dirname(dest_path)
    if parent and not os.path.exists(parent):
        return {"ok": False,
                "error": f"parent directory does not exist: {parent}"}

    requests = _require_requests()
    session = requests.Session()
    abstract_url = f"https://ieeexplore.ieee.org/document/{arn}/"
    pdf_url = (f"https://ieeexplore.ieee.org/stampPDF/getPDF.jsp"
               f"?tp=&arnumber={arn}&ref=")

    # Step 1: visit abstract page to seed session cookies
    try:
        r1 = session.get(abstract_url, headers=_HTML_HEADERS, timeout=30)
    except Exception as e:
        return {"ok": False, "error": f"abstract fetch error: {e}"}
    if r1.status_code != 200:
        return {
            "ok": False,
            "error": f"abstract HTTP {r1.status_code} for arnumber {arn}",
            "url": abstract_url,
        }

    # Step 2: download PDF with Referer header
    pdf_headers = dict(_PDF_HEADERS)
    pdf_headers["Referer"] = abstract_url
    try:
        r2 = session.get(pdf_url, headers=pdf_headers,
                         timeout=120, stream=True)
    except Exception as e:
        return {"ok": False, "error": f"PDF fetch error: {e}"}
    if r2.status_code != 200:
        return {
            "ok": False,
            "error": f"PDF HTTP {r2.status_code} for arnumber {arn}",
            "url": pdf_url,
        }

    with open(dest_path, "wb") as f:
        for chunk in r2.iter_content(chunk_size=65536):
            if chunk:
                f.write(chunk)

    verify = _verify_pdf(dest_path)
    if not verify["ok"]:
        # PDF didn't validate — likely an HTML anti-bot page disguised
        # as PDF (Cloudflare or "Temporarily Unavailable")
        return {
            "ok": False,
            "error": "downloaded file is not a valid PDF",
            "dest_path": dest_path,
            "verify": verify,
            "hint": "IEEE may have served an anti-bot HTML page. "
                    "Check your IP has institutional access. "
                    "If this works in a browser but not here, the "
                    "session cookie pattern may need updating.",
        }

    return {
        "ok": True,
        "dest_path": dest_path,
        "arnumber": arn,
        "doi": doi,
        **verify,
    }


def paper_writing_sciencedirect_download_pdf(
    article_pdf_url: str,
    dest_path: str,
    article_landing_url: Optional[str] = None,
    overwrite: bool = False,
) -> dict:
    """Download an Elsevier ScienceDirect PDF.

    Same cookie-seeded session pattern as IEEE / Emerald.  Requires
    the caller's IP to have institutional access (i.e. a subscribing
    institution's ScienceDirect entitlement).

    Args:
        article_pdf_url:     The article PDF URL.  Find via the article
                              landing page → "Download PDF" button.
                              Typical pattern:
                              https://www.sciencedirect.com/science/article/pii/<PII>/pdfft?md5=...&pid=...-main.pdf
        dest_path:           Where to save the PDF.
        article_landing_url: Optional landing URL to seed cookies.
                              If None, derived from article_pdf_url by
                              stripping "/pdfft?..." → "/abs/<PII>".
        overwrite:           If False, refuses to overwrite existing.

    Returns:
        Same shape as paper_writing_ieee_download_pdf.
    """
    if os.path.exists(dest_path) and not overwrite:
        return {
            "ok": False,
            "error": f"dest_path already exists: {dest_path}. "
                     f"Pass overwrite=True to replace.",
        }
    parent = os.path.dirname(dest_path)
    if parent and not os.path.exists(parent):
        return {"ok": False,
                "error": f"parent directory does not exist: {parent}"}

    # Derive landing URL if not given (strip /pdfft?...)
    if not article_landing_url:
        # Extract PII from URL like /science/article/pii/<PII>/pdfft?...
        m = re.search(r"/science/article/pii/([A-Z0-9]+)", article_pdf_url)
        if m:
            article_landing_url = (
                f"https://www.sciencedirect.com/science/article/pii/{m.group(1)}"
            )
        else:
            article_landing_url = article_pdf_url  # fallback

    requests = _require_requests()
    session = requests.Session()
    try:
        r1 = session.get(article_landing_url, headers=_HTML_HEADERS,
                         timeout=30)
    except Exception as e:
        return {"ok": False, "error": f"landing-page fetch error: {e}"}

    pdf_headers = dict(_PDF_HEADERS)
    pdf_headers["Referer"] = article_landing_url
    try:
        r2 = session.get(article_pdf_url, headers=pdf_headers,
                         timeout=120, stream=True, allow_redirects=True)
    except Exception as e:
        return {"ok": False, "error": f"PDF fetch error: {e}"}
    if r2.status_code != 200:
        return {"ok": False,
                "error": f"PDF HTTP {r2.status_code}",
                "url": article_pdf_url}

    with open(dest_path, "wb") as f:
        for chunk in r2.iter_content(chunk_size=65536):
            if chunk:
                f.write(chunk)

    verify = _verify_pdf(dest_path)
    if not verify["ok"]:
        return {
            "ok": False,
            "error": "downloaded file is not a valid PDF",
            "dest_path": dest_path,
            "verify": verify,
            "hint": "ScienceDirect may require an active institutional "
                    "session.  If the article-pdf URL requires a token "
                    "(md5/pid query params), it may have expired -- "
                    "re-fetch from the landing page first.",
        }

    return {"ok": True, "dest_path": dest_path, **verify}


def _make_citation_key(authors: list, year: int, title: str) -> str:
    """Build the same lowercase author-year-word key as bibliography tools."""
    from ..bibliography._bibparse import BibEntry, make_cite_key

    entry = BibEntry(
        kind="article",
        key="",
        fields={
            "author": " and ".join(authors),
            "year": str(year) if year else "",
            "title": title,
        },
    )
    return make_cite_key(entry)


def paper_writing_doi_to_bibtex(doi: str,
                                  citation_key: Optional[str] = None) -> dict:
    """Generate a BibTeX entry for a DOI via Crossref metadata.

    Args:
        doi:          DOI string (with or without https://doi.org/ prefix).
        citation_key: Override the auto-generated AuthorYearWord key.

    Returns:
        {ok, doi, citation_key, bibtex} on success.
    """
    meta = paper_writing_resolve_doi(doi)
    if not meta["ok"]:
        return meta

    if not citation_key:
        citation_key = _make_citation_key(
            meta["authors"], meta["year"], meta["title"]
        )

    # Build BibTeX entry; entry type depends on Crossref 'type'
    ctype = meta.get("type", "journal-article")
    if "journal" in ctype:
        bib_type = "article"
    elif "proceedings" in ctype:
        bib_type = "inproceedings"
    elif "book" in ctype:
        bib_type = "book"
    else:
        bib_type = "misc"

    # Authors → "Last, First and Last, First and ..."
    authors_bib = " and ".join(meta["authors"])
    # Strip HTML/Crossref italic markers from title
    title_clean = re.sub(r"</?[a-zA-Z]+>", "", meta["title"])
    for _ in range(3):
        decoded = html.unescape(title_clean)
        if decoded == title_clean:
            break
        title_clean = decoded

    lines = [f"@{bib_type}{{{citation_key},",
              f"  title   = {{{title_clean}}},",
              f"  author  = {{{authors_bib}}},"]
    if meta.get("year") is not None:
        lines.append(f"  year    = {{{meta['year']}}},")
    if meta.get("journal"):
        if bib_type == "article":
            lines.append(f"  journal = {{{meta['journal']}}},")
        else:
            lines.append(f"  booktitle = {{{meta['journal']}}},")
    lines.append(f"  doi     = {{{meta['doi']}}},")
    lines.append(f"  url     = {{{meta['url']}}},")
    lines.append("}")
    bibtex = "\n".join(lines)

    return {
        "ok": True,
        "doi": meta["doi"],
        "citation_key": citation_key,
        "bibtex": bibtex,
        "metadata": meta,
    }


def paper_writing_fetch_and_cite(
    doi: str,
    dest_dir: str,
    filename: Optional[str] = None,
    citation_key: Optional[str] = None,
    overwrite: bool = False,
) -> dict:
    """One-shot: download IEEE PDF + generate BibTeX entry from a DOI.

    For non-IEEE DOIs, the BibTeX is still generated but no PDF is
    downloaded (returns pdf_path=None).

    Args:
        doi:          DOI string.
        dest_dir:     Directory to save the PDF.  Must exist.
        filename:     PDF filename.  If None, derived from BibTeX key.
        citation_key: Override the auto-generated BibTeX key.
        overwrite:    Refuse to overwrite existing PDF unless True.

    Returns:
        {ok, pdf_path (or None), bibtex, citation_key, metadata}.
    """
    # Step 1: resolve metadata + build BibTeX
    bib = paper_writing_doi_to_bibtex(doi, citation_key=citation_key)
    if not bib["ok"]:
        return bib

    citation_key = bib["citation_key"]
    meta = bib["metadata"]

    # Step 2: determine if IEEE → can download
    is_ieee = meta["doi"].lower().startswith("10.1109/")
    if not is_ieee:
        return {
            "ok": True,
            "pdf_path": None,
            "bibtex": bib["bibtex"],
            "citation_key": citation_key,
            "metadata": meta,
            "note": (f"non-IEEE DOI ({meta['doi']}) — BibTeX only, "
                     f"no PDF download attempted.  Use the "
                     f"paper_writing_emerald_download_pdf or "
                     f"_sciencedirect_download_pdf for those publishers."),
        }

    # Step 3: download PDF
    if not filename:
        # Sanitize title for filename: keep ASCII alphanumeric + spaces
        title_safe = re.sub(r"[^A-Za-z0-9 \-]", "", meta["title"]).strip()
        title_safe = re.sub(r"\s+", " ", title_safe)[:120]
        filename = f"{title_safe}.pdf" if title_safe else f"{citation_key}.pdf"

    if not os.path.exists(dest_dir):
        return {"ok": False, "error": f"dest_dir does not exist: {dest_dir}"}
    pdf_path = os.path.join(dest_dir, filename)

    dl = paper_writing_ieee_download_pdf(
        doi_or_arnumber=meta["doi"],
        dest_path=pdf_path,
        overwrite=overwrite,
    )
    if not dl["ok"]:
        return {
            "ok": False,
            "error": dl["error"],
            "stage": "download",
            "bibtex": bib["bibtex"],
            "citation_key": citation_key,
            "metadata": meta,
        }

    return {
        "ok": True,
        "pdf_path": pdf_path,
        "bibtex": bib["bibtex"],
        "citation_key": citation_key,
        "metadata": meta,
        "download": dl,
    }


def paper_writing_emerald_download_pdf(
    article_pdf_url: str,
    dest_path: str,
    article_landing_url: Optional[str] = None,
    overwrite: bool = False,
) -> dict:
    """Download an Emerald Publishing PDF (e.g. COMPEL journal).

    Emerald uses Cloudflare anti-bot.  The cookie-seeded session
    pattern (same as IEEE) works.

    Args:
        article_pdf_url:     The article-pdf URL ending in .pdf, e.g.
                              "https://www.emerald.com/compel/article-pdf/44/5/832/.../compel-12-2024-0520en.pdf"
        dest_path:           Where to save the PDF.
        article_landing_url: Optional article landing URL to seed cookies
                              from.  If None, derived from article_pdf_url
                              by removing "/article-pdf/" → "/article/".
        overwrite:           If False (default), refuses to overwrite.

    Returns:
        {ok, dest_path, size_bytes, page_count, title} on success,
        {ok: False, error: ...} on failure.
    """
    if os.path.exists(dest_path) and not overwrite:
        return {
            "ok": False,
            "error": f"dest_path already exists: {dest_path}. "
                     f"Pass overwrite=True to replace.",
        }
    parent = os.path.dirname(dest_path)
    if parent and not os.path.exists(parent):
        return {"ok": False,
                "error": f"parent directory does not exist: {parent}"}

    # Derive landing URL if not given
    if not article_landing_url:
        article_landing_url = article_pdf_url.replace(
            "/article-pdf/", "/article/"
        )
        # Trim the trailing .pdf and filename portion
        article_landing_url = re.sub(r"/[^/]+\.pdf$", "", article_landing_url)

    requests = _require_requests()
    session = requests.Session()
    try:
        r1 = session.get(article_landing_url, headers=_HTML_HEADERS,
                         timeout=30)
    except Exception as e:
        return {"ok": False, "error": f"landing-page fetch error: {e}"}

    pdf_headers = dict(_PDF_HEADERS)
    pdf_headers["Referer"] = article_landing_url
    try:
        r2 = session.get(article_pdf_url, headers=pdf_headers,
                         timeout=120, stream=True)
    except Exception as e:
        return {"ok": False, "error": f"PDF fetch error: {e}"}
    if r2.status_code != 200:
        return {"ok": False,
                "error": f"PDF HTTP {r2.status_code}",
                "url": article_pdf_url}

    with open(dest_path, "wb") as f:
        for chunk in r2.iter_content(chunk_size=65536):
            if chunk:
                f.write(chunk)

    verify = _verify_pdf(dest_path)
    if not verify["ok"]:
        return {
            "ok": False,
            "error": "downloaded file is not a valid PDF",
            "dest_path": dest_path,
            "verify": verify,
            "hint": "Emerald Cloudflare may have served a challenge page. "
                    "Retry, or check if you need to access from a different IP.",
        }

    return {
        "ok": True,
        "dest_path": dest_path,
        **verify,
    }
