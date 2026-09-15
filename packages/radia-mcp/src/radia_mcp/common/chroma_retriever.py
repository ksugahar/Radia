"""ChromaDB + sentence-transformers RAG infrastructure.

Adapted from wjc9011/COMSOL_Multiphysics_MCP src/knowledge/retriever.py
which uses the same stack to index 75 COMSOL official PDFs (~500MB)
into a persistent ChromaDB store, then expose semantic search to LLMs.

For radia_mcp the target corpus is public-safe curated corpus
(3,889 files / 51 GB) plus the freshly-organized 04_機械学習と最適化/.

Optional dependencies (graceful degradation if missing):
  - chromadb
  - sentence-transformers (defaults to all-MiniLM-L6-v2, ~80 MB;
    pass ``embedding_model="paraphrase-multilingual-MiniLM-L12-v2"``
    for better Japanese / Chinese retrieval on bilingual lab corpus)
  - pymupdf (for PDF text extraction)

Install:
    pip install chromadb sentence-transformers pymupdf

Usage:
    from radia_mcp.common.chroma_retriever import ChromaRetriever

    rag = ChromaRetriever(
        db_dir="C:/temp/radia_lit_chroma",
        collection="lit_em",
        embedding_model="all-MiniLM-L6-v2",
    )

    # Index a batch of chunks (each chunk = dict with at least
    # 'id', 'text', and a metadata dict).
    rag.add_chunks([
        {"id": "doc1_p1", "text": "Maxwell equations ...",
         "metadata": {"source": "MaxwellBook.pdf", "page": 1,
                       "language": "en"}},
        ...
    ])

    # Semantic search.
    hits = rag.search("rotational core loss tester", n_results=5)
    # -> [{"text": ..., "metadata": {...}, "score": 0.87}, ...]

    # Restrict to Japanese-language hits only (requires chunks to
    # have been tagged with language= at indexing time).
    hits = rag.search("ヒステリシス", n_results=5, language_filter="ja")

Design choice: this module does NOT depend on radia_mcp.literature_index
-- it is a generic utility usable by any subpackage. literature_index
adopts it as the engine for `literature_semantic_search`.

Multilingual support (2026-05-25, ported from wjc9011 COMSOL fork):
  - ``detect_filename_language(filename)`` -- CJK Unicode-range heuristic
  - ``extract_pdf_chunks(..., default_language=, auto_detect_language=)``
  - ``ChromaRetriever.search(..., language_filter="ja")`` filter
  - Multilingual chapter regex (JA: 第N章, 第N節; ZH: kanji numerals)
"""

from __future__ import annotations
import logging
import re
from pathlib import Path
from typing import Any, Optional, Sequence

logger = logging.getLogger(__name__)

DEFAULT_EMBEDDING_MODEL = "all-MiniLM-L6-v2"

# Multilingual chapter detection patterns -- callers that want chapter-
# aware splitting can run `find_chapters(text)` from this module.
# Patterns ordered: English (most common in lab PDFs) -> Japanese -> Chinese.
CHAPTER_PATTERNS: list[str] = [
    # English
    r"^Chapter\s+(\d+)[:\s]+(.+)$",
    r"^(\d+(?:\.\d+)*)\s+(.+)$",
    r"^([A-Z][A-Za-z\s]+)$",
    # Japanese: 第1章, 第10章
    r"^第\s*(\d+)\s*章[:：\s]*(.*)$",
    # Japanese: 1章, 2章 (章 prefix omitted)
    r"^(\d+)\s*章[:：\s]*(.*)$",
    # Japanese: 第1節, 1.1節
    r"^第\s*(\d+(?:\.\d+)*)\s*節[:：\s]*(.*)$",
    # Chinese: 第一章, 第二章 (kanji number)
    r"^第([一二三四五六七八九十百千]+)章[:：\s]*(.*)$",
]


def detect_filename_language(filename: str) -> Optional[str]:
    """Heuristic: guess language from PDF filename characters.

    Returns "ja" if hiragana / katakana present (definitively Japanese),
    "zh" if only CJK ideographs are present (likely Chinese -- no kana),
    "en" if only ASCII letters / digits / common punctuation,
    None if undetermined (mixed or empty).

    Cheap heuristic for auto-tagging at index-build time when the
    user has not passed an explicit ``default_language`` flag. For
    accurate detection on document BODY text use langdetect or
    fasttext (heavier dependencies).

    Ported from wjc9011/COMSOL_Multiphysics_MCP fork
    (src/knowledge/pdf_processor.py::detect_filename_language).
    """
    if not filename:
        return None
    name = Path(filename).stem
    has_hiragana_or_katakana = False
    has_cjk_ideograph = False
    has_ascii_letter = False
    for ch in name:
        cp = ord(ch)
        if 0x3040 <= cp <= 0x30FF:
            has_hiragana_or_katakana = True
        elif 0x4E00 <= cp <= 0x9FFF:
            has_cjk_ideograph = True
        elif ch.isascii() and ch.isalpha():
            has_ascii_letter = True
    if has_hiragana_or_katakana:
        return "ja"
    if has_cjk_ideograph:
        return "zh"  # ideograph-only filename more likely Chinese
    if has_ascii_letter:
        return "en"
    return None


def find_chapters(text: str) -> list[tuple[str, int, int]]:
    """Detect chapter boundaries in text (multilingual).

    Returns list of (chapter_title, line_start, line_end). Useful when
    a caller wants to attach a per-chunk ``chapter`` metadata field.

    Recognizes English, Japanese (第N章, 第N節, N章), and Chinese
    (kanji-numeral chapters) headings via CHAPTER_PATTERNS.
    """
    chapters: list[tuple[str, int, int]] = []
    lines = text.split("\n")
    current_chapter = "Introduction"
    chapter_start = 0
    for i, raw in enumerate(lines):
        line = raw.strip()
        if not line:
            continue
        for pattern in CHAPTER_PATTERNS:
            match = re.match(pattern, line)
            if match:
                if chapter_start < i:
                    chapters.append((current_chapter, chapter_start, i))
                if len(match.groups()) >= 2:
                    current_chapter = match.group(2).strip() or current_chapter
                else:
                    current_chapter = match.group(1).strip()
                chapter_start = i
                break
    if chapter_start < len(lines):
        chapters.append((current_chapter, chapter_start, len(lines)))
    return chapters


def _import_chromadb():
    """Lazy import of chromadb with a clear error message."""
    try:
        import chromadb
        from chromadb.utils import embedding_functions
        return chromadb, embedding_functions
    except ImportError as e:
        raise ImportError(
            "ChromaDB RAG requires the optional dependencies:\n"
            "  pip install chromadb sentence-transformers pymupdf\n"
            f"Underlying error: {e}"
        )


def _find_local_st_cache(model_name: str) -> Optional[str]:
    """Look for a locally-cached sentence-transformers model.

    Returns the snapshot path if found (offline mode), else None.
    """
    import os
    import glob
    base = os.path.expanduser("~/.cache/huggingface/hub")
    pattern = os.path.join(
        base, f"models--sentence-transformers--{model_name}",
        "snapshots", "*")
    matches = glob.glob(pattern)
    return matches[0] if matches else None


class ChromaRetriever:
    """Wrapper around ChromaDB persistent collection.

    Designed to be:
    - SAFE TO IMPORT even without chromadb installed (lazy load)
    - PERSISTENT across MCP server restarts (PersistentClient)
    - OFFLINE-CAPABLE (uses local sentence-transformers cache if present)
    """

    def __init__(
        self,
        db_dir: str | Path,
        collection: str = "default",
        embedding_model: str = DEFAULT_EMBEDDING_MODEL,
    ):
        self.db_dir = Path(db_dir)
        self.collection_name = collection
        self.embedding_model = embedding_model
        self._client = None
        self._collection = None
        self._embedding_fn = None

    @property
    def is_initialized(self) -> bool:
        return self._collection is not None

    def _get_embedding_fn(self):
        if self._embedding_fn is not None:
            return self._embedding_fn
        _, embedding_functions = _import_chromadb()

        # Prefer local cache if available (works offline)
        local_path = _find_local_st_cache(self.embedding_model)
        model_name = local_path if local_path else self.embedding_model
        if local_path:
            logger.info(f"Using cached embedding model: {local_path}")

        self._embedding_fn = embedding_functions.SentenceTransformerEmbeddingFunction(
            model_name=model_name)
        return self._embedding_fn

    def initialize(self) -> bool:
        """Create / open the ChromaDB collection.

        Returns True on success, False on error (with log).
        """
        try:
            chromadb, _ = _import_chromadb()
            self.db_dir.mkdir(parents=True, exist_ok=True)
            self._client = chromadb.PersistentClient(path=str(self.db_dir))
            self._collection = self._client.get_or_create_collection(
                name=self.collection_name,
                embedding_function=self._get_embedding_fn(),
            )
            logger.info(
                f"ChromaDB collection '{self.collection_name}' open "
                f"({self._collection.count()} documents).")
            return True
        except ImportError:
            raise
        except Exception as e:
            logger.error(f"Failed to initialize ChromaDB: {e}")
            return False

    def add_chunks(
        self,
        chunks: Sequence[dict],
        batch_size: int = 100,
    ) -> int:
        """Add chunks to the collection.

        Each chunk must be a dict with:
            - "id": unique string
            - "text": str (the content to embed)
            - "metadata": dict of simple types (str/int/float/bool)

        Returns the number of chunks added.
        """
        if not self.is_initialized and not self.initialize():
            return 0
        if not chunks:
            return 0

        ids = [c["id"] for c in chunks]
        docs = [c["text"] for c in chunks]
        metas = [_flatten_metadata(c.get("metadata", {})) for c in chunks]

        added = 0
        for i in range(0, len(chunks), batch_size):
            self._collection.add(
                ids=ids[i:i + batch_size],
                documents=docs[i:i + batch_size],
                metadatas=metas[i:i + batch_size],
            )
            added += min(batch_size, len(chunks) - i)
        logger.info(f"Added {added} chunks to {self.collection_name}")
        return added

    def search(
        self,
        query: str,
        n_results: int = 5,
        where: Optional[dict] = None,
        language_filter: Optional[str] = None,
    ) -> list[dict]:
        """Semantic search.

        Args:
            query: free-text query (embedded with same model as docs)
            n_results: number of top hits to return
            where: optional ChromaDB metadata filter, e.g.
                   {"source": {"$eq": "MaxwellBook.pdf"}}
            language_filter: optional ISO 639-1 language tag (e.g.
                   "en", "ja", "zh"). ANDed with ``where`` if both are
                   given. Only effective when chunks were tagged with
                   ``language`` at index-build time -- see
                   ``extract_pdf_chunks(default_language=, auto_detect_language=)``.

        Returns:
            List of dicts: [{"id", "text", "metadata", "score"}, ...]
            (score in [0, 1] -- higher is more similar)
        """
        if not self.is_initialized and not self.initialize():
            return []

        # Build effective `where` filter combining caller-supplied where
        # with language_filter. ChromaDB requires multi-clause filters
        # to be wrapped in `$and`.
        eff_where: Optional[dict] = None
        if where and language_filter:
            eff_where = {"$and": [where, {"language": language_filter}]}
        elif where:
            eff_where = where
        elif language_filter:
            eff_where = {"language": language_filter}

        try:
            res = self._collection.query(
                query_texts=[query],
                n_results=n_results,
                where=eff_where,
                include=["documents", "metadatas", "distances"],
            )
        except Exception as e:
            logger.error(f"Search failed: {e}")
            return []

        out: list[dict] = []
        if res and res.get("documents") and res["documents"][0]:
            docs = res["documents"][0]
            metas = res.get("metadatas", [[]])[0] or [{}] * len(docs)
            dists = res.get("distances", [[]])[0] or [0] * len(docs)
            ids = res.get("ids", [[]])[0] or [""] * len(docs)
            for i, doc in enumerate(docs):
                # Convert L2 distance to similarity score in [0, 1]
                score = max(0.0, 1.0 - min(dists[i], 1.0))
                out.append({
                    "id": ids[i],
                    "text": doc,
                    "metadata": metas[i] or {},
                    "score": score,
                })
        return out

    def stats(self) -> dict:
        """Collection statistics."""
        if not self.is_initialized and not self.initialize():
            return {"initialized": False}
        n = self._collection.count()
        return {
            "initialized": True,
            "collection": self.collection_name,
            "db_dir": str(self.db_dir),
            "embedding_model": self.embedding_model,
            "n_docs": n,
        }

    def clear(self) -> bool:
        """Delete + recreate the collection."""
        if self._client is None and not self.initialize():
            return False
        try:
            self._client.delete_collection(self.collection_name)
            self._collection = None
            return self.initialize()
        except Exception as e:
            logger.error(f"Clear failed: {e}")
            return False


def _flatten_metadata(meta: dict) -> dict:
    """ChromaDB metadata must be flat (str/int/float/bool only)."""
    out = {}
    for k, v in meta.items():
        if v is None:
            continue
        if isinstance(v, (str, int, float, bool)):
            out[k] = v
        else:
            out[k] = str(v)
    return out


# ============================================================
# PDF helpers (optional pymupdf dependency)
# ============================================================

def _text_readability(s: str) -> float:
    """Fraction of chars that are alphanumeric (unicode/CJK-aware) or whitespace.

    A CORRUPT PDF text layer (CID font with a missing/wrong ToUnicode CMap --
    e.g. the Zienkiewicz FEM volumes) DISPLAYS fine but ``get_text()`` extracts
    control-char garbage -> very low ratio.  Clean text in ANY script (Latin,
    Japanese, ...) -> high ratio.  Used as the quality gate that triggers OCR.
    """
    if not s:
        return 0.0
    return sum(1 for c in s if c.isalnum() or c.isspace()) / len(s)


def chunk_garble_fraction(chunks: list[dict],
                          readability_threshold: float = 0.55) -> float:
    """Fraction of ``chunks`` whose text reads as garble (corrupt text layer).

    Post-extraction quality gate: a born-digital PDF with a corrupt CID /
    ToUnicode text layer (the Zienkiewicz FEM class) extracts text that is
    non-empty but meaningless -- ``_text_readability`` near zero. Computing the
    fraction of such chunks lets the ingestion path FLAG these books (so they
    can be re-OCR'd) instead of silently storing garbage. Returns 0.0 for empty
    input. Each chunk is a dict with a ``"text"`` key (as produced by
    :func:`extract_pdf_chunks`).
    """
    if not chunks:
        return 0.0
    g = sum(1 for c in chunks
            if _text_readability(c.get("text", "")) < readability_threshold)
    return g / len(chunks)


_EASYOCR_READER = None
_EASYOCR_LANGS = None


def _ocr_page_easyocr(page, languages, dpi: int) -> str:
    """OCR a rendered PDF page with easyocr (in-process; GPU if torch.cuda).
    Reader cached module-wide (construction loads models)."""
    global _EASYOCR_READER, _EASYOCR_LANGS
    import fitz  # pymupdf
    import numpy as np
    langs = tuple(languages)
    if _EASYOCR_READER is None or _EASYOCR_LANGS != langs:
        import easyocr
        # verbose=False: easyocr's U+2588 progress bar crashes a cp932 console.
        _EASYOCR_READER = easyocr.Reader(list(langs), verbose=False)
        _EASYOCR_LANGS = langs
    pix = page.get_pixmap(dpi=dpi, colorspace=fitz.csRGB)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)
    return "\n".join(_EASYOCR_READER.readtext(img, detail=0, paragraph=True))


def _ocr_page_ndlocr(page, dpi: int) -> str:
    """OCR a rendered PDF page with NDLOCR-Lite (the LAB'S PREFERRED OCR:
    layout-aware reading order, CPU-fast ~1.6s/page, handles Latin + Japanese --
    e.g. recovers the Zienkiewicz FEM volumes' garbled text layer to clean
    English).  Runs the ``ndlocr-lite`` CLI on a temp PNG.  Point NDLOCR_LITE_DIR
    at its ``src`` dir (default ``C:\\Tools\\ndlocr-lite\\src``)."""
    import fitz  # pymupdf
    import os
    import sys
    import shutil
    import tempfile
    import subprocess
    src = os.environ.get("NDLOCR_LITE_DIR", r"C:\Tools\ndlocr-lite\src")
    if not os.path.isfile(os.path.join(src, "ocr.py")):
        raise FileNotFoundError(
            f"NDLOCR-Lite not found at {src!r} (set env NDLOCR_LITE_DIR)")
    tmp = tempfile.mkdtemp(prefix="ndlocr_")
    try:
        png = os.path.join(tmp, "page.png")
        out = os.path.join(tmp, "out")
        os.makedirs(out, exist_ok=True)
        page.get_pixmap(dpi=dpi, colorspace=fitz.csRGB).save(png)
        subprocess.run([sys.executable, "ocr.py", "--sourceimg", png, "--output", out],
                       cwd=src, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        txt = os.path.join(out, "page.txt")
        if os.path.isfile(txt):
            with open(txt, encoding="utf-8", errors="replace") as fh:
                return fh.read()
        return ""
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _ocr_pdf_page(page, languages, dpi: int, engine: str = "ndlocr") -> str:
    """Rasterize + OCR a PDF page -> text.  engine='ndlocr' (default; the lab
    standard, layout-aware) or 'easyocr' (in-process, GPU)."""
    if engine == "easyocr":
        return _ocr_page_easyocr(page, languages, dpi)
    return _ocr_page_ndlocr(page, dpi)


def extract_pdf_chunks(
    pdf_path: str | Path,
    chunk_size: int = 1000,
    overlap: int = 100,
    default_language: Optional[str] = None,
    auto_detect_language: bool = False,
    ocr_fallback: bool = False,
    ocr_engine: str = "ndlocr",
    ocr_languages: tuple = ("en",),
    ocr_dpi: int = 200,
    ocr_min_readability: float = 0.6,
) -> list[dict]:
    """Extract text chunks from a PDF, ready to feed `add_chunks`.

    Simple page-based chunking with overlap. Returns chunks with:
        - id       = <stem>__p<page>__<chunk_in_page>
        - text     = the chunk content
        - metadata = {"source", "stem", "page",
                       "language" (optional)}

    Args:
        pdf_path: PDF file path
        chunk_size: characters per chunk (default 1000)
        overlap: chunk overlap (default 100)
        default_language: ISO 639-1 tag attached to every chunk
            (e.g. "en", "ja"). Use when an entire directory is known
            to be one language. None = no language metadata written.
        auto_detect_language: when True AND default_language is None,
            fall back to ``detect_filename_language(pdf_path.name)``.
            Useful for mixed-language collections (the radia-mcp lab
            corpus is roughly 50/50 English papers + Japanese
            textbooks side by side under public-safe curated corpus).

    Requires `pymupdf` (a.k.a. `fitz`).

    Ported from wjc9011/COMSOL_Multiphysics_MCP fork
    (src/knowledge/pdf_processor.py::PDFProcessor.process_pdf).
    """
    try:
        import fitz  # pymupdf
    except ImportError:
        raise ImportError(
            "PDF extraction requires pymupdf: pip install pymupdf")
    pdf_path = Path(pdf_path)

    # Resolve per-file language tag.
    if default_language:
        per_file_language: Optional[str] = default_language
    elif auto_detect_language:
        per_file_language = detect_filename_language(pdf_path.name)
    else:
        per_file_language = None

    chunks: list[dict] = []
    doc = fitz.open(str(pdf_path))
    stem = pdf_path.stem
    for page_idx in range(len(doc)):
        page = doc[page_idx]
        text = page.get_text() or ""
        ocr_used = False
        if ocr_fallback:
            has_text = bool(text.strip())
            readable = _text_readability(text) if has_text else 0.0
            # OCR when the text layer is GARBAGE (corrupt CID font) or is absent
            # on an image page; keep the OCR text only if it is actually cleaner.
            need_ocr = (has_text and readable < ocr_min_readability) or \
                       (not has_text and bool(page.get_images()))
            if need_ocr:
                try:
                    ocr_text = _ocr_pdf_page(page, ocr_languages, ocr_dpi, ocr_engine)
                    if (len(ocr_text.strip()) > 50
                            and _text_readability(ocr_text) > max(readable, ocr_min_readability)):
                        text = ocr_text
                        ocr_used = True
                except Exception:
                    pass  # OCR is best-effort; fall back to the extracted text
        if not text.strip():
            continue
        # Sliding-window chunks per page
        start = 0
        chunk_in_page = 0
        while start < len(text):
            end = min(start + chunk_size, len(text))
            chunk_text = text[start:end].strip()
            if len(chunk_text) > 50:   # skip near-empty chunks
                meta: dict = {
                    "source": pdf_path.name,
                    "stem": stem,
                    "page": page_idx + 1,
                }
                if per_file_language:
                    meta["language"] = per_file_language
                if ocr_used:
                    meta["ocr"] = True
                chunks.append({
                    "id": f"{stem}__p{page_idx + 1}__{chunk_in_page}",
                    "text": chunk_text,
                    "metadata": meta,
                })
                chunk_in_page += 1
            if end >= len(text):
                break
            start = end - overlap
    doc.close()
    return chunks
