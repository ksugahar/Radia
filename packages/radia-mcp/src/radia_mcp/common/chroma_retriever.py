"""ChromaDB + sentence-transformers RAG infrastructure.

Adapted from wjc9011/COMSOL_Multiphysics_MCP src/knowledge/retriever.py
which uses the same stack to index 75 COMSOL official PDFs (~500MB)
into a persistent ChromaDB store, then expose semantic search to LLMs.

For radia_mcp the target corpus is W:/03_文献・論文/00_電磁界解析/
(3,889 files / 51 GB) plus the freshly-organized 04_機械学習と最適化/.

Optional dependencies (graceful degradation if missing):
  - chromadb
  - sentence-transformers (defaults to all-MiniLM-L6-v2, ~80 MB)
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
         "metadata": {"source": "MaxwellBook.pdf", "page": 1}},
        ...
    ])

    # Semantic search.
    hits = rag.search("rotational core loss tester", n_results=5)
    # → [{"text": ..., "metadata": {...}, "score": 0.87}, ...]

Design choice: this module does NOT depend on radia_mcp.literature_index
— it is a generic utility usable by any subpackage. literature_index
adopts it as the engine for `literature_semantic_search`.
"""

from __future__ import annotations
import logging
from pathlib import Path
from typing import Any, Optional, Sequence

logger = logging.getLogger(__name__)

DEFAULT_EMBEDDING_MODEL = "all-MiniLM-L6-v2"


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
    ) -> list[dict]:
        """Semantic search.

        Args:
            query: free-text query (embedded with same model as docs)
            n_results: number of top hits to return
            where: optional ChromaDB metadata filter, e.g.
                   {"source": {"$eq": "MaxwellBook.pdf"}}

        Returns:
            List of dicts: [{"id", "text", "metadata", "score"}, ...]
            (score in [0, 1] — higher is more similar)
        """
        if not self.is_initialized and not self.initialize():
            return []
        try:
            res = self._collection.query(
                query_texts=[query],
                n_results=n_results,
                where=where,
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

def extract_pdf_chunks(
    pdf_path: str | Path,
    chunk_size: int = 1000,
    overlap: int = 100,
) -> list[dict]:
    """Extract text chunks from a PDF, ready to feed `add_chunks`.

    Simple page-based chunking with overlap. Returns chunks with:
        - id = <stem>__p<page>__<chunk_in_page>
        - text = the chunk content
        - metadata = {"source": filename, "page": int, "stem": stem}

    Requires `pymupdf` (a.k.a. `fitz`).
    """
    try:
        import fitz  # pymupdf
    except ImportError:
        raise ImportError(
            "PDF extraction requires pymupdf: pip install pymupdf")
    pdf_path = Path(pdf_path)
    chunks: list[dict] = []
    doc = fitz.open(str(pdf_path))
    stem = pdf_path.stem
    for page_idx in range(len(doc)):
        page = doc[page_idx]
        text = page.get_text() or ""
        if not text.strip():
            continue
        # Sliding-window chunks per page
        start = 0
        chunk_in_page = 0
        while start < len(text):
            end = min(start + chunk_size, len(text))
            chunk_text = text[start:end].strip()
            if len(chunk_text) > 50:   # skip near-empty chunks
                chunks.append({
                    "id": f"{stem}__p{page_idx + 1}__{chunk_in_page}",
                    "text": chunk_text,
                    "metadata": {
                        "source": pdf_path.name,
                        "stem": stem,
                        "page": page_idx + 1,
                    },
                })
                chunk_in_page += 1
            if end >= len(text):
                break
            start = end - overlap
    doc.close()
    return chunks
