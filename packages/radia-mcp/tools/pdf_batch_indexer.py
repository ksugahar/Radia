"""
PDF batch indexer for radia-mcp.

Walks a literature directory (e.g. W:\\03_文献・論文), extracts metadata
from each PDF (title, first-page abstract, suggested keywords), classifies
each PDF into a radia-mcp subpackage by keyword matching, and emits a
structured JSON catalog.

The JSON catalog is then ingested by `bibliography_index.py` knowledge
modules in each subpackage.

Usage:
    python pdf_batch_indexer.py \\
        --root "W:/03_文献・論文/00_電磁界解析" \\
        --output catalog.json \\
        --max-files 100   # batch size for one run
        --skip-existing   # don't re-index PDFs already in the catalog

Designed for incremental / batched processing — re-run with --skip-existing
to add new PDFs without re-processing the entire library.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

try:
    import fitz  # PyMuPDF
except ImportError:
    print("Error: PyMuPDF not installed. pip install pymupdf")
    raise


# ============================================================
# Classification rules: filename / abstract keyword → subpackage
# ============================================================
# Order matters: first match wins. Put specific topics before general.
CLASSIFICATION_RULES = [
    # --- CAE-AI / ML / optimization ---
    (r"physics[- ]informed.*neural|pinn", "pinn"),
    (r"physics[- ]informed.*gaussian|gp.+pde", "pinn"),
    (r"neural network|deep learning|machine learning|cnn|rnn|gnn", "pinn"),
    (r"bayesian opt|gaussian process|gp regression|optuna", "pinn"),
    (r"topology optim|shape optim|sensitivity analysis", "topology_optimization"),
    (r"shape derivative|topological derivative", "topology_optimization"),

    # --- MOR / model reduction ---
    (r"cauer ladder|cln method", "mor"),
    (r"model order reduction|reduced order model|prima|lanczos|arnoldi|krylov", "mor"),
    (r"pod method|proper orthogonal|reduced[- ]basis", "mor"),
    (r"hybrid twin|digital twin.*ml", "mor"),

    # --- Forces / Mechanics ---
    (r"electromagnetic force|maxwell stress|magnetic force|nodal force|local force", "differential_forms"),
    (r"magnetostrict|ferromagnetic force|electromechan", "differential_forms"),
    (r"eggshell method|arkkio|coulomb.+method", "differential_forms"),
    (r"force density|jacobian.+force|virtual work", "differential_forms"),

    # --- Differential forms / Whitney / FEEC / cohomology ---
    (r"whitney|edge element|differential form|exterior calculus|de rham", "differential_forms"),
    (r"hodge|kelvin transform|maxwell house|tonti diagram", "differential_forms"),
    (r"feec|finite element exterior calc|cochain", "differential_forms"),
    (r"tree.cotree|gauge.*fix|null.space", "differential_forms"),

    # --- Accelerator ---
    (r"accelerator|synchrotron|cyclotron|insertion device|undulator|wiggler", "accelerator"),
    (r"multipole|sextupole|quadrupole|end[- ]?pole|chamfer", "accelerator"),
    (r"rotating coil|harmonic analysis.+magnet", "accelerator"),

    # --- Application areas (build out as new subpackages added) ---
    (r"induction heating|workpiece|sibc|skin depth.+heating", "ih"),
    (r"peec|partial element|loop[- ]star|fasthenry", "peec"),
    (r"litz wire|stranded conductor|proximity loss", "peec"),
    (r"wireless power|wpt|inductive coupling|kerl coil", "peec"),
    # --- Motor / electric machine (rotating machines, IPM, SynRM, IM, SRM) ---
    (r"motor|machine|stator|rotor|ipm|synrm|reluctance.+machine", "motor"),
    (r"hameyer|henrotte|lange.+linearization", "motor"),
    (r"field.circuit coupling|moving band|sliding band|air.gap.+coupling", "motor"),
    (r"switched reluctance|wound.field|squirrel cage|park.+transform", "motor"),
    (r"cogging|back.?emf|torque ripple|d.q axis|incremental permeability", "motor"),
    # --- Static / DC / accelerator electromagnets (not rotating) ---
    (r"permanent magnet|halbach|magnetic levitation|maglev", "electromagnet"),

    # --- High-frequency / BEM / scattering ---
    (r"h.matrix|hacapk|fast multipole|fmm", "radia_ngsolve"),
    (r"moment method|method of moments|boundary integral|bem", "radia_ngsolve"),
    (r"scattering|microwave|antenna|radar|rcs", "radia_ngsolve"),
    (r"pml|perfectly matched|absorbing.+boundary", "radia_ngsolve"),
    (r"nonconfirming|non[- ]conformal", "radia_ngsolve"),

    # --- Magnetic materials / hysteresis ---
    (r"hysteresis|play model|preisach|jiles[- ]atherton|stoner[- ]wohlfarth", "electromagnet"),
    (r"iron loss|eddy current loss|core loss|joule loss", "ih"),
    (r"saturation|bh curve|nonlinear magneto", "electromagnet"),

    # --- Mesh / FEM ---
    (r"mesh generation|cubit|gmsh|netgen|tetra.+mesh|hex.+mesh", "cubit"),
    (r"finite element|fem formulation", "radia_ngsolve"),

    # --- Math / theory ---
    (r"topology|cohomology|homology|simplicial", "differential_forms"),
    (r"sobolev|hilbert|functional analysis", "differential_forms"),
    (r"fourier|spectral|wavelet", "radia_ngsolve"),

    # --- Specialized ---
    (r"nondestructive|ndt|defect detection|crack detection", "ndt"),  # future subpackage
    (r"metamaterial|negative index", "metamaterial"),  # future
    (r"transmission line|tem|tlm method", "peec"),
    (r"superconduct|hts|low.tc|nb3sn", "accelerator"),
    (r"plasma|tokamak|stellarator|fusion", "accelerator"),
    (r"nmr|mri|spectrometer", "electromagnet"),
    (r"meta[- ]heuristic|particle swarm|evolutionary|genetic algorithm", "pinn"),
    (r"surrogate|response surface|kriging", "pinn"),

    # --- Japanese-language broad EM literature (catch-all by language) ---
    (r"電磁界数値解析|解析積分公式|公式集", "radia_ngsolve"),
    (r"電磁界解析|電磁場解析|電磁気解析", "radia_ngsolve"),
    (r"電磁界|電磁場.+解析|electromagnetic.+analysis", "radia_ngsolve"),
    (r"四面体|六面体|デローニ|delaunay|メッシュ生成", "cubit"),
    (r"四角形要素|形状関数|shape function|形状最適化", "differential_forms"),
    (r"高精度化|高精度|精度向上|精度改善", "radia_ngsolve"),
    (r"インダクタ|inductor", "peec"),
    (r"フェライト|ferrite|磁心|core loss", "electromagnet"),
    (r"トランス|transformer|変圧器", "electromagnet"),
    (r"磁気回路|magnetic circuit|reluctance network", "electromagnet"),
    (r"修士|博士|学位論文|thesis|dissertation", "radia_ngsolve"),
    (r"electromagnetic compatibility|emc|emi|electromagnetic interfer", "radia_ngsolve"),
    (r"surface impedance|nonlinear effective surface", "ih"),
    (r"レイル|rail|magnetic levitation|maglev", "electromagnet"),
    (r"power converter|inverter|switching|電力変換", "peec"),
    (r"flyback|bridgeless|dc.dc converter", "peec"),
    (r"plasma|fusion|tokamak", "accelerator"),
    (r"高速.+計算|高速.+技術|fast.+algorithm|gpu.+computing", "radia_ngsolve"),
    (r"先進.+解析|先進電磁|高度化", "radia_ngsolve"),
    (r"committee report|専門委員会|研究会|tech.*report", "radia_ngsolve"),
    (r"ieee trans.*mag|ieee transactions on magnetics", "radia_ngsolve"),
    (r"arxiv|preprint.*math|preprint.*na", "radia_ngsolve"),
    (r"journal of research|nbs.*standard|bureau.*standard", "radia_ngsolve"),

    # --- Final round: more specific rules to reduce _misc ---
    (r"magnetic resonance|nmr|mri|spectrometer|relax.*time", "nmr_mri"),
    (r"demagnetiz|demagnetisation|stoner.*demag", "electromagnet"),
    (r"medical|particle.*therapy|粒子線治療|cancer.*treatment|がん", "accelerator"),
    (r"team problem|problem[ -]?[0-9]+|felix", "differential_forms"),
    (r"演習|講義|lecture|tutorial|exercise|note|教科書", "_textbooks"),
    (r"多芯線|stranded.*wire|hybrid.*coil|wound.*coil", "peec"),
    (r"emp.*interact|electromagnetic.*pulse|lightning|落雷", "radia_ngsolve"),
    (r"shielding|遮蔽|シールド", "radia_ngsolve"),
    (r"benchmark|comparison|validation|verification", "radia_ngsolve"),
    (r"hall.*effect|magnetoresist|spin.*orbit|spintronic", "electromagnet"),
    (r"levitation|magnetic.*bearing|bearingless", "electromagnet"),
    (r"saito|kindai|kindai.*univ|学長|email|invoice|order", "_admin"),

    # --- Final catch-all fallback ---
    (r".*", "_misc"),
]


def classify(text: str) -> str:
    """Return the radia-mcp subpackage name based on text content."""
    text_lower = text.lower()
    for pattern, subpackage in CLASSIFICATION_RULES:
        if re.search(pattern, text_lower):
            return subpackage
    return "_misc"


def extract_metadata(pdf_path: Path) -> dict[str, Any]:
    """Extract title, authors, abstract from a PDF.  Returns a dict with
    fields { 'title', 'authors', 'abstract_excerpt', 'page_count',
              'filename', 'is_scanned' }."""
    doc = fitz.open(pdf_path)
    page_count = len(doc)

    # Get the first 2-3 pages of text for title/abstract extraction
    first_text = ""
    for i in range(min(3, page_count)):
        try:
            first_text += doc[i].get_text("text") + "\n"
        except Exception:
            pass
    doc.close()

    # Detect scanned (1 font, lots of images)
    is_scanned = False
    if len(first_text.strip()) < 100 and page_count > 0:
        is_scanned = True

    # Heuristic title extraction: longest line in first 30 lines
    lines = [l.strip() for l in first_text.split("\n") if l.strip()][:30]
    title_candidates = [l for l in lines if 10 < len(l) < 200]
    title = title_candidates[0] if title_candidates else pdf_path.stem

    # Heuristic abstract: first 500-1000 chars after "Abstract"
    abstract_excerpt = ""
    m = re.search(r"(?i)abstract\s*[—-]?\s*(.{50,1500})", first_text)
    if m:
        abstract_excerpt = m.group(1).strip()[:1000]
    else:
        # Fallback: first 500 chars after the title
        idx = first_text.find(title)
        if idx >= 0:
            abstract_excerpt = first_text[idx + len(title):idx + len(title) + 800].strip()

    # Heuristic authors: line with multiple commas + capitalized words,
    # right after title
    authors = ""
    title_idx = first_text.find(title)
    if title_idx >= 0:
        after_title = first_text[title_idx + len(title):title_idx + len(title) + 500]
        author_line_m = re.search(r"^\s*([A-Z][^\n]{5,200})$", after_title, re.MULTILINE)
        if author_line_m:
            authors = author_line_m.group(1).strip()

    return {
        "filename": pdf_path.name,
        "relative_path": str(pdf_path),
        "page_count": page_count,
        "is_scanned": is_scanned,
        "title": title,
        "authors": authors,
        "abstract_excerpt": abstract_excerpt,
        "subpackage": classify(title + " " + abstract_excerpt + " " + pdf_path.name),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, help="Literature root directory")
    parser.add_argument("--output", required=True, help="Output JSON catalog")
    parser.add_argument("--max-files", type=int, default=0,
                        help="Process at most this many new files (0 = no limit)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip PDFs already in the output catalog")
    args = parser.parse_args()

    root = Path(args.root)
    output_path = Path(args.output)

    # Load existing catalog
    existing: dict[str, dict[str, Any]] = {}
    if args.skip_existing and output_path.exists():
        try:
            with output_path.open("r", encoding="utf-8") as f:
                existing = {e["relative_path"]: e for e in json.load(f)}
            print(f"Loaded {len(existing)} existing entries from {output_path}")
        except Exception as e:
            print(f"Warning: could not load existing catalog ({e})")

    # Walk PDFs
    pdfs = list(root.rglob("*.pdf"))
    print(f"Found {len(pdfs)} PDFs under {root}")

    def safe_print(s: str) -> None:
        """Print, falling back to ASCII when terminal is cp932."""
        try:
            print(s)
        except UnicodeEncodeError:
            print(s.encode("ascii", errors="replace").decode("ascii"))

    def write_catalog(all_entries: list[dict[str, Any]]) -> None:
        with output_path.open("w", encoding="utf-8") as f:
            json.dump(all_entries, f, indent=2, ensure_ascii=False)

    new_entries = []
    n_done = 0
    n_skipped = 0
    n_failed = 0

    for pdf_path in pdfs:
        rel = str(pdf_path)
        if rel in existing:
            n_skipped += 1
            continue
        if args.max_files and n_done >= args.max_files:
            break
        try:
            entry = extract_metadata(pdf_path)
            new_entries.append(entry)
            n_done += 1
            if n_done % 10 == 0:
                safe_print(f"  [{n_done}] {entry['subpackage']:20s} | "
                           f"{entry['filename'][:60]}")
            # Incremental save every 100 entries
            if n_done % 100 == 0:
                write_catalog(list(existing.values()) + new_entries)
        except Exception as e:
            n_failed += 1
            safe_print(f"  FAIL: {pdf_path.name}: {e}")

    # Final write
    all_entries = list(existing.values()) + new_entries
    write_catalog(all_entries)

    print()
    print(f"Done: {n_done} new, {n_skipped} skipped, {n_failed} failed")
    print(f"Total entries: {len(all_entries)} → {output_path}")
    print()

    # Per-subpackage counts
    counts: dict[str, int] = {}
    for e in all_entries:
        counts[e["subpackage"]] = counts.get(e["subpackage"], 0) + 1
    print("Per-subpackage counts:")
    for sub, n in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {n:5d}  {sub}")


if __name__ == "__main__":
    main()
