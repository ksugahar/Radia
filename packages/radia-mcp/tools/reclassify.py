"""
Re-classify an existing catalog.json using the latest CLASSIFICATION_RULES
from pdf_batch_indexer.  Useful after improving the classifier.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from pdf_batch_indexer import classify  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", required=True)
    args = parser.parse_args()

    cat_path = Path(args.catalog)
    with cat_path.open("r", encoding="utf-8") as f:
        entries = json.load(f)

    changed = 0
    for e in entries:
        text = " ".join([
            str(e.get("title", "")),
            str(e.get("authors", "")),
            str(e.get("abstract_excerpt", "")),
            str(e.get("filename", "")),
        ])
        new_sub = classify(text)
        if new_sub != e.get("subpackage"):
            e["subpackage"] = new_sub
            changed += 1

    print(f"Re-classified {changed} of {len(entries)} entries.")

    with cat_path.open("w", encoding="utf-8") as f:
        json.dump(entries, f, indent=2, ensure_ascii=False)

    c = Counter(e["subpackage"] for e in entries)
    print()
    print("New distribution:")
    for sub, n in c.most_common():
        print(f"  {n:5d}  {sub}")


if __name__ == "__main__":
    main()
