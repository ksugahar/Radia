"""Check tracked index text blobs for CRLF/mixed endings, honoring attributes.

Git's index EOL classifier detects binary data; explicit -text/binary evidence
is preserved. Check the index, not an autocrlf working copy. No files/config are
modified. UTF-8/LF is the authoring policy; this guard checks line endings only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


def check_index(root: Path) -> dict:
    result = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--eol", "-z"],
        check=True, capture_output=True)
    violations = []
    checked = preserved = binary = 0
    for record in result.stdout.split(b"\0"):
        if not record:
            continue
        metadata, path = record.split(b"\t", 1)
        fields = metadata.decode("ascii").split()
        indexed = fields[0].removeprefix("i/")
        # Attribute override wins over automatic text classification. This
        # includes historical hash-bound JSON and source-form evidence.
        if "attr/-text" in fields:
            preserved += 1
            continue
        if indexed == "-text":
            binary += 1
            continue
        checked += 1
        if indexed in {"crlf", "mixed"}:
            violations.append({"path": path.decode("utf-8"), "index_eol": indexed})
    return {"schema": "radia.tracked-text-eol.v1", "passed": not violations,
            "checked_text_files": checked, "preserved_attribute_files": preserved,
            "auto_binary_files": binary, "violations": violations}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args(argv)
    result = check_index(args.root.resolve())
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
