"""Read citation decisions already made by TeX, without executing aux commands."""
import hashlib
from pathlib import Path
import re

from ..paper_writing._tex_lex import mask_tex_noncode


def read_compiled_aux(path: Path):
    keys, styles, snapshots = [], [], []
    seen, active = set(), set()
    root = path.resolve().parent

    def visit(current, depth=0):
        current = current.resolve()
        if current in active:
            raise ValueError("cyclic compiled aux chain")
        if current in seen:
            return
        if depth > 16:
            raise ValueError("compiled aux chain exceeds depth limit")
        if current.suffix.casefold() != ".aux":
            raise ValueError("compiled input must be an aux file")
        raw = current.read_bytes()
        text = mask_tex_noncode(raw.decode("utf-8", errors="strict"))
        snapshots.append({"path": str(current), "sha256": hashlib.sha256(raw).hexdigest()})
        active.add(current)
        records = list(re.finditer(r"(?<!\\)\\(citation|bibstyle|@input)\s*\{([^{}]*)\}", text))
        heads = list(re.finditer(r"(?<!\\)\\(?:citation|bibstyle|@input)(?![A-Za-z@])", text))
        if len(records) != len(heads):
            raise ValueError("malformed compiled aux citation/input/style record")
        for match in records:
            command, value = match.groups()
            if command == "@input":
                visit(root / value.strip(), depth + 1)
            elif command == "bibstyle":
                styles.append(value.strip())
            else:
                for key in value.split(","):
                    key = key.strip()
                    if key and key not in keys:
                        keys.append(key)
        active.remove(current)
        seen.add(current)

    visit(path)
    if len(set(styles)) > 1:
        raise ValueError("conflicting bibliography styles in compiled aux")
    if not keys:
        raise ValueError("compiled aux contains no citation records (biber/bcf is not BibTeX aux)")
    return keys, styles[0] if styles else "", snapshots
