r"""Convert a directory of .eqn files and score the result.

The unit tests say the editor still behaves; they say nothing about the
documents that already worked, because the MTEF passes are shared -- every
document goes through them.  This is the check that does: run a corpus before
and after a change, and look at what actually moved.

    python convert_corpus.py <dir> --out before.json
    ... make the change, rebuild ...
    python convert_corpus.py <dir> --out after.json
    python convert_corpus.py --diff before.json after.json
    python convert_corpus.py --health after.json

The diff normalises away the size markers, because a change to those touches
every file at once and would bury the structural differences that can be
regressions.

No corpus is named here.  The lab's own .eqn collection is private, and only
the derived LaTeX is written, to wherever --out says.
"""

from __future__ import annotations

import argparse
import io
import json
import os
import re
import sys
from collections import Counter
from pathlib import Path

# Run against the tree this file lives in, before any installed copy.
for _p in Path(__file__).resolve().parents:
    if (_p / "src" / "radia").exists():
        sys.path.insert(0, str(_p / "src"))
        break

B = chr(92)

# Shapes that cannot be right: an empty slot, a bracket round nothing, a
# command that reached the unknown-command path.  A document carrying one is a
# document still to fix.
CHECKS = [
    # \left, ONE delimiter -- a command word or a single character -- then
    # nothing but space before \right.  Reading it as "up to twelve of
    # anything" counted \left( 1+x_{ij} \right) as empty, and nineteen of the
    # twenty documents it flagged were fine.
    ("empty fence", re.compile(re.escape(B + "left")
                               + r"(?:" + re.escape(B) + r"[a-zA-Z]+\s?|.)"
                               + r"\s*" + re.escape(B + "right"))),
    ("empty fraction slot", re.compile(re.escape(B + "dfrac{}")
                                       + "|" + re.escape(B + "frac{}"))),
    ("unknown command", re.compile(re.escape(B + "text{") + re.escape(B))),
    ("stray style marker", re.compile(re.escape(B + "scriptstyle"))),
    ("empty matrix row", re.compile(re.escape(B + "begin{matrix}")
                                    + r"\s*" + re.escape(B + B))),
    ("replacement char", re.compile("�")),
]

STYLE = re.compile(re.escape(B)
                   + r"(displaystyle|textstyle|scriptstyle|scriptscriptstyle)\s*")


def convert(root, out_path):
    import radia.equation as eq

    files = []
    for base, _dirs, names in os.walk(root):
        for nm in names:
            if nm.lower().endswith(".eqn"):
                files.append(os.path.join(base, nm))
    files.sort()

    result = {}
    errors = 0
    for p in files:
        key = os.path.relpath(p, root).replace(os.sep, "/")
        try:
            result[key] = {"latex": eq.mtef_to_latex(eq.read_eqn(p))}
        except Exception as exc:                      # noqa: BLE001
            result[key] = {"error": str(exc)}
            errors += 1
    io.open(out_path, "w", encoding="utf-8").write(
        json.dumps(result, ensure_ascii=False, indent=1))
    print("files: %d  errors: %d  -> %s" % (len(files), errors, out_path))
    return 1 if errors else 0


def text_of(entry):
    return entry.get("latex", entry.get("error", ""))


def diff(a_path, b_path, show):
    a = json.load(io.open(a_path, encoding="utf-8"))
    b = json.load(io.open(b_path, encoding="utf-8"))

    def norm(s):
        return STYLE.sub("", s).replace(" ", "")

    same = style_only = 0
    changed = []
    for k in sorted(a):
        if k not in b:
            changed.append(k)
            continue
        x, y = text_of(a[k]), text_of(b[k])
        if x == y:
            same += 1
        elif norm(x) == norm(y):
            style_only += 1
        else:
            changed.append(k)

    print("files                    :", len(a))
    print("identical                :", same)
    print("size-marker/spacing only :", style_only)
    print("structurally different   :", len(changed))
    for k in changed[:show]:
        print("=" * 72)
        print(k)
        print("  before:", text_of(a[k])[:260])
        print("  after :", text_of(b.get(k, {}))[:260])
    if len(changed) > show:
        print("... and", len(changed) - show, "more")
    return 0


def health(path):
    data = json.load(io.open(path, encoding="utf-8"))
    counts = Counter()
    examples = {}
    clean = 0
    for k, v in sorted(data.items()):
        tex = text_of(v)
        hits = [name for name, rx in CHECKS if rx.search(tex)]
        if not hits:
            clean += 1
        for h in hits:
            counts[h] += 1
            examples.setdefault(h, k)
    total = max(1, len(data))
    print("documents        :", len(data))
    print("no defect marker : %d (%.1f%%)" % (clean, 100.0 * clean / total))
    for name, c in counts.most_common():
        print("  %-22s %4d   e.g. %s" % (name, c, examples[name]))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("root", nargs="?", help="directory of .eqn files")
    ap.add_argument("--out", help="where to write the converted LaTeX")
    ap.add_argument("--diff", nargs=2, metavar=("BEFORE", "AFTER"))
    ap.add_argument("--health", metavar="JSON")
    ap.add_argument("--show", type=int, default=10,
                    help="how many structural differences to print")
    args = ap.parse_args(argv)

    if args.diff:
        return diff(args.diff[0], args.diff[1], args.show)
    if args.health:
        return health(args.health)
    if not args.root or not args.out:
        ap.error("give a directory and --out, or --diff, or --health")
    return convert(args.root, args.out)


if __name__ == "__main__":
    sys.exit(main())
