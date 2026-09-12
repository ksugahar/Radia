"""One model/font operation followed by process teardown; disposable CI only.

An external observer must monitor the font host and then launch the fixed EXE.
Successful child exit is not a claim that the host survived.
"""
import argparse
import os
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("load-only", "metrics-ab", "svg-ab", "svg-math"))
    args = parser.parse_args()
    if (os.environ.get("GITHUB_ACTIONS") != "true" or
            os.environ.get("EQNEDIT64_ISOLATED_TEST_SESSION") != "1"):
        print("FAIL: disposable diagnostic CI required", flush=True)
        return 90
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "build"))
    import eqnedit_core as core
    if not core.math_font_loaded():
        return 1
    print(f"OPERATION_BEGIN {args.case}", flush=True)
    if args.case == "metrics-ab":
        equation = core.Equation()
        equation.insert_text("ab")
        width, height, _baseline = equation.metrics()
        if width <= 0 or height <= 0:
            return 2
    elif args.case in ("svg-ab", "svg-math"):
        source = "ab" if args.case == "svg-ab" else r"\frac{a}{b}+\sqrt{x^{2}}+\sum_{i=1}^{n}x_i"
        if "<svg" not in core.tex_to_svg(source):
            return 3
    print(f"OPERATION_END {args.case}", flush=True)
    print("DIAGNOSTIC_ONLY: operation completed; external host observation required", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
