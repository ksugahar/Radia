"""One model/font operation followed by process teardown; disposable CI only.

An external observer must monitor the font host and then launch the fixed EXE.
Successful child exit is not a claim that the host survived.
"""
import argparse
import importlib.util
import os
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("load-only", "metrics-ab", "svg-ab", "svg-math",
                                        "edit-400", "edit-1000", "edit-1525", "edit-all",
                                        "deep-once", "deep-sequence"))
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
    elif args.case.startswith("deep-"):
        # Reproduce the depth fixture without any unrelated editing suite.
        depths = (200,) if args.case == "deep-once" else range(1, 201)
        for depth in depths:
            equation = core.Equation()
            source = r"\sqrt{" * depth + "}" * depth
            if not equation.load_latex(source):
                return 4
            equation.metrics()
        if "<svg" not in equation.svg():
            return 5
    elif args.case.startswith("edit-"):
        # A diagnostic prefix of the pinned source, not a shortened acceptance test.
        path = Path(__file__).with_name("test_edit.py")
        spec = importlib.util.spec_from_file_location("font_probe_edit", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        limit = None if args.case == "edit-all" else int(args.case.split("-")[1])

        class PrefixComplete(Exception):
            pass

        def trace(frame, event, _arg):
            if event == "line" and frame.f_code is module.main.__code__:
                if frame.f_lineno >= limit:
                    print(f"EDIT_PREFIX_STOP before_line={frame.f_lineno}", flush=True)
                    raise PrefixComplete()
            return trace

        try:
            if limit is not None:
                sys.settrace(trace)
            result = module.main()
            if result not in (None, 0):
                return int(result)
        except PrefixComplete:
            pass
        finally:
            sys.settrace(None)
    print(f"OPERATION_END {args.case}", flush=True)
    print("DIAGNOSTIC_ONLY: operation completed; external host observation required", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
