"""Run the headless Eqnedit64 model suites in one Python process.

Each process loads the embedded math font privately.  Keeping these suites in
one host avoids needless register/unregister churn while preserving the same
individual ``main()`` entry points for focused debugging.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import sys


TESTS = (
    "test_edit.py",
    "test_symbols.py",
    "test_operations_fuzz.py",
    "test_tex_fuzz.py",
    "test_layout.py",
    "test_usability_trace.py",
)


def main() -> int:
    if os.environ.get("EQNEDIT64_ISOLATED_TEST_SESSION") != "1":
        print(
            "FAIL  model tests use a private math font; run them in a "
            "disposable CI/VM/user session with "
            "EQNEDIT64_ISOLATED_TEST_SESSION=1"
        )
        return 90
    test_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(test_dir.parent / "build"))
    from eqnedit_core import math_font_loaded

    if not math_font_loaded():
        print("FAIL  embedded math font unavailable before model tests")
        return 1
    for index, name in enumerate(TESTS):
        path = test_dir / name
        spec = importlib.util.spec_from_file_location(f"eqnedit64_model_{index}", path)
        if spec is None or spec.loader is None:
            print(f"FAIL  cannot load {path}")
            return 1
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        result = module.main()
        if result not in (None, 0):
            print(f"FAIL  {name} returned {result}")
            return int(result)
    print(f"PASS: {len(TESTS)} headless model suites in one process")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
