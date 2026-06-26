"""Skip the Cubit validation suite unless explicitly opted in.

The validation_test/cubit/* modules exercise Coreform Cubit 2025.12's APREPRO
export plugin and may import `cubit` / call `cubit.init([...])`. Outside a
Cubit-enabled interpreter that import either fails (ModuleNotFound) or --
worse -- `cubit.init()` blocks on the Coreform engine/license, hanging pytest
collection. So these tests must be run only in a deliberate Cubit context,
gated by an explicit env var.

This conftest makes a bare `pytest validation_test/` skip the Cubit subtree
cleanly -- without importing cubit -- instead of erroring or hanging.

To RUN the Cubit tests (inside a Cubit-enabled Python), set:  RADIA_RUN_CUBIT_TESTS=1
"""
import os

if not os.environ.get("RADIA_RUN_CUBIT_TESTS"):
    # Skip collecting the whole directory.  We must NOT import cubit here to decide --
    # `import cubit` / `cubit.init()` can block, which is exactly the hang we avoid.
    collect_ignore_glob = ["*"]
