#!python
"""Compatibility startup shim for legacy .cubit entries.

New installs do not play this file directly.  ``radia.install_panels``
generates a startup shim under ProgramData or LocalAppData with the
absolute ``register_toolbar.py`` path baked in.  This file remains only
for older .cubit blocks that still point into the package.
"""

import glob
import os
import sys
import traceback


def _add_cubit_site_packages():
    cubit_bin = (
        sys._cubit_bin if hasattr(sys, "_cubit_bin")
        else os.path.dirname(
            os.path.abspath(os.path.join(os.path.dirname(sys.executable), "cubit.py"))
        )
    )
    candidates = glob.glob(os.path.join(cubit_bin, "python*", "lib", "site-packages"))
    candidates += glob.glob(
        os.path.join(cubit_bin, "python*", "lib", "python*", "site-packages")
    )
    if candidates and candidates[0] not in sys.path:
        sys.path.insert(0, candidates[0])


def _register_toolbar_path():
    if "__file__" not in globals():
        raise RuntimeError(
            "Cubit did not expose __file__ for startup.py; re-run "
            "cubit-plugin-install to regenerate the ProgramData startup shim."
        )
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "register_toolbar.py")
    if not os.path.isfile(path):
        raise RuntimeError(f"register_toolbar.py not found next to startup.py: {path}")
    return path


try:
    _add_cubit_site_packages()
    register_path = _register_toolbar_path()
    __file__ = register_path
    with open(register_path, encoding="utf-8") as f:
        exec(f.read())
except Exception:
    traceback.print_exc()
