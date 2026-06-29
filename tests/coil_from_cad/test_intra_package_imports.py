"""Regression: intra-package imports must be package-qualified, not bare.

Several library modules under ``src/radia/`` historically imported their
sibling modules with a *bare* top-level name, e.g.::

    # coil_from_step.py (inside to_coil_builder)
    from coil_builder import CoilBuilder

That only resolves when ``src/radia`` itself is on ``sys.path`` -- the
legacy dev convention.  For a normal ``pip install radia`` / ``import
radia`` environment (only the package *parent* ``src`` -- or
site-packages -- is on the path), the bare form raises
``ModuleNotFoundError: No module named 'coil_builder'`` the moment the
function is called.

Concretely, ``filaments_from_step(...)`` -> ``_filaments_via_coil_builder``
-> ``to_coil_builder`` hit ``coil_from_step.py`` line 737 and blew up for
installed users.  The fix is to use ``from radia.X import Y`` everywhere
(the convention documented in ``tests/conftest.py``: src/radia is
deliberately NOT on sys.path, "All tests now use ``from radia.X import
Y`` exclusively").

These tests lock that convention:

1. ``test_fixed_modules_import_only_src_on_path`` spawns a child with
   ONLY ``src`` on the path, asserts a bare sibling import fails (proving
   the config), then imports every previously-broken module via
   ``radia.X`` -- a bare-sibling ``ModuleNotFoundError`` is a hard fail
   (optional-dep absence, e.g. ngsolve, is a skip).

2. ``test_filaments_from_step_reaches_coil_builder`` runs the real
   ``filaments_from_step`` end-to-end on the keiko fixture STEP, which
   exercises ``to_coil_builder`` -- catching a re-introduced bare import
   inside the function body that the import-time test would miss.  Runs
   under the conftest only-src config.
"""
from __future__ import annotations

import os
import subprocess
import sys
import textwrap

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(REPO_ROOT, "src")
FIXTURE_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
KEIKO_STEP = os.path.join(FIXTURE_DIR, "keiko_outsideline.step")

# Modules that previously contained an unguarded bare sibling import.
# (Guarded `try: from .X except ImportError: from X` modules and
# `if __name__ == "__main__"` demo blocks are intentionally excluded --
# they already work for installed users.)
FIXED_MODULES = [
    "radia.coil_from_step",      # from coil_builder import CoilBuilder
    "radia.kelvin_source",       # from biot_savart import ...
    "radia.kelvin_solver",       # from kelvin_material import ...
    "radia.kelvin_validate",     # from kelvin_geometry / kelvin_solver import ...
    "radia.ngsbem_interface",    # from peec_topology import ...
    "radia.streamfunction_volume",  # from stream_function import ...
    "radia.em_material",         # from esim_cell_problem import ...
    "radia.peec_matrices",       # from peec_matrices import ... (C++ stub)
]


def test_fixed_modules_import_only_src_on_path():
    """Every fixed module must import with ONLY `src` on sys.path.

    Run in a subprocess so the parent's already-imported modules cannot
    mask a regression, and so we fully control sys.path.
    """
    child = textwrap.dedent(
        """
        import os, sys, importlib

        SRC = sys.argv[1]
        # Keep ONLY `src` of the repo on the path (drop any stray src/radia
        # the parent env may have injected); the editable/site-packages
        # `radia` install still resolves the package itself.
        radia_pkg = os.path.normcase(os.path.join(SRC, "radia"))
        sys.path[:] = [p for p in sys.path
                       if os.path.normcase(os.path.abspath(p or ".")) != radia_pkg]
        if SRC not in sys.path:
            sys.path.insert(0, SRC)

        mkl = os.path.join(sys.prefix, "Library", "bin")
        if os.path.isdir(mkl) and hasattr(os, "add_dll_directory"):
            os.add_dll_directory(mkl)

        import radia  # noqa: F401

        # Sanity: the config genuinely has src/radia OFF the path, so a
        # bare sibling import MUST fail.  If it succeeds, the test is not
        # actually exercising the installed-user scenario.
        try:
            import coil_builder  # noqa: F401
            print("CONFIG-ERROR: bare 'import coil_builder' unexpectedly succeeded")
            sys.exit(3)
        except ModuleNotFoundError:
            pass

        # The primary fix: the CoilBuilder symbol used by to_coil_builder.
        from radia.coil_builder import CoilBuilder  # noqa: F401

        SIBLINGS = {
            "coil_builder", "biot_savart", "kelvin_material", "kelvin_geometry",
            "kelvin_solver", "peec_topology", "stream_function",
            "esim_cell_problem", "peec_matrices", "peec_bundle",
        }
        mods = sys.argv[2:]
        failures = []
        for m in mods:
            try:
                importlib.import_module(m)
            except ModuleNotFoundError as e:
                missing = (e.name or "").split(".")[-1]
                if missing in SIBLINGS:
                    failures.append(f"{m}: bare-sibling ModuleNotFoundError -> {e}")
                # else: optional dependency (ngsolve, netgen, ...) -> ignore
            except Exception:
                # A non-import error (e.g. heavy module-level work) is not an
                # import-resolution regression; ignore here.
                pass

        if failures:
            print("IMPORT-REGRESSION:\\n" + "\\n".join(failures))
            sys.exit(1)
        print("OK")
        sys.exit(0)
        """
    )
    env = dict(os.environ)
    env["PYTHONPATH"] = SRC  # only the package parent, never src/radia
    proc = subprocess.run(
        [sys.executable, "-c", child, SRC, *FIXED_MODULES],
        capture_output=True, text=True, env=env, cwd=REPO_ROOT,
    )
    out = (proc.stdout or "") + (proc.stderr or "")
    assert proc.returncode == 0, (
        f"intra-package import regression (rc={proc.returncode}):\n{out}"
    )


@pytest.mark.slow
def test_filaments_from_step_reaches_coil_builder():
    """End-to-end: filaments_from_step -> to_coil_builder -> coil_builder.

    This is the exact path from the bug report.  Runs under the conftest
    only-src config (src/radia is NOT on sys.path), so a bare
    `from coil_builder import CoilBuilder` would raise ModuleNotFoundError
    here.
    """
    pytest.importorskip("build123d")
    pytest.importorskip("netgen")
    assert os.path.isfile(KEIKO_STEP), KEIKO_STEP

    from radia.coil_from_cad import filaments_from_step

    res = filaments_from_step(
        KEIKO_STEP, nwinc=2, nhinc=2, cad_units_per_meter=1000.0,
    )
    # CoilBuilder path returns a dict including the built coil_builder.
    assert "coil_builder" in res
    assert res.get("filament_paths"), "no filament paths reconstructed"
