"""Standalone GMSH viewer launcher used by the Radia analysis windows.

The IH/EM/PCB panels invoke this script as a detached subprocess after a
Run completes successfully::

    python open_gmsh.py <path-to-geo-or-msh>

Why a separate script and not ``python -c "..."``?

    - Quoting ``r'<path>'`` inside a ``-c`` string is fragile, especially
      when the user is on a UNC share (``//192.168.11.100/work/...``)
      where forward / backslash conversions can break the literal.
    - ``-c`` failures vanish silently because the panel runs the launcher
      with ``CREATE_NO_WINDOW``.  A real script can wrap the gmsh calls
      in ``try / except`` and append the traceback to the shared panel
      debug log so the failure is visible to the user (and to Claude).

The script is intentionally minimal: ``gmsh.initialize(['-noconfig'])``
skips ``%APPDATA%\\gmsh-options`` so only the exact companion option file
controls display settings. Prefer ``case.geo`` + ``case.geo.opt`` for the
normal launch path; ``case.msh`` + ``case.msh.opt`` is raw inspection.
"""

from __future__ import annotations

import os
import sys


def _import_panel_log():
    """Best-effort ``panel_log`` import.

    The script lives in ``src/radia/panels/`` next to ``panel_log.py``,
    so the directory of ``__file__`` should always be on ``sys.path``.
    Falls back to no-op shims if the import fails for any reason --
    we never want a logging failure to suppress GMSH itself.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    try:
        from panel_log import (init_panel_log, panel_log,
                               panel_log_exception)
        init_panel_log("gmsh-launch", truncate=False, banner=False)
        return panel_log, panel_log_exception
    except Exception:
        return (lambda _msg: None), (lambda _prefix="": None)


def main(argv):
    panel_log, panel_log_exception = _import_panel_log()

    if len(argv) < 2:
        panel_log("open_gmsh.py: missing .geo/.msh argument")
        return 2

    input_path = argv[1]
    panel_log(f"open_gmsh.py: input={input_path}")

    if not os.path.isfile(input_path):
        panel_log(f"open_gmsh.py: ERROR file does not exist: {input_path}")
        return 3

    open_path = input_path
    stem, ext = os.path.splitext(input_path)
    sibling_merge_stem = stem if ext.lower() == ".msh" else ""
    if ext.lower() == ".msh":
        geo_path = stem + ".geo"
        if os.path.isfile(geo_path):
            open_path = geo_path
            panel_log(f"open_gmsh.py: using geo companion={geo_path}")

    try:
        import gmsh
    except ImportError:
        panel_log_exception("open_gmsh.py: gmsh import failed")
        return 4

    try:
        gmsh.initialize(["-noconfig"])
        gmsh.open(open_path)
        panel_log("open_gmsh.py: gmsh.open OK")

        # Auto-merge sibling overlay files written by calc_*.py:
        #   <stem>_filaments.msh  : PEEC filament 1D lines + |I|
        #                            (calc_inductance.py --coil-solver peec)
        #   <stem>_wp_surface.msh : workpiece surface mesh from wp BEM solve
        #   <stem>_wp.msh         : ditto, alternate naming
        #   <stem>_surface.msh    : generic surface mesh overlay
        # This lets a single Open GMSH click show all relevant overlays
        # (B field on bbox + filament currents + wp surface) in one
        # window without the user having to "File -> Merge..." manually.
        if sibling_merge_stem or os.path.splitext(open_path)[1].lower() != ".geo":
            stem = sibling_merge_stem or os.path.splitext(open_path)[0]
            merge_suffixes = ("_filaments.msh", "_wp_surface.msh",
                              "_wp.msh", "_surface.msh")
            for suffix in merge_suffixes:
                sibling = stem + suffix
                if os.path.isfile(sibling):
                    try:
                        gmsh.merge(sibling)
                        panel_log(f"open_gmsh.py: merged "
                                  f"{os.path.basename(sibling)}")
                    except Exception:
                        panel_log_exception(
                            f"open_gmsh.py: merge {os.path.basename(sibling)}")

        panel_log("open_gmsh.py: entering fltk.run")
        gmsh.fltk.run()
        gmsh.finalize()
        panel_log("open_gmsh.py: gmsh.finalize OK (window closed)")
        return 0
    except Exception:
        panel_log_exception("open_gmsh.py: gmsh raised")
        try:
            gmsh.finalize()
        except Exception:
            pass
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
