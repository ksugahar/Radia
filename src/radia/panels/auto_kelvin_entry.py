"""Cubit entry point for Auto-Kelvin addition.

Invoked by the Radia-NGSolve C++ launcher (RadiaComp.cpp::
launch_radia_ngsolve) via::

    play "<panels_dir>/auto_kelvin_entry.py"

immediately before the ``radia_export netgen`` call.  The launcher's
"Add Kelvin open boundary (auto)" checkbox gates whether this file is
played.

All the work lives in ``add_kelvin.auto_add_kelvin_from_current_model``;
this entry script only handles sys.path wiring so the import works no
matter which working directory Cubit is in.
"""
import os
import sys

# This script lives in radia/panels/ so its own directory is exactly
# where add_kelvin.py also lives.  Insert at the front so a local
# add_kelvin.py shadows any stale site-packages copy during development.
_here = os.path.dirname(os.path.abspath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)

from add_kelvin import auto_add_kelvin_from_current_model

try:
    _info = auto_add_kelvin_from_current_model()
    if _info is not None:
        ox, oy, oz = _info["center"]
        print(f"[auto_kelvin_entry] Kelvin added at "
              f"({ox:.3f}, {oy:.3f}, {oz:.3f}), R={_info['R']:.4f}, "
              f"symmetry={_info['symmetry']}")
except Exception as _e:
    print(f"[auto_kelvin_entry] ERROR: {_e}")
    # Do not re-raise — the launcher should proceed even if Kelvin
    # auto-addition fails (user will get Dirichlet truncation).
