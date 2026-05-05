"""Cubit entry point for Auto-Kelvin addition (argument-driven).

Invoked by ``radia_export netgen ... add_kelvin ...`` (C++) via::

    play "<plugin_dir>/cubit_helpers/auto_kelvin_entry.py"

just before the actual Netgen .vol export.  The C++ command writes a
JSON config file first and exports its path via the
``RADIA_LAUNCHER_CONFIG`` environment variable; this script reads the
JSON and dispatches to ``add_kelvin.auto_add_kelvin_from_current_model``.

Config schema (all keys optional; defaults in parens)::

    {
      "add_kelvin":           bool         (default True)
      "kelvin_air_block":     str          (default "air")
      "kelvin_block_name":    str          (default "kelvin")
      "kelvin_mesh_size":     float|null   (default null = auto from air)
      "kelvin_reduction":     dict|null    (default null = mesh-seam mode)
    }

``kelvin_mesh_size`` is the tet edge length [m] imposed on the Kelvin
exterior sphere.  ``null`` (default) lets add_kelvin_cubit inherit the
size from the air outer surface via copy-mesh.

``kelvin_reduction`` is a JSON object keyed by axis ("x"|"y"|"z") with
value in {"bn=0", "ht=0"}.  When present, it triggers domain-reduction
mode in add_kelvin_cubit.  ``null`` (default) runs the existing
symmetry auto-detection (mesh-seam mode).  Example::

    "kelvin_reduction": {"x": "ht=0", "z": "bn=0"}   // 1/4 xz model

If ``RADIA_LAUNCHER_CONFIG`` is not set or the file is missing, all
defaults apply.
"""
import json
import os
import sys

# Locate add_kelvin.py.  This file and add_kelvin.py live in the same
# directory (cubit_helpers/).  Cubit's `play` exec's .py via
# `exec(compile(source, "<string>", "exec"))` and does NOT bind
# __file__ in the exec globals, so `os.path.abspath(__file__)` raises
# NameError when invoked through `play`.  Three fallbacks in order:
#   1. __file__              -- available when imported normally
#   2. CUBIT_HELPERS_DIR env -- set by the C++ caller alongside
#                                RADIA_LAUNCHER_CONFIG
#   3. sys.path search       -- find add_kelvin.py on any entry
_here = None
if "__file__" in globals():
    _here = os.path.dirname(os.path.abspath(__file__))
if not _here:
    _here = os.environ.get("CUBIT_HELPERS_DIR")
if not _here or not os.path.isfile(os.path.join(_here, "add_kelvin.py")):
    # Back-compat with the pre-2026-05-05 RADIA_PANELS_DIR name; old
    # code paths and tests may still set it.
    legacy = os.environ.get("RADIA_PANELS_DIR")
    if legacy and os.path.isfile(os.path.join(legacy, "add_kelvin.py")):
        _here = legacy
if not _here or not os.path.isfile(os.path.join(_here, "add_kelvin.py")):
    for _p in sys.path:
        if _p and os.path.isfile(os.path.join(_p, "add_kelvin.py")):
            _here = _p
            break
if not _here or not os.path.isfile(os.path.join(_here, "add_kelvin.py")):
    raise RuntimeError(
        "auto_kelvin_entry.py: cannot locate cubit_helpers dir "
        "(no __file__, no CUBIT_HELPERS_DIR env, add_kelvin.py not in "
        "sys.path).  The C++ caller must set CUBIT_HELPERS_DIR or "
        "prepend the cubit_helpers directory to sys.path before `play`.")
if _here not in sys.path:
    sys.path.insert(0, _here)

from add_kelvin import auto_add_kelvin_from_current_model


def _load_config():
    """Return the dict config.  Defaults apply to missing keys."""
    defaults = {
        "add_kelvin":        True,
        "kelvin_air_block":  "air",
        "kelvin_block_name": "kelvin",
        "kelvin_mesh_size":  None,
        "kelvin_reduction":  None,
    }
    cfg_path = os.environ.get("RADIA_LAUNCHER_CONFIG", "")
    if not cfg_path:
        print("[auto_kelvin_entry] no RADIA_LAUNCHER_CONFIG set -- "
              "using defaults (add_kelvin=True)")
        return defaults
    if not os.path.isfile(cfg_path):
        print(f"[auto_kelvin_entry] WARN: config file not found: "
              f"{cfg_path} -- using defaults")
        return defaults
    try:
        with open(cfg_path, "r", encoding="utf-8") as f:
            loaded = json.load(f)
        out = dict(defaults)
        out.update({k: v for k, v in loaded.items() if k in defaults})
        print(f"[auto_kelvin_entry] loaded config from {cfg_path}: {out}")
        return out
    except Exception as e:
        print(f"[auto_kelvin_entry] WARN: failed to parse {cfg_path}: "
              f"{e} -- using defaults")
        return defaults


cfg = _load_config()

if not cfg["add_kelvin"]:
    print("[auto_kelvin_entry] add_kelvin=False -- skipping "
          "(explicit opt-out via config).")
else:
    try:
        _info = auto_add_kelvin_from_current_model(
            air_block=cfg["kelvin_air_block"],
            kelvin_block=cfg["kelvin_block_name"],
            mesh_size=cfg["kelvin_mesh_size"],
            reduction=cfg["kelvin_reduction"])
        if _info is not None:
            ox, oy, oz = _info["center"]
            print(f"[auto_kelvin_entry] Kelvin added at "
                  f"({ox:.3f}, {oy:.3f}, {oz:.3f}), "
                  f"R={_info['R']:.4f}, symmetry={_info['symmetry']}")
        else:
            print("[auto_kelvin_entry] auto_add_kelvin returned None "
                  "(already present or no air block) -- noop")
    except Exception as _e:
        # Never re-raise: the export should proceed even if Kelvin
        # fails (user gets Dirichlet truncation + a console warning).
        print(f"[auto_kelvin_entry] ERROR: {_e}")
