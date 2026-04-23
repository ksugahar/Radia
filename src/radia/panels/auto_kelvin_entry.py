"""Cubit entry point for Auto-Kelvin addition (argument-driven).

Invoked by the Radia-NGSolve C++ launcher via::

    play "<panels_dir>/auto_kelvin_entry.py"

just before the `radia_export netgen` call.  The launcher writes a JSON
config first and exports its path via the ``RADIA_LAUNCHER_CONFIG``
environment variable; this script reads the JSON and dispatches.

Config schema (all keys optional; defaults in parens)::

    {
      "add_kelvin":           bool         (default True)
      "kelvin_air_block":     str          (default "air")
      "kelvin_block_name":    str          (default "kelvin")
      "kelvin_mesh_size":     float|null   (default null = auto from air)
    }

``kelvin_mesh_size`` is the tet edge length [m] imposed on the Kelvin
exterior sphere.  ``null`` (default) lets add_kelvin_cubit inherit the
size from the air outer surface via copy-mesh — usually fine.  For
large models you may want a COARSER Kelvin mesh (e.g. 2-3x the air
surface size); specify an explicit float to override.

If ``RADIA_LAUNCHER_CONFIG`` is not set or the file is missing, all
defaults apply — i.e. the pre-argument-driven behaviour is preserved,
so existing samples / tests keep working.

Argument-driven rationale (2026-04-23, CLAUDE.md § Panel Design
Workflow Policy):  the Cubit panel's only job becomes collecting
dialog state and writing this JSON.  Every knob that matters for
computation is in the config file, so `pytest` can exercise every
combination without driving the GUI.  Per the policy, a panel is
Stage-3-ready only after Stage-2 (this CLI / config layer) is
golden-locked.
"""
import json
import os
import sys

# Locate add_kelvin: this file and add_kelvin.py live in the same
# directory (radia/panels/), so its own dir is on sys.path.
_here = os.path.dirname(os.path.abspath(__file__))
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
            mesh_size=cfg["kelvin_mesh_size"])
        if _info is not None:
            ox, oy, oz = _info["center"]
            print(f"[auto_kelvin_entry] Kelvin added at "
                  f"({ox:.3f}, {oy:.3f}, {oz:.3f}), "
                  f"R={_info['R']:.4f}, symmetry={_info['symmetry']}")
        else:
            print("[auto_kelvin_entry] auto_add_kelvin returned None "
                  "(already present or no air block) -- noop")
    except Exception as _e:
        # Never re-raise: the launcher should proceed even if Kelvin
        # fails (user gets Dirichlet truncation + a console warning).
        print(f"[auto_kelvin_entry] ERROR: {_e}")
