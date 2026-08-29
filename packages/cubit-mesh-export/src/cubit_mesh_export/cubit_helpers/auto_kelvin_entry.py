"""Cubit entry point for Auto-Kelvin addition (argument-driven).

Invoked by ``export netgen ... add_kelvin ...`` (C++) via::

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
import importlib.util
import os
import sys


def _process_environment(name):
    """Read the live process environment, bypassing embedded-Python caches."""
    if os.name != "nt":
        return os.environ.get(name)
    import ctypes
    kernel32 = ctypes.windll.kernel32
    size = kernel32.GetEnvironmentVariableW(name, None, 0)
    if not size:
        return None
    buffer = ctypes.create_unicode_buffer(size)
    if not kernel32.GetEnvironmentVariableW(name, buffer, size):
        return None
    return buffer.value

# Locate add_kelvin.py.  This file and add_kelvin.py live in the same
# directory (cubit_helpers/).  Cubit's `play` exec's .py via
# `exec(compile(source, "<string>", "exec"))` and does NOT bind
# __file__ in the exec globals, so `os.path.abspath(__file__)` raises
# NameError when invoked through `play`.  Two authoritative locations:
#   1. __file__              -- available when imported normally
#   2. CUBIT_HELPERS_DIR env -- set by the C++ caller alongside
#                                RADIA_LAUNCHER_CONFIG
# ``play`` reuses Cubit's global Python dictionary, including a stale
# ``__file__`` left by a previously executed Radia script.  The C++ command's
# per-call environment variable is therefore authoritative.
_here = _process_environment("CUBIT_HELPERS_DIR")
if not _here:
    _exe_dir = os.path.dirname(os.path.abspath(sys.executable))
    for _candidate in (os.path.join(_exe_dir, "plugins", "cubit_helpers"), os.path.join(_exe_dir, "bin", "plugins", "cubit_helpers")):
        if os.path.isfile(os.path.join(_candidate, "add_kelvin.py")):
            _here = _candidate
            break
if not _here and "__file__" in globals():
    _here = os.path.dirname(os.path.abspath(__file__))
_add_kelvin_path = os.path.join(_here or "", "add_kelvin.py")
if not _here or not os.path.isfile(_add_kelvin_path):
    raise RuntimeError("auto_kelvin_entry.py: cannot locate cubit_helpers dir (no __file__ and no valid CUBIT_HELPERS_DIR env). The C++ caller must set CUBIT_HELPERS_DIR before `play`.")
print(f"[auto_kelvin_entry] helper source: {_add_kelvin_path}")

# Cubit's long-lived Python interpreter may already contain an unrelated
# ``sys.modules['add_kelvin']`` from an earlier Radia installation.  A normal
# import would silently reuse it even after CUBIT_HELPERS_DIR is updated.  Load
# the deployed sibling by exact path and a private name so the export command
# always uses the helper shipped beside this entry point.
_spec = importlib.util.spec_from_file_location("_cubit_mesh_export_auto_kelvin_impl", _add_kelvin_path)
if _spec is None or _spec.loader is None:
    raise RuntimeError(f"auto_kelvin_entry.py: cannot load helper: {_add_kelvin_path}")
_add_kelvin_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_add_kelvin_module)
auto_add_kelvin_from_current_model = _add_kelvin_module.auto_add_kelvin_from_current_model


def _load_config():
    """Return the dict config.  Defaults apply to missing keys."""
    defaults = {
        "add_kelvin":        True,
        "kelvin_air_block":  "air",
        "kelvin_block_name": "kelvin",
        "kelvin_mesh_size":  None,
        "kelvin_reduction":  None,
    }
    cfg_path = _process_environment("RADIA_LAUNCHER_CONFIG") or ""
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
    _info = auto_add_kelvin_from_current_model(air_block=cfg["kelvin_air_block"], kelvin_block=cfg["kelvin_block_name"], mesh_size=cfg["kelvin_mesh_size"], reduction=cfg["kelvin_reduction"])
    if _info is not None:
        ox, oy, oz = _info["center"]
        print(f"[auto_kelvin_entry] Kelvin added at ({ox:.3f}, {oy:.3f}, {oz:.3f}), R={_info['R']:.4f}, symmetry={_info['symmetry']}")
    else:
        print("[auto_kelvin_entry] auto_add_kelvin returned None "
              "(already present or no eligible air block) -- noop")
