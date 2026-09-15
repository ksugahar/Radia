"""Compatibility bridge; cubit-mesh-export owns the implementation."""
from pathlib import Path
import runpy

_bridge = runpy.run_path(str(Path(__file__).resolve().parent.parent / '_cubit_gui_compat.py'))
_bridge['load_exporter_file']('cubit_gui/radia_export_menu.py', globals())
