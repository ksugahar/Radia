"""Compatibility bridge; cubit-mesh-export owns the implementation."""
from pathlib import Path
import runpy

_bridge = runpy.run_path(str(Path(__file__).resolve().parent / '_cubit_gui_compat.py'))
_bridge['load_exporter_file']('toolbar_install.py', globals())
