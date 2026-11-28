"""Register the export menu inside Cubit's private Python runtime."""
import os
import sys

_directory = os.path.dirname(os.path.abspath(__file__))
if _directory not in sys.path:
    sys.path.insert(0, _directory)

import cubit_export_menu

if cubit_export_menu.install_menu():
    print("[radia-mcp] Export menu installed on the menu bar")
else:
    print("[radia-mcp] Cubit main window unavailable (batch mode)")
