"""Register the export menu inside Cubit's private Python runtime."""
import os
import sys

_directory = os.path.dirname(os.path.abspath(__file__))
if _directory not in sys.path:
    sys.path.insert(0, _directory)

import radia_export_menu

if radia_export_menu.install_menu():
    print("[cubit-mesh-export] Export menu installed through Claro")
else:
    print("[cubit-mesh-export] Claro GUI unavailable (batch mode)")
