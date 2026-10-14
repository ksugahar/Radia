#!python
import os
import sys

# Cubit's play route can retain a previous __file__. WorkflowToolbar sets
# workingdir to this self-contained scripts directory before invoking us.
_here = os.getcwd()
if not os.path.isfile(os.path.join(_here, 'cubit_export_menu.py')) and '__file__' in globals():
    _here = os.path.dirname(os.path.abspath(__file__))
if not os.path.isfile(os.path.join(_here, 'cubit_export_menu.py')):
    raise RuntimeError('Toolbar working directory is invalid; run cubit-plugin-install to repair the toolbar')
if _here not in sys.path:
    sys.path.insert(0, _here)
from cubit_export_menu import launch_export

launch_export("gmsh")
