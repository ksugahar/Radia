#!python
import os
import sys

# Cubit's official play route does not define __file__. WorkflowToolbar sets
# workingdir to this self-contained scripts directory before invoking us.
_here = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
if not os.path.isfile(os.path.join(_here, 'cubit_export_menu.py')):
    raise RuntimeError('Toolbar working directory is invalid; reimport the exporter toolbar package')
if _here not in sys.path:
    sys.path.insert(0, _here)
from cubit_export_menu import launch_export

launch_export("netgen_vol")
