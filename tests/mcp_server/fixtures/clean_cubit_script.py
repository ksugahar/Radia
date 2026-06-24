"""Clean Cubit script: should produce zero lint findings.

Demonstrates correct patterns for all common Cubit operations.
"""

import os
import sys
import cubit

# Correct: cubit.init() called
cubit.init(["cubit", "-nojournal", "-nographics"])

# Correct: relative path usage
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src"))

# Correct: mesh command present
cubit.cmd("mesh volume all")

# Correct: block registration with names
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "iron"')
cubit.cmd("block 2 add tri all")
cubit.cmd('block 2 name "boundary"')

# Correct: element type set AFTER adding elements
cubit.cmd("block 1 element type tetra4")

# Correct: get_expanded_connectivity for 2nd-order (if needed)
# nodes = cubit.get_expanded_connectivity("tet", 1)

# Correct: export with matching extension and proper setup
import cubit_mesh_export as cme
cme.export_Gmsh_ver4(cubit, "output.msh")
cubit.cmd('export netgen "output.vol" order 2 overwrite')
