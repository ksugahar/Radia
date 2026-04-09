"""Clean GMSH script: should produce zero lint findings.

Demonstrates correct patterns for GMSH post-processing workflow.
"""

import os
from pathlib import Path

# Correct: no gmsh import (use Cubit or OCC for meshing)
from ngsolve import Mesh

# Correct: load .vol file (not .msh via ReadGmsh)
mesh = Mesh("model.vol")

# Correct: no meshio import
# Correct: no removed builder class reference

# Correct: post-processing only, no mesh generation
print(f"Number of elements: {mesh.ne}")
