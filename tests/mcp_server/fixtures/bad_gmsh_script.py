"""Bad GMSH script: intentionally triggers gmsh lint rules.

Expected findings:
- gmsh-mesh-generation (CRITICAL): GMSH API used for mesh generation
- gmsh-builder-removed (CRITICAL): GmshBuilder reference
- readgmsh-deprecated (HIGH): ReadGmsh usage
- pip-gmsh-import (HIGH): import gmsh in regular script
- meshio-removed (HIGH): import meshio
- numsubedges-missing (MODERATE): high-order elements without NumSubEdges
"""

# pip-gmsh-import: import gmsh in regular script
import gmsh

# gmsh-mesh-generation: GMSH Python API for mesh generation
gmsh.initialize()
gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.finalize()

# gmsh-builder-removed: reference to removed GmshBuilder class
from radia import GmshBuilder

# readgmsh-deprecated: deprecated ReadGmsh function
from ngsolve import Mesh
mesh = Mesh("model.msh")  # not the trigger
from netgen.read_gmsh import ReadGmsh
vol_mesh = ReadGmsh("model.msh")

# meshio-removed: import meshio (removed from project)
import meshio

# numsubedges-missing: high-order + GmshPostExport without NumSubEdges
mesh.Curve(3)
from radia.gmsh_post import GmshPostExport
GmshPostExport(mesh, "output.msh")
