"""Lint fixture: intentionally violates every GMSH policy rule.

NEVER executed -- only read as text by the gmsh lint self-test and the
pytest fixture lock. Each line below exists to trip one specific rule.
"""
import gmsh  # pip-gmsh-import
import meshio  # meshio-removed
from radia.gmsh_builder import GmshBuilder  # gmsh-builder-removed

gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)  # gmsh-mesh-generation
gmsh.model.mesh.generate(3)  # gmsh-mesh-generation
gmsh.option.setNumber("Mesh.Volumes", 1)  # invalid-gmsh-option
gmsh.option.setNumber("General.GraphicsSizeX", 800)  # invalid-gmsh-option

mesh = ReadGmsh("legacy.msh")  # readgmsh-deprecated
