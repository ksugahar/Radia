import os, sys
import numpy as np
sys.path.append("C:/Program Files/Coreform Cubit 2025.3/bin")

import cubit
cubit.init(['cubit','-nojournal','-batch'])

with open('York.jou','r', encoding='utf8') as fid:
	strLines = fid.readlines()
	for n in range(len(strLines)):
		cubit.cmd(strLines[n])

import cubit_mesh_export
FileName = 'York'
cubit_mesh_export.export_Gmsh_ver2(cubit, FileName + '.msh')
cubit_mesh_export.export_Nastran(cubit, FileName + '.bdf', DIM='3D', PYRAM=False)
cubit_mesh_export.export_vtk(cubit, FileName + '.vtk', ORDER="1st")

