import os, sys
from numpy import *
import radia as rad
from my_module import *

sys.path.append("C:/Program Files/Coreform Cubit 2023.11/bin")
import cubit
cubit.init(["cubit", "-nobanner", "-nographics", "-noecho", "-nojournal", "-information", "off", "-warning", "off"])

FileName = 'Sphere'
with open(FileName + '.jou','r') as fid:
	strLines = fid.readlines()
	for n in range(len(strLines)):
		cubit.silent_cmd(strLines[n])

for block_id in cubit.get_block_id_list():
	name = cubit.get_exodus_entity_name("block",block_id)
	print(f'block name = "{name}"')
	surface_list = cubit.get_block_surfaces(block_id)
	for surface_id in surface_list:
		tri_list = cubit.get_surface_tris(surface_id)
		if len(tri_list)>0:
			print(f'number of tris = {len(tri_list)}')
			node_list = cubit.get_surface_nodes(surface_id)
			coords = []
			for node_id in node_list:
				coord = cubit.get_nodal_coordinates(node_id)
				coords.append(list(coord))
			faces = []
			for tri_id in tri_list:
				connectivity_list = cubit.get_connectivity("tri", tri_id)
				face = [connectivity_list[0], connectivity_list[2], connectivity_list[1]]
				faces.append(face)

M = [0,0,0]
g1 = rad.ObjPolyhdr(coords,faces,M)
mat1 = rad.MatLin([999,999], [0, 1, 1])
rad.MatApl(g1,mat1)

g2 = rad.ObjRecMag([0,0,100],[1,1,1],[0,0,1])
rad.MatApl(g2,mat1)

rad.ObjDrwOpenGL(g1)
g3 = rad.ObjBckg([0,1,0])

g = rad.ObjCnt([ g1, g3])
#g = rad.ObjCnt([g2])
#res = rad.RlxPre(g, 3, 1000, 3)

#print(rad.UtiDmp(2))
print(rad.UtiVer())

res = rad.Solve(g, 0.1, 100, 4)
print(res)

exportGeometryToVTK(g, 'RADIA')

x = linspace(-2.0, 2.0, 201)
y = linspace(-2.0, 2.0, 201)
[xx,yy] = meshgrid(x,y)
hz = 0
Bx = zeros((shape(xx)))
By = zeros((shape(xx)))
Bz = zeros((shape(xx)))
for ny in range(len(y)):
	for nx in range(len(x)):
		Bx[ny,nx] = rad.Fld(g,"Bx",[x[nx],y[ny],z])
		By[ny,nx] = rad.Fld(g,"By",[x[nx],y[ny],z])
		Bz[ny,nx] = rad.Fld(g,"Bz",[x[nx],y[ny],z])

from scipy.io import savemat
data = {'xx': xx, 'yy': yy, 'Bx':Bx, 'By':By, 'Bz':Bz}
FileName = f"{os.path.splitext(__file__)[0]}.mat"
savemat(FileName,data)

