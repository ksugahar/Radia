import os,  sys
import radia as rad

from numpy import *
from COIL import cCOIL

COIL = []
I = 1265.0

V = array([[1,0,0], [0,1,0], [0,0,1]]).T
W = 122
H = 122
L = 16.43186645*2
x0 = array([(48+170/1),-L/2,-20-W/2])
R = []
phi = 0
tilt = 0
COIL.append(cCOIL('GCE',I,x0,V,W,H,R,phi,tilt,L))

x0 = COIL[-1].x1
W  = COIL[-1].W
H  = COIL[-1].H
V  = COIL[-1].V1
R = 60+61
phi = 64.59228189
tilt = 90
L = 0.0
COIL.append(cCOIL('ARC',I,x0,V,W,H,R,phi,tilt,L))

x0 = COIL[-1].x1
W  = COIL[-1].W
H  = COIL[-1].H
V  = COIL[-1].V1
R = []
phi = 0
tilt = 90
L = 1018.51313197
COIL.append(cCOIL('GCE',I,x0,V,W,H,R,phi,tilt,L))

x0 = COIL[-1].x1
W  = COIL[-1].W
H  = COIL[-1].H
V  = COIL[-1].V1
R = 60+61
phi = 115.40771811
tilt = -90
L = 0.0
COIL.append(cCOIL('ARC',I,x0,V,W,H,R,phi,tilt,L))

x0 = COIL[-1].x1
W  = COIL[-1].W
H  = COIL[-1].H
V  = COIL[-1].V1
R = []
phi = 0
tilt = 90
L = 453.43186645*2
COIL.append(cCOIL('GCE',I,x0,V,W,H,R,phi,tilt,L))

x0 = COIL[-1].x1
W  = COIL[-1].W
H  = COIL[-1].H
V  = COIL[-1].V1
R = 60+61
phi = 115.40771811
tilt = -90
L = 0.0
COIL.append(cCOIL('ARC',I,x0,V,W,H,R,phi,tilt,L))

x0 = COIL[-1].x1
W  = COIL[-1].W
H  = COIL[-1].H
V  = COIL[-1].V1
R = []
phi = 0
tilt = 90
L = 1018.51313197
COIL.append(cCOIL('GCE',I,x0,V,W,H,R,phi,tilt,L))

x0 = COIL[-1].x1
W  = COIL[-1].W
H  = COIL[-1].H
V  = COIL[-1].V1
R = 60+61
phi = 64.59228189
tilt = -90
L = 0.0
COIL.append(cCOIL('ARC',I,x0,V,W,H,R,phi,tilt,L))

coils = []
for n in range(len(COIL)):
	if (COIL[n].type=='ARC'):
		phi1 = COIL[n].rot[0]/180*pi
		if (phi1 <= 0):
			phi1 = phi1 + 2*pi
		phi2 = (COIL[n].rot[0]+COIL[n].phi)/180*pi
		if (phi1 > phi2) or (phi2 <= 0):
			phi2 = phi2 + 2*pi
		coil = rad.ObjArcCur([0,0,0],[COIL[n].R-COIL[n].W/2,COIL[n].R+COIL[n].W/2], [phi1,phi2], COIL[n].H, 10, COIL[n].I/COIL[n].W/COIL[n].H, "auto")
		trf = rad.TrfRot([0,0,0], [0,0,1], COIL[n].rot[2]/180*pi)
		trf = rad.TrfCmbR(trf,rad.TrfRot([0,0,0], [1,0,0], COIL[n].rot[1]/180*pi))
	else:
		coil = rad.ObjRecCur([0,0,0],[COIL[n].W, COIL[n].L, COIL[n].H], [0,COIL[n].I/COIL[n].W/COIL[n].H,0])
		trf = rad.TrfRot([0,0,0], [0,0,1], COIL[n].rot[2]/180*pi)
		trf = rad.TrfCmbR(trf,rad.TrfRot([0,0,0], [1,0,0], COIL[n].rot[1]/180*pi))
		trf = rad.TrfCmbR(trf,rad.TrfRot([0,0,0], [0,0,1], COIL[n].rot[0]/180*pi))
	trf = rad.TrfCmbL(trf,rad.TrfTrsl([COIL[n].c[0],COIL[n].c[1],COIL[n].c[2]]))
	coil = rad.TrfOrnt(coil,trf)
	coils.append(coil)

coils = rad.ObjCnt(coils)
rad.ObjDrwAtr(coils, [1,0,0], 0.1)
rad.TrfZerPara(coils, [0,0,0], [1,0,0])

sys.path.append("C:/Program Files/Coreform Cubit 2023.11/bin")
import cubit
cubit.init(['-batch','-nojournal','-noecho','warning on','-nobanner'])

FileName = '2024_01_31_york'
#FileName = 'test'
with open(FileName + '.jou','r', encoding="utf-8") as fid:
	strLines = fid.readlines()
	for n in range(len(strLines)):
		cubit.silent_cmd(strLines[n])

#rad.FldLenTol(2e-3, 2e-3, 2e-3)

hex = []
for block_id in cubit.get_block_id_list():
	name = cubit.get_exodus_entity_name("block",block_id)
	print(f'block name = "{name}"')
	volume_list = cubit.get_block_volumes(block_id)
	for volume_id in volume_list:
		hex_list = cubit.get_volume_hexes(volume_id)
		if len(hex_list)>0:
			print(f'number of hexahedron = {len(hex_list)}')
			for hex_id in hex_list:
				connectivity_list = cubit.get_connectivity("hex", hex_id)
				coords = []
				for node_id in connectivity_list:
					coord = cubit.get_nodal_coordinates(node_id)
					coords.append(list(coord))

#				v21 = array(coords[1])-array(coords[0])
#				v31 = array(coords[2])-array(coords[0])
#				v41 = array(coords[3])-array(coords[0])
#				v21xv31 = cross(v21,v31)
#				print(v21xv31)
#				print(dot(v21xv31,v41))

#				face = [[1,4,3,2],[1,2,6,5],[2,3,7,6],[3,4,8,7],[4,1,5,8],[5,6,7,8]];
				face = [[1,2,3,4],[1,5,6,2],[2,6,7,3],[3,7,8,4],[1,4,8,5],[5,8,7,6]];

				M = [0, 0, 0.5]
				hex.append(rad.ObjPolyhdr(coords,face,M))


print(len(hex))

york = rad.ObjCnt(hex)
rad.ObjDrwAtr(york, [0,0,1], 0.1)

mat = rad.MatLin([1000,1000], [0, 0, 1])
rad.MatApl(york, mat)
rad.TrfZerPerp(york, [0,0,0], [0,1,0])
rad.TrfZerPara(york, [0,0,0], [1,0,0])

g = rad.ObjCnt([coils, york])

rad.ObjDrwOpenGL(g)
from my_module import *
exportGeometryToVTK(g, 'RADIA')

res = rad.Solve(g, 0.01, 1000, 4)
print(res)

B = rad.Fld(g, "BxByBz", [0,0,500])
print(B)

