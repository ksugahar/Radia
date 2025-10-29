import os,  sys
import radia as rad

from numpy import *
from COIL import cCOIL

for I_ in [130,245,460,620,700,805,900,1000,1150,1200,1265]:
	COIL = []
	I = 64*I_

	V = array([[1,0,0], [0,1,0], [0,0,1]]).T
	W = 122
	H = 122
	L = 16.43186645*2
	x0 = array([(48+170)/2,-L/2,-20-W/2])
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
	cubit.init(['-batch','-nojournal','-noecho','-warning','on','-nobanner'])

#	FileName = '2024_01_31_york'
#	FileName = '2024_01_31_york_tetmesh'
	FileName = '2024_02_02_york_tetmesh'
	with open(FileName + '.jou','r', encoding="utf-8") as fid:
		strLines = fid.readlines()
		for n in range(len(strLines)):
			cubit.silent_cmd(strLines[n])

#rad.FldLenTol(2e-3, 2e-3, 2e-3)

	york = []
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
					face = [[1,2,3,4],[1,5,6,2],[2,6,7,3],[3,7,8,4],[1,4,8,5],[5,8,7,6]];
					M = [0, 0, 0.5]
					hex.append(rad.ObjPolyhdr(coords,face,M))
			tet_list = cubit.get_volume_tets(volume_id)
			if len(tet_list)>0:
				print(f'number of tetrahedron = {len(tet_list)}')
				for tet_id in tet_list:
					connectivity_list = cubit.get_connectivity("tet", tet_id)
					coords = []
					for node_id in connectivity_list:
						coord = cubit.get_nodal_coordinates(node_id)
						coords.append(list(coord))
					face = [[1,3,2],[1,2,4],[2,3,4],[3,1,4]];
					M = [0, 0, 0.5]
					york.append(rad.ObjPolyhdr(coords,face,M))

#					v21 = array(coords[1])-array(coords[0])
#					v31 = array(coords[2])-array(coords[0])
#					v41 = array(coords[3])-array(coords[0])
#					v21xv31 = cross(v21,v31)
#					print(v21xv31)
#					print(dot(v21xv31,v41))
#					face = [[1,4,3,2],[1,2,6,5],[2,3,7,6],[3,4,8,7],[4,1,5,8],[5,6,7,8]];

	york = rad.ObjCnt(york)
	rad.ObjDrwAtr(york, [0,0,1], 0.1)

	if True:
		mu0 = 4.*3.141592653589793e-07
		H = [10, 16, 20, 25, 30, 34, 40, 47, 52, 59, 62, 66, 72, 76, 80, 86, 90, 100, 110, 120, 130, 140, 150, 180, 220, 250, 300, 400, 520, 600, 700, 900, 1300, 1700, 2500, 3000, 4200, 5000, 6600, 9000, 10000, 15000, 22000, 30000, 42000, 66000, 100000, 640599, 6207043]
		B = [0.0156, 0.02925, 0.039, 0.0585, 0.078, 0.0975, 0.1365, 0.195, 0.2535, 0.351, 0.4095, 0.4875, 0.60938, 0.6825, 0.741, 0.819, 0.858, 0.9555, 1.0335, 1.092, 1.131001, 1.170003, 1.189506, 1.248049, 1.306743, 1.326433, 1.366402, 1.408496, 1.456667, 1.485215, 1.517785, 1.554001, 1.628569, 1.662114, 1.71874, 1.74328, 1.788677, 1.810448, 1.853225, 1.911445, 1.927287, 1.980139, 2.043564, 2.085955, 2.139015, 2.192476, 2.246281, 3.01226, 10.57212]
		M = list(array(B) - mu0*array(H))
		mat = rad.MatSatIsoTab([[mu0*H[i], M[i]] for i in range(len(H))])
	else:
		mat = rad.MatLin([999,999], [0, 0, 1])

	rad.MatApl(york, mat)
	rad.TrfZerPerp(york, [0,0,0], [0,1,0])
	rad.TrfZerPara(york, [0,0,0], [1,0,0])

	g = rad.ObjCnt([coils, york])

	rad.ObjDrwOpenGL(g)
	from my_module import *
	exportGeometryToVTK(g, 'RADIA')

	import time
	start_time = time.time()
	res = rad.Solve(g, 0.01, 5000, 4)
	print(f'number of iterations = {res[3]:.0f}')
	print(f'elapsed_time = {time.time()-start_time:.2f} sec')

	B = rad.Fld(g, "BxByBz", [0,0,600])
	print(B)

	z = linspace(-350,1370,173)
	Bx = rad.FldLst(g,"Bx", [0,0,z[0]], [0,0,z[-1]], int(len(z)), "noarg")

	from scipy.io import savemat
	data = {"z": z, "Bx": Bx}
	savemat(f'O:/Radia_tetra_I={I_}.mat',data)

sys.exit()

Nx = 31
Ny = 31
Nz = 31
x = linspace(-500,500,Nx)
y = linspace(-500,500,Ny)
z = linspace(-500,500,Nz)
dx = (x[1]-x[0])
dy = (y[1]-y[0])
dz = (z[1]-z[0])
with open('MagneticFluxDensity.vtk','w') as fid:
	fid.write(f'# vtk DataFile Version 2.0\n')
	fid.write(f'Generated by ksugahar\n')
	fid.write(f'ASCII\n')
	fid.write(f'DATASET STRUCTURED_POINTS\n')
	fid.write(f'DIMENSIONS {Nx} {Ny} {Nz}\n')
	fid.write(f'ORIGIN {x[0]} {y[0]} {z[0]}\n')
	fid.write(f'SPACING {dx} {dy} {dz}\n')
	fid.write(f'POINT_DATA {Nx*Ny*Nz}\n')
	fid.write(f'SCALARS B_magnitude float\n')
	fid.write(f'LOOKUP_TABLE default\n')
	for nz in range(Nz):
		for ny in range(Ny):
			for nx in range(Nx):
				bx = rad.Fld(g,"Bx",[x[nx],y[ny],z[nz]])
				by = rad.Fld(g,"By",[x[nx],y[ny],z[nz]])
				bz = rad.Fld(g,"Bz",[x[nx],y[ny],z[nz]])
				b = sqrt(bx**2+by**2+bz**2)
				fid.write(f'{b}\n');
	fid.write(f'VECTORS B_vectors float\n')
	for nz in range(Nz):
		for ny in range(Ny):
			for nx in range(Nx):
				bx = rad.Fld(g,"Bx",[x[nx],y[ny],z[nz]])
				by = rad.Fld(g,"By",[x[nx],y[ny],z[nz]])
				bz = rad.Fld(g,"Bz",[x[nx],y[ny],z[nz]])
				fid.write(f'{bx} {by} {bz}\n');
	fid.close()

print("end")
