import os, sys
import radia as rad
from numpy import *
from my_module import *

meshsize = 10

sys.path.append("C:/Program Files/Coreform Cubit 2023.11/bin")
import cubit
cubit.init(["cubit","-nojournal"])
cubit.cmd('reset')
cubit.cmd('create Cylinder height 20 radius 60')
cubit.cmd('move Volume 1 z 10 include_merged ')
cubit.cmd('webcut volume all with plane zplane offset 10 ')
cubit.cmd('mearge all')
cubit.cmd('imprint all')
cubit.cmd('volume all size {meshsize}')
cubit.cmd('mesh volume all')
cubit.cmd('block 1 add vol 1')
cubit.cmd('block 1 name "magnet"')
cubit.cmd('block 2 add vol 2')
cubit.cmd('block 2 name "iron"')


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
				M = [0,0,0.5]
				g_ = rad.ObjPolyhdr(coords,face,M)
				if block_id == 1:
					if "g1" in locals():
						g1 = rad.ObjCnt([g1, g_])
					else:
						g1 = g_
				if block_id == 2:
					if "g2" in locals():
						g2 = rad.ObjCnt([g2, g_])
					else:
						g2 = g_

g = rad.ObjCnt([g1, g2])

exportGeometryToVTK(g, 'RADIA')

print("end")

	#磁化された多面体を作成する
	#points： {{x1,y1,z1}, {x2,y2,z2},...}: 多面体の頂点のデカルト座標 (デフォルトでは mm) を指定する 3 つの実数のリストのリスト
	#face：{{f1n1,f1n2,...}, {f2n1,f2n2,...},...}: 多面体の面を形成する前のリストの頂点のインデックスを指定する整数のリストのリスト
	#M：{mx,my,mz}: 多面体の内部の磁化ベクトルを指定する 3 つの実数のリスト (デフォルトではテスラ)
#	g = rad.ObjPolyhdr(points.tolist(),face,M)


