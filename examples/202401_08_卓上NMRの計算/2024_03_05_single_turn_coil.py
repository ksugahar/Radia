import sys
import radia as rad
from numpy import *
from my_module import *

mu0 = 4*pi*1e-7

g1 = rad.ObjArcCur([0, 0, 0], [600, 700], [0*pi, 1*pi], 100, 12, 1.0, "auto")
g2 = rad.ObjArcCur([0, 0, 0], [600, 700], [1*pi, 2*pi], 100, 12, 1.0, "auto")
g = rad.ObjCnt([g1, g2])
rad.ObjDrwOpenGL(g)

Bz = rad.Fld(g,"Bz",[0, 0, 0])

print(f'Bz={Bz}T')

I = 1.0*100*100
a = 0.25
print(f'電磁気の教科書に載っている答えは、{mu0*I/(2*a)}T')

N = 101
x = linspace(-1000.0, 1000.0, N)
y = linspace(-1000.0, 1000.0, N)
z = linspace(-1000.0, 1000.0, N)

with open('radia_single_turn.csv', 'w') as fid:
	for nx in range(N):
		for ny in range(N):
			for nz in range(N):
				Bx = rad.Fld(g,"Bx",[x[nx],y[ny],z[nz]])
				By = rad.Fld(g,"By",[x[nx],y[ny],z[nz]])
				Bz = rad.Fld(g,"Bz",[x[nx],y[ny],z[nz]])
				if isinf(Bx):
					Bx = 0
				if isinf(By):
					By = 0
				if isinf(Bz):
					Bz = 0
				Bx = nan_to_num(Bx)
				By = nan_to_num(By)
				Bz = nan_to_num(Bz)
				fid.write(f'{x[nx]/1000}, {y[ny]/1000}, {z[nz]/1000}, {Bx/mu0}, {By/mu0}, {Bz/mu0}\n')
