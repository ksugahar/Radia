import os, sys
import radia as rad
from numpy import *
os.system('taskkill.exe /f /im i_view64.exe')

p = []
p.append(rad.ObjArcCur([0,0,0], [95, 105], [0,2*pi], 10, 36, 1e6/10/10, "man"))
g1 = rad.ObjCnt(p)

##rad.TrfZerPara(g, [0,0,0], [1,0,0])
#rad.TrfZerPerp(g, [0,0,0], [1,0,0])
#rad.TrfZerPerp(g, [0,0,0], [0,1,0])
#rad.ObjDrwOpenGL(g)

M = [0,0,0]
g2 = rad.ObjRecMag( [0,0,0], [10,10,5], [0,0,1])
mat2 = rad.MatLin([1,1], [0, 0, 1])
rad.MatApl(g2, mat2)

g3 = []
z = arange(5, 20, 5)
for n in range(len(z)):
	g3.append(rad.ObjRecMag([0,0,z[n]],[10,10,5],[0,0,1]))

print(g3)

g3 = rad.ObjCnt(g3)

rad.TrfZerPerp(g3, [0,0,0], [0,0,1])

rad.ObjDrwAtr(g3, [0,1,1], 0.1)
mat3 = rad.MatLin([1,1], [0, 0, 1])
rad.MatApl(g3, mat3)

print(rad.UtiDmp(g2))
g = rad.ObjCnt([ g1, g2, g3])

print(g)
rad.ObjDrwOpenGL(g)

from my_module import *
exportGeometryToVTK(g, 'RADIA')

res = rad.Solve(g, 0.01, 1000, 4)
print(res)

x1 =    0
y1 =    0
z1 = -100
x2 =    0
y2 =    0
z2 =  100
np =  100
z = linspace(z1,z2,np);

import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('mathtext', **{'rm':'serif','it':'serif:italic','bf':'serif:bold','fontset':'cm'})

plt.figure(figsize=(3,3),dpi=400)
ax = plt.axes([0.20,0.13,0.77,0.82])


Bz = rad.FldLst(g,"Bz", [x1,y1,z1], [x2,y2,z2], int(np), "noarg")

ax.plot(z, Bz, 'r-', linewidth=0.5)

mu0 = 4*pi*1e-7*1e6
a = 0.1
I = 1
z = z/1000
Bz = mu0*I*a**2/2/(a**2+z**2)**(3/2)
ax.plot(z*1000, Bz, 'k--', linewidth=0.5)

#ax.set_xlim(-1.2,1.2)
#ax.set_ylim(-1.5,0.5)
ax.set_xlabel('${\\it x}$ (m)',fontname='times new roman')
ax.set_ylabel('${\\it B_y}$ (A/m)',fontname='times new roman')

ax.legend(['Radia', 'Analytic'], frameon=False, bbox_to_anchor=(0.20,0.85),loc='center',framealpha=0.0,fancybox='none',edgecolor='none',prop={"family":'Times New Roman',"size":6})
ax.xaxis.grid(True, which='major', linestyle=':', linewidth=0.2)
ax.xaxis.grid(True, which='minor', linestyle=':', linewidth=0.1)
ax.yaxis.grid(True, which='major', linestyle=':', linewidth=0.2)
ax.yaxis.grid(True, which='minor', linestyle=':', linewidth=0.1)

ax.tick_params(left=True, right=True, top=True, bottom=True, which="major", direction='in', length=2 ,width=0.4, pad=3)

plt.setp(ax.get_xticklabels(),fontname='times new roman')
plt.setp(ax.get_yticklabels(),fontname='times new roman')

FileName = os.path.basename(__file__).replace(".py",".png")
plt.savefig(FileName)
os.startfile(FileName)

