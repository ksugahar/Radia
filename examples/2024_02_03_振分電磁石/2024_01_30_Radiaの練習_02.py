import os, sys
import radia as rad
from numpy import *
os.system('taskkill.exe /f /im i_view64.exe')

p = []
p.append(rad.ObjArcCur([0,0,0], [95, 105], [0,2*pi], 2000, 36, 1e6/2000/10, "man"))
g1 = rad.ObjCnt(p)

print(rad.UtiDmp(g1))
g = rad.ObjCnt([ g1])

print(g)
rad.ObjDrwOpenGL(g)


from my_module import *
exportGeometryToVTK(g, 'RADIA')

#res = rad.Solve(g, 0.01, 1000, 4)
#print(res)

x1 =  -50
y1 =    0
z1 =    0
x2 =   50
y2 =    0
z2 =    0 
np =  100

x = linspace(x1,x2,np);
y = linspace(y1,y2,np);
z = linspace(z1,z2,np);

import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('mathtext', **{'rm':'serif','it':'serif:italic','bf':'serif:bold','fontset':'cm'})

plt.figure(figsize=(3,3),dpi=400)
ax = plt.axes([0.20,0.13,0.77,0.82])

Bz = rad.FldLst(g,"Bz", [x1,y1,z1], [x2,y2,z2], int(np), "noarg")

#rad.FldCmpCrt(0.01, 0.001, 0.001, 0.01, 0.01, 0.001)
rad.FldCmpPrc('PrcB->1e-10, PrcCoord->1e-12')

r =  rad.FldPtcTrj(g, 0.002, [10.0, 0.0, 1.0, 0.0], [ 0.0, 10.0],  50)

r = array(r)
y  = r[:,0]
x  = r[:,1]
xp = r[:,2]
z  = r[:,3]
zp = r[:,4]
ax.plot(x,y, 'r.-', linewidth=0.5)

#ax.set_xlim(-1.2,1.2)
#ax.set_ylim(-1.5,0.5)
ax.set_xlabel('${\\it x}$ (mm)',fontname='times new roman')
ax.set_ylabel('${\\it y}$ (mm)',fontname='times new roman')
ax.axis('equal')

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

