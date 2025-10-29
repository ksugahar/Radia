import os, sys
import scipy.io
from numpy import *

import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('mathtext', **{'rm':'serif','it':'serif:italic','bf':'serif:bold','fontset':'cm'})

plt.figure(figsize=(3,3),dpi=400)
ax = plt.axes([0.20,0.13,0.77,0.82])

data = scipy.io.loadmat("2024_01_28_Sphere.mat")
xx = data['xx']
By = data['By']-1

print(By)

ax.plot(xx[100,:], By[100,:], 'k-', linewidth=0.5)

mur = 1000
x = xx[100,:]
Hy = zeros_like(x)
a = 0.5
Hy[abs(x) < a] = -1 + 3 / (mur + 2)
t = -(mur - 1) / (mur + 2) * (a / abs(x)) ** 3
Hy[abs(x) >= a] = t[abs(x) >= a]
ax.plot(x, Hy, 'r:', linewidth=0.5)

ax.set_xlim(-1.2,1.2)
ax.set_ylim(-1.5,0.5)
ax.set_xlabel('${\\it x}$ (m)',fontname='times new roman')
ax.set_ylabel('${\\it B_y}$ (A/m)',fontname='times new roman')

ax.legend(['Radia', 'Analytic'], frameon=False, bbox_to_anchor=(0.35,0.85),loc='center',framealpha=0.0,fancybox='none',edgecolor='none',prop={"family":'Times New Roman',"size":8})
ax.xaxis.grid(True, which='major', linestyle=':', linewidth=0.2)
ax.xaxis.grid(True, which='minor', linestyle=':', linewidth=0.1)
ax.yaxis.grid(True, which='major', linestyle=':', linewidth=0.2)
ax.yaxis.grid(True, which='minor', linestyle=':', linewidth=0.1)
ax.tick_params(left=True, right=True, top=True, bottom=True, which="major", direction='in', length=2 ,width=0.4, pad=3)

plt.setp(ax.get_xticklabels(),fontname='times new roman')
plt.setp(ax.get_yticklabels(),fontname='times new roman')

FileName = f"{os.path.splitext(__file__)[0]}.png"
plt.savefig(FileName)
os.startfile(FileName)

