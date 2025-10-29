import os, sys
import scipy.io
from numpy import *

import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('mathtext', **{'rm':'serif','it':'serif:italic','bf':'serif:bold','fontset':'cm'})

plt.figure(figsize=(3,3),dpi=400)
ax = plt.axes([0.17,0.13,0.77,0.82])

data = scipy.io.loadmat("2024_01_28_Sphere.mat")
xx = data['xx']
yy = data['yy']
Bx = data['Bx']
By = data['By']-1
B = sqrt(Bx**2+By**2)

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
matplotlib.rc('mathtext', **{'rm':'serif','it':'serif:italic','bf':'serif:bold','fontset':'cm'})

fig = plt.figure(figsize=(3,3),dpi=200)
ax1 = plt.axes([0.2, 0.2, 0.77, 0.77])

lv = arange(0.0, 2.0, 0.01)
ax1.contourf(xx, yy, B, levels=lv,cmap='jet',extend='both')

ax1.set_xlim(-1.5,1.5)
ax1.set_ylim(-1.5,1.5)
plt.setp(ax1.get_xticklabels(),fontname='Times New Roman',fontsize=10)
ax1.set_xlabel('${\it x}$ (m)',fontname='Times New Roman',fontsize=10)
plt.setp(ax1.get_yticklabels(),fontname='Times New Roman',fontsize=10)
ax1.set_ylabel('${\it y}$ (m)',fontname='Times New Roman',fontsize=10)
ax1.minorticks_on()
ax1.tick_params(which = 'major', direction="in",top=True)
ax1.tick_params(which = 'minor', direction="in",top=True)
ax1.grid(axis='x', which = 'major', c='gainsboro', linestyle=':', linewidth=0.2)
ax1.grid(axis='x', which = 'minor', c='gainsboro', linestyle='--' , linewidth=0.1)
ax1.grid(axis='y', which = 'major', c='gainsboro', linestyle=':', linewidth=0.2)
ax1.grid(axis='y', which = 'minor', c='gainsboro', linestyle='--' , linewidth=0.1)

FileName = f"{os.path.splitext(__file__)[0]}.png"
plt.savefig(FileName)
os.startfile(FileName)



