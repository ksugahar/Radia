#!/usr/bin/env python
# coding: utf-8

# In[12]:


from netgen.occ import *
from netgen.webgui import Draw as DrawGeo
from ngsolve import *
from ngsolve.webgui import Draw
from math import pi
from numpy import *

import scipy.sparse as sp
import matplotlib.pylab as plt
from scipy.io import savemat

import os, sys
sys.path.append('../../ICCG/JP-MARs/SparseSolv')

import SparseSolvPy
import scipy.sparse as sp


# In[13]:


#tetrahedron
cir = Circle((0, 0), 0.01).Face()
cir.edges.name = 'outer'
geo = OCCGeometry(cir, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=0.001)).Curve(3)
Draw(mesh);

mesh.GetMaterials(), mesh.GetBoundaries()


# In[26]:


Rn=[]
Ln=[]

order = 1
fes = H1(mesh, order=order, dirichlet="outer")

# u and v refer to trial and test-functions in the definition of forms below
u = fes.TrialFunction()
v = fes.TestFunction()


r = 0.01
sigma = 1e6
mu = 4*pi*1e-7
E = CoefficientFunction(1)
J = sigma*E

for nStage in range(2):
    print(nStage+1, "-stage")
    R = 1/Integrate(J*J/sigma*dx, mesh)
    Rn.append(R)
    R_theory = (2*nStage+1)/(pi*r*r*sigma)
    print("R_theory[",2*nStage,"]:", R_theory)
    print("       R[",2*nStage,"]:", R)
    R_err = abs(R - R_theory)/abs(R_theory)
    print("     R_err[",nStage+1,"]:",R_err) 

    #静磁界計算
    a = BilinearForm(fes)
    a += 1/mu*grad(u)*grad(v)*dx
    Precon = Preconditioner(a, type="bddc")
    f = LinearForm(fes)
    f += v*J * dx

    gfA = GridFunction(fes)
    a.Assemble()
    f.Assemble()

    solvers.CG(sol=gfA.vec, rhs=f.vec, mat=a.mat, pre=Precon.mat, tol=1e-8, printrates=True, maxsteps=10000)


    if nStage == 0:
        Apotential = R*(gfA)
        A1 = Apotential
        B = R*grad(gfA)
        B1 = B
    else:
        Apotential = Apotential + R*(gfA)
        A3 = Apotential
        B = B + R*grad(gfA)
        B3 = B

    L = Integrate(R*J*Apotential*dx, mesh)
    Ln.append(L)
    L_theory = (mu)/(4*2*(nStage+1)*pi)
    print("L_theory[",2*nStage+1,"]:", L_theory)
    print("       L[",2*nStage+1,"]:", L)
    L_err = abs(L - L_theory)/abs(L_theory)
    print("     L_err[",2*nStage+1,"]:",L_err) 

    J = J - sigma*Apotential/L

    if nStage == 0:    
        J2 = J   
    else:
        J4 = J

"""    
data = {'Rn': Rn, 'Ln': Ln}
FileName = f"2D_{order}.mat"  
savemat(FileName,data) 
"""

print("A_InnerProduct", Integrate(InnerProduct(A1, A3), mesh))
print("B1_H3_InnerProduct", Integrate(InnerProduct(B1, B3/mu), mesh))
print("J2_J2_InnerProduct", Integrate(InnerProduct(J2, J2), mesh))
print("J2_E4_InnerProduct", Integrate(InnerProduct(J2, J4/sigma), mesh))
print("J_InnerProduct", Integrate(J2*J4, mesh))
print("R2", (1/Integrate(J2*J2/sigma, mesh)))


# In[ ]:




