#!/usr/bin/env python
# coding: utf-8

# In[3]:


from netgen.occ import *
from netgen.webgui import Draw as DrawGeo
from ngsolve import *
from ngsolve.webgui import Draw
from math import pi
from numpy import *

from ngsolve.krylovspace import CGSolver

import scipy.sparse as sp
import matplotlib.pylab as plt
from scipy.io import savemat
import scipy.sparse as sp

import os, sys
sys.path.append('../../../ICCG/JP-MARs/SparseSolv')

import SparseSolvPy
#dir(SparseSolvPy.MatSolvers)


# In[4]:


#tetrahedron
cyl = Cylinder(Pnt(0,0,0), Z, r=1, h=0.2).bc("S_side")
cyl.mat("cyl")

cyl.faces.Min(Z).name="S_bottom"
cyl.faces.Max(Z).name="S_top"

geo = OCCGeometry(cyl)
mesh = Mesh(geo.GenerateMesh(maxh=0.1))
mesh.Curve(5)
Draw (mesh)#, clipping=True);
mesh.GetMaterials(), mesh.GetBoundaries()

#help(mesh)
edges_no=mesh.nedge
faces_no=mesh.nface
elements_no=mesh.nv
calcd_ndof = (mesh.nedge*3)+( mesh.nface*6)+( mesh.ne*3)
print('#nv=', mesh.nv)
print('#nedge=', mesh.nedge)
print('#nface=', mesh.nface)
print('#ne=', mesh.ne)
print(calcd_ndof)


# In[5]:


cyl = Cylinder ((0,0,0), Z, r=1, h=0.2, bottom="bot", top="top").bc("S_side")
cyl.faces.Min(Z).name="S_bottom"
cyl.faces.Max(Z).name="S_top"
cyl.faces.Min(Z).Identify(cyl.faces.Max(Z), "bot-top", type=IdentificationType.CLOSESURFACES)
Draw(cyl);
ngmesh = OCCGeometry(cyl).GenerateMesh(maxh=0.1)#, quad_dominated=True)
ngmesh.ZRefine("bot-top", [0.25,0.5,0.75])
mesh = Mesh(ngmesh)
mesh.Curve(5)
Draw(mesh);
mesh.GetMaterials(), mesh.GetBoundaries()


# In[1]:


Rn=[]
Ln=[]
log1min=[]
ndof=[]
nstage=[]

method = "accICCG"
mesh_type_size="tet_0.1"
order=1
#fes = HCurl(mesh, order=order, dirichlet="S_side|S_bottom|S_top", nograds = True)
fes = HCurl(mesh, order=order, dirichlet="S_side|S_bottom|S_top", type1 = True)
print ("Hcurl_ndof =", fes.ndof)
gauge = H1(mesh, order=order, dirichlet="S_side|S_bottom|S_top")
print ("H1_ndof =", gauge.ndof)
print ("ndof =", fes.ndof+gauge.ndof)
ndof.append(fes.ndof+gauge.ndof)


# u and v refer to trial and test-functions in the definition of forms below
u = fes.TrialFunction()
v = fes.TestFunction()

uu = gauge.TrialFunction()
vv = gauge.TestFunction()

r = 1
h = 0.2
sigma = 3.4e7
mu = 4*pi*1e-7
E = CoefficientFunction((0,0,1))
J = sigma*E

tol = 1e-12
max_iter = 10000
acc = [1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01]



# Historical field recurrence and dependent outputs retired.
# Mesh, initial field and independent solver preparation above are retained.

# Independent closed-form scalar reference expressions retained from this version.
# nStage / nn denote the scalar index; these are not a field-update algorithm.
# R_theory = (2*nStage+1)/(pi*r*r*sigma*h)
# L_theory = (mu)/(4*2*(nStage+1)*pi*h)
# R0 = (2*nn+1)/pi/a**2/sig/h
# L0 =mur*mu0/(4*(2*nn+2))/pi/h
