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


# Historical field recurrence and dependent outputs retired.
# Mesh, initial field and independent solver preparation above are retained.

# Independent closed-form scalar reference expressions retained from this version.
# nStage / nn denote the scalar index; these are not a field-update algorithm.
# R_theory = (2*nStage+1)/(pi*r*r*sigma)
# L_theory = (mu)/(4*2*(nStage+1)*pi)
