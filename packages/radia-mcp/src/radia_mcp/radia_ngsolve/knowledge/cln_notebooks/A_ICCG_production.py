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


for nStage in range(10):
    print(nStage+1, "-stage")
    nstage.append(nStage+1)
    R = 1/Integrate(J*J/sigma*dx, mesh)
    Rn.append(R)
    R_theory = (2*nStage+1)/(pi*r*r*sigma*h)
    print("R_theory[",2*nStage,"]:", R_theory)
    print("       R[",2*nStage,"]:", R)
    R_err = abs(R - R_theory)/abs(R_theory)
    print("     R_err[",nStage+1,"]:",R_err) 

    #静磁界計算
    a = BilinearForm(fes)
    a += 1/mu*curl(u)*curl(v)*dx
    #a += 1e-6/mu*u*v*dx

    f = LinearForm(fes)
    f += v*J * dx

    gfA = GridFunction(fes)
    with TaskManager():
        a.Assemble()
        f.Assemble()

    A = sp.csr_matrix (a.mat.CSR())
    #print(A)
    Acut = A[:,fes.FreeDofs()][fes.FreeDofs(),:]
    fcut = array(f.vec.FV())[fes.FreeDofs()]
    ucut = array(f.vec.FV(), copy=True)[fes.FreeDofs()]

    rows, cols = Acut.nonzero()
    vals = Acut[rows, cols]
    vals = ravel(vals)
    mat = SparseSolvPy.SparseMat(len(rows), rows, cols, vals)

    solver = SparseSolvPy.MatSolvers()
    solver.setSaveBest(True)
    solver.setSaveLog(True)
    solver.setDiagScale(False)
    solver.setDirvegeType(1)
    solver.setBadDivCount(1000)
    solver.setBadDivVal(100.0)
    solver.solveICCG_py(len(fcut), tol, max_iter, acc[nStage], mat, fcut, ucut, True)


    log1 = solver.getResidualLog_py()
    print(log)

    plt.plot(range(len(log1)), log1)    
    plt.yscale('log')
    plt.show(block=False)  

    array(gfA.vec.FV(), copy=False)[fes.FreeDofs()] = ucut
    print("min:", min(log1))
    #log1min.append(min(log1))

    result = Acut.dot(ucut) - fcut
    norm = linalg.norm(result)/linalg.norm(fcut)
    print("結果のノルム:", norm)


    #補正計算 
    aa = BilinearForm(gauge)
    aa += grad(uu)*grad(vv)*dx
    cc = Preconditioner(aa, "local")
    #cc = Preconditioner(aa, "bddc", coarsetype="h1amg")
    ff = LinearForm(gauge)
    ff += grad(vv)*gfA*dx
    aa.Assemble()
    ff.Assemble()        
    gfu = GridFunction(gauge)
    solvers.CG(sol=gfu.vec, rhs=-ff.vec, mat=aa.mat, pre=cc.mat, tol=1e-12, printrates='\r', maxsteps=10000)
    """
    solver = CGSolver(aa.mat, cc.mat, tol=1e-12, maxiter=10000)
    gfu.vec.data = solver * -ff.vec
    #print("iteration count =", solver.iterations)
    #print("convergence history = ", solver.errors)
    plt.plot(range(len(solver.errors)), solver.errors) 
    plt.yscale('log')
    plt.show(block=False)     
    #print("min:", min(solver.errors))
    """

    if nStage == 0:
        Apotential = R*(gfA+grad(gfu))
        #Apotential = R*gfA
        B = R*curl(gfA)
    else:
        Apotential = Apotential + R*(gfA+grad(gfu))
        #Apotential = Apotential + R*gfA
        B = B + R*curl(gfA)

    L = Integrate(R*J*Apotential*dx, mesh)
    Ln.append(L)
    L_theory = (mu)/(4*2*(nStage+1)*pi*h)
    print("L_theory[",2*nStage+1,"]:", L_theory)
    print("       L[",2*nStage+1,"]:", L)
    L_err = abs(L - L_theory)/abs(L_theory)
    print("     L_err[",nStage+1,"]:",L_err) 

    J = J - sigma*Apotential/L

    """
    Draw (Apotential*h, mesh, "A-field", draw_surf=False, \
          clipping = { "z" : -0.01, "function":True}, \
          vectors = { "grid_size":50});        
    Draw (J, mesh, "J-field", draw_surf=False, \
          clipping = { "z" : -0.01, "function":True}, \
          vectors = { "grid_size":50});    
    Draw (B*h, mesh, "B-field", draw_surf=False, \
          clipping = { "z" : -0.01, "function":True}, \
          vectors = { "grid_size":50});   

    """
    N = 1000
    xpnts = [2*i/N-1 for i in range(N+1)]
    Ah = Apotential*h
    vals_A = [Ah(mesh(xi,0,0.1)) for xi in xpnts]
    plt.plot(xpnts,vals_A)
    plt.show(block=False)

    vals_J = [J(mesh(xi,0,0.1)) for xi in xpnts]
    plt.plot(xpnts,vals_J)
    plt.show(block=False)

    Bh = B*h
    vals_B = [Bh(mesh(xi,0,0.1)) for xi in xpnts]
    plt.plot(xpnts,vals_B)
    plt.show(block=False)

    """
    data = {'Rn': Rn, 'Ln': Ln}
    #FileName = f"A_{order}_{mesh_type_size}_{method}_nograds.mat"
    FileName = f"A_{order}_{mesh_type_size}_{method}_type1.mat"  
    savemat(FileName,data)
    """


# In[4]:


a = 1
h = 0.2
sig = 3.4e7
mur = 4*pi*1e-7
mu0 = 1

#print(Rn)
N = shape(Rn)[0]
#print(N)
nn = array(list(range(N)))
#print(nn)
R0 = (2*nn+1)/pi/a**2/sig/h
L0 =mur*mu0/(4*(2*nn+2))/pi/h
#print((Rn-R0)/Rn)
#print((Ln-L0)/Ln)
plt.plot(2*nn, abs(Rn-R0)/R0, ".:", linewidth=0.5)
plt.plot(2*nn+1, abs(Ln-L0)/L0, ".:", linewidth=0.5)
#plt.plot([0,20],[1e-0,1e-0],'k--',linewidth=1)
plt.yscale('log')
plt.xticks(range(0, 2*N, 1))
plt.show(block=False)


# In[ ]:




