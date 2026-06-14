# -*- coding: utf-8 -*-
# DEMO (cc): "Is the condensed Kelvin operator still FEM, or has it become BEM?" -- the taxonomy,
# grounded by one measurement. ANSWER: it is reduced/condensed FEM (static condensation = macro-element
# = superelement; SBFEM family), NOT BEM. The single defining criterion that separates BEM from the
# FEM family is: DOES IT USE THE FUNDAMENTAL SOLUTION (Green's function)?  BEM does (boundary integral
# of G -> dense, surface-only mesh, singular quadrature, G required). The Kelvin route does NOT (weak
# form, local polynomial basis, sparse volume assembly) -- it only LOOKS like BEM AFTER condensation,
# because both produce the SAME boundary object: the exterior Dirichlet-to-Neumann (Steklov-Poincare)
# operator, which is method-agnostic (demo_r/demo_s: they converge to one operator). The literature
# already named this convergence: Wolf & Song's SBFEM is literally subtitled "a FUNDAMENTAL-SOLUTION-
# LESS boundary-element method" (wolfsong2001) -- the paradoxical name IS this very observation.
#
# This script measures the generator to make the distinction concrete on one sphere Gamma:
#   (1) the Kelvin-FEM VOLUME matrix A is SPARSE (local grad-grad, ~tens nnz/row, no G, SPD);
#   (2) static condensation (Schur) onto Gamma turns it DENSE -- and that dense block IS a DtN matrix
#       whose dipole eigenvalue equals the analytic 2/R, i.e. the same operator BEM would assemble;
#   (3) NO fundamental solution and NO singular quadrature were ever used.
# => same OPERATOR as BEM, completely different GENERATOR. "FEM vs BEM" is the wrong axis; the right one
#    is "two routes to the exterior DtN: the Green-function route (BEM) vs the transformed sparse-volume
#    route (Kelvin-FEM)" -- and only the latter needs no G, stays SPD, and works for arbitrary mu(x).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import eigh
import ngsolve as ng
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry, Glue

a, R_out, offset = 0.5, 1.0, 3.0
order, maxh = 3, 0.30

inner = Sphere(Pnt(0, 0, 0), a); outer = Sphere(Pnt(0, 0, 0), R_out)
for f in inner.faces: f.name = "inner"
for f in outer.faces: f.name = "kelvin_int"
ext = (outer - inner); ext.mat("vac")
kball = Sphere(Pnt(offset, 0, 0), R_out)
for f in kball.faces: f.name = "kelvin_ext"
kball.mat("kelvin")
gnd = occ.Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"
tr = occ.gp_Trsf.Translation(Vec(offset, 0, 0))
[f for f in ext.faces if f.name == "kelvin_int"][0].Identify(
    [f for f in kball.faces if f.name == "kelvin_ext"][0], "kelvin", IdentificationType.PERIODIC, tr)
mesh = ng.Mesh(OCCGeometry(Glue([ext, kball, gnd])).GenerateMesh(maxh=maxh)).Curve(order + 1)
x, y, z = ng.x, ng.y, ng.z; rp2 = (x - offset) ** 2 + y * y + z * z + 1e-20
mu = mesh.MaterialCF({"vac": 1.0, "kelvin": R_out ** 2 / rp2}, default=1.0)
fes = ng.Periodic(ng.H1(mesh, order=order, dirichlet="GND")); u, v = fes.TnT()
A = ng.BilinearForm(mu * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=8)); A.Assemble()   # LOCAL grad-grad, no G
Mf = ng.BilinearForm(u * v * ng.ds("inner", bonus_intorder=8)); Mf.Assemble()

nd = fes.ndof; free = fes.FreeDofs(); gb = fes.GetDofs(mesh.Boundaries("inner"))
gamma = [i for i in range(nd) if gb[i] and free[i]]
rest = [i for i in range(nd) if free[i] and not gb[i]]


def csr(m):
    r, c_, val = m.COO()
    return sp.csr_matrix((np.array(val), (np.array(r), np.array(c_))), shape=(nd, nd))
Asp = csr(A.mat)
nnz_vol = Asp.nnz; fill_vol = nnz_vol / float(nd * nd); nnzrow = nnz_vol / nd
Agg = Asp[gamma][:, gamma].toarray(); Agr = Asp[gamma][:, rest]
Arg = Asp[rest][:, gamma]; Arr = Asp[rest][:, rest].tocsc()
Lam = Agg - (Agr @ spla.spsolve(Arr, Arg.toarray())); Lam = 0.5 * (Lam + Lam.T)
ng_ = len(gamma)
fill_Lam = np.count_nonzero(np.abs(Lam) > 1e-12 * np.max(np.abs(Lam))) / float(ng_ * ng_)
Mgg = csr(Mf.mat)[gamma][:, gamma].toarray(); Mgg = 0.5 * (Mgg + Mgg.T)
spec = np.sort(eigh(Lam, Mgg, eigvals_only=True))
dipole = spec[1:4].mean()

print("ONE sphere Gamma (R=%.2f): the GENERATOR, measured\n" % a)
print(" Kelvin-FEM VOLUME matrix A : %d DoF, %.1f nnz/row, fill %.2f%%  -> SPARSE, SPD, LOCAL grad-grad, NO Green fn"
      % (nd, nnzrow, 100 * fill_vol))
print(" Schur condensation onto Gamma -> dense DtN matrix Lambda : %dx%d, fill %.0f%%  -> DENSE (looks BEM-like)"
      % (ng_, ng_, 100 * fill_Lam))
print(" Lambda dipole eigenvalue = %.4f   vs analytic (n+1)/R = %.4f   rel %.1e   (SAME operator as BEM)"
      % (dipole, 2.0 / a, abs(dipole - 2.0 / a) / (2.0 / a)))

print("\n THE TAXONOMY (what actually separates the methods)")
print("   feature                 |  BEM                      |  reduced/condensed Kelvin-FEM (this)")
print("   ------------------------+---------------------------+-------------------------------------")
print("   fundamental solution G  |  REQUIRED (kernel)        |  NEVER used")
print("   what is meshed          |  surface Gamma only       |  volume (then interior eliminated)")
print("   raw system matrix       |  dense (nnz=DoF^2)        |  SPARSE (%.0f nnz/row), then condensed" % nnzrow)
print("   quadrature              |  singular (1/r kernels)   |  regular, local polynomial")
print("   definiteness            |  often indefinite/complex |  SPD")
print("   exterior material mu(x) |  needs layered/Sommerfeld G or a volume IE | arbitrary mu(x), native (demo_t/bb)")
print("   the boundary product    |  the SAME exterior DtN (Steklov-Poincare) operator -- method-agnostic")

print("\n=> it is STILL FEM: a condensed FE substructure (static condensation = macro-element/superelement;")
print("   SBFEM = 'fundamental-solution-less BEM', wolfsong2001). It is NOT BEM -- no Green's function is")
print("   ever formed. The dense boundary block is the SAME DtN OPERATOR BEM targets, produced by a")
print("   DIFFERENT generator. So 'reduced-FEM vs BEM' blurs only at the OUTPUT (the shared DtN); they")
print("   stay distinct at the GENERATOR (Green-function integral vs sparse weak-form volume + condensation).")
