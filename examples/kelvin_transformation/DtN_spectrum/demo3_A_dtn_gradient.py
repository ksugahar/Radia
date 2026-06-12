# -*- coding: utf-8 -*-
# DEMO 3 (verified PARTIAL): A-formulation DtN gradient block = -(n+1)/R (formulation-independence).
# Part 1: scalar exterior DtN ladder -(n+1)/R (dipole -2/R).
# Part 2a: A-side B.n = curl_Gamma(A_t) of a magnetic dipole == n=1 zonal harmonic.
# Part 2b: that n=1 harmonic fed through the scalar exterior DtN returns -2/R.
# (Assembling a full independent vector H(curl) exterior DtN operator is a TODO.)
import numpy as np
import ngsolve as ng
from ngsolve import GridFunction, BilinearForm, LinearForm, ds, Integrate, specialcf
from netgen.occ import Sphere, Pnt, OCCGeometry
from radia_mcp.radia_ngsolve.bem_integral import (
    exterior_dtn_spectrum, laplace_exterior_dtn, _dense_mat,
)
R = 1.0; m = 1.0

print("PART 1 -- scalar exterior Laplace DtN spectrum  lambda_n = -(n+1)/R")
spec = exterior_dtn_spectrum(R=R, maxh=0.5, order=1, intorder=10, nmax=3)
for md in spec["modes"]:
    print(f"  n={md['n']} exact={md['lambda_exact']:.6f} bem={md['lambda_mean']:.6f} rel_err={md['rel_err']:.3e}")

print("\nPART 2a -- A-side: B.n of dipole = curl_Gamma(A_t)")
mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=0.3)).Curve(4)
x, y, z = ng.x, ng.y, ng.z
r = ng.sqrt(x * x + y * y + z * z)
pref = m / (4.0 * np.pi * r ** 3)
A_vec = ng.CoefficientFunction((-pref * y, pref * x, 0.0))   # A = m x r/(4pi r^3), m=(0,0,m)
nrm = specialcf.normal(3)
A_t_cf = A_vec - (A_vec * nrm) * nrm
fesH = ng.HCurl(mesh, order=4, definedon=mesh.Boundaries(".*"))
gfA = GridFunction(fesH); uH, vH = fesH.TnT()
massH = BilinearForm(uH.Trace() * vH.Trace() * ds(bonus_intorder=10)).Assemble().mat
rhsH = LinearForm(A_t_cf * vH.Trace() * ds(bonus_intorder=10)).Assemble()
gfA.vec.data = massH.Inverse(freedofs=fesH.FreeDofs()) * rhsH.vec
Bn_curl = ng.InnerProduct(ng.curl(gfA), nrm)
Bn_exact = 2.0 * m * (z / r) / (4.0 * np.pi * R ** 3)        # n=1 zonal harmonic
err = Integrate((Bn_curl - Bn_exact) ** 2 * ds(bonus_intorder=10), mesh).real ** 0.5
ref = Integrate(Bn_exact ** 2 * ds(bonus_intorder=10), mesh).real ** 0.5
print(f"  ||curl_Gamma(A_t)-B.n||={err:.3e}  ||B.n||={ref:.3e}  rel={err/ref:.3e}")

print("\nPART 2b -- A-side n=1 harmonic -> scalar exterior DtN eigenvalue (-> -2/R)")
mesh2 = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=0.4)).Curve(3)
fesS = ng.SurfaceL2(mesh2, order=1, dual_mapping=True); uS, vS = fesS.TnT()
dtn = laplace_exterior_dtn(fesS, intorder=10)
massS = BilinearForm(uS.Trace() * vS.Trace() * ds(bonus_intorder=6)).Assemble().mat
M_d = _dense_mat(massS, fesS.ndof)
rhs_h = LinearForm((ng.z / R) * vS.Trace() * ds(bonus_intorder=6)).Assemble()
gf_h = GridFunction(fesS)
gf_h.vec.data = massS.Inverse(freedofs=fesS.FreeDofs()) * rhs_h.vec
h_vec = gf_h.vec.FV().NumPy()[:]
sigma_vec = dtn @ h_vec
lam = float(np.real(sigma_vec.conj() @ M_d @ h_vec) / np.real(h_vec.conj() @ M_d @ h_vec))
print(f"  scalar DtN eig on n=1 harmonic = {lam:.6f}  (exact -2/R)  rel_err={abs(lam+2.0)/2.0:.3e}")
