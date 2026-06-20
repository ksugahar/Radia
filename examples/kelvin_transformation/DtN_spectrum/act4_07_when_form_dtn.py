# -*- coding: utf-8 -*-
# DEMO (dd): "If you're going to build the DtN, why not just solve the FEM directly?" -- the honest
# cost answer. For SOLVING ONE configuration the questioner is RIGHT: do NOT form the dense Schur
# complement Lambda -- just solve the monolithic sparse Kelvin-FEM (one factor + one back-solve).
# Forming the dense exterior DtN matrix costs ~N_Gamma back-solves (one per boundary DoF) + dense
# storage O(N_Gamma^2): far more than a single solve. This script MEASURES that gap and states the
# narrow conditions under which forming Lambda actually pays.
#
# WHEN forming Lambda IS worth it (the operator is the DELIVERABLE, not the field of one problem):
#   (1) ANALYSIS / CERTIFICATION: the SPECTRUM of Lambda (the -(n+1)/R ladder) is a predictive
#       datasheet -- you compute a few eigenvalues, you do NOT solve. (This is the paper's actual point;
#       act1_02_multipole_spectrum..s use the spectrum, not a solve.)
#   (2) PORTABILITY: hand the small dense Lambda to ANOTHER code (a commercial FEM/BEM/interior solver)
#       that has NO open-boundary capability -> an exact, Green-function-free open BC as a drop-in block.
#   (3) SUBSTRUCTURING REUSE: a FIXED exterior with a CHANGING interior MATRIX (design sweep, moving
#       part) -> form Lambda once, refactor only the small interior each time. (NB: if only the RHS
#       changes, the MONOLITHIC factorization already amortizes -- Lambda gives nothing. And in Kelvin
#       the exterior is a cheap thin ball, so the reuse win is small unless the interior dominates.)
#   (4) SURROGATE: Lambda IS the sparse-generated stand-in for the dense layered/Sommerfeld Green
#       operator -- valuable when something downstream NEEDS that operator (not when you just want u).
# So the paper must NOT claim "form the DtN to solve faster than FEM" (false). It claims the SPECTRUM
# (datasheet) and the OPERATOR-as-artifact (portable, material-aware, Green-function-free surrogate).
import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import ngsolve as ng
import netgen.occ as occ
from netgen.occ import Sphere, Pnt, Vec, IdentificationType, OCCGeometry, Glue

a, R_out, offset = 0.5, 1.0, 3.0
order, maxh = 3, 0.28

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
A = ng.BilinearForm(mu * ng.grad(u) * ng.grad(v) * ng.dx(bonus_intorder=6)); A.Assemble()

nd = fes.ndof; free = fes.FreeDofs()
gb = fes.GetDofs(mesh.Boundaries("inner"))
gamma = [i for i in range(nd) if gb[i] and free[i]]
rest = [i for i in range(nd) if free[i] and not gb[i]]
freelist = [i for i in range(nd) if free[i]]
r_, c_, val = A.mat.COO()
Asp = sp.csr_matrix((np.array(val), (np.array(r_), np.array(c_))), shape=(nd, nd))
Aff = Asp[freelist][:, freelist].tocsc()
Arr = Asp[rest][:, rest].tocsc()
Arg = Asp[rest][:, gamma].toarray()
Agg = Asp[gamma][:, gamma].toarray()
Agr = Asp[gamma][:, rest]
NG = len(gamma)

# (A) MONOLITHIC: one factor + one back-solve (what you do to SOLVE one configuration) -------------
rhs = np.ones(Aff.shape[0])
t0 = time.time(); lu = spla.splu(Aff); t_fac_mono = time.time() - t0
t0 = time.time(); _ = lu.solve(rhs); t_solve_mono = time.time() - t0
t_mono = t_fac_mono + t_solve_mono

# (B) FORM the dense exterior DtN Lambda (the OPERATOR) -- factor A_rest + N_Gamma back-solves -------
t0 = time.time(); lu_r = spla.splu(Arr); t_fac_L = time.time() - t0
t0 = time.time(); X = lu_r.solve(Arg); Lam = Agg - (Agr @ X); t_back_L = time.time() - t0
t_formL = t_fac_L + t_back_L

print("Cost: SOLVE one configuration (monolithic sparse FEM)  vs  FORM the dense DtN operator Lambda\n")
print("   free DoF = %d,  N_Gamma (boundary DoF) = %d,  Lambda is %dx%d dense (%.1f MB)"
      % (Aff.shape[0], NG, NG, NG, NG * NG * 8 / 1e6))
print("   (A) monolithic solve : factor %.3fs + back-solve %.4fs = %.3fs   <- one open-boundary solve"
      % (t_fac_mono, t_solve_mono, t_mono))
print("   (B) form Lambda      : factor %.3fs + %d back-solves+GEMM %.3fs = %.3fs"
      % (t_fac_L, NG, t_back_L, t_formL))
print("   -> forming the operator costs %.0fx a single solve (the N_Gamma=%d back-solves dominate)\n"
      % (t_formL / max(t_mono, 1e-9), NG))

# when (if ever) the cost of forming Lambda amortizes -- two honest scenarios ----------------------
print("Does forming Lambda ever amortize? two cases, honestly:")
print("   (i)  only the RHS/source changes  -> the MONOLITHIC factorization is reused too (back-solve")
print("        %.4fs each); Lambda gives NOTHING and wastes the %.2fs build. Do NOT form Lambda."
      % (t_solve_mono, t_formL))
print("   (ii) the INTERIOR matrix changes, exterior FIXED (substructuring) -> monolithic must refactor")
print("        the whole system (~%.2fs) each config; condensed refactors only the (smaller) interior +" % t_mono)
print("        reuses the fixed dense Lambda. Pays off IF the exterior factor you avoid re-doing exceeds")
print("        the per-config interior refactor over enough configs. In KELVIN the exterior is a cheap")
print("        thin ball, so that saving is modest -- the reuse case is real but rarely the main reason.")

print("\n=> You are right for SOLVING: do not form the DtN to solve one problem -- the monolithic sparse")
print("   Kelvin-FEM is the cheaper, simpler path (act5_02_layered_green_function/x/aa solve exactly this way). Form Lambda ONLY")
print("   when the OPERATOR is the deliverable: (1) its SPECTRUM as a predictive datasheet (the paper's")
print("   real point -- analysis, not a solve); (2) a PORTABLE open-boundary block for a code that lacks")
print("   one; (3) substructuring reuse with a changing interior; (4) a sparse SURROGATE for the dense")
print("   layered/Sommerfeld Green operator. 'Form the DtN to solve faster than FEM' is NOT a claim.")
