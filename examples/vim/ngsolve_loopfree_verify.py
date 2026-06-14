"""
ngsolve_loopfree_verify.py -- REAL-NGSolve confirmation that the FEEC magnetization
"loops" are FIELD-NULL by construction on DISTORTED hex elements (the actual elements a
C++ VIM would interface with -- not symbolic Mathematica).

This is the Python companion to
  packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/vim_loopfree.wls
which proves the same result symbolically.  Run with the lab's Python 3.12 + NGSolve
6.2.2604:   python ngsolve_loopfree_verify.py

ANSWER to "does NGSolve's basis solve the loop-mode problem?":
  YES for the loop-mode itself (the distorted-element + high-mu_r FIELD defect).
  A magnetization M is FIELD-NULL (a "loop", in ker N) IFF it is CHARGE-FREE:
        rho = -div M = 0   AND   sigma = M.n = 0.
  The FEEC loop  M = curl(interior H(curl))  (interior = zero tangential trace) has BOTH:
        div(curl W) = 0                 -> no volume charge
        (curl W).n = curlS(W_t) = 0     -> no surface charge   (W interior => W_t = 0)
  => charge-free => field-null EVERYWHERE, at ANY mu_r, on ANY (distorted) element, because
  the contravariant Piola map preserves both div and normal-trace (de Rham commuting diagram).

  This does NOT address the SEPARATE high-mu_r CONDITIONING problem (iters ~ mu_r^1.5,
  the >15^3 H-ILU factor wall): that is the near-null SPECTRUM of A=(1/chi)I-N -> -N,
  formulation-independent (hex ~ tet), not a kernel/loop issue.  See memory/.
"""
import json, os
import numpy as np
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
from math import sin, pi

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))
results = {}

# genuinely distorted (non-affine) hex mesh on ~[0,1]^3 -- the de Rham defect regime
def distort(x, y, z):
    return (x + 0.12 * sin(pi * y) * z, y + 0.10 * x * z, z + 0.08 * x * sin(pi * y))

with ng.TaskManager():
    mesh = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3, mapping=distort)
    print("mesh: ne =", mesh.ne, "(distorted hexes), nv =", mesh.nv)

    order = 3
    fes = ng.HCurl(mesh, order=order)
    gfw = ng.GridFunction(fes)
    free, bnd = fes.FreeDofs(), fes.GetDofs(mesh.Boundaries(".*"))
    interior = [i for i in range(fes.ndof) if free[i] and not bnd[i]]
    print(f"H(curl) order {order}: ndof={fes.ndof}, interior dofs={len(interior)}")
    rng = np.random.default_rng(0)
    gfw.vec[:] = 0.0
    for i in interior:
        gfw.vec[i] = rng.standard_normal()       # a generic interior H(curl) field
    M = ng.curl(gfw)                              # loop candidate (in H(div))

    fesH = ng.HDiv(mesh, order=order)
    gfM = ng.GridFunction(fesH); gfM.Set(M)
    vol = ng.Integrate(ng.CF(1.0), mesh)
    Mnorm = np.sqrt(ng.Integrate(M * M, mesh) / vol)
    div_rel = np.sqrt(ng.Integrate(ng.div(gfM) ** 2, mesh) / vol) / Mnorm
    nn = ng.specialcf.normal(mesh.dim)
    area = ng.Integrate(ng.CF(1.0), mesh, definedon=mesh.Boundaries(".*"))
    nt_rel = np.sqrt(ng.Integrate((M * nn) ** 2, mesh, definedon=mesh.Boundaries(".*")) / area) / Mnorm

    # field-null: CHARGE FORM (no cancellation; the field comes only from the ~0 charges)
    rho = -ng.div(gfM)
    def Hcharge(obs, io=10):
        R = ng.CF((obs[0]-ng.x, obs[1]-ng.y, obs[2]-ng.z)); Rn = ng.sqrt(R*R)
        vt = np.array([ng.Integrate(rho*R[k]/Rn**3, mesh, order=io) for k in range(3)])
        st = np.array([ng.Integrate((gfM*nn)*R[k]/Rn**3, mesh,
                       definedon=mesh.Boundaries(".*"), order=io) for k in range(3)])
        return (vt + st) / (4*pi)
    ext = [(2.0,0.3,0.4),(0.5,2.0,-0.6),(-1.5,0.5,1.5)]
    ref = 1.0/(4*pi*1.5**3)
    Hmax = max(np.linalg.norm(Hcharge(p)) for p in ext)

    results.update(div_rel=float(div_rel), normaltrace_rel=float(nt_rel),
                   Hcharge_max=float(Hmax), Hcharge_over_ref=float(Hmax/ref))
    print(f"(i)   ||div M||/||M||      = {div_rel:.3e}   (volume charge rho)")
    print(f"(ii)  ||M.n||_bnd/||M||    = {nt_rel:.3e}   (surface charge sigma)")
    print(f"(iii) charge-form |H|_max  = {Hmax:.3e}  (ref dipole d=1.5 ~ {ref:.2e}) -> field-null")

verdict = results["div_rel"] < 1e-7 and results["normaltrace_rel"] < 1e-7 and results["Hcharge_over_ref"] < 1e-6
results["verdict_loop_field_null_on_distorted"] = bool(verdict)
print(f"\n=> FEEC loop CHARGE-FREE => FIELD-NULL on DISTORTED hexes: {'YES' if verdict else 'NO'}")
with open(os.path.join(HERE, "ngsolve_loopfree_verify.json"), "w") as f:
    json.dump(results, f, indent=2)
print("saved", os.path.join(HERE, "ngsolve_loopfree_verify.json"))
