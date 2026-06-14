"""RETIRED field-recurrence experiment. Only independent setup remains executable.
The following is historical context, not current functionality:

Step 0a v3: A-form Kameari with Tanimoto canonical L = Integrate(B^2/mu dx)
over FULL mesh (magnetic energy).

Difference from v2:
- v2 used L = Integrate(R*J*Apot*dx, conductor) — works for closed PEC,
  fails for air-box because of conductor-surface boundary term.
- v3 uses L = Integrate(B*B/mu*dx, full_mesh) — physical magnetic energy,
  always positive, valid for any geometry.

Reference: CLN_AT.ipynb (Tanimoto thesis, 修論).

In closed PEC: int_par_part B^2 = J*A so v2's formula = v3's formula.
In air-box: gauss-bound integration introduces conductor-surface boundary
term: int_full B^2/mu = int_cond J*A + boundary_term. The boundary term
matters and is what v3 captures.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
print("Starting Step 0a v3 (Tanimoto B^2/mu energy)...", flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, H1, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, grad, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from math import pi
import json, time, argparse
from pathlib import Path

mu0 = 4 * pi * 1e-7
NU_0 = 1.0 / mu0
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

H_COND = 0.5e-3
H_AIR = 2.0e-3
ORDER = 2
N_STAGES = 8
BONUS_INT = 8

ELF_TAU_LEAD = 11.51e-6
TAU_PEC_LIMIT = mu0 * sigma_Cu * ax**2 * ay**2 / (pi**2 * (ax**2 + ay**2))
R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


def build_geo(R_outer):
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.name = "conductor"
    cuboid.maxh = H_COND

    sphere_outer = Sphere(Pnt(0, 0, 0), R_outer)
    for f in sphere_outer.faces:
        f.name = "outer"
    sphere_outer.maxh = H_AIR

    air = sphere_outer - cuboid
    air.name = "air"

    geo = Glue([cuboid, air])
    return OCCGeometry(geo)


def main(R_outer, order=ORDER, n_stages=N_STAGES):
    ngsglobals.msg_level = 0
    print(f"=== Step 0a v3: A-form + air-box + B^2 energy, R_outer = {R_outer*1000:.1f} mm ===", flush=True)
    print(f"  Closed-PEC analytic tau_lead     = {TAU_PEC_LIMIT*1e6:.3f} us", flush=True)
    print(f"  ELF (vacuum) reference tau_lead  = {ELF_TAU_LEAD*1e6:.2f} us", flush=True)
    print(f"  N_STAGES = {n_stages}, ORDER = {order}\n", flush=True)

    geo = build_geo(R_outer)
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(order + 1)
    print(f"  ne = {mesh.ne}  ({time.time()-t0:.1f}s)\n", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0})

    fes = HCurl(mesh, order=order, dirichlet="outer", nograds=True)
    gauge = H1(mesh, order=order, dirichlet="outer")
    print(f"  HCurl ndof = {fes.ndof}, H1 ndof = {gauge.ndof}\n", flush=True)

    u, v = fes.TnT()
    uu, vv = gauge.TnT()

    a_HC = BilinearForm(fes)
    a_HC += NU_0 * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    a_HH = BilinearForm(gauge)
    a_HH += grad(uu) * grad(vv) * dx(bonus_intorder=BONUS_INT)

    print("  Assembling+factor...", flush=True)
    t0 = time.time()
    with TaskManager():
        a_HC.Assemble()
        a_HH.Assemble()
        try:
            inv_HC = a_HC.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_HC = a_HC.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        try:
            inv_HH = a_HH.mat.Inverse(gauge.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_HH = a_HH.mat.Inverse(gauge.FreeDofs(), inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    A_s_cf = CoefficientFunction((-y, x, 0)) * 0.5
    J = sigma_cf * A_s_cf

    R0_check = float(Integrate(J * J * sigma_inv_cf
                               * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    print(f"  R_0 init: {1.0/R0_check:.6e} Ohm vs analytic {R0_anal:.6e}\n", flush=True)

    diag = []
    Apot = None   # CF accumulator for A
    B_pot = None  # CF accumulator for B = curl(A)
    J_history = []  # GridFunctions for orthogonality diagnostic

    # Historical field-recurrence driver retired; assembly above is retained.
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--R_outer_mm", type=float, default=25.0)
    parser.add_argument("--order", type=int, default=ORDER)
    parser.add_argument("--n_stages", type=int, default=N_STAGES)
    args = parser.parse_args()
    main(R_outer=args.R_outer_mm * 1e-3, order=args.order, n_stages=args.n_stages)
