"""Step 0a: A-formulation Kameari Cauer ladder for cuboid 5x2x1 + uniform B_z,
finite air-box with Dirichlet outer BC (NO Kelvin).

Sanity check before T-Omega + Kelvin (Step 1):
- Closed PEC analytic limit (R_outer -> 0): tau_lead = 25.46 us
  (TE_z(1,1,0) of cuboid surrounded by PEC; v3.py Case B reproduces this
   as max tau_n over 12 stages = 25.33 us)
- Vacuum limit (R_outer -> inf, ELF reference): tau_lead = 11.51 us
- This Step 0a sweeps R_outer to see where tau_lead ends up.

Drive: J_imp_0 = sigma * A_s where A_s = (B_0/2)(-y, x, 0) (uniform B_z=1).

Goal: verify the formulation is sane before adding Kelvin / switching to T-Omega.

Usage:
    python cuboid_521_step0a_airbox_no_kelvin.py --R_outer_mm 25
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
sys.path.insert(0, r"S:/Radia/01_GitHub/src/radia")
print("Starting Step 0a (A-form + finite air-box, no Kelvin)...", flush=True)

from netgen.occ import Box, Sphere, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals, BitArray,
)
from collections import deque
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
GAUGE_EPS = 1e-8

ELF_TAU_LEAD = 11.51e-6   # seconds (4-stage Cauer reference)
TAU_PEC_LIMIT = mu0 * sigma_Cu * ax**2 * ay**2 / (pi**2 * (ax**2 + ay**2))
R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


def build_spanning_tree(mesh):
    nv = mesh.nv
    visited = [False] * nv
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        adj[v0].append((v1, ed.nr))
        adj[v1].append((v0, ed.nr))
    visited[0] = True
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for vn, edn in adj[v]:
            if not visited[vn]:
                visited[vn] = True
                tree_edges.append(edn)
                queue.append(vn)
    return tree_edges


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


def main(R_outer):
    ngsglobals.msg_level = 0
    print(f"=== Step 0a: A-form + air-box (NO Kelvin), R_outer = {R_outer*1000:.1f} mm ===", flush=True)
    print(f"  Closed-PEC analytic tau_lead   = {TAU_PEC_LIMIT*1e6:.3f} us", flush=True)
    print(f"  ELF (vacuum) reference tau_lead = {ELF_TAU_LEAD*1e6:.2f} us", flush=True)
    print(f"  Analytic R_0 (Case B)          = {R0_anal:.6e} Ohm", flush=True)
    print(f"  N_STAGES = {N_STAGES}, ORDER = {ORDER}\n", flush=True)

    geo = build_geo(R_outer)
    print("Generating mesh...", flush=True)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER + 1)
    print(f"  ne = {mesh.ne}  ({time.time()-t0:.1f}s)", flush=True)
    print(f"  materials = {mesh.GetMaterials()}", flush=True)
    print(f"  boundaries = {mesh.GetBoundaries()}\n", flush=True)

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0})

    fes = HCurl(mesh, order=ORDER, dirichlet="outer", nograds=True)
    print(f"  HCurl ndof = {fes.ndof}", flush=True)

    tree_edges = build_spanning_tree(mesh)
    fd = BitArray(fes.FreeDofs())
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  tree-cotree masked {masked} edges, active = {sum(fd)}\n", flush=True)

    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += NU_0 * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)
    a += GAUGE_EPS * NU_0 * u * v * dx(bonus_intorder=BONUS_INT)

    print("  Assembling+factor...", flush=True)
    t0 = time.time()
    with TaskManager():
        a.Assemble()
        try:
            inv = a.mat.Inverse(fd, inverse="pardiso")
        except Exception:
            inv = a.mat.Inverse(fd, inverse="sparsecholesky")
    print(f"    {time.time()-t0:.1f}s\n", flush=True)

    A_s_cf = CoefficientFunction((-y, x, 0)) * 0.5
    J_imp_cf = sigma_cf * A_s_cf

    R0_check = float(Integrate(sigma_cf * A_s_cf * A_s_cf
                               * dx("conductor", bonus_intorder=BONUS_INT), mesh))
    R0_init = 1.0 / R0_check if R0_check > 1e-30 else float('inf')
    print(f"  R_0 init: {R0_init:.6e} Ohm vs analytic {R0_anal:.6e}, "
          f"ratio {R0_init/R0_anal:.6f}\n", flush=True)

    diag = []
    gfApot = GridFunction(fes)
    gfApot.vec[:] = 0.0

    for n in range(N_STAGES):
        t_stage = time.time()

        f = LinearForm(fes)
        f += J_imp_cf * v * dx("conductor", bonus_intorder=BONUS_INT)
        with TaskManager():
            f.Assemble()
        gfA_n = GridFunction(fes)
        with TaskManager():
            gfA_n.vec.data = inv * f.vec

        R_inv = float(Integrate(J_imp_cf * J_imp_cf * sigma_inv_cf
                                * dx("conductor", bonus_intorder=BONUS_INT),
                                mesh))
        if R_inv < 1e-30:
            print(f"  [{n}] R_inv too small ({R_inv:.2e}), STOP", flush=True)
            break
        R_n = 1.0 / R_inv

        gfApot.vec.data += R_n * gfA_n.vec

        Ln_int = float(Integrate(J_imp_cf * gfApot
                                 * dx("conductor", bonus_intorder=BONUS_INT),
                                 mesh))
        L_n = R_n * Ln_int

        signL = "+" if L_n > 0 else ("-" if L_n < 0 else "0")
        tau_us = L_n / R_n * 1e6 if R_n > 0 and L_n > 0 else None
        elapsed = time.time() - t_stage
        tau_str = f"{tau_us:.4f}" if tau_us else "N/A"
        print(f"  [{n}] R={R_n:.4e}, L={L_n:.4e}({signL}), "
              f"tau=L/R={tau_str}us  ({elapsed:.1f}s)", flush=True)
        diag.append({"n": n, "R_n": R_n, "L_n": L_n, "tau_us": tau_us,
                     "energy_norm_inv": R_inv, "sign_L": signL})

        if L_n <= 0:
            print(f"      ** SIGN FLIP at n={n} **", flush=True)
            break

        J_imp_cf = J_imp_cf - sigma_cf * gfApot * (1.0 / L_n)

    print("\n" + "=" * 78, flush=True)
    print(f"STEP 0a (A-form + air-box, R_outer={R_outer*1000:.1f} mm) RESULTS", flush=True)
    print("=" * 78, flush=True)
    for d in diag:
        tau_s = f"{d['tau_us']:.4f}" if d['tau_us'] else "N/A"
        print(f"  n={d['n']}: R={d['R_n']:.4e}, L={d['L_n']:.4e}({d['sign_L']}), "
              f"tau={tau_s} us", flush=True)

    valid_taus = [d['tau_us'] for d in diag if d['tau_us'] is not None]
    if valid_taus:
        tau_max = max(valid_taus)
        ratio_pec = tau_max / (TAU_PEC_LIMIT * 1e6)
        ratio_elf = tau_max / (ELF_TAU_LEAD * 1e6)
        print(f"\n  max tau_n = {tau_max:.4f} us", flush=True)
        print(f"    vs closed-PEC limit {TAU_PEC_LIMIT*1e6:.2f} us  ratio {ratio_pec:.3f}", flush=True)
        print(f"    vs ELF (vacuum)     {ELF_TAU_LEAD*1e6:.2f} us  ratio {ratio_elf:.3f}", flush=True)

    out = {"method": "A-form + finite air-box (NO Kelvin), Case B drive",
           "R_outer_mm": R_outer*1000,
           "h_cond_mm": H_COND*1000, "h_air_mm": H_AIR*1000,
           "order": ORDER, "ne": mesh.ne, "ndof": fes.ndof,
           "stages": diag, "R0_analytic": R0_anal,
           "tau_pec_limit_us": TAU_PEC_LIMIT*1e6,
           "tau_elf_us": ELF_TAU_LEAD*1e6}
    out_path = Path(__file__).parent / f"cuboid_521_step0a_Router_{int(R_outer*1000):03d}mm.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--R_outer_mm", type=float, default=25.0)
    parser.add_argument("--order", type=int, default=ORDER)
    parser.add_argument("--n_stages", type=int, default=N_STAGES)
    args = parser.parse_args()
    ORDER = args.order
    N_STAGES = args.n_stages
    main(R_outer=args.R_outer_mm * 1e-3)
