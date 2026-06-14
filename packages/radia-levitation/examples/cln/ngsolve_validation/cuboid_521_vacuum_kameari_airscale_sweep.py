"""Air-box scale dependence of Kameari Stage 0 tau.

Hypothesis: tau_0 of Kameari with air-box truncation depends on AIR_SCALE
(grows with box size due to PEC cavity modes), so it does NOT converge to
ELF MAGIC BEM reference tau_lead = 11.51 us.

Run only Stage 0 (cheap) for AIR_SCALE in [3, 5, 8, 12, 16].
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from netgen.occ import Box, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from collections import deque
from math import pi
import json, time
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az
H_COND = 0.30e-3
ORDER = 2

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


def build_geo(air_scale, h_air):
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = H_COND
    Bx, By, Bz = air_scale * ax / 2, air_scale * ay / 2, air_scale * az / 2
    air = Box(Pnt(-Bx, -By, -Bz), Pnt(Bx, By, Bz))
    air.mat("air").bc("outer_box")
    air.maxh = h_air
    return OCCGeometry(Glue([air - cuboid, cuboid]))


def stage0_tau(air_scale, h_air):
    print(f"\n=== AIR_SCALE = {air_scale}, h_air = {h_air*1000} mm ===")
    geo = build_geo(air_scale, h_air)
    mesh = Mesh(geo.GenerateMesh(maxh=h_air))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}")

    fes = HCurl(mesh, order=ORDER, dirichlet="outer_box", nograds=True)
    print(f"  HCurl ndof = {fes.ndof}")

    tree_edges = build_spanning_tree(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0})
    u, v = fes.TnT()
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    a_form = BilinearForm(fes)
    a_form += (1.0/mu0) * curl(u) * curl(v) * dx
    t0 = time.time()
    with TaskManager():
        a_form.Assemble()
        inv = a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    print(f"  factor time = {time.time()-t0:.1f}s")

    J = sigma_cf * A_ext
    f = LinearForm(fes)
    f += J * v * dx("conductor")
    with TaskManager():
        f.Assemble()
    gfA = GridFunction(fes)
    with TaskManager():
        gfA.vec.data = inv * f.vec

    R_inv = float(Integrate(J * J * sigma_inv_cf * dx("conductor"), mesh))
    R0 = 1.0 / R_inv
    L0 = R0 * R0 * float(Integrate(J * gfA * dx("conductor"), mesh))
    tau0 = L0 / R0 * 1e6
    print(f"  R_0 = {R0:.4e}, L_0 = {L0:.4e}, tau_0 = {tau0:.3f} us")
    return {"AIR_SCALE": air_scale, "ne": mesh.ne, "ndof": fes.ndof,
            "R0": R0, "L0": L0, "tau0_us": tau0,
            "ratio_to_ELF": tau0 / 11.51}


def main():
    ngsglobals.msg_level = 0
    print("=== Air-box scale dependence of Kameari stage-0 tau ===")
    print(f"  Conductor: 5x2x1 mm Cu")
    print(f"  ELF Foster reference tau_lead = 11.51 us\n")

    results = []
    # h_air scales with air box (keep similar element count)
    for air_scale, h_air in [(3, 1.0e-3), (5, 1.5e-3), (8, 2.5e-3), (12, 4.0e-3)]:
        r = stage0_tau(air_scale, h_air)
        results.append(r)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"{'AIR_SCALE':>10} {'ne':>6} {'tau_0 [us]':>12} {'tau_0/tau_ELF':>14}")
    for r in results:
        print(f"{r['AIR_SCALE']:>10} {r['ne']:>6} {r['tau0_us']:>12.3f} {r['ratio_to_ELF']:>14.4f}")

    out = {
        "method": "Air-box scale sweep, Kameari stage-0",
        "elf_foster_tau_lead_us": 11.51,
        "results": results,
    }
    out_path = Path(__file__).parent / "cuboid_521_vacuum_kameari_airscale_sweep.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
