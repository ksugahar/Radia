"""3D cuboid in vacuum, frequency sweep of chi(s) for Stage 4 (beta') validation.

For each real s value, solve the eddy current problem:
  (1/mu_0) curl A_ind . curl v + s sigma A_ind . v = - s sigma A_ext . v   (in conductor)
  (1/mu_0) curl A_ind . curl v                       = 0                   (in air)
  A_ind tangential = 0 on outer air box                                    (PEC truncation)

Then m_z = (1/2) integral_cond (r x J_total)_z dV with J_total = -s sigma (A_ext + A_ind)
     chi_V(s) = m_z mu_0 / (V B_0)   (using B_0 = 1, H_0 = 1/mu_0)

Reference: Mathematica Stage 4 (beta') sphere benchmark with R_eff = R_iz.
Output: chi_cub_NGSolve(s) at s in {0, 1, 10, 1e2, ..., 1e8}, plus eta_2(s) = ratio.
"""
from netgen.occ import Box, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z,
    Integrate, TaskManager, ngsglobals,
)
from collections import deque
from math import pi
import json
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

# Conductor (centered at origin)
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

# Air box scale
AIR_SCALE = 5
H_COND = 0.30e-3
H_AIR = 1.5e-3
ORDER = 2

# s-values to sweep (matching Mathematica reference table)
S_VALUES = [1.0, 10.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]

# Reference (closed form)
ALPHA_LF_CUB = -mu0 * sigma_Cu * (ax**2 + ay**2) / 48.0
R_EFF = (5 * (ax**2 + ay**2) / 24)**0.5


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


def build_geometry():
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = H_COND
    Bx, By, Bz = AIR_SCALE*ax/2, AIR_SCALE*ay/2, AIR_SCALE*az/2
    air = Box(Pnt(-Bx, -By, -Bz), Pnt(Bx, By, Bz))
    air.mat("air").bc("outer_box")
    air.maxh = H_AIR
    return OCCGeometry(Glue([air - cuboid, cuboid]))


def chi_at_s(mesh, fes, free_dofs, s_val, sigma_cf, A_ext, sigma_inv_cf):
    """Solve eddy current at real s, return chi_V(s)."""
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (1.0/mu0) * curl(u) * curl(v) * dx
    a += s_val * sigma_cf * u * v * dx

    f = LinearForm(fes)
    f += -s_val * sigma_cf * A_ext * v * dx("conductor")

    with TaskManager():
        a.Assemble()
        f.Assemble()
        try:
            inv = a.mat.Inverse(free_dofs, inverse="sparsecholesky")
        except Exception as e:
            print(f"  s={s_val:.2e}: solver failed: {e}")
            return None

    gf = GridFunction(fes)
    with TaskManager():
        gf.vec.data = inv * f.vec

    A_total_x = A_ext[0] + gf[0]
    A_total_y = A_ext[1] + gf[1]
    # m_z = (1/2) integral_cond (x J_total_y - y J_total_x) dV
    #     = (1/2) (-s sigma) integral_cond (x A_total_y - y A_total_x) dV
    integrand = (x * A_total_y - y * A_total_x)
    integral = Integrate(integrand * dx("conductor"), mesh)
    m_z = 0.5 * (-s_val * sigma_Cu) * float(integral)
    # B_0 = 1, H_0 = B_0/mu_0 = 1/mu_0
    # chi_V = m_z / (V_cond H_0) = m_z mu_0 / V_cond
    chi = m_z * mu0 / V_cond
    return chi


def main():
    ngsglobals.msg_level = 0
    print("3D cuboid in vacuum, chi(s) frequency sweep (Stage 4 beta')")
    print(f"  Conductor: {ax*1000}x{ay*1000}x{az*1000} mm Cu")
    print(f"  Air box: {AIR_SCALE}x conductor")
    print(f"  R_eff (alpha_LF match) = {R_EFF*1000:.6f} mm")
    print(f"  alpha_LF^cub = {ALPHA_LF_CUB:.6e}")
    print()

    geo = build_geometry()
    print("Generating mesh...")
    mesh = Mesh(geo.GenerateMesh(maxh=H_AIR))
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
    print(f"  Tree-cotree masked {masked} Whitney edges, active DoFs = {sum(fd)}")

    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0})
    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0})
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    # Sphere reference (Mittag-Leffler fast form)
    def chi_sphere(R, s, n_max=80):
        if s == 0:
            return 0.0
        total = 0.0
        msr2 = mu0 * sigma_Cu * R * R
        for n in range(1, n_max + 1):
            tau_n = msr2 / (n * pi) ** 2
            total += 1.0 / ((n * pi) ** 4 * (1 + s * tau_n))
        return -9.0 * s * msr2 * total

    print()
    print("Sweeping s values...")
    print(f"{'s':>12} | {'chi_NGSolve':>20} | {'chi_sphere(R_eff)':>20} | {'eta_2 = ratio':>14}")
    print("-" * 78)
    results = []
    for s_val in S_VALUES:
        chi_ng = chi_at_s(mesh, fes, fd, s_val, sigma_cf, A_ext, sigma_inv_cf)
        chi_sph = chi_sphere(R_EFF, s_val)
        eta = chi_ng / chi_sph if chi_sph != 0 else float('nan')
        results.append({"s": s_val, "chi_ngsolve": chi_ng, "chi_sphere": chi_sph, "eta_2": eta})
        print(f"{s_val:12.3e} | {chi_ng:20.10e} | {chi_sph:20.10e} | {eta:14.6f}")

    # Low-s consistency: chi/s -> alpha_LF^cub
    print()
    print("Low-frequency consistency: chi(s)/s should -> alpha_LF^cub = "
          f"{ALPHA_LF_CUB:.6e}")
    for r in results[:4]:
        s_val = r["s"]
        ratio = r["chi_ngsolve"] / s_val
        rel_err = abs(ratio - ALPHA_LF_CUB) / abs(ALPHA_LF_CUB)
        print(f"  s = {s_val:.0e}:  chi/s = {ratio:.6e}  (rel err vs alpha_LF: {rel_err:.3e})")

    # High-s screening: chi -> shape-dependent constant (sphere -> -3/2)
    print()
    print("High-frequency screening (sphere limit -3/2):")
    for r in results[-3:]:
        s_val = r["s"]
        print(f"  s = {s_val:.0e}:  chi_cub = {r['chi_ngsolve']:.6f}, "
              f"chi_sph = {r['chi_sphere']:.6f}, eta = {r['eta_2']:.6f}")

    out = {
        "method": "Cuboid vacuum chi(s) frequency sweep (Stage 4 beta')",
        "geometry_mm": [ax*1000, ay*1000, az*1000],
        "R_eff_mm": R_EFF * 1000,
        "alpha_LF_cub": ALPHA_LF_CUB,
        "air_scale": AIR_SCALE,
        "h_cond_mm": H_COND*1000,
        "h_air_mm": H_AIR*1000,
        "order": ORDER,
        "ne": mesh.ne,
        "results": results,
    }
    out_path = Path(__file__).parent / "cuboid_521_vacuum_chi_sweep_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
