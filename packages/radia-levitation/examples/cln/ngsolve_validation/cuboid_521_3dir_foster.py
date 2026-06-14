"""5x2x1 mm Cu cuboid - 3-direction independent Foster (B_x, B_y, B_z).

Full 3D HCurl with interior PEC (n x A = 0 on conductor surface).
For each eigenmode phi_n, compute 3 coupling coefficients:
  beta_n^(k) = <phi_n, A_ext^(k)> for k=x,y,z

External excitation A fields giving uniform B_k = 1:
  B_x:  A_ext = (0, -z/2, y/2)
  B_y:  A_ext = (z/2, 0, -x/2)
  B_z:  A_ext = (-y/2, x/2, 0)

Foster series per direction:
  Y_kk(s) = sigma * sum_n beta_n^(k)^2 * tau_n / (1 + s tau_n)
  R_n^(k) = 1/(sigma * beta_n^(k)^2),  L_n^(k) = R_n^(k) * tau_n

Verification:
  B_z leading tau: TE_z(1,1,0) = mu0 sigma / (pi^2 (1/a^2 + 1/b^2)) = 25.46 us
  B_y leading tau: TE_y(1,0,1) = mu0 sigma / (pi^2 (1/a^2 + 1/c^2)) =  7.10 us
  B_x leading tau: TE_x(0,1,1) = mu0 sigma / (pi^2 (1/b^2 + 1/c^2)) =  5.91 us

Date: 2026-05-03
"""
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, GridFunction, Preconditioner,
    CoefficientFunction, curl, dx, x, y, z, BND,
    Integrate, TaskManager, ngsglobals,
)
from ngsolve.solvers import PINVIT
from collections import deque
from math import pi
import time
import numpy as np


def classify_vertices(mesh):
    boundary_v = set()
    for el in mesh.Elements(BND):
        for v in el.vertices:
            boundary_v.add(v.nr)
    return set(range(mesh.nv)) - boundary_v, boundary_v


def build_spanning_tree_interior(mesh):
    interior_v, _ = classify_vertices(mesh)
    nv = mesh.nv
    visited = [False] * nv
    for v in range(nv):
        if v not in interior_v:
            visited[v] = True
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        if v0 in interior_v and v1 in interior_v:
            adj[v0].append((v1, ed.nr))
            adj[v1].append((v0, ed.nr))
    for start in interior_v:
        if visited[start]:
            continue
        visited[start] = True
        queue = deque([start])
        while queue:
            v = queue.popleft()
            for vn, edn in adj[v]:
                if not visited[vn]:
                    visited[vn] = True
                    tree_edges.append(edn)
                    queue.append(vn)
    return tree_edges

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
H_COND = 0.20e-3
ORDER = 1
N_EIGEN = 20
PINVIT_MAXIT = 1500

# Analytical reference (interior PEC TE leading per direction)
def tau_analytic(L1, L2):
    """Leading tau for TE mode k^2 = pi^2(1/L1^2 + 1/L2^2)."""
    k2 = pi**2 * (1.0/L1**2 + 1.0/L2**2)
    return mu0 * sigma_Cu / k2

TAU_REF = {
    'x': tau_analytic(ay, az),  # 2x1: 5.91 us
    'y': tau_analytic(ax, az),  # 5x1: 7.10 us
    'z': tau_analytic(ax, ay),  # 5x2: 25.46 us
}


def build_cuboid():
    box = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    box.mat("conductor").bc("conductor_surface")
    box.maxh = H_COND
    return OCCGeometry(box)


def main():
    ngsglobals.msg_level = 0
    print("=" * 72)
    print("Cuboid 5x2x1 Cu - 3-direction Foster (interior PEC, full 3D HCurl)")
    print("=" * 72)
    print(f"  h = {H_COND*1000} mm, order = {ORDER}, N_eigen = {N_EIGEN}")
    print(f"  Analytical refs:")
    for k, tau in TAU_REF.items():
        print(f"    B_{k} leading tau = {tau*1e6:.4f} us")

    mesh = Mesh(build_cuboid().GenerateMesh(maxh=H_COND))
    print(f"\n  ne = {mesh.ne}, nv = {mesh.nv}")

    fes = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)
    print(f"  HCurl ndof = {fes.ndof}")

    # Tree-cotree mask in-place to remove residual gradient kernel beyond nograds
    tree_edges = build_spanning_tree_interior(mesh)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  Tree-cotree masked {masked} DoFs, active = {sum(fd)}")

    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += (1.0/mu0) * curl(u) * curl(v) * dx
    m = BilinearForm(fes)
    m += sigma_Cu * u * v * dx

    print("  Assembling...")
    pre = Preconditioner(a, type="local")
    with TaskManager():
        a.Assemble()
        m.Assemble()

    print(f"  Running PINVIT for {N_EIGEN} smallest eigenpairs...")
    t0 = time.time()
    gf_modes = GridFunction(fes, multidim=N_EIGEN, name="modes")

    with TaskManager():
        evals = PINVIT(a.mat, m.mat, pre=pre.mat, num=N_EIGEN,
                       maxit=PINVIT_MAXIT, printrates=False)
    print(f"  PINVIT done in {time.time()-t0:.1f}s")

    if isinstance(evals, tuple):
        evals_vec = evals[0]
        evec_multi = evals[1]
    else:
        evals_vec = evals
        evec_multi = gf_modes.vecs

    if hasattr(evals_vec, "NumPy"):
        evals_arr = np.array(evals_vec.NumPy())
    elif hasattr(evals_vec, "FV"):
        evals_arr = np.array(list(evals_vec.FV()))
    else:
        s = str(evals_vec).strip()
        evals_arr = np.array([float(xv.strip()) for xv in s.split("\n") if xv.strip()])

    # ev = lambda such that S v = lambda M v
    # S = (1/mu0)curl-curl, M = sigma*I  =>  tau_n = 1/lambda_n
    print(f"\n  Raw evals[:8]: {evals_arr[:8]}")

    A_ext = {
        'x': CoefficientFunction((0, -z/2, y/2)),
        'y': CoefficientFunction((z/2, 0, -x/2)),
        'z': CoefficientFunction((-y/2, x/2, 0)),
    }

    print(f"\n=== Per-mode coupling (beta^(k)^2) ===")
    print(f"{'mode':>4} {'tau (us)':>12}  {'beta^2_x':>12}  {'beta^2_y':>12}  {'beta^2_z':>12}  dominant")

    foster = {'x': [], 'y': [], 'z': []}   # list of (tau, beta2, R, L)
    EV_MIN = 1e-3
    for k in range(N_EIGEN):
        ev = float(evals_arr[k])
        if ev < EV_MIN:
            continue
        gf = GridFunction(fes)
        gf.vec.data = evec_multi[k]
        # M-normalize
        nrm2 = float(Integrate(sigma_Cu * gf * gf * dx, mesh))
        if nrm2 < 1e-30:
            continue
        gf.vec.data = (1.0/nrm2**0.5) * gf.vec
        tau_n = 1.0 / ev
        betas = {}
        for kdir in ['x', 'y', 'z']:
            b = float(Integrate(sigma_Cu * gf * A_ext[kdir] * dx, mesh))
            betas[kdir] = b
        b2 = {kdir: betas[kdir]**2 for kdir in betas}
        b2_sum = sum(b2.values())
        if b2_sum < 1e-30:
            continue
        dominant = max(b2, key=lambda kk: b2[kk])
        print(f"{k:>4} {tau_n*1e6:>12.4f}  {b2['x']:>12.4e}  {b2['y']:>12.4e}  {b2['z']:>12.4e}  B_{dominant}")
        for kdir in ['x', 'y', 'z']:
            if b2[kdir] / b2_sum > 1e-6:   # significant coupling
                Rn = 1.0 / (sigma_Cu * b2[kdir])
                Ln = Rn * tau_n
                foster[kdir].append((tau_n, b2[kdir], Rn, Ln))

    print(f"\n=== Foster decomposition per direction (sorted by tau desc) ===")
    for kdir in ['x', 'y', 'z']:
        modes = sorted(foster[kdir], key=lambda t: -t[0])
        if not modes:
            print(f"\n[B_{kdir}] no coupled modes found")
            continue
        tau_lead = modes[0][0] * 1e6
        tau_ref = TAU_REF[kdir] * 1e6
        err_pct = (tau_lead - tau_ref) / tau_ref * 100
        print(f"\n[B_{kdir}] {len(modes)} coupled modes, leading tau = {tau_lead:.4f} us  "
              f"(analytical {tau_ref:.4f} us, err = {err_pct:+.2f}%)")
        for j, (tau, b2, Rn, Ln) in enumerate(modes[:6]):
            print(f"  rank {j}: tau = {tau*1e6:>10.4f} us, beta^2 = {b2:.4e}, "
                  f"R = {Rn:.4e}, L = {Ln*1e9:.4e} nH")

    # DC sum check: sum tau_n beta_n^2 = chi(0) for each direction
    print(f"\n=== DC consistency: chi_kk(0) = sum_n tau_n beta_n^2 ===")
    for kdir in ['x', 'y', 'z']:
        modes = foster[kdir]
        chi0 = sum(tau * b2 for tau, b2, _, _ in modes)
        print(f"  chi_{kdir}{kdir}(0) = {chi0:.6e}")


if __name__ == "__main__":
    main()
