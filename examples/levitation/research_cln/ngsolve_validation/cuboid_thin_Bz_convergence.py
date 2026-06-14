"""z-thin cuboid Bz convergence study: 5x2xc Cu, c sweep, B_z only.

Goal: verify FEM Foster pipeline against the rigorous 2D limit.

In z-thin limit (c << a,b), the eddy current induced by uniform B_z
is purely in xy plane (J_z = 0, no z-dependence). This reduces to the
2D rectangular plate problem on a x b rectangle:

   leading tau = mu0 * sigma / (pi^2 * (1/a^2 + 1/b^2))

For 5x2 mm:  tau_TE_z(1,1) = 25.4648 us  (independent of c)

By making c small, all TM_z (J_z type) modes get pushed to high freq
(~ k_z^2 contribution dominates), leaving the vortex modes (which
couple to B_z) cleanly isolated near TE_z(m,n,0) family.

Sweep c = {1.0, 0.5, 0.2, 0.1} mm. Report leading tau (=TE_z(1,1,0))
and verify it stays at 25.46 us (no c-dependence in true 2D limit).

Also provides axis-permutation argument:
  B_x of (a,b,c) cuboid = B_z of (b,c,a) cuboid
  B_y of (a,b,c) cuboid = B_z of (c,a,b) cuboid

So one z-aligned validation = all-direction validation by relabeling.

Date: 2026-05-03
"""
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z, BND,
    Integrate, TaskManager, ngsglobals,
)
from collections import deque
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import numpy as np
from math import pi
import time

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay = 5e-3, 2e-3
ORDER = 2

# Analytical TE_z(m,n,0) leading tau (independent of c)
def tau_te_z(L1, L2, m=1, n=1):
    k2 = pi**2 * ((m / L1) ** 2 + (n / L2) ** 2)
    return mu0 * sigma_Cu / k2

TAU_REF = tau_te_z(ax, ay) * 1e6  # us
print(f"Analytical: TE_z(1,1,0) leading tau = {TAU_REF:.6f} us  (independent of c)")


def classify_vertices(mesh):
    boundary_v = set()
    for el in mesh.Elements(BND):
        for v in el.vertices:
            boundary_v.add(v.nr)
    return set(range(mesh.nv)) - boundary_v, boundary_v


def build_spanning_tree_interior(mesh):
    interior_v, _ = classify_vertices(mesh)
    nv = mesh.nv
    visited = [v not in interior_v for v in range(nv)]
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


def to_csr(mat, free_dofs_mask):
    rows, cols, vals = mat.COO()
    n = mat.height
    csr = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    keep = np.array([free_dofs_mask[i] for i in range(n)], dtype=bool)
    return csr[keep][:, keep], keep


def run_one(az, h_xy, h_z, n_modes=10):
    """Solve 5x2xaz cuboid with B_z external. Return [(tau_us, beta2)] sorted."""
    box = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    box.mat("conductor").bc("conductor_surface")
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=min(h_xy, h_z)))
    fes = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)

    tree_edges = build_spanning_tree_interior(mesh)
    fd = fes.FreeDofs()
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False

    u, v = fes.TnT()
    a_form = BilinearForm(fes); a_form += (1.0/mu0) * curl(u) * curl(v) * dx
    m_form = BilinearForm(fes); m_form += sigma_Cu * u * v * dx
    with TaskManager():
        a_form.Assemble()
        m_form.Assemble()
    Ssp, keep = to_csr(a_form.mat, fd)
    Msp, _ = to_csr(m_form.mat, fd)
    Ssp = (Ssp + Ssp.T) / 2.0
    Msp = (Msp + Msp.T) / 2.0

    A_ext = CoefficientFunction((-y/2, x/2, 0))

    # Multi-target shift-invert: TE_z(odd,odd,0) modes for B_z
    targets_te = []
    for m in range(1, 12, 2):
        for n in range(1, 12, 2):
            targets_te.append((tau_te_z(ax, ay, m, n), m, n))
    targets_te.sort(reverse=True)
    targets_te = targets_te[:n_modes]

    foster = []
    for tau, m, n in targets_te:
        target_lam = 1.0 / tau
        try:
            evals_t, evecs_t = spla.eigsh(
                Ssp, k=8, M=Msp, sigma=target_lam, which='LM',
                tol=1e-9, maxiter=3000
            )
        except Exception as e:
            print(f"   ({m},{n}) FAILED: {e}")
            continue
        # Pick max-β² candidate
        best = None
        for j in range(len(evals_t)):
            ev = float(evals_t[j])
            if ev < 1e-3:
                continue
            full_vec = np.zeros(fes.ndof)
            full_vec[keep] = evecs_t[:, j]
            gf = GridFunction(fes)
            gf.vec.FV().NumPy()[:] = full_vec
            nrm2 = float(Integrate(sigma_Cu * gf * gf * dx, mesh))
            if nrm2 < 1e-30:
                continue
            gf.vec.FV().NumPy()[:] /= nrm2**0.5
            beta = float(Integrate(sigma_Cu * gf * A_ext * dx, mesh))
            b2 = beta * beta
            if best is None or b2 > best[1]:
                best = (1.0/ev * 1e6, b2, m, n, tau*1e6)
        if best:
            foster.append(best)

    # Dedupe by tau (rel-tol 0.5%)
    foster.sort()
    unique = []
    for f in foster:
        if unique and abs(f[0] - unique[-1][0])/max(f[0], 1e-30) < 0.005:
            continue
        unique.append(f)
    unique.sort(key=lambda f: -f[0])
    return unique, fes.ndof, mesh.ne


def main():
    ngsglobals.msg_level = 0
    print("=" * 78)
    print(" z-thin 5x2xc Cu cuboid: B_z Foster convergence (c sweep)")
    print("=" * 78)
    print(f" Geometry xy: {ax*1000} x {ay*1000} mm  (fixed)")
    print(f" Analytical TE_z(1,1,0) tau = {TAU_REF:.4f} us (independent of c)")

    cases = [
        # az_mm, h_xy_mm, h_z_mm
        (1.0,    0.20, 0.20),    # bulk
        (0.5,    0.20, 0.10),    # half thickness
        (0.2,    0.15, 0.05),    # thin
        (0.1,    0.10, 0.025),   # very thin
    ]
    print()
    print(f" {'c (mm)':>8} {'h_xy':>6} {'h_z':>6} {'ndof':>8} {'ne':>6}  {'leading tau':>14} {'err vs TE(1,1,0)':>18}")
    print(" " + "-" * 76)
    results = []
    for az_mm, h_xy_mm, h_z_mm in cases:
        az = az_mm * 1e-3
        h_xy = h_xy_mm * 1e-3
        h_z = h_z_mm * 1e-3
        t0 = time.time()
        modes, ndof, ne = run_one(az, h_xy, h_z, n_modes=10)
        if not modes:
            print(f" {az_mm:>8.2f} {h_xy_mm:>6.2f} {h_z_mm:>6.2f}  NO MODES")
            continue
        lead_tau, lead_b2, m, n, te_tau = modes[0]
        err = (lead_tau - TAU_REF) / TAU_REF * 100
        print(f" {az_mm:>8.2f} {h_xy_mm:>6.2f} {h_z_mm:>6.2f} {ndof:>8d} {ne:>6d}  {lead_tau:>14.4f} {err:>+15.3f}%  ({time.time()-t0:.1f}s)")
        results.append({
            "c_mm": az_mm,
            "ndof": ndof,
            "leading_tau_us": lead_tau,
            "err_pct": err,
            "n_modes": len(modes),
            "modes": modes,
        })

    print()
    print("=" * 78)
    print(" Top 5 modes per case (each row: c=...mm)")
    print("=" * 78)
    for r in results:
        print(f"\n c = {r['c_mm']} mm: ({r['n_modes']} modes)")
        for k, (tau, b2, m, n, te_tau) in enumerate(r['modes'][:6]):
            err = (tau - te_tau) / te_tau * 100
            print(f"   rank {k}: tau={tau:>10.4f} us, b2={b2:.3e},  TE_z({m},{n},0) ref={te_tau:.4f} us, err={err:+.2f}%")


if __name__ == "__main__":
    main()
