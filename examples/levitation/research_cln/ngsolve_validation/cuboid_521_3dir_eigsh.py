"""5x2x1 mm Cu cuboid - 3-direction Foster via scipy eigsh shift-invert.

PINVIT failed to converge to lowest physical modes for HCurl with full Dirichlet.
Switch to scipy.sparse.linalg.eigsh with shift-invert at analytical target.

Strategy:
  1. Assemble S (curl-curl/mu0) and M (sigma * I) on conductor
  2. Use shift-invert at lambda_target ~ k^2/(mu0 sigma) for each direction
  3. Verify leading tau matches analytical TE mode

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
ax, ay, az = 5e-3, 2e-3, 1e-3
H_COND = 0.20e-3
ORDER = 2
N_PER_DIR = 12   # eigvecs near each direction's target

def tau_analytic(L1, L2):
    k2 = pi**2 * (1.0/L1**2 + 1.0/L2**2)
    return mu0 * sigma_Cu / k2

TAU_REF = {'x': tau_analytic(ay, az), 'y': tau_analytic(ax, az), 'z': tau_analytic(ax, ay)}
LAM_REF = {k: 1.0/v for k, v in TAU_REF.items()}


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
    """Convert NGSolve matrix to scipy CSR, restricted to active rows/cols."""
    rows, cols, vals = mat.COO()
    n_full = mat.height
    csr_full = sp.csr_matrix((vals, (rows, cols)), shape=(n_full, n_full))
    # Restrict to active rows/cols
    keep = np.array([free_dofs_mask[i] for i in range(n_full)], dtype=bool)
    return csr_full[keep][:, keep], keep


def main():
    ngsglobals.msg_level = 0
    print("=" * 72)
    print("Cuboid 5x2x1 Cu - 3-direction Foster via scipy shift-invert")
    print("=" * 72)
    print(f"  h={H_COND*1000} mm, order={ORDER}")
    print(f"  Analytical (interior PEC):")
    for k, tau in TAU_REF.items():
        print(f"    B_{k} TE leading tau = {tau*1e6:.4f} us  (lambda = {LAM_REF[k]:.2e})")

    box = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    box.mat("conductor").bc("conductor_surface")
    box.maxh = H_COND
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=H_COND))
    print(f"  ne={mesh.ne}, nv={mesh.nv}")

    fes = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)
    print(f"  HCurl ndof = {fes.ndof}")

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
    with TaskManager():
        a.Assemble()
        m.Assemble()

    Ssp, keep = to_csr(a.mat, fd)
    Msp, _ = to_csr(m.mat, fd)
    Ssp = (Ssp + Ssp.T) / 2.0
    Msp = (Msp + Msp.T) / 2.0
    print(f"  Active matrices: {Ssp.shape}, S nnz={Ssp.nnz}, M nnz={Msp.nnz}")

    # Per-direction shift-invert
    A_ext = {
        'x': CoefficientFunction((0, -z/2, y/2)),
        'y': CoefficientFunction((z/2, 0, -x/2)),
        'z': CoefficientFunction((-y/2, x/2, 0)),
    }

    foster_results = {}
    for kdir in ['z', 'y', 'x']:
        target = LAM_REF[kdir]
        print(f"\n=== B_{kdir}: shift-invert near lambda={target:.4e} ===")
        t0 = time.time()
        try:
            evals, evecs = spla.eigsh(
                Ssp, k=N_PER_DIR, M=Msp, sigma=target, which='LM',
                tol=1e-9, maxiter=3000
            )
        except Exception as e:
            print(f"  eigsh failed: {e}")
            continue
        order = np.argsort(evals)
        evals = evals[order]
        evecs = evecs[:, order]
        print(f"  eigsh done in {time.time()-t0:.1f}s")
        print(f"  eigenvalues: {evals}")

        modes = []
        for j in range(len(evals)):
            ev = float(evals[j])
            if ev < 1e-3:
                continue
            # Reconstruct full DOF vector and put into GridFunction
            full_vec = np.zeros(fes.ndof)
            full_vec[keep] = evecs[:, j]
            gf = GridFunction(fes)
            gf.vec.FV().NumPy()[:] = full_vec
            nrm2 = float(Integrate(sigma_Cu * gf * gf * dx, mesh))
            if nrm2 < 1e-30:
                continue
            gf.vec.FV().NumPy()[:] /= nrm2**0.5
            beta = {kk: float(Integrate(sigma_Cu * gf * A_ext[kk] * dx, mesh))
                    for kk in ['x', 'y', 'z']}
            b2 = {kk: beta[kk]**2 for kk in beta}
            tau_n = 1.0 / ev
            modes.append((tau_n, ev, b2))
            print(f"  rank {j}: lam={ev:.4e}, tau={tau_n*1e6:.4f} us, "
                  f"b2_x={b2['x']:.3e}, b2_y={b2['y']:.3e}, b2_z={b2['z']:.3e}")
        foster_results[kdir] = modes

    print(f"\n=== Summary: leading B_k coupled tau ===")
    for kdir in ['x', 'y', 'z']:
        modes = foster_results.get(kdir, [])
        # Find mode with largest beta^2 in this direction
        best = max(modes, key=lambda mm: mm[2][kdir]) if modes else None
        if best:
            tau_lead = best[0] * 1e6
            tau_ref = TAU_REF[kdir] * 1e6
            err = (tau_lead - tau_ref) / tau_ref * 100
            print(f"  B_{kdir}: leading tau = {tau_lead:.4f} us  "
                  f"(analytical {tau_ref:.4f}, err={err:+.3f}%, beta^2={best[2][kdir]:.3e})")
        else:
            print(f"  B_{kdir}: no modes found")


if __name__ == "__main__":
    main()
