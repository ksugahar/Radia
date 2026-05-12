"""5x2x1 mm Cu cuboid - B_z mesh convergence study (verification finalization).

By axis symmetry, B_x and B_y are obtained from this by relabeling
(a,b,c) -> (b,c,a) etc., so B_z verification suffices.

For each mesh size h, run shift-invert at TE_z(1,1,0) target and confirm
the FEM leading tau converges to the analytical value 25.4648 us as h -> 0.

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
ORDER = 2

# Analytical TE_z(1,1,0): tau = mu0*sigma / [pi^2 (1/a^2 + 1/b^2)]
def tau_te_z(L1, L2, m=1, n=1):
    k2 = pi**2 * ((m / L1) ** 2 + (n / L2) ** 2)
    return mu0 * sigma_Cu / k2

TAU_REF = tau_te_z(ax, ay) * 1e6


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


def run_one(h_mm, n_modes_top=6):
    h = h_mm * 1e-3
    box = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    box.mat("conductor").bc("conductor_surface")
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=h))
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

    # Top n_modes_top TE_z(odd,odd,0) targets sorted by tau desc
    targets = []
    for m in range(1, 12, 2):
        for n in range(1, 12, 2):
            targets.append((tau_te_z(ax, ay, m, n), m, n))
    targets.sort(reverse=True)
    targets = targets[:n_modes_top]

    foster = []
    for tau, m, n in targets:
        target_lam = 1.0 / tau
        try:
            evals_t, evecs_t = spla.eigsh(
                Ssp, k=8, M=Msp, sigma=target_lam, which='LM',
                tol=1e-9, maxiter=3000
            )
        except Exception:
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
    print(" 5x2x1 mm Cu cuboid: B_z FEM convergence vs analytical TE_z(1,1,0)")
    print("=" * 78)
    print(f" Analytical: tau_TE_z(1,1,0) = {TAU_REF:.4f} us")
    print()
    print(f" {'h (mm)':>7} {'ndof':>7} {'ne':>6}  {'leading tau':>14} {'err':>10}  {'time':>6}")
    print(" " + "-" * 60)

    cases = [0.40, 0.30, 0.25, 0.20, 0.15, 0.12]
    results = []
    for h in cases:
        t0 = time.time()
        modes, ndof, ne = run_one(h, n_modes_top=6)
        if not modes:
            print(f" {h:>7.3f}  failed")
            continue
        lead_tau = modes[0][0]
        err = (lead_tau - TAU_REF) / TAU_REF * 100
        elapsed = time.time() - t0
        print(f" {h:>7.3f} {ndof:>7d} {ne:>6d}  {lead_tau:>14.4f} {err:>+9.3f}%  {elapsed:>5.1f}s")
        results.append({"h": h, "ndof": ndof, "lead_tau": lead_tau, "err": err, "modes": modes})

    print()
    print(" Top 6 TE_z(odd,odd,0) modes per h (mesh convergence per mode):")
    print(" {:>7} | ".format("h(mm)") + " | ".join([f"({m},{n})" for _, _, m, n, _ in results[0]["modes"][:6]]))
    print(" " + "-" * 78)
    for r in results:
        line = f" {r['h']:>7.3f} | "
        line += " | ".join([f"{m[0]:>6.3f} ({(m[0]-m[4])/m[4]*100:+5.2f}%)" for m in r["modes"][:6]])
        print(line)

    print()
    print(" Convergence trend (Richardson-like): ratio of consecutive errors")
    if len(results) >= 2:
        for i in range(1, len(results)):
            r1, r2 = results[i-1], results[i]
            ratio_h = r1["h"] / r2["h"]
            ratio_e = abs(r1["err"]) / abs(r2["err"]) if r2["err"] != 0 else float("inf")
            order_observed = np.log(ratio_e) / np.log(ratio_h) if ratio_h > 1 else 0
            print(f"  h={r1['h']:.3f} -> {r2['h']:.3f}  err: {r1['err']:+.3f}% -> {r2['err']:+.3f}%  observed order: {order_observed:.2f}")


if __name__ == "__main__":
    main()
