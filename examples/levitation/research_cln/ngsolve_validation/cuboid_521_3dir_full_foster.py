"""5x2x1 mm Cu cuboid - full Foster series per direction (25-30 modes each).

Multi-target shift-invert: for each direction, run scipy.eigsh at multiple
analytical TE-mode targets to populate the Foster spectrum. Classify each
mode by dominant beta^2 direction, dedupe across targets, sort by tau desc.

Output:
  cuboid_521_foster.json
    {
      "B_x": [{"rank":1,"tau_us":...,"lambda":...,"beta2":...,"R":...,"L":...}, ...],
      "B_y": [...],
      "B_z": [...],
    }

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
import json
from itertools import product

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
ax, ay, az = 5e-3, 2e-3, 1e-3
H_COND = 0.20e-3
ORDER = 2
N_PER_TARGET = 8     # eigvecs per shift-invert call
DEDUPE_REL_TOL = 1e-3   # initial dedupe: identical eigenvalues from same Krylov pass
CLUSTER_REL_TOL = 0.05  # final cluster: λ within 5% → keep max-β² (handles TE+TM degeneracy + FEM noise)
EV_MIN = 1e-3        # filter spurious near-zero
TAU_MIN_US = 0.3     # discard modes with tau < 0.3 us (noise)


def tau_te(L1, L2, m, n):
    """Analytical TE mode tau (k^2 = pi^2(m^2/L1^2 + n^2/L2^2))."""
    k2 = pi**2 * ((m / L1) ** 2 + (n / L2) ** 2)
    return mu0 * sigma_Cu / k2


def gen_targets(L1, L2, n_modes=10):
    """Generate analytical TE odd-odd targets for one direction, sorted by tau desc."""
    pairs = []
    for m, n in product(range(1, 12, 2), range(1, 12, 2)):  # odd 1,3,5,...,11
        tau = tau_te(L1, L2, m, n)
        pairs.append((tau, m, n))
    pairs.sort(reverse=True)
    return pairs[:n_modes]


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
    n_full = mat.height
    csr = sp.csr_matrix((vals, (rows, cols)), shape=(n_full, n_full))
    keep = np.array([free_dofs_mask[i] for i in range(n_full)], dtype=bool)
    return csr[keep][:, keep], keep


def main():
    ngsglobals.msg_level = 0
    print("=" * 78)
    print(" Cuboid 5x2x1 Cu - Full Foster series per direction (multi-target shift-invert) ")
    print("=" * 78)
    print(f" h={H_COND*1000} mm, order={ORDER}, N_per_target={N_PER_TARGET}")

    # Analytical targets per direction
    targets = {
        'x': gen_targets(ay, az),  # B_x: J in y-z plane, dimensions b=2, c=1
        'y': gen_targets(ax, az),  # B_y: J in x-z plane, dimensions a=5, c=1
        'z': gen_targets(ax, ay),  # B_z: J in x-y plane, dimensions a=5, b=2
    }
    print("\n Analytical TE targets (odd-odd) per direction:")
    for kdir in ('x', 'y', 'z'):
        line = f"  B_{kdir}: "
        for tau, m, n in targets[kdir][:8]:
            line += f"({m},{n})={tau*1e6:.2f}us  "
        print(line)

    # Build mesh + matrices
    print("\n Building mesh + assembling...")
    box = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    box.mat("conductor").bc("conductor_surface")
    box.maxh = H_COND
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=H_COND))
    fes = HCurl(mesh, order=ORDER, dirichlet="conductor_surface", nograds=True)
    print(f"  ne={mesh.ne}, nv={mesh.nv}, HCurl ndof={fes.ndof}")

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
    a_form = BilinearForm(fes)
    a_form += (1.0/mu0) * curl(u) * curl(v) * dx
    m_form = BilinearForm(fes)
    m_form += sigma_Cu * u * v * dx
    with TaskManager():
        a_form.Assemble()
        m_form.Assemble()

    Ssp, keep = to_csr(a_form.mat, fd)
    Msp, _ = to_csr(m_form.mat, fd)
    Ssp = (Ssp + Ssp.T) / 2.0
    Msp = (Msp + Msp.T) / 2.0
    print(f"  Active matrices: {Ssp.shape}, S nnz={Ssp.nnz}")

    A_ext = {
        'x': CoefficientFunction((0, -z/2, y/2)),
        'y': CoefficientFunction((z/2, 0, -x/2)),
        'z': CoefficientFunction((-y/2, x/2, 0)),
    }

    # Per-direction: for each analytical TE target, run shift-invert and keep
    # the FEM mode with MAX β² for this direction (the physical TE/TM-vortex pair).
    # → guarantees one Foster mode per analytical TE target. Robust against
    # TE/TM degeneracy + FEM discretization spread (no clustering needed).
    foster = {}
    for kdir in ('x', 'y', 'z'):
        print(f"\n=== B_{kdir} direction ===")
        modes = []
        t_total = time.time()
        for tau, m, n in targets[kdir]:
            target = 1.0 / tau
            t0 = time.time()
            try:
                evals_t, evecs_t = spla.eigsh(
                    Ssp, k=N_PER_TARGET, M=Msp, sigma=target, which='LM',
                    tol=1e-9, maxiter=3000
                )
            except Exception as e:
                print(f"  ({m},{n}) target tau={tau*1e6:.3f}us FAILED: {e}")
                continue
            # Among 8 candidates, find the one with MAX β² for this direction
            best = None
            for j in range(len(evals_t)):
                ev = float(evals_t[j])
                if ev < EV_MIN:
                    continue
                full_vec = np.zeros(fes.ndof)
                full_vec[keep] = evecs_t[:, j]
                gf = GridFunction(fes)
                gf.vec.FV().NumPy()[:] = full_vec
                nrm2 = float(Integrate(sigma_Cu * gf * gf * dx, mesh))
                if nrm2 < 1e-30:
                    continue
                gf.vec.FV().NumPy()[:] /= nrm2**0.5
                beta = float(Integrate(sigma_Cu * gf * A_ext[kdir] * dx, mesh))
                b2 = beta * beta
                if best is None or b2 > best["beta2"]:
                    best = {"lambda": ev, "tau_us": 1.0/ev * 1e6, "beta2": b2,
                            "te_target": (m, n), "te_tau_us": tau * 1e6}
            if best is None:
                print(f"  ({m},{n}) target tau={tau*1e6:.3f}us: no valid mode")
                continue
            err = (best["tau_us"] - best["te_tau_us"]) / best["te_tau_us"] * 100
            print(f"  ({m},{n}) target tau={tau*1e6:.3f}us -> FEM tau={best['tau_us']:.3f}us "
                  f"(err {err:+.2f}%), b2={best['beta2']:.3e}  ({time.time()-t0:.1f}s)")
            modes.append(best)

        # Dedupe (same FEM eigenvalue picked from multiple analytical targets)
        modes.sort(key=lambda mm: mm["lambda"])
        unique = []
        for mm in modes:
            if unique and abs(mm["lambda"] - unique[-1]["lambda"]) / max(mm["lambda"], 1e-30) < DEDUPE_REL_TOL:
                # Same FEM mode; keep the one with closer-to-analytical-target match
                if abs(mm["te_tau_us"] - mm["tau_us"]) < abs(unique[-1]["te_tau_us"] - unique[-1]["tau_us"]):
                    unique[-1] = mm
                continue
            unique.append(mm)
        if len(unique) < len(modes):
            print(f"  Dedupe: {len(modes)} → {len(unique)} (same FEM mode selected by multiple targets)")
        modes = unique

        # Sort by tau desc, assign rank, compute R/L
        modes.sort(key=lambda d: -d["tau_us"])
        for j, mm in enumerate(modes):
            mm["rank"] = j
            mm["R"] = 1.0 / (sigma_Cu * mm["beta2"])
            mm["L"] = mm["R"] * mm["tau_us"] * 1e-6
        foster[kdir] = modes
        print(f"  Final: {len(modes)} Foster modes ({time.time()-t_total:.1f}s)")

    # Print per-direction summary
    print("\n" + "=" * 78)
    print(" FOSTER SERIES per direction (sorted by tau desc, classified by dominant beta^2)")
    print("=" * 78)
    for kdir in ('x', 'y', 'z'):
        n = len(foster[kdir])
        if n == 0:
            print(f"\n B_{kdir}: NO MODES (empty)")
            continue
        tau_lead = foster[kdir][0]["tau_us"]
        ana_lead = targets[kdir][0][0] * 1e6
        err = (tau_lead - ana_lead) / ana_lead * 100
        print(f"\n B_{kdir}: {n} coupled modes  (leading tau={tau_lead:.4f} us, analytic={ana_lead:.4f}, err={err:+.2f}%)")
        print(f"  {'rank':>4} {'tau (us)':>12} {'beta^2':>14} {'R (Ohm)':>14} {'L (nH)':>14}")
        for mm in foster[kdir][:30]:
            print(f"  {mm['rank']:>4} {mm['tau_us']:>12.4f} {mm['beta2']:>14.4e} {mm['R']:>14.4e} {mm['L']*1e9:>14.4e}")
        if n > 30:
            print(f"  ... ({n - 30} more, see JSON output)")

    # DC moment chi(0) = sum tau * beta^2
    print("\n" + "=" * 78)
    print(" DC moment per direction: chi_kk(0) = sum_n tau_n * beta_n^2")
    print("=" * 78)
    for kdir in ('x', 'y', 'z'):
        chi0 = sum(mm["tau_us"] * 1e-6 * mm["beta2"] for mm in foster[kdir])
        print(f"  chi_{kdir}{kdir}(0) = {chi0:.6e}")

    # Save JSON
    out_path = "cuboid_521_foster.json"
    json_out = {
        f"B_{kdir}": [{
            "rank": mm["rank"],
            "tau_us": mm["tau_us"],
            "lambda": mm["lambda"],
            "beta2": mm["beta2"],
            "R": mm["R"],
            "L": mm["L"],
        } for mm in foster[kdir]] for kdir in ('x', 'y', 'z')
    }
    json_out["meta"] = {
        "geometry_mm": [ax*1000, ay*1000, az*1000],
        "sigma": sigma_Cu,
        "mu0": mu0,
        "mesh_h_mm": H_COND*1000,
        "fes_order": ORDER,
        "ndof": fes.ndof,
        "n_targets_per_dir": len(targets['x']),
    }
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(json_out, f, ensure_ascii=False, indent=2)
    print(f"\n Wrote {out_path} ({sum(len(foster[k]) for k in 'xyz')} modes total)")


if __name__ == "__main__":
    main()
