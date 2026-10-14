"""Compare global K, M for a 2-triangle mesh sharing an edge.

Two triangles forming a quad, all conductor. Compare:
  Python: AxifemP2Triangle per-tri + global assembly into 9-DOF (4 vertex +
          5 mid-edge) array.
  C++: H1Henrotte + AxiHenrotteStiffnessBFI/SigmaMassBFI assembly.
Permute DOFs to match; print max abs diff. This isolates whether C++ correctly
handles the SHARED mid-edge DOF between adjacent triangles.
"""
from __future__ import annotations
import sys
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import scipy.sparse as sp
from math import pi

import netgen.meshing as ng_meshing
from netgen.meshing import (
    Mesh as NgMesh, EdgeDescriptor, Element1D, Element2D, FaceDescriptor,
    MeshPoint, Pnt,
)
from ngsolve import (
    Mesh, BilinearForm, CoefficientFunction, TaskManager, ngsglobals,
)
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)
from axifem_p2_triangle import AxifemP2Triangle

MU0 = 4 * pi * 1e-7
SIGMA = 5.8e7

# Mesh: 4 vertices forming a quad, split into 2 triangles via diagonal v0-v2.
# Vertex coords (r, z), all off-axis to keep Vandermonde well-conditioned.
V = np.array([
    [1e-3, 0.0],
    [5e-3, 0.0],
    [5e-3, 2e-3],
    [1e-3, 2e-3],
])

# Triangles: tri_a = [v0, v1, v2], tri_b = [v0, v2, v3].
TRI = [(0, 1, 2), (0, 2, 3)]


def build_cpp_mesh():
    m = NgMesh()
    m.dim = 2
    m.SetMaterial(1, "conductor")
    m.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    m.SetBCName(0, "outer")
    edge = EdgeDescriptor()
    edge.edgenr = 1
    edge.surfnr = (1, -1)
    edge.domin = 1
    edge.domout = 0
    edge.name = "outer"
    m.Add(edge)
    pids = [m.Add(MeshPoint(Pnt(V[i, 0], V[i, 1], 0))) for i in range(4)]
    for t in TRI:
        m.Add(Element2D(1, [pids[t[0]], pids[t[1]], pids[t[2]]]))
    # Boundary edges (perimeter of quad).
    m.Add(Element1D([pids[0], pids[1]], index=1))
    m.Add(Element1D([pids[1], pids[2]], index=1))
    m.Add(Element1D([pids[2], pids[3]], index=1))
    m.Add(Element1D([pids[3], pids[0]], index=1))
    return Mesh(m), pids


def py_assemble():
    """Build P2 nodes list with edge dedup, assemble into global K, M."""
    nodes_p1 = V.tolist()
    nodes_p2 = list(nodes_p1)
    edge_to_mid = {}
    triangles_p2 = []
    for tri in TRI:
        v0, v1, v2 = tri
        def get_mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in edge_to_mid:
                mid = 0.5 * (V[a] + V[b])
                edge_to_mid[key] = len(nodes_p2)
                nodes_p2.append(mid.tolist())
            return edge_to_mid[key]
        m01 = get_mid(v0, v1)
        m12 = get_mid(v1, v2)
        m20 = get_mid(v2, v0)
        triangles_p2.append([v0, v1, v2, m01, m12, m20])
    nodes_p2 = np.array(nodes_p2)
    triangles_p2 = np.array(triangles_p2)
    N = nodes_p2.shape[0]
    K = np.zeros((N, N))
    M = np.zeros((N, N))
    for tri in triangles_p2:
        rn = nodes_p2[tri, 0]
        zn = nodes_p2[tri, 1]
        t = AxifemP2Triangle(rn, zn)
        Ke = t.stiffness(mu_r=MU0)
        Me = t.sigma_mass(sigma=SIGMA)
        for il in range(6):
            for jl in range(6):
                K[tri[il], tri[jl]] += Ke[il, jl]
                M[tri[il], tri[jl]] += Me[il, jl]
    return K, M, nodes_p2, edge_to_mid


def cpp_assemble():
    ngsglobals.msg_level = 0
    mesh, pids = build_cpp_mesh()
    fes = H1Henrotte(mesh, order=2, dirichlet=None)
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += AxiHenrotteStiffnessBFI(CoefficientFunction(MU0))
    with TaskManager(): a.Assemble()
    sigma_cf = mesh.MaterialCF({"conductor": SIGMA}, default=0.0)
    mbf = BilinearForm(fes, symmetric=True, check_unused=False)
    mbf += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager(): mbf.Assemble()
    def to_dense(bf, n):
        rs, cs, vs = bf.COO()
        return sp.coo_matrix((np.array(vs), (np.array(rs), np.array(cs))),
                             shape=(n, n)).toarray()
    return to_dense(a.mat, fes.ndof), to_dense(mbf.mat, fes.ndof), fes, mesh


def main():
    Kp, Mp, nodes_p2, e2m = py_assemble()
    print(f"Python: {len(nodes_p2)} nodes (4 vertex + {len(e2m)} mid-edge)")
    print(f"Python node coords:")
    for i, (r, z) in enumerate(nodes_p2):
        kind = "v" if i < 4 else "m"
        print(f"  [{i}] {kind} ({r:.4e}, {z:.4e})")

    Kc, Mc, fes, mesh = cpp_assemble()
    nv = mesh.nv
    ne = fes.ndof - nv - mesh.ne  # NEdges = ndof - nv - NE_faces
    print(f"\nC++: ndof={fes.ndof}  nv={nv}  ne={ne}  NE_tri={mesh.ne}")
    # Get C++ DOF -> physical coord mapping for verification.
    # Vertex DOFs 0..nv-1; edge DOFs nv..nv+ne-1; face DOFs at the end.
    print(f"C++ K (full):  max={np.max(np.abs(Kc)):.3e}  nnz_rows={np.sum(np.any(np.abs(Kc) > 1e-15, axis=1))}")
    print(f"C++ M (full):  max={np.max(np.abs(Mc)):.3e}")
    print(f"Python K (full): max={np.max(np.abs(Kp)):.3e}  nnz_rows={np.sum(np.any(np.abs(Kp) > 1e-15, axis=1))}")
    print(f"Python M (full): max={np.max(np.abs(Mp)):.3e}")

    # Find permutation by matching coords. Python orders [v0..v3, m01, m12, m20, m23, m30].
    # C++ orders vertex DOFs by netgen vertex order, then edge DOFs by edge index.
    # We need to identify which C++ DOF corresponds to which Python DOF.
    # Approach: solve K_c v = K_p v with constraint that nv vertex DOFs map first.
    # Since both have 4 vertex DOFs (in some order), iterate vertex permutations,
    # then compute edge mapping via mid-edge coords (we need NGSolve internal edge
    # midpoints to match).

    # Simpler: re-build Python's K, M as 9x9 (vertices + 5 mid-edges) and compare
    # via sums: if both assemblies are correct, trace(K), trace(M), sum(K),
    # sum(M) should agree exactly.
    Kc_useful = Kc[:nv + ne, :nv + ne]
    Mc_useful = Mc[:nv + ne, :nv + ne]
    print(f"\nC++ K (useful {nv + ne}x{nv + ne}): trace={np.trace(Kc_useful):.6e}  "
          f"sum={Kc_useful.sum():.6e}")
    print(f"Python K (9x9): trace={np.trace(Kp):.6e}  sum={Kp.sum():.6e}")
    print(f"C++ M (useful): trace={np.trace(Mc_useful):.6e}  sum={Mc_useful.sum():.6e}")
    print(f"Python M (9x9): trace={np.trace(Mp):.6e}  sum={Mp.sum():.6e}")

    # Diff via sum of eigenvalues
    Kc_sym = 0.5 * (Kc_useful + Kc_useful.T)
    Kp_sym = 0.5 * (Kp + Kp.T)
    ev_c = np.linalg.eigvalsh(Kc_sym)
    ev_p = np.linalg.eigvalsh(Kp_sym)
    print(f"\nC++ K eigvalsh: {ev_c}")
    print(f"Python K eigvalsh: {ev_p}")

    # Match DOFs by physical position. C++ DOF i corresponds to NGSolve's
    # i-th DOF. For vertex DOFs i<nv, it's the i-th mesh vertex. For edge
    # DOFs nv+e, it's the e-th mesh edge midpoint.
    # We need C++ -> Python mapping.
    # Vertex mapping: C++ vertex i is at the coordinate of mesh.vertices[i].
    cpp_dof_coord = np.zeros((nv + ne, 2))
    for vi, v in enumerate(mesh.vertices):
        cpp_dof_coord[vi] = (v.point[0], v.point[1])
    # Edge DOFs: midpoint of edge's two vertices.
    edges_list = []
    for ei, e in enumerate(mesh.edges):
        ev = list(e.vertices)
        v0 = mesh.vertices[ev[0].nr]
        v1 = mesh.vertices[ev[1].nr]
        mid = ((v0.point[0] + v1.point[0]) / 2, (v0.point[1] + v1.point[1]) / 2)
        cpp_dof_coord[nv + ei] = mid
        edges_list.append((ev[0].nr, ev[1].nr, mid))
    print(f"\nC++ DOFs (vertex + edge):")
    for k in range(nv + ne):
        print(f"  [{k}] ({cpp_dof_coord[k, 0]:.4e}, {cpp_dof_coord[k, 1]:.4e})")

    # Build permutation P such that nodes_p2[P[i]] == cpp_dof_coord[i].
    perm = []
    for i in range(nv + ne):
        target = cpp_dof_coord[i]
        diffs = np.linalg.norm(nodes_p2 - target, axis=1)
        idx = np.argmin(diffs)
        if diffs[idx] > 1e-12:
            raise RuntimeError(f"C++ DOF {i} at {target} has no Python match")
        perm.append(idx)
    perm = np.array(perm)
    print(f"\nPermutation C++ -> Python: {perm}")

    # Reorder Python K to match C++ ordering.
    Kp_perm = Kp[np.ix_(perm, perm)]
    Mp_perm = Mp[np.ix_(perm, perm)]
    print(f"\nKp(reordered) max abs - Kc max abs = "
          f"{np.max(np.abs(Kp_perm)):.3e} - {np.max(np.abs(Kc_useful)):.3e}")
    print(f"Diff |Kp_perm - Kc| max abs: {np.max(np.abs(Kp_perm - Kc_useful)):.3e}")
    print(f"Diff |Mp_perm - Mc| max abs: {np.max(np.abs(Mp_perm - Mc_useful)):.3e}")
    print(f"\nDiff matrix (Kc - Kp_perm) for k_largest mag entries:")
    diff = Kc_useful - Kp_perm
    abs_diff = np.abs(diff)
    flat_top = np.argsort(abs_diff.ravel())[::-1][:8]
    for f in flat_top:
        i, j = f // diff.shape[1], f % diff.shape[1]
        print(f"  ({i}, {j}) Kc={Kc_useful[i, j]:+.4e}  Kp={Kp_perm[i, j]:+.4e}  "
              f"diff={diff[i, j]:+.4e}")


if __name__ == "__main__":
    main()
