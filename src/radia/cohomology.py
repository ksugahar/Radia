"""Pure-Python (gmsh-free) cohomology of a multiply-connected NGSolve/Netgen mesh.

NGSolve/Netgen do not ship a homology/cohomology solver, and the legacy ``cohomology_cut.py`` path leaned on
Gmsh's ``computeHomology`` plus a ``.msh`` -> ``.vol`` transfer.  This module computes the first cohomology of
a tetrahedral mesh -- the "loops"/"cuts" of a holed magnetic body -- WITHOUT Gmsh, straight from the mesh.

The discrete de Rham complex of the tet mesh (lowest-order Whitney forms):

    C0 (vertices) --d0=G--> C1 (edges) --d1=Curl--> C2 (faces)

is built with the CANONICAL sorted-global-vertex orientation (edge oriented lo->hi; triangle a<b<c maps to
+e(a,b) +e(b,c) -e(a,c)), which guarantees the complex property ``Curl . G = 0``.  The first cohomology

    H^1 = ker(d1) / im(d0)   (closed 1-cochains modulo exact),    dim = b1 = number of independent loops,

is realised by the HARMONIC 1-cochains = the kernel of the combinatorial Hodge Laplacian (identity Hodge
stars)

    L1 = G G^T + Curl^T Curl       (E x E),     nullity(L1) = b1,     ker(L1) = harmonic 1-cochains = H^1.

The dimension b1 is topological (metric-independent), so the identity-star Laplacian gives the correct b1 and
a valid closed-non-exact representative in each class -- exactly what a T-Omega scalar-potential "cut" or a
loop basis needs (the period/circulation is what matters, not the metric shape).

The harmonic cochains are realised as HCurl(order=0) GridFunctions (``cohomology_basis``) and normalised to
UNIT CIRCULATION around a homology-generator basis (so each h_k carries one ampere-turn of one loop), using a
spanning-tree gauge + cotree-period construction -- again no Gmsh.

Validated (see ``__main__`` self-test): solid b1=0, washer (genus 1) b1=1, two-hole plate b1=2; the washer
basis function is curl-free (||curl h||/||h|| ~ 1e-13), has unit circulation around the hole (~+1) and ~0 on
a contractible loop; the two-hole period matrix is a non-singular "one loop per hole" permutation.

Public API:
    chain_complex(mesh)                 -> (G, Curl, eidx, EI, V, E, F)   sparse d0, d1 + edge tables
    betti_numbers(mesh)                 -> (b0, b1)
    cohomology(mesh)                    -> (b1, harm, ctx, vals)          raw harmonic 1-cochains (E x b1)
    cohomology_basis(mesh, ...)         -> (basis, b1, fes, ctx, loops)   HCurl(order=0) unit-circulation h_k
    circulation(h, mesh, cx, cy, rho)   -> oint h.dl around a test circle (validation helper)

This is the order-0 (Whitney) realisation; it is the topology engine the T-Omega ``CohomologyCutSolver`` is
being migrated onto (dropping its Gmsh dependency).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla
import ngsolve as ng

__all__ = ["chain_complex", "betti_numbers", "cohomology", "cohomology_basis", "circulation"]


def chain_complex(mesh):
    """Build d0=G (E x V) and d1=Curl (F x E) of the tet mesh in the canonical sorted-vertex orientation.

    Edges (6 per tet) and faces (4 per tet) are enumerated FROM THE TETS and deduped by sorted-vertex key --
    version-robust (only needs ``el.vertices``), and the sorted-vertex orientation makes ``Curl @ G == 0`` by
    construction.  Returns ``(G, Curl, eidx, EI, V, E, F)`` where ``eidx[(lo,hi)] -> edge index`` and
    ``EI[edge index] -> (lo,hi)``.
    """
    V = mesh.nv
    eidx, EI = {}, []
    fidx, FI = {}, []
    for el in mesh.Elements(ng.VOL):
        vs = sorted(int(v.nr) for v in el.vertices)            # tet: 4 vertices
        for i in range(4):
            for j in range(i + 1, 4):
                key = (vs[i], vs[j])
                if key not in eidx:
                    eidx[key] = len(EI); EI.append(key)
        for k in range(4):                                     # 4 triangle faces (drop vertex k); pre-sorted
            tri = tuple(vs[m] for m in range(4) if m != k)
            if tri not in fidx:
                fidx[tri] = len(FI); FI.append(tri)
    E, F = len(EI), len(FI)
    r, c, d = [], [], []
    for ei, (lo, hi) in enumerate(EI):                         # d0 = G : (G f)[e] = f(hi) - f(lo)
        r += [ei, ei]; c += [hi, lo]; d += [1.0, -1.0]
    G = sp.csr_matrix((d, (r, c)), shape=(E, V))
    r, c, d = [], [], []
    for fi, (a, b, cc) in enumerate(FI):                       # d1 = Curl : boundary(a<b<c) = +ab +bc -ac
        for (p, q), s in (((a, b), 1.0), ((b, cc), 1.0), ((a, cc), -1.0)):
            r.append(fi); c.append(eidx[(p, q)]); d.append(s)
    Curl = sp.csr_matrix((d, (r, c)), shape=(F, E))
    return G, Curl, eidx, EI, V, E, F


def cohomology(mesh, maxgen=24, tol=1e-4):
    """First cohomology via the combinatorial Hodge-Laplacian nullspace.

    Returns ``(b1, harm, ctx, vals)``: ``harm`` is ``(E x b1)`` -- columns are harmonic 1-cochains (edge
    values in the canonical orientation), ``ctx = (G, Curl, eidx, EI, V, E, F)``, ``vals`` the smallest
    eigenvalues (for inspecting the spectral gap).  The CALLER opens ``with ng.TaskManager():`` if desired.
    """
    G, Curl, eidx, EI, V, E, F = chain_complex(mesh)
    assert abs(Curl @ G).sum() < 1e-9, "Curl.G != 0 -- orientation bug in the chain complex"
    L1 = (G @ G.T + Curl.T @ Curl).tocsc()
    k = max(1, min(maxgen, E - 2))
    # shift-invert just BELOW the (PSD) spectrum: (L1 + 1e-6 I)^-1 is well-conditioned and surfaces the
    # smallest eigenvalues; the b1 harmonic ones come out ~0, well separated from the first nonzero.
    vals, vecs = spla.eigsh(L1, k=k, sigma=-1e-6, which="LM")
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    sel = vals < tol
    b1 = int(np.sum(sel))
    harm = vecs[:, sel] if b1 else np.zeros((E, 0))
    return b1, harm, (G, Curl, eidx, EI, V, E, F), vals


def betti_numbers(mesh):
    """Return ``(b0, b1)`` -- number of connected components and number of independent loops."""
    G, Curl, eidx, EI, V, E, F = chain_complex(mesh)
    parent = list(range(V))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a

    for (lo, hi) in EI:
        ra, rb = find(lo), find(hi)
        if ra != rb:
            parent[ra] = rb
    b0 = len({find(v) for v in range(V)})
    b1, _, _, _ = cohomology(mesh)
    return b0, b1


# ----------------------------------------------------------------------- unit-circulation normalisation
def _spanning_tree(V, EI):
    """BFS spanning forest of the 1-skeleton -> (order, parent_edge, cotree)."""
    from collections import deque
    adj = [[] for _ in range(V)]
    for ei, (lo, hi) in enumerate(EI):
        adj[lo].append((hi, ei)); adj[hi].append((lo, ei))
    parent_edge = [None] * V
    visited = [False] * V
    intree, order = set(), []
    for s in range(V):
        if visited[s]:
            continue
        visited[s] = True; dq = deque([s])
        while dq:
            u = dq.popleft(); order.append(u)
            for (w, ei) in adj[u]:
                if not visited[w]:
                    visited[w] = True; parent_edge[w] = ei; intree.add(ei); dq.append(w)
    cotree = [ei for ei in range(len(EI)) if ei not in intree]
    return order, parent_edge, cotree


def _periods(harm, V, EI):
    """Tree-gauge each closed cochain (subtract a gradient so it vanishes on tree edges); its values on the
    COTREE edges are then its periods (circulations around the fundamental loops).  Returns ``(Pi, cotree)``
    with ``Pi`` (ncotree x b1); ``rank(Pi)=b1`` and b1 independent cotree rows form a homology generator basis.
    """
    order, parent_edge, cotree = _spanning_tree(V, EI)
    b1 = harm.shape[1]
    phi = np.zeros((V, b1))
    for w in order:                                            # BFS order: parent set before child
        ei = parent_edge[w]
        if ei is None:
            continue
        lo, hi = EI[ei]                                        # want phi[hi]-phi[lo] = harm[ei] on tree edges
        phi[w] = phi[lo] + harm[ei] if w == hi else phi[hi] - harm[ei]
    Pi = np.array([harm[ei] - (phi[hi] - phi[lo]) for ei in cotree
                   for lo, hi in [EI[ei]]]).reshape(len(cotree), b1)
    return Pi, cotree


def cohomology_basis(mesh, unit_circulation=True):
    """Realise the harmonic 1-cochains as HCurl(order=0) GridFunctions ``h_k`` -- the cohomology basis the
    T-Omega scalar-potential cut needs.  They are curl-free (the sorted-vertex orientation matches NGSolve's
    HCurl order-0 lo->hi edge convention, so the discrete curl vanishes).  With ``unit_circulation`` the basis
    is normalised so the periods around a homology-generator basis are the IDENTITY (each ``h_k`` carries one
    ampere-turn of loop k).  Returns ``(basis, b1, fes, ctx, loops)``; ``loops`` are the generator cotree edge
    indices (or ``None`` if not normalised).  The CALLER opens ``with ng.TaskManager():``.
    """
    b1, harm, ctx, vals = cohomology(mesh)
    V, EI = ctx[4], ctx[3]
    loops = None
    if unit_circulation and b1 > 0:
        Pi, cotree = _periods(harm, V, EI)
        _, _, piv = sla.qr(Pi.T, pivoting=True, mode="economic")   # most-independent cotree rows first
        sel = piv[:b1]
        harm = harm @ np.linalg.inv(Pi[sel])                       # periods on the selected loops -> identity
        loops = [cotree[s] for s in sel]
    fes = ng.HCurl(mesh, order=0)
    ngdof = {}
    for e in mesh.edges:                                           # NGSolve edge -> its order-0 HCurl DOF
        vs = sorted(int(v.nr) for v in e.vertices)
        ngdof[(vs[0], vs[1])] = list(fes.GetDofNrs(e))[0]
    basis = []
    for k in range(b1):
        gf = ng.GridFunction(fes)
        vec = gf.vec.FV().NumPy(); vec[:] = 0.0
        for i, (lo, hi) in enumerate(EI):
            vec[ngdof[(lo, hi)]] = harm[i, k]
        basis.append(gf)
    return basis, b1, fes, ctx, loops


def circulation(h, mesh, cx, cy, rho, z=0.0, M=240):
    """oint h.dl around a circle (centre (cx,cy,z), radius rho in the xy-plane) -- numeric line integral.
    Validation helper: ~1 around a hole that h_k threads, ~0 on a contractible loop."""
    s = 0.0
    for t in np.linspace(0, 2 * np.pi, M, endpoint=False):
        hv = h(mesh(cx + rho * np.cos(t), cy + rho * np.sin(t), z))
        s += (-hv[0] * np.sin(t) + hv[1] * np.cos(t)) * (2 * np.pi * rho / M)
    return s


if __name__ == "__main__":
    from netgen.occ import Box, Cylinder, Pnt, Dir, OCCGeometry

    ax = Dir(0, 0, 1)
    cyl = lambda z0, r, h: Cylinder(Pnt(0, 0, z0), ax, r=r, h=h)
    cases = [("solid_cyl  (b1=0)", cyl(-0.015, 0.05, 0.03), 0),
             ("washer     (b1=1)", cyl(-0.015, 0.05, 0.03) - cyl(-0.025, 0.02, 0.05), 1)]
    twohole = (Box(Pnt(-0.08, -0.04, -0.015), Pnt(0.08, 0.04, 0.015))
               - cyl(-0.025, 0.015, 0.05).Move((-0.04, 0, 0)) - cyl(-0.025, 0.015, 0.05).Move((0.04, 0, 0)))
    cases.append(("two_hole   (b1=2)", twohole, 2))
    print(f"{'geom':22s} {'V':>5s} {'E':>6s} {'F':>6s} {'b1':>3s} {'exp':>4s}  smallest eigs")
    allok = True
    for name, shape, exp in cases:
        mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=0.02))
        b1, harm, ctx, vals = cohomology(mesh)   # caller-wraps policy: self-test runs serial
        V, E, F = ctx[4], ctx[5], ctx[6]
        allok &= (b1 == exp)
        print(f"{name:22s} {V:5d} {E:6d} {F:6d} {b1:3d} {exp:4d}  " + " ".join(f"{v:+.1e}" for v in vals[:4])
              + ("  OK" if b1 == exp else "  *** MISMATCH ***"))

    print("\n--- HCurl basis, UNIT-normalised (washer) ---")
    mesh = ng.Mesh(OCCGeometry(cyl(-0.015, 0.05, 0.03) - cyl(-0.025, 0.02, 0.05)).GenerateMesh(maxh=0.02))
    basis, b1, fes, ctx, loops = cohomology_basis(mesh)
    h = basis[0]
    curl_rel = (np.sqrt(ng.Integrate(ng.InnerProduct(ng.curl(h), ng.curl(h)), mesh))
                / np.sqrt(ng.Integrate(ng.InnerProduct(h, h), mesh)))
    c_hole = circulation(h, mesh, 0, 0, 0.035)
    c_contr = circulation(h, mesh, 0.035, 0, 0.006)
    print(f"   ||curl h||/||h||={curl_rel:.1e}  oint_hole={c_hole:+.4f} (unit~+-1)  oint_contractible={c_contr:+.1e}")
    ok1 = curl_rel < 1e-6 and abs(abs(c_hole) - 1) < 0.05 and abs(c_contr) < 1e-2

    print("--- two-hole period matrix (one loop per hole?) ---")
    mesh2 = ng.Mesh(OCCGeometry(twohole).GenerateMesh(maxh=0.02))
    basis2, b12, _, _, _ = cohomology_basis(mesh2)
    P = np.array([[circulation(basis2[k], mesh2, cx, 0.0, 0.026) for cx in (-0.04, 0.04)]
                  for k in range(b12)])
    print(f"   det(period matrix) = {np.linalg.det(P):+.4f} (non-singular separates the holes)")
    ok2 = abs(np.linalg.det(P)) > 0.1
    print(f"\n   {'SELF-TEST PASS' if (allok and ok1 and ok2) else 'SELF-TEST FAIL'} -- gmsh-free cohomology.")
