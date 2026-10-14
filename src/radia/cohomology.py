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
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla
import ngsolve as ng

__all__ = ["chain_complex", "betti_numbers", "cohomology", "cohomology_basis", "circulation",
           "surface_chain_complex", "surface_cohomology", "surface_fundamental_cycles",
           "surface_homology_loops"]


def _assemble_d0_d1(V, EI, FI, eidx):
    """Assemble d0=G (E x V) and d1=Curl (F x E) from the canonical
    sorted-vertex edge/face tables (shared by the tet and the triangle-
    surface complexes; ``Curl @ G == 0`` by the sorted orientation)."""
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
    return G, Curl


def _harmonic_cochains(ctx, maxgen, tol):
    """Nullspace of the combinatorial Hodge Laplacian L1 = G G^T + Curl^T Curl
    (shared by ``cohomology`` and ``surface_cohomology``).  Returns
    ``(b1, harm, vals)``."""
    G, Curl, eidx, EI, V, E, F = ctx
    assert abs(Curl @ G).sum() < 1e-9, "Curl.G != 0 -- orientation bug in the chain complex"
    L1 = (G @ G.T + Curl.T @ Curl).tocsc()
    k = max(1, min(maxgen, E - 2))
    # shift-invert just BELOW the (PSD) spectrum: (L1 + 1e-6 I)^-1 is well-conditioned and surfaces the
    # smallest eigenvalues; the b1 harmonic ones come out ~0, well separated from the first nonzero.
    # DETERMINISM CONTRACT: fix the Lanczos start vector -- with a random
    # v0 the (degenerate, b1-dimensional) nullspace basis rotates from
    # call to call, so downstream generator selection (QR pivoting on
    # the period matrix) would pick a different -- equivalent but not
    # identical -- cut each run.
    v0 = np.full(E, 1.0 / math.sqrt(E))
    vals, vecs = spla.eigsh(L1, k=k, sigma=-1e-6, which="LM", v0=v0)
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    sel = vals < tol
    b1 = int(np.sum(sel))
    harm = vecs[:, sel] if b1 else np.zeros((E, 0))
    return b1, harm, vals


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
    G, Curl = _assemble_d0_d1(V, EI, FI, eidx)
    return G, Curl, eidx, EI, V, len(EI), len(FI)


def surface_chain_complex(tris, nv=None):
    """Build d0=G and d1=Curl of a TRIANGULATED SURFACE (pure arrays, no NGSolve mesh).

    The triangle-surface counterpart of ``chain_complex``: C0 = vertices, C1 = edges, C2 = the triangles
    themselves, in the same canonical sorted-vertex orientation (so ``Curl @ G == 0`` by construction).
    The complex is combinatorial -- only connectivity enters, no coordinates.  For a CLOSED orientable
    surface, dim H^1 = b1 = 2 * genus.

    Args:
        tris: (nt, 3) vertex indices (winding irrelevant -- faces are canonicalised by sorting).
        nv: vertex-id space size (default ``tris.max() + 1``); pass the mesh vertex count when isolated
            vertices must keep their ids.

    Returns ``(G, Curl, eidx, EI, V, E, F)`` exactly like ``chain_complex``.
    """
    tris = np.asarray(tris, dtype=np.int64)
    V = int(nv) if nv is not None else int(tris.max()) + 1
    eidx, EI = {}, []
    fidx, FI = {}, []
    for t in tris:
        vs = sorted(int(v) for v in t)
        for i in range(3):
            for j in range(i + 1, 3):
                key = (vs[i], vs[j])
                if key not in eidx:
                    eidx[key] = len(EI); EI.append(key)
        tri = tuple(vs)
        if tri not in fidx:
            fidx[tri] = len(FI); FI.append(tri)
    G, Curl = _assemble_d0_d1(V, EI, FI, eidx)
    return G, Curl, eidx, EI, V, len(EI), len(FI)


def cohomology(mesh, maxgen=24, tol=1e-4):
    """First cohomology via the combinatorial Hodge-Laplacian nullspace.

    Returns ``(b1, harm, ctx, vals)``: ``harm`` is ``(E x b1)`` -- columns are harmonic 1-cochains (edge
    values in the canonical orientation), ``ctx = (G, Curl, eidx, EI, V, E, F)``, ``vals`` the smallest
    eigenvalues (for inspecting the spectral gap).  The CALLER opens ``with ng.TaskManager():`` if desired.
    """
    ctx = chain_complex(mesh)
    b1, harm, vals = _harmonic_cochains(ctx, maxgen, tol)
    return b1, harm, ctx, vals


def surface_cohomology(tris, nv=None, maxgen=24, tol=1e-4):
    """First cohomology of a triangulated surface (b1 = 2*genus for closed orientable) -- the surface
    counterpart of ``cohomology``.  Pure numpy/scipy: same combinatorial Hodge-Laplacian nullspace, on the
    ``surface_chain_complex``.  Returns ``(b1, harm, ctx, vals)``.
    """
    ctx = surface_chain_complex(tris, nv)
    b1, harm, vals = _harmonic_cochains(ctx, maxgen, tol)
    return b1, harm, ctx, vals


def surface_fundamental_cycles(tris, nv=None, maxgen=24, tol=1e-4):
    """Harmonic periods over ALL cotree fundamental cycles of a triangulated surface, plus an expander.

    Every cotree edge closes a FUNDAMENTAL CYCLE through the spanning tree -- always a simple vertex loop --
    and its homology class is encoded by its period vector against the harmonic 1-cochain basis (the class
    map H_1 -> R^b1 is linear, so ANY class-valued function of a cycle -- e.g. geometric winding numbers --
    is a fixed linear image of its period vector).  This lets a caller classify thousands of candidate
    cycles in O(b1) each and pick a representative with prescribed class AND good geometry, which a fixed
    b1-sized generator basis cannot offer (its representatives may all be mixed-class).

    Returns ``(b1, Pi, expand, cotree)``:
      * ``b1``      : dim H^1 (= 2*genus for a closed orientable surface).
      * ``Pi``      : (ncotree x b1) periods of the harmonic basis over each fundamental cycle.
      * ``expand(k) -> [v0, v1, ...]`` : vertex loop (implicitly closed) of fundamental cycle ``k``.
      * ``cotree``  : the cycles' cotree edge indices into the complex's ``EI``.
    """
    b1, harm, ctx, _vals = surface_cohomology(tris, nv, maxgen, tol)
    V, EI = ctx[4], ctx[3]
    if b1 == 0:
        return 0, np.zeros((0, 0)), None, []
    Pi, cotree = _periods(harm, V, EI)

    # Vertex-parent tree from the (deterministic BFS) spanning forest.
    _order, parent_edge, _cotree2 = _spanning_tree(V, EI)
    parent = [None] * V
    for w in range(V):
        ei = parent_edge[w]
        if ei is None:
            continue
        lo, hi = EI[ei]
        parent[w] = lo if w == hi else hi

    def _tree_path(a, b):
        anc = []
        u = a
        while u is not None:
            anc.append(u)
            u = parent[u]
        s = set(anc)
        pb = []
        u = b
        while u not in s:
            pb.append(u)
            u = parent[u]
        return anc[:anc.index(u) + 1] + pb[::-1]

    def expand(k):
        lo, hi = EI[cotree[int(k)]]
        return _tree_path(lo, hi)

    return b1, Pi, expand, cotree


def surface_homology_loops(tris, nv=None):
    """Vertex-cycle representatives of a homology-generator basis of a triangulated surface.

    Returns a list of b1 (= 2*genus, closed orientable) vertex loops, each a list of vertex ids that is
    IMPLICITLY closed (last -> first is a mesh edge).  Generator selection goes through the PERIOD MATRIX of
    the harmonic 1-cochains (tree gauge + QR pivoting on the cotree periods, ``_periods``) -- independence of
    the loops is certified by a nonsingular b1 x b1 period matrix, not by a combinatorial dual-tree argument.

    NOTE: the representatives are only guaranteed INDEPENDENT, not class-pure (a genus-1 basis may come out
    as {toroidal+poloidal, poloidal}).  A caller that needs a PRESCRIBED class (e.g. the pure toroidal cut of
    the ``radia.bem_loop_extension`` loop DOF) should classify all fundamental cycles via
    ``surface_fundamental_cycles`` instead.  This is the surface H^1 sibling of the volume
    ``cohomology_basis`` the T-Omega cut consumes.
    """
    b1, Pi, expand, _cotree = surface_fundamental_cycles(tris, nv)
    if b1 == 0:
        return []
    _, _, piv = sla.qr(Pi.T, pivoting=True, mode="economic")   # most-independent cotree rows first
    return [expand(k) for k in piv[:b1]]


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
