"""Tetrahedral volume PEEC of a conductor with lowest-order face currents.

NGSolve owns the space: ``HDiv(order=0, RT=True)`` supplies the orientation and
scaling of every face current, the resistance ``int u.v/sigma``, the discrete
divergence against cell constants and the drive.  The native kernel
(``peec_matrices.tet_inductance``, C++ ``rad_peec_tet``) assembles the Galerkin
partial inductance and ``peec_matrices.tet_kkt_solve`` solves the dense
circuit ``[R + j omega L, D^T; D, 0] [I; phi] = [b; 0]``.  With
``operator="hmatrix"`` the same inductance is compressed by HACApK
(``radia._radia_pybind.TetPEECHMatrix``) and the circuit is solved by GMRES,
right-preconditioned by the exact constraint preconditioner of
``Z0 = diag(R) + j omega diag(L)``, whose cell Schur complement is factored by
NGSolve's sparsecholesky.

A loop drive applies the voltage across a named internal cut surface; faces on
the conductor's outer surface carry no current.  Flat tetrahedra only.
Normalization: per 1 A peak terminal current, ``exp(+j omega t)``.
"""

from __future__ import annotations

import time

import numpy as np

from .peec_matrices import tet_inductance, tet_kkt_solve


def _face_basis(mesh, fes):
    """Per (tetrahedron, local face k): global DOF and c in w = c (x - v_k), from NGSolve's basis."""
    import ngsolve as ng
    if mesh.GetCurveOrder() > 1 or mesh.deformation is not None:
        raise ValueError("peec_tet: curved or deformed meshes are not supported (flat tetrahedra only)")
    els = list(mesh.Elements(ng.VOL))
    if any(el.type != ng.ET.TET for el in els):
        raise ValueError("peec_tet: every volume element must be a tetrahedron")
    pts = np.array([mesh[v].point for v in mesh.vertices], dtype=float)
    tets = np.array([[v.nr for v in el.vertices] for el in els], dtype=np.int32)
    dofs = np.array([list(fes.GetDofNrs(el)) for el in els], dtype=np.int64)
    if dofs.shape[1] != 4:
        raise ValueError("peec_tet: expected one lowest-order DOF per face")
    V = pts[tets]                                               # (m, 4, 3)
    # The DOF of the face opposite vertex k: the one whose face omits vertex k (from the face numbers).
    face_of_dof = {}
    for el in els:
        for f, d in zip(el.faces, fes.GetDofNrs(el)):
            face_of_dof[d] = f.nr
    face_verts = {f.nr: frozenset(v.nr for v in mesh[f].vertices) for f in mesh.faces}
    order = np.empty_like(dofs)
    for e, el in enumerate(els):
        for d in dofs[e]:
            missing = set(tets[e]) - face_verts[face_of_dof[int(d)]]
            k = list(tets[e]).index(missing.pop())
            order[e, k] = d
    # One random coefficient vector: w(x) = sum_k a_k (x - v_k) in each element, a_k = r_d c_k.
    rng = np.random.default_rng(20261116)
    gf = ng.GridFunction(fes)
    r = rng.uniform(0.5, 1.5, fes.ndof)
    gf.vec.FV().NumPy()[:] = r
    lam = np.array([[0.25, 0.25, 0.25, 0.25], [0.4, 0.2, 0.2, 0.2], [0.1, 0.2, 0.3, 0.4]])
    X = np.einsum("pk,mkj->mpj", lam, V)                        # (m, 3, 3) interior points
    flat = X.reshape(-1, 3)
    mp = mesh(flat[:, 0], flat[:, 1], flat[:, 2])
    if not np.array_equal(np.asarray(mp["nr"]).reshape(-1, 3), np.repeat(np.arange(len(els))[:, None], 3, 1)):
        raise RuntimeError("peec_tet: interior evaluation points were located in another element")
    W = np.asarray(gf(mp)).reshape(-1, 3, 3)
    A = np.einsum("mj,mj->m", W[:, 1] - W[:, 0], X[:, 1] - X[:, 0]) / np.einsum(
        "mj,mj->m", X[:, 1] - X[:, 0], X[:, 1] - X[:, 0])
    b = W[:, 0] - A[:, None]*X[:, 0]                            # w = A x + b, b = -sum a_k v_k
    M = np.concatenate([np.transpose(V, (0, 2, 1)), np.ones((len(els), 1, 4))], axis=1)   # (m, 4, 4)
    a = np.linalg.solve(M, np.concatenate([-b, A[:, None]], axis=1)[..., None])[..., 0]
    pred = A[:, None]*X[:, 2] + b
    scale = np.abs(W).max(axis=(1, 2))                         # per element
    if np.any(np.abs(pred - W[:, 2]).max(axis=1) > 1e-9*scale):
        raise ValueError("peec_tet: the face basis is not affine on some element (curved geometry?)")
    coef = a / r[order]
    return pts, tets, order, coef


def _check_one_body(order, interior, ne):
    """Refuse meshes whose tetrahedra do not form one face-connected body (one cell is grounded)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    e, _ = np.nonzero(np.isin(order, interior))
    d = order[np.isin(order, interior)]
    pairs = np.argsort(d, kind="stable")
    a, b = e[pairs][0::2], e[pairs][1::2]                       # each interior face joins two cells
    n_comp, _ = connected_components(coo_matrix((np.ones(len(a)), (a, b)), shape=(ne, ne)), directed=False)
    if n_comp != 1:
        raise ValueError(f"peec_tet: the conductor mesh has {n_comp} face-connected bodies; one is supported")


def _schur_inverse(mesh, s_full):
    """sparsecholesky of the cell matrix ``s_full`` (scipy, ne x ne, cell 0 grounded) in NGSolve storage."""
    import ngsolve as ng
    l2 = ng.L2(mesh, order=0, dgjumps=True, complex=True)
    u, v = l2.TnT()
    bf = ng.BilinearForm(l2, symmetric=False)
    bf += u*v*ng.dx
    bf += (u - u.Other())*(v - v.Other())*ng.dx(skeleton=True)     # pattern: cells sharing a face
    bf.Assemble()
    r, c, _ = bf.mat.COO()
    r, c = np.asarray(r), np.asarray(c)
    vals = np.asarray(s_full[r, c]).ravel()
    if not np.isclose(np.abs(vals).sum(), abs(s_full).sum(), rtol=1e-12, atol=0.0):
        raise RuntimeError("peec_tet: the Schur complement leaves the face-neighbour pattern")
    bf.mat.AsVector().FV().NumPy()[:] = vals
    free = ng.BitArray(mesh.ne)
    free.Set()
    free.Clear(0)
    inv = bf.mat.Inverse(free, inverse="sparsecholesky")
    work_in, work_out = bf.mat.CreateColVector(), bf.mat.CreateColVector()

    def solve(t):
        work_in.FV().NumPy()[0] = 0.0
        work_in.FV().NumPy()[1:] = t
        work_out.data = inv*work_in
        return work_out.FV().NumPy()[1:].copy()
    return solve


def _solve_kkt_hmatrix(mesh, H, n, m, rr, rc, rv, dr, dc, dv, g, omega, rtol, max_iter, restart):
    """GMRES on [R + j omega L, D^T; D, 0], right-preconditioned by the exact constraint preconditioner of
    Z0 = diag(R) + j omega diag(L)."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import LinearOperator, gmres
    Rs = sp.csr_matrix((rv, (rr, rc)), shape=(n, n))
    Ds = sp.csr_matrix((dv, (dr, dc)), shape=(m, n))
    DsT = Ds.T.tocsr()
    z0 = Rs.diagonal() + 1j*omega*np.asarray(H.diagonal())
    if not np.all(z0.real > 0):
        raise RuntimeError("peec_tet: nonpositive resistance diagonal")
    S = (Ds @ sp.diags(1.0/z0) @ DsT).tocoo()
    s_full = sp.csr_matrix((np.concatenate([S.data, [1.0]]),
                            (np.concatenate([S.row + 1, [0]]), np.concatenate([S.col + 1, [0]]))),
                           shape=(m + 1, m + 1))
    schur = _schur_inverse(mesh, s_full)

    def apply_A(x):
        u, p = x[:n], x[n:]
        return np.concatenate([Rs @ u + 1j*omega*np.asarray(H.matvec(u)) + DsT @ p, Ds @ u])

    def apply_P(y):
        f, q = y[:n], y[n:]
        p = schur(Ds @ (f/z0) - q)
        return np.concatenate([(f - DsT @ p)/z0, p])

    N = n + m
    b = np.concatenate([g.astype(complex), np.zeros(m, dtype=complex)])
    count = [0]

    def cb(_):
        count[0] += 1
    AP = LinearOperator((N, N), matvec=lambda y: apply_A(apply_P(y)), dtype=complex)
    y, info = gmres(AP, b, rtol=rtol, atol=0.0, restart=restart, maxiter=max_iter, callback=cb,
                    callback_type="pr_norm")
    x = apply_P(y)
    res = np.linalg.norm(b - apply_A(x))/np.linalg.norm(b)
    if info != 0 or not res <= 1e-6:
        raise RuntimeError(f"peec_tet: GMRES did not converge (info={info}, true relative residual {res:.3e}, "
                           f"{count[0]} iterations)")
    return x[:n], float(res), count[0]


def solve_tet_peec(mesh, sigma, frequency_hz, cut, region=None, operator="dense", near_factor=3.0,
                   moderate_factor=10.0, near_order=3, max_dense_gib=16.0, aca_eps=1e-8, leaf_size=64, eta=2.0,
                   max_rank=400, rtol=1e-9, max_iter=50, restart=200):
    """Impedance and losses of a conductor driven by a loop voltage across the internal surface ``cut``.

    Args:
        mesh: NGSolve mesh of the conductor only (flat, undeformed tetrahedra forming one face-connected
            body; checked).
        sigma: conductivity [S/m], number or real scalar CoefficientFunction; it must be finite and
            positive at the resistance quadrature points (checked).
        frequency_hz: drive frequency [Hz] (> 0), or a sequence of them sharing one assembled L.
        cut: boundary name of the internal surface carrying the voltage gap.
        region: optional dict name -> volume region regex; losses are integrated per entry.
        operator: ``"dense"`` (dense L, LU) or ``"hmatrix"`` (HACApK L, preconditioned GMRES).
        max_dense_gib: memory budget of the dense route (L, the circuit and its factor); larger
            problems raise before assembly.
        aca_eps, leaf_size, eta, max_rank: HACApK parameters of the hmatrix route.
        rtol, max_iter, restart: GMRES tolerance, restart cycles and restart length (hmatrix route);
            the true relative residual of the circuit must reach 1e-6.

    Returns:
        dict (a list of dicts for a sequence of frequencies) with frequency_hz, R_ohm, L_H, impedance_ohm,
        power_balance, region_loss_W_per_A2, total_loss_W_per_A2, true_relative_residual, unknowns, tets,
        assembly_s, solve_s, operator, iterations and, for the hmatrix route, the H-matrix statistics.
    """
    import ngsolve as ng
    freqs = np.atleast_1d(np.asarray(frequency_hz, dtype=float))
    if freqs.ndim != 1 or freqs.size == 0 or not np.all(freqs > 0):
        raise ValueError("frequency_hz must be positive (a number or a sequence)")
    if cut not in mesh.GetBoundaries():
        raise ValueError(f"cut surface {cut!r} not among the mesh boundaries")
    ir = ng.IntegrationRule(ng.ET.TET, 2)                       # exact for the affine face currents
    sig = ng.CoefficientFunction(sigma)
    if sig.dim != 1 or sig.is_complex:
        raise ValueError("sigma must be a real scalar conductivity")
    sig_q = np.asarray(sig(mesh.MapToAllElements(ir, ng.VOL))).ravel()
    if not (np.all(np.isfinite(sig_q)) and np.all(sig_q > 0)):
        raise ValueError(f"sigma must be finite and positive at the resistance quadrature points "
                         f"(min {np.min(sig_q):.6g})")
    fes = ng.HDiv(mesh, order=0, RT=True)
    pts, tets, order, coef = _face_basis(mesh, fes)
    counts = np.bincount(order.ravel(), minlength=fes.ndof)
    interior = np.flatnonzero(counts == 2)                      # outer-surface faces carry no current
    compact = -np.ones(fes.ndof, dtype=np.int64)
    compact[interior] = np.arange(len(interior))
    face_dof = compact[order].astype(np.int32)
    _check_one_body(order, interior, mesh.ne)
    if operator not in ("dense", "hmatrix"):
        raise ValueError(f"operator must be 'dense' or 'hmatrix', not {operator!r}")
    n, N = len(interior), len(interior) + mesh.ne - 1
    need = (2*8.0*n*n + 2*16.0*N*N)/2**30                      # L twice (native + array), circuit and factor
    if operator == "dense" and not need <= max_dense_gib:
        raise MemoryError(f"peec_tet: dense route needs {need:.1f} GiB for {n} face unknowns, above "
                          f"max_dense_gib={max_dense_gib}; use operator='hmatrix'")
    t0 = time.perf_counter()
    if operator == "dense":
        L = tet_inductance(pts, tets, coef, face_dof, n, near_factor, moderate_factor, near_order)
    else:
        from ._radia_pybind import TetPEECHMatrix
        H = TetPEECHMatrix(pts, tets, coef, face_dof, n, near_factor, moderate_factor, near_order,
                           aca_eps, leaf_size, eta, max_rank)
    t1 = time.perf_counter()
    u, v = fes.TnT()
    l2 = ng.L2(mesh, order=0)
    q = l2.TestFunction()
    Rf = ng.BilinearForm(fes)
    Rf += (1/sig)*u*v*ng.dx(intrules={ng.ET.TET: ir})
    Rf.Assemble()
    Df = ng.BilinearForm(trialspace=fes, testspace=l2)
    Df += ng.div(u)*q*ng.dx
    Df.Assemble()
    gform = ng.LinearForm(fes)
    gform += v.Trace()*ng.specialcf.normal(3)*ng.ds(cut)
    gform.Assemble()
    g = gform.vec.FV().NumPy()[interior].copy()
    rr, rc, rv = Rf.mat.COO()
    keep = (compact[np.asarray(rr)] >= 0) & (compact[np.asarray(rc)] >= 0)
    dr, dc, dv = Df.mat.COO()
    dr, dc, dv = np.asarray(dr), np.asarray(dc), np.asarray(dv)
    dkeep = (compact[dc] >= 0) & (dr != 0)                     # ground cell 0
    r_row, r_col = compact[np.asarray(rr)][keep].astype(np.int32), compact[np.asarray(rc)][keep].astype(np.int32)
    d_row, d_col = (dr[dkeep] - 1).astype(np.int32), compact[dc[dkeep]].astype(np.int32)
    results = []
    for f_hz in freqs:
        results.append(_solve_frequency(
            mesh, fes, sig, interior, n, g, operator, L if operator == "dense" else H, r_row, r_col,
            np.asarray(rv)[keep], d_row, d_col, dv[dkeep], float(f_hz), region, max_dense_gib, rtol, max_iter,
            restart, t1 - t0))
    return results if np.ndim(frequency_hz) else results[0]


def _solve_frequency(mesh, fes, sig, interior, n, g, operator, Lop, r_row, r_col, r_val, d_row, d_col, d_val,
                     frequency_hz, region, max_dense_gib, rtol, max_iter, restart, assembly_s):
    import ngsolve as ng
    t1 = time.perf_counter()
    omega = 2*np.pi*frequency_hz
    extra = {}
    if operator == "dense":
        sol = tet_kkt_solve(Lop, r_row, r_col, r_val, mesh.ne - 1, d_row, d_col, d_val,
                            g.astype(complex), omega, float(max_dense_gib))
        current, residual, iterations = np.asarray(sol["current"]), float(sol["true_relative_residual"]), 0
    else:
        current, residual, iterations = _solve_kkt_hmatrix(
            mesh, Lop, n, mesh.ne - 1, r_row, r_col, r_val, d_row, d_col, d_val, g, omega, rtol, max_iter, restart)
        extra["hmatrix"] = dict(Lop.stats())
    t2 = time.perf_counter()
    total = complex(g @ current)                                # terminal current at V = 1
    z = 1.0/total
    full = np.zeros(fes.ndof, dtype=complex)
    full[interior] = current/total                              # per 1 A terminal current
    gr, gi = ng.GridFunction(fes), ng.GridFunction(fes)
    gr.vec.FV().NumPy()[:] = full.real
    gi.vec.FV().NumPy()[:] = full.imag
    dens = 0.5*(gr*gr + gi*gi)/sig
    total_loss = ng.Integrate(dens, mesh)
    regional = None if region is None else {
        k: ng.Integrate(dens*mesh.MaterialCF({rx: 1.0}, default=0.0), mesh) for k, rx in region.items()}
    return dict(frequency_hz=frequency_hz, R_ohm=z.real, L_H=z.imag/omega, impedance_ohm=z,
                power_balance=abs(2*total_loss - z.real)/z.real, total_loss_W_per_A2=total_loss,
                region_loss_W_per_A2=regional, true_relative_residual=residual,
                unknowns=int(n), tets=int(mesh.ne), assembly_s=assembly_s, solve_s=t2 - t1,
                operator=operator, iterations=int(iterations), **extra)
