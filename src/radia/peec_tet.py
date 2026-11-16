"""Tetrahedral volume PEEC of a conductor with lowest-order face currents.

NGSolve owns the space: ``HDiv(order=0, RT=True)`` supplies the orientation and
scaling of every face current, the resistance ``int u.v/sigma``, the discrete
divergence against cell constants and the drive.  The native kernel
(``peec_matrices.tet_inductance``, C++ ``rad_peec_tet``) assembles the Galerkin
partial inductance and ``peec_matrices.tet_kkt_solve`` solves the dense
circuit ``[R + j omega L, D^T; D, 0] [I; phi] = [b; 0]``.

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
    scale = np.abs(W).max()
    if np.abs(pred - W[:, 2]).max() > 1e-9*scale:
        raise ValueError("peec_tet: the face basis is not affine on some element (curved geometry?)")
    coef = a / r[order]
    return pts, tets, order, coef


def solve_tet_peec(mesh, sigma, frequency_hz, cut, region=None, near_factor=3.0, moderate_factor=10.0,
                   near_order=3, max_dense_gib=16.0):
    """Impedance and losses of a conductor driven by a loop voltage across the internal surface ``cut``.

    Args:
        mesh: NGSolve mesh of the conductor only (flat tetrahedra, one connected body).
        sigma: conductivity [S/m], number or NGSolve CoefficientFunction.
        frequency_hz: drive frequency [Hz] (> 0).
        cut: boundary name of the internal surface carrying the voltage gap.
        region: optional dict name -> volume region regex; losses are integrated per entry.
        max_dense_gib: memory budget of the dense route (L, the circuit and its factor); larger
            problems raise before assembly.

    Returns:
        dict with R_ohm, L_H, impedance_ohm, power_balance, region_loss_W_per_A2, total_loss_W_per_A2,
        true_relative_residual, unknowns, tets, assembly_s, solve_s.
    """
    import ngsolve as ng
    if not frequency_hz > 0:
        raise ValueError("frequency_hz must be positive")
    if cut not in mesh.GetBoundaries():
        raise ValueError(f"cut surface {cut!r} not among the mesh boundaries")
    sig = ng.CoefficientFunction(sigma)
    fes = ng.HDiv(mesh, order=0, RT=True)
    pts, tets, order, coef = _face_basis(mesh, fes)
    counts = np.bincount(order.ravel(), minlength=fes.ndof)
    interior = np.flatnonzero(counts == 2)                      # outer-surface faces carry no current
    compact = -np.ones(fes.ndof, dtype=np.int64)
    compact[interior] = np.arange(len(interior))
    face_dof = compact[order].astype(np.int32)
    n, N = len(interior), len(interior) + mesh.ne - 1
    need = (2*8.0*n*n + 2*16.0*N*N)/2**30                      # L twice (native + array), circuit and factor
    if not need <= max_dense_gib:
        raise MemoryError(f"peec_tet: dense route needs {need:.1f} GiB for {n} face unknowns, above "
                          f"max_dense_gib={max_dense_gib}; a compressed operator is required")
    t0 = time.perf_counter()
    L = tet_inductance(pts, tets, coef, face_dof, len(interior), near_factor, moderate_factor, near_order)
    t1 = time.perf_counter()
    u, v = fes.TnT()
    l2 = ng.L2(mesh, order=0)
    q = l2.TestFunction()
    Rf = ng.BilinearForm(fes)
    Rf += (1/sig)*u*v*ng.dx
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
    sol = tet_kkt_solve(L, compact[np.asarray(rr)][keep].astype(np.int32),
                        compact[np.asarray(rc)][keep].astype(np.int32), np.asarray(rv)[keep],
                        mesh.ne - 1, (dr[dkeep] - 1).astype(np.int32),
                        compact[dc[dkeep]].astype(np.int32), dv[dkeep], g.astype(complex),
                        2*np.pi*frequency_hz, float(max_dense_gib))
    t2 = time.perf_counter()
    current = np.asarray(sol["current"])
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
    omega = 2*np.pi*frequency_hz
    return dict(R_ohm=z.real, L_H=z.imag/omega, impedance_ohm=z,
                power_balance=abs(2*total_loss - z.real)/z.real, total_loss_W_per_A2=total_loss,
                region_loss_W_per_A2=regional, true_relative_residual=float(sol["true_relative_residual"]),
                unknowns=int(len(interior)), tets=int(mesh.ne), assembly_s=t1 - t0, solve_s=t2 - t1)
