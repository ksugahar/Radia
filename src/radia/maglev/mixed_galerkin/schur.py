"""schur.py -- mixed Galerkin admittance in Schur (coupled) form on a box.

The bulk eigenmodes phi_n (Dirichlet Laplacian, vanish on the boundary) and
the surface envelope psi (vanishes on the boundary, tends to -1 inside)
share ONE Galerkin space.  For the scalar diffusion problem

    (-Laplacian + t^2) v = -t^2 f,   v = 0 on dOmega,   t^2 = s mu sigma,

with driving function f, the trial function v = sum_n xi_n phi_n + xi_s f psi
gives the block system

    [ K_bb  K_bs ] [xi_b]        [b_b]
    [ K_sb  K_ss ] [xi_s]  = -t^2 [b_s],

    K_bb = diag(kappa_n + t^2)                 (M-orthonormal eigenmodes)
    K_bs = (kappa_n + t^2) int phi_n f psi dV  (Green: psi f = 0 on dOmega)
    K_ss = int |grad(f psi)|^2 + t^2 int (f psi)^2 dV
    b_b  = int phi_n f dV,  b_s = int f^2 psi dV,

and the admittance Y(s) = sigma [ int f^2 dV + xi_b . b_b + xi_s b_s ].
Eliminating xi_b is the Schur complement S(s) = K_ss - K_sb K_bb^-1 K_bs,
the Gram-Schmidt residual of the envelope against the bulk modes: the
crossover between the Foster ladder and the sqrt(s) tail is decided by the
projection, nothing is fitted.  Y is exact at DC (xi -> 0) and expels the
field at high frequency through the "-1" of the envelope.

With P driving functions f_p the surface family is {f_p psi} and every
block gains port indices; Y(s) is then the P x P admittance matrix of
`alpha.bulk_foster_matrix_via_eigen` / `alpha.Y_matrix_mixed`, but built by
the projection instead of the additive Y_bulk + K/sqrt(s) + c_1/s.

Scope: a BOX conductor [x0, x0+Lx] x [y0, y0+Ly] x [z0, z0+Lz].  The tensor
envelope psi = f_x(x) f_y(y) f_z(z), f_i(u) = cosh(t(u - L_i/2))/cosh(t L_i/2) - 1,
is the exact boundary-layer solution of the right-angle wedge and corner, so
the box is the one polyhedron where a single envelope carries face, edge and
corner at once.  General polyhedra need per-face, per-edge and per-corner
layer quadrature and are not implemented here.

Layer integrals.  The eigenmodes are finite-element functions that the mesh
resolves; the envelope varies on the skin depth 1/|t|, which the mesh does
not.  The coupling integrals are therefore computed on a tensor quadrature
grid graded geometrically towards each face (down to a_min = 1e-5 L, well
below the smallest skin depth of the intended band), on which the
eigenmodes are sampled ONCE by vectorised point evaluation; per frequency
only the separable envelope weights change.  The 1-D rule is checked
against the closed forms of `box_layer_1d` in the tests.

Validation: validation_test/mixed_galerkin/cube3d (closed-form K_ss,
analytic-sine bulk, NGSolve ground truth) and the golden test
tests/test_maglev_mixed_galerkin_golden.py.
"""
from __future__ import annotations

import math

import numpy as np

from .alpha import _dirichlet_eigenmodes


# ---------------------------------------------------------------------------
# 1-D layer function and its closed-form integrals
# ---------------------------------------------------------------------------
def box_layer_1d(t, L):
    """Closed-form integrals of the 1-D layer function on [0, L].

    f(u) = c(u) - 1,  c(u) = (exp(-t u) + exp(-t (L - u))) / (1 + exp(-t L)).

    Returns (F0, F2, D2) = (int f du, int f^2 du, int f'^2 du).  Written
    with q = exp(-t L) so that Re t >= 0 never overflows.
    """
    t = complex(t)
    q = np.exp(-t * L)
    C1 = (2.0 / t) * (1.0 - q) / (1.0 + q)                   # int c
    C2 = ((1.0 - q * q) / t + 2.0 * L * q) / (1.0 + q) ** 2   # int c^2
    D2 = t * t * ((1.0 - q * q) / t - 2.0 * L * q) / (1.0 + q) ** 2  # int c'^2
    F0 = C1 - L
    F2 = C2 - 2.0 * C1 + L
    return F0, F2, D2


def layer_values(u, t, L):
    """f(u) and f'(u) of the 1-D layer function at the points u (array)."""
    u = np.asarray(u, dtype=float)
    q = np.exp(-t * L)
    e0 = np.exp(-t * u)
    e1 = np.exp(-t * (L - u))
    f = (e0 + e1) / (1.0 + q) - 1.0
    df = t * (-e0 + e1) / (1.0 + q)
    return f, df


# ---------------------------------------------------------------------------
# graded tensor quadrature
# ---------------------------------------------------------------------------
def graded_axis_rule(L, *, a_min_frac=1e-5, ratio=2.0, max_panel_frac=0.125,
                     nodes_per_panel=3):
    """Gauss-Legendre panels on [0, L], graded geometrically towards both ends.

    Panel widths grow by `ratio` from `a_min_frac * L` at each wall until
    they reach `max_panel_frac * L`, then stay uniform; the rule is mirror
    symmetric about L/2.  Returns (nodes, weights); no node lies on a wall.

    The defaults give 108 nodes per axis and integrate the layer function
    f, f^2 and f'^2 of `box_layer_1d` to 3e-6, 5e-6 and 2.5e-4 for every
    |t| L between 0.1 and 3400 (the copper cube from 1 Hz to 1 GHz);
    `a_min_frac` must stay below about 0.1 / (|t| L) at the top of the band.
    (ratio=2.5 gives 84 nodes and 1.4e-3 on f'^2; nodes_per_panel=4 with
    ratio=2.0 gives 144 nodes and 1.7e-5.)
    """
    edges = [0.0, a_min_frac * L]
    while True:
        step = min(edges[-1] * (ratio - 1.0), max_panel_frac * L)
        nxt = edges[-1] + step
        if nxt >= 0.5 * L - 1e-12 * L:
            edges.append(0.5 * L)
            break
        edges.append(nxt)
    # a sliver at the centre is merged into its neighbour
    if len(edges) > 3 and (edges[-1] - edges[-2]) < 0.25 * (edges[-2] - edges[-3]):
        del edges[-2]
    half = np.array(edges)
    full = np.concatenate([half, L - half[-2::-1]])
    gx, gw = np.polynomial.legendre.leggauss(nodes_per_panel)
    nodes, weights = [], []
    for a, b in zip(full[:-1], full[1:]):
        nodes.append(0.5 * (b - a) * gx + 0.5 * (a + b))
        weights.append(0.5 * (b - a) * gw)
    return np.concatenate(nodes), np.concatenate(weights)


def check_axis_rule(L, t, nodes, weights):
    """Relative errors of the rule on (int f, int f^2, int f'^2) at this t."""
    f, df = layer_values(nodes, t, L)
    F0, F2, D2 = box_layer_1d(t, L)
    got = (weights @ f, weights @ (f * f), weights @ (df * df))
    return tuple(abs(g - e) / abs(e) for g, e in zip(got, (F0, F2, D2)))


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
def box_of_mesh(mesh, rel_tol=1e-6):
    """(origin, dims) of the axis-aligned box that the mesh fills.

    Raises if a boundary vertex is not on a face of the bounding box: the
    tensor envelope is only right for a box.
    """
    from ngsolve import BND
    pts = np.array([v.point for v in mesh.vertices], dtype=float)
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    dims = hi - lo
    if np.any(dims <= 0):
        raise ValueError("mesh is not three-dimensional: bounding box %r" % (dims,))
    bnd = set()
    for el in mesh.Elements(BND):
        for v in el.vertices:
            bnd.add(v.nr)
    tol = rel_tol * dims.max()
    for nr in bnd:
        p = pts[nr]
        on_face = np.any(np.abs(p - lo) < tol) or np.any(np.abs(p - hi) < tol)
        if not on_face:
            raise ValueError(
                "boundary vertex %r is not on the bounding box %r..%r: the "
                "tensor envelope needs a box conductor" % (p, lo, hi))
    return lo, dims


# ---------------------------------------------------------------------------
# the coupled model
# ---------------------------------------------------------------------------
class BoxMixedGalerkin:
    """Mixed Galerkin (Schur form) eddy-current admittance of a box conductor.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        Tetrahedral mesh of the box conductor; the whole boundary carries
        `dirichlet_label`.
    sigma, mu : float
        Conductivity (S/m) and permeability (H/m).
    drive_cfs : sequence of CoefficientFunction or float
        Driving functions f_p (default [1.0], the monopole port).
    n_eigen : int
        Number of bulk Dirichlet eigenmodes.
    dirichlet_label : str
        Boundary label of the conductor surface (default "outer").
    rule : dict, optional
        Keyword arguments of `graded_axis_rule`.

    Memory: the eigenmodes are stored on the grid, 8 bytes per mode and grid
    point (n_eigen = 40 on the default 108^3 grid is 400 MB); pass a coarser
    `rule` or fewer modes when that is too much.  Setup samples the modes
    once (a few seconds); each Y(s) is a tensor contraction plus an
    (n_eigen + P) x (n_eigen + P) solve.

    Attributes
    ----------
    lam, tau, V, B, M, origin, dims, nodes, weights, phi, F, G, W3
    """

    def __init__(self, mesh, sigma, mu, drive_cfs=(1.0,), n_eigen=40,
                 dirichlet_label="outer", rule=None):
        from ngsolve import CoefficientFunction, GridFunction, Integrate, x, y, z

        self.sigma = float(sigma)
        self.mu = float(mu)
        self.origin, self.dims = box_of_mesh(mesh)

        lam, vecs, Mmass, free, fes, V = _dirichlet_eigenmodes(mesh, n_eigen, dirichlet_label)
        self.lam = lam
        self.tau = self.mu * self.sigma / lam
        self.V = V

        cfs = [CoefficientFunction(c) for c in drive_cfs]
        P = len(cfs)
        self.P = P

        # drive Gram M_pq = int f_p f_q dV and projections B_np = int phi_n f_p dV
        self.M = np.zeros((P, P))
        for p in range(P):
            for q in range(p, P):
                val = float(Integrate(cfs[p] * cfs[q], mesh, order=4))
                self.M[p, q] = self.M[q, p] = val
        N = vecs.shape[1]
        self.B = np.zeros((N, P))
        gfu = GridFunction(fes)
        for p, cf in enumerate(cfs):
            gfu.Set(cf)
            full = np.array(gfu.vec.FV().NumPy())
            self.B[:, p] = vecs.T @ (Mmass[free, :] @ full)

        # graded tensor grid
        rule = dict(rule or {})
        axes = [graded_axis_rule(float(L), **rule) for L in self.dims]
        self.nodes = [a[0] for a in axes]
        self.weights = [a[1] for a in axes]
        X, Y, Z = np.meshgrid(self.nodes[0] + self.origin[0],
                              self.nodes[1] + self.origin[1],
                              self.nodes[2] + self.origin[2], indexing="ij")
        shape = X.shape
        mips = mesh(X.ravel(), Y.ravel(), Z.ravel())

        # eigenmodes sampled once on the grid
        self.phi = np.empty((N,) + shape)
        gf = GridFunction(fes)
        full = np.zeros(fes.ndof)
        for n in range(N):
            full[:] = 0.0
            full[free] = vecs[:, n]
            gf.vec.FV().NumPy()[:] = full
            self.phi[n] = np.asarray(gf(mips)).reshape(shape)

        # drives and their gradients on the grid
        self.F = np.empty((P,) + shape)
        self.G = np.empty((P, 3) + shape)
        for p, cf in enumerate(cfs):
            self.F[p] = np.asarray(cf(mips)).reshape(shape)
            for c, xc in enumerate((x, y, z)):
                self.G[p, c] = np.asarray(cf.Diff(xc)(mips)).reshape(shape)

        self.W3 = (self.weights[0][:, None, None] * self.weights[1][None, :, None]
                   * self.weights[2][None, None, :])

    # -- per-frequency pieces -------------------------------------------------
    def _t(self, s):
        t = complex(np.sqrt(complex(s) * self.mu * self.sigma))
        if t.real <= 0.0:
            raise ValueError(
                "s = %r gives t = sqrt(s mu sigma) = %r with Re t <= 0; the "
                "envelope needs Re t > 0 (s = 0 is handled as the exact DC "
                "limit)" % (s, t))
        return t

    def blocks(self, s):
        """(K_bs, K_ss, b_s, I) for complex s with Re t > 0.

        K_bs : (N, P) coupling, K_ss : (P, P) surface block,
        b_s : (P, P) with b_s[p, q] = int f_p f_q psi dV,
        I : (N, P) overlaps int phi_n f_p psi dV.
        """
        t = self._t(s)
        t2 = t * t
        f1, d1 = [], []
        for L, u in zip(self.dims, self.nodes):
            f, df = layer_values(u, t, float(L))
            f1.append(f)
            d1.append(df)
        psi = f1[0][:, None, None] * f1[1][None, :, None] * f1[2][None, None, :]
        gpsi = np.stack([
            d1[0][:, None, None] * f1[1][None, :, None] * f1[2][None, None, :],
            f1[0][:, None, None] * d1[1][None, :, None] * f1[2][None, None, :],
            f1[0][:, None, None] * f1[1][None, :, None] * d1[2][None, None, :],
        ])
        Wpsi = self.W3 * psi
        P, N = self.P, self.phi.shape[0]

        I = np.empty((N, P), dtype=complex)
        for p in range(P):
            I[:, p] = np.tensordot(self.phi, Wpsi * self.F[p], axes=([1, 2, 3], [0, 1, 2]))
        K_bs = (self.lam[:, None] + t2) * I

        # surface trial functions f_p psi and their gradients on the grid
        Fp = self.F
        grad_fp_psi = Fp[:, None] * gpsi[None] + self.G * psi[None, None]   # (P, 3, grid)
        K_ss = np.empty((P, P), dtype=complex)
        b_s = np.empty((P, P), dtype=complex)
        for p in range(P):
            for q in range(p, P):
                dot = np.sum(grad_fp_psi[p] * grad_fp_psi[q], axis=0)
                val = np.sum(self.W3 * (dot + t2 * Fp[p] * Fp[q] * psi * psi))
                K_ss[p, q] = K_ss[q, p] = val
                bval = np.sum(Wpsi * Fp[p] * Fp[q])
                b_s[p, q] = b_s[q, p] = bval
        return K_bs, K_ss, b_s, I

    # -- admittance -----------------------------------------------------------
    def Y(self, s):
        """P x P admittance matrix at complex s (Re s >= 0); scalar for P = 1."""
        if s == 0:
            Y = self.sigma * self.M
            return Y[0, 0] if self.P == 1 else Y
        t = self._t(s)
        t2 = t * t
        K_bs, K_ss, b_s, _ = self.blocks(s)
        N, P = K_bs.shape
        A = np.zeros((N + P, N + P), dtype=complex)
        A[np.arange(N), np.arange(N)] = self.lam + t2
        A[:N, N:] = K_bs
        A[N:, :N] = K_bs.T
        A[N:, N:] = K_ss
        rhs = np.empty((N + P, P), dtype=complex)
        rhs[:N] = -t2 * self.B
        rhs[N:] = -t2 * b_s
        # symmetric diagonal scaling keeps the envelope block, which scales
        # like t^12 at low frequency, on the same footing as the bulk
        d = 1.0 / np.sqrt(np.abs(np.diag(A)))
        xi = d[:, None] * np.linalg.solve(d[:, None] * A * d[None, :], d[:, None] * rhs)
        Y = self.sigma * (self.M + self.B.T @ xi[:N] + b_s.T @ xi[N:])
        return Y[0, 0] if P == 1 else Y

    def Y_bulk(self, s):
        """Bulk modes alone (same truncation, no envelope); exact at DC."""
        if s == 0:
            Y = self.sigma * self.M
        else:
            t2 = self._t(s) ** 2
            Y = self.sigma * (self.M - t2 * (self.B.T @ (self.B / (self.lam[:, None] + t2))))
        return Y[0, 0] if self.P == 1 else Y

    def alpha(self, s):
        """Polarizability alpha(s) = V I - Y(s) / sigma."""
        Y = self.Y(s)
        if self.P == 1:
            return self.V - Y / self.sigma
        return self.V * np.eye(self.P) - Y / self.sigma
