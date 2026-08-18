"""CanonicalHCurl: tracking-specialized vacuum HCurl subspace on a beam chain.

The space (named 2026-08-17) is the structured HCurl(p) subspace used by the
EarlyTimes Lie/A-map routes.  Per element (one element per cross-section,
plain rectangular section, chain along the design orbit):

  * ``a_x = 0`` identically (axial-gauge contract, no DOFs),
  * median-plane parity: ``a_s`` even in ``y``, ``a_y`` odd in ``y``,
  * gauge rigidity: ``a_y`` carries no ``x^0`` modes and ``a_s`` no
    ``x^0 y^0`` modes, so the residual ``chi(y,s)`` gauge family is removed
    STRUCTURALLY and ``A`` vanishes on the design orbit,
  * vacuum: ``curl curl A = 0`` in the curved chart ``g = 1 + htilde(s) x``
    as denominator-cleared polynomial identities.

Verified structure (exact over Q by the Mathematica oracle, float-verified
here as construction self-tests):

  * dimension law: ``dim = p_x * (p_s + 1)`` exactly -- the element is
    parameterized by the midplane multipole profiles ``b_m(s)``,
  * ``p_y`` is generated, not chosen: the null space saturates at
    ``p_y <= p_x + 2`` (internal bound, never a user knob),
  * for ``htilde != 0`` the STRICT kernel collapses (h-continuations of the
    top families need out-of-truncation x-degrees), so the curved space is
    DEFINED as the fixed-dimension least-defect subspace: the
    ``p_x*(p_s+1)`` smallest right singular vectors of the constraint,
  * s-interface contract (graded ``L1``): the full tangential+normal trace
    match degenerates to analytic continuation; the production contract is
    the ``a_y`` trace (HCurl conformity, kills sheet currents / ``p_y``
    kicks) plus the ``b_m`` VALUE continuity (midplane ``a_s`` trace, kills
    phase kicks).  Chained production therefore needs ``p_s >= 2``
    (``p_s = 1`` chains degenerate to a globally linear ``b_m``).

Projection is a FULL-VOLUME 3-component B least squares (B is
gauge-invariant and curl is injective on the gauge-rigid space, so the fit
determines the canonical-gauge ``A`` uniquely); midplane-only estimation is
banned -- its own residual under-reports the true aperture error by an
order of magnitude.  The fitted chain feeds the Lie kernel through
per-segment transverse coefficient arrays (covariant ``a_s``) and the
canonical A-map RK through the same polynomials.

This module is deliberately NumPy-only: the space is a reference-element
construction plus frame bookkeeping, not an NGSolve FE space.  The NGSolve
HCurl bridge (Z-matrix) and the ring-closure variant (periodic splines plus
one cohomology DOF for the linked flux) are follow-ups.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "CanonicalHCurlElement",
    "CanonicalHCurlChain",
    "CanonicalHCurlFit",
    "graded_breaks",
]


def graded_breaks(s_nodes_m, weight, element_count, *, strength=1.0):
    """Element boundaries equidistributing the monitor ``1 + strength*w``.

    ``weight`` are non-negative monitor samples at ``s_nodes_m`` (e.g.
    ``|dB_y/ds|`` on the design orbit, normalized or not -- it is rescaled
    to unit maximum).  ``strength`` sets how much the fringe attracts
    elements: 0 reproduces the uniform split, larger values shrink elements
    where the monitor is large.  Returns ``element_count + 1`` strictly
    increasing breaks spanning ``[s_nodes[0], s_nodes[-1]]``.
    """
    s = np.asarray(s_nodes_m, dtype=float).reshape(-1)
    w = np.asarray(weight, dtype=float).reshape(-1)
    if s.size != w.size or s.size < 2 or not np.all(np.diff(s) > 0.0):
        raise ValueError("s_nodes_m must be strictly increasing and match w")
    if np.any(w < 0.0) or not np.all(np.isfinite(w)):
        raise ValueError("weight must be finite and non-negative")
    count = int(element_count)
    if count < 1:
        raise ValueError("element_count must be positive")
    peak = float(np.max(w))
    density = 1.0 + float(strength) * (w / peak if peak > 0.0 else w)
    cdf = np.concatenate(([0.0], np.cumsum(
        0.5 * (density[1:] + density[:-1]) * np.diff(s))))
    targets = np.linspace(0.0, cdf[-1], count + 1)
    breaks = np.interp(targets, cdf, s)
    breaks[0] = s[0]
    breaks[-1] = s[-1]
    if not np.all(np.diff(breaks) > 0.0):
        raise ValueError("monitor produced degenerate breaks; lower strength")
    return breaks


def _monomials(order_x, order_y, order_s, parity):
    """Exponent triples ``(i, j, k)`` with the requested y-parity."""
    out = []
    for i in range(order_x + 1):
        for j in range(order_y + 1):
            if j % 2 != parity:
                continue
            for k in range(order_s + 1):
                out.append((i, j, k))
    return out


class _DofPolynomial:
    """Sparse polynomial with DOF-linear coefficients: {exps: {dof: coef}}."""

    __slots__ = ("terms",)

    def __init__(self):
        self.terms = {}

    def add(self, exp, dof, coef):
        row = self.terms.setdefault(exp, {})
        row[dof] = row.get(dof, 0.0) + coef

    def diff(self, axis, inverse_scale):
        out = _DofPolynomial()
        for exp, row in self.terms.items():
            n = exp[axis]
            if n == 0:
                continue
            ne = list(exp)
            ne[axis] = n - 1
            ne = tuple(ne)
            for dof, c in row.items():
                out.add(ne, dof, c * n * inverse_scale)
        return out

    def plus(self, other):
        out = _DofPolynomial()
        for source in (self, other):
            for exp, row in source.terms.items():
                for dof, c in row.items():
                    out.add(exp, dof, c)
        return out

    def minus(self, other):
        out = _DofPolynomial()
        for exp, row in self.terms.items():
            for dof, c in row.items():
                out.add(exp, dof, c)
        for exp, row in other.terms.items():
            for dof, c in row.items():
                out.add(exp, dof, -c)
        return out

    def scaled(self, factor):
        out = _DofPolynomial()
        for exp, row in self.terms.items():
            for dof, c in row.items():
                out.add(exp, dof, c * factor)
        return out

    def times_monomial_sum(self, shifts_and_coefs):
        """Multiply by ``sum_t coef_t * xi^dx_t * eta^dy_t * zeta^ds_t``."""
        out = _DofPolynomial()
        for (dx, dy, ds), factor in shifts_and_coefs:
            if factor == 0.0:
                continue
            for exp, row in self.terms.items():
                ne = (exp[0] + dx, exp[1] + dy, exp[2] + ds)
                for dof, c in row.items():
                    out.add(ne, dof, c * factor)
        return out


@dataclass(frozen=True)
class _ElementOperators:
    """Constraint matrix and the flat-operator reachable exponent sets."""

    constraint: np.ndarray
    reachable_rows: int
    dropped_rows: int


@dataclass(frozen=True)
class CanonicalHCurlElement:
    """One reference element of the CanonicalHCurl space.

    Parameters
    ----------
    order_x, order_s:
        The only resolution knobs ``(p_x, p_s)``.  ``p_y`` is internal.
    half_width_m, half_height_m, half_length_m:
        Physical half-apertures ``(a_x, a_y, a_s)``; the basis lives on
        normalized coordinates ``xi = x/a_x`` etc. so design columns stay
        O(1) (raw-meter monomials are a measured Vandermonde failure).
    curvature_poly_per_m:
        Coefficients of ``htilde(zeta) = sum_k c_k zeta^k`` (physical 1/m,
        ``zeta`` the normalized s), degree at most ``order_s``.  The metric
        is ``g = 1 + htilde(zeta) * x`` and the stored ``a_s`` DOF is the
        COVARIANT component ``g * A_s_physical``.
    """

    order_x: int
    order_s: int
    half_width_m: float
    half_height_m: float
    half_length_m: float
    curvature_poly_per_m: tuple = (0.0,)
    ay_exponents: tuple = field(init=False, repr=False, default=())
    as_exponents: tuple = field(init=False, repr=False, default=())
    basis: np.ndarray = field(init=False, repr=False, default=None)
    vacuum_defects: np.ndarray = field(init=False, repr=False, default=None)
    vacuum_defect_scale: float = field(init=False, repr=False, default=0.0)

    def __post_init__(self):
        order_x = int(self.order_x)
        order_s = int(self.order_s)
        if order_x < 1 or order_s < 0:
            raise ValueError("order_x >= 1 and order_s >= 0 are required")
        scales = (float(self.half_width_m), float(self.half_height_m),
                  float(self.half_length_m))
        if not all(np.isfinite(s) and s > 0.0 for s in scales):
            raise ValueError("all half-apertures must be finite and positive")
        curvature = np.asarray(self.curvature_poly_per_m, dtype=float).reshape(-1)
        if curvature.size == 0 or curvature.size > order_s + 1 or not np.all(
                np.isfinite(curvature)):
            raise ValueError(
                "curvature_poly_per_m needs 1..order_s+1 finite coefficients")
        object.__setattr__(self, "order_x", order_x)
        object.__setattr__(self, "order_s", order_s)
        object.__setattr__(self, "curvature_poly_per_m",
                           tuple(float(c) for c in curvature))

        order_y = order_x + 2  # generated-series saturation bound (measured)
        ay_exps = tuple(e for e in _monomials(order_x, order_y, order_s, 1)
                        if e[0] >= 1)                      # gauge rigidity
        as_exps = tuple(e for e in _monomials(order_x, order_y, order_s, 0)
                        if not (e[0] == 0 and e[1] == 0))  # centre gauge
        object.__setattr__(self, "ay_exponents", ay_exps)
        object.__setattr__(self, "as_exponents", as_exps)

        operators = self._build_operators(ay_exps, as_exps, scales, curvature)
        target_dim = order_x * (order_s + 1)
        basis, defects = self._least_defect_basis(operators.constraint,
                                                  target_dim)
        scale = (float(np.linalg.norm(operators.constraint, 2))
                 if operators.constraint.size else 1.0)
        object.__setattr__(self, "basis", basis)
        object.__setattr__(self, "vacuum_defects", defects)
        object.__setattr__(self, "vacuum_defect_scale", scale)
        if basis.shape[1] != target_dim:
            raise AssertionError(
                "dimension law violated: expected "
                f"{target_dim}, built {basis.shape[1]}")

    # -- construction ----------------------------------------------------
    @property
    def dimension(self) -> int:
        return int(self.order_x) * (int(self.order_s) + 1)

    @property
    def scales(self):
        return (float(self.half_width_m), float(self.half_height_m),
                float(self.half_length_m))

    def _component_polynomials(self, ay_exps, as_exps):
        ay = _DofPolynomial()
        for idx, exp in enumerate(ay_exps):
            ay.add(exp, idx, 1.0)
        a_s = _DofPolynomial()
        for idx, exp in enumerate(as_exps):
            a_s.add(exp, len(ay_exps) + idx, 1.0)
        return ay, a_s

    def _build_operators(self, ay_exps, as_exps, scales, curvature):
        inv = tuple(1.0 / s for s in scales)
        ay, a_s = self._component_polynomials(ay_exps, as_exps)

        def d(poly, axis):
            return poly.diff(axis, inv[axis])

        # g = 1 + x*htilde(zeta): normalized-coordinate monomial sum.
        g_terms = [((0, 0, 0), 1.0)] + [
            ((1, 0, k), float(c) * scales[0])
            for k, c in enumerate(curvature)
        ]

        def g_times(poly):
            return poly.times_monomial_sum(g_terms)

        def h_times(poly):
            return poly.times_monomial_sum(
                [((0, 0, k), float(c)) for k, c in enumerate(curvature)])

        N = d(a_s, 1).minus(d(ay, 2))      # g*Bx
        M = d(a_s, 0)                      # -g*By
        P = d(ay, 0)                       # Bs
        gP = g_times(P)
        c1 = g_times(d(gP, 1)).plus(d(M, 2))
        c2 = d(N, 2).minus(g_times(d(gP, 0)))
        c3 = g_times(d(M, 0)).minus(h_times(M)).plus(g_times(d(N, 1)))

        flat = (
            d(P, 1).plus(d(M, 2)),
            d(N, 2).minus(d(P, 0)),
            d(M, 0).plus(d(N, 1)),
        )
        n = len(ay_exps) + len(as_exps)
        rows, dropped = [], 0
        for comp, flat_comp in zip((c1, c2, c3), flat):
            reachable = set(flat_comp.terms.keys())
            for exp, row in sorted(comp.terms.items()):
                if exp not in reachable:
                    dropped += 1
                    continue
                dense = np.zeros(n)
                for dof, c in row.items():
                    dense[dof] = c
                if np.any(dense != 0.0):
                    rows.append(dense)
        constraint = (np.asarray(rows) if rows else np.zeros((0, n)))
        return _ElementOperators(constraint=constraint,
                                 reachable_rows=len(rows),
                                 dropped_rows=dropped)

    @staticmethod
    def _least_defect_basis(constraint, dim):
        n = constraint.shape[1]
        if constraint.shape[0] == 0:
            return np.eye(n)[:, :dim], np.zeros(dim)
        _, s, vt = np.linalg.svd(constraint, full_matrices=True)
        s_full = np.concatenate((s, np.zeros(max(0, n - len(s)))))
        order = np.argsort(s_full)
        keep = order[:dim]
        return np.ascontiguousarray(vt[keep].T), s_full[keep]

    # -- evaluation ------------------------------------------------------
    def _metric(self, xi, zeta):
        h = np.zeros_like(np.asarray(zeta, dtype=float))
        for k, c in enumerate(self.curvature_poly_per_m):
            h = h + c * np.asarray(zeta, dtype=float) ** k
        return 1.0 + h * (np.asarray(xi, dtype=float) * self.half_width_m), h

    def _exponent_arrays(self):
        """Cached integer exponent arrays for vectorized evaluation."""
        cached = getattr(self, "_exponents_cache", None)
        if cached is None:
            cached = tuple(
                tuple(np.asarray(axis, dtype=np.int64) for axis in zip(*exps))
                for exps in (self.ay_exponents, self.as_exponents)
            )
            object.__setattr__(self, "_exponents_cache", cached)
        return cached

    @staticmethod
    def _monomials_and_derivative(coordinate, exponent):
        """``coordinate**exponent`` and its derivative, vectorized.

        Returns ``(value, derivative)`` of shapes ``(n_points, n_terms)``.
        The per-term Python loops these matrices replace were the measured
        bottleneck of every chain evaluation (A-RK rhs, fit design build).
        """
        maximum = int(exponent.max(initial=0))
        powers = np.empty((coordinate.size, maximum + 1))
        powers[:, 0] = 1.0
        for degree in range(1, maximum + 1):
            powers[:, degree] = powers[:, degree - 1] * coordinate
        value = powers[:, exponent]
        derivative = np.zeros_like(value)
        positive = exponent > 0
        derivative[:, positive] = (
            exponent[positive]
            * powers[:, np.maximum(exponent[positive] - 1, 0)]
        )
        return value, derivative

    def _basis_blocks(self):
        n_ay = len(self.ay_exponents)
        return self.basis[:n_ay], self.basis[n_ay:]

    def component_columns(self, xi, eta, zeta):
        """Value columns of (a_y, covariant a_s) at normalized points."""
        xi = np.asarray(xi, dtype=float).reshape(-1)
        eta = np.asarray(eta, dtype=float).reshape(-1)
        zeta = np.asarray(zeta, dtype=float).reshape(-1)
        (ay_i, ay_j, ay_k), (as_i, as_j, as_k) = self._exponent_arrays()
        basis_ay, basis_as = self._basis_blocks()
        ay_monomials = (
            self._monomials_and_derivative(xi, ay_i)[0]
            * self._monomials_and_derivative(eta, ay_j)[0]
            * self._monomials_and_derivative(zeta, ay_k)[0]
        )
        as_monomials = (
            self._monomials_and_derivative(xi, as_i)[0]
            * self._monomials_and_derivative(eta, as_j)[0]
            * self._monomials_and_derivative(zeta, as_k)[0]
        )
        return ay_monomials @ basis_ay, as_monomials @ basis_as

    def b_row_columns(self, xi, eta, zeta):
        """Columns of the POLYNOMIAL row quantities (g*Bx, g*By, Bs).

        The fit contract multiplies the sampled physical ``Bx, By`` by the
        metric ``g`` pointwise so both sides stay polynomial.
        """
        xi = np.asarray(xi, dtype=float).reshape(-1)
        eta = np.asarray(eta, dtype=float).reshape(-1)
        zeta = np.asarray(zeta, dtype=float).reshape(-1)
        ax, ay_scale, as_scale = self.scales
        (ay_i, ay_j, ay_k), (as_i, as_j, as_k) = self._exponent_arrays()
        basis_ay, basis_as = self._basis_blocks()
        ay_x, ay_dx = self._monomials_and_derivative(xi, ay_i)
        ay_y, _ = self._monomials_and_derivative(eta, ay_j)
        ay_z, ay_dz = self._monomials_and_derivative(zeta, ay_k)
        as_x, as_dx = self._monomials_and_derivative(xi, as_i)
        as_y, as_dy = self._monomials_and_derivative(eta, as_j)
        as_z, _ = self._monomials_and_derivative(zeta, as_k)
        gbx = (as_x * as_dy * as_z / ay_scale) @ basis_as \
            - (ay_x * ay_y * ay_dz / as_scale) @ basis_ay
        gby = -(as_dx * as_y * as_z / ax) @ basis_as
        bs = (ay_dx * ay_y * ay_z / ax) @ basis_ay
        return gbx, gby, bs

    # -- interface traces (graded L1 contract) ---------------------------
    def ay_trace_matrix(self, side):
        """Tangential-trace map: rows = face (i, j) monomials of ``a_y``."""
        zeta_value = 1.0 if side > 0 else -1.0
        pairs = sorted({(i, j) for (i, j, _) in self.ay_exponents})
        index = {p: r for r, p in enumerate(pairs)}
        T = np.zeros((len(pairs), self.dimension))
        for col in range(self.dimension):
            for idx, (i, j, k) in enumerate(self.ay_exponents):
                T[index[(i, j)], col] += self.basis[idx, col] * zeta_value**k
        return T

    def b_value_trace_matrix(self, side):
        """Midplane ``a_s`` face trace: rows = (i, 0) monomials (b_m values)."""
        zeta_value = 1.0 if side > 0 else -1.0
        n_ay = len(self.ay_exponents)
        rows = sorted({i for (i, j, _) in self.as_exponents if j == 0})
        index = {i: r for r, i in enumerate(rows)}
        T = np.zeros((len(rows), self.dimension))
        for col in range(self.dimension):
            for idx, (i, j, k) in enumerate(self.as_exponents):
                if j == 0:
                    T[index[i], col] += self.basis[n_ay + idx, col] \
                        * zeta_value**k
        return T

    # -- Lie / A-map consumption ----------------------------------------
    def transverse_coefficients(self, zeta, degree=5):
        """Physical ``(d+1, d+1)`` coefficient maps at fixed ``zeta``.

        Returns ``(Ay_map, As_map)``: matrices of shape ``(d+1, d+1,
        dimension)`` mapping modal coefficients to the physical transverse
        polynomial arrays ``coef[i, j] * x^i * y^j`` in T*m units.  ``As`` is
        the COVARIANT component, matching the Lie kernel's
        ``longitudinal_component='covariant'`` contract; the design-orbit
        gauge zeros are structural.
        """
        degree = int(degree)
        if degree < 1:
            raise ValueError("degree must be at least 1")
        zeta_value = float(zeta)
        ax, ay_scale, _ = self.scales
        n_ay = len(self.ay_exponents)
        Ay_map = np.zeros((degree + 1, degree + 1, self.dimension))
        As_map = np.zeros((degree + 1, degree + 1, self.dimension))
        for col in range(self.dimension):
            coeffs = self.basis[:, col]
            for idx, (i, j, k) in enumerate(self.ay_exponents):
                if i <= degree and j <= degree and i + j <= degree:
                    Ay_map[i, j, col] += coeffs[idx] * zeta_value**k \
                        / (ax**i * ay_scale**j)
            for idx, (i, j, k) in enumerate(self.as_exponents):
                if i <= degree and j <= degree and i + j <= degree:
                    As_map[i, j, col] += coeffs[n_ay + idx] * zeta_value**k \
                        / (ax**i * ay_scale**j)
        return Ay_map, As_map


@dataclass(frozen=True)
class CanonicalHCurlFit:
    """Fit result and its honesty certificates."""

    coefficients: np.ndarray
    maximum_residual_t: float
    field_scale_t: float
    maximum_interface_ay_jump: float
    maximum_interface_b_value_jump: float
    sample_count: int

    @property
    def relative_residual(self) -> float:
        return self.maximum_residual_t / self.field_scale_t


class CanonicalHCurlChain:
    """An open chain of CanonicalHCurl elements along the arc length.

    ``s_breaks_m`` are the element boundaries (monotone, ``E+1`` values for
    ``E`` elements).  All elements share the transverse half-apertures.
    ``curvature_per_m`` is a callable ``h(s)`` (physical design-orbit
    curvature with the metric convention ``g = 1 + htilde x``); it is
    sampled per element and represented by a degree-``order_s`` polynomial
    (curve order = ``p_s``).

    Chained production requires ``order_s >= 2``: the measured L1 interface
    contract (a_y trace + b_m values) collapses ``order_s = 1`` chains to a
    globally linear ``b_m``.  A single element accepts any ``order_s``.

    ``periodic=True`` closes the chain into a ring: the L1 contract also
    couples the last element's exit face to the first element's entrance
    face, and the fixed-dimension target becomes the periodic spline count
    ``p_x * E`` per the maximal-smoothness spline dimension.  The ring's
    cohomology obstruction -- the design orbit links a flux, so the
    on-orbit circulation of ``a_s`` cannot vanish although the centre gauge
    forces it to -- is carried by ONE explicit global DOF: a constant
    covariant ``a_s`` (a zero-field harmonic 1-form) set through
    :meth:`set_ring_circulation` from the gauge-invariant loop integral.
    It never enters the Lie/A-RK dynamics (a constant in ``H`` exerts no
    force); it only restores the flux observable.
    """

    def __init__(self, s_breaks_m, half_width_m, half_height_m, *,
                 order_x, order_s, curvature_per_m=None, periodic=False):
        breaks = np.asarray(s_breaks_m, dtype=float).reshape(-1)
        if breaks.size < 2 or not np.all(np.diff(breaks) > 0.0):
            raise ValueError("s_breaks_m must be strictly increasing, >= 2")
        self.periodic = bool(periodic)
        if int(order_s) < 2 and (breaks.size > 2 or self.periodic):
            raise ValueError(
                "chained CanonicalHCurl requires order_s >= 2 (the L1 "
                "interface contract degenerates order_s=1 chains)")
        if self.periodic and breaks.size < 3:
            raise ValueError("a periodic chain needs at least 2 elements")
        self.s_breaks_m = breaks
        self.elements = []
        for e in range(breaks.size - 1):
            s0, s1 = float(breaks[e]), float(breaks[e + 1])
            half_length = 0.5 * (s1 - s0)
            if curvature_per_m is None:
                poly = (0.0,)
            else:
                nodes = np.linspace(-1.0, 1.0, int(order_s) + 1)
                values = np.asarray([
                    float(curvature_per_m(0.5 * (s0 + s1)
                                          + half_length * node))
                    for node in nodes
                ])
                poly = tuple(np.polynomial.polynomial.polyfit(
                    nodes, values, int(order_s)))
            self.elements.append(CanonicalHCurlElement(
                order_x=int(order_x), order_s=int(order_s),
                half_width_m=float(half_width_m),
                half_height_m=float(half_height_m),
                half_length_m=half_length,
                curvature_poly_per_m=poly,
            ))
        self.interface_defects = np.zeros(0)
        self.interface_defect_scale = 1.0
        self._reduced = self._interface_null_space()
        self._fit = None
        self.ring_circulation_t_m = 0.0

    # -- chain structure -------------------------------------------------
    @property
    def element_count(self) -> int:
        return len(self.elements)

    @property
    def total_dimension(self) -> int:
        return sum(el.dimension for el in self.elements)

    @property
    def chain_dimension(self) -> int:
        return self._reduced.shape[1]

    def _interface_pairs(self):
        pairs = [(e, e + 1) for e in range(self.element_count - 1)]
        if self.periodic:
            pairs.append((self.element_count - 1, 0))
        return pairs

    def _interface_constraint_rows(self):
        offsets = np.concatenate(([0], np.cumsum(
            [el.dimension for el in self.elements])))
        total = offsets[-1]
        rows = []
        for e_left, e_right in self._interface_pairs():
            left, right = self.elements[e_left], self.elements[e_right]
            for maker in ("ay_trace_matrix", "b_value_trace_matrix"):
                Tp = getattr(left, maker)(+1)
                Tm = getattr(right, maker)(-1)
                if Tp.shape[0] != Tm.shape[0]:
                    raise AssertionError("trace row layouts must match")
                block = np.zeros((Tp.shape[0], total))
                block[:, offsets[e_left]:offsets[e_left + 1]] = Tp
                block[:, offsets[e_right]:offsets[e_right + 1]] = -Tm
                rows.append(block)
        return (np.vstack(rows) if rows else np.zeros((0, total)))

    def _interface_null_space(self):
        """Fixed-dimension least-defect reduction of the interface contract.

        With per-element curvature polynomials the neighbours' generated
        traces lose the flat-case alignment, so the STRICT joint kernel
        collapses (the same truncation phenomenon as the curved element
        kernel).  The chain is therefore DEFINED by the measured spline law:
        the ``p_x * (E + p_s)`` smallest singular directions of the L1
        contract (exact strict kernel at h=0; interface jumps O(h * defect)
        otherwise, reported by the fit certificate).
        """
        constraint = self._interface_constraint_rows()
        total = constraint.shape[1]
        order_x = self.elements[0].order_x
        order_s = self.elements[0].order_s
        if self.periodic:
            # Maximal-smoothness periodic spline dimension per multipole.
            target = order_x * self.element_count
        else:
            target = order_x * (self.element_count + order_s)
        if constraint.shape[0] == 0:
            return np.eye(total)[:, :target]
        _, s, vt = np.linalg.svd(constraint, full_matrices=True)
        s_full = np.concatenate((s, np.zeros(max(0, total - len(s)))))
        order = np.argsort(s_full)
        keep = order[:target]
        self.interface_defects = s_full[keep]
        self.interface_defect_scale = float(s[0]) if s.size else 1.0
        return np.ascontiguousarray(vt[keep].T)

    def _locate(self, s_m):
        s = np.asarray(s_m, dtype=float).reshape(-1)
        if np.any(s < self.s_breaks_m[0] - 1e-12) or np.any(
                s > self.s_breaks_m[-1] + 1e-12):
            raise ValueError("s outside the chain range")
        index = np.clip(np.searchsorted(self.s_breaks_m, s, side="right") - 1,
                        0, self.element_count - 1)
        s0 = self.s_breaks_m[index]
        s1 = self.s_breaks_m[index + 1]
        zeta = np.clip((2.0 * s - s0 - s1) / (s1 - s0), -1.0, 1.0)
        return index, zeta

    # -- projection ------------------------------------------------------
    def fit_frame_samples(self, x_m, y_m, s_m, bx_t, by_t, bs_t):
        """Full-volume B fit in frame components (the production projection).

        ``bx, by, bs`` are the physical frame components of ``B`` sampled at
        ``(x, y, s)``; the rows fitted are the polynomial quantities
        ``(g bx, g by, bs)``.  Midplane-only clouds are rejected: the
        vacuum-constrained fit "solves" from midplane data while
        under-reporting the true aperture error by an order of magnitude.
        """
        x = np.asarray(x_m, dtype=float).reshape(-1)
        y = np.asarray(y_m, dtype=float).reshape(-1)
        s = np.asarray(s_m, dtype=float).reshape(-1)
        bx = np.asarray(bx_t, dtype=float).reshape(-1)
        by = np.asarray(by_t, dtype=float).reshape(-1)
        bs = np.asarray(bs_t, dtype=float).reshape(-1)
        if not (x.size == y.size == s.size == bx.size == by.size == bs.size):
            raise ValueError("sample arrays must share one length")
        if x.size < 4 * self.chain_dimension:
            raise ValueError("need at least 4 samples per chain DOF")
        half_height = self.elements[0].half_height_m
        if float(np.max(np.abs(y), initial=0.0)) < 0.2 * half_height:
            raise ValueError(
                "midplane-only sample cloud rejected: off-plane B data is "
                "required for an honest aperture certificate")

        index, zeta = self._locate(s)
        counts = np.bincount(index, minlength=self.element_count)
        starved = [int(e) for e in range(self.element_count)
                   if counts[e] < self.elements[e].dimension]
        if starved:
            raise ValueError(
                "sample starvation: elements "
                f"{starved} hold fewer samples than their dimension "
                f"({self.elements[0].dimension}); distribute the cloud per "
                "element (thin graded elements need their own stations)")
        offsets = np.concatenate(([0], np.cumsum(
            [el.dimension for el in self.elements])))
        design = np.zeros((3 * x.size, offsets[-1]))
        rhs = np.zeros(3 * x.size)
        n = x.size
        for e, element in enumerate(self.elements):
            mask = index == e
            if not np.any(mask):
                continue
            xi = x[mask] / element.half_width_m
            eta = y[mask] / element.half_height_m
            gbx, gby, bs_cols = element.b_row_columns(xi, eta, zeta[mask])
            g, _ = element._metric(xi, zeta[mask])
            where = np.flatnonzero(mask)
            sl = slice(offsets[e], offsets[e + 1])
            design[where, sl] = gbx
            design[n + where, sl] = gby
            design[2 * n + where, sl] = bs_cols
            rhs[where] = g * bx[mask]
            rhs[n + where] = g * by[mask]
            rhs[2 * n + where] = bs[mask]
        reduced = design @ self._reduced
        solution, *_ = np.linalg.lstsq(reduced, rhs, rcond=None)
        coefficients = self._reduced @ solution
        residual = reduced @ solution - rhs
        per_point = np.sqrt(residual[:n]**2 + residual[n:2 * n]**2
                            + residual[2 * n:]**2)
        fit = CanonicalHCurlFit(
            coefficients=coefficients,
            maximum_residual_t=float(np.max(per_point)),
            field_scale_t=float(np.max(np.sqrt(bx**2 + by**2 + bs**2))),
            maximum_interface_ay_jump=self._interface_jump(
                coefficients, "ay_trace_matrix"),
            maximum_interface_b_value_jump=self._interface_jump(
                coefficients, "b_value_trace_matrix"),
            sample_count=int(n),
        )
        self._fit = fit
        return fit

    def _interface_jump(self, coefficients, maker):
        offsets = np.concatenate(([0], np.cumsum(
            [el.dimension for el in self.elements])))
        worst = 0.0
        for e_left, e_right in self._interface_pairs():
            left = getattr(self.elements[e_left], maker)(+1) @ \
                coefficients[offsets[e_left]:offsets[e_left + 1]]
            right = getattr(self.elements[e_right], maker)(-1) @ \
                coefficients[offsets[e_right]:offsets[e_right + 1]]
            if left.size:
                worst = max(worst, float(np.max(np.abs(left - right))))
        return worst

    def set_ring_circulation(self, circulation_t_m2):
        """Set the ring's cohomology DOF from the gauge-invariant loop flux.

        ``circulation_t_m2`` is the design-orbit loop integral of the source
        vector potential (= the flux linked by the closed orbit, in T*m^2);
        the stored constant covariant ``a_s`` is that value divided by the
        circumference.  Zero-field mode: it never affects B, the fit, or
        the tracking dynamics -- only the ``loop-integral A`` observable.
        """
        if not self.periodic:
            raise ValueError("ring circulation applies to periodic chains")
        length = float(self.s_breaks_m[-1] - self.s_breaks_m[0])
        self.ring_circulation_t_m = float(circulation_t_m2) / length

    def _require_fit(self):
        if self._fit is None:
            raise RuntimeError("call fit_frame_samples before evaluation")
        return self._fit.coefficients

    # -- evaluation ------------------------------------------------------
    def vector_potential_frame(self, x_m, y_m, s_m):
        """``(a_x, a_y, a_s_covariant)`` at frame points, in T*m."""
        coefficients = self._require_fit()
        x = np.asarray(x_m, dtype=float).reshape(-1)
        y = np.asarray(y_m, dtype=float).reshape(-1)
        s = np.asarray(s_m, dtype=float).reshape(-1)
        index, zeta = self._locate(s)
        offsets = np.concatenate(([0], np.cumsum(
            [el.dimension for el in self.elements])))
        out = np.zeros((x.size, 3))
        for e, element in enumerate(self.elements):
            mask = index == e
            if not np.any(mask):
                continue
            ay_cols, as_cols = element.component_columns(
                x[mask] / element.half_width_m,
                y[mask] / element.half_height_m, zeta[mask])
            c = coefficients[offsets[e]:offsets[e + 1]]
            out[mask, 1] = ay_cols @ c
            out[mask, 2] = as_cols @ c
        out[:, 2] += self.ring_circulation_t_m
        return out

    def vector_potential_and_gradient_frame(self, x_m, y_m, s_m):
        """``(a, da/d(x,y))`` for the canonical A-map RK.

        Returns ``(a, gradient)`` with ``a`` of shape ``(n, 3)`` as in
        :meth:`vector_potential_frame` (covariant ``a_s``) and ``gradient``
        of shape ``(n, 3, 2)`` holding ``d a_i / d(x, y)`` in T*m/m --
        exactly the layout of
        ``accelerator_lie_topopt.canonical_vector_potential_hamiltonian_rhs``.
        """
        coefficients = self._require_fit()
        x = np.asarray(x_m, dtype=float).reshape(-1)
        y = np.asarray(y_m, dtype=float).reshape(-1)
        s = np.asarray(s_m, dtype=float).reshape(-1)
        index, zeta = self._locate(s)
        offsets = np.concatenate(([0], np.cumsum(
            [el.dimension for el in self.elements])))
        a = np.zeros((x.size, 3))
        gradient = np.zeros((x.size, 3, 2))
        for e, element in enumerate(self.elements):
            mask = index == e
            if not np.any(mask):
                continue
            xi = x[mask] / element.half_width_m
            eta = y[mask] / element.half_height_m
            zl = zeta[mask]
            c = coefficients[offsets[e]:offsets[e + 1]]
            ax_scale, ay_scale, _ = element.scales
            basis_ay, basis_as = element._basis_blocks()
            modal_ay = basis_ay @ c
            modal_as = basis_as @ c
            (ay_i, ay_j, ay_k), (as_i, as_j, as_k) = \
                element._exponent_arrays()
            ay_x, ay_dx = element._monomials_and_derivative(xi, ay_i)
            ay_y, ay_dy = element._monomials_and_derivative(eta, ay_j)
            ay_z, _ = element._monomials_and_derivative(zl, ay_k)
            as_x, as_dx = element._monomials_and_derivative(xi, as_i)
            as_y, as_dy = element._monomials_and_derivative(eta, as_j)
            as_z, _ = element._monomials_and_derivative(zl, as_k)
            a[mask, 1] = (ay_x * ay_y * ay_z) @ modal_ay
            a[mask, 2] = (as_x * as_y * as_z) @ modal_as
            gradient[mask, 1, 0] = (ay_dx * ay_y * ay_z) @ modal_ay / ax_scale
            gradient[mask, 1, 1] = (ay_x * ay_dy * ay_z) @ modal_ay / ay_scale
            gradient[mask, 2, 0] = (as_dx * as_y * as_z) @ modal_as / ax_scale
            gradient[mask, 2, 1] = (as_x * as_dy * as_z) @ modal_as / ay_scale
        a[:, 2] += self.ring_circulation_t_m
        return a, gradient

    def magnetic_flux_density_frame(self, x_m, y_m, s_m):
        """Physical frame ``(Bx, By, Bs)`` of the fitted field, in T."""
        coefficients = self._require_fit()
        x = np.asarray(x_m, dtype=float).reshape(-1)
        y = np.asarray(y_m, dtype=float).reshape(-1)
        s = np.asarray(s_m, dtype=float).reshape(-1)
        index, zeta = self._locate(s)
        offsets = np.concatenate(([0], np.cumsum(
            [el.dimension for el in self.elements])))
        out = np.zeros((x.size, 3))
        for e, element in enumerate(self.elements):
            mask = index == e
            if not np.any(mask):
                continue
            xi = x[mask] / element.half_width_m
            gbx, gby, bs_cols = element.b_row_columns(
                xi, y[mask] / element.half_height_m, zeta[mask])
            g, _ = element._metric(xi, zeta[mask])
            c = coefficients[offsets[e]:offsets[e + 1]]
            out[mask, 0] = (gbx @ c) / g
            out[mask, 1] = (gby @ c) / g
            out[mask, 2] = bs_cols @ c
        return out

    def lie_element_spoly_arrays(self, degree=5):
        """Per-ELEMENT s-polynomial transverse arrays for the Lie kernel.

        Returns ``(Ay, As_covariant, lengths, curvature_polys)`` with
        shapes ``(E, p_s+1, d+1, d+1)`` (entry ``[e, k, i, j]`` multiplies
        ``zeta**k * x**i * y**j`` on element ``e``'s normalized
        ``zeta in [-1, 1]``), ``(E,)`` and ``(E, p_s+1)``.  One Lie
        segment per chain element then consumes the ``p_s`` dependence
        DIRECTLY through the kernel's nonautonomous stage-jet RK4 --
        the in-segment alternative to midpoint staging.
        """
        coefficients = self._require_fit()
        degree = int(degree)
        if degree < 1:
            raise ValueError("degree must be at least 1")
        order_s = self.elements[0].order_s
        count = self.element_count
        offsets = np.concatenate(([0], np.cumsum(
            [el.dimension for el in self.elements])))
        Ay = np.zeros((count, order_s + 1, degree + 1, degree + 1))
        As = np.zeros((count, order_s + 1, degree + 1, degree + 1))
        curvature_polys = np.zeros((count, order_s + 1))
        for e, element in enumerate(self.elements):
            c = coefficients[offsets[e]:offsets[e + 1]]
            modal = element.basis @ c
            n_ay = len(element.ay_exponents)
            ax, ay_scale, _ = element.scales
            for idx, (i, j, k) in enumerate(element.ay_exponents):
                if i <= degree and j <= degree and i + j <= degree:
                    Ay[e, k, i, j] += modal[idx] / (ax**i * ay_scale**j)
            for idx, (i, j, k) in enumerate(element.as_exponents):
                if i <= degree and j <= degree and i + j <= degree:
                    As[e, k, i, j] += modal[n_ay + idx] \
                        / (ax**i * ay_scale**j)
            poly = element.curvature_poly_per_m
            curvature_polys[e, :len(poly)] = poly
        lengths = np.diff(self.s_breaks_m)
        return Ay, As, lengths, curvature_polys

    def lie_segment_arrays(self, segment_count, degree=5):
        """Per-segment transverse coefficient arrays for the Lie kernel.

        Splits the chain into ``segment_count`` equal arc-length segments,
        evaluates the transverse polynomials at each segment midpoint
        (s-order staging: midpoint evaluation now, in-segment s-polynomial
        flow later), and returns ``(Ay, As_covariant, lengths, curvatures)``
        with shapes ``(n, d+1, d+1)``, ``(n,)``, ``(n,)`` matching
        ``longitudinal_component='covariant'``.
        """
        coefficients = self._require_fit()
        count = int(segment_count)
        if count < 1:
            raise ValueError("segment_count must be positive")
        edges = np.linspace(self.s_breaks_m[0], self.s_breaks_m[-1],
                            count + 1)
        mids = 0.5 * (edges[:-1] + edges[1:])
        lengths = np.diff(edges)
        index, zeta = self._locate(mids)
        offsets = np.concatenate(([0], np.cumsum(
            [el.dimension for el in self.elements])))
        Ay = np.zeros((count, degree + 1, degree + 1))
        As = np.zeros((count, degree + 1, degree + 1))
        curvatures = np.zeros(count)
        for seg in range(count):
            element = self.elements[int(index[seg])]
            Ay_map, As_map = element.transverse_coefficients(
                zeta[seg], degree=degree)
            c = coefficients[offsets[index[seg]]:offsets[index[seg] + 1]]
            Ay[seg] = Ay_map @ c
            As[seg] = As_map @ c
            h = 0.0
            for k, ck in enumerate(element.curvature_poly_per_m):
                h += ck * float(zeta[seg]) ** k
            curvatures[seg] = h
        return Ay, As, lengths, curvatures
