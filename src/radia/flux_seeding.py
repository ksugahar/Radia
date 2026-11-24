"""Equal-flux magnetic field-line seeds on a parameterized transverse surface.

A one-to-one chart F(xi, eta) maps a coordinate rectangle onto the seed
surface. The oriented normal flux density pulled back to the chart is

    g(xi, eta) = B(F(xi, eta)) . (F_xi x F_eta)    [T m^2 per chart area]

With rho(eta) = int g dxi and M = int rho deta, the normalized marginal
V(eta) = int_0^eta rho / M and the conditional
U(xi | eta) = int_0^xi g(., eta) / rho(eta) define flux coordinates (u, v)
whose Jacobian is g / M. Each cell of a uniform (u, v) grid therefore carries
the flux M du dv. A seed sits at the image of each cell centre: V is inverted
first, then U is inverted at that exact eta.

Sampling model. The marginal uses a table of g on xi_nodes x eta_nodes, with
g piecewise linear in xi on each row and rho piecewise linear in eta, so V is
piecewise quadratic and is inverted in closed form. The conditional is
re-sampled at every seed row eta in one vectorized call of the flux density
on xi_nodes, never read from the nearest table row, and is inverted the same
way. The model is exact when g is bilinear in (xi, eta).

Supported geometry. Any one-to-one chart whose pulled-back density has one
sign in the interior of the coordinate rectangle. planar_patch builds the
regular planar parallelogram. The density may vanish on the rectangle edges
(a coordinate singularity at an endpoint) but not at an interior sample.
Surfaces with holes, periodic seams or a sign-changing normal flux must be
split by the caller. Seed generation is separate from trajectory
integration; no field-line closure is implied.
"""

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

METHOD = "two-stage inverse CDF, piecewise-linear table, exact-eta conditional"


@dataclass(frozen=True)
class EqualFluxSeeds:
    """Seeds of a uniform n_u x n_v flux-coordinate grid.

    ids run row-major: id = j * n_u + i for flux row j and column i.
    chart_coords are (xi, eta); flux_coords are the cell centres (u, v) in
    [0, 1]. Fluxes in Wb carry the chart orientation sign.
    """

    ids: NDArray[np.int64]
    points_m: NDArray[np.float64]
    chart_coords: NDArray[np.float64]
    flux_coords: NDArray[np.float64]
    n_u: int
    n_v: int
    orientation: int
    total_flux_Wb: float
    target_cell_flux_Wb: float
    diagnostics: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)


@dataclass(frozen=True)
class PlanarPatch:
    """Regular planar parallelogram F(xi, eta) = origin + xi e_xi + eta e_eta."""

    origin_m: NDArray[np.float64]
    e_xi_m: NDArray[np.float64]
    e_eta_m: NDArray[np.float64]

    @property
    def area_vector_m2(self) -> NDArray[np.float64]:
        return np.cross(self.e_xi_m, self.e_eta_m)

    def chart(self, xi, eta):
        xi, eta = np.broadcast_arrays(np.asarray(xi, float), np.asarray(eta, float))
        return (self.origin_m + xi[..., None] * self.e_xi_m
                + eta[..., None] * self.e_eta_m)

    def flux_density(self, field_T: Callable) -> Callable:
        """Return g(xi, eta) = B . (e_xi x e_eta) for field_T(points) -> (..., 3)."""
        normal = self.area_vector_m2

        def g(xi, eta):
            b = np.asarray(field_T(self.chart(xi, eta)), dtype=float)
            return b @ normal

        return g


def planar_patch(origin_m: ArrayLike, e_xi_m: ArrayLike, e_eta_m: ArrayLike) -> PlanarPatch:
    """Build a planar chart on the unit square; reject a degenerate one."""
    vectors = [np.asarray(v, dtype=float) for v in (origin_m, e_xi_m, e_eta_m)]
    if any(v.shape != (3,) or not np.all(np.isfinite(v)) for v in vectors):
        raise ValueError("origin_m, e_xi_m and e_eta_m must be finite 3-vectors")
    patch = PlanarPatch(*vectors)
    scale = np.linalg.norm(vectors[1]) * np.linalg.norm(vectors[2])
    if not scale > 0 or np.linalg.norm(patch.area_vector_m2) <= 1e-12 * scale:
        raise ValueError("planar chart is not invertible: e_xi and e_eta are parallel or zero")
    return patch


def _axis(values: ArrayLike, name: str) -> NDArray[np.float64]:
    axis = np.asarray(values, dtype=float)
    if (axis.ndim != 1 or len(axis) < 2 or not np.all(np.isfinite(axis))
            or not np.all(np.diff(axis) > 0)):
        raise ValueError(f"{name} must be a finite, strictly increasing 1-D grid")
    return axis


def _check_sign(values: NDArray, where: str) -> int:
    if not np.all(np.isfinite(values)):
        raise ValueError(f"flux density must be finite ({where})")
    positive, negative = np.any(values > 0), np.any(values < 0)
    if positive and negative:
        raise ValueError(f"flux density changes sign ({where}); split the surface")
    if not (positive or negative):
        raise ValueError(f"flux density is identically zero ({where})")
    return 1 if positive else -1


def _cumulative(x: NDArray, density: NDArray) -> NDArray:
    """Trapezoidal cumulative integral along the last axis, starting at 0."""
    steps = 0.5 * (density[..., 1:] + density[..., :-1]) * np.diff(x)
    return np.concatenate((np.zeros(density.shape[:-1] + (1,)), np.cumsum(steps, axis=-1)),
                          axis=-1)


def _invert(x: NDArray, density: NDArray, cumulative: NDArray, targets: NDArray) -> NDArray:
    """Invert the piecewise-quadratic cumulative of a piecewise-linear density.

    density and cumulative have shape (rows, len(x)); targets (rows, k) lie in
    [0, cumulative[:, -1]]. Uses the cancellation-free quadratic root.
    """
    rows = np.arange(density.shape[0])[:, None]
    seg = np.empty(targets.shape, dtype=np.int64)
    for r in range(density.shape[0]):
        seg[r] = np.searchsorted(cumulative[r], targets[r], side="right") - 1
    seg = np.clip(seg, 0, len(x) - 2)
    h = np.diff(x)[seg]
    d0, d1 = density[rows, seg], density[rows, seg + 1]
    rest = np.maximum(targets - cumulative[rows, seg], 0.0)
    a = (d1 - d0) / (2 * h)
    disc = np.sqrt(np.maximum(d0 * d0 + 4 * a * rest, 0.0))
    t = 2 * rest / (d0 + disc)
    return x[seg] + np.minimum(t, h)


def _evaluate(x: NDArray, density: NDArray, cumulative: NDArray, at: NDArray) -> NDArray:
    """Evaluate the same piecewise-quadratic cumulative at points `at`."""
    rows = np.arange(density.shape[0])[:, None]
    seg = np.empty(at.shape, dtype=np.int64)
    for r in range(density.shape[0]):
        seg[r] = np.searchsorted(x, at[r], side="right") - 1
    seg = np.clip(seg, 0, len(x) - 2)
    h = np.diff(x)[seg]
    t = at - x[seg]
    d0, d1 = density[rows, seg], density[rows, seg + 1]
    return cumulative[rows, seg] + d0 * t + (d1 - d0) * t * t / (2 * h)


def _row_density(flux_density: Callable, xi: NDArray, eta: NDArray) -> NDArray:
    g = np.asarray(flux_density(xi[None, :], eta[:, None]), dtype=float)
    return np.broadcast_to(g, (len(eta), len(xi))).copy()


def _check_rows(g: NDArray, orientation: int, where: str, *, edge_rows: bool = False) -> NDArray:
    """Return orientation * g after checking finiteness and interior sign.

    With edge_rows the first and last rows lie on the rectangle edge and may
    vanish there; only interior rows must be nonzero at interior xi.
    """
    if not np.all(np.isfinite(g)):
        raise ValueError(f"flux density must be finite ({where})")
    oriented = orientation * g
    if np.any(oriented < 0):
        raise ValueError(f"flux density changes sign ({where}); split the surface")
    interior = oriented[1:-1, 1:-1] if edge_rows else oriented[:, 1:-1]
    if np.any(interior == 0):
        raise ValueError(f"flux density vanishes at an interior xi sample ({where})")
    return oriented


def equal_flux_seeds(
    flux_density: Callable,
    *,
    xi_nodes: ArrayLike,
    eta_nodes: ArrayLike,
    n_u: int,
    n_v: int,
    chart: Callable,
) -> EqualFluxSeeds:
    """Place n_u * n_v equal-flux seeds by the two-stage inverse.

    flux_density(xi, eta) broadcasts and returns g in Wb per unit chart area
    (T m^2 when xi, eta are dimensionless); chart(xi, eta) returns physical
    points in metres with shape (..., 3). The table rows xi_nodes x eta_nodes
    must resolve g; tables are a sampling choice, not a default accuracy.

    Raises ValueError for nonfinite samples, a sign change, an interior zero,
    a zero total flux, a flat marginal or conditional segment, or a chart that
    returns nonfinite points.
    """
    xi, eta = _axis(xi_nodes, "xi_nodes"), _axis(eta_nodes, "eta_nodes")
    n_u, n_v = int(n_u), int(n_v)
    if n_u < 1 or n_v < 1:
        raise ValueError("n_u and n_v must be positive")

    table = _row_density(flux_density, xi, eta)
    orientation = _check_sign(table, "table")
    table = _check_rows(table, orientation, "table", edge_rows=True)

    rho = _cumulative(xi, table)[:, -1]
    if np.any(rho[1:-1] <= 0):
        raise ValueError("an interior eta row carries no flux")
    marginal = _cumulative(eta, rho[None, :])
    total = float(marginal[0, -1])
    if not np.isfinite(total) or total <= 0:
        raise ValueError("total flux must be finite and nonzero")
    if np.any(np.diff(marginal[0]) <= 0):
        raise ValueError("marginal flux is not strictly increasing; a table segment is flat")

    v_c = (np.arange(n_v) + 0.5) / n_v
    u_c = (np.arange(n_u) + 0.5) / n_u
    eta_c = _invert(eta, rho[None, :], marginal, total * v_c[None, :])[0]

    rows = _check_rows(_row_density(flux_density, xi, eta_c), orientation, "seed rows")
    row_cum = _cumulative(xi, rows)
    row_rho = row_cum[:, -1]
    if np.any(row_rho <= 0) or np.any(np.diff(row_cum, axis=1) <= 0):
        raise ValueError("a seed row has a flat conditional flux segment")
    targets = row_rho[:, None] * u_c[None, :]
    xi_c = _invert(xi, rows, row_cum, targets)

    v_check = _evaluate(eta, rho[None, :], marginal, eta_c[None, :])[0] / total
    u_check = _evaluate(xi, rows, row_cum, xi_c) / row_rho[:, None]

    chart_coords = np.column_stack((xi_c.ravel(), np.repeat(eta_c, n_u)))
    points = np.asarray(chart(chart_coords[:, 0], chart_coords[:, 1]), dtype=float)
    if points.shape != (n_u * n_v, 3) or not np.all(np.isfinite(points)):
        raise ValueError("chart must return finite 3-D points for every seed")
    flux_coords = np.column_stack((np.tile(u_c, n_v), np.repeat(v_c, n_u)))

    signed_total = orientation * total
    return EqualFluxSeeds(
        ids=np.arange(n_u * n_v, dtype=np.int64),
        points_m=points,
        chart_coords=chart_coords,
        flux_coords=flux_coords,
        n_u=n_u,
        n_v=n_v,
        orientation=orientation,
        total_flux_Wb=signed_total,
        target_cell_flux_Wb=signed_total / (n_u * n_v),
        diagnostics={
            "forward_residual_v_max": float(np.max(np.abs(v_check - v_c))),
            "forward_residual_u_max": float(np.max(np.abs(u_check - u_c[None, :]))),
            "table_shape": [len(eta), len(xi)],
            "seed_rows_resampled": int(n_v),
        },
        provenance={"method": METHOD, "module": __name__},
    )


def _gauss(order: int):
    nodes, weights = np.polynomial.legendre.leggauss(int(order))
    return 0.5 * (nodes + 1), 0.5 * weights


def cell_flux_quadrature(
    flux_density: Callable,
    seeds: EqualFluxSeeds,
    *,
    xi_nodes: ArrayLike,
    eta_nodes: ArrayLike,
    order: int = 8,
) -> dict:
    """Integrate the original density over each placement cell by Gauss rules.

    Cells are the model's own boundaries: eta edges from the table marginal,
    curved xi edges from the conditional re-sampled at each quadrature eta.
    The integrand is the caller's flux_density, never the table, so the
    result measures the equal-flux deviation of the placement. The total is
    an independent composite Gauss integral over the whole rectangle.
    Returns fluxes in Wb with the chart orientation sign.
    """
    xi, eta = _axis(xi_nodes, "xi_nodes"), _axis(eta_nodes, "eta_nodes")
    n_u, n_v, sign = seeds.n_u, seeds.n_v, seeds.orientation
    table = _check_rows(_row_density(flux_density, xi, eta), sign, "table", edge_rows=True)
    rho = _cumulative(xi, table)[:, -1]
    marginal = _cumulative(eta, rho[None, :])
    edges_v = np.linspace(0, 1, n_v + 1)
    eta_edges = _invert(eta, rho[None, :], marginal, marginal[0, -1] * edges_v[None, :])[0]
    eta_edges[0], eta_edges[-1] = eta[0], eta[-1]
    t, w = _gauss(order)

    cell = np.zeros((n_v, n_u))
    edges_u = np.linspace(0, 1, n_u + 1)
    for j in range(n_v):
        d_eta = eta_edges[j + 1] - eta_edges[j]
        eta_q = eta_edges[j] + d_eta * t
        rows = _check_rows(_row_density(flux_density, xi, eta_q), sign, "quadrature rows")
        cum = _cumulative(xi, rows)
        xi_edges = _invert(xi, rows, cum, cum[:, -1:] * edges_u[None, :])
        xi_edges[:, 0], xi_edges[:, -1] = xi[0], xi[-1]
        width = np.diff(xi_edges, axis=1)
        xi_q = xi_edges[:, :-1, None] + width[..., None] * t
        g = np.asarray(flux_density(xi_q, eta_q[:, None, None]), dtype=float)
        inner = sign * np.sum(g * w, axis=-1) * width
        cell[j] = d_eta * np.sum(w[:, None] * inner, axis=0)

    xi_t = (xi[:-1, None] + np.diff(xi)[:, None] * t).ravel()
    xi_w = (np.diff(xi)[:, None] * w).ravel()
    eta_t = (eta[:-1, None] + np.diff(eta)[:, None] * t).ravel()
    eta_w = (np.diff(eta)[:, None] * w).ravel()
    total = sign * float(eta_w @ (np.asarray(flux_density(xi_t[None, :], eta_t[:, None]),
                                             dtype=float) @ xi_w))
    target = total / (n_u * n_v)
    return {
        "cell_flux_Wb": sign * cell.ravel(),
        "independent_total_flux_Wb": sign * total,
        "max_relative_cell_deviation": float(np.max(np.abs(cell / target - 1))),
        "coverage": float(cell.sum() / total),
        "order": int(order),
    }
