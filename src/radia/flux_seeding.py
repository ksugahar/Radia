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

Supported geometry. Any one-to-one chart. planar_patch builds the planar
parallelogram, annular_patch a planar annulus or disk sector in polar
coordinates, cylindrical_patch a circular cylinder sector. Each chart also
gives its metric (physical length per unit chart coordinate) and locate(),
the signed offset from the surface and the chart coordinates of a point.

Two layouts. equal_flux_seeds is the n_u x n_v grid above; its density must
keep one sign in the interior (it may vanish on the rectangle edges, a
coordinate singularity). balanced_equal_flux_seeds places any number of
seeds by recursive equal-flux bisection: every split cuts the longer
physical side of a chart rectangle so that each part carries flux in
proportion to its seed count (a one-dimensional inverse of the marginal),
so cells stay near square where the grid gives thin strips. It seeds one
sign of a sign-changing surface: in a field with closed lines a surface is
usually crossed once in each direction, so seeding both signs would draw
each flux tube twice. Seed generation is separate from trajectory
integration; no field-line closure is implied.
"""

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

METHOD = "two-stage inverse CDF, piecewise-linear table, exact-eta conditional"
BALANCED_METHOD = ("recursive equal-flux bisection of the longer physical side, "
                   "bilinear table, one-sided flux")
SIGNS = ("auto", "positive", "negative")


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
class BalancedFluxSeeds:
    """Seeds of the recursive equal-flux bisection, one per chart cell.

    cells are chart rectangles (xi0, xi1, eta0, eta1); ids follow the
    bisection order. Only the selected sign is seeded: total_flux_Wb is the
    signed flux of that sign, positive_flux_Wb and negative_flux_Wb (<= 0)
    are both parts of the surface flux.
    """

    ids: NDArray[np.int64]
    points_m: NDArray[np.float64]
    chart_coords: NDArray[np.float64]
    cells: NDArray[np.float64]
    n_lines: int
    orientation: int
    total_flux_Wb: float
    target_cell_flux_Wb: float
    positive_flux_Wb: float
    negative_flux_Wb: float
    diagnostics: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)


def _grid(xi, eta):
    return np.broadcast_arrays(np.asarray(xi, float), np.asarray(eta, float))


class _Chart:
    """Common chart interface: chart, area_vector, metric, locate, flux_density."""

    def monotone_points(self, points_m: NDArray) -> tuple[NDArray, NDArray]:
        """Insert points so the signed offset is monotone along every segment.

        Returns the points and the new index of each original vertex. The
        offset of a planar chart is linear along a segment: nothing to add.
        """
        return points_m, np.arange(len(points_m))

    def flux_density(self, field_T: Callable) -> Callable:
        """Return g(xi, eta) = B . (F_xi x F_eta) for field_T(points) -> (..., 3)."""

        def g(xi, eta):
            xi, eta = _grid(xi, eta)
            b = np.asarray(field_T(self.chart(xi, eta)), dtype=float)
            return np.einsum("...k,...k->...", b, self.area_vector(xi, eta))

        return g


@dataclass(frozen=True)
class PlanarPatch(_Chart):
    """Regular planar parallelogram F(xi, eta) = origin + xi e_xi + eta e_eta."""

    origin_m: NDArray[np.float64]
    e_xi_m: NDArray[np.float64]
    e_eta_m: NDArray[np.float64]

    @property
    def area_vector_m2(self) -> NDArray[np.float64]:
        return np.cross(self.e_xi_m, self.e_eta_m)

    def chart(self, xi, eta):
        xi, eta = _grid(xi, eta)
        return (self.origin_m + xi[..., None] * self.e_xi_m
                + eta[..., None] * self.e_eta_m)

    def area_vector(self, xi, eta):
        xi, _ = _grid(xi, eta)
        return np.broadcast_to(self.area_vector_m2, xi.shape + (3,))

    def metric(self, xi, eta):
        xi, _ = _grid(xi, eta)
        return (np.full(xi.shape, np.linalg.norm(self.e_xi_m)),
                np.full(xi.shape, np.linalg.norm(self.e_eta_m)))

    def locate(self, points_m):
        """Signed offset along the area normal and chart coordinates of points."""
        d = np.asarray(points_m, dtype=float) - self.origin_m
        normal = self.area_vector_m2 / np.linalg.norm(self.area_vector_m2)
        s = d @ normal
        a, b = self.e_xi_m, self.e_eta_m
        gram = np.array([[a @ a, a @ b], [a @ b, b @ b]])
        coords = np.linalg.solve(gram, np.stack((d @ a, d @ b), axis=-1)[..., None])[..., 0]
        return s, coords[..., 0], coords[..., 1]


def _frame(normal, reference):
    """Unit normal n and in-plane axes e1, e2 = n x e1 (e1 from reference)."""
    n = np.asarray(normal, dtype=float)
    if n.shape != (3,) or not np.all(np.isfinite(n)) or np.linalg.norm(n) == 0:
        raise ValueError("the axis or normal must be a finite nonzero 3-vector")
    n = n / np.linalg.norm(n)
    if reference is None:
        reference = np.eye(3)[int(np.argmin(np.abs(n)))]
    r = np.asarray(reference, dtype=float)
    if r.shape != (3,) or not np.all(np.isfinite(r)):
        raise ValueError("reference must be a finite 3-vector")
    e1 = r - (r @ n) * n
    if np.linalg.norm(e1) <= 1e-12 * max(np.linalg.norm(r), 1e-300):
        raise ValueError("reference must not be parallel to the axis or normal")
    e1 = e1 / np.linalg.norm(e1)
    return n, e1, np.cross(n, e1)


def _theta_range(theta_range):
    t = np.asarray(theta_range, dtype=float)
    if (t.shape != (2,) or not np.all(np.isfinite(t))
            or not 0 < t[1] - t[0] <= 2 * np.pi + 1e-12):
        raise ValueError("theta_range must be [theta0, theta1] radians with "
                         "0 < theta1 - theta0 <= 2 pi")
    return float(t[0]), float(min(t[1] - t[0], 2 * np.pi))


@dataclass(frozen=True)
class AnnularPatch(_Chart):
    """Planar annulus sector: r = r_inner + xi (r_outer - r_inner), theta = theta0 + eta dtheta.

    The area normal is the given normal; theta is measured from the
    reference direction projected into the plane, toward normal x reference.
    """

    center_m: NDArray[np.float64]
    normal: NDArray[np.float64]
    e1: NDArray[np.float64]
    e2: NDArray[np.float64]
    r_inner_m: float
    r_outer_m: float
    theta0: float
    dtheta: float

    def _polar(self, xi, eta):
        xi, eta = _grid(xi, eta)
        return (self.r_inner_m + xi * (self.r_outer_m - self.r_inner_m),
                self.theta0 + eta * self.dtheta)

    def chart(self, xi, eta):
        r, th = self._polar(xi, eta)
        return (self.center_m + (r * np.cos(th))[..., None] * self.e1
                + (r * np.sin(th))[..., None] * self.e2)

    def area_vector(self, xi, eta):
        r, _ = self._polar(xi, eta)
        return ((self.r_outer_m - self.r_inner_m) * self.dtheta * r)[..., None] * self.normal

    def metric(self, xi, eta):
        r, _ = self._polar(xi, eta)
        return np.full(r.shape, self.r_outer_m - self.r_inner_m), r * self.dtheta

    def locate(self, points_m):
        d = np.asarray(points_m, dtype=float) - self.center_m
        x, y = d @ self.e1, d @ self.e2
        xi = (np.hypot(x, y) - self.r_inner_m) / (self.r_outer_m - self.r_inner_m)
        eta = np.mod(np.arctan2(y, x) - self.theta0, 2 * np.pi) / self.dtheta
        return d @ self.normal, xi, eta


@dataclass(frozen=True)
class CylindricalPatch(_Chart):
    """Circular cylinder sector: theta = theta0 + xi dtheta, height eta * |axis|.

    The area normal points away from the axis.
    """

    axis_origin_m: NDArray[np.float64]
    axis_unit: NDArray[np.float64]
    e1: NDArray[np.float64]
    e2: NDArray[np.float64]
    height_m: float
    radius_m: float
    theta0: float
    dtheta: float

    def chart(self, xi, eta):
        xi, eta = _grid(xi, eta)
        th = self.theta0 + xi * self.dtheta
        return (self.axis_origin_m + (self.radius_m * np.cos(th))[..., None] * self.e1
                + (self.radius_m * np.sin(th))[..., None] * self.e2
                + (eta * self.height_m)[..., None] * self.axis_unit)

    def area_vector(self, xi, eta):
        xi, _ = _grid(xi, eta)
        th = self.theta0 + xi * self.dtheta
        radial = np.cos(th)[..., None] * self.e1 + np.sin(th)[..., None] * self.e2
        return self.radius_m * self.dtheta * self.height_m * radial

    def metric(self, xi, eta):
        xi, _ = _grid(xi, eta)
        return (np.full(xi.shape, self.radius_m * self.dtheta),
                np.full(xi.shape, self.height_m))

    def locate(self, points_m):
        d = np.asarray(points_m, dtype=float) - self.axis_origin_m
        x, y = d @ self.e1, d @ self.e2
        xi = np.mod(np.arctan2(y, x) - self.theta0, 2 * np.pi) / self.dtheta
        return np.hypot(x, y) - self.radius_m, xi, (d @ self.axis_unit) / self.height_m

    def monotone_points(self, points_m):
        """Add each segment's point nearest the axis: the radius is convex along
        a segment, so a chord entering and leaving the cylinder gets both roots."""
        p = np.asarray(points_m, dtype=float)
        d = p - self.axis_origin_m
        perp = d - np.outer(d @ self.axis_unit, self.axis_unit)
        step = np.diff(perp, axis=0)
        length2 = np.einsum("ij,ij->i", step, step)
        t = np.where(length2 > 0, -np.einsum("ij,ij->i", perp[:-1], step)
                     / np.where(length2 > 0, length2, 1.0), -1.0)
        inner = np.nonzero((t > 1e-12) & (t < 1 - 1e-12))[0]
        if not len(inner):
            return p, np.arange(len(p))
        feet = p[inner] + t[inner, None] * (p[inner + 1] - p[inner])
        index = np.arange(len(p)) + np.searchsorted(inner, np.arange(len(p)), side="left")
        out = np.empty((len(p) + len(inner), 3))
        out[index] = p
        out[inner + np.arange(1, len(inner) + 1)] = feet
        return out, index


def surface_crossings(points_m: ArrayLike, chart, *, seed_vertex: int,
                      closed: bool) -> list[int]:
    """Signs (+1 along the area normal) of a polyline's crossings of a seed surface.

    The polyline is first made monotone in the signed offset (see
    monotone_points). A crossing is a change of sign between consecutive
    vertices off the surface; vertices on it (|offset| <= 1e-12 of the
    polyline's coordinate scale) do not count, so a line that touches the
    surface and returns is no crossing and one that passes through a vertex
    on it is one. The crossing must lie inside the chart's [0, 1]^2. The
    crossing at the seed vertex, and for a closed line at its first and last
    vertex (the seed), is not counted.
    """
    original = np.asarray(points_m, dtype=float)
    if len(original) < 2:
        return []
    pts, index = chart.monotone_points(original)
    s, _, _ = chart.locate(pts)
    tol = 1e-12 * max(float(np.max(np.abs(pts))), 1e-300)
    side = np.where(np.abs(s) <= tol, 0, np.sign(s)).astype(int)
    excluded = {int(index[seed_vertex])}
    if closed:
        excluded.update({0, len(pts) - 1})
    off = np.nonzero(side)[0]
    signs = []
    for i, j in zip(off[:-1], off[1:]):
        if side[i] == side[j] or any(i <= e <= j for e in excluded):
            continue
        if j == i + 1:
            at = pts[i] + s[i] / (s[i] - s[j]) * (pts[j] - pts[i])
        else:
            at = pts[i + 1]  # the first vertex on the surface
        _, xi, eta = chart.locate(at)
        if -1e-9 <= xi <= 1 + 1e-9 and -1e-9 <= eta <= 1 + 1e-9:
            signs.append(int(side[j]))
    return signs


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


def annular_patch(center_m: ArrayLike, normal: ArrayLike, r_outer_m: float, *,
                  r_inner_m: float = 0.0, theta_range=(0.0, 2 * np.pi),
                  reference: ArrayLike | None = None) -> AnnularPatch:
    """Build a planar annulus (or disk, r_inner_m = 0) sector chart."""
    c = np.asarray(center_m, dtype=float)
    if c.shape != (3,) or not np.all(np.isfinite(c)):
        raise ValueError("center_m must be a finite 3-vector")
    n, e1, e2 = _frame(normal, reference)
    r0, r1 = float(r_inner_m), float(r_outer_m)
    if not (np.isfinite(r0) and np.isfinite(r1) and 0 <= r0 < r1):
        raise ValueError("radii must satisfy 0 <= r_inner_m < r_outer_m")
    return AnnularPatch(c, n, e1, e2, r0, r1, *_theta_range(theta_range))


def cylindrical_patch(axis_origin_m: ArrayLike, axis_m: ArrayLike, radius_m: float, *,
                      theta_range=(0.0, 2 * np.pi),
                      reference: ArrayLike | None = None) -> CylindricalPatch:
    """Build a cylinder sector chart; axis_m spans the height from axis_origin_m."""
    o = np.asarray(axis_origin_m, dtype=float)
    if o.shape != (3,) or not np.all(np.isfinite(o)):
        raise ValueError("axis_origin_m must be a finite 3-vector")
    a, e1, e2 = _frame(axis_m, reference)
    radius = float(radius_m)
    if not (np.isfinite(radius) and radius > 0):
        raise ValueError("radius_m must be positive")
    return CylindricalPatch(o, a, e1, e2, float(np.linalg.norm(np.asarray(axis_m, float))),
                            radius, *_theta_range(theta_range))


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
    denominator = d0 + disc
    # A zero-density segment start with nothing left to place stays at its start.
    t = np.where(denominator > 0, 2 * rest / np.where(denominator > 0, denominator, 1.0), 0.0)
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
    metric: Callable | None = None,
) -> EqualFluxSeeds:
    """Place n_u * n_v equal-flux seeds by the two-stage inverse.

    flux_density(xi, eta) broadcasts and returns g in Wb per unit chart area
    (T m^2 when xi, eta are dimensionless); chart(xi, eta) returns physical
    points in metres with shape (..., 3). With metric(xi, eta) -> (|F_xi|,
    |F_eta|) the diagnostics also report the physical cell aspect ratios. The table rows xi_nodes x eta_nodes
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
    diagnostics = {
        "forward_residual_v_max": float(np.max(np.abs(v_check - v_c))),
        "forward_residual_u_max": float(np.max(np.abs(u_check - u_c[None, :]))),
        "table_shape": [len(eta), len(xi)],
        "seed_rows_resampled": int(n_v),
    }
    if metric is not None:
        # Cell edges of the same model: eta edges from the marginal, xi edges
        # from the conditional on each seed row.
        eta_edges = _invert(eta, rho[None, :], marginal,
                            total * np.linspace(0, 1, n_v + 1)[None, :])[0]
        eta_edges[0], eta_edges[-1] = eta[0], eta[-1]
        xi_edges = _invert(xi, rows, row_cum,
                           row_rho[:, None] * np.linspace(0, 1, n_u + 1)[None, :])
        xi_edges[:, 0], xi_edges[:, -1] = xi[0], xi[-1]
        lx, ly = metric(xi_c, np.broadcast_to(eta_c[:, None], xi_c.shape))
        width = lx * np.diff(xi_edges, axis=1)
        height = ly * np.diff(eta_edges)[:, None]
        diagnostics.update(_aspect_summary(width, height))
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
        diagnostics=diagnostics,
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


def _aspect_summary(width: NDArray, height: NDArray) -> dict:
    """Physical cell aspect ratios (>= 1) of cells with the given side lengths."""
    width, height = np.abs(np.ravel(width)), np.abs(np.ravel(height))
    aspect = np.maximum(width, height) / np.maximum(np.minimum(width, height), 1e-300)
    return {"max_cell_aspect": float(np.max(aspect)),
            "median_cell_aspect": float(np.median(aspect))}


def _child_aspect(a: float, b: float, k: int) -> float:
    """Best attainable cell aspect of an a x b rectangle cut into k cells.

    Elongation up to k can be absorbed by laying the k cells along the long
    side; the rest remains in every cell.
    """
    long_side, short_side = max(a, b), max(min(a, b), 1e-300)
    return max(1.0, long_side / short_side / k)


def _pl_eval(nodes: NDArray, values: NDArray, at: NDArray) -> NDArray:
    """Linear interpolation of values (..., len(nodes)) at points `at`."""
    at = np.asarray(at, dtype=float)
    k = np.clip(np.searchsorted(nodes, at, side="right") - 1, 0, len(nodes) - 2)
    w = (at - nodes[k]) / (nodes[k + 1] - nodes[k])
    return values[..., k] * (1 - w) + values[..., k + 1] * w


def _breaks(nodes: NDArray, a: float, b: float) -> NDArray:
    return np.concatenate(([a], nodes[(nodes > a) & (nodes < b)], [b]))


def _pl_integral(nodes: NDArray, values: NDArray, a: float, b: float) -> NDArray:
    """Exact integral over [a, b] of the piecewise-linear interpolant."""
    pts = _breaks(nodes, a, b)
    vals = _pl_eval(nodes, values, pts)
    return np.sum(0.5 * (vals[..., 1:] + vals[..., :-1]) * np.diff(pts), axis=-1)


def _marginal(table: NDArray, a_nodes: NDArray, b_nodes: NDArray, a0, a1, b0, b1):
    """Marginal flux density along a of the bilinear table over a rectangle.

    table has shape (len(b_nodes), len(a_nodes)). The marginal is exactly
    piecewise linear in a with breaks at the returned points.
    """
    pts = _breaks(a_nodes, a0, a1)
    columns = _pl_eval(a_nodes, table, pts)          # (len(b_nodes), len(pts))
    return pts, _pl_integral(b_nodes, columns.T, b0, b1)


def _split(pts: NDArray, density: NDArray, fraction: float) -> tuple[float, float]:
    """Point where the cumulative of a piecewise-linear density reaches fraction."""
    cumulative = _cumulative(pts, density[None, :])
    total = float(cumulative[0, -1])
    at = _invert(pts, density[None, :], cumulative, np.array([[total * fraction]]))[0, 0]
    return float(at), total


def _mixed_cells(table: NDArray) -> NDArray:
    """Table cells whose corner values take both strict signs."""
    corners = np.stack((table[:-1, :-1], table[:-1, 1:], table[1:, :-1], table[1:, 1:]))
    return np.any(corners > 0, axis=0) & np.any(corners < 0, axis=0)


def _edge_zeros(nodes: NDArray, values: NDArray) -> NDArray:
    """Linear zero crossings along the last axis between strictly opposite signs."""
    a, b = values[..., :-1], values[..., 1:]
    hit = a * b < 0
    rows, k = np.nonzero(hit.reshape(-1, hit.shape[-1]))
    a, b = a.reshape(-1, a.shape[-1])[rows, k], b.reshape(-1, b.shape[-1])[rows, k]
    return nodes[k] + a / (a - b) * (nodes[k + 1] - nodes[k])


def _resolve_zero_line(flux_density: Callable, xi: NDArray, eta: NDArray, table: NDArray):
    """Insert table nodes where cell edges cross zero, and resample the density.

    The clipped bilinear interpolant puts support on both sides of a zero line
    that runs through a cell. Adding the edge crossings of the mixed cells as
    nodes (one pass) makes a zero line that is straight along either chart
    axis lie on table lines; a curved one leaves smaller mixed cells, whose
    share of the flux is reported.
    """
    mixed = _mixed_cells(table)
    if not np.any(mixed):
        return xi, eta, table
    rows, cols = np.nonzero(mixed)
    # Edges of mixed cells only: horizontal edges (constant eta) give xi zeros,
    # vertical edges (constant xi) give eta zeros.
    new_xi, new_eta = [], []
    for r, c in zip(rows, cols):
        for rr in (r, r + 1):
            new_xi.append(_edge_zeros(xi[c:c + 2], table[rr, c:c + 2]))
        for cc in (c, c + 1):
            new_eta.append(_edge_zeros(eta[r:r + 2], table[r:r + 2, cc]))

    def merged(nodes, extra):
        extra = np.concatenate(extra) if extra else np.empty(0)
        out = np.unique(np.concatenate((nodes, extra)))
        keep = np.concatenate(([True], np.diff(out) > 1e-12 * (nodes[-1] - nodes[0])))
        return out[keep]

    xi, eta = merged(xi, new_xi), merged(eta, new_eta)
    table = _row_density(flux_density, xi, eta)
    if not np.all(np.isfinite(table)):
        raise ValueError("flux density must be finite (refined table)")
    return xi, eta, table


def balanced_equal_flux_seeds(
    flux_density: Callable,
    *,
    xi_nodes: ArrayLike,
    eta_nodes: ArrayLike,
    n_lines: int,
    chart: Callable,
    metric: Callable,
    sign: str = "auto",
) -> BalancedFluxSeeds:
    """Place n_lines equal-flux seeds by recursive bisection of one flux sign.

    flux_density, chart and the table nodes are as for equal_flux_seeds;
    metric(xi, eta) -> (|F_xi|, |F_eta|) gives physical lengths per unit
    chart coordinate. The density is the bilinear interpolant of the table,
    clipped to the selected sign at the nodes. Where the normal flux changes
    sign inside a table cell, the edge zero crossings of that cell are first
    added as table nodes and the density resampled there, so clipping does
    not spread the support across the zero line; diagnostics report the
    added nodes and the share of the seeded flux still in cells whose
    corners take both signs. A rectangle holding n seeds
    is cut across its longer physical side where the marginal flux reaches
    floor(n / 2) / n of its own, until each cell holds one seed; the seed is
    the flux median of its cell (marginal along the longer side, then the
    conditional at that point). Every cell then carries total / n_lines of
    the model flux exactly.

    sign "auto" seeds the only sign present and raises when the normal flux
    changes sign; "positive" or "negative" seeds that part. Field lines that
    cross the surface only with the other sign are not seeded.

    Raises ValueError for nonfinite samples, a zero or absent selected flux,
    a bad sign or n_lines, nonfinite chart points, or a seed where the
    caller's flux_density has the other sign (refine the table there).
    """
    xi, eta = _axis(xi_nodes, "xi_nodes"), _axis(eta_nodes, "eta_nodes")
    n_lines = int(n_lines)
    if n_lines < 1:
        raise ValueError("n_lines must be positive")
    if sign not in SIGNS:
        raise ValueError(f"sign must be one of {list(SIGNS)}")
    table = _row_density(flux_density, xi, eta)
    if not np.all(np.isfinite(table)):
        raise ValueError("flux density must be finite (table)")
    given_shape = [len(eta), len(xi)]
    xi, eta, table = _resolve_zero_line(flux_density, xi, eta, table)

    def flux_of(part):
        return float(_pl_integral(eta, _pl_integral(xi, part, xi[0], xi[-1]), eta[0], eta[-1]))

    positive, negative = flux_of(np.maximum(table, 0)), flux_of(np.maximum(-table, 0))
    if positive == 0 and negative == 0:
        raise ValueError("flux density is identically zero (table)")
    if sign == "auto":
        if positive > 0 and negative > 0:
            raise ValueError("normal flux changes sign on the surface; choose sign='positive' "
                             "or 'negative' (lines crossing only with the other sign are "
                             "then not seeded)")
        orientation = 1 if positive > 0 else -1
    else:
        orientation = 1 if sign == "positive" else -1
    density = np.maximum(orientation * table, 0)
    total = positive if orientation > 0 else negative
    if not total > 0:
        raise ValueError(f"no {'positive' if orientation > 0 else 'negative'} flux "
                         "through the surface")

    cells, seeds_chart = [], []
    stack = [((xi[0], xi[-1], eta[0], eta[-1]), n_lines)]
    while stack:
        (x0, x1, y0, y1), n = stack.pop()
        lx, ly = metric(np.array(0.5 * (x0 + x1)), np.array(0.5 * (y0 + y1)))
        along_xi = float(lx) * (x1 - x0) >= float(ly) * (y1 - y0)
        if along_xi:
            pts, m = _marginal(density, xi, eta, x0, x1, y0, y1)
        else:
            pts, m = _marginal(density.T, eta, xi, y0, y1, x0, x1)
        if n == 1:
            at, _ = _split(pts, m, 0.5)
            if along_xi:
                ys = _breaks(eta, y0, y1)
                column = _pl_eval(eta, _pl_eval(xi, density, np.array([at]))[:, 0], ys)
                other, row_total = _split(ys, column, 0.5)
                point = (at, other)
            else:
                xs = _breaks(xi, x0, x1)
                row = _pl_eval(xi, _pl_eval(eta, density.T, np.array([at]))[:, 0], xs)
                other, row_total = _split(xs, row, 0.5)
                point = (other, at)
            if not row_total > 0:
                raise ValueError("a seed cell has no flux on its median line; refine the table")
            cells.append((x0, x1, y0, y1))
            seeds_chart.append(point)
            continue
        left = n // 2
        width, height = float(lx) * (x1 - x0), float(ly) * (y1 - y0)
        options = []
        for cut_xi in (along_xi, not along_xi):  # longer side first wins ties
            if cut_xi == along_xi:
                cut_pts, cut_m = pts, m
            elif cut_xi:
                cut_pts, cut_m = _marginal(density, xi, eta, x0, x1, y0, y1)
            else:
                cut_pts, cut_m = _marginal(density.T, eta, xi, y0, y1, x0, x1)
            cut, _ = _split(cut_pts, cut_m, left / n)
            if cut_xi:
                sides = [(float(lx) * (cut - x0), height), (float(lx) * (x1 - cut), height)]
                parts = ((x0, cut, y0, y1), (cut, x1, y0, y1))
            else:
                sides = [(width, float(ly) * (cut - y0)), (width, float(ly) * (y1 - cut))]
                parts = ((x0, x1, y0, cut), (x0, x1, cut, y1))
            # Judge a cut by its children's worst aspect ratio, each child seen
            # as the square-ish cells it will be divided into.
            worst = max(_child_aspect(a, b, k) for (a, b), k in zip(sides, (left, n - left)))
            options.append((worst, parts))
        parts = min(options, key=lambda option: option[0])[1]
        stack.append((parts[1], n - left))
        stack.append((parts[0], left))

    cells = np.asarray(cells, dtype=float)
    chart_coords = np.asarray(seeds_chart, dtype=float)
    at_seeds = orientation * np.asarray(flux_density(chart_coords[:, 0], chart_coords[:, 1]),
                                        dtype=float)
    wrong = np.nonzero(~(at_seeds >= 0))[0]
    if len(wrong):
        raise ValueError(f"{len(wrong)} of {n_lines} seeds lie where the normal flux has the "
                         f"other sign (first at chart {chart_coords[wrong[0]].tolist()}); "
                         "the table does not resolve the zero line, refine it")
    corner_mean = 0.25 * (density[:-1, :-1] + density[:-1, 1:] + density[1:, :-1]
                          + density[1:, 1:]) * np.outer(np.diff(eta), np.diff(xi))
    mixed_share = float(corner_mean[_mixed_cells(table)].sum() / total)
    points = np.asarray(chart(chart_coords[:, 0], chart_coords[:, 1]), dtype=float)
    if points.shape != (n_lines, 3) or not np.all(np.isfinite(points)):
        raise ValueError("chart must return finite 3-D points for every seed")
    centre_x, centre_y = 0.5 * (cells[:, 0] + cells[:, 1]), 0.5 * (cells[:, 2] + cells[:, 3])
    lx, ly = metric(centre_x, centre_y)
    signed_total = orientation * total
    return BalancedFluxSeeds(
        ids=np.arange(n_lines, dtype=np.int64),
        points_m=points,
        chart_coords=chart_coords,
        cells=cells,
        n_lines=n_lines,
        orientation=orientation,
        total_flux_Wb=signed_total,
        target_cell_flux_Wb=signed_total / n_lines,
        positive_flux_Wb=positive,
        negative_flux_Wb=-negative,
        diagnostics={"table_shape": [len(eta), len(xi)],
                     "given_table_shape": given_shape,
                     "mixed_cell_flux_share": mixed_share,
                     **_aspect_summary(lx * (cells[:, 1] - cells[:, 0]),
                                       ly * (cells[:, 3] - cells[:, 2]))},
        provenance={"method": BALANCED_METHOD, "module": __name__},
    )


def balanced_cell_flux_quadrature(
    flux_density: Callable,
    seeds: BalancedFluxSeeds,
    *,
    xi_nodes: ArrayLike,
    eta_nodes: ArrayLike,
    order: int = 6,
) -> dict:
    """Integrate the original one-signed density over each bisection cell.

    Each cell is split at the table nodes and integrated by tensor Gauss
    rules of the caller's flux_density clipped to the seeded sign, never the
    table; the total is an independent composite Gauss integral over the
    table cells. Returns fluxes in Wb with the seeded sign.
    """
    xi, eta = _axis(xi_nodes, "xi_nodes"), _axis(eta_nodes, "eta_nodes")
    t, w = _gauss(order)
    xs_q, ys_q, ws_q, owner = [], [], [], []
    for index, (x0, x1, y0, y1) in enumerate(seeds.cells):
        xb, yb = _breaks(xi, x0, x1), _breaks(eta, y0, y1)
        qx = (xb[:-1, None] + np.diff(xb)[:, None] * t).ravel()
        wx = (np.diff(xb)[:, None] * w).ravel()
        qy = (yb[:-1, None] + np.diff(yb)[:, None] * t).ravel()
        wy = (np.diff(yb)[:, None] * w).ravel()
        gx, gy = np.meshgrid(qx, qy, indexing="xy")
        xs_q.append(gx.ravel())
        ys_q.append(gy.ravel())
        ws_q.append(np.outer(wy, wx).ravel())
        owner.append(np.full(gx.size, index))
    xs_q, ys_q = np.concatenate(xs_q), np.concatenate(ys_q)
    ws_q, owner = np.concatenate(ws_q), np.concatenate(owner)
    sign = seeds.orientation
    g = np.maximum(sign * np.asarray(flux_density(xs_q, ys_q), dtype=float), 0)
    cell = np.bincount(owner, weights=g * ws_q, minlength=len(seeds.cells))

    xi_t = (xi[:-1, None] + np.diff(xi)[:, None] * t).ravel()
    xi_w = (np.diff(xi)[:, None] * w).ravel()
    eta_t = (eta[:-1, None] + np.diff(eta)[:, None] * t).ravel()
    eta_w = (np.diff(eta)[:, None] * w).ravel()
    whole = np.maximum(sign * np.asarray(flux_density(xi_t[None, :], eta_t[:, None]),
                                         dtype=float), 0)
    total = float(eta_w @ whole @ xi_w)
    target = total / seeds.n_lines
    return {
        "cell_flux_Wb": sign * cell,
        "independent_total_flux_Wb": sign * total,
        "max_relative_cell_deviation": float(np.max(np.abs(cell / target - 1))),
        "coverage": float(cell.sum() / total),
        "order": int(order),
    }
