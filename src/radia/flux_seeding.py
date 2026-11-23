"""Equal-flux field-line seeds on a prescribed oriented surface chart."""

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import cumulative_trapezoid, trapezoid


@dataclass(frozen=True)
class EqualFluxSeeds:
    """Physical centres, matching flux-coordinate cells and modelled fluxes.

    The cells_uv columns are u_lo, u_hi, v_lo, v_hi. Fluxes in Wb are
    signed with respect to the chart orientation. Coverage is the
    complete-cell flux divided by the total sampled patch flux.
    """

    points_m: NDArray[np.float64]
    cells_uv: NDArray[np.float64]
    signed_cell_flux_Wb: float
    signed_patch_flux_Wb: float

    @property
    def coverage(self) -> float:
        return len(self.points_m) * (self.signed_cell_flux_Wb / self.signed_patch_flux_Wb)


def _axis(values: ArrayLike, name: str) -> NDArray[np.float64]:
    axis = np.asarray(values, dtype=float)
    if (axis.ndim != 1 or len(axis) < 2 or not np.all(np.isfinite(axis))
            or not np.all(np.diff(axis) > 0)):
        raise ValueError(f"{name} must be a finite, strictly increasing 1-D grid")
    return axis


def equal_flux_seeds(
    chart: Callable,
    flux_density: Callable,
    *,
    s_nodes: ArrayLike,
    v_nodes: ArrayLike,
    u_edges: ArrayLike,
    v_edges: ArrayLike,
    s_reference: float,
) -> EqualFluxSeeds:
    """Invert the cumulative flux at fixed v and seed each complete cell.

    chart(s, v) maps broadcast arrays to physical points in metres with
    shape (..., 3). flux_density(s, v) returns the signed coefficient of
    the pulled-back magnetic flux form, B dot (X_s cross X_v). Supply a
    one-to-one oriented chart and a field of one sign on the patch.
    Sign and non-vanishing are checked at the sampled interior points.
    The coefficient may vanish on the s boundaries and at isolated
    boundary points. A v slice on which it vanishes identically, such as
    a v boundary with no normal flux, has no strictly monotone cumulative
    flux and is rejected.

    Integrate in s with the trapezoidal rule and set u=0 at s_reference.
    The sampled cumulative integral is piecewise linear in s and v. A
    cell is retained if it fits the interpolated patch at both v edges
    and every intervening table knot. Equal spacing in u and v gives
    equal flux weights. An empty complete-cell set is valid.

    The reported fluxes are those of the sampled map. To verify against
    the original field, integrate it over the mapped cells.
    """
    s = _axis(s_nodes, "s_nodes")
    v = _axis(v_nodes, "v_nodes")
    u_grid = _axis(u_edges, "u_edges")
    v_grid = _axis(v_edges, "v_edges")
    if not s[0] <= s_reference <= s[-1]:
        raise ValueError("s_reference must lie inside the sampled s interval")
    if v_grid[0] < v[0] or v_grid[-1] > v[-1]:
        raise ValueError("v_edges must lie inside the sampled v interval")

    du, dv = np.diff(u_grid), np.diff(v_grid)
    for steps, grid, name in ((du, u_grid, "u_edges"), (dv, v_grid, "v_edges")):
        roundoff = 8 * np.finfo(float).eps * np.max(np.abs(grid))
        if not np.all(np.abs(steps - steps[0]) <= roundoff):
            raise ValueError(f"{name} must be uniformly spaced")

    g = np.broadcast_to(np.asarray(flux_density(s[:, None], v[None, :]), dtype=float),
                        (len(s), len(v)))
    if not np.all(np.isfinite(g)):
        raise ValueError("flux density samples must be finite")
    positive = np.any(g > 0)
    if positive == np.any(g < 0) or np.any(g[1:-1, 1:-1] == 0):
        raise ValueError("flux density must have one nonzero sign on the patch interior")
    sign = 1 if positive else -1
    u = cumulative_trapezoid(g, s, axis=0, initial=0)
    ref_row = min(np.searchsorted(s, s_reference, side="right") - 1, len(s) - 2)
    ref_fraction = (s_reference - s[ref_row]) / (s[ref_row + 1] - s[ref_row])
    u -= (1 - ref_fraction) * u[ref_row] + ref_fraction * u[ref_row + 1]
    if not np.all(np.isfinite(u)) or np.any(sign * np.diff(u, axis=0) <= 0):
        raise ValueError("cumulative flux must be finite and strictly monotone on each v slice")

    if sign < 0:
        u = u[::-1]
        s = s[::-1]
    total = float(trapezoid(u[-1] - u[0], v))
    if not np.isfinite(total) or total <= 0:
        raise ValueError("patch flux magnitude must be finite and positive")

    cells, seed_s, seed_v = [], [], []
    for v_lo, v_hi in zip(v_grid[:-1], v_grid[1:]):
        knots = v[(v > v_lo) & (v < v_hi)]
        samples = np.r_[v_lo, knots, v_hi]
        left = np.interp(samples, v, u[0]).max()
        right = np.interp(samples, v, u[-1]).min()
        indices = np.flatnonzero((u_grid[:-1] >= left) & (u_grid[1:] <= right))
        if not len(indices):
            continue
        lo, hi = u_grid[indices], u_grid[indices + 1]
        v_mid = v_lo + (v_hi - v_lo) / 2
        j = max(np.searchsorted(v, v_mid) - 1, 0)
        frac = (v_mid - v[j]) / (v[j + 1] - v[j])
        u_slice = (1 - frac) * u[:, j] + frac * u[:, j + 1]
        seed_s.extend(np.interp(lo + (hi - lo) / 2, u_slice, s))
        seed_v.extend([v_mid] * len(lo))
        cells.extend(zip(lo, hi, [v_lo] * len(lo), [v_hi] * len(lo)))

    points = (np.asarray(chart(np.asarray(seed_s), np.asarray(seed_v)), dtype=float)
              if cells else np.empty((0, 3)))
    if points.shape != (len(cells), 3) or not np.all(np.isfinite(points)):
        raise ValueError("chart must return finite 3-D points for every seed")
    return EqualFluxSeeds(points, np.asarray(cells, dtype=float).reshape(-1, 4),
                         float(sign * du[0] * dv[0]), sign * total)
