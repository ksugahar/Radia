"""coil_profile.py

Cross-section profile abstractions for filament-approximation PEEC.

A Profile describes a 2D cross-section in the conductor's local frame
(u = x_hat axis = width direction, v = z_hat axis = height direction).
Each profile exposes a canonical normalized parametrization (alpha, beta)
in [0, 1]^2 that is AREA-PRESERVING: the Jacobian of (alpha, beta) ->
(u, v) is constant and equals total_area(). This means:

    - filament k at (alpha_k, beta_k) has a consistent position across
      profiles of different types when evaluated via sample_at();
    - cell areas A(alpha, beta) are uniform, so a uniform (n_a x n_b)
      grid on [0, 1]^2 gives cells of equal area regardless of profile
      shape;
    - cross-type loft (e.g. rect -> circle) between two area-preserving
      profiles preserves DC-current uniformity: each filament carries
      the same integral int(1 / (sigma * A(s))) |dpath/ds| ds because
      A(s) is the same for all filaments.

Stage 3a: RectProfile.
Stage 3b (this file): CircleProfile, AnnularProfile, InterpolatedProfile.
Stage 3c: OCCProfile (from build123d Face via UV parametrization).
"""

from abc import ABC, abstractmethod
import numpy as np


class Profile(ABC):
    """Abstract 2D cross-section profile with area-preserving (alpha, beta)
    parametrization.

    Subclasses MUST override:
        sample_at(alpha, beta) -- local-frame (u, v)
        total_area()           -- area of the profile
        bounding_wh()          -- axis-aligned bounding (w, h)
        _interpolate_same(other, t) -- same-type loft (may return NotImpl)

    Defaults provided for:
        sample(n_a, n_b)       -- uniform (n_a x n_b) grid on [0,1]^2
                                  pushed through sample_at
        cell_wh(n_a, n_b)      -- isotropic (sqrt(A), sqrt(A)); subclasses
                                  can override for anisotropic cells
        interpolate(other, t)  -- same-type via _interpolate_same,
                                  cross-type via InterpolatedProfile
    """

    @abstractmethod
    def sample_at(self, alpha, beta):
        """(u, v) at normalized (alpha, beta) in [0,1]^2. Supports arrays."""

    @abstractmethod
    def sample_perimeter(self, n):
        """Return (u, v) arrays of length n, equally spaced by arc length
        along the outer boundary of the cross-section.

        Used for perimeter-only filament placement (thin-skin limit, IH
        typical use).  Each sample is cell-centered on an arc-length bin:
        s_k = (k + 0.5) * L_perimeter / n for k in [0, n).

        The effective cross-section area per filament is taken as
        ``total_area() / n`` (all perimeter filaments share the solid
        cross-section equally, matching DC resistance of the bundle to
        the solid conductor).
        """

    @abstractmethod
    def total_area(self):
        """Total cross-section area."""

    @abstractmethod
    def bounding_wh(self):
        """(w, h) bounding box."""

    @abstractmethod
    def _interpolate_same(self, other, t):
        """Same-type interpolation. Return an instance of the same Profile
        subclass at path parameter t in [0, 1]. Raise NotImplementedError
        if not supported."""

    def sample(self, n_a, n_b):
        """Cell-center (u, v, area) on uniform (n_a, n_b) grid in [0,1]^2.

        alpha_k = (i + 0.5) / n_a, beta_k = (j + 0.5) / n_b for indices
        i in [0, n_a), j in [0, n_b); flatten index k = i * n_b + j.
        """
        n_a = int(n_a)
        n_b = int(n_b)
        i_idx, j_idx = np.meshgrid(np.arange(n_a), np.arange(n_b),
                                    indexing='ij')
        alphas = ((i_idx + 0.5) / n_a).ravel()
        betas = ((j_idx + 0.5) / n_b).ravel()
        u, v = self.sample_at(alphas, betas)
        A = self.total_area() / (n_a * n_b)
        area = np.full(u.size, A)
        return u, v, area

    def cell_wh(self, n_a, n_b):
        """Equivalent-square sub-filament dimensions (dw, dh). Isotropic by
        default: sqrt(A_cell) for both. Anisotropic profiles (e.g. Rect)
        can override."""
        A = self.total_area() / (int(n_a) * int(n_b))
        side = float(np.sqrt(A))
        return side, side

    def interpolate(self, other, t):
        """Profile at path parameter t in [0, 1] between self (t=0) and
        other (t=1). Same-type lofts use the subclass _interpolate_same()
        fast path; cross-type falls back to InterpolatedProfile."""
        if type(self) is type(other):
            try:
                return self._interpolate_same(other, t)
            except NotImplementedError:
                pass
        return InterpolatedProfile(self, other, t)


class RectProfile(Profile):
    """Rectangular cross-section of total size (w, h).

    Canonical parametrization (area-preserving on a rectangle is trivial):
        u = (alpha - 0.5) * w     alpha in [0, 1] -> u in [-w/2, w/2]
        v = (beta - 0.5) * h      beta in [0, 1]  -> v in [-h/2, h/2]
    """

    def __init__(self, w, h):
        self.w = float(w)
        self.h = float(h)

    def sample_at(self, alpha, beta):
        alpha = np.asarray(alpha, dtype=float)
        beta = np.asarray(beta, dtype=float)
        u = (alpha - 0.5) * self.w
        v = (beta - 0.5) * self.h
        return u, v

    def sample_perimeter(self, n):
        n = int(n)
        if n < 1:
            raise ValueError(f"n_peri must be >= 1, got {n}")
        L = 2.0 * (self.w + self.h)
        s = (np.arange(n) + 0.5) * (L / n)
        u = np.empty(n)
        v = np.empty(n)
        # CCW starting from lower-left corner (-w/2, -h/2):
        #   [0, w)        bottom:   u = -w/2 + s,            v = -h/2
        #   [w, w+h)      right:    u =  w/2,                v = -h/2 + (s - w)
        #   [w+h, 2w+h)   top:      u =  w/2 - (s - w - h),  v =  h/2
        #   [2w+h, 2w+2h) left:     u = -w/2,                v =  h/2 - (s - 2w - h)
        w, h = self.w, self.h
        m1 = s < w
        m2 = (s >= w) & (s < w + h)
        m3 = (s >= w + h) & (s < 2 * w + h)
        m4 = s >= 2 * w + h
        u[m1] = -0.5 * w + s[m1];         v[m1] = -0.5 * h
        u[m2] =  0.5 * w;                 v[m2] = -0.5 * h + (s[m2] - w)
        u[m3] =  0.5 * w - (s[m3] - w - h); v[m3] =  0.5 * h
        u[m4] = -0.5 * w;                 v[m4] =  0.5 * h - (s[m4] - 2 * w - h)
        return u, v

    def total_area(self):
        return self.w * self.h

    def bounding_wh(self):
        return self.w, self.h

    def cell_wh(self, n_a, n_b):
        """Anisotropic: sub-filament cell is exactly (w/n_a, h/n_b)."""
        return self.w / int(n_a), self.h / int(n_b)

    def _interpolate_same(self, other, t):
        t = float(t)
        return RectProfile((1.0 - t) * self.w + t * other.w,
                            (1.0 - t) * self.h + t * other.h)

    def __repr__(self):
        return f"RectProfile(w={self.w:g}, h={self.h:g})"


class CircleProfile(Profile):
    """Circular cross-section of radius r.

    Area-preserving polar parametrization:
        r(alpha) = r_outer * sqrt(alpha)   alpha in [0, 1]
        theta(beta) = 2*pi * beta           beta in [0, 1]
        u = r(alpha) * cos(theta(beta))
        v = r(alpha) * sin(theta(beta))

    Jacobian = pi * r_outer^2 = total_area (constant over [0,1]^2), so a
    uniform grid on (alpha, beta) yields cells of equal area.

    alpha is a radial coordinate (0 at the disk center, 1 at the outer
    boundary); beta is angular (0 = +u axis, 0.25 = +v axis, etc).
    """

    def __init__(self, r):
        self.r = float(r)

    def sample_at(self, alpha, beta):
        alpha = np.asarray(alpha, dtype=float)
        beta = np.asarray(beta, dtype=float)
        r = self.r * np.sqrt(np.clip(alpha, 0.0, 1.0))
        theta = 2.0 * np.pi * beta
        return r * np.cos(theta), r * np.sin(theta)

    def sample_perimeter(self, n):
        n = int(n)
        if n < 1:
            raise ValueError(f"n_peri must be >= 1, got {n}")
        theta = 2.0 * np.pi * (np.arange(n) + 0.5) / n
        return self.r * np.cos(theta), self.r * np.sin(theta)

    def total_area(self):
        return np.pi * self.r * self.r

    def bounding_wh(self):
        return 2.0 * self.r, 2.0 * self.r

    def _interpolate_same(self, other, t):
        t = float(t)
        return CircleProfile((1.0 - t) * self.r + t * other.r)

    def __repr__(self):
        return f"CircleProfile(r={self.r:g})"


class AnnularProfile(Profile):
    """Hollow annular cross-section (r_in <= r <= r_out).

    Area-preserving polar parametrization:
        r(alpha) = sqrt(r_in^2 + alpha * (r_out^2 - r_in^2))
        theta(beta) = 2 * pi * beta
        u = r(alpha) * cos(theta)
        v = r(alpha) * sin(theta)

    alpha = 0 maps to the inner boundary, alpha = 1 to the outer boundary.
    Jacobian = pi * (r_out^2 - r_in^2) constant, so cell areas are equal.

    Useful for water-cooled IH copper tubes.
    """

    def __init__(self, r_in, r_out):
        self.r_in = float(r_in)
        self.r_out = float(r_out)
        if not (0.0 <= self.r_in < self.r_out):
            raise ValueError(
                f"AnnularProfile requires 0 <= r_in < r_out; got "
                f"r_in={self.r_in}, r_out={self.r_out}")

    def sample_at(self, alpha, beta):
        alpha = np.asarray(alpha, dtype=float)
        beta = np.asarray(beta, dtype=float)
        rr2 = self.r_in * self.r_in \
              + np.clip(alpha, 0.0, 1.0) * (self.r_out * self.r_out
                                             - self.r_in * self.r_in)
        r = np.sqrt(rr2)
        theta = 2.0 * np.pi * beta
        return r * np.cos(theta), r * np.sin(theta)

    def sample_perimeter(self, n):
        """Outer boundary only.  Inner boundary filaments are not placed;
        for hollow tubes where both surfaces carry skin current, use a
        larger n on the outer or switch to a volume-grid model."""
        n = int(n)
        if n < 1:
            raise ValueError(f"n_peri must be >= 1, got {n}")
        theta = 2.0 * np.pi * (np.arange(n) + 0.5) / n
        return self.r_out * np.cos(theta), self.r_out * np.sin(theta)

    def total_area(self):
        return np.pi * (self.r_out * self.r_out - self.r_in * self.r_in)

    def bounding_wh(self):
        return 2.0 * self.r_out, 2.0 * self.r_out

    def _interpolate_same(self, other, t):
        t = float(t)
        return AnnularProfile(
            (1.0 - t) * self.r_in + t * other.r_in,
            (1.0 - t) * self.r_out + t * other.r_out)

    def __repr__(self):
        return (f"AnnularProfile(r_in={self.r_in:g}, r_out={self.r_out:g})")


class PolygonProfile(Profile):
    """Cross-section defined by an arbitrary (simply-connected, convex or
    star-shaped from centroid) polygon in the local (u, v) plane.

    Vertices are specified counter-clockwise; the polygon is implicitly
    closed (last vertex connects back to the first).

    Parametrization (approximate-area-preserving, exact for centered
    regular polygons):
        beta in [0, 1] -> parametric point B(beta) on the boundary
            (linear interpolation along perimeter parameter);
        alpha in [0, 1] -> sqrt-scaled radial blend between centroid C
            and boundary B:
                (u, v) = C + sqrt(alpha) * (B(beta) - C)

    The sqrt scaling gives Jacobian |J| = pi * d(beta)^2 which is
    exactly area-preserving for rotationally symmetric polygons (e.g.
    regular n-gon from centroid) and a reasonable approximation for
    other simply-connected shapes. For a rectangle, this gives a small
    (~20%) anisotropy in cell area between the short and long axes, so
    DC uniformity is not guaranteed to the same precision as
    RectProfile. The abstraction is meant for CAD-authored shapes that
    don't fit the named-primitive catalogue (RectProfile, CircleProfile,
    AnnularProfile); for named shapes, prefer the dedicated class.

    total_area() is computed exactly from the shoelace formula.
    """

    def __init__(self, vertices, centroid=None):
        """Args:
            vertices: (n, 2) array-like of (u, v) coordinates, CCW ordered.
                Closing edge is implicit.
            centroid: optional explicit (u_c, v_c). Defaults to the
                area-weighted centroid computed from vertices.
        """
        V = np.asarray(vertices, dtype=float)
        if V.ndim != 2 or V.shape[1] != 2 or V.shape[0] < 3:
            raise ValueError(
                f"PolygonProfile vertices must be (n>=3, 2); got {V.shape}")
        self.vertices = V
        self._area = self._signed_area(V)
        if self._area <= 0:
            # Ensure CCW (positive area). Flip if the input was CW.
            self.vertices = V[::-1]
            self._area = -self._area

        if centroid is None:
            centroid = self._area_centroid(self.vertices, self._area)
        self.centroid = np.asarray(centroid, dtype=float)

        # Precompute perimeter and cumulative-arc-length parametrization
        edges = np.roll(self.vertices, -1, axis=0) - self.vertices
        self._edge_lens = np.linalg.norm(edges, axis=1)
        self._perimeter = float(np.sum(self._edge_lens))
        self._cum_arc = np.concatenate(([0.0], np.cumsum(self._edge_lens)))

        # Cached bounding box
        self._bounds = (float(self.vertices[:, 0].min()),
                         float(self.vertices[:, 0].max()),
                         float(self.vertices[:, 1].min()),
                         float(self.vertices[:, 1].max()))

    @staticmethod
    def _signed_area(V):
        x = V[:, 0]
        y = V[:, 1]
        return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))

    @staticmethod
    def _area_centroid(V, A):
        x = V[:, 0]
        y = V[:, 1]
        x_next = np.roll(x, -1)
        y_next = np.roll(y, -1)
        cross = x * y_next - x_next * y
        cx = float(np.sum((x + x_next) * cross)) / (6.0 * A)
        cy = float(np.sum((y + y_next) * cross)) / (6.0 * A)
        return np.array([cx, cy])

    def _boundary_point(self, beta):
        """Interpolated boundary point B(beta) in global (u, v).

        beta in [0, 1] is parametrized by arc length along the closed
        polygon boundary. beta=0 -> V[0], beta=1 -> V[0] (wraps).
        Supports scalar or array inputs.
        """
        beta = np.asarray(beta, dtype=float)
        s = np.mod(beta, 1.0) * self._perimeter
        # Locate segment index via searchsorted on cumulative arc length
        idx = np.searchsorted(self._cum_arc, s, side='right') - 1
        idx = np.clip(idx, 0, self._edge_lens.size - 1)
        s0 = self._cum_arc[idx]
        t = (s - s0) / np.maximum(self._edge_lens[idx], 1e-30)
        V0 = self.vertices[idx]
        V1 = self.vertices[(idx + 1) % self.vertices.shape[0]]
        return V0 + t[..., None] * (V1 - V0)

    def sample_at(self, alpha, beta):
        alpha = np.asarray(alpha, dtype=float)
        beta = np.asarray(beta, dtype=float)
        B = self._boundary_point(beta)
        # sqrt-scaled radial blend
        s = np.sqrt(np.clip(alpha, 0.0, 1.0))
        u = self.centroid[0] + s * (B[..., 0] - self.centroid[0])
        v = self.centroid[1] + s * (B[..., 1] - self.centroid[1])
        return u, v

    def sample_perimeter(self, n):
        n = int(n)
        if n < 1:
            raise ValueError(f"n_peri must be >= 1, got {n}")
        betas = (np.arange(n) + 0.5) / n
        B = self._boundary_point(betas)
        return B[..., 0], B[..., 1]

    def total_area(self):
        return self._area

    def bounding_wh(self):
        umin, umax, vmin, vmax = self._bounds
        return umax - umin, vmax - vmin

    def _interpolate_same(self, other, t):
        """Same-type PolygonProfile interpolation works only when vertex
        counts match. Otherwise fall back to InterpolatedProfile."""
        if self.vertices.shape != other.vertices.shape:
            raise NotImplementedError(
                f"PolygonProfile._interpolate_same requires equal vertex "
                f"counts; got {self.vertices.shape[0]} vs "
                f"{other.vertices.shape[0]}")
        t = float(t)
        V_t = (1.0 - t) * self.vertices + t * other.vertices
        return PolygonProfile(V_t)

    def __repr__(self):
        return (f"PolygonProfile(n_vertices={self.vertices.shape[0]}, "
                f"area={self._area:g})")


class InterpolatedProfile(Profile):
    """Cross-type loft: linear interpolation between two profiles of
    potentially different types at path parameter t in [0, 1].

    sample_at(alpha, beta) = (1-t)*a.sample_at(..) + t*b.sample_at(..)

    total_area(t) linearly blends. This is NOT in general area-preserving
    at intermediate t (the intermediate shape can have a slightly different
    area than the blend of a.area and b.area), but the deviation is
    typically <5% for physically meaningful cross-type lofts (rect <->
    circle, etc.) and vanishes at the endpoints.

    For the PEEC filament model this approximation means:
        - at t=0 and t=1: exact profile, exact area
        - at intermediate t: approximate area used for R_k and cell (w, h)
        - the filament COORDINATES (u, v) are well-defined at every t
          because (alpha, beta) is shared.
    """

    def __init__(self, a, b, t):
        self.a = a
        self.b = b
        self.t = float(t)

    def sample_at(self, alpha, beta):
        ua, va = self.a.sample_at(alpha, beta)
        ub, vb = self.b.sample_at(alpha, beta)
        t = self.t
        return (1.0 - t) * ua + t * ub, (1.0 - t) * va + t * vb

    def sample_perimeter(self, n):
        ua, va = self.a.sample_perimeter(n)
        ub, vb = self.b.sample_perimeter(n)
        t = self.t
        return (1.0 - t) * ua + t * ub, (1.0 - t) * va + t * vb

    def total_area(self):
        t = self.t
        return (1.0 - t) * self.a.total_area() + t * self.b.total_area()

    def bounding_wh(self):
        wa, ha = self.a.bounding_wh()
        wb, hb = self.b.bounding_wh()
        t = self.t
        return ((1.0 - t) * wa + t * wb, (1.0 - t) * ha + t * hb)

    def _interpolate_same(self, other, t):
        raise NotImplementedError(
            "Cannot nest InterpolatedProfile; "
            "build a fresh InterpolatedProfile from the base profiles.")

    def __repr__(self):
        return (f"InterpolatedProfile(a={self.a!r}, b={self.b!r}, "
                 f"t={self.t:g})")
