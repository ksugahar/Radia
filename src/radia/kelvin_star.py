"""Kelvin inversion of an exterior bounded by a star-shaped surface.

The map ``y = T(x) = rho(n)^2 x / r^2`` about the centre, with ``r = |x|``, ``n = x/r`` and ``rho(n)`` the boundary
radius along ``n``, fixes the boundary pointwise and is its own inverse. It is the spherical inversion when ``rho`` is
constant. With the boundary normal ``v``, set ``t = n - v/(v.n) = r grad(rho)/rho``, which is tangent (``n.t = 0``,
``|t| = tan(angle(n, v))``), and ``w = n - t``. Then

    J = (rho/r)^2 A,   A = I - 2 n w^T,   A^2 = I,   det A = -1,

so ``T`` reverses orientation and the exterior material is a tensor, not a scalar:

    nu' = nu_0 J^T J / |det J| = nu_0 (r/rho)^2 (I - 2 (n t^T + t n^T) + 4 t t^T)     (curl-curl form)
    mu' = mu_0 J^-1 J^-T |det J| = mu_0 (rho/r)^2 (I + 2 (n t^T + t n^T) + 4 |t|^2 n n^T)   (scalar potential form)

The two are exact inverses. The scalar potential is a twisted 0-form, so with ``s = sign(det J) = -1``

    Omega'(x) = s Omega(T x),   H'(x) = s J^T H(T x) = -(rho/r)^2 (I - 2 w n^T) H(T x).

Across the identified boundary ``t.H' = -t.H`` for a tangent ``t`` and ``v.mu'.H' = v.H`` with ``v`` outward on both
copies, the conditions of the periodic scalar-potential space, whose exterior unknown is ``Omega o T``. The condition number of ``mu'`` and ``nu'`` is
``(tan a + sqrt(1 + tan^2 a))^4``, about ``16 tan^4 a``, where ``a`` is the angle between ``n`` and ``v``: the boundary
shape sets the conditioning of the exterior block.

Each boundary class gives ``radius_and_normal`` for coefficient functions and for NumPy arrays of unit directions, so
one expression serves the assembly and the checks. The boundary must be star-shaped about the centre and C1.
"""
from __future__ import annotations

import math

import numpy as np

from .kelvin_material import MU_0, NU_0, _is_kelvin_material


def _checked_positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite; got {value!r}")
    return value


class SphereBoundary:
    """The sphere of radius ``radius``: the tensors reduce to the isotropic Kelvin factors."""

    def __init__(self, radius):
        self.radius = _checked_positive(radius, "radius")

    def radius_and_normal(self, direction, module=None):
        if module is None:
            return self.radius*np.ones(np.shape(direction)[:-1]), np.asarray(direction, dtype=float)
        return module.CF(self.radius), direction


class RoundedDiskBoundary:
    """Two flat faces ``z = +-edge_radius`` for ``hypot(x, y) <= flat_radius`` joined by a half-torus edge of ``edge_radius``."""

    def __init__(self, flat_radius, edge_radius):
        self.flat_radius = _checked_positive(flat_radius, "flat_radius")
        self.edge_radius = _checked_positive(edge_radius, "edge_radius")

    @property
    def maximum_tangent(self):
        """Largest ``|t|``, at the rim where the flat face meets the edge."""
        return self.flat_radius/self.edge_radius

    def radius_and_normal(self, direction, module=None):
        flat, edge = self.flat_radius, self.edge_radius
        if module is None:
            n = np.asarray(direction, dtype=float)
            s = np.hypot(n[..., 0], n[..., 1])
            c = np.abs(n[..., 2])
            edge_rho = flat*s + np.sqrt(np.maximum(flat**2*s**2 - flat**2 + edge**2, 0.0))
            on_face = flat*c - edge*s > 0.0
            plane = np.stack([n[..., 0], n[..., 1], np.zeros_like(s)], axis=-1)/np.where(s > 1e-12, s, 1.0)[..., None]
            edge_normal = (edge_rho[..., None]*n - flat*plane)/edge
            face_normal = np.stack([np.zeros_like(s), np.zeros_like(s), np.sign(n[..., 2])], axis=-1)
            rho = np.where(on_face, edge/np.where(c > 1e-12, c, 1.0), edge_rho)
            return rho, np.where(on_face[..., None], face_normal, edge_normal)
        ng = module
        s = ng.sqrt(direction[0]*direction[0] + direction[1]*direction[1])
        c = ng.sqrt(direction[2]*direction[2])
        safe_s = ng.IfPos(s - 1e-12, s, 1.0)
        safe_c = ng.IfPos(c - 1e-12, c, 1.0)
        discriminant = flat**2*s*s - flat**2 + edge**2
        edge_rho = flat*s + ng.sqrt(ng.IfPos(discriminant, discriminant, 0.0))
        on_face = flat*c - edge*s
        plane = ng.CF((direction[0]/safe_s, direction[1]/safe_s, 0.0))
        edge_normal = (edge_rho*direction - flat*plane)/edge
        face_normal = ng.CF((0.0, 0.0, direction[2]/safe_c))
        return ng.IfPos(on_face, edge/safe_c, edge_rho), ng.IfPos(on_face, face_normal, edge_normal)


def star_kelvin_map(points, boundary, center):
    """``T`` and its own inverse: an exterior physical point to its computational point and back."""
    center = np.asarray(center, dtype=float)
    relative = np.asarray(points, dtype=float) - center
    r = np.linalg.norm(relative, axis=-1, keepdims=True)
    if np.any(r <= 0.0):
        raise ValueError("The Kelvin map is singular at its centre")
    rho, _ = boundary.radius_and_normal(relative/r)
    return center + (np.asarray(rho)[..., None]/r)**2*relative


def _frame(relative, boundary):
    r = np.linalg.norm(relative, axis=-1, keepdims=True)
    if np.any(r <= 0.0):
        raise ValueError("The Kelvin frame is singular at its centre")
    n = relative/r
    rho, normal = boundary.radius_and_normal(n)
    tangent = n - normal/np.sum(normal*n, axis=-1, keepdims=True)
    return np.asarray(rho)[..., None]/r, n, tangent


def star_kelvin_jacobian(points, boundary, center):
    """``J = (rho/r)^2 (I - 2 n w^T)`` at computational points, shape ``(..., 3, 3)``."""
    ratio, n, tangent = _frame(np.asarray(points, dtype=float) - np.asarray(center, dtype=float), boundary)
    return (ratio**2)[..., None]*(np.eye(3) - 2.0*n[..., None]*(n - tangent)[..., None, :])


def star_kelvin_nu(points, boundary, center, nu_0=NU_0):
    """The curl-curl tensor ``nu' = nu_0 J^T J / |det J|``, shape ``(..., 3, 3)``."""
    ratio, n, t = _frame(np.asarray(points, dtype=float) - np.asarray(center, dtype=float), boundary)
    cross = n[..., None]*t[..., None, :]
    return nu_0*(ratio**-2)[..., None]*(np.eye(3) - 2.0*(cross + np.swapaxes(cross, -1, -2)) + 4.0*t[..., None]*t[..., None, :])


def star_kelvin_mu(points, boundary, center, mu_0=MU_0):
    """The scalar-potential tensor ``mu' = mu_0 J^-1 J^-T |det J|``, shape ``(..., 3, 3)``."""
    ratio, n, t = _frame(np.asarray(points, dtype=float) - np.asarray(center, dtype=float), boundary)
    cross = n[..., None]*t[..., None, :]
    outer = np.sum(t*t, axis=-1)[..., None, None]*n[..., None]*n[..., None, :]
    return mu_0*(ratio**2)[..., None]*(np.eye(3) + 2.0*(cross + np.swapaxes(cross, -1, -2)) + 4.0*outer)


def star_kelvin_physical_h(h_computational, points, boundary, center):
    """The physical ``H(T x)`` from the computational ``H'(x) = s J^T H(T x)`` at the computational ``points``."""
    jacobian = star_kelvin_jacobian(points, boundary, center)
    return -np.linalg.solve(np.swapaxes(jacobian, -1, -2), np.asarray(h_computational, dtype=float)[..., None])[..., 0]


def star_kelvin_mu_cf(mesh, boundary, center, mu_0=MU_0, kelvin_mats=("kelvin",), mu_r_by_material=None, kelvin_match_exact=False):
    """Fully Kelvin-aware permeability tensor for the ``mu_cf`` argument of the mixed Omega solvers.

    Kelvin materials take ``mu'`` about ``center``; the others take ``mu_0 mu_r I``. Unknown material names and a
    relative permeability on a Kelvin material refuse, as in ``make_kelvin_mu_cf``.
    """
    import ngsolve as ng
    ox, oy, oz = (float(value) for value in center)
    relative = ng.CF((ng.x - ox, ng.y - oy, ng.z - oz))
    r_squared = relative*relative + 1e-24
    n = relative/ng.sqrt(r_squared)
    rho, normal = boundary.radius_and_normal(n, ng)
    t = n - normal/(normal*n)
    cross = ng.OuterProduct(n, t)
    exterior = mu_0*(rho*rho/r_squared)*(ng.Id(3) + 2.0*(cross + cross.trans) + 4.0*(t*t)*ng.OuterProduct(n, n))
    names = tuple(mesh.GetMaterials())
    relative_mu = dict(mu_r_by_material or {})
    unknown = sorted(set(relative_mu) - set(names))
    if unknown:
        raise ValueError(f"mu_r_by_material names are not mesh materials: {unknown}; available: {sorted(set(names))}")
    values = {}
    for name in names:
        if _is_kelvin_material(name, kelvin_mats, exact=kelvin_match_exact):
            if name in relative_mu:
                raise ValueError(f"mu_r_by_material must not override Kelvin material {name!r}")
            values[name] = exterior
            continue
        values[name] = mu_0*_checked_positive(relative_mu.get(name, 1.0), f"mu_r_by_material[{name!r}]")*ng.Id(3)
    return mesh.MaterialCF(values, default=mu_0*ng.Id(3))


def star_kelvin_source_h_cf(field, boundary, center):
    """``H'(x) = s J^T H(T x)`` of an analytic source for a reduced exterior.

    ``field(X, Y, Z)`` returns the physical H as a 3-vector coefficient function of its coordinate arguments, so the
    composition with ``T`` stays exact on the quadrature points.
    """
    import ngsolve as ng
    ox, oy, oz = (float(value) for value in center)
    relative = ng.CF((ng.x - ox, ng.y - oy, ng.z - oz))
    r_squared = relative*relative + 1e-24
    n = relative/ng.sqrt(r_squared)
    rho, normal = boundary.radius_and_normal(n, ng)
    w = normal/(normal*n)
    ratio = rho*rho/r_squared
    image = field(ox + ratio*relative[0], oy + ratio*relative[1], oz + ratio*relative[2])
    return -ratio*(ng.Id(3) - 2.0*ng.OuterProduct(w, n))*image
