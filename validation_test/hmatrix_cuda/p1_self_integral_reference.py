"""Independent star-simplex quadrature for P1 charge self integrals.

No Radia potentials or native modules are imported. This is a candidate
reference, not an adopted oracle. Refine inner and outer rules separately.
"""
import numpy as np


def interval(order):
    if not isinstance(order, int) or order < 2:
        raise ValueError('Quadrature order must be an integer >= 2')
    x, w = np.polynomial.legendre.leggauss(order)
    return (x + 1) / 2, w / 2


def simplex_rule(order, dimension):
    q, w = interval(order)
    if dimension == 2:
        s, t = np.meshgrid(q, q, indexing='ij')
        ws, wt = np.meshgrid(w, w, indexing='ij')
        return np.stack((s, (1-s)*t), axis=-1).reshape(-1, 2), (ws*wt*(1-s)).ravel()
    if dimension == 3:
        s, t, u = np.meshgrid(q, q, q, indexing='ij')
        ws, wt, wu = np.meshgrid(w, w, w, indexing='ij')
        return np.stack((s, (1-s)*t, (1-s)*(1-t)*u), axis=-1).reshape(-1, 3), (ws*wt*wu*(1-s)**2*(1-t)).ravel()
    raise ValueError('Only triangles and tetrahedra are supported')


def self_integrals(vertices, outer_order, inner_order):
    """Return constant TET or TRI (1, eta, xi) self energies / (4*pi).

    For each outer x, partition the source simplex around x. The radial
    variable cancels the 1/r singularity: TRI integrates an affine charge
    to (m(x)+m(edge))/2; constant TET integrates radial u to 1/2.
    All remaining integrals use open Gauss nodes, with no cutoff/drop rule.
    """
    vertices = np.asarray(vertices, dtype=float)
    if vertices.shape not in ((3, 3), (4, 3)) or not np.isfinite(vertices).all():
        raise ValueError('Expected finite TRI or TET vertices')
    dim = len(vertices)-1
    edges = vertices[1:] - vertices[0]
    measure = (np.linalg.norm(np.cross(*edges)) if dim == 2
               else abs(np.linalg.det(edges)))
    if not np.isfinite(measure) or measure <= 0:
        raise ValueError('Degenerate simplex')
    outer, ow = simplex_rule(outer_order, dim)
    x = vertices[0] + outer @ edges
    result = np.zeros(3 if dim == 2 else 1)
    if dim == 2:
        v, vw = interval(inner_order)
        corner = np.array([[1., 0., 0.], [1., 0., 1.], [1., 1., 0.]])
        modes = np.column_stack((np.ones(len(x)), outer[:, 1], outer[:, 0]))
        for i, j in ((0, 1), (1, 2), (2, 0)):
            a, b = vertices[i]-x, vertices[j]-x
            r = a[:, None, :]*(1-v)[None, :, None] + b[:, None, :]*v[None, :, None]
            distance = np.linalg.norm(r, axis=-1)
            if np.any(distance <= 0) or not np.isfinite(distance).all():
                raise ValueError('Invalid transformed distance')
            edge_modes = corner[i]*(1-v[:, None]) + corner[j]*v[:, None]
            inner = np.einsum('nv,nvk,v->nk', 1/distance,
                              (modes[:, None, :]+edge_modes[None, :, :])/2, vw)
            result += np.sum((ow*np.linalg.norm(np.cross(a, b), axis=1))[:, None]*modes*inner, axis=0)
    else:
        uv, iw = simplex_rule(inner_order, 2)
        for start in range(0, len(x), 32):
            batch = x[start:start+32]
            potential = np.zeros(len(batch))
            for i, j, k in ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)):
                a, b, c = vertices[i]-batch, vertices[j]-batch, vertices[k]-batch
                r = a[:, None, :] + uv[None, :, :1]*(b-a)[:, None, :] + uv[None, :, 1:]*(c-a)[:, None, :]
                distance = np.linalg.norm(r, axis=-1)
                if np.any(distance <= 0) or not np.isfinite(distance).all():
                    raise ValueError('Invalid transformed distance')
                determinant = np.abs(np.einsum('ij,ij->i', a, np.cross(b, c)))
                potential += determinant/2 * ((1/distance) @ iw)
            result[0] += ow[start:start+32] @ potential
    result *= measure/(4*np.pi)
    if not np.isfinite(result).all() or np.any(result <= 0):
        raise ValueError('Invalid self integral')
    return result
