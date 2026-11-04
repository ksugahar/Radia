"""Coulomb integrals of polynomial charges on a box, by a route that shares nothing with the charge kernels.

A charge is a product of one monomial per axis, on an interval or at a point. The kernel 1/r depends on differences
only, so each axis reduces to the weight W(u) = integral a(x) b(x - u) dx, and the pair integral to integrals of
W1(u) W2(v) W3(w)/r over boxes. On a box with the origin at a corner the map (u, u s, u t) per dominant axis
(Duffy, SIAM J. Numer. Anal. 19, 1982) removes the singularity: the integrand is analytic and Gauss-Legendre converges
geometrically. Closed forms of the box integrals of Bailey, Borwein and Crandall (J. Comput. Appl. Math. 206, 2007)
check the reduction, and a uniform magnetization of a cube has the exact demagnetizing factor 1/3.

A polynomial magnetization has no jump across interior faces, so on any mesh of the cube its charge is the volume
term -div M and the surface term M.n of the cube itself. The rational pyramid charge is outside this reference.
"""
from itertools import product

import numpy as np

GAUSS = 24


def _gauss01(n):
    nodes, weights = np.polynomial.legendre.leggauss(n)
    return .5*(nodes + 1), .5*weights


def _split(low, high, weight):
    if low < 0. < high:
        return [(low, 0., weight), (0., high, weight)]
    return [(low, high, weight)]


def axis_weight(a, b):
    """Pieces (low, high, W) of u = x - x'; low == high marks a point mass. A factor is (lo, hi, exponent)."""
    (a_lo, a_hi, a_e), (b_lo, b_hi, b_e) = a, b
    if a_lo == a_hi and b_lo == b_hi:
        return [(a_lo - b_lo, a_lo - b_lo, lambda u: np.ones_like(u))]
    if b_lo == b_hi:
        return _split(a_lo - b_lo, a_hi - b_lo, lambda u: (u + b_lo)**a_e)
    if a_lo == a_hi:
        return _split(a_lo - b_hi, a_lo - b_lo, lambda u: (a_lo - u)**b_e)
    nodes, weights = _gauss01((a_e + b_e)//2 + 1)

    def weight(u):
        low = np.maximum(a_lo, b_lo + u)[..., None]
        length = np.minimum(a_hi, b_hi + u)[..., None] - low
        x = low + length*nodes
        return (length*weights*x**a_e*(x - u[..., None])**b_e).sum(-1)
    breaks = sorted({a_lo - b_hi, a_lo - b_lo, a_hi - b_hi, a_hi - b_lo, 0.})
    breaks = [value for value in breaks if a_lo - b_hi <= value <= a_hi - b_lo]
    return [(low, high, weight) for low, high in zip(breaks, breaks[1:])]


def _corner_box(pieces, free, n):
    """Origin at a corner of the box: one transformed pyramid per free axis."""
    x, w = _gauss01(n)
    total = 0.
    lengths = {i: pieces[i][1] - pieces[i][0] for i in free}
    signs = {i: (1. if pieces[i][0] == 0. else -1.) for i in free}
    for dominant in free:
        others = [i for i in free if i != dominant]
        grids = np.meshgrid(*([x]*len(free)), indexing="ij")
        weight = np.ones_like(grids[0])
        for axis in range(len(free)):
            shape = [1]*len(free)
            shape[axis] = n
            weight = weight*w.reshape(shape)
        xi = lengths[dominant]*grids[0]
        coordinates = {dominant: xi}
        radial = np.ones_like(xi)
        scale = lengths[dominant]
        for i, s in zip(others, grids[1:]):
            ratio = lengths[i]/lengths[dominant]
            coordinates[i] = ratio*xi*s
            radial = radial + (ratio*s)**2
            scale = scale*ratio
        value = weight*scale*xi**(len(free) - 2)/np.sqrt(radial)
        for axis, (low, _, axis_function) in enumerate(pieces):
            argument = signs[axis]*coordinates[axis] if axis in free else np.full_like(xi, low)
            value = value*axis_function(argument)
        total += float(value.sum())
    return total


def _tensor_box(pieces, free, n):
    x, w = _gauss01(n)
    grids = np.meshgrid(*[pieces[i][0] + (pieces[i][1] - pieces[i][0])*x for i in free], indexing="ij")
    weight = np.ones_like(grids[0])
    for axis, i in enumerate(free):
        shape = [1]*len(free)
        shape[axis] = n
        weight = weight*((pieces[i][1] - pieces[i][0])*w).reshape(shape)
    coordinates = {i: grid for i, grid in zip(free, grids)}
    squared = np.zeros_like(weight)
    value = weight
    for axis, (low, _, axis_function) in enumerate(pieces):
        argument = coordinates[axis] if axis in free else np.full_like(weight, low)
        squared = squared + argument**2
        value = value*axis_function(argument)
    return float((value/np.sqrt(squared)).sum())


def coulomb(a, b, n=GAUSS):
    """Integral of a(x) b(y)/|x - y| for two charges, each a tuple of three factors."""
    total = 0.
    for pieces in product(*(axis_weight(x, y) for x, y in zip(a, b))):
        free = [i for i, (low, high, _) in enumerate(pieces) if high > low]
        if all(low == 0. or high == 0. for low, high, _ in pieces):
            if len(free) < 2:
                raise ValueError("Line charges meet on a line: the Coulomb integral diverges")
            total += _corner_box(pieces, free, n)
            continue
        nearest = np.sqrt(sum(min(abs(low), abs(high))**2 if low*high > 0 or low == high else 0. for low, high, _ in pieces))
        size = np.sqrt(sum((high - low)**2 for low, high, _ in pieces))
        if nearest < .2*size:
            raise ValueError(f"Box is near the singularity without touching it: distance {nearest}, size {size}")
        total += _tensor_box(pieces, free, n)
    return total


def magnetization_modes(order):
    """Monomials of total degree <= order in the centered coordinates of the cube, per component."""
    return [(component, (a, b, c)) for component in range(3) for a in range(order + 1)
            for b in range(order + 1 - a) for c in range(order + 1 - a - b)]


def mode_charges(mode, half=.5):
    """Terms (coefficient, factors) of the volume charge -div M and the surface charge M.n."""
    component, exponents = mode
    terms = []
    if exponents[component]:
        lowered = list(exponents)
        lowered[component] -= 1
        terms.append((-float(exponents[component]), tuple((-half, half, e) for e in lowered)))
    for side in (-1., 1.):
        factors = [(-half, half, e) for e in exponents]
        factors[component] = (side*half, side*half, 0)
        terms.append((side*(side*half)**exponents[component], tuple(factors)))
    return terms


def mode_energies(modes, half=.5, n=GAUSS):
    """Q_ij = integral of rho_i(x) rho_j(y)/(4 pi |x - y|) over all charges of two modes."""
    charges = [mode_charges(mode, half) for mode in modes]
    cache = {}

    def integral(a, b):
        key = (a, b) if a <= b else (b, a)
        if key not in cache:
            cache[key] = coulomb(key[0], key[1], n)
        return cache[key]
    energies = np.zeros((len(modes), len(modes)))
    for i, first in enumerate(charges):
        for j in range(i, len(modes)):
            energies[i, j] = energies[j, i] = sum(c*d*integral(a, b) for c, a in first for d, b in charges[j])/(4*np.pi)
    return energies


def check_reference(n=GAUSS, limit=1e-12):
    """Actual against closed-form values of the reduction; raises when one leaves the limit."""
    interval = (-.5, .5, 0)
    face = (.5, .5, 0)
    root2, root3 = np.sqrt(2), np.sqrt(3)
    closed = dict(cube=(.4*(1 + root2 - 2*root3) - 2*np.pi/3 + 2*np.log(1 + root2) + 4*np.log((1 + root3)/root2),
                        coulomb((interval,)*3, (interval,)*3, n)),
                  square=(4*np.log(1 + root2) - (4/3)*(root2 - 1), coulomb((face, interval, interval), (face, interval, interval), n)),
                  demagnetizing_factor=(1/3, float(mode_energies([(0, (0, 0, 0))], n=n)[0, 0])))
    for name, (expected, actual) in closed.items():
        if abs(actual - expected) > limit*abs(expected):
            raise RuntimeError(f"Polynomial reference {name}: actual {actual!r}, expected {expected!r}, limit {limit}")
    return {name: dict(expected=float(expected), actual=float(actual)) for name, (expected, actual) in closed.items()}
