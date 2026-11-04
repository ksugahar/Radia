"""Coulomb energies of the rational charge x y/(1-z) on one pyramid, by a reduction independent of the charge Gram.

On the pyramid with the reference coordinates as physical ones, y = (t a, t b, 1-t) with (a, b) in the unit square and
measure t**2 dt da db, so the charge is t a b. Writing the lower of two heights t' = s t removes t in closed form, and
the two axes separate into weights W(u; s) = integral of f(u + s c) g(c) dc over the c that keep u + s c in [0, 1]:

    E = integral of s**k ds, integral of W1(u) W1(v) / sqrt(u**2 + v**2 + (1-s)**2) du dv, times a power of t's integral.

W is piecewise polynomial in u with breakpoints -s, 0, 1-s, 1. The rest is Gauss quadrature on panels shrinking to the
peak at u = v = 0, whose width is 1-s, and on s graded to 1. This route shares no code or rule with the charge Gram.
"""
import numpy as np

CONSTANT = 0
LINEAR = 1


def _gauss01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return .5*(x + 1), .5*w


def _weight(u, s, outer, inner, rule=_gauss01(4)):
    """W(u; s) with f(x) = x**outer and g(c) = c**inner; four points are exact on these polynomials."""
    low, high = np.maximum(0., -u/s), np.minimum(1., (1. - u)/s)
    width = np.maximum(high - low, 0.)
    c = low[:, None] + width[:, None]*rule[0]
    return width*np.sum(rule[1]*(u[:, None] + s*c)**outer*c**inner, axis=1)


def _axis(s, panel_rule):
    """Points and weights of one axis on [-s, 1]: panels split at the breakpoints and doubling away from the peak at 0."""
    gap = 1. - s
    edges = {-s, 0., gap, 1.}
    size = gap
    while size < 1.:
        edges |= {size} | ({-size} if size < s else set())
        size *= 2.
    edges = np.array(sorted(edges))
    x = (edges[:-1, None] + np.diff(edges)[:, None]*panel_rule[0]).ravel()
    w = (np.diff(edges)[:, None]*panel_rule[1]).ravel()
    return x, w


def height_integral(power, weights_of_s, *, heights=40, panel=14):
    """Integral over s in [0, 1] of s**power times the separated double integral, graded toward s = 1."""
    r, rw = _gauss01(heights)
    s_values, s_weights = 1. - r**4, rw*4.*r**3
    panel_rule = _gauss01(panel)
    total = 0.
    for s, weight in zip(s_values, s_weights):
        x, w = _axis(s, panel_rule)
        first, second = weights_of_s(x, s)
        gap2 = (1. - s)**2
        kernel = 1./np.sqrt(x[:, None]**2 + x[None, :]**2 + gap2)
        total += weight*s**power*float((w*first) @ kernel @ (w*second))
    return total


def energies(*, heights=40, panel=14):
    """Double integrals of the constant and rational charges, without the 1/(4 pi) of the Gram."""
    both = lambda outer, inner: (lambda x, s: (_weight(x, s, outer, inner), _weight(x, s, outer, inner)))
    constant = 2./5.*height_integral(2, both(0, 0), heights=heights, panel=panel)
    rational = 2./7.*height_integral(3, both(1, 1), heights=heights, panel=panel)
    # Rational at the higher point, then at the lower one.
    mixed = (height_integral(2, both(1, 0), heights=heights, panel=panel)
             + height_integral(3, both(0, 1), heights=heights, panel=panel))/6.
    return dict(constant=constant, rational=rational, mixed=mixed)


def separated_integral(s, *, panel=14):
    """The inner double integral of the rational pair at one s, for a check against a four-dimensional integration."""
    x, w = _axis(s, _gauss01(panel))
    first = _weight(x, s, 1, 1)
    kernel = 1./np.sqrt(x[:, None]**2 + x[None, :]**2 + (1. - s)**2)
    return float((w*first) @ kernel @ (w*first))
