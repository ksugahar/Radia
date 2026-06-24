"""Readable P1 tetrahedron FEM element formulas for scalar 3D Poisson problems.

These local matrices are the 3D counterpart of :mod:`scalar_fem2d`: constant
shape-function gradients, exact stiffness, consistent mass, and constant-source
load for first-order tetrahedra.  They are intentionally small and dependency
free so the same formulas can be mirrored in teaching MATLAB prototypes.
"""


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _solve3(rows, rhs):
    """Solve a 3x3 system with row-major ``rows`` using Cramer's rule."""
    c0 = (rows[0][0], rows[1][0], rows[2][0])
    c1 = (rows[0][1], rows[1][1], rows[2][1])
    c2 = (rows[0][2], rows[1][2], rows[2][2])
    det = _dot(c0, _cross(c1, c2))
    if det == 0.0:
        raise ValueError("singular 3x3 system")
    return (
        _dot(rhs, _cross(c1, c2)) / det,
        _dot(c0, _cross(rhs, c2)) / det,
        _dot(c0, _cross(c1, rhs)) / det,
    )


def p1_tetrahedron_geometry(vertices):
    """Volume and constant P1 shape-function gradients for a tetrahedron.

    ``vertices`` is ``[(x0,y0,z0), ..., (x3,y3,z3)]``.  Returns
    ``{"volume": V, "gradients": [(dN0/dx,dN0/dy,dN0/dz), ...]}``.
    The returned gradient order follows the supplied local node order.
    """
    if len(vertices) != 4:
        raise ValueError("a tetrahedron needs exactly four vertices")
    v = [tuple(float(x) for x in p) for p in vertices]
    if any(len(p) != 3 for p in v):
        raise ValueError("tetrahedron vertices must be 3D points")

    j0 = _sub(v[1], v[0])
    j1 = _sub(v[2], v[0])
    j2 = _sub(v[3], v[0])
    detj = _dot(j0, _cross(j1, j2))
    if detj == 0.0:
        raise ValueError("degenerate tetrahedron")

    jt_rows = (j0, j1, j2)
    g1 = _solve3(jt_rows, (1.0, 0.0, 0.0))
    g2 = _solve3(jt_rows, (0.0, 1.0, 0.0))
    g3 = _solve3(jt_rows, (0.0, 0.0, 1.0))
    g0 = (-(g1[0] + g2[0] + g3[0]),
          -(g1[1] + g2[1] + g3[1]),
          -(g1[2] + g2[2] + g3[2]))
    return {"volume": abs(detj) / 6.0, "gradients": [g0, g1, g2, g3]}


def p1_tetrahedron_stiffness(vertices, coeff=1.0):
    """Local scalar P1 stiffness matrix for ``-div(coeff grad u)`` on a tetrahedron."""
    g = p1_tetrahedron_geometry(vertices)
    volume = g["volume"]
    grads = g["gradients"]
    c = float(coeff)
    return [[c * volume * _dot(gi, gj) for gj in grads] for gi in grads]


def p1_tetrahedron_mass(vertices, density=1.0):
    """Consistent scalar P1 mass matrix ``int density N_i N_j dV`` on a tetrahedron."""
    volume = p1_tetrahedron_geometry(vertices)["volume"]
    d = float(density)
    return [[d * volume * (2.0 if i == j else 1.0) / 20.0 for j in range(4)] for i in range(4)]


def p1_tetrahedron_constant_load(vertices, source=1.0):
    """Local P1 load vector for a constant source term over a tetrahedron."""
    volume = p1_tetrahedron_geometry(vertices)["volume"]
    return [float(source) * volume / 4.0] * 4
