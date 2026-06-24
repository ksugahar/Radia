"""Readable P1 tetrahedron FEM element formulas for scalar 3D Poisson problems.

These local matrices are the 3D counterpart of :mod:`scalar_fem2d`: constant
shape-function gradients, exact stiffness, consistent mass, and constant-source
load for first-order tetrahedra.  They are intentionally small and dependency
free so the same formulas can be mirrored in teaching MATLAB prototypes.
"""

import math


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


def _scale(a, s):
    return (s * a[0], s * a[1], s * a[2])


def _norm(a):
    return math.sqrt(_dot(a, a))


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


def p1_tetrahedron_gradient(vertices, nodal_values):
    """Constant gradient of a scalar P1 tetrahedron field.

    ``nodal_values`` follows the same local node order as ``vertices``.  The
    result is ``sum_i u_i grad(N_i)``.  This is the local post-processing block
    students usually need immediately after solving a readable P1 system.
    """
    values = [float(v) for v in nodal_values]
    if len(values) != 4:
        raise ValueError("nodal_values must contain four values")
    grads = p1_tetrahedron_geometry(vertices)["gradients"]
    return tuple(sum(values[i] * grads[i][axis] for i in range(4)) for axis in range(3))


def p1_tetrahedron_flux(vertices, nodal_values, coeff=1.0):
    """Constant physical flux ``-coeff grad(u)`` for a scalar P1 tetrahedron."""
    grad_u = p1_tetrahedron_gradient(vertices, nodal_values)
    c = float(coeff)
    return tuple(-c * value for value in grad_u)


def p1_tetrahedron_boundary_fluxes(vertices, nodal_values, coeff=1.0):
    """Outward Neumann flux rows on the four faces of a P1 tetrahedron.

    Each row reports the local face opposite one tetrahedron node, its outward
    area vector, flux density ``q.n`` and integrated flux ``int_face q.n dS``
    for ``q=-coeff grad(u)``.  The face-node ids are zero-based local ids so
    they can be mapped directly to a local element or to one-based ``.vol``
    nodes by the caller.
    """
    v = [tuple(float(x) for x in p) for p in vertices]
    if len(v) != 4 or any(len(p) != 3 for p in v):
        raise ValueError("a tetrahedron needs exactly four 3D vertices")
    q = p1_tetrahedron_flux(v, nodal_values, coeff=coeff)
    rows = []
    for opposite in range(4):
        face = tuple(i for i in range(4) if i != opposite)
        a, b, cpt = (v[i] for i in face)
        area_vector = _scale(_cross(_sub(b, a), _sub(cpt, a)), 0.5)
        face_centroid = tuple((a[axis] + b[axis] + cpt[axis]) / 3.0 for axis in range(3))
        to_opposite = _sub(v[opposite], face_centroid)
        if _dot(area_vector, to_opposite) > 0.0:
            area_vector = _scale(area_vector, -1.0)
        area_vector = tuple(0.0 if component == 0.0 else component for component in area_vector)
        area = _norm(area_vector)
        if area == 0.0:
            raise ValueError("degenerate tetrahedron face")
        integrated = _dot(q, area_vector)
        rows.append({
            "opposite_local_node": opposite,
            "face_local_nodes": face,
            "area": area,
            "outward_area_vector": area_vector,
            "normal_flux_density": integrated / area,
            "integrated_flux": integrated,
        })
    return rows


def p1_surface_triangle_geometry(vertices):
    """Area, oriented normal, and surface gradients for a 3D P1 triangle.

    ``vertices`` is ``[(x0,y0,z0), (x1,y1,z1), (x2,y2,z2)]``.  Returns
    ``{"area", "area_vector", "unit_normal", "gradients"}``, where each
    gradient is the constant tangential gradient of the corresponding P1 shape
    function.  This is the boundary-triangle counterpart of the tetrahedron
    formulas and is the small clean-room block behind readable SurfaceL2/BEM
    assembly.
    """
    if len(vertices) != 3:
        raise ValueError("a surface triangle needs exactly three vertices")
    a, b, c = [tuple(float(x) for x in p) for p in vertices]
    if any(len(p) != 3 for p in (a, b, c)):
        raise ValueError("surface triangle vertices must be 3D points")

    normal2 = _cross(_sub(b, a), _sub(c, a))
    double_area = _norm(normal2)
    if double_area == 0.0:
        raise ValueError("degenerate surface triangle")
    area = 0.5 * double_area
    unit_normal = _scale(normal2, 1.0 / double_area)
    gradients = [
        _scale(_cross(unit_normal, _sub(c, b)), 1.0 / double_area),
        _scale(_cross(unit_normal, _sub(a, c)), 1.0 / double_area),
        _scale(_cross(unit_normal, _sub(b, a)), 1.0 / double_area),
    ]
    return {
        "area": area,
        "area_vector": _scale(normal2, 0.5),
        "unit_normal": unit_normal,
        "gradients": gradients,
    }


def p1_surface_triangle_stiffness(vertices, coeff=1.0):
    """Local P1 surface stiffness ``int coeff grad_s N_i . grad_s N_j dS``."""
    g = p1_surface_triangle_geometry(vertices)
    area = g["area"]
    grads = g["gradients"]
    c = float(coeff)
    return [[c * area * _dot(gi, gj) for gj in grads] for gi in grads]


def p1_surface_triangle_mass(vertices, density=1.0):
    """Consistent P1 surface mass matrix ``int density N_i N_j dS``."""
    area = p1_surface_triangle_geometry(vertices)["area"]
    d = float(density)
    return [[d * area * (2.0 if i == j else 1.0) / 12.0 for j in range(3)] for i in range(3)]


def p1_surface_triangle_constant_load(vertices, source=1.0):
    """Local P1 surface load vector for a constant source over a triangle."""
    area = p1_surface_triangle_geometry(vertices)["area"]
    return [float(source) * area / 3.0] * 3


def assemble_p1_tet_robin_system(points, tetrahedra, surface_triangles,
                                  volume_coeff=1.0, source=0.0,
                                  robin_coeff=0.0, boundary_flux=0.0):
    """Assemble a readable dense P1 tet system with boundary Robin terms.

    This is the small clean-room volume/surface trace system that MATLAB or
    Gypsilab-style teaching scripts can mirror directly.  Inputs use one-based
    element node ids, so a parsed Netgen ``.vol`` mesh can be passed without
    renumbering:

    ``tetrahedra`` may contain plain ``(n1,n2,n3,n4)`` tuples or records with
    ``nodes`` and optional ``matnr`` attributes. ``surface_triangles`` may
    contain plain ``(n1,n2,n3)`` tuples or records with ``nodes`` and optional
    ``bcnr`` attributes.

    The weak form is

        int_Omega k grad(u).grad(v) dV + int_Gamma r u v dS
        = int_Omega f v dV + int_Gamma g v dS.

    ``volume_coeff`` and ``source`` may be scalars or dictionaries keyed by
    material number.  ``robin_coeff`` and ``boundary_flux`` may be scalars or
    dictionaries keyed by boundary number.  The returned matrix is dense for
    readability; ``matrix_triplets`` gives one-based sparse triplets for MATLAB.
    """
    pts = [tuple(float(x) for x in point) for point in points]
    if any(len(point) != 3 for point in pts):
        raise ValueError("points must be 3D coordinates")
    n = len(pts)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    rhs = [0.0 for _ in range(n)]
    volume_by_material: dict[int, float] = {}
    boundary_area_by_number: dict[int, float] = {}
    robin_area_weight = 0.0
    flux_area_weight = 0.0

    for tet in tetrahedra:
        nodes = _one_based_nodes(tet, 4)
        matnr = int(getattr(tet, "matnr", 0))
        ids = _zero_based_ids(nodes, n)
        coords = [pts[i] for i in ids]
        coeff = _coefficient_value(volume_coeff, matnr, default=1.0)
        src = _coefficient_value(source, matnr, default=0.0)
        kloc = p1_tetrahedron_stiffness(coords, coeff=coeff)
        floc = p1_tetrahedron_constant_load(coords, source=src)
        volume = p1_tetrahedron_geometry(coords)["volume"]
        volume_by_material[matnr] = volume_by_material.get(matnr, 0.0) + volume
        for a, ia in enumerate(ids):
            rhs[ia] += floc[a]
            for b, ib in enumerate(ids):
                matrix[ia][ib] += kloc[a][b]

    for tri in surface_triangles:
        nodes = _one_based_nodes(tri, 3)
        bcnr = int(getattr(tri, "bcnr", 0))
        ids = _zero_based_ids(nodes, n)
        coords = [pts[i] for i in ids]
        geom = p1_surface_triangle_geometry(coords)
        area = geom["area"]
        rcoeff = _coefficient_value(robin_coeff, bcnr, default=0.0)
        flux = _coefficient_value(boundary_flux, bcnr, default=0.0)
        boundary_area_by_number[bcnr] = boundary_area_by_number.get(bcnr, 0.0) + area
        robin_area_weight += rcoeff * area
        flux_area_weight += flux * area
        if rcoeff:
            mloc = p1_surface_triangle_mass(coords, density=rcoeff)
            for a, ia in enumerate(ids):
                for b, ib in enumerate(ids):
                    matrix[ia][ib] += mloc[a][b]
        if flux:
            floc = p1_surface_triangle_constant_load(coords, source=flux)
            for a, ia in enumerate(ids):
                rhs[ia] += floc[a]

    return {
        "matrix": matrix,
        "rhs": rhs,
        "matrix_triplets": _dense_triplets(matrix),
        "node_count": n,
        "volume_by_material": dict(sorted(volume_by_material.items())),
        "boundary_area_by_number": dict(sorted(boundary_area_by_number.items())),
        "robin_area_weight": robin_area_weight,
        "flux_area_weight": flux_area_weight,
        "policy": "dense_readable_p1_tet_with_p1_boundary_robin_trace",
    }


def _one_based_nodes(record, expected_count):
    nodes = getattr(record, "nodes", record)
    nodes = tuple(int(node) for node in nodes)
    if len(nodes) != expected_count:
        raise ValueError(f"expected {expected_count} one-based nodes, got {len(nodes)}")
    return nodes


def _zero_based_ids(nodes, point_count):
    ids = []
    for node in nodes:
        if node < 1 or node > point_count:
            raise ValueError(f"node id {node} outside 1..{point_count}")
        ids.append(node - 1)
    return ids


def _coefficient_value(value, key, default):
    if value is None:
        return float(default)
    if isinstance(value, dict):
        return float(value.get(key, default))
    return float(value)


def _dense_triplets(matrix):
    triplets = []
    for i, row in enumerate(matrix, start=1):
        for j, value in enumerate(row, start=1):
            if value != 0.0:
                triplets.append({"row": i, "col": j, "value": value})
    return triplets
