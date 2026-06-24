"""Executable electromagnetic force / energy extractors for the radia-ngsolve
A-formulation FEM path.

These are analytic- and regression-validated methods as RUNNABLE code (not just
the theory in ``differential_forms``). Each function takes ``B`` -- the magnetic
flux density CoefficientFunction, ``B = curl(gfu)`` for an HCurl GridFunction
``gfu`` -- plus the NGSolve mesh, and returns an SI quantity (force [N],
energy [J], inductance [H]).

Validated against independent references and analytics; see the
``force_validation`` MCP tool for the agreement table (sphere 0.11 %,
coil+iron force ~3 %, self-inductance 0.01 %, ...). The regression tests in
``validation/force/validate_force_xval.py`` assert these keep matching.

#25 lesson baked in: for a HIGH-permeability body do NOT carve a separate nested
"shell" material around it -- the nested-sphere interface isolates the body and
zeroes its interior B. Put the body directly in the surrounding air;
``eggshell_force`` integrates a radial weight band over that plain air, so no
extra material region is required.
"""
import math

from ngsolve import (CoefficientFunction, InnerProduct, sqrt, dx, ds, Integrate,
                     IfPos, specialcf, Conj, x, y, z)

from .scalar_fem3d import p1_surface_triangle_geometry, p1_tetrahedron_geometry

MU0 = 4.0e-7 * math.pi
EPS0 = 8.8541878128e-12
C0 = 299792458.0
ETA0 = MU0 * C0


def _float_vector(values, name):
    vec = [float(value) for value in values]
    if len(vec) not in (2, 3):
        raise ValueError(f"{name} must have length 2 or 3")
    return vec


def _complex_vector(values, name):
    vec = [complex(value) for value in values]
    if len(vec) not in (2, 3):
        raise ValueError(f"{name} must have length 2 or 3")
    return vec


def _unit_vector(values, name):
    vec = _float_vector(values, name)
    norm = math.sqrt(sum(value * value for value in vec))
    if norm <= 0.0:
        raise ValueError(f"{name} must be nonzero")
    return [value / norm for value in vec]


def _xy_point(values, name):
    point = [float(value) for value in values]
    if len(point) != 2:
        raise ValueError(f"{name} must have length 2")
    return point


def _phasor_average_factor(amplitude):
    if amplitude == "peak":
        return 0.5
    if amplitude == "rms":
        return 1.0
    raise ValueError("amplitude must be 'peak' or 'rms'")


def maxwell_stress_tensor_air(B, mu=MU0):
    """Pointwise magnetic Maxwell stress tensor in air.

    ``B`` is a 2- or 3-component flux-density vector [T].  The returned nested
    list is

        T_ij = (B_i B_j - 0.5 |B|^2 delta_ij) / mu

    in pascals.  This dependency-free helper mirrors the integrand used by the
    surface and weighted-stress FEM extractors, so examples can teach the local
    traction identity before moving to mesh integrals.
    """

    mu = float(mu)
    if mu <= 0.0:
        raise ValueError("mu must be > 0")
    b = _float_vector(B, "B")
    b2 = sum(value * value for value in b)
    dim = len(b)
    return [
        [
            (b[i] * b[j] - (0.5 * b2 if i == j else 0.0)) / mu
            for j in range(dim)
        ]
        for i in range(dim)
    ]


def maxwell_traction_air(B, normal, mu=MU0):
    """Maxwell traction vector ``T n`` in air for a unit surface normal.

    ``normal`` is normalised internally; it must have the same length as ``B``.
    For a uniform normal field this returns ``p n`` with
    ``p = B^2/(2 mu)``, which is exactly :func:`air_gap_maxwell_pressure`.
    A purely tangential field gives ``-p n`` (magnetic tension).
    """

    b = _float_vector(B, "B")
    n = _unit_vector(normal, "normal")
    if len(b) != len(n):
        raise ValueError("B and normal must have the same length")
    tensor = maxwell_stress_tensor_air(b, mu=mu)
    return [
        sum(tensor[i][j] * n[j] for j in range(len(n)))
        for i in range(len(n))
    ]


def maxwell_traction_summary(B, normal, area_m2=1.0, mu=MU0):
    """JSON-friendly Maxwell traction decomposition for one surface patch.

    The normal component is

        traction . n = (B_n^2 - |B_t|^2) / (2 mu)

    and the tangential component has magnitude ``|B_n B_t| / mu``.  ``area_m2``
    scales the traction to a force vector for simple patch/air-gap examples.
    """

    area = float(area_m2)
    if area < 0.0:
        raise ValueError("area_m2 must be >= 0")
    mu = float(mu)
    if mu <= 0.0:
        raise ValueError("mu must be > 0")
    b = _float_vector(B, "B")
    n = _unit_vector(normal, "normal")
    if len(b) != len(n):
        raise ValueError("B and normal must have the same length")
    traction = maxwell_traction_air(b, n, mu=mu)
    b_normal = sum(bi * ni for bi, ni in zip(b, n))
    b2 = sum(bi * bi for bi in b)
    b_tangent2 = max(0.0, b2 - b_normal * b_normal)
    normal_traction = sum(ti * ni for ti, ni in zip(traction, n))
    tangential_traction = [
        ti - normal_traction * ni
        for ti, ni in zip(traction, n)
    ]
    return {
        "B": b,
        "normal": n,
        "mu": mu,
        "area_m2": area,
        "B_normal_T": b_normal,
        "B_tangent_T": math.sqrt(b_tangent2),
        "traction_Pa": traction,
        "normal_traction_Pa": normal_traction,
        "normal_traction_identity_Pa": (b_normal * b_normal - b_tangent2) / (2.0 * mu),
        "tangential_traction_Pa": tangential_traction,
        "tangential_traction_magnitude_Pa": math.sqrt(
            sum(value * value for value in tangential_traction)
        ),
        "force_N": [area * value for value in traction],
    }


def electrostatic_stress_tensor(E, eps=EPS0):
    """Pointwise electrostatic Maxwell stress tensor.

    ``E`` is a 2- or 3-component electric-field vector [V/m].  The returned
    tensor is

        T_ij = eps (E_i E_j - 0.5 |E|^2 delta_ij)

    in pascals, the electric counterpart of
    :func:`maxwell_stress_tensor_air`.
    """

    eps = float(eps)
    if eps <= 0.0:
        raise ValueError("eps must be > 0")
    e = _float_vector(E, "E")
    e2 = sum(value * value for value in e)
    dim = len(e)
    return [
        [
            eps * (e[i] * e[j] - (0.5 * e2 if i == j else 0.0))
            for j in range(dim)
        ]
        for i in range(dim)
    ]


def electrostatic_traction(E, normal, eps=EPS0):
    """Electrostatic Maxwell traction vector ``T n`` for a unit normal."""

    e = _float_vector(E, "E")
    n = _unit_vector(normal, "normal")
    if len(e) != len(n):
        raise ValueError("E and normal must have the same length")
    tensor = electrostatic_stress_tensor(e, eps=eps)
    return [
        sum(tensor[i][j] * n[j] for j in range(len(n)))
        for i in range(len(n))
    ]


def electrostatic_traction_summary(E, normal, area_m2=1.0, eps=EPS0):
    """JSON-friendly electrostatic Maxwell traction decomposition."""

    area = float(area_m2)
    if area < 0.0:
        raise ValueError("area_m2 must be >= 0")
    e = _float_vector(E, "E")
    n = _unit_vector(normal, "normal")
    if len(e) != len(n):
        raise ValueError("E and normal must have the same length")
    traction = electrostatic_traction(e, n, eps=eps)
    e_normal = sum(ei * ni for ei, ni in zip(e, n))
    e2 = sum(ei * ei for ei in e)
    e_tangent2 = max(0.0, e2 - e_normal * e_normal)
    normal_traction = sum(ti * ni for ti, ni in zip(traction, n))
    tangential_traction = [
        ti - normal_traction * ni
        for ti, ni in zip(traction, n)
    ]
    return {
        "E_V_per_m": e,
        "normal": n,
        "eps": float(eps),
        "area_m2": area,
        "E_normal_V_per_m": e_normal,
        "E_tangent_V_per_m": math.sqrt(e_tangent2),
        "traction_Pa": traction,
        "normal_traction_Pa": normal_traction,
        "normal_traction_identity_Pa": float(eps) * (e_normal * e_normal - e_tangent2) / 2.0,
        "tangential_traction_Pa": tangential_traction,
        "tangential_traction_magnitude_Pa": math.sqrt(
            sum(value * value for value in tangential_traction)
        ),
        "force_N": [area * value for value in traction],
    }


def maxwell_line_segment_force_2d(p0, p1, B, mu=MU0, normal_side="right"):
    """Maxwell-stress force on one 2D contour segment, per unit depth.

    ``p0`` and ``p1`` are the segment endpoints in the planar model.  ``B`` is
    the local 2D flux density in tesla.  The returned force has units ``N/m``:
    a Maxwell traction ``T n`` [N/m^2] integrated over contour length ``ds``
    and one metre of out-of-plane depth.

    For a counter-clockwise contour around a body, the outward normal is on the
    segment's right side.  For a clockwise contour, use ``normal_side="left"``.
    """

    a = _xy_point(p0, "p0")
    b = _xy_point(p1, "p1")
    dx_seg = b[0] - a[0]
    dy_seg = b[1] - a[1]
    length = math.hypot(dx_seg, dy_seg)
    if length <= 0.0:
        raise ValueError("contour segment length must be > 0")
    tangent = [dx_seg / length, dy_seg / length]
    if normal_side == "right":
        normal = [tangent[1], -tangent[0]]
    elif normal_side == "left":
        normal = [-tangent[1], tangent[0]]
    else:
        raise ValueError("normal_side must be 'right' or 'left'")
    field = _float_vector(B, "B")
    if len(field) != 2:
        raise ValueError("B must have length 2 for a 2D contour segment")
    traction = maxwell_traction_air(field, normal, mu=mu)
    force = [length * value for value in traction]
    normal_force = sum(value * normal[i] for i, value in enumerate(force))
    tangential_force = [
        force[i] - normal_force * normal[i]
        for i in range(2)
    ]
    return {
        "p0": a,
        "p1": b,
        "midpoint": [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5],
        "length_m": length,
        "tangent": tangent,
        "normal_side": normal_side,
        "unit_normal": normal,
        "B_T": field,
        "traction_N_per_m2": traction,
        "force_per_depth_N_per_m": force,
        "normal_force_per_depth_N_per_m": normal_force,
        "tangential_force_per_depth_N_per_m": tangential_force,
    }


def _try_float_vector(values):
    try:
        vector = [float(value) for value in values]
    except (TypeError, ValueError):
        return None
    return vector if len(vector) in (2, 3) else None


def _contour_field_at(B, index, midpoint, n_segments):
    if callable(B):
        return B(midpoint)
    values = list(B)
    vector = _try_float_vector(values)
    if vector is not None:
        return vector
    if len(values) != n_segments:
        raise ValueError("B must be a constant vector, a callable, or one vector per segment")
    return values[index]


def _polygon_signed_area(points):
    if len(points) < 3:
        return 0.0
    return 0.5 * sum(
        points[i][0] * points[(i + 1) % len(points)][1]
        - points[(i + 1) % len(points)][0] * points[i][1]
        for i in range(len(points))
    )


def maxwell_contour_force_2d(vertices, B, mu=MU0, orientation="ccw", closed=True):
    """Integrate 2D Maxwell stress around a polyline contour.

    The result is a JSON-friendly summary of ``integral T n ds`` in ``N/m``
    (force per out-of-plane depth).  ``B`` can be a constant 2-vector, a
    callable taking a segment midpoint, or a list of one 2-vector per segment.

    ``orientation="ccw"`` means the body is on the left as the vertices are
    traversed, so the outward normal is the segment's right normal.  Set
    ``orientation="cw"`` for clockwise vertex order.
    """

    points = [_xy_point(point, f"vertices[{index}]") for index, point in enumerate(vertices)]
    if closed and len(points) >= 2 and points[0] == points[-1]:
        points = points[:-1]
    min_points = 3 if closed else 2
    if len(points) < min_points:
        raise ValueError(f"contour needs at least {min_points} vertices")
    if orientation == "ccw":
        normal_side = "right"
    elif orientation == "cw":
        normal_side = "left"
    else:
        raise ValueError("orientation must be 'ccw' or 'cw'")

    pairs = [
        (points[i], points[(i + 1) % len(points)])
        for i in range(len(points) if closed else len(points) - 1)
    ]
    segments = []
    for index, (p_start, p_end) in enumerate(pairs):
        midpoint = [(p_start[0] + p_end[0]) * 0.5, (p_start[1] + p_end[1]) * 0.5]
        field = _contour_field_at(B, index, midpoint, len(pairs))
        row = maxwell_line_segment_force_2d(
            p_start,
            p_end,
            field,
            mu=mu,
            normal_side=normal_side,
        )
        row["index"] = index + 1
        segments.append(row)

    total = [
        sum(row["force_per_depth_N_per_m"][axis] for row in segments)
        for axis in range(2)
    ]
    scalar_length = sum(row["length_m"] for row in segments)
    normal_abs = sum(abs(row["normal_force_per_depth_N_per_m"]) for row in segments)
    tangential_abs = sum(
        math.hypot(*row["tangential_force_per_depth_N_per_m"])
        for row in segments
    )
    return {
        "closed": bool(closed),
        "orientation": orientation,
        "normal_side": normal_side,
        "n_segments": len(segments),
        "contour_length_m": scalar_length,
        "polygon_signed_area_m2": _polygon_signed_area(points) if closed else None,
        "segments": segments,
        "total_force_per_depth_N_per_m": total,
        "total_force_magnitude_per_depth_N_per_m": math.hypot(total[0], total[1]),
        "sum_abs_normal_force_per_depth_N_per_m": normal_abs,
        "sum_abs_tangential_force_per_depth_N_per_m": tangential_abs,
    }


def time_average_maxwell_stress_tensor(E, H, eps=EPS0, mu=MU0, amplitude="peak"):
    """Time-average Maxwell stress tensor [Pa] for complex harmonic phasors.

    For peak phasors, the average dyadic factor is ``1/2``:

        <T_ij> = 1/2 Re(eps E_i E_j* + mu H_i H_j*)
                 - 1/4 (eps |E|^2 + mu |H|^2) delta_ij

    For RMS phasors, set ``amplitude="rms"`` and the dyadic factor becomes
    one.  The function returns a real nested list in pascals.  It is the local
    time-harmonic counterpart of :func:`maxwell_stress_tensor_air`.
    """

    eps = float(eps)
    mu = float(mu)
    if eps <= 0.0:
        raise ValueError("eps must be > 0")
    if mu <= 0.0:
        raise ValueError("mu must be > 0")
    e = _complex_vector(E, "E")
    h = _complex_vector(H, "H")
    if len(e) != len(h):
        raise ValueError("E and H must have the same length")
    factor = _phasor_average_factor(amplitude)
    e2 = sum(abs(value) ** 2 for value in e)
    h2 = sum(abs(value) ** 2 for value in h)
    scalar = 0.5 * (eps * e2 + mu * h2)
    dim = len(e)
    return [
        [
            factor * (
                eps * (e[i] * e[j].conjugate()).real
                + mu * (h[i] * h[j].conjugate()).real
                - (scalar if i == j else 0.0)
            )
            for j in range(dim)
        ]
        for i in range(dim)
    ]


def time_average_maxwell_traction(E, H, normal, eps=EPS0, mu=MU0, amplitude="peak"):
    """Time-average Maxwell traction vector ``<T> n`` for harmonic fields."""

    e = _complex_vector(E, "E")
    h = _complex_vector(H, "H")
    n = _unit_vector(normal, "normal")
    if len(e) != len(n) or len(h) != len(n):
        raise ValueError("E, H, and normal must have the same length")
    tensor = time_average_maxwell_stress_tensor(e, h, eps=eps, mu=mu, amplitude=amplitude)
    return [
        sum(tensor[i][j] * n[j] for j in range(len(n)))
        for i in range(len(n))
    ]


def time_average_maxwell_traction_summary(
    E,
    H,
    normal,
    area_m2=1.0,
    eps=EPS0,
    mu=MU0,
    amplitude="peak",
):
    """JSON-friendly time-harmonic Maxwell traction summary for one patch."""

    area = float(area_m2)
    if area < 0.0:
        raise ValueError("area_m2 must be >= 0")
    n = _unit_vector(normal, "normal")
    e = _complex_vector(E, "E")
    h = _complex_vector(H, "H")
    if len(e) != len(n) or len(h) != len(n):
        raise ValueError("E, H, and normal must have the same length")
    factor = _phasor_average_factor(amplitude)
    tensor = time_average_maxwell_stress_tensor(e, h, eps=eps, mu=mu, amplitude=amplitude)
    traction = [sum(tensor[i][j] * n[j] for j in range(len(n))) for i in range(len(n))]
    normal_traction = sum(ti * ni for ti, ni in zip(traction, n))
    tangential_traction = [ti - normal_traction * ni for ti, ni in zip(traction, n)]
    energy_density = factor * 0.5 * (
        float(eps) * sum(abs(value) ** 2 for value in e)
        + float(mu) * sum(abs(value) ** 2 for value in h)
    )
    return {
        "E": [[value.real, value.imag] for value in e],
        "H": [[value.real, value.imag] for value in h],
        "normal": n,
        "eps": float(eps),
        "mu": float(mu),
        "amplitude": amplitude,
        "area_m2": area,
        "stress_tensor_Pa": tensor,
        "traction_Pa": traction,
        "normal_traction_Pa": normal_traction,
        "tangential_traction_Pa": tangential_traction,
        "tangential_traction_magnitude_Pa": math.sqrt(
            sum(value * value for value in tangential_traction)
        ),
        "force_N": [area * value for value in traction],
        "average_energy_density_J_per_m3": energy_density,
    }


def surface_triangle_constant_traction_load_summary(vertices, traction_N_per_m2, pivot_m=None):
    """P1 equivalent nodal loads for a constant vector traction on a triangle.

    This is the tiny readable FEM/BEM boundary-load block: for one flat P1
    surface triangle with area ``A`` and constant traction ``t`` [N/m2],

        F_e = A t,        f_i = F_e / 3  for i=1,2,3.

    The equal nodal distribution preserves both the integrated force and the
    moment of the constant traction, because the triangle centroid is the mean
    of its vertices.  That makes this helper a clean bridge between MATLAB
    teaching scripts, `.vol` boundary triangles, and the Maxwell-traction
    special cases below.
    """

    verts = [_float_vector(vertex, f"vertices[{index}]") for index, vertex in enumerate(vertices)]
    if len(verts) != 3:
        raise ValueError("a surface triangle needs exactly three vertices")
    if any(len(vertex) != 3 for vertex in verts):
        raise ValueError("surface triangle vertices must be 3D points")
    traction = _float_vector(traction_N_per_m2, "traction_N_per_m2")
    if len(traction) != 3:
        raise ValueError("traction_N_per_m2 must have length 3")

    geom = p1_surface_triangle_geometry(verts)
    area = geom["area"]
    centroid = [
        sum(vertex[axis] for vertex in verts) / 3.0
        for axis in range(3)
    ]
    integrated_force = [area * value for value in traction]
    nodal_loads = [[value / 3.0 for value in integrated_force] for _ in range(3)]
    patch_resultant = force_moment_resultant_summary(
        [centroid],
        [integrated_force],
        pivot_m=pivot_m,
    )
    nodal_resultant = force_moment_resultant_summary(
        verts,
        nodal_loads,
        pivot_m=pivot_m,
    )
    moment_errors = [
        abs(nodal_resultant["total_moment"][axis] - patch_resultant["total_moment"][axis])
        for axis in range(3)
    ]
    force_errors = [
        abs(nodal_resultant["total_force"][axis] - integrated_force[axis])
        for axis in range(3)
    ]
    return {
        "area": area,
        "centroid_m": centroid,
        "unit_normal": list(geom["unit_normal"]),
        "area_vector": list(geom["area_vector"]),
        "traction_N_per_m2": traction,
        "integrated_force_N": integrated_force,
        "nodal_force_loads_N": nodal_loads,
        "p1_shape_function_integral_m2": area / 3.0,
        "pivot_m": patch_resultant["pivot_m"],
        "patch_resultant": patch_resultant,
        "nodal_resultant": nodal_resultant,
        "force_preservation_abs_errors_N": force_errors,
        "moment_preservation_abs_errors_Nm": moment_errors,
    }


def surface_triangle_maxwell_traction_summary(vertices, B, mu=MU0):
    """P1 surface-triangle Maxwell traction and equivalent nodal force load.

    ``vertices`` is one 3D boundary triangle, using its stored orientation.
    The returned force is ``area * (T n)`` for the triangle unit normal.  The
    constant traction is also distributed as a P1 equivalent nodal load:
    each of the three nodes receives one third of the integrated force.

    This is intentionally small and explicit so the same block can be mirrored
    in first-order teaching scripts before using a full FEM surface
    integral over a curved or high-order mesh.
    """

    geom = p1_surface_triangle_geometry(vertices)
    traction = maxwell_traction_air(B, geom["unit_normal"], mu=mu)
    force = [geom["area"] * value for value in traction]
    return {
        "area": geom["area"],
        "unit_normal": list(geom["unit_normal"]),
        "area_vector": list(geom["area_vector"]),
        "B": _float_vector(B, "B"),
        "mu": float(mu),
        "traction_Pa": traction,
        "integrated_force_N": force,
        "nodal_force_loads_N": [[value / 3.0 for value in force] for _ in range(3)],
    }


def tetrahedron_lorentz_force_summary(vertices, current_density, B):
    """P1 tetrahedron equivalent nodal load for constant Lorentz force density.

    ``current_density`` is ``J`` [A/m2] and ``B`` is flux density [T].  For a
    constant field over a first-order tetrahedron,

        f = J x B,    F_e = volume * f

    and the consistent P1 body-force load gives one quarter of ``F_e`` to each
    vertex.  This is the volume counterpart of
    :func:`surface_triangle_maxwell_traction_summary`.
    """

    geom = p1_tetrahedron_geometry(vertices)
    volume = geom["volume"]
    j = _float_vector(current_density, "current_density")
    b = _float_vector(B, "B")
    if len(j) != 3 or len(b) != 3:
        raise ValueError("current_density and B must have length 3")
    force_density = [
        j[1] * b[2] - j[2] * b[1],
        j[2] * b[0] - j[0] * b[2],
        j[0] * b[1] - j[1] * b[0],
    ]
    force = [volume * value for value in force_density]
    return {
        "volume_m3": volume,
        "current_density_A_per_m2": j,
        "B_T": b,
        "force_density_N_per_m3": force_density,
        "integrated_force_N": force,
        "nodal_force_loads_N": [[value / 4.0 for value in force] for _ in range(4)],
        "p1_shape_function_integral_m3": volume / 4.0,
    }


def planar_lorentz_force_summary(Jz_A_per_m2, B_xy_T, area_m2=1.0):
    """2D planar Lorentz force from uniform out-of-plane current density.

    For magnetostatic ``A_z`` formulations and planar block integrals,
    ``J = (0, 0, Jz)`` and ``B = (Bx, By, 0)``.  The body-force density is

        J x B = (-Jz * By, Jz * Bx, 0).

    Integrating over cross-section area gives force per out-of-plane depth
    [N/m].  This dependency-free helper mirrors the local integrand of
    :func:`lorentz_force_2d`.
    """

    jz = float(Jz_A_per_m2)
    b = _float_vector(B_xy_T, "B_xy_T")
    if len(b) != 2:
        raise ValueError("B_xy_T must have length 2")
    area = float(area_m2)
    if area < 0.0:
        raise ValueError("area_m2 must be >= 0")
    force_density = [-jz * b[1], jz * b[0]]
    force = [area * value for value in force_density]
    return {
        "Jz_A_per_m2": jz,
        "B_xy_T": b,
        "area_m2": area,
        "current_A": jz * area,
        "force_density_N_per_m3": force_density,
        "force_per_depth_N_per_m": force,
        "force_magnitude_per_depth_N_per_m": math.hypot(force[0], force[1]),
    }


def _xy_separation_vector(separation_xy_m):
    try:
        values = [float(value) for value in separation_xy_m]
    except TypeError:
        distance = float(separation_xy_m)
        if distance <= 0.0:
            raise ValueError("scalar separation must be > 0")
        values = [distance, 0.0]
    if len(values) != 2:
        raise ValueError("separation_xy_m must be a scalar distance or a 2D vector")
    distance = math.hypot(values[0], values[1])
    if distance <= 0.0:
        raise ValueError("separation_xy_m must be nonzero")
    return values, distance


def parallel_wire_lorentz_force_summary(current1_A, current2_A, separation_xy_m):
    """Ampere two-wire force as a signed 2D Lorentz-force summary.

    Wire 1 is at the origin, wire 2 is displaced by ``separation_xy_m``, and
    positive current flows along ``+z``.  The magnetic field from wire 1 at wire
    2 is evaluated with the right-hand rule, then the force on wire 2 is
    ``I2 zhat x B1``.  Like currents attract, so for a scalar positive
    separation the force on wire 2 points in ``-x``.
    """

    i1 = float(current1_A)
    i2 = float(current2_A)
    separation, distance = _xy_separation_vector(separation_xy_m)
    unit = [separation[0] / distance, separation[1] / distance]
    tangent = [-unit[1], unit[0]]
    b_mag = MU0 * i1 / (2.0 * math.pi * distance)
    field_at_wire2 = [b_mag * tangent[0], b_mag * tangent[1]]
    signed_ampere_force = MU0 * i1 * i2 / (2.0 * math.pi * distance)
    force_on_wire2 = [-signed_ampere_force * unit[0], -signed_ampere_force * unit[1]]
    force_on_wire1 = [signed_ampere_force * unit[0], signed_ampere_force * unit[1]]
    if signed_ampere_force > 0.0:
        interaction = "attraction"
    elif signed_ampere_force < 0.0:
        interaction = "repulsion"
    else:
        interaction = "zero"
    return {
        "current1_A": i1,
        "current2_A": i2,
        "separation_xy_m": separation,
        "separation_m": distance,
        "unit_from_wire1_to_wire2": unit,
        "right_hand_tangent_at_wire2": tangent,
        "field_from_wire1_at_wire2_T": field_at_wire2,
        "signed_ampere_force_per_length_N_per_m": signed_ampere_force,
        "force_magnitude_per_length_N_per_m": abs(signed_ampere_force),
        "force_on_wire1_N_per_m": force_on_wire1,
        "force_on_wire2_N_per_m": force_on_wire2,
        "interaction": interaction,
    }


def force_moment_resultant_summary(points_m, forces, pivot_m=None):
    """Resultant force and torque from discrete force rows.

    ``points_m`` and ``forces`` are matching 2D or 3D vectors.  The returned
    moment is

        M_p = sum_i (r_i - p) x F_i

    about ``pivot_m``.  For 2D inputs the cross product is the scalar out-of-
    plane moment ``Mz``; for 3D inputs it is a 3-vector.  This dependency-free
    helper is the common final post-processing step for Maxwell-stress patches,
    Lorentz body-force elements, nodal loads, and mesh boundary pressure rows.
    """

    points = [_float_vector(point, f"points_m[{index}]") for index, point in enumerate(points_m)]
    force_rows = [_float_vector(force, f"forces[{index}]") for index, force in enumerate(forces)]
    if not points:
        raise ValueError("at least one force row is required")
    if len(points) != len(force_rows):
        raise ValueError("points_m and forces must have the same length")
    dim = len(points[0])
    if dim not in (2, 3):
        raise ValueError("force rows must be 2D or 3D")
    if any(len(point) != dim for point in points):
        raise ValueError("all points_m rows must have the same dimension")
    if any(len(force) != dim for force in force_rows):
        raise ValueError("all forces rows must match the point dimension")
    pivot = [0.0] * dim if pivot_m is None else _float_vector(pivot_m, "pivot_m")
    if len(pivot) != dim:
        raise ValueError("pivot_m must match the point dimension")

    total_force = [sum(force[axis] for force in force_rows) for axis in range(dim)]
    rows = []
    if dim == 2:
        total_moment = 0.0
        for index, (point, force) in enumerate(zip(points, force_rows), start=1):
            lever = [point[0] - pivot[0], point[1] - pivot[1]]
            moment = lever[0] * force[1] - lever[1] * force[0]
            total_moment += moment
            rows.append({
                "index": index,
                "point_m": point,
                "force": force,
                "lever_arm_m": lever,
                "moment_z": moment,
            })
        moment_magnitude = abs(total_moment)
    else:
        total_moment = [0.0, 0.0, 0.0]
        for index, (point, force) in enumerate(zip(points, force_rows), start=1):
            lever = [point[axis] - pivot[axis] for axis in range(3)]
            moment = [
                lever[1] * force[2] - lever[2] * force[1],
                lever[2] * force[0] - lever[0] * force[2],
                lever[0] * force[1] - lever[1] * force[0],
            ]
            total_moment = [total_moment[axis] + moment[axis] for axis in range(3)]
            rows.append({
                "index": index,
                "point_m": point,
                "force": force,
                "lever_arm_m": lever,
                "moment": moment,
            })
        moment_magnitude = math.sqrt(sum(value * value for value in total_moment))

    return {
        "dimension": dim,
        "n_rows": len(rows),
        "pivot_m": pivot,
        "rows": rows,
        "total_force": total_force,
        "total_force_magnitude": math.sqrt(sum(value * value for value in total_force)),
        "total_moment": total_moment,
        "total_moment_magnitude": moment_magnitude,
    }


def air_gap_maxwell_pressure(B_T, mu=MU0):
    """Magnetic pressure [Pa] for a normal flux density in an air gap.

    For a locally uniform normal field at an iron/air interface, the Maxwell
    stress gives

        p = B^2 / (2 mu)

    with ``mu=mu0`` for air.  The same value is the magnetic energy density in
    the gap.  This tiny helper is intentionally dependency-free so magnetic
    circuit examples can turn a solved ``B_T`` directly into a holding-force
    estimate before running a full weighted-stress FEM extraction.
    """

    B = float(B_T)
    mu = float(mu)
    if mu <= 0.0:
        raise ValueError("mu must be > 0")
    return B * B / (2.0 * mu)


def air_gap_holding_force(B_T, area_m2, faces=1, mu=MU0):
    """Uniform-gap holding force [N] from flux density and active pole area.

    ``faces`` is the number of active, equal pole faces/gaps contributing the
    same pressure.  Use ``faces=2`` for a symmetric two-pole yoke with two equal
    gaps; keep ``faces=1`` for a single plunger or one pole face.
    """

    area = float(area_m2)
    faces = int(faces)
    if area < 0.0:
        raise ValueError("area_m2 must be >= 0")
    if faces < 1:
        raise ValueError("faces must be >= 1")
    return air_gap_maxwell_pressure(B_T, mu=mu) * area * faces


def air_gap_force_summary(B_T, area_m2, faces=1, mu=MU0):
    """Readable JSON-friendly air-gap force summary."""

    pressure = air_gap_maxwell_pressure(B_T, mu=mu)
    area = float(area_m2)
    faces = int(faces)
    force = air_gap_holding_force(B_T, area, faces=faces, mu=mu)
    return {
        "B_T": float(B_T),
        "mu": float(mu),
        "area_m2": area,
        "faces": faces,
        "pressure_Pa": pressure,
        "energy_density_J_per_m3": pressure,
        "force_N": force,
        "force_per_area_N_per_m2": force / (area * faces) if area > 0.0 else math.inf,
    }


def air_gap_shear_stress(B_radial_T, B_tangential_T, mu=MU0):
    """Tangential Maxwell shear stress [Pa] in a cylindrical air gap.

    For radial and tangential flux-density components ``Br`` and ``Bt`` on a
    cylindrical integration surface, the tangential traction is

        tau = Br Bt / mu

    with sign set by ``Bt``.  This is the local stress used by air-gap motor
    torque estimates and by FE Maxwell-stress post-processing.
    """

    mu = float(mu)
    if mu <= 0.0:
        raise ValueError("mu must be > 0")
    return float(B_radial_T) * float(B_tangential_T) / mu


def air_gap_shear_torque(
    B_radial_T,
    B_tangential_T,
    radius_m,
    axial_length_m=1.0,
    angle_rad=2.0 * math.pi,
    mu=MU0,
):
    """Torque [N.m] from uniform air-gap Maxwell shear stress.

    For a cylindrical surface patch, ``area = radius * angle * axial_length``
    and ``torque = radius * tau * area``.  Use ``angle_rad=2*pi`` for a full
    machine, or a sector angle for a symmetry-sector result before multiplying
    by the sector count.
    """

    radius = float(radius_m)
    length = float(axial_length_m)
    angle = float(angle_rad)
    if radius < 0.0:
        raise ValueError("radius_m must be >= 0")
    if length < 0.0:
        raise ValueError("axial_length_m must be >= 0")
    if angle < 0.0:
        raise ValueError("angle_rad must be >= 0")
    shear = air_gap_shear_stress(B_radial_T, B_tangential_T, mu=mu)
    return shear * radius * radius * angle * length


def air_gap_shear_torque_summary(
    B_radial_T,
    B_tangential_T,
    radius_m,
    axial_length_m=1.0,
    angle_rad=2.0 * math.pi,
    mu=MU0,
):
    """JSON-friendly air-gap shear-stress torque summary."""

    radius = float(radius_m)
    length = float(axial_length_m)
    angle = float(angle_rad)
    if radius < 0.0:
        raise ValueError("radius_m must be >= 0")
    if length < 0.0:
        raise ValueError("axial_length_m must be >= 0")
    if angle < 0.0:
        raise ValueError("angle_rad must be >= 0")
    shear = air_gap_shear_stress(B_radial_T, B_tangential_T, mu=mu)
    area = radius * angle * length
    force = shear * area
    torque = force * radius
    return {
        "B_radial_T": float(B_radial_T),
        "B_tangential_T": float(B_tangential_T),
        "mu": float(mu),
        "radius_m": radius,
        "axial_length_m": length,
        "angle_rad": angle,
        "surface_area_m2": area,
        "shear_stress_Pa": shear,
        "tangential_force_N": force,
        "torque_Nm": torque,
        "torque_per_axial_length_N": torque / length if length > 0.0 else math.inf,
    }


def air_gap_shear_torque_from_angle_samples(
    angles_rad,
    B_radial_T,
    B_tangential_T,
    radius_m,
    axial_length_m=1.0,
    periodic=True,
    period_rad=2.0 * math.pi,
    mu=MU0,
):
    """Integrate sampled air-gap Maxwell shear stress into machine torque.

    Electric-machine solvers often export air-gap samples around a cylindrical
    contour.  For each angle sample,

        tau(theta) = Br(theta) Bt(theta) / mu

    and the torque is

        T = r^2 L integral tau(theta) dtheta.

    If ``periodic=True`` the last segment wraps from the last sample to the
    first sample plus ``period_rad``; omit a duplicate endpoint.  If
    ``periodic=False`` only the provided angular span is integrated.  The
    returned rows are segment contributions, useful for teaching and for
    checking sector-model scaling before comparing whole-machine torque.
    """

    angles = [float(value) for value in angles_rad]
    br = [float(value) for value in B_radial_T]
    bt = [float(value) for value in B_tangential_T]
    if len(angles) != len(br) or len(angles) != len(bt):
        raise ValueError("angles_rad, B_radial_T, and B_tangential_T must have the same length")
    if len(angles) < 2:
        raise ValueError("at least two angle samples are required")
    if any(angles[i + 1] <= angles[i] for i in range(len(angles) - 1)):
        raise ValueError("angles_rad must be strictly increasing")
    radius = float(radius_m)
    length = float(axial_length_m)
    if radius < 0.0:
        raise ValueError("radius_m must be >= 0")
    if length < 0.0:
        raise ValueError("axial_length_m must be >= 0")
    mu = float(mu)
    if mu <= 0.0:
        raise ValueError("mu must be > 0")
    period = float(period_rad)
    if periodic and period <= 0.0:
        raise ValueError("period_rad must be > 0")

    shear = [
        air_gap_shear_stress(bri, bti, mu=mu)
        for bri, bti in zip(br, bt)
    ]
    n = len(angles)
    segment_count = n if periodic else n - 1
    rows = []
    integral_shear = 0.0
    for i in range(segment_count):
        j = (i + 1) % n
        theta0 = angles[i]
        theta1 = angles[j]
        if periodic and j == 0:
            theta1 += period
        dtheta = theta1 - theta0
        if dtheta <= 0.0:
            raise ValueError("angle segment width must be > 0")
        shear_avg = 0.5 * (shear[i] + shear[j])
        tangential_force = shear_avg * radius * length * dtheta
        torque = tangential_force * radius
        integral_shear += shear_avg * dtheta
        rows.append({
            "segment_index": i + 1,
            "angle_start_rad": theta0,
            "angle_end_rad": theta1,
            "angle_width_rad": dtheta,
            "B_radial_start_T": br[i],
            "B_radial_end_T": br[j],
            "B_tangential_start_T": bt[i],
            "B_tangential_end_T": bt[j],
            "shear_start_Pa": shear[i],
            "shear_end_Pa": shear[j],
            "shear_average_Pa": shear_avg,
            "tangential_force_N": tangential_force,
            "torque_Nm": torque,
        })

    torque_total = radius * radius * length * integral_shear
    force_total = radius * length * integral_shear
    integrated_angle = sum(row["angle_width_rad"] for row in rows)
    return {
        "n_samples": n,
        "n_segments": len(rows),
        "periodic": bool(periodic),
        "period_rad": period,
        "radius_m": radius,
        "axial_length_m": length,
        "mu": mu,
        "integrated_angle_rad": integrated_angle,
        "integral_shear_dtheta_Pa_rad": integral_shear,
        "average_shear_stress_Pa": integral_shear / integrated_angle if integrated_angle > 0.0 else math.nan,
        "tangential_force_N": force_total,
        "torque_Nm": torque_total,
        "torque_per_axial_length_N": torque_total / length if length > 0.0 else math.inf,
        "rows": rows,
    }


def coenergy_torque_from_angle_samples(
    angles_rad,
    coenergy_J,
    periodic=False,
    period_rad=2.0 * math.pi,
):
    """Differentiate a coenergy-vs-angle table into torque samples.

    For fixed currents, the virtual-work torque is

        T(theta) = dW'(theta) / dtheta

    where ``W'`` is magnetic coenergy.  This helper turns an angle sweep into a
    readable finite-difference table, matching the way motor solvers often
    report torque from a sequence of rotor-position solves.

    If ``periodic=True``, the first and last samples are wrapped across
    ``period_rad``; provide one sample per angle and omit the duplicate endpoint.
    """

    angles = [float(value) for value in angles_rad]
    values = [float(value) for value in coenergy_J]
    if len(angles) != len(values):
        raise ValueError("angles_rad and coenergy_J must have the same length")
    if len(angles) < 3:
        raise ValueError("at least three samples are required")
    if any(angles[i + 1] <= angles[i] for i in range(len(angles) - 1)):
        raise ValueError("angles_rad must be strictly increasing")
    period = float(period_rad)
    if periodic and period <= 0.0:
        raise ValueError("period_rad must be > 0")

    rows = []
    n = len(angles)
    for i, (angle, value) in enumerate(zip(angles, values)):
        if periodic:
            im = (i - 1) % n
            ip = (i + 1) % n
            angle_minus = angles[im]
            angle_plus = angles[ip]
            if im > i:
                angle_minus -= period
            if ip < i:
                angle_plus += period
            stencil = "central_periodic"
        elif i == 0:
            im, ip = 0, 1
            angle_minus = angles[im]
            angle_plus = angles[ip]
            stencil = "forward"
        elif i == n - 1:
            im, ip = n - 2, n - 1
            angle_minus = angles[im]
            angle_plus = angles[ip]
            stencil = "backward"
        else:
            im, ip = i - 1, i + 1
            angle_minus = angles[im]
            angle_plus = angles[ip]
            stencil = "central"

        denom = angle_plus - angle_minus
        if denom <= 0.0:
            raise ValueError("finite-difference angle denominator must be > 0")
        torque = (values[ip] - values[im]) / denom
        rows.append({
            "index": i + 1,
            "angle_rad": angle,
            "coenergy_J": value,
            "torque_Nm": torque,
            "stencil": stencil,
            "angle_minus_rad": angle_minus,
            "angle_plus_rad": angle_plus,
        })
    return tuple(rows)


def coenergy_torque_summary(
    angles_rad,
    coenergy_J,
    periodic=False,
    period_rad=2.0 * math.pi,
):
    """JSON-friendly summary for coenergy-derived torque samples."""

    rows = list(coenergy_torque_from_angle_samples(
        angles_rad,
        coenergy_J,
        periodic=periodic,
        period_rad=period_rad,
    ))
    torques = [row["torque_Nm"] for row in rows]
    return {
        "n_samples": len(rows),
        "periodic": bool(periodic),
        "period_rad": float(period_rad),
        "torque_min_Nm": min(torques),
        "torque_max_Nm": max(torques),
        "torque_peak_abs_Nm": max(abs(value) for value in torques),
        "torque_mean_Nm": sum(torques) / len(torques),
        "rows": rows,
    }


def _virtual_work_energy_sign(energy_kind):
    kind = str(energy_kind).lower().replace("-", "_").replace(" ", "_")
    if kind in ("coenergy", "magnetic_coenergy", "wco", "w_prime", "constant_current"):
        return "coenergy", 1.0, "F = dW_co/dx at fixed current"
    if kind in ("stored_energy", "field_energy", "energy", "magnetic_energy", "constant_flux"):
        return "stored_energy", -1.0, "F = -dW/dx at fixed flux/source-free displacement"
    raise ValueError(
        "energy_kind must be 'coenergy'/'constant_current' or "
        "'stored_energy'/'field_energy'"
    )


def virtual_work_force_from_displacement_samples(
    positions_m,
    energy_J,
    energy_kind="coenergy",
):
    """Differentiate an energy/coenergy-vs-displacement table into force samples.

    This is the straight-line counterpart of
    :func:`coenergy_torque_from_angle_samples`.  It makes the sign convention
    explicit for validation sweeps and solver cross-checks:

    * fixed current / magnetic coenergy: ``F = dW_co/dx``
    * stored field energy at fixed flux: ``F = -dW/dx``

    The derivative is reported in newtons because ``J/m = N``.  End points use
    one-sided differences; interior rows use central differences.  Use matched
    meshes or a deliberately stable remeshing recipe when the samples come from
    separate FEM solves.
    """

    positions = [float(value) for value in positions_m]
    values = [float(value) for value in energy_J]
    if len(positions) != len(values):
        raise ValueError("positions_m and energy_J must have the same length")
    if len(positions) < 3:
        raise ValueError("at least three samples are required")
    if any(positions[i + 1] <= positions[i] for i in range(len(positions) - 1)):
        raise ValueError("positions_m must be strictly increasing")
    normalized_kind, sign, identity = _virtual_work_energy_sign(energy_kind)

    rows = []
    n = len(positions)
    for i, (position, value) in enumerate(zip(positions, values)):
        if i == 0:
            im, ip = 0, 1
            stencil = "forward"
        elif i == n - 1:
            im, ip = n - 2, n - 1
            stencil = "backward"
        else:
            im, ip = i - 1, i + 1
            stencil = "central"

        denom = positions[ip] - positions[im]
        if denom <= 0.0:
            raise ValueError("finite-difference displacement denominator must be > 0")
        derivative = (values[ip] - values[im]) / denom
        force = sign * derivative
        rows.append({
            "index": i + 1,
            "position_m": position,
            "energy_J": value,
            "energy_kind": normalized_kind,
            "virtual_work_identity": identity,
            "denergy_dx_N": derivative,
            "force_N": force,
            "stencil": stencil,
            "position_minus_m": positions[im],
            "position_plus_m": positions[ip],
        })
    return tuple(rows)


def virtual_work_force_summary(
    positions_m,
    energy_J,
    energy_kind="coenergy",
):
    """JSON-friendly summary for virtual-work force samples."""

    normalized_kind, sign, identity = _virtual_work_energy_sign(energy_kind)
    rows = list(virtual_work_force_from_displacement_samples(
        positions_m,
        energy_J,
        energy_kind=normalized_kind,
    ))
    forces = [row["force_N"] for row in rows]
    return {
        "n_samples": len(rows),
        "energy_kind": normalized_kind,
        "energy_to_force_sign": sign,
        "virtual_work_identity": identity,
        "force_min_N": min(forces),
        "force_max_N": max(forces),
        "force_peak_abs_N": max(abs(value) for value in forces),
        "force_mean_N": sum(forces) / len(forces),
        "rows": rows,
    }


def virtual_work_symmetric_pair_force_summary(
    displacement_m,
    energy_minus_J,
    energy_plus_J,
    energy_kind="coenergy",
    center_position_m=0.0,
    energy_center_J=None,
):
    """Virtual-work force from a matched ``x0 +/- h`` energy pair.

    This is the compact two-solve central-difference gate often used for
    magnetostatic force sweeps when the displaced geometries are deliberately
    paired:

        F = +/- (E(x0 + h) - E(x0 - h)) / (2 h)

    The sign is positive for fixed-current coenergy and negative for fixed-flux
    stored energy, matching :func:`virtual_work_force_summary`.  If
    ``energy_center_J`` is supplied, the even residual
    ``0.5 * (E_+ + E_-) - E_0`` is reported as a curvature/noise indicator.
    """

    h = float(displacement_m)
    if h <= 0.0:
        raise ValueError("displacement_m must be > 0")
    normalized_kind, sign, identity = _virtual_work_energy_sign(energy_kind)
    e_minus = float(energy_minus_J)
    e_plus = float(energy_plus_J)
    derivative = (e_plus - e_minus) / (2.0 * h)
    force = sign * derivative
    center = float(center_position_m)
    even_residual = None
    if energy_center_J is not None:
        even_residual = 0.5 * (e_plus + e_minus) - float(energy_center_J)

    return {
        "stencil": "symmetric_central_pair",
        "center_position_m": center,
        "displacement_m": h,
        "position_minus_m": center - h,
        "position_plus_m": center + h,
        "energy_minus_J": e_minus,
        "energy_plus_J": e_plus,
        "energy_center_J": None if energy_center_J is None else float(energy_center_J),
        "energy_kind": normalized_kind,
        "energy_to_force_sign": sign,
        "virtual_work_identity": identity,
        "denergy_dx_N": derivative,
        "force_N": force,
        "even_energy_residual_J": even_residual,
        "even_energy_residual_abs_J": None if even_residual is None else abs(even_residual),
    }


def _validate_optical_coefficients(absorptance, reflectance):
    absorptance = float(absorptance)
    reflectance = float(reflectance)
    if absorptance < 0.0:
        raise ValueError("absorptance must be >= 0")
    if reflectance < 0.0:
        raise ValueError("reflectance must be >= 0")
    if absorptance + reflectance > 1.0 + 1e-15:
        raise ValueError("absorptance + reflectance must be <= 1")
    return absorptance, reflectance


def _validate_reflectance_transmittance(reflectance, transmittance):
    reflectance = float(reflectance)
    transmittance = float(transmittance)
    if reflectance < 0.0:
        raise ValueError("reflectance must be >= 0")
    if transmittance < 0.0:
        raise ValueError("transmittance must be >= 0")
    if reflectance + transmittance > 1.0 + 1e-15:
        raise ValueError("reflectance + transmittance must be <= 1")
    absorptance = max(0.0, 1.0 - reflectance - transmittance)
    return reflectance, transmittance, absorptance


def plane_wave_intensity_from_electric_field(
    electric_field_V_per_m,
    impedance_ohm=ETA0,
    amplitude="rms",
):
    """Time-average plane-wave intensity [W/m2] from electric-field amplitude.

    ``amplitude="rms"`` uses ``I = E_rms^2 / eta``.  ``amplitude="peak"`` uses
    ``I = E_peak^2 / (2 eta)``.  This small helper makes RF power-flow examples
    use the same RMS/peak convention before converting Poynting flux to force.
    """

    field = float(electric_field_V_per_m)
    impedance = float(impedance_ohm)
    if impedance <= 0.0:
        raise ValueError("impedance_ohm must be > 0")
    if field < 0.0:
        raise ValueError("electric_field_V_per_m must be >= 0")
    if amplitude == "rms":
        return field * field / impedance
    if amplitude == "peak":
        return field * field / (2.0 * impedance)
    raise ValueError("amplitude must be 'rms' or 'peak'")


def radiation_pressure_from_intensity(
    intensity_W_per_m2,
    absorptance=1.0,
    reflectance=0.0,
    speed=C0,
):
    """Normal-incidence radiation pressure [Pa] from time-average intensity.

    For a normally incident wave with intensity ``I``, the normal momentum flux
    is ``I / c``.  A perfectly absorbing surface takes one unit of photon
    momentum, while a perfect mirror reverses it:

        p = (absorptance + 2 reflectance) I / c

    The remaining fraction is transmitted and contributes no force on the
    surface.  This is the RF/time-harmonic counterpart of the static Maxwell
    stress helpers above.
    """

    intensity = float(intensity_W_per_m2)
    speed = float(speed)
    if intensity < 0.0:
        raise ValueError("intensity_W_per_m2 must be >= 0")
    if speed <= 0.0:
        raise ValueError("speed must be > 0")
    absorptance, reflectance = _validate_optical_coefficients(absorptance, reflectance)
    return (absorptance + 2.0 * reflectance) * intensity / speed


def radiation_force_from_power(
    power_W,
    absorptance=1.0,
    reflectance=0.0,
    speed=C0,
):
    """Normal force [N] from normally incident time-average RF/optical power."""

    power = float(power_W)
    speed = float(speed)
    if power < 0.0:
        raise ValueError("power_W must be >= 0")
    if speed <= 0.0:
        raise ValueError("speed must be > 0")
    absorptance, reflectance = _validate_optical_coefficients(absorptance, reflectance)
    return (absorptance + 2.0 * reflectance) * power / speed


def radiation_force_from_normal_scattering(
    power_incident_W,
    reflectance,
    transmittance=0.0,
    speed=C0,
):
    """Normal force [N] from one-sided normal-incidence scattering powers.

    Use ``reflectance=|S11|^2`` and ``transmittance=|S21|^2`` when the input and
    output ports have the same power normalization.  The momentum balance on the
    scatterer is

        F = (1 + R - T) P_inc / c

    which is identical to ``(A + 2R) P_inc/c`` with
    ``A = 1 - R - T``.  The limits are: absorber ``P/c``, mirror ``2P/c``,
    and transparent through-line ``0``.
    """

    power = float(power_incident_W)
    speed = float(speed)
    if power < 0.0:
        raise ValueError("power_incident_W must be >= 0")
    if speed <= 0.0:
        raise ValueError("speed must be > 0")
    reflectance, transmittance, _absorptance = _validate_reflectance_transmittance(
        reflectance,
        transmittance,
    )
    return (1.0 + reflectance - transmittance) * power / speed


def radiation_scattering_force_summary(
    power_incident_W,
    reflectance,
    transmittance=0.0,
    speed=C0,
):
    """JSON-friendly normal-incidence scattering radiation-force summary."""

    power = float(power_incident_W)
    speed = float(speed)
    if power < 0.0:
        raise ValueError("power_incident_W must be >= 0")
    if speed <= 0.0:
        raise ValueError("speed must be > 0")
    reflectance, transmittance, absorptance = _validate_reflectance_transmittance(
        reflectance,
        transmittance,
    )
    factor = 1.0 + reflectance - transmittance
    force = factor * power / speed
    return {
        "power_incident_W": power,
        "reflectance": reflectance,
        "transmittance": transmittance,
        "absorptance": absorptance,
        "power_reflected_W": reflectance * power,
        "power_transmitted_W": transmittance * power,
        "power_absorbed_W": absorptance * power,
        "speed_m_per_s": speed,
        "momentum_transfer_factor": factor,
        "absorber_reflector_equivalent_factor": absorptance + 2.0 * reflectance,
        "force_N": force,
        "force_from_absorptance_reflectance_N": radiation_force_from_power(
            power,
            absorptance=absorptance,
            reflectance=reflectance,
            speed=speed,
        ),
    }


def two_port_scattering_momentum_force_summary(
    power_incident_W,
    incident_direction,
    reflectance,
    transmittance=0.0,
    transmitted_direction=None,
    speed=C0,
):
    """Vector force [N] from a one-sided two-port scattering momentum balance.

    ``incident_direction`` and ``transmitted_direction`` are propagation
    directions of the incident and transmitted power flows.  Reflection is
    assumed to leave back toward the source along ``-incident_direction``.
    The force on the scatterer is the incoming field momentum minus outgoing
    field momentum:

        F = P/c * ((1 + R) k_inc - T k_out)

    For a straight through-line ``k_out = k_inc`` this reduces to
    :func:`radiation_force_from_normal_scattering`.  For a lossless 90-degree
    bend with no reflection it gives ``P/c * (k_inc - k_out)``, the readable
    vector gate for RF waveguide momentum bookkeeping.
    """

    power = float(power_incident_W)
    speed = float(speed)
    if power < 0.0:
        raise ValueError("power_incident_W must be >= 0")
    if speed <= 0.0:
        raise ValueError("speed must be > 0")
    inc = _unit_vector(incident_direction, "incident_direction")
    if transmitted_direction is None:
        out = list(inc)
    else:
        out = _unit_vector(transmitted_direction, "transmitted_direction")
        if len(out) != len(inc):
            raise ValueError("incident_direction and transmitted_direction must have the same length")
    reflectance, transmittance, absorptance = _validate_reflectance_transmittance(
        reflectance,
        transmittance,
    )

    momentum_scale = power / speed
    incident_momentum = [momentum_scale * value for value in inc]
    reflected_momentum = [-reflectance * momentum_scale * value for value in inc]
    transmitted_momentum = [transmittance * momentum_scale * value for value in out]
    force = [
        incident_momentum[i] - reflected_momentum[i] - transmitted_momentum[i]
        for i in range(len(inc))
    ]
    axial_force = sum(fi * ki for fi, ki in zip(force, inc))
    return {
        "power_incident_W": power,
        "reflectance": reflectance,
        "transmittance": transmittance,
        "absorptance": absorptance,
        "power_reflected_W": reflectance * power,
        "power_transmitted_W": transmittance * power,
        "power_absorbed_W": absorptance * power,
        "speed_m_per_s": speed,
        "incident_direction": inc,
        "transmitted_direction": out,
        "incident_momentum_flow_N": incident_momentum,
        "reflected_momentum_flow_N": reflected_momentum,
        "transmitted_momentum_flow_N": transmitted_momentum,
        "force_N": force,
        "force_magnitude_N": math.sqrt(sum(value * value for value in force)),
        "axial_force_along_incident_direction_N": axial_force,
        "straight_through_equivalent_force_N": radiation_force_from_normal_scattering(
            power,
            reflectance,
            transmittance,
            speed=speed,
        ),
    }


def radiation_pressure_summary(
    intensity_W_per_m2,
    area_m2=1.0,
    absorptance=1.0,
    reflectance=0.0,
    speed=C0,
):
    """JSON-friendly radiation-pressure and force summary for a flat patch."""

    intensity = float(intensity_W_per_m2)
    area = float(area_m2)
    if intensity < 0.0:
        raise ValueError("intensity_W_per_m2 must be >= 0")
    if area < 0.0:
        raise ValueError("area_m2 must be >= 0")
    pressure = radiation_pressure_from_intensity(
        intensity,
        absorptance=absorptance,
        reflectance=reflectance,
        speed=speed,
    )
    return {
        "intensity_W_per_m2": intensity,
        "area_m2": area,
        "absorptance": float(absorptance),
        "reflectance": float(reflectance),
        "transmittance": 1.0 - float(absorptance) - float(reflectance),
        "speed_m_per_s": float(speed),
        "momentum_transfer_factor": float(absorptance) + 2.0 * float(reflectance),
        "incident_power_W": intensity * area,
        "pressure_Pa": pressure,
        "force_N": pressure * area,
    }


def oblique_radiation_pressure_summary(
    intensity_W_per_m2,
    incidence_angle_rad,
    area_m2=1.0,
    absorptance=1.0,
    reflectance=0.0,
    speed=C0,
):
    """Radiation force on a flat patch for oblique plane-wave incidence.

    ``incidence_angle_rad`` is measured from the surface normal: zero is normal
    incidence and ``pi/2`` is grazing.  The normal force scales with
    ``cos(angle)^2`` because both intercepted power and normal momentum carry a
    cosine factor.  Absorption also transfers tangential momentum; specular
    reflection reverses only the normal momentum component.
    """

    intensity = float(intensity_W_per_m2)
    angle = float(incidence_angle_rad)
    area = float(area_m2)
    speed = float(speed)
    if intensity < 0.0:
        raise ValueError("intensity_W_per_m2 must be >= 0")
    if angle < 0.0 or angle > 0.5 * math.pi:
        raise ValueError("incidence_angle_rad must be in [0, pi/2]")
    if area < 0.0:
        raise ValueError("area_m2 must be >= 0")
    if speed <= 0.0:
        raise ValueError("speed must be > 0")
    absorptance, reflectance = _validate_optical_coefficients(absorptance, reflectance)
    c = math.cos(angle)
    s = math.sin(angle)
    incident_power = intensity * area * c
    normal_pressure = (absorptance + 2.0 * reflectance) * intensity * c * c / speed
    tangential_pressure = absorptance * intensity * s * c / speed
    return {
        "intensity_W_per_m2": intensity,
        "incidence_angle_rad": angle,
        "incidence_angle_deg": math.degrees(angle),
        "area_m2": area,
        "absorptance": absorptance,
        "reflectance": reflectance,
        "transmittance": 1.0 - absorptance - reflectance,
        "speed_m_per_s": speed,
        "incident_power_on_patch_W": incident_power,
        "normal_pressure_Pa": normal_pressure,
        "tangential_pressure_Pa": tangential_pressure,
        "normal_force_N": normal_pressure * area,
        "tangential_force_N": tangential_pressure * area,
        "force_components_N": {
            "tangent": tangential_pressure * area,
            "normal": normal_pressure * area,
        },
        "normal_momentum_factor": absorptance + 2.0 * reflectance,
        "tangential_momentum_factor": absorptance,
    }


def poynting_patch_force_summary(
    poynting_W_per_m2,
    surface_normal,
    area_m2=1.0,
    absorptance=1.0,
    reflectance=0.0,
    speed=C0,
):
    """Vector radiation force on a flat patch from a Poynting vector.

    ``poynting_W_per_m2`` points along propagation.  ``surface_normal`` points
    out of the illuminated side of the patch, so an incident wave has
    ``dot(k, normal) < 0``.  The intercepted power is

        P_inc = |S| area max(0, -k.n).

    Absorption transfers momentum along ``k``; specular reflection reverses the
    normal component.  The returned force vector is in newtons.
    """

    s_vec = _float_vector(poynting_W_per_m2, "poynting_W_per_m2")
    normal = _unit_vector(surface_normal, "surface_normal")
    if len(s_vec) != len(normal):
        raise ValueError("poynting_W_per_m2 and surface_normal must have the same length")
    area = float(area_m2)
    speed = float(speed)
    if area < 0.0:
        raise ValueError("area_m2 must be >= 0")
    if speed <= 0.0:
        raise ValueError("speed must be > 0")
    absorptance, reflectance = _validate_optical_coefficients(absorptance, reflectance)

    s_mag = math.sqrt(sum(value * value for value in s_vec))
    if s_mag > 0.0:
        k = [value / s_mag for value in s_vec]
        cos_incidence = max(0.0, -sum(ki * ni for ki, ni in zip(k, normal)))
    else:
        k = [0.0 for _ in s_vec]
        cos_incidence = 0.0
    incident_intensity = s_mag * cos_incidence
    incident_power = incident_intensity * area
    absorption_scale = absorptance * incident_power / speed
    reflection_scale = -2.0 * reflectance * incident_power * cos_incidence / speed
    force = [
        absorption_scale * k[i] + reflection_scale * normal[i]
        for i in range(len(s_vec))
    ]
    normal_force_into_surface = -sum(fi * ni for fi, ni in zip(force, normal))
    tangential_force = [
        force[i] + normal_force_into_surface * normal[i]
        for i in range(len(force))
    ]
    return {
        "poynting_W_per_m2": s_vec,
        "poynting_magnitude_W_per_m2": s_mag,
        "propagation_direction": k,
        "surface_normal": normal,
        "area_m2": area,
        "absorptance": absorptance,
        "reflectance": reflectance,
        "transmittance": 1.0 - absorptance - reflectance,
        "speed_m_per_s": speed,
        "cos_incidence": cos_incidence,
        "incident_intensity_W_per_m2": incident_intensity,
        "incident_power_on_patch_W": incident_power,
        "force_N": force,
        "force_magnitude_N": math.sqrt(sum(value * value for value in force)),
        "normal_force_into_surface_N": normal_force_into_surface,
        "tangential_force_N": tangential_force,
        "tangential_force_magnitude_N": math.sqrt(
            sum(value * value for value in tangential_force)
        ),
    }


def electrostatic_eggshell_force(E, mesh, gradg, air_region="air"):
    """Weighted Maxwell-stress ("eggshell") ELECTROSTATIC force -- the electric twin
    of :func:`eggshell_force` (ε0 E in place of B/μ0). ``E`` is the electric field
    CoefficientFunction (``E = -grad(gfV)``); ``gradg`` = grad(g) of a smooth weight g
    (=1 on the body side, 0 on the far side of a band that lies in air), a vector CF
    nonzero only inside the band:

        F_k = - int_air [ eps0 E_k (E.gradg) - (eps0/2) |E|^2 d_k g ] dV

    Returns (Fx, Fy, Fz) in newtons. For a compact body use a spherical band
    (gradg = band * (r-center)/|r-center| * -1/(r_outer-r_inner)); for a plate/gap use
    an axis-aligned ramp (e.g. gradg = (0,0,g'(z)) across the gap).
    """
    Edg = InnerProduct(E, gradg)
    E2 = InnerProduct(E, E)
    region = dx(definedon=mesh.Materials(air_region))
    return tuple(-Integrate((EPS0 * E[k] * Edg - 0.5 * EPS0 * E2 * gradg[k]) * region, mesh)
                 for k in range(3))


def electrostatic_eggshell_force_2d(E, mesh, gradg, air_region="air"):
    """2D weighted Maxwell-stress ("eggshell") ELECTROSTATIC force [N/m] -- the 2D
    twin of :func:`electrostatic_eggshell_force` (range 2, per unit out-of-plane
    depth). ``E`` = ``-grad(V)`` (a 2-vector CF); ``gradg`` = grad of a smooth weight
    ``g`` (=1 on the body side, 0 on the far side across a band lying in the
    dielectric), a 2-vector CF nonzero only inside the band:

        F_k = - int [ eps0 E_k (E.gradg) - (eps0/2) |E|^2 d_k g ] dA .

    For a plate/gap use an axis-aligned ramp ``gradg = (0, g'(y))`` (a horizontal
    band in the gap); for a compact body a radial band. Returns (Fx, Fy) in N/m.

    WHY a volume band, not a boundary trace of ``grad(V)``: on a CONDUCTOR face
    ``V`` is constant, so the boundary trace of ``grad(V)`` (in a ``ds`` integral)
    keeps only the TANGENTIAL derivative (=0) and DROPS the normal field -- the
    surface-stress integral then reads ~0. The volume gradient is the true field, so
    the eggshell band (integrating in the dielectric AROUND the conductor) is the
    correct, robust extractor. (DEAD END recorded in ngsolve_usage("electro_mechanical").)
    """
    Edg = InnerProduct(E, gradg)
    E2 = InnerProduct(E, E)
    region = dx(definedon=mesh.Materials(air_region))
    return tuple(-Integrate((EPS0 * E[k] * Edg - 0.5 * EPS0 * E2 * gradg[k]) * region, mesh)
                 for k in range(2))


def magnetic_energy_2d(B, mesh, region=None):
    """2D magnetic field energy PER UNIT LENGTH [J/m]:  W = int |B|^2/(2 mu0) dA  over
    ``region`` (a material name; whole mesh if None).  For a current-driven problem the
    inductance is ``L = 2 W / I^2`` [H/m].  ``B`` = ``CF((grad(A)[1], -grad(A)[0]))``.

    NOTE: for a 1/r field (e.g. a coaxial line) the |B|^2 integrand is sharply peaked at
    small r, so refine the mesh near the inner radius -- the coax energy converged
    1.8 % -> 0.5 % under such refinement (examples/comsol_class/coax_line.py)."""
    dom = dx if region is None else dx(definedon=mesh.Materials(region))
    return Integrate(InnerProduct(B, B) / (2.0 * MU0) * dom, mesh)


def eggshell_force(B, mesh, center, r_inner, r_outer, air_region="air"):
    """Weighted Maxwell-stress ("eggshell") force on the body inside ``r_inner``.

    Robust on unstructured meshes: replaces the sharp surface integral by a
    volume integral over a radial band ``r_inner < |r-center| < r_outer`` with a
    smooth weight g (=1 at r_inner, 0 at r_outer):

        F_k = - int_band [ (1/mu0) B_k (B.grad g) - (1/2mu0) |B|^2 d_k g ] dV

    The band must lie in the air surrounding the body (choose ``r_inner`` >= the
    body radius). Returns (Fx, Fy, Fz) in newtons.
    """
    cx, cy, cz = center
    rho = sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
    band = IfPos(rho - r_inner, IfPos(r_outer - rho, 1.0, 0.0), 0.0)
    gscale = -1.0 / (r_outer - r_inner)
    gradg = band * CoefficientFunction((x - cx, y - cy, z - cz)) / rho * gscale
    Bdg = InnerProduct(B, gradg)
    B2 = InnerProduct(B, B)
    region = dx(definedon=mesh.Materials(air_region))
    F = []
    for k in range(3):
        integ = (1.0 / MU0) * B[k] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[k]
        F.append(-Integrate(integ * region, mesh))
    return tuple(F)


def eggshell_torque(B, mesh, center, r_inner, r_outer, pivot=(0.0, 0.0, 0.0),
                    air_region="air"):
    """3D weighted Maxwell-stress ("eggshell") TORQUE [N m] about ``pivot`` on
    the body inside ``r_inner``.  3D analogue of :func:`eggshell_torque_2d`:

        tau = - int_band  r' x S  dV,   r' = r - pivot,
        S_k = (1/mu0) B_k (B.grad g) - (1/2mu0) |B|^2 d_k g,

    same radial weight band (g=1 at r_inner, 0 at r_outer) in the air around the
    body. Returns (Tx, Ty, Tz). Use the same band as :func:`eggshell_force`;
    validated on a magnetised cylinder in a uniform field (tau = m x B0,
    examples/comsol_class/motor_torque.py)."""
    cx, cy, cz = center
    px, py, pz = pivot
    rho = sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
    band = IfPos(rho - r_inner, IfPos(r_outer - rho, 1.0, 0.0), 0.0)
    gscale = -1.0 / (r_outer - r_inner)
    gradg = band * CoefficientFunction((x - cx, y - cy, z - cz)) / rho * gscale
    Bdg = InnerProduct(B, gradg)
    B2 = InnerProduct(B, B)
    S = [(1.0 / MU0) * B[k] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[k] for k in range(3)]
    rp = (x - px, y - py, z - pz)
    cross = ((rp[1] * S[2] - rp[2] * S[1]),     # r' x S
             (rp[2] * S[0] - rp[0] * S[2]),
             (rp[0] * S[1] - rp[1] * S[0]))
    region = dx(definedon=mesh.Materials(air_region))
    return tuple(-Integrate(c * region, mesh) for c in cross)


def eggshell_force_2d(B, mesh, center, r_inner, r_outer, air_region="air"):
    """2D PLANAR weighted Maxwell-stress ("eggshell") force [N/m] on the body
    inside ``r_inner`` (force PER UNIT LENGTH in the out-of-plane direction).

    Same identity as :func:`eggshell_force` but in the (x, y) plane with the
    2-vector ``B = (Bx, By)`` (e.g. from the A_z solver,
    ``B = CF((grad(gfu)[1], -grad(gfu)[0]))``):

        F_k = - int_band [ (1/mu0) B_k (B.grad g) - (1/2mu0) |B|^2 d_k g ] dA

    The radial band ``r_inner < |r-center| < r_outer`` must lie in the air
    surrounding the body and enclose it. Validated on the two-parallel-wire
    benchmark F/L = mu0 I1 I2 / (2 pi d) to ~1 % (see tests/test_planar_force.py).
    Returns (Fx, Fy) in N/m.
    """
    cx, cy = center
    rho = sqrt((x - cx)**2 + (y - cy)**2)
    band = IfPos(rho - r_inner, IfPos(r_outer - rho, 1.0, 0.0), 0.0)
    gscale = -1.0 / (r_outer - r_inner)
    gradg = band * CoefficientFunction((x - cx, y - cy)) / rho * gscale
    Bdg = InnerProduct(B, gradg)
    B2 = InnerProduct(B, B)
    region = dx(definedon=mesh.Materials(air_region))
    F = []
    for k in range(2):
        integ = (1.0 / MU0) * B[k] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[k]
        F.append(-Integrate(integ * region, mesh))
    return tuple(F)


def eggshell_torque_2d(B, mesh, center, r_inner, r_outer, pivot=(0.0, 0.0),
                       air_region="air"):
    """2D PLANAR weighted Maxwell-stress torque [N] (per unit length) about
    ``pivot``, on the body inside ``r_inner``:

        tau_z = - int_band [ x' S_y - y' S_x ] dA,
        S_k = (1/mu0) B_k (B.grad g) - (1/2mu0)|B|^2 d_k g,   r' = r - pivot.

    Same eggshell band convention as :func:`eggshell_force_2d`. The lever-arm
    weighting is validated to reproduce r' x F of the validated force.
    """
    cx, cy = center
    px, py = pivot
    rho = sqrt((x - cx)**2 + (y - cy)**2)
    band = IfPos(rho - r_inner, IfPos(r_outer - rho, 1.0, 0.0), 0.0)
    gscale = -1.0 / (r_outer - r_inner)
    gradg = band * CoefficientFunction((x - cx, y - cy)) / rho * gscale
    Bdg = InnerProduct(B, gradg)
    B2 = InnerProduct(B, B)
    Sx = (1.0 / MU0) * B[0] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[0]
    Sy = (1.0 / MU0) * B[1] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[1]
    integ = (x - px) * Sy - (y - py) * Sx
    return -Integrate(integ * dx(definedon=mesh.Materials(air_region)), mesh)


def eggshell_force_axi(B, mesh, center, r_inner, r_outer, air_region="air"):
    """AXISYMMETRIC net AXIAL force [N] (full 3D torus) on the body whose
    meridional cross-section lies inside ``r_inner`` of ``center=(rc, zc)``,
    via the weighted Maxwell-stress ("eggshell") band in the (r, z) half-plane.

    ``B`` = ``CF((B_r, B_z))`` is the meridional flux density
    (``B_r = -grad(u)[1]``, ``B_z = grad(u)[0] + u/r`` from
    :func:`solve_axi_magnetostatic`).  The band ``r_inner < rho < r_outer``,
    ``rho = |(r, z) - center|``, must lie in the air enclosing the body's
    cross-section.  Only the AXIAL force is returned -- an axisymmetric body has
    no net radial force.  The ``2*pi*r`` toroidal weight makes this a 3D force:

        F_z = -2 pi int_band [ (1/mu0) B_z (B.grad g)
                               - (1/2mu0)|B|^2 d_z g ] r dr dz,

    g = 1 at r_inner, 0 at r_outer.  Validated on two coaxial loops against the
    exact mutual-inductance force I1 I2 dM/dz (M via elliptic integrals)
    (tests/test_axi_force.py)."""
    rc, zc = center
    rho = sqrt((x - rc)**2 + (y - zc)**2)
    band = IfPos(rho - r_inner, IfPos(r_outer - rho, 1.0, 0.0), 0.0)
    gscale = -1.0 / (r_outer - r_inner)
    gradg = band * CoefficientFunction((x - rc, y - zc)) / rho * gscale
    Bdg = InnerProduct(B, gradg)          # B.grad g
    B2 = InnerProduct(B, B)
    integ_z = (1.0 / MU0) * B[1] * Bdg - (1.0 / (2.0 * MU0)) * B2 * gradg[1]
    return -2.0 * math.pi * Integrate(
        integ_z * x * dx(definedon=mesh.Materials(air_region)), mesh)


def maxwell_surface_force(B, mesh, surface):
    """Maxwell-stress force as a surface integral over a named closed boundary
    ``surface`` in air enclosing the body (outward normal):

        F_k = oint (1/mu0) [ B_k (B.n) - 1/2 |B|^2 n_k ] dS

    Use ``eggshell_force`` instead unless you have a clean meshed surface --
    point/surface traces are noisier than the volume-band method. (Fx, Fy, Fz) [N].
    """
    n = specialcf.normal(mesh.dim)
    Bn = InnerProduct(B, n)
    B2 = InnerProduct(B, B)
    bnd = ds(definedon=mesh.Boundaries(surface))
    F = []
    for k in range(3):
        integ = (1.0 / MU0) * B[k] * Bn - (1.0 / (2.0 * MU0)) * B2 * n[k]
        F.append(Integrate(integ * bnd, mesh))
    return tuple(F)


def maxwell_surface_force_harmonic(B, mesh, surface):
    """TIME-AVERAGED Maxwell-stress force [N] over a closed boundary ``surface``
    for a COMPLEX time-harmonic flux density ``B`` (phasor).  The time-average of
    the quadratic Maxwell stress is

        <F_k> = oint_S [ (1/(2 mu0)) Re( B_k conj(B.n) )
                          - (1/(4 mu0)) |B|^2 n_k ] dS ,

    where the extra 1/2 vs the static :func:`maxwell_surface_force` is the
    time-average of cos^2(wt).  For a field oscillating with REAL peak amplitude
    B0 this equals exactly ``0.5 * maxwell_surface_force(B0)``.

    Intended for the time-harmonic eddy force on a body enclosed by an air
    surface -- in particular the ``"sibc"`` hole boundary of
    ``calc_fem_kelvin`` (the workpiece is a HOLE, so the only force handle is the
    Maxwell stress over its surface; ``B = curl(gfu)`` from the SIBC A-solve).
    Returns (Fx, Fy, Fz) in newtons.  Validated by reduction to the
    static :func:`maxwell_surface_force` reference (tests/test_maxwell_surface_harmonic.py).
    """
    n = specialcf.normal(mesh.dim)
    Bn = sum(B[k] * n[k] for k in range(3))                  # B . n  (n real)
    B2 = sum((B[k] * Conj(B[k])).real for k in range(3))     # |B|^2 (real)
    bnd = ds(definedon=mesh.Boundaries(surface))
    F = []
    for k in range(3):
        integ = (0.5 / MU0) * (B[k] * Conj(Bn)).real - (0.25 / MU0) * B2 * n[k]
        F.append(Integrate(integ * bnd, mesh))
    return tuple(F)


def ohmic_loss_2d(Ez, mesh, sigma, region=None):
    """Time-averaged ohmic loss PER UNIT LENGTH [W/m] for a 2D planar harmonic
    eddy problem:  P = 1/2 int sigma |E_z|^2 dA  (E_z = -j w A_z (+ Vc)).

    With a current-driven conductor (net current I), the AC resistance per length
    is ``Rac = 2 P / |I|^2``. Validated on the round-wire skin effect (Rac/Rdc vs
    Kelvin functions, 0.07 %; see tests/test_planar_eddy.py)."""
    integrand = 0.5 * sigma * (Ez * Conj(Ez)).real
    dom = dx if region is None else dx(definedon=mesh.Materials(region))
    return Integrate(integrand * dom, mesh)


def lorentz_force_2d(Jz, B, mesh, region):
    """2D PLANAR Lorentz force PER UNIT LENGTH [N/m] on the conductor ``region``
    carrying out-of-plane current density ``Jz`` [A/m^2] in flux density
    ``B = (Bx, By)``:

        F = int J x B dA = int Jz (zhat x B) dA   ->   Fx = -int Jz By,  Fy = int Jz Bx.

    Pass the TOTAL field B (the conductor's own self-field exerts ZERO net force by
    symmetry, so it drops out). The direct current-source twin of the Maxwell-stress
    :func:`eggshell_force_2d`. Returns ``(Fx, Fy)`` [N/m]. Validated on parallel
    busbars: |F| = mu0 I1 I2/(2 pi d) (examples/comsol_class/busbar_force.py)."""
    dom = dx(definedon=mesh.Materials(region))
    Fx = -Integrate(Jz * B[1] * dom, mesh)
    Fy = Integrate(Jz * B[0] * dom, mesh)
    return Fx, Fy


def magnetic_energy(B, mesh, region=None):
    """Field energy  W = 1/2 integral |B|^2 / mu0 dV  [J].
    ``region=None`` integrates the whole domain; else a material name."""
    integrand = 0.5 * InnerProduct(B, B) / MU0
    if region is None:
        return Integrate(integrand * dx, mesh)
    return Integrate(integrand * dx(definedon=mesh.Materials(region)), mesh)


def self_inductance(B, mesh, i_terminal):
    """Self-inductance  L = 2W / I^2  [H] from the field energy (I = terminal
    current = ampere-turns / N)."""
    return 2.0 * magnetic_energy(B, mesh) / (i_terminal * i_terminal)


def inductance_2d(B, mesh, nu, current, region=None):
    """2D PLANAR inductance PER UNIT LENGTH [H/m] via the energy method
    L = 2W/I^2,  W = 1/2 int nu |B|^2 dA  (``nu`` = reluctivity CF, B in-plane).

    ``region=None`` integrates the whole domain (total L); pass a material name
    for a partial energy (e.g. a conductor's internal inductance). Validated:
    round-wire internal inductance L_int = mu0/(8 pi) = 5.0e-8 H/m, radius-
    independent, to 0.06 % (tests/test_planar_inductance.py)."""
    dom = dx if region is None else dx(definedon=mesh.Materials(region))
    W = 0.5 * Integrate(nu * InnerProduct(B, B) * dom, mesh)
    return 2.0 * W / (current * current)


def inductance_axi(B, mesh, nu, current, region=None):
    """3D INDUCTANCE [H] of an axisymmetric coil via the energy method.

    L = 2 W_3D / I^2,  W_3D = 2*pi * int_half-plane  nu/2 |B|^2 r dr dz

    where the factor 2*pi converts the meridional-plane energy to the full
    toroidal (3D) energy.  B = (B_z, B_r) = (grad(u)[0]+u/r, -grad(u)[1])
    from ``solve_axi_magnetostatic``; ``nu`` = reluctivity; ``current`` = total
    current through the coil cross-section [A].

    ``region=None`` uses the whole domain; pass a material name to restrict."""
    dom = dx if region is None else dx(definedon=mesh.Materials(region))
    W_half = 0.5 * Integrate(nu * InnerProduct(B, B) * x * dom, mesh)
    W_3D = 2.0 * math.pi * W_half
    return 2.0 * W_3D / (current * current)


def ohmic_loss_axi(E_phi, mesh, sigma, region=None):
    """Time-averaged ohmic loss [W] (full 3D torus) for an axisymmetric
    time-harmonic eddy problem:

        P = 2*pi * int_half-plane  sigma/2 |E_phi|^2 r dr dz

    E_phi = -j*omega*A_phi (+ Vc/r for a driven conductor with NumberSpace Vc).
    AC resistance of a ring conductor: Rac = 2*P / |I|^2.

    Analogous to ``ohmic_loss_2d`` but accounts for the toroidal (2*pi*r) volume.
    """
    integrand = 0.5 * sigma * (E_phi * Conj(E_phi)).real * x
    dom = dx if region is None else dx(definedon=mesh.Materials(region))
    return 2.0 * math.pi * Integrate(integrand * dom, mesh)
