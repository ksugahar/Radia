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
import re
from datetime import datetime, timezone

from ngsolve import (CoefficientFunction, InnerProduct, sqrt, dx, ds, Integrate,
                     IfPos, specialcf, Conj, x, y, z)

from .scalar_fem3d import p1_surface_triangle_geometry, p1_tetrahedron_geometry

MU0 = 4.0e-7 * math.pi
EPS0 = 8.8541878128e-12
C0 = 299792458.0
ETA0 = MU0 * C0


def _parse_utc_like_datetime(value):
    if value in (None, ""):
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


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


def _metadata_truthy(value):
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    text = str(value).strip().lower()
    return text in {"1", "true", "yes", "y", "ok", "loaded", "solved", "pass", "passed"}


def electromagnetic_force_method_selection_gate(
    target_kind,
    requested_method,
    *,
    relative_permeability=1.0,
    weighted_stress_available=False,
    virtual_work_samples_available=False,
    contour_clearance_mesh_layers=0,
):
    """Select a robust electromagnetic-force extraction method.

    Unit-permeability current conductors prefer a Lorentz body-force integral.
    Magnetic bodies prefer a weighted-stress volume and fall back to coenergy
    virtual work when a displacement sweep exists. A contour Maxwell-stress
    integral remains a sensitivity diagnostic, not the primary acceptance
    method.
    """

    kind = str(target_kind or "").strip().lower().replace("-", "_").replace(" ", "_")
    method = str(requested_method or "").strip().lower().replace("-", "_").replace(" ", "_")
    mu_r = float(relative_permeability)
    contour_layers = int(contour_clearance_mesh_layers)
    conductor_kinds = {"conductor", "current_conductor", "coil", "busbar"}
    magnetic_kinds = {"magnetic_body", "ferromagnetic_body", "magnet", "iron"}

    if kind in conductor_kinds and abs(mu_r - 1.0) <= 1.0e-12:
        recommended = "lorentz_body_force"
        reason = "unit-mu current conductor supports a direct J cross B volume integral"
    elif kind in magnetic_kinds and bool(weighted_stress_available):
        recommended = "weighted_stress_volume"
        reason = "magnetic-body force should use a mesh-robust weighted stress volume"
    elif kind in magnetic_kinds and bool(virtual_work_samples_available):
        recommended = "coenergy_virtual_work"
        reason = "displacement samples support a coenergy derivative fallback"
    else:
        recommended = "none"
        reason = "no robust primary force evidence is available"

    contour_requested = method in {
        "contour_maxwell_stress",
        "maxwell_stress_contour",
        "line_maxwell_stress",
    }
    checks = {
        "target_kind_supported": kind in conductor_kinds | magnetic_kinds,
        "relative_permeability_positive": mu_r > 0.0,
        "robust_primary_method_available": recommended != "none",
        "requested_method_matches_recommendation": method == recommended,
        "contour_not_used_as_primary": not contour_requested,
        "contour_clearance_recorded_when_requested": not contour_requested or contour_layers > 0,
        "lorentz_restricted_to_unit_mu_conductor": (
            method != "lorentz_body_force"
            or (kind in conductor_kinds and abs(mu_r - 1.0) <= 1.0e-12)
        ),
    }
    return {
        "policy": "electromagnetic_force_method_selection_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "target_kind": kind,
        "relative_permeability": mu_r,
        "requested_method": method,
        "recommended_method": recommended,
        "recommendation_reason": reason,
        "weighted_stress_available": bool(weighted_stress_available),
        "virtual_work_samples_available": bool(virtual_work_samples_available),
        "contour_clearance_mesh_layers": contour_layers,
        "checks": checks,
        "notes": [
            "Use contour Maxwell stress as a sensitivity diagnostic, not the sole primary acceptance result.",
            "Cross-check magnetic-body forces with weighted stress and coenergy when both are available.",
        ],
    }


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


def maxwell_contour_segment_balance_summary_2d(
    vertices,
    B,
    mu=MU0,
    orientation="ccw",
    closed=True,
    expected_force_per_depth_N_per_m=None,
    force_abs_tolerance_N_per_m=1.0e-9,
):
    """Return a teaching-oriented balance audit for a 2D stress contour.

    The underlying contour integral is :func:`maxwell_contour_force_2d`.  This
    wrapper keeps the per-segment rows but adds a compact reference comparison,
    cancellation ratio, and orientation check so rectangular sanity cases can
    be read as a table rather than only as a net force.
    """

    tolerance = float(force_abs_tolerance_N_per_m)
    if tolerance < 0.0:
        raise ValueError("force_abs_tolerance_N_per_m must be >= 0")
    contour = maxwell_contour_force_2d(
        vertices,
        B,
        mu=mu,
        orientation=orientation,
        closed=closed,
    )
    segment_rows = []
    for row in contour["segments"]:
        tangent_force = row["tangential_force_per_depth_N_per_m"]
        tangent_magnitude = math.hypot(tangent_force[0], tangent_force[1])
        force = row["force_per_depth_N_per_m"]
        force_magnitude = math.hypot(force[0], force[1])
        normal_force = float(row["normal_force_per_depth_N_per_m"])
        segment_rows.append({
            "index": row["index"],
            "p0": row["p0"],
            "p1": row["p1"],
            "midpoint": row["midpoint"],
            "length_m": row["length_m"],
            "unit_normal": row["unit_normal"],
            "B_T": row["B_T"],
            "force_per_depth_N_per_m": force,
            "force_magnitude_per_depth_N_per_m": force_magnitude,
            "normal_force_per_depth_N_per_m": normal_force,
            "tangential_force_per_depth_N_per_m": tangent_force,
            "tangential_force_magnitude_per_depth_N_per_m": tangent_magnitude,
            "dominant_contribution": "normal" if abs(normal_force) >= tangent_magnitude else "tangential",
        })

    total = contour["total_force_per_depth_N_per_m"]
    total_magnitude = contour["total_force_magnitude_per_depth_N_per_m"]
    contribution_scale = (
        contour["sum_abs_normal_force_per_depth_N_per_m"]
        + contour["sum_abs_tangential_force_per_depth_N_per_m"]
    )
    cancellation_ratio = total_magnitude / contribution_scale if contribution_scale > 0.0 else 0.0

    expected = None
    abs_error = None
    max_abs_error = None
    reference_pass = None
    if expected_force_per_depth_N_per_m is not None:
        expected = _float_vector(expected_force_per_depth_N_per_m, "expected_force_per_depth_N_per_m")
        if len(expected) != 2:
            raise ValueError("expected_force_per_depth_N_per_m must have length 2")
        abs_error = [abs(total[i] - expected[i]) for i in range(2)]
        max_abs_error = max(abs_error)
        reference_pass = max_abs_error <= tolerance

    signed_area = contour["polygon_signed_area_m2"]
    orientation_consistent = None
    if closed and signed_area is not None:
        if orientation == "ccw":
            orientation_consistent = signed_area > tolerance
        else:
            orientation_consistent = signed_area < -tolerance

    issues = []
    if reference_pass is False:
        issues.append("net contour force differs from the supplied reference")
    if orientation_consistent is False:
        issues.append("vertex order sign does not match the requested orientation")

    dominant = max(
        segment_rows,
        key=lambda row: row["force_magnitude_per_depth_N_per_m"],
        default=None,
    )
    return {
        "policy": "maxwell_contour_segment_balance_2d",
        "status": "ok" if not issues else "needs_attention",
        "issues": issues,
        "closed": contour["closed"],
        "orientation": contour["orientation"],
        "orientation_consistent": orientation_consistent,
        "normal_side": contour["normal_side"],
        "n_segments": contour["n_segments"],
        "contour_length_m": contour["contour_length_m"],
        "polygon_signed_area_m2": signed_area,
        "total_force_per_depth_N_per_m": total,
        "total_force_magnitude_per_depth_N_per_m": total_magnitude,
        "sum_abs_normal_force_per_depth_N_per_m": contour["sum_abs_normal_force_per_depth_N_per_m"],
        "sum_abs_tangential_force_per_depth_N_per_m": contour["sum_abs_tangential_force_per_depth_N_per_m"],
        "cancellation_ratio": cancellation_ratio,
        "expected_force_per_depth_N_per_m": expected,
        "force_abs_tolerance_N_per_m": tolerance,
        "reference_force_abs_error_N_per_m": abs_error,
        "max_reference_force_abs_error_N_per_m": max_abs_error,
        "reference_pass": reference_pass,
        "dominant_segment_index": None if dominant is None else dominant["index"],
        "segment_rows": segment_rows,
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


def parallel_wire_virtual_work_force_summary(
    current1_A,
    current2_A,
    separation_xy_m,
    displacement_step_m=None,
    reference_separation_m=None,
):
    """Compare two-wire Lorentz force with fixed-current virtual work.

    The mutual magnetic coenergy per unit length has the separation-dependent
    part

        W'(d) = -mu0 I1 I2 / (2 pi) log(d / d_ref).

    Differentiating it with respect to the scalar separation ``d`` gives the
    radial force per unit length.  Like currents have a negative radial force
    in the increasing-separation coordinate, which is the same attraction
    direction reported by :func:`parallel_wire_lorentz_force_summary`.
    """

    pair = parallel_wire_lorentz_force_summary(current1_A, current2_A, separation_xy_m)
    distance = pair["separation_m"]
    if displacement_step_m is None:
        h = 1.0e-4 * distance
    else:
        h = float(displacement_step_m)
    if h <= 0.0:
        raise ValueError("displacement_step_m must be > 0")
    if h >= distance:
        raise ValueError("displacement_step_m must be smaller than the wire separation")
    if reference_separation_m is None:
        d_ref = distance
    else:
        d_ref = float(reference_separation_m)
        if d_ref <= 0.0:
            raise ValueError("reference_separation_m must be > 0")

    i1 = float(current1_A)
    i2 = float(current2_A)
    coefficient = MU0 * i1 * i2 / (2.0 * math.pi)

    def coenergy_per_length(d):
        return -coefficient * math.log(d / d_ref)

    e_minus = coenergy_per_length(distance - h)
    e_center = coenergy_per_length(distance)
    e_plus = coenergy_per_length(distance + h)
    virtual = virtual_work_symmetric_pair_force_summary(
        h,
        e_minus,
        e_plus,
        energy_kind="coenergy",
        center_position_m=distance,
        energy_center_J=e_center,
    )

    radial_force = virtual["force_N"]
    analytic_radial = -pair["signed_ampere_force_per_length_N_per_m"]
    unit = pair["unit_from_wire1_to_wire2"]
    virtual_force_on_wire2 = [radial_force * unit[0], radial_force * unit[1]]
    vector_error = [
        virtual_force_on_wire2[axis] - pair["force_on_wire2_N_per_m"][axis]
        for axis in range(2)
    ]
    vector_abs_error = math.hypot(vector_error[0], vector_error[1])
    reference_force = max(pair["force_magnitude_per_length_N_per_m"], 1.0e-300)
    radial_abs_error = abs(radial_force - analytic_radial)
    return {
        "current1_A": i1,
        "current2_A": i2,
        "separation_xy_m": pair["separation_xy_m"],
        "separation_m": distance,
        "reference_separation_m": d_ref,
        "displacement_step_m": h,
        "coenergy_formula": "-mu0*I1*I2/(2*pi)*log(d/d_ref) per unit length",
        "coenergy_minus_J_per_m": e_minus,
        "coenergy_center_J_per_m": e_center,
        "coenergy_plus_J_per_m": e_plus,
        "virtual_work_units_note": "energy samples are J/m, so the differentiated force is N/m",
        "virtual_work": virtual,
        "lorentz": pair,
        "virtual_work_radial_force_per_length_N_per_m": radial_force,
        "analytic_radial_force_per_length_N_per_m": analytic_radial,
        "virtual_work_force_on_wire2_N_per_m": virtual_force_on_wire2,
        "lorentz_force_on_wire2_N_per_m": pair["force_on_wire2_N_per_m"],
        "force_vector_abs_error_N_per_m": vector_abs_error,
        "radial_force_abs_error_N_per_m": radial_abs_error,
        "force_rel_error": vector_abs_error / reference_force,
        "interaction": pair["interaction"],
    }


def parallel_wire_force_result_package_gate(
    row,
    *,
    expected_model_id=None,
    expected_operating_point_id=None,
    expected_artifact_id=None,
    expected_result_set_id=None,
    expected_parameter_set_artifact_id=None,
    expected_parameter_set_digest=None,
    expected_parameter_set_path=None,
    expected_model_input_artifact_id=None,
    expected_model_input_digest=None,
    expected_model_input_path=None,
    expected_solution_artifact_id=None,
    expected_block_label_artifact_id=None,
    expected_source_tool=None,
    expected_source_group_id=None,
    expected_target_group_id=None,
    expected_source_center_xy_m=None,
    expected_target_center_xy_m=None,
    expected_source_region=None,
    expected_target_region=None,
    expected_source_material=None,
    expected_target_material=None,
    expected_postprocess_trace_id=None,
    expected_postprocess_command_digest=None,
    expected_postprocess_output_artifact_id=None,
    expected_postprocess_output_digest=None,
    expected_postprocess_output_schema_id=None,
    expected_postprocess_output_columns=None,
    expected_postprocess_output_units=None,
    expected_postprocess_row_convention_schema_id=None,
    expected_postprocess_script_artifact_id=None,
    expected_postprocess_script_digest=None,
    expected_postprocess_script_path=None,
    expected_force_observable_id=None,
    expected_force_observable_family=None,
    expected_force_convention_schema_id=None,
    expected_force_component_basis_schema_id=None,
    expected_force_unit_basis_schema_id=None,
    expected_objective_observable_id=None,
    expected_objective_observable_family=None,
    expected_force_component_frame=None,
    expected_radial_projection_axis=None,
    expected_force_sign_convention=None,
    expected_force_extraction_method=None,
    expected_block_integral_types=None,
    expected_current_source_artifact_id=None,
    expected_current_definition_method=None,
    expected_problem_type=None,
    expected_length_unit=None,
    expected_frequency_hz=None,
    expected_solver_precision=None,
    max_solver_precision=None,
    expected_min_angle_deg=None,
    expected_created_at_utc=None,
    expected_run_timestamp_utc=None,
    expected_solver_version=None,
    expected_radia_mcp_version=None,
    max_created_run_skew_s=None,
    min_timing_sections=4,
    require_solution_loaded=False,
    require_selection_clear=False,
    require_postprocess_command_trace=False,
    require_postprocess_output_artifact=False,
    require_postprocess_output_schema=False,
    require_postprocess_row_convention_schema=False,
    require_postprocess_script_artifact=False,
    require_force_convention_schema=False,
    require_force_component_basis_schema=False,
    require_force_unit_basis_schema=False,
    require_model_input_artifact=False,
    require_parameter_set_artifact=False,
    require_execution_metadata=False,
    require_timing_breakdown=False,
    rtol=1.0e-6,
):
    """Check a two-parallel-wire force row against Ampere's force law.

    This is the result-package form used when a FEMM, NGSolve, or notebook
    postprocess table reports a two-wire force.  The row must say which model
    and operating point it came from, record SI units per length, identify the
    result artifact/table, and preserve the sign convention.  Like currents
    attract: for a positive scalar separation the force on wire 2 points
    toward wire 1.
    """

    data = dict(row)
    tol = float(rtol)
    if tol < 0.0:
        raise ValueError("rtol must be >= 0")

    def _normalize_problem_type(value):
        text = str(value).strip().lower().replace("-", "_").replace(" ", "_")
        aliases = {
            "planar_2d": "planar",
            "2d_planar": "planar",
            "axi": "axisymmetric",
            "axisym": "axisymmetric",
            "axisymmetric_2d": "axisymmetric",
        }
        return aliases.get(text, text)

    def _normalize_length_unit(value):
        text = str(value).strip().lower().replace(" ", "").replace("_", "")
        aliases = {
            "m": "meters",
            "meter": "meters",
            "meters": "meters",
            "metre": "meters",
            "metres": "meters",
            "mm": "millimeters",
            "millimeter": "millimeters",
            "millimeters": "millimeters",
            "millimetre": "millimeters",
            "millimetres": "millimeters",
            "cm": "centimeters",
            "centimeter": "centimeters",
            "centimeters": "centimeters",
            "centimetre": "centimeters",
            "centimetres": "centimeters",
            "inch": "inches",
            "inches": "inches",
            "in": "inches",
        }
        return aliases.get(text, text)

    def _normalize_method(value):
        return str(value).strip().lower().replace("-", "_").replace(" ", "_")

    def _xy_pair(value):
        if value in (None, ""):
            return None
        if isinstance(value, dict):
            if "x" in value and "y" in value:
                return (float(value["x"]), float(value["y"]))
            if "X" in value and "Y" in value:
                return (float(value["X"]), float(value["Y"]))
        pair = list(value) if not isinstance(value, str) else [
            item for item in re.split(r"[,;\s]+", value.strip()) if item
        ]
        if len(pair) != 2:
            raise ValueError("wire center coordinates must contain exactly two values")
        return (float(pair[0]), float(pair[1]))

    model_id = str(data.get("model_id", "")).strip()
    operating_point_id = str(data.get("operating_point_id", "")).strip()
    artifact_id = str(data.get("artifact_id", data.get("case_artifact_id", ""))).strip()
    result_set_id = str(
        data.get("result_set_id", data.get("table_id", data.get("run_id", "")))
    ).strip()
    parameter_set_artifact_id = str(
        data.get(
            "parameter_set_artifact_id",
            data.get(
                "design_parameter_set_artifact_id",
                data.get(
                    "force_parameter_set_artifact_id",
                    data.get("operating_point_parameter_set_artifact_id", ""),
                ),
            ),
        )
    ).strip()
    parameter_set_digest = str(
        data.get(
            "parameter_set_digest",
            data.get(
                "parameter_set_sha256",
                data.get(
                    "design_parameter_set_digest",
                    data.get("force_parameter_set_digest", ""),
                ),
            ),
        )
    ).strip()
    parameter_set_path = str(
        data.get(
            "parameter_set_path",
            data.get(
                "parameter_set_file",
                data.get(
                    "design_parameter_set_path",
                    data.get("force_parameter_set_path", ""),
                ),
            ),
        )
    ).strip()
    model_input_artifact_id = str(
        data.get(
            "model_input_artifact_id",
            data.get(
                "fem_artifact_id",
                data.get("input_model_artifact_id", data.get("model_artifact_id", "")),
            ),
        )
    ).strip()
    model_input_digest = str(
        data.get(
            "model_input_digest",
            data.get(
                "fem_digest",
                data.get("input_model_digest", data.get("model_digest", "")),
            ),
        )
    ).strip()
    model_input_path = str(
        data.get(
            "model_input_path",
            data.get(
                "fem_path",
                data.get("input_model_path", data.get("model_path", "")),
            ),
        )
    ).strip()
    solution_artifact_id = str(
        data.get(
            "solution_artifact_id",
            data.get("ans_artifact_id", data.get("loaded_solution_artifact_id", "")),
        )
    ).strip()
    block_label_artifact_id = str(
        data.get(
            "block_label_artifact_id",
            data.get(
                "block_labels_artifact_id",
                data.get("source_contract_artifact_id", ""),
            ),
        )
    ).strip()
    solution_loaded_raw = data.get(
        "solution_loaded",
        data.get("postprocessor_solution_loaded", data.get("mo_solution_loaded")),
    )
    solution_loaded_recorded = solution_loaded_raw not in (None, "")
    solution_loaded = _metadata_truthy(solution_loaded_raw)
    source_tool = str(data.get("source_tool", "")).strip()
    source_function = str(data.get("source_function", "")).strip()
    source_group_id = str(
        data.get(
            "source_group_id",
            data.get("current_group_id", data.get("wire1_group_id", "")),
        )
    ).strip()
    target_group_id = str(
        data.get(
            "target_group_id",
            data.get(
                "force_group_id",
                data.get("block_integral_group_id", data.get("wire2_group_id", "")),
            ),
        )
    ).strip()
    source_center_xy_m = _xy_pair(
        data.get(
            "source_center_xy_m",
            data.get("wire1_center_xy_m", data.get("source_center_m")),
        )
    )
    target_center_xy_m = _xy_pair(
        data.get(
            "target_center_xy_m",
            data.get("wire2_center_xy_m", data.get("target_center_m")),
        )
    )
    selection_function = str(
        data.get("selection_function", data.get("block_integral_selection", ""))
    ).strip()
    selection_lower = selection_function.lower()
    postprocess_trace_id = str(
        data.get(
            "postprocess_trace_id",
            data.get("postprocess_command_trace_id", data.get("command_trace_id", "")),
        )
    ).strip()
    postprocess_command_digest = str(
        data.get(
            "postprocess_command_digest",
            data.get("command_trace_sha256", data.get("selection_command_digest", "")),
        )
    ).strip()
    postprocess_output_artifact_id = str(
        data.get(
            "postprocess_output_artifact_id",
            data.get(
                "postprocess_table_artifact_id",
                data.get("postprocess_result_artifact_id", data.get("output_artifact_id", "")),
            ),
        )
    ).strip()
    postprocess_output_digest = str(
        data.get(
            "postprocess_output_digest",
            data.get(
                "postprocess_output_sha256",
                data.get("postprocess_table_sha256", data.get("output_digest", "")),
            ),
        )
    ).strip()
    postprocess_output_path = str(
        data.get(
            "postprocess_output_path",
            data.get("postprocess_table_path", data.get("result_table_path", "")),
        )
    ).strip()
    postprocess_output_schema_id = str(
        data.get(
            "postprocess_output_schema_id",
            data.get(
                "postprocess_table_schema_id",
                data.get("force_table_schema_id", data.get("output_schema_id", "")),
            ),
        )
    ).strip()
    postprocess_row_convention_schema_id = str(
        data.get(
            "postprocess_row_convention_schema_id",
            data.get(
                "postprocess_convention_schema_id",
                data.get(
                    "force_postprocess_convention_schema_id",
                    data.get("row_convention_schema_id", ""),
                ),
            ),
        )
    ).strip()

    def _normalize_output_columns(value):
        if value in (None, ""):
            return []
        if isinstance(value, str):
            return [item.strip() for item in re.split(r"[,;\s]+", value) if item.strip()]
        return [str(item).strip() for item in value if str(item).strip()]

    def _normalize_output_units(value):
        if value in (None, ""):
            return {}
        if isinstance(value, dict):
            return {
                str(key).strip(): str(val).strip()
                for key, val in value.items()
                if str(key).strip()
            }
        rows = list(value)
        return {
            str(index): str(val).strip()
            for index, val in enumerate(rows)
            if str(val).strip()
        }

    postprocess_output_columns = _normalize_output_columns(
        data.get(
            "postprocess_output_columns",
            data.get(
                "postprocess_table_columns",
                data.get("force_table_columns", data.get("output_columns", data.get("columns"))),
            ),
        )
    )
    postprocess_output_units = _normalize_output_units(
        data.get(
            "postprocess_output_units",
            data.get(
                "postprocess_table_units",
                data.get("force_table_units", data.get("output_units", data.get("column_units"))),
            ),
        )
    )
    postprocess_script_artifact_id = str(
        data.get(
            "postprocess_script_artifact_id",
            data.get(
                "postprocess_script_id",
                data.get("script_artifact_id", data.get("postprocess_file_artifact_id", "")),
            ),
        )
    ).strip()
    postprocess_script_digest = str(
        data.get(
            "postprocess_script_digest",
            data.get(
                "postprocess_script_sha256",
                data.get("script_digest", data.get("script_sha256", data.get("postprocess_file_digest", ""))),
            ),
        )
    ).strip()
    postprocess_script_path = str(
        data.get(
            "postprocess_script_path",
            data.get("postprocess_script_file", data.get("script_path", data.get("postprocess_file_path", ""))),
        )
    ).strip()
    force_observable_id = str(
        data.get(
            "force_observable_id",
            data.get("force_integral_observable_id", data.get("observable_id", "")),
        )
    ).strip()
    force_observable_family = str(
        data.get(
            "force_observable_family",
            data.get("force_integral_family", data.get("observable_family", "")),
        )
    ).strip()
    force_convention_schema_id = str(
        data.get(
            "force_convention_schema_id",
            data.get(
                "force_physics_convention_schema_id",
                data.get(
                    "source_material_convention_schema_id",
                    data.get("physics_convention_schema_id", ""),
                ),
            ),
        )
    ).strip()
    force_component_basis_schema_id = str(
        data.get(
            "force_component_basis_schema_id",
            data.get(
                "component_basis_schema_id",
                data.get(
                    "force_component_convention_schema_id",
                    data.get("force_component_row_convention_schema_id", ""),
                ),
            ),
        )
    ).strip()
    force_unit_basis_schema_id = str(
        data.get(
            "force_unit_basis_schema_id",
            data.get(
                "unit_basis_schema_id",
                data.get(
                    "femm_force_unit_basis_schema_id",
                    data.get("planar_force_unit_basis_schema_id", ""),
                ),
            ),
        )
    ).strip()
    objective_observable_id = str(
        data.get(
            "objective_observable_id",
            data.get(
                "force_objective_observable_id",
                data.get("objective_id", data.get("objective_function_id", "")),
            ),
        )
    ).strip()
    objective_observable_family = str(
        data.get(
            "objective_observable_family",
            data.get(
                "force_objective_observable_family",
                data.get("objective_family", data.get("objective_kind", "")),
            ),
        )
    ).strip()
    command_sequence_source = data.get(
        "postprocess_commands",
        data.get("postprocess_command_sequence", data.get("command_sequence")),
    )
    if command_sequence_source is None:
        postprocess_commands = []
    elif isinstance(command_sequence_source, str):
        postprocess_commands = [
            item.strip()
            for item in command_sequence_source.replace("\n", ";").split(";")
            if item.strip()
        ]
    else:
        postprocess_commands = [str(item).strip() for item in command_sequence_source if str(item).strip()]
    command_trace_text = "; ".join(postprocess_commands).lower()
    source_region = str(
        data.get(
            "source_region",
            data.get("source_region_name", data.get("wire1_region", "")),
        )
    ).strip()
    target_region = str(
        data.get(
            "target_region",
            data.get("target_region_name", data.get("wire2_region", "")),
        )
    ).strip()
    source_material = str(
        data.get(
            "source_material",
            data.get("source_material_name", data.get("wire1_material", "")),
        )
    ).strip()
    target_material = str(
        data.get(
            "target_material",
            data.get("target_material_name", data.get("wire2_material", "")),
        )
    ).strip()
    force_component_frame = str(
        data.get("force_component_frame", data.get("component_frame", ""))
    ).strip()
    radial_projection_axis = str(
        data.get(
            "radial_projection_axis",
            data.get("radial_axis", data.get("projection_axis", "")),
        )
    ).strip()
    force_sign_convention = str(
        data.get(
            "force_sign_convention",
            data.get("sign_convention", data.get("radial_force_sign_convention", "")),
        )
    ).strip()
    force_extraction_method = str(
        data.get(
            "force_extraction_method",
            data.get(
                "extraction_method",
                data.get("postprocess_method", data.get("force_method", "")),
            ),
        )
    ).strip()
    force_extraction_method_normalized = (
        _normalize_method(force_extraction_method) if force_extraction_method else ""
    )
    current_source_artifact_id = str(
        data.get(
            "current_source_artifact_id",
            data.get(
                "current_snapshot_artifact_id",
                data.get("current_definition_artifact_id", data.get("current_artifact_id", "")),
            ),
        )
    ).strip()
    current_definition_method = str(
        data.get(
            "current_definition_method",
            data.get("current_method", data.get("current_kind", "")),
        )
    ).strip()
    current_definition_method_normalized = (
        _normalize_method(current_definition_method) if current_definition_method else ""
    )
    block_integral_source = data.get(
        "block_integral_types",
        data.get("block_integral_type", data.get("integral_types")),
    )
    if block_integral_source is None:
        block_integral_text = f"{selection_function} {source_function}"
        block_integral_types = [
            int(match)
            for match in re.findall(r"blockintegral\s*\(\s*(\d+)\s*\)", block_integral_text, flags=re.I)
        ]
    elif isinstance(block_integral_source, (list, tuple, set)):
        block_integral_types = [int(value) for value in block_integral_source]
    else:
        block_integral_types = [
            int(value)
            for value in re.findall(r"\d+", str(block_integral_source))
        ]
    units = str(data.get("force_units", data.get("unit", ""))).strip()
    units_compact = units.lower().replace(" ", "")
    raw_basis = str(data.get("force_unit_basis", data.get("unit_basis", ""))).strip()
    basis_key = raw_basis.lower().replace("-", "_").replace(" ", "_").replace("/", "_")
    if not basis_key:
        if "n/m" in units_compact or "n_per_m" in units_compact:
            force_unit_basis = "per_length"
        elif units_compact == "n":
            force_unit_basis = "depth_integrated"
        else:
            force_unit_basis = ""
    elif basis_key in {
        "per_length",
        "per_unit_length",
        "per_depth",
        "per_unit_depth",
        "force_per_length",
        "n_per_m",
    }:
        force_unit_basis = "per_length"
    elif basis_key in {
        "depth_integrated",
        "depth_integrated_total",
        "total_depth_integrated",
        "total_force",
        "total",
        "n",
    }:
        force_unit_basis = "depth_integrated"
    else:
        force_unit_basis = raw_basis
    planar_depth_source = data.get(
        "problem_depth_m",
        data.get("planar_depth_m", data.get("femm_problem_depth_m", data.get("depth_m"))),
    )
    planar_depth_m = None
    if planar_depth_source not in (None, ""):
        planar_depth_m = float(planar_depth_source)
    if "separation_xy_m" in data:
        separation = data["separation_xy_m"]
    else:
        separation = data.get("separation_m")
    problem_type = str(
        data.get(
            "problem_type",
            data.get("femm_problem_type", data.get("analysis_type", "")),
        )
    ).strip()
    problem_type_normalized = _normalize_problem_type(problem_type) if problem_type else ""
    length_unit = str(
        data.get(
            "length_unit",
            data.get("length_units", data.get("problem_length_unit", data.get("femm_length_unit", ""))),
        )
    ).strip()
    length_unit_normalized = _normalize_length_unit(length_unit) if length_unit else ""
    frequency_source = data.get(
        "frequency_hz",
        data.get("problem_frequency_hz", data.get("femm_frequency_hz", data.get("freq_hz"))),
    )
    frequency_hz = None
    if frequency_source not in (None, ""):
        frequency_hz = float(frequency_source)
    solver_precision_source = data.get(
        "solver_precision",
        data.get(
            "problem_precision",
            data.get("femm_solver_precision", data.get("mi_probdef_precision")),
        ),
    )
    solver_precision = None
    if solver_precision_source not in (None, ""):
        solver_precision = float(solver_precision_source)
    min_angle_source = data.get(
        "min_angle_deg",
        data.get(
            "minangle_deg",
            data.get("triangle_min_angle_deg", data.get("femm_minangle_deg")),
        ),
    )
    min_angle_deg = None
    if min_angle_source not in (None, ""):
        min_angle_deg = float(min_angle_source)
    created_at_utc = str(
        data.get(
            "created_at_utc",
            data.get("artifact_created_at_utc", data.get("created_at", "")),
        )
    ).strip()
    run_timestamp_utc = str(
        data.get(
            "run_timestamp_utc",
            data.get(
                "executed_at_utc",
                data.get("run_date_utc", data.get("date_utc", data.get("run_date", ""))),
            ),
        )
    ).strip()
    solver_version = str(
        data.get(
            "solver_version",
            data.get("femm_version", data.get("source_tool_version", "")),
        )
    ).strip()
    radia_mcp_version = str(
        data.get(
            "radia_mcp_version",
            data.get("radia_ngsolve_version", data.get("mcp_server_version", "")),
        )
    ).strip()
    run_duration_source = data.get(
        "run_duration_s",
        data.get("elapsed_s", data.get("runtime_s", data.get("wall_time_s"))),
    )
    run_duration_s = None
    if run_duration_source not in (None, ""):
        run_duration_s = float(run_duration_source)

    def _timing_duration(value):
        if isinstance(value, dict):
            for key in ("duration_s", "elapsed_s", "runtime_s", "seconds", "s"):
                if value.get(key) not in (None, ""):
                    return float(value[key])
            return None
        if value in (None, ""):
            return None
        return float(value)

    def _timing_rows(value):
        if value in (None, ""):
            return []
        rows = []
        if isinstance(value, dict):
            iterable = value.items()
        else:
            iterable = enumerate(value)
        for key, entry in iterable:
            if isinstance(entry, dict):
                name = str(
                    entry.get(
                        "name",
                        entry.get("stage", entry.get("phase", entry.get("label", key))),
                    )
                ).strip()
                duration = _timing_duration(entry)
            elif isinstance(entry, (list, tuple)) and len(entry) >= 2:
                name = str(entry[0]).strip()
                duration = _timing_duration(entry[1])
            else:
                name = str(key).strip()
                duration = _timing_duration(entry)
            if name and duration is not None:
                rows.append({"name": name, "duration_s": duration})
        return rows

    timing_breakdown_rows = _timing_rows(
        data.get(
            "timing_breakdown_s",
            data.get("timing_breakdown", data.get("timing_breakdown_rows", data.get("timings"))),
        )
    )
    timing_durations = [row["duration_s"] for row in timing_breakdown_rows]
    timing_total_s = sum(timing_durations) if timing_durations else None
    timing_sections_min = int(min_timing_sections or 0)
    top_durations = timing_durations[: max(timing_sections_min, 1)]
    timing_top_sections_descending = all(
        top_durations[index] >= top_durations[index + 1]
        for index in range(len(top_durations) - 1)
    )
    analytic = parallel_wire_lorentz_force_summary(
        data.get("current1_A"),
        data.get("current2_A"),
        separation,
    )
    expected_vector = analytic["force_on_wire2_N_per_m"]
    expected_radial = -analytic["signed_ampere_force_per_length_N_per_m"]

    measured_vector = data.get("force_on_wire2_N_per_m")
    measured_radial = data.get("radial_force_on_wire2_N_per_m")
    vector_error = math.nan
    radial_error = math.nan
    vector_ok = True
    radial_ok = True
    reference_force = max(analytic["force_magnitude_per_length_N_per_m"], 1.0e-300)
    if measured_vector is not None:
        vector = [float(value) for value in measured_vector]
        if len(vector) != 2:
            raise ValueError("force_on_wire2_N_per_m must have length 2")
        vector_error = math.hypot(
            vector[0] - expected_vector[0],
            vector[1] - expected_vector[1],
        )
        vector_ok = vector_error <= tol * reference_force
    if measured_radial is not None:
        radial = float(measured_radial)
        radial_error = abs(radial - expected_radial)
        radial_ok = radial_error <= tol * reference_force
    interaction = str(data.get("interaction", "")).strip().lower()
    expected_model = None if expected_model_id is None else str(expected_model_id).strip()
    expected_op = (
        None if expected_operating_point_id is None else str(expected_operating_point_id).strip()
    )
    expected_artifact = (
        None if expected_artifact_id is None else str(expected_artifact_id).strip()
    )
    expected_result_set = (
        None if expected_result_set_id is None else str(expected_result_set_id).strip()
    )
    expected_parameter_set_artifact = (
        None
        if expected_parameter_set_artifact_id is None
        else str(expected_parameter_set_artifact_id).strip()
    )
    expected_parameter_set_digest_text = (
        None
        if expected_parameter_set_digest is None
        else str(expected_parameter_set_digest).strip()
    )
    expected_parameter_set_path_text = (
        None
        if expected_parameter_set_path is None
        else str(expected_parameter_set_path).strip()
    )
    expected_model_input_artifact = (
        None
        if expected_model_input_artifact_id is None
        else str(expected_model_input_artifact_id).strip()
    )
    expected_model_input_digest_text = (
        None
        if expected_model_input_digest is None
        else str(expected_model_input_digest).strip()
    )
    expected_model_input_path_text = (
        None
        if expected_model_input_path is None
        else str(expected_model_input_path).strip()
    )
    expected_solution_artifact = (
        None
        if expected_solution_artifact_id is None
        else str(expected_solution_artifact_id).strip()
    )
    expected_block_label_artifact = (
        None
        if expected_block_label_artifact_id is None
        else str(expected_block_label_artifact_id).strip()
    )
    expected_tool = None if expected_source_tool is None else str(expected_source_tool).strip()
    expected_source_group = (
        None if expected_source_group_id is None else str(expected_source_group_id).strip()
    )
    expected_target_group = (
        None if expected_target_group_id is None else str(expected_target_group_id).strip()
    )
    expected_source_center = _xy_pair(expected_source_center_xy_m)
    expected_target_center = _xy_pair(expected_target_center_xy_m)
    expected_source_region_text = (
        None if expected_source_region is None else str(expected_source_region).strip()
    )
    expected_target_region_text = (
        None if expected_target_region is None else str(expected_target_region).strip()
    )
    expected_source_material_text = (
        None if expected_source_material is None else str(expected_source_material).strip()
    )
    expected_target_material_text = (
        None if expected_target_material is None else str(expected_target_material).strip()
    )
    expected_trace_id = (
        None if expected_postprocess_trace_id is None else str(expected_postprocess_trace_id).strip()
    )
    expected_command_digest = (
        None
        if expected_postprocess_command_digest is None
        else str(expected_postprocess_command_digest).strip()
    )
    expected_output_artifact = (
        None
        if expected_postprocess_output_artifact_id is None
        else str(expected_postprocess_output_artifact_id).strip()
    )
    expected_output_digest = (
        None
        if expected_postprocess_output_digest is None
        else str(expected_postprocess_output_digest).strip()
    )
    expected_output_schema_id = (
        None
        if expected_postprocess_output_schema_id is None
        else str(expected_postprocess_output_schema_id).strip()
    )
    expected_output_columns = (
        None
        if expected_postprocess_output_columns is None
        else _normalize_output_columns(expected_postprocess_output_columns)
    )
    expected_output_units = (
        None
        if expected_postprocess_output_units is None
        else _normalize_output_units(expected_postprocess_output_units)
    )
    expected_postprocess_row_convention_schema = (
        None
        if expected_postprocess_row_convention_schema_id is None
        else str(expected_postprocess_row_convention_schema_id).strip()
    )
    expected_script_artifact = (
        None
        if expected_postprocess_script_artifact_id is None
        else str(expected_postprocess_script_artifact_id).strip()
    )
    expected_script_digest = (
        None
        if expected_postprocess_script_digest is None
        else str(expected_postprocess_script_digest).strip()
    )
    expected_script_path = (
        None
        if expected_postprocess_script_path is None
        else str(expected_postprocess_script_path).strip()
    )
    expected_observable_id = (
        None if expected_force_observable_id is None else str(expected_force_observable_id).strip()
    )
    expected_observable_family = (
        None if expected_force_observable_family is None else str(expected_force_observable_family).strip()
    )
    expected_force_convention_schema = (
        None
        if expected_force_convention_schema_id is None
        else str(expected_force_convention_schema_id).strip()
    )
    expected_force_component_basis_schema = (
        None
        if expected_force_component_basis_schema_id is None
        else str(expected_force_component_basis_schema_id).strip()
    )
    expected_force_unit_basis_schema = (
        None
        if expected_force_unit_basis_schema_id is None
        else str(expected_force_unit_basis_schema_id).strip()
    )
    expected_objective_id = (
        None
        if expected_objective_observable_id is None
        else str(expected_objective_observable_id).strip()
    )
    expected_objective_family = (
        None
        if expected_objective_observable_family is None
        else str(expected_objective_observable_family).strip()
    )
    expected_component_frame = (
        None if expected_force_component_frame is None else str(expected_force_component_frame).strip()
    )
    expected_projection_axis = (
        None if expected_radial_projection_axis is None else str(expected_radial_projection_axis).strip()
    )
    expected_sign_convention = (
        None if expected_force_sign_convention is None else str(expected_force_sign_convention).strip()
    )
    expected_extraction_method = (
        None
        if expected_force_extraction_method is None
        else _normalize_method(expected_force_extraction_method)
    )
    expected_integral_types = (
        None
        if expected_block_integral_types is None
        else sorted(int(value) for value in expected_block_integral_types)
    )
    expected_current_artifact = (
        None
        if expected_current_source_artifact_id is None
        else str(expected_current_source_artifact_id).strip()
    )
    expected_current_method = (
        None
        if expected_current_definition_method is None
        else _normalize_method(expected_current_definition_method)
    )
    expected_problem = (
        None if expected_problem_type is None else _normalize_problem_type(expected_problem_type)
    )
    expected_length = (
        None if expected_length_unit is None else _normalize_length_unit(expected_length_unit)
    )
    expected_frequency = (
        None if expected_frequency_hz is None else float(expected_frequency_hz)
    )
    expected_precision = (
        None if expected_solver_precision is None else float(expected_solver_precision)
    )
    max_precision = None if max_solver_precision is None else float(max_solver_precision)
    expected_min_angle = (
        None if expected_min_angle_deg is None else float(expected_min_angle_deg)
    )
    expected_created_at = (
        None
        if expected_created_at_utc is None
        else str(expected_created_at_utc).strip()
    )
    expected_run_timestamp = (
        None
        if expected_run_timestamp_utc is None
        else str(expected_run_timestamp_utc).strip()
    )
    expected_solver_version_text = (
        None if expected_solver_version is None else str(expected_solver_version).strip()
    )
    expected_radia_mcp_version_text = (
        None if expected_radia_mcp_version is None else str(expected_radia_mcp_version).strip()
    )
    frequency_tolerance = (
        None
        if expected_frequency is None
        else max(1.0e-12, tol * max(abs(expected_frequency), 1.0))
    )
    precision_tolerance = (
        None
        if expected_precision is None
        else max(1.0e-18, tol * max(abs(expected_precision), 1.0))
    )
    min_angle_tolerance = (
        None
        if expected_min_angle is None
        else max(1.0e-12, tol * max(abs(expected_min_angle), 1.0))
    )
    max_created_run_skew = (
        None if max_created_run_skew_s is None else float(max_created_run_skew_s)
    )
    created_at_dt = _parse_utc_like_datetime(created_at_utc)
    run_timestamp_dt = _parse_utc_like_datetime(run_timestamp_utc)
    created_run_skew_s = None
    if created_at_dt is not None and run_timestamp_dt is not None:
        created_run_skew_s = abs((created_at_dt - run_timestamp_dt).total_seconds())
    separation_float = float(analytic["separation_m"])
    center_distance_m = None
    center_distance_error_m = None
    center_tolerance_m = max(1.0e-12, tol * max(abs(separation_float), 1.0))
    if source_center_xy_m is not None and target_center_xy_m is not None:
        dx = target_center_xy_m[0] - source_center_xy_m[0]
        dy = target_center_xy_m[1] - source_center_xy_m[1]
        center_distance_m = math.hypot(dx, dy)
        center_distance_error_m = abs(center_distance_m - separation_float)

    def _center_matches(actual, expected):
        if expected is None:
            return True
        if actual is None:
            return False
        return (
            abs(actual[0] - expected[0]) <= center_tolerance_m
            and abs(actual[1] - expected[1]) <= center_tolerance_m
        )

    solution_loaded_required = bool(require_solution_loaded)
    selection_clear_required = bool(require_selection_clear)
    trace_required = bool(require_postprocess_command_trace)
    execution_metadata_required = (
        bool(require_execution_metadata)
        or expected_created_at is not None
        or expected_run_timestamp is not None
        or expected_solver_version_text is not None
        or expected_radia_mcp_version_text is not None
        or max_created_run_skew is not None
    )
    created_at_required = expected_created_at is not None or max_created_run_skew is not None
    timing_breakdown_required = bool(require_timing_breakdown)
    parameter_set_artifact_required = (
        bool(require_parameter_set_artifact)
        or expected_parameter_set_artifact is not None
        or expected_parameter_set_digest_text is not None
        or expected_parameter_set_path_text is not None
    )
    parameter_set_digest_required = (
        bool(require_parameter_set_artifact)
        or expected_parameter_set_digest_text is not None
    )
    parameter_set_path_required = (
        bool(require_parameter_set_artifact)
        or expected_parameter_set_path_text is not None
    )
    output_artifact_required = (
        bool(require_postprocess_output_artifact)
        or expected_output_artifact is not None
        or expected_output_digest is not None
    )
    output_schema_required = (
        bool(require_postprocess_output_schema)
        or expected_output_schema_id is not None
        or expected_output_columns is not None
        or expected_output_units is not None
    )
    postprocess_row_convention_schema_required = (
        bool(require_postprocess_row_convention_schema)
        or expected_postprocess_row_convention_schema is not None
    )
    force_convention_schema_required = (
        bool(require_force_convention_schema)
        or expected_force_convention_schema is not None
    )
    force_component_basis_schema_required = (
        bool(require_force_component_basis_schema)
        or expected_force_component_basis_schema is not None
    )
    force_unit_basis_schema_required = (
        bool(require_force_unit_basis_schema)
        or expected_force_unit_basis_schema is not None
    )
    script_artifact_required = (
        bool(require_postprocess_script_artifact)
        or expected_script_artifact is not None
        or expected_script_digest is not None
        or expected_script_path is not None
    )
    script_digest_required = (
        bool(require_postprocess_script_artifact)
        or expected_script_digest is not None
    )
    script_path_required = (
        bool(require_postprocess_script_artifact)
        or expected_script_path is not None
    )
    model_input_artifact_required = (
        bool(require_model_input_artifact)
        or expected_model_input_artifact is not None
        or expected_model_input_digest_text is not None
        or expected_model_input_path_text is not None
    )
    model_input_digest_required = (
        bool(require_model_input_artifact)
        or expected_model_input_digest_text is not None
    )
    model_input_path_required = (
        bool(require_model_input_artifact)
        or expected_model_input_path_text is not None
    )
    output_digest_required = bool(require_postprocess_output_artifact) or expected_output_digest is not None
    trace_id_required = trace_required or expected_trace_id is not None
    command_digest_required = trace_required or expected_command_digest is not None
    source_group_required = expected_source_group is not None
    target_group_required = expected_target_group is not None
    current_artifact_required = expected_current_artifact is not None
    current_method_required = expected_current_method is not None
    selection_mentions_target_group = True
    if expected_target_group is not None:
        pattern = rf"(?<![A-Za-z0-9_]){re.escape(expected_target_group)}(?![A-Za-z0-9_])"
        selection_mentions_target_group = bool(
            selection_function and re.search(pattern, selection_function)
        )
    clear_index = selection_lower.find("mo_clearblock")
    select_index = selection_lower.find("mo_groupselectblock")
    selection_clear_before_select = (
        clear_index >= 0
        and select_index >= 0
        and clear_index < select_index
    )
    trace_clear_index = command_trace_text.find("mo_clearblock")
    trace_select_index = command_trace_text.find("mo_groupselectblock")
    trace_integral18_index = command_trace_text.find("mo_blockintegral(18")
    trace_integral19_index = command_trace_text.find("mo_blockintegral(19")
    trace_has_clear_select_force_xy = (
        trace_clear_index >= 0
        and trace_select_index >= 0
        and trace_integral18_index >= 0
        and trace_integral19_index >= 0
        and trace_clear_index < trace_select_index
        and trace_select_index < min(trace_integral18_index, trace_integral19_index)
    )
    trace_mentions_target_group = True
    if expected_target_group is not None and postprocess_commands:
        trace_mentions_target_group = bool(
            re.search(pattern, command_trace_text)
        )
    has_force_measurement = measured_vector is not None or measured_radial is not None
    femm_source = source_tool.lower().startswith("femm")
    radial_axis_text = radial_projection_axis.lower()
    weighted_stress_extraction = "weighted_stress" in force_extraction_method_normalized
    checks = {
        "model_id_recorded": bool(model_id),
        "operating_point_id_recorded": bool(operating_point_id),
        "artifact_id_recorded": expected_artifact is None or bool(artifact_id),
        "result_set_id_recorded": expected_result_set is None or bool(result_set_id),
        "parameter_set_artifact_id_recorded": not parameter_set_artifact_required
        or bool(parameter_set_artifact_id),
        "parameter_set_digest_recorded": not parameter_set_digest_required
        or bool(parameter_set_digest),
        "parameter_set_path_recorded": not parameter_set_path_required
        or bool(parameter_set_path),
        "model_input_artifact_id_recorded": not model_input_artifact_required
        or bool(model_input_artifact_id),
        "model_input_digest_recorded": not model_input_digest_required
        or bool(model_input_digest),
        "model_input_path_recorded": not model_input_path_required
        or bool(model_input_path),
        "solution_artifact_id_recorded": expected_solution_artifact is None
        or bool(solution_artifact_id),
        "block_label_artifact_id_recorded": expected_block_label_artifact is None
        or bool(block_label_artifact_id),
        "expected_model_id_matches": expected_model is None or model_id == expected_model,
        "expected_operating_point_id_matches": expected_op is None or operating_point_id == expected_op,
        "expected_artifact_id_matches": expected_artifact is None or artifact_id == expected_artifact,
        "expected_result_set_id_matches": expected_result_set is None
        or result_set_id == expected_result_set,
        "expected_parameter_set_artifact_id_matches": (
            expected_parameter_set_artifact is None
            or parameter_set_artifact_id == expected_parameter_set_artifact
        ),
        "expected_parameter_set_digest_matches": (
            expected_parameter_set_digest_text is None
            or parameter_set_digest == expected_parameter_set_digest_text
        ),
        "expected_parameter_set_path_matches": (
            expected_parameter_set_path_text is None
            or parameter_set_path == expected_parameter_set_path_text
        ),
        "expected_model_input_artifact_id_matches": (
            expected_model_input_artifact is None
            or model_input_artifact_id == expected_model_input_artifact
        ),
        "expected_model_input_digest_matches": (
            expected_model_input_digest_text is None
            or model_input_digest == expected_model_input_digest_text
        ),
        "expected_model_input_path_matches": (
            expected_model_input_path_text is None
            or model_input_path == expected_model_input_path_text
        ),
        "expected_solution_artifact_id_matches": expected_solution_artifact is None
        or solution_artifact_id == expected_solution_artifact,
        "expected_block_label_artifact_id_matches": expected_block_label_artifact is None
        or block_label_artifact_id == expected_block_label_artifact,
        "solution_loaded_recorded": not solution_loaded_required or solution_loaded_recorded,
        "solution_loaded_before_postprocess": not solution_loaded_required or solution_loaded,
        "source_group_id_recorded": not source_group_required or bool(source_group_id),
        "target_group_id_recorded": not target_group_required or bool(target_group_id),
        "expected_source_group_id_matches": expected_source_group is None
        or source_group_id == expected_source_group,
        "expected_target_group_id_matches": expected_target_group is None
        or target_group_id == expected_target_group,
        "source_center_xy_recorded_when_expected": expected_source_center is None
        or source_center_xy_m is not None,
        "target_center_xy_recorded_when_expected": expected_target_center is None
        or target_center_xy_m is not None,
        "expected_source_center_xy_matches": _center_matches(source_center_xy_m, expected_source_center),
        "expected_target_center_xy_matches": _center_matches(target_center_xy_m, expected_target_center),
        "wire_center_separation_matches_separation_m": (
            center_distance_error_m is None
            or center_distance_error_m <= center_tolerance_m
        ),
        "source_region_recorded": expected_source_region_text is None or bool(source_region),
        "target_region_recorded": expected_target_region_text is None or bool(target_region),
        "source_material_recorded": expected_source_material_text is None or bool(source_material),
        "target_material_recorded": expected_target_material_text is None or bool(target_material),
        "expected_source_region_matches": expected_source_region_text is None
        or source_region == expected_source_region_text,
        "expected_target_region_matches": expected_target_region_text is None
        or target_region == expected_target_region_text,
        "expected_source_material_matches": expected_source_material_text is None
        or source_material == expected_source_material_text,
        "expected_target_material_matches": expected_target_material_text is None
        or target_material == expected_target_material_text,
        "postprocess_trace_id_recorded": not trace_id_required
        or bool(postprocess_trace_id),
        "expected_postprocess_trace_id_matches": expected_trace_id is None
        or postprocess_trace_id == expected_trace_id,
        "postprocess_command_digest_recorded": not command_digest_required
        or bool(postprocess_command_digest),
        "expected_postprocess_command_digest_matches": expected_command_digest is None
        or postprocess_command_digest == expected_command_digest,
        "postprocess_output_artifact_id_recorded": not output_artifact_required
        or bool(postprocess_output_artifact_id),
        "expected_postprocess_output_artifact_id_matches": expected_output_artifact is None
        or postprocess_output_artifact_id == expected_output_artifact,
        "postprocess_output_digest_recorded": not output_digest_required
        or bool(postprocess_output_digest),
        "expected_postprocess_output_digest_matches": expected_output_digest is None
        or postprocess_output_digest == expected_output_digest,
        "postprocess_output_schema_id_recorded": not output_schema_required
        or bool(postprocess_output_schema_id),
        "expected_postprocess_output_schema_id_matches": expected_output_schema_id is None
        or postprocess_output_schema_id == expected_output_schema_id,
        "postprocess_output_columns_recorded": not output_schema_required
        or bool(postprocess_output_columns),
        "expected_postprocess_output_columns_match": expected_output_columns is None
        or postprocess_output_columns == expected_output_columns,
        "postprocess_output_units_recorded": not output_schema_required
        or bool(postprocess_output_units),
        "expected_postprocess_output_units_match": expected_output_units is None
        or postprocess_output_units == expected_output_units,
        "postprocess_row_convention_schema_id_recorded": (
            not postprocess_row_convention_schema_required
            or bool(postprocess_row_convention_schema_id)
        ),
        "expected_postprocess_row_convention_schema_id_matches": (
            expected_postprocess_row_convention_schema is None
            or postprocess_row_convention_schema_id == expected_postprocess_row_convention_schema
        ),
        "postprocess_script_artifact_id_recorded": not script_artifact_required
        or bool(postprocess_script_artifact_id),
        "postprocess_script_digest_recorded": not script_digest_required
        or bool(postprocess_script_digest),
        "postprocess_script_path_recorded": not script_path_required
        or bool(postprocess_script_path),
        "expected_postprocess_script_artifact_id_matches": expected_script_artifact is None
        or postprocess_script_artifact_id == expected_script_artifact,
        "expected_postprocess_script_digest_matches": expected_script_digest is None
        or postprocess_script_digest == expected_script_digest,
        "expected_postprocess_script_path_matches": expected_script_path is None
        or postprocess_script_path == expected_script_path,
        "force_observable_id_recorded": expected_observable_id is None
        or bool(force_observable_id),
        "expected_force_observable_id_matches": expected_observable_id is None
        or force_observable_id == expected_observable_id,
        "force_observable_family_recorded": expected_observable_family is None
        or bool(force_observable_family),
        "expected_force_observable_family_matches": expected_observable_family is None
        or force_observable_family == expected_observable_family,
        "force_convention_schema_id_recorded": not force_convention_schema_required
        or bool(force_convention_schema_id),
        "expected_force_convention_schema_id_matches": (
            expected_force_convention_schema is None
            or force_convention_schema_id == expected_force_convention_schema
        ),
        "force_component_basis_schema_id_recorded": not force_component_basis_schema_required
        or bool(force_component_basis_schema_id),
        "expected_force_component_basis_schema_id_matches": (
            expected_force_component_basis_schema is None
            or force_component_basis_schema_id == expected_force_component_basis_schema
        ),
        "force_unit_basis_schema_id_recorded": not force_unit_basis_schema_required
        or bool(force_unit_basis_schema_id),
        "expected_force_unit_basis_schema_id_matches": (
            expected_force_unit_basis_schema is None
            or force_unit_basis_schema_id == expected_force_unit_basis_schema
        ),
        "objective_observable_id_recorded": expected_objective_id is None
        or bool(objective_observable_id),
        "expected_objective_observable_id_matches": expected_objective_id is None
        or objective_observable_id == expected_objective_id,
        "objective_observable_family_recorded": expected_objective_family is None
        or bool(objective_observable_family),
        "expected_objective_observable_family_matches": expected_objective_family is None
        or objective_observable_family == expected_objective_family,
        "force_component_frame_recorded_when_expected": expected_component_frame is None
        or bool(force_component_frame),
        "expected_force_component_frame_matches": expected_component_frame is None
        or force_component_frame == expected_component_frame,
        "radial_projection_axis_recorded_when_expected": expected_projection_axis is None
        or bool(radial_projection_axis),
        "expected_radial_projection_axis_matches": expected_projection_axis is None
        or radial_projection_axis == expected_projection_axis,
        "force_sign_convention_recorded_when_expected": expected_sign_convention is None
        or bool(force_sign_convention),
        "expected_force_sign_convention_matches": expected_sign_convention is None
        or force_sign_convention == expected_sign_convention,
        "force_extraction_method_recorded_when_expected": expected_extraction_method is None
        or bool(force_extraction_method_normalized),
        "expected_force_extraction_method_matches": expected_extraction_method is None
        or force_extraction_method_normalized == expected_extraction_method,
        "current_source_artifact_id_recorded_when_expected": not current_artifact_required
        or bool(current_source_artifact_id),
        "expected_current_source_artifact_id_matches": expected_current_artifact is None
        or current_source_artifact_id == expected_current_artifact,
        "current_definition_method_recorded_when_expected": not current_method_required
        or bool(current_definition_method_normalized),
        "expected_current_definition_method_matches": expected_current_method is None
        or current_definition_method_normalized == expected_current_method,
        "weighted_stress_extraction_uses_force_xy_integrals": (
            not weighted_stress_extraction or set(block_integral_types) == {18, 19}
        ),
        "postprocess_commands_recorded": not trace_required or bool(postprocess_commands),
        "postprocess_commands_clear_select_force_xy": (
            not trace_required or trace_has_clear_select_force_xy
        ),
        "postprocess_commands_mention_target_group": trace_mentions_target_group,
        "source_target_groups_distinct": not (source_group_id and target_group_id)
        or source_group_id != target_group_id,
        "selection_function_recorded": not target_group_required or bool(selection_function),
        "selection_mentions_target_group": selection_mentions_target_group,
        "selection_clear_before_groupselect": (
            not selection_clear_required or selection_clear_before_select
        ),
        "source_tool_recorded": bool(source_tool),
        "expected_source_tool_matches": expected_tool is None or source_tool == expected_tool,
        "source_function_recorded": bool(source_function),
        "force_component_frame_recorded": not femm_source or bool(force_component_frame),
        "radial_projection_axis_recorded": not femm_source
        or measured_radial is None
        or bool(radial_projection_axis),
        "radial_projection_axis_names_wire_pair": not femm_source
        or measured_radial is None
        or (
            ("wire1" in radial_axis_text and "wire2" in radial_axis_text)
            or "separation" in radial_axis_text
        ),
        "force_units_are_per_length": "n/m" in units_compact or "n_per_m" in units_compact,
        "force_unit_basis_is_per_length": force_unit_basis == "per_length",
        "depth_integrated_force_not_used_for_per_length_gate": force_unit_basis != "depth_integrated",
        "femm_planar_depth_recorded": not femm_source or planar_depth_m is not None,
        "femm_planar_depth_positive": planar_depth_m is None or planar_depth_m > 0.0,
        "problem_type_recorded": expected_problem is None or bool(problem_type_normalized),
        "expected_problem_type_matches": expected_problem is None
        or problem_type_normalized == expected_problem,
        "length_unit_recorded": expected_length is None or bool(length_unit_normalized),
        "expected_length_unit_matches": expected_length is None
        or length_unit_normalized == expected_length,
        "frequency_hz_recorded": expected_frequency is None or frequency_hz is not None,
        "expected_frequency_hz_matches": expected_frequency is None
        or (
            frequency_hz is not None
            and abs(frequency_hz - expected_frequency) <= frequency_tolerance
        ),
        "solver_precision_recorded": (
            expected_precision is None and max_precision is None
        )
        or solver_precision is not None,
        "expected_solver_precision_matches": expected_precision is None
        or (
            solver_precision is not None
            and abs(solver_precision - expected_precision) <= precision_tolerance
        ),
        "solver_precision_within_max": max_precision is None
        or (solver_precision is not None and solver_precision <= max_precision),
        "min_angle_deg_recorded": expected_min_angle is None or min_angle_deg is not None,
        "expected_min_angle_deg_matches": expected_min_angle is None
        or (
            min_angle_deg is not None
            and abs(min_angle_deg - expected_min_angle) <= min_angle_tolerance
        ),
        "min_angle_deg_positive": min_angle_deg is None or min_angle_deg > 0.0,
        "created_at_utc_recorded": not created_at_required or bool(created_at_utc),
        "created_at_utc_parseable": not created_at_utc or created_at_dt is not None,
        "expected_created_at_utc_matches": expected_created_at is None
        or created_at_utc == expected_created_at,
        "run_timestamp_utc_recorded": not execution_metadata_required
        or bool(run_timestamp_utc),
        "run_timestamp_utc_parseable": not run_timestamp_utc
        or run_timestamp_dt is not None,
        "expected_run_timestamp_utc_matches": expected_run_timestamp is None
        or run_timestamp_utc == expected_run_timestamp,
        "created_run_timestamp_skew_within_limit": max_created_run_skew is None
        or (
            created_run_skew_s is not None
            and created_run_skew_s <= max_created_run_skew
        ),
        "solver_version_recorded": not execution_metadata_required
        or bool(solver_version),
        "expected_solver_version_matches": expected_solver_version_text is None
        or solver_version == expected_solver_version_text,
        "radia_mcp_version_recorded": not execution_metadata_required
        or bool(radia_mcp_version),
        "expected_radia_mcp_version_matches": expected_radia_mcp_version_text is None
        or radia_mcp_version == expected_radia_mcp_version_text,
        "run_duration_s_recorded": not (
            execution_metadata_required or timing_breakdown_required
        )
        or run_duration_s is not None,
        "run_duration_s_positive": run_duration_s is None or run_duration_s > 0.0,
        "timing_breakdown_recorded": not timing_breakdown_required
        or bool(timing_breakdown_rows),
        "timing_breakdown_has_required_sections": not timing_breakdown_required
        or len(timing_breakdown_rows) >= timing_sections_min,
        "timing_breakdown_sections_named": not timing_breakdown_required
        or all(bool(row["name"]) for row in timing_breakdown_rows),
        "timing_breakdown_values_nonnegative": not timing_breakdown_required
        or all(value >= 0.0 for value in timing_durations),
        "timing_breakdown_top_sections_descending": not timing_breakdown_required
        or timing_top_sections_descending,
        "timing_breakdown_total_within_run_duration": (
            run_duration_s is None
            or timing_total_s is None
            or timing_total_s <= run_duration_s * (1.0 + tol) + 1.0e-12
        ),
        "force_measurement_present": has_force_measurement,
        "vector_force_matches_ampere": vector_ok,
        "radial_force_matches_ampere": radial_ok,
        "interaction_matches_current_sign": not interaction
        or interaction == analytic["interaction"],
    }
    if expected_integral_types is not None or block_integral_types:
        sorted_integral_types = sorted(block_integral_types)
        checks["block_integral_types_recorded"] = bool(sorted_integral_types)
        checks["block_integral_types_match_expected"] = (
            expected_integral_types is None
            or sorted_integral_types == expected_integral_types
        )
        checks["block_integral_types_are_force_xy"] = set(sorted_integral_types) == {18, 19}
        checks["block_integral_types_exclude_torque"] = 22 not in sorted_integral_types
    else:
        sorted_integral_types = []
    return {
        "policy": "parallel_wire_force_result_package_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "model_id": model_id,
        "operating_point_id": operating_point_id,
        "artifact_id": artifact_id,
        "result_set_id": result_set_id,
        "parameter_set_artifact_id": parameter_set_artifact_id,
        "parameter_set_digest": parameter_set_digest,
        "parameter_set_path": parameter_set_path,
        "model_input_artifact_id": model_input_artifact_id,
        "model_input_digest": model_input_digest,
        "model_input_path": model_input_path,
        "solution_artifact_id": solution_artifact_id,
        "block_label_artifact_id": block_label_artifact_id,
        "source_group_id": source_group_id,
        "target_group_id": target_group_id,
        "source_center_xy_m": list(source_center_xy_m) if source_center_xy_m is not None else None,
        "target_center_xy_m": list(target_center_xy_m) if target_center_xy_m is not None else None,
        "expected_source_center_xy_m": list(expected_source_center) if expected_source_center is not None else None,
        "expected_target_center_xy_m": list(expected_target_center) if expected_target_center is not None else None,
        "center_distance_m": center_distance_m,
        "center_distance_error_m": center_distance_error_m,
        "source_region": source_region,
        "target_region": target_region,
        "source_material": source_material,
        "target_material": target_material,
        "postprocess_trace_id": postprocess_trace_id,
        "postprocess_command_digest": postprocess_command_digest,
        "postprocess_output_artifact_id": postprocess_output_artifact_id,
        "postprocess_output_digest": postprocess_output_digest,
        "postprocess_output_path": postprocess_output_path,
        "postprocess_output_schema_id": postprocess_output_schema_id,
        "postprocess_output_columns": postprocess_output_columns,
        "postprocess_output_units": postprocess_output_units,
        "postprocess_row_convention_schema_id": postprocess_row_convention_schema_id,
        "postprocess_script_artifact_id": postprocess_script_artifact_id,
        "postprocess_script_digest": postprocess_script_digest,
        "postprocess_script_path": postprocess_script_path,
        "force_observable_id": force_observable_id,
        "force_observable_family": force_observable_family,
        "force_convention_schema_id": force_convention_schema_id,
        "force_component_basis_schema_id": force_component_basis_schema_id,
        "force_unit_basis_schema_id": force_unit_basis_schema_id,
        "objective_observable_id": objective_observable_id,
        "objective_observable_family": objective_observable_family,
        "postprocess_commands": postprocess_commands,
        "expected_artifact_id": expected_artifact,
        "expected_result_set_id": expected_result_set,
        "expected_parameter_set_artifact_id": expected_parameter_set_artifact,
        "expected_parameter_set_digest": expected_parameter_set_digest_text,
        "expected_parameter_set_path": expected_parameter_set_path_text,
        "expected_model_input_artifact_id": expected_model_input_artifact,
        "expected_model_input_digest": expected_model_input_digest_text,
        "expected_model_input_path": expected_model_input_path_text,
        "expected_solution_artifact_id": expected_solution_artifact,
        "expected_block_label_artifact_id": expected_block_label_artifact,
        "solution_loaded": solution_loaded if solution_loaded_recorded else None,
        "solution_loaded_required": solution_loaded_required,
        "expected_source_tool": expected_tool,
        "expected_source_group_id": expected_source_group,
        "expected_target_group_id": expected_target_group,
        "expected_source_region": expected_source_region_text,
        "expected_target_region": expected_target_region_text,
        "expected_source_material": expected_source_material_text,
        "expected_target_material": expected_target_material_text,
        "expected_postprocess_trace_id": expected_trace_id,
        "expected_postprocess_command_digest": expected_command_digest,
        "expected_postprocess_output_artifact_id": expected_output_artifact,
        "expected_postprocess_output_digest": expected_output_digest,
        "expected_postprocess_output_schema_id": expected_output_schema_id,
        "expected_postprocess_output_columns": expected_output_columns,
        "expected_postprocess_output_units": expected_output_units,
        "expected_postprocess_row_convention_schema_id": expected_postprocess_row_convention_schema,
        "expected_postprocess_script_artifact_id": expected_script_artifact,
        "expected_postprocess_script_digest": expected_script_digest,
        "expected_postprocess_script_path": expected_script_path,
        "expected_force_observable_id": expected_observable_id,
        "expected_force_observable_family": expected_observable_family,
        "expected_force_convention_schema_id": expected_force_convention_schema,
        "expected_force_component_basis_schema_id": expected_force_component_basis_schema,
        "expected_force_unit_basis_schema_id": expected_force_unit_basis_schema,
        "expected_objective_observable_id": expected_objective_id,
        "expected_objective_observable_family": expected_objective_family,
        "expected_force_component_frame": expected_component_frame,
        "expected_radial_projection_axis": expected_projection_axis,
        "expected_force_sign_convention": expected_sign_convention,
        "force_extraction_method": force_extraction_method_normalized or None,
        "raw_force_extraction_method": force_extraction_method,
        "expected_force_extraction_method": expected_extraction_method,
        "current_source_artifact_id": current_source_artifact_id,
        "current_definition_method": current_definition_method_normalized or None,
        "raw_current_definition_method": current_definition_method,
        "expected_current_source_artifact_id": expected_current_artifact,
        "expected_current_definition_method": expected_current_method,
        "selection_function": selection_function,
        "force_component_frame": force_component_frame,
        "radial_projection_axis": radial_projection_axis,
        "force_sign_convention": force_sign_convention,
        "block_integral_types": sorted_integral_types,
        "expected_block_integral_types": expected_integral_types,
        "selection_clear_required": selection_clear_required,
        "postprocess_command_trace_required": trace_required,
        "postprocess_output_artifact_required": output_artifact_required,
        "postprocess_output_schema_required": output_schema_required,
        "postprocess_row_convention_schema_required": postprocess_row_convention_schema_required,
        "postprocess_script_artifact_required": script_artifact_required,
        "force_convention_schema_required": force_convention_schema_required,
        "force_component_basis_schema_required": force_component_basis_schema_required,
        "force_unit_basis_schema_required": force_unit_basis_schema_required,
        "parameter_set_artifact_required": parameter_set_artifact_required,
        "model_input_artifact_required": model_input_artifact_required,
        "source_tool": source_tool,
        "source_function": source_function,
        "force_units": units,
        "force_unit_basis": force_unit_basis,
        "problem_depth_m": planar_depth_m,
        "problem_type": problem_type_normalized,
        "raw_problem_type": problem_type,
        "expected_problem_type": expected_problem,
        "length_unit": length_unit_normalized,
        "raw_length_unit": length_unit,
        "expected_length_unit": expected_length,
        "frequency_hz": frequency_hz,
        "expected_frequency_hz": expected_frequency,
        "solver_precision": solver_precision,
        "expected_solver_precision": expected_precision,
        "max_solver_precision": max_precision,
        "min_angle_deg": min_angle_deg,
        "expected_min_angle_deg": expected_min_angle,
        "created_at_utc": created_at_utc,
        "expected_created_at_utc": expected_created_at,
        "run_timestamp_utc": run_timestamp_utc,
        "expected_run_timestamp_utc": expected_run_timestamp,
        "created_run_timestamp_skew_s": created_run_skew_s,
        "max_created_run_skew_s": max_created_run_skew,
        "solver_version": solver_version,
        "expected_solver_version": expected_solver_version_text,
        "radia_mcp_version": radia_mcp_version,
        "expected_radia_mcp_version": expected_radia_mcp_version_text,
        "run_duration_s": run_duration_s,
        "timing_breakdown_s": timing_breakdown_rows,
        "timing_total_s": timing_total_s,
        "min_timing_sections": timing_sections_min,
        "execution_metadata_required": execution_metadata_required,
        "timing_breakdown_required": timing_breakdown_required,
        "analytic": analytic,
        "expected_force_on_wire2_N_per_m": expected_vector,
        "expected_radial_force_on_wire2_N_per_m": expected_radial,
        "vector_force_abs_error_N_per_m": vector_error,
        "radial_force_abs_error_N_per_m": radial_error,
        "rtol": tol,
        "checks": checks,
        "notes": [
            "Use this before comparing FEMM block-integral force rows with radia-ngsolve force rows.",
            "The sign convention is radial separation increasing away from wire 1; like currents give negative radial force on wire 2.",
            "Store force per unit length for 2D planar wire benchmarks, not total force unless depth is explicit.",
            "For FEMM planar rows, archive the mi_probdef problem depth but keep this gate on the N/m comparison basis.",
            "Bind the force unit-basis schema separately from the force component basis so per-length N/m rows, stored problem depth, and depth-integrated total-force rows cannot be silently interchanged.",
            "For FEMM planar magnetic force rows, keep mo_blockintegral(18/19) separate from torque mo_blockintegral(22).",
            "Record the force extraction method, e.g. weighted_stress_block_integral_xy, so N/m force rows are not confused with Lorentz, contour, gap-harmonic, or torque-derived observables.",
            "Bind mo_blockintegral(18/19) rows to a named force observable id and family so weighted-stress block-force evidence is not confused with Lorentz, contour, or torque postprocessing.",
            "Bind recovered rows to artifact_id and result_set_id so stale tables from another run cannot pass by value alone.",
            "Bind force rows to the current-source artifact and current-definition method so correct force values are not joined to a stale current table or RMS/peak convention.",
            "Bind optimization/notebook force rows to the parameter-set artifact id/digest/path and objective observable id/family so a correct force value is not reused for a stale design point or scalar objective.",
            "Bind FEMM solver rows to the input .fem model artifact id/digest/path so a copied .ans or force table cannot be joined to a stale geometry/source definition.",
            "For FEMM solver-ready rows, bind the loaded .ans solution artifact and require mi_loadsolution()/postprocessor-loaded state before block integrals.",
            "Archive mi_probdef/mo_getprobleminfo problem_type, length_unit, and frequency_hz so planar-static N/m rows are not mixed with axisymmetric, mm-scaled, or AC runs.",
            "Archive mi_probdef precision and Triangle minangle so force rows from coarse solver or mesh settings do not pass by value alone.",
            "Archive the FEMM block-label/source-contract artifact plus source/target region and material names so a correct force value from the right numeric group cannot be joined to stale block labels.",
            "Archive the vector component frame and radial projection axis before comparing scalar radial forces.",
            "When promoting FEMM force rows, match the expected component frame and projection axis explicitly so local/global or wire1/wire2 sign-convention drift cannot pass by value alone.",
            "Bind the force component-basis schema separately from the physics convention so Fx/Fy ordering, global-vs-local basis, and radial projection basis cannot silently drift while values remain plausible.",
            "For FEMM rows promoted to solver-ready evidence, require mo_clearblock() before mo_groupselectblock(<target>) so stale block selections cannot leak into the force integral.",
            "For repeatable FEMM postprocessing, archive the command trace id, digest, and replay command sequence that cleared selection, selected the target group, and called mo_blockintegral(18/19).",
            "Bind the postprocess command trace to the CSV/JSON output artifact id and digest so a correct replay script cannot be joined to a stale exported result table.",
            "Bind the postprocess output schema id, columns, and column units so a correct output file cannot be read through a stale force-table layout.",
            "For notebook/crossval reuse, archive run_timestamp_utc, solver_version, radia_mcp_version, run_duration_s, and the top timing_breakdown_s sections in the result JSON.",
            "When created_at_utc is available, keep it close to run_timestamp_utc so a copied result artifact is not mistaken for a freshly executed force row.",
        ],
    }


def magnetic_field_probe_result_package_gate(
    row,
    *,
    expected_model_id=None,
    expected_operating_point_id=None,
    expected_artifact_id=None,
    expected_solution_artifact_id=None,
    expected_solution_digest=None,
    expected_solution_path=None,
    expected_source_tool=None,
    expected_probe_id=None,
    expected_probe_point_xy_m=None,
    expected_problem_length_unit=None,
    expected_probe_point_input_unit=None,
    expected_coordinate_scale_to_m=None,
    expected_field_component_frame=None,
    expected_field_units=None,
    expected_field_probe_method=None,
    expected_postprocess_trace_id=None,
    expected_postprocess_command_digest=None,
    expected_probe_output_artifact_id=None,
    expected_probe_output_digest=None,
    require_solution_artifact=False,
    require_solution_loaded=False,
    require_postprocess_command_trace=False,
    require_probe_output_artifact=False,
    require_probe_coordinate_scale=False,
    point_atol_m=1.0e-9,
):
    """Check a point magnetic-field probe result row before solver comparison.

    This is the small result-package gate for FEMM ``mo_getb(px, py)`` rows,
    NGSolve point samples, or notebook-exported B-field probes.  It does not
    judge the field value against physics; it prevents a plausible Bx/By value
    from travelling with a stale solution, wrong probe point, wrong component
    frame, or stale output table.
    """

    data = dict(row)
    point_tol = float(point_atol_m)
    if point_tol < 0.0:
        raise ValueError("point_atol_m must be >= 0")

    def _normalize_label(value):
        return str(value).strip().lower().replace("-", "_").replace(" ", "_")

    def _point(value):
        if value in (None, ""):
            return None
        if isinstance(value, dict):
            if "x" in value and "y" in value:
                return (float(value["x"]), float(value["y"]))
            if "X" in value and "Y" in value:
                return (float(value["X"]), float(value["Y"]))
        parts = list(value) if not isinstance(value, str) else [
            item for item in re.split(r"[,;\s]+", value.strip()) if item
        ]
        if len(parts) == 3 and abs(float(parts[2])) <= point_tol:
            parts = parts[:2]
        if len(parts) != 2:
            raise ValueError("field probe point must contain x/y or x/y/z with z=0")
        return (float(parts[0]), float(parts[1]))

    def _field_vector():
        for key in ("B_T", "b_T", "field_T", "magnetic_flux_density_T"):
            if key in data and data[key] is not None:
                values = list(data[key])
                if len(values) == 3 and abs(float(values[2])) <= 1.0e-300:
                    values = values[:2]
                if len(values) != 2:
                    raise ValueError("field vector must contain two planar components")
                return (float(values[0]), float(values[1]))
        bx = data.get("Bx_T", data.get("bx_T", data.get("B_x_T")))
        by = data.get("By_T", data.get("by_T", data.get("B_y_T")))
        if bx is None or by is None:
            return None
        return (float(bx), float(by))

    model_id = str(data.get("model_id", "")).strip()
    operating_point_id = str(data.get("operating_point_id", "")).strip()
    artifact_id = str(data.get("artifact_id", data.get("case_artifact_id", ""))).strip()
    solution_artifact_id = str(
        data.get(
            "solution_artifact_id",
            data.get("ans_artifact_id", data.get("loaded_solution_artifact_id", "")),
        )
    ).strip()
    solution_digest = str(
        data.get(
            "solution_digest",
            data.get(
                "solution_sha256",
                data.get("ans_digest", data.get("ans_sha256", data.get("loaded_solution_digest", ""))),
            ),
        )
    ).strip()
    solution_path = str(
        data.get(
            "solution_path",
            data.get("ans_path", data.get("loaded_solution_path", "")),
        )
    ).strip()
    solution_loaded_raw = data.get(
        "solution_loaded",
        data.get("postprocessor_solution_loaded", data.get("mo_solution_loaded")),
    )
    solution_loaded_recorded = solution_loaded_raw not in (None, "")
    solution_loaded = _metadata_truthy(solution_loaded_raw)
    source_tool = str(data.get("source_tool", "")).strip()
    source_function = str(data.get("source_function", data.get("probe_function", ""))).strip()
    probe_id = str(
        data.get("field_probe_id", data.get("probe_id", data.get("observable_id", "")))
    ).strip()
    probe_point = _point(
        data.get(
            "probe_point_xy_m",
            data.get(
                "field_probe_point_xy_m",
                data.get("sample_point_xy_m", data.get("field_probe_point_xyz_m")),
            ),
        )
    )
    problem_length_unit = str(
        data.get("problem_length_unit", data.get("length_unit", data.get("femm_length_unit", "")))
    ).strip()
    probe_point_input = _point(
        data.get(
            "probe_point_input_xy",
            data.get("probe_point_xy_native", data.get("mo_getb_input_xy")),
        )
    )
    probe_point_input_unit = str(
        data.get(
            "probe_point_input_unit",
            data.get("probe_coordinate_unit", data.get("mo_getb_coordinate_unit", "")),
        )
    ).strip()
    coordinate_scale_raw = data.get(
        "coordinate_scale_to_m",
        data.get("length_unit_scale_to_m", data.get("probe_coordinate_scale_to_m")),
    )
    coordinate_scale_to_m = None if coordinate_scale_raw in (None, "") else float(coordinate_scale_raw)
    field_vector = _field_vector()
    field_units = str(
        data.get("field_units", data.get("B_units", data.get("unit", "")))
    ).strip()
    component_frame = str(
        data.get("field_component_frame", data.get("component_frame", ""))
    ).strip()
    probe_method = str(
        data.get("field_probe_method", data.get("probe_method", data.get("postprocess_method", "")))
    ).strip()
    probe_method_normalized = _normalize_label(probe_method) if probe_method else ""
    postprocess_trace_id = str(
        data.get(
            "postprocess_trace_id",
            data.get("postprocess_command_trace_id", data.get("command_trace_id", "")),
        )
    ).strip()
    postprocess_command_digest = str(
        data.get(
            "postprocess_command_digest",
            data.get("command_trace_sha256", data.get("probe_command_digest", "")),
        )
    ).strip()
    probe_output_artifact_id = str(
        data.get(
            "field_probe_output_artifact_id",
            data.get("probe_output_artifact_id", data.get("output_artifact_id", "")),
        )
    ).strip()
    probe_output_digest = str(
        data.get(
            "field_probe_output_digest",
            data.get("probe_output_sha256", data.get("output_digest", "")),
        )
    ).strip()
    probe_output_path = str(
        data.get(
            "field_probe_output_path",
            data.get("probe_output_path", data.get("output_path", "")),
        )
    ).strip()
    command_sequence_source = data.get(
        "postprocess_commands",
        data.get("postprocess_command_sequence", data.get("command_sequence")),
    )
    if command_sequence_source is None:
        postprocess_commands = []
    elif isinstance(command_sequence_source, str):
        postprocess_commands = [
            item.strip()
            for item in command_sequence_source.replace("\n", ";").split(";")
            if item.strip()
        ]
    else:
        postprocess_commands = [str(item).strip() for item in command_sequence_source if str(item).strip()]
    command_text = f"{source_function}; {'; '.join(postprocess_commands)}".lower()

    expected_model = None if expected_model_id is None else str(expected_model_id).strip()
    expected_op = (
        None if expected_operating_point_id is None else str(expected_operating_point_id).strip()
    )
    expected_artifact = (
        None if expected_artifact_id is None else str(expected_artifact_id).strip()
    )
    expected_solution = (
        None if expected_solution_artifact_id is None else str(expected_solution_artifact_id).strip()
    )
    expected_solution_hash = (
        None if expected_solution_digest is None else str(expected_solution_digest).strip()
    )
    expected_solution_file = (
        None if expected_solution_path is None else str(expected_solution_path).strip()
    )
    expected_tool = None if expected_source_tool is None else str(expected_source_tool).strip()
    expected_probe = None if expected_probe_id is None else str(expected_probe_id).strip()
    expected_point = _point(expected_probe_point_xy_m)
    expected_length_unit = (
        None if expected_problem_length_unit is None else str(expected_problem_length_unit).strip()
    )
    expected_input_unit = (
        None if expected_probe_point_input_unit is None else str(expected_probe_point_input_unit).strip()
    )
    expected_scale = (
        None if expected_coordinate_scale_to_m is None else float(expected_coordinate_scale_to_m)
    )
    expected_frame = (
        None if expected_field_component_frame is None else str(expected_field_component_frame).strip()
    )
    expected_units = None if expected_field_units is None else str(expected_field_units).strip()
    expected_method = (
        None if expected_field_probe_method is None else _normalize_label(expected_field_probe_method)
    )
    expected_trace = (
        None if expected_postprocess_trace_id is None else str(expected_postprocess_trace_id).strip()
    )
    expected_command_digest = (
        None
        if expected_postprocess_command_digest is None
        else str(expected_postprocess_command_digest).strip()
    )
    expected_output_artifact = (
        None
        if expected_probe_output_artifact_id is None
        else str(expected_probe_output_artifact_id).strip()
    )
    expected_output_digest = (
        None if expected_probe_output_digest is None else str(expected_probe_output_digest).strip()
    )

    field_finite = (
        field_vector is not None
        and all(math.isfinite(value) for value in field_vector)
    )
    point_matches = True
    if expected_point is not None:
        point_matches = (
            probe_point is not None
            and all(abs(actual - expected) <= point_tol for actual, expected in zip(probe_point, expected_point))
        )
    trace_required = bool(require_postprocess_command_trace)
    output_required = bool(require_probe_output_artifact or expected_output_artifact or expected_output_digest)
    solution_artifact_required = bool(
        require_solution_artifact
        or expected_solution is not None
        or expected_solution_hash is not None
        or expected_solution_file is not None
    )
    solution_required = bool(require_solution_loaded)
    coordinate_scale_required = bool(
        require_probe_coordinate_scale
        or expected_length_unit is not None
        or expected_input_unit is not None
        or expected_scale is not None
    )
    scale_matches = True
    if expected_scale is not None:
        scale_matches = (
            coordinate_scale_to_m is not None
            and math.isclose(coordinate_scale_to_m, expected_scale, rel_tol=1.0e-12, abs_tol=point_tol)
        )
    input_scale_matches = True
    if coordinate_scale_required:
        input_scale_matches = (
            probe_point is not None
            and probe_point_input is not None
            and coordinate_scale_to_m is not None
            and all(
                abs(actual * coordinate_scale_to_m - expected) <= point_tol
                for actual, expected in zip(probe_point_input, probe_point)
            )
        )
    checks = {
        "model_id_recorded": expected_model is None or bool(model_id),
        "expected_model_id_matches": expected_model is None or model_id == expected_model,
        "operating_point_id_recorded": expected_op is None or bool(operating_point_id),
        "expected_operating_point_id_matches": expected_op is None or operating_point_id == expected_op,
        "artifact_id_recorded": expected_artifact is None or bool(artifact_id),
        "expected_artifact_id_matches": expected_artifact is None or artifact_id == expected_artifact,
        "solution_artifact_id_recorded": expected_solution is None or bool(solution_artifact_id),
        "solution_artifact_id_recorded_when_required": (
            not solution_artifact_required or bool(solution_artifact_id)
        ),
        "expected_solution_artifact_id_matches": expected_solution is None or solution_artifact_id == expected_solution,
        "solution_digest_recorded_when_required": (
            not solution_artifact_required or bool(solution_digest)
        ),
        "expected_solution_digest_matches": (
            expected_solution_hash is None or solution_digest == expected_solution_hash
        ),
        "solution_path_recorded_when_required": (
            not solution_artifact_required or bool(solution_path)
        ),
        "expected_solution_path_matches": (
            expected_solution_file is None or solution_path == expected_solution_file
        ),
        "solution_loaded_recorded_when_required": not solution_required or solution_loaded_recorded,
        "solution_loaded_true_when_required": not solution_required or solution_loaded,
        "source_tool_recorded": bool(source_tool),
        "expected_source_tool_matches": expected_tool is None or source_tool == expected_tool,
        "source_function_recorded": bool(source_function),
        "probe_id_recorded": expected_probe is None or bool(probe_id),
        "expected_probe_id_matches": expected_probe is None or probe_id == expected_probe,
        "probe_point_recorded": expected_point is None or probe_point is not None,
        "expected_probe_point_xy_matches": point_matches,
        "problem_length_unit_recorded": not coordinate_scale_required or bool(problem_length_unit),
        "expected_problem_length_unit_matches": (
            expected_length_unit is None
            or problem_length_unit.lower() == expected_length_unit.lower()
        ),
        "probe_point_input_xy_recorded": not coordinate_scale_required or probe_point_input is not None,
        "probe_point_input_unit_recorded": not coordinate_scale_required or bool(probe_point_input_unit),
        "expected_probe_point_input_unit_matches": (
            expected_input_unit is None
            or probe_point_input_unit.lower() == expected_input_unit.lower()
        ),
        "coordinate_scale_to_m_recorded": not coordinate_scale_required or coordinate_scale_to_m is not None,
        "expected_coordinate_scale_to_m_matches": scale_matches,
        "probe_point_input_scale_matches_probe_point_xy_m": input_scale_matches,
        "field_vector_present": field_vector is not None,
        "field_vector_finite": field_finite,
        "field_units_recorded_when_expected": expected_units is None or bool(field_units),
        "expected_field_units_matches": expected_units is None or field_units == expected_units,
        "component_frame_recorded_when_expected": expected_frame is None or bool(component_frame),
        "expected_field_component_frame_matches": expected_frame is None or component_frame == expected_frame,
        "probe_method_recorded_when_expected": expected_method is None or bool(probe_method_normalized),
        "expected_field_probe_method_matches": expected_method is None or probe_method_normalized == expected_method,
        "postprocess_trace_id_recorded": not trace_required or bool(postprocess_trace_id),
        "expected_postprocess_trace_id_matches": expected_trace is None or postprocess_trace_id == expected_trace,
        "postprocess_command_digest_recorded": not trace_required or bool(postprocess_command_digest),
        "expected_postprocess_command_digest_matches": (
            expected_command_digest is None or postprocess_command_digest == expected_command_digest
        ),
        "postprocess_commands_recorded": not trace_required or bool(postprocess_commands),
        "postprocess_commands_include_mo_getb": not trace_required or "mo_getb" in command_text,
        "probe_output_artifact_id_recorded": not output_required or bool(probe_output_artifact_id),
        "expected_probe_output_artifact_id_matches": (
            expected_output_artifact is None or probe_output_artifact_id == expected_output_artifact
        ),
        "probe_output_digest_recorded": not output_required or bool(probe_output_digest),
        "expected_probe_output_digest_matches": (
            expected_output_digest is None or probe_output_digest == expected_output_digest
        ),
        "probe_output_path_recorded": not output_required or bool(probe_output_path),
    }
    return {
        "policy": "magnetic_field_probe_result_package_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "model_id": model_id,
        "operating_point_id": operating_point_id,
        "artifact_id": artifact_id,
        "solution_artifact_id": solution_artifact_id,
        "solution_digest": solution_digest or None,
        "solution_path": solution_path or None,
        "solution_loaded": solution_loaded if solution_loaded_recorded else None,
        "source_tool": source_tool,
        "source_function": source_function,
        "field_probe_id": probe_id,
        "probe_point_xy_m": list(probe_point) if probe_point is not None else None,
        "expected_probe_point_xy_m": list(expected_point) if expected_point is not None else None,
        "problem_length_unit": problem_length_unit,
        "probe_point_input_xy": list(probe_point_input) if probe_point_input is not None else None,
        "probe_point_input_unit": probe_point_input_unit,
        "coordinate_scale_to_m": coordinate_scale_to_m,
        "B_T": list(field_vector) if field_vector is not None else None,
        "field_units": field_units,
        "field_component_frame": component_frame,
        "field_probe_method": probe_method_normalized or None,
        "raw_field_probe_method": probe_method,
        "postprocess_trace_id": postprocess_trace_id,
        "postprocess_command_digest": postprocess_command_digest,
        "postprocess_commands": postprocess_commands,
        "probe_output_artifact_id": probe_output_artifact_id,
        "probe_output_digest": probe_output_digest,
        "probe_output_path": probe_output_path,
        "expected_model_id": expected_model,
        "expected_operating_point_id": expected_op,
        "expected_artifact_id": expected_artifact,
        "expected_solution_artifact_id": expected_solution,
        "expected_solution_digest": expected_solution_hash,
        "expected_solution_path": expected_solution_file,
        "expected_source_tool": expected_tool,
        "expected_probe_id": expected_probe,
        "expected_problem_length_unit": expected_length_unit,
        "expected_probe_point_input_unit": expected_input_unit,
        "expected_coordinate_scale_to_m": expected_scale,
        "expected_field_component_frame": expected_frame,
        "expected_field_units": expected_units,
        "expected_field_probe_method": expected_method,
        "expected_postprocess_trace_id": expected_trace,
        "expected_postprocess_command_digest": expected_command_digest,
        "expected_probe_output_artifact_id": expected_output_artifact,
        "expected_probe_output_digest": expected_output_digest,
        "solution_loaded_required": solution_required,
        "solution_artifact_required": solution_artifact_required,
        "postprocess_command_trace_required": trace_required,
        "probe_output_artifact_required": output_required,
        "probe_coordinate_scale_required": coordinate_scale_required,
        "point_atol_m": point_tol,
        "checks": checks,
        "notes": [
            "Use this before comparing FEMM mo_getb, NGSolve, or notebook B-field probe rows.",
            "Bind point B samples to the loaded solution artifact and probe point; value agreement alone is not reusable evidence.",
            "For FEMM .ans reuse, bind solution artifact id, digest, and path before trusting point B rows.",
            "For FEMM, record mi_loadsolution state and a postprocess command trace containing mo_getb(px, py).",
            "When comparing across solvers, record FEMM problem length units, the raw mo_getb input coordinate, and the scale used to express probe_point_xy_m in meters.",
        ],
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

    from radia.force import air_gap_maxwell_pressure as pressure

    return pressure(B_T, permeability_H_per_m=mu)


def air_gap_holding_force(B_T, area_m2, faces=1, mu=MU0):
    """Uniform-gap holding force [N] from flux density and active pole area.

    ``faces`` is the number of active, equal pole faces/gaps contributing the
    same pressure.  Use ``faces=2`` for a symmetric two-pole yoke with two equal
    gaps; keep ``faces=1`` for a single plunger or one pole face.
    """

    from radia.force import air_gap_holding_force as holding_force

    return holding_force(
        B_T,
        area_m2,
        faces=faces,
        permeability_H_per_m=mu,
    )


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

    from radia.force import air_gap_shear_stress as shear_stress

    return shear_stress(
        B_radial_T,
        B_tangential_T,
        permeability_H_per_m=mu,
    )


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

    from radia.force import air_gap_shear_torque as shear_torque

    return shear_torque(
        B_radial_T,
        B_tangential_T,
        radius_m,
        axial_length_m=axial_length_m,
        angle_rad=angle_rad,
        permeability_H_per_m=mu,
    )


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
    """Compatibility adapter to the canonical :mod:`radia.force` kernel."""

    from radia.force import air_gap_shear_torque_from_angle_samples as integrate

    result = integrate(
        angles_rad,
        B_radial_T,
        B_tangential_T,
        radius_m,
        axial_length_m=axial_length_m,
        periodic=periodic,
        period_rad=period_rad,
        permeability_H_per_m=mu,
    )
    result["mu"] = result["permeability_H_per_m"]
    return result


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

    from radia.force import coenergy_torque_from_angle_samples as differentiate

    torque_values = differentiate(
        angles,
        values,
        periodic=periodic,
        period_rad=period,
    )

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
        torque = float(torque_values[i])
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


def coenergy_torque_table_consistency_summary(
    angles_rad,
    coenergy_J,
    torque_Nm,
    periodic=False,
    period_rad=2.0 * math.pi,
    torque_abs_tolerance_Nm=0.0,
    torque_rel_tolerance=1.0e-6,
    comparison_stencils=None,
):
    """Compare a torque-angle table with ``d(coenergy)/dtheta``.

    The summary preserves the finite-difference torque inferred from coenergy,
    the supplied torque table, selected-stencil errors, and basic angle-step
    diagnostics.  For nonperiodic sweeps the default reference check uses only
    central rows; for periodic sweeps it uses all ``central_periodic`` rows.

    A nonzero mean torque over a full angle sweep implies that the coenergy
    table contains a work term such as ``T_mean * theta`` and should be treated
    as nonperiodic.  A purely periodic coenergy table has zero mean derivative
    over its period.
    """

    angles = [float(value) for value in angles_rad]
    values = [float(value) for value in coenergy_J]
    references = [float(value) for value in torque_Nm]
    if len(angles) != len(references):
        raise ValueError("torque_Nm must have the same length as angles_rad")
    rows = [
        dict(row)
        for row in coenergy_torque_from_angle_samples(
            angles,
            values,
            periodic=periodic,
            period_rad=period_rad,
        )
    ]
    if comparison_stencils is None:
        stencil_filter = {"central_periodic"} if periodic else {"central"}
    elif isinstance(comparison_stencils, str):
        token = comparison_stencils.lower().strip()
        stencil_filter = None if token in ("all", "*") else {token}
    else:
        stencil_filter = {str(value).lower().strip() for value in comparison_stencils}
        if "all" in stencil_filter or "*" in stencil_filter:
            stencil_filter = None
    if stencil_filter == set():
        raise ValueError("comparison_stencils must not be empty")

    abs_tol = float(torque_abs_tolerance_Nm)
    rel_tol = float(torque_rel_tolerance)
    if abs_tol < 0.0:
        raise ValueError("torque_abs_tolerance_Nm must be >= 0")
    if rel_tol < 0.0:
        raise ValueError("torque_rel_tolerance must be >= 0")

    selected_rows = []
    for row, reference in zip(rows, references):
        error = row["torque_Nm"] - reference
        abs_error = abs(error)
        rel_error = abs_error / max(abs(reference), 1.0e-300)
        selected = stencil_filter is None or row["stencil"].lower() in stencil_filter
        row["reference_torque_Nm"] = reference
        row["torque_error_Nm"] = error
        row["torque_abs_error_Nm"] = abs_error
        row["torque_rel_error"] = rel_error
        row["selected_for_reference_check"] = selected
        if selected:
            selected_rows.append(row)
    if not selected_rows:
        raise ValueError("comparison_stencils selected no rows")

    step_rows = [angles[i + 1] - angles[i] for i in range(len(angles) - 1)]
    reference_pass = all(
        row["torque_abs_error_Nm"]
        <= abs_tol + rel_tol * max(abs(row["reference_torque_Nm"]), 1.0e-300)
        for row in selected_rows
    )
    inferred = [row["torque_Nm"] for row in rows]
    errors = [row["torque_error_Nm"] for row in selected_rows]
    abs_errors = [abs(value) for value in errors]
    rel_errors = [row["torque_rel_error"] for row in selected_rows]
    trapezoid_work = sum(
        0.5 * (references[i] + references[i + 1]) * step_rows[i]
        for i in range(len(step_rows))
    )
    coenergy_delta = values[-1] - values[0]
    period = float(period_rad)
    angle_span = angles[-1] - angles[0]
    periodic_gap = period - angle_span if periodic else None
    return {
        "policy": "coenergy_torque_angle_table_consistency",
        "n_samples": len(rows),
        "periodic": bool(periodic),
        "period_rad": period,
        "angle_min_rad": angles[0],
        "angle_max_rad": angles[-1],
        "angle_span_rad": angle_span,
        "angle_step_min_rad": min(step_rows) if step_rows else 0.0,
        "angle_step_max_rad": max(step_rows) if step_rows else 0.0,
        "periodic_gap_rad": periodic_gap,
        "coenergy_delta_J": coenergy_delta,
        "reference_torque_trapezoid_work_J": trapezoid_work,
        "reference_work_minus_coenergy_delta_J": trapezoid_work - coenergy_delta,
        "comparison_stencils": (
            "all" if stencil_filter is None else sorted(stencil_filter)
        ),
        "reference_checked_count": len(selected_rows),
        "torque_abs_tolerance_Nm": abs_tol,
        "torque_rel_tolerance": rel_tol,
        "reference_pass": reference_pass,
        "max_torque_abs_error_Nm": max(abs_errors),
        "max_torque_rel_error": max(rel_errors),
        "mean_torque_error_Nm": sum(errors) / len(errors),
        "rms_torque_error_Nm": math.sqrt(sum(value * value for value in errors) / len(errors)),
        "inferred_torque_min_Nm": min(inferred),
        "inferred_torque_max_Nm": max(inferred),
        "reference_torque_min_Nm": min(references),
        "reference_torque_max_Nm": max(references),
        "status": "ok" if reference_pass else "needs_attention",
        "ok_for_torque_table": reference_pass,
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

    from radia.force import virtual_work_force_from_displacement_samples as differentiate

    force_values = differentiate(
        positions,
        values,
        energy_kind=normalized_kind,
    )
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
        force = float(force_values[i])
        derivative = force / sign
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


def virtual_work_force_sweep_audit_summary(
    positions_m,
    energy_J,
    energy_kind="coenergy",
    reference_force_N=None,
    force_abs_tolerance_N=0.0,
    force_rel_tolerance=1.0e-6,
    comparison_stencils=("central",),
):
    """Audit a virtual-work force sweep against an optional reference table.

    This wraps :func:`virtual_work_force_from_displacement_samples` with the
    extra bookkeeping needed for validation-class sweeps: second differences
    of the energy table, force-gradient estimates, optional reference-force
    errors, and pass/fail tolerances.  By default only central-difference rows
    are used for the reference check, because endpoint one-sided derivatives
    are expected to be lower order.
    """

    positions = [float(value) for value in positions_m]
    values = [float(value) for value in energy_J]
    normalized_kind, sign, identity = _virtual_work_energy_sign(energy_kind)
    rows = [
        dict(row)
        for row in virtual_work_force_from_displacement_samples(
            positions,
            values,
            energy_kind=normalized_kind,
        )
    ]
    n = len(rows)
    for i, row in enumerate(rows):
        if 0 < i < n - 1:
            dx_left = positions[i] - positions[i - 1]
            dx_right = positions[i + 1] - positions[i]
            left_slope = (values[i] - values[i - 1]) / dx_left
            right_slope = (values[i + 1] - values[i]) / dx_right
            second = 2.0 * (right_slope - left_slope) / (positions[i + 1] - positions[i - 1])
            row["energy_second_derivative_J_per_m2"] = second
            row["force_gradient_N_per_m"] = sign * second
        else:
            row["energy_second_derivative_J_per_m2"] = None
            row["force_gradient_N_per_m"] = None

    if comparison_stencils is None:
        stencil_filter = None
    elif isinstance(comparison_stencils, str):
        token = comparison_stencils.lower().strip()
        stencil_filter = None if token in ("all", "*") else {token}
    else:
        stencil_filter = {str(value).lower().strip() for value in comparison_stencils}
        if "all" in stencil_filter or "*" in stencil_filter:
            stencil_filter = None
    if stencil_filter == set():
        raise ValueError("comparison_stencils must not be empty")

    abs_tol = float(force_abs_tolerance_N)
    rel_tol = float(force_rel_tolerance)
    if abs_tol < 0.0:
        raise ValueError("force_abs_tolerance_N must be >= 0")
    if rel_tol < 0.0:
        raise ValueError("force_rel_tolerance must be >= 0")

    reference_pass = None
    reference_checked_count = 0
    selected_error_rows = []
    if reference_force_N is not None:
        references = [float(value) for value in reference_force_N]
        if len(references) != n:
            raise ValueError("reference_force_N must have the same length as positions_m")
        for row, reference in zip(rows, references):
            error = row["force_N"] - reference
            abs_error = abs(error)
            rel_error = abs_error / max(abs(reference), 1.0e-300)
            row["reference_force_N"] = reference
            row["force_error_N"] = error
            row["force_abs_error_N"] = abs_error
            row["force_rel_error"] = rel_error
            selected = stencil_filter is None or row["stencil"].lower() in stencil_filter
            row["selected_for_reference_check"] = selected
            if selected:
                selected_error_rows.append(row)
        if not selected_error_rows:
            raise ValueError("comparison_stencils selected no rows")
        reference_checked_count = len(selected_error_rows)
        reference_pass = all(
            row["force_abs_error_N"]
            <= abs_tol + rel_tol * max(abs(row["reference_force_N"]), 1.0e-300)
            for row in selected_error_rows
        )
    else:
        for row in rows:
            row["selected_for_reference_check"] = False

    forces = [row["force_N"] for row in rows]
    gradients = [
        row["force_gradient_N_per_m"]
        for row in rows
        if row["force_gradient_N_per_m"] is not None
    ]
    selected_abs_errors = [row["force_abs_error_N"] for row in selected_error_rows]
    selected_rel_errors = [row["force_rel_error"] for row in selected_error_rows]
    ok = True if reference_pass is None else bool(reference_pass)
    return {
        "policy": "virtual_work_force_sweep_energy_table_audit",
        "n_samples": n,
        "energy_kind": normalized_kind,
        "energy_to_force_sign": sign,
        "virtual_work_identity": identity,
        "force_min_N": min(forces),
        "force_max_N": max(forces),
        "force_peak_abs_N": max(abs(value) for value in forces),
        "force_span_N": max(forces) - min(forces),
        "force_mean_N": sum(forces) / n,
        "max_abs_force_gradient_N_per_m": (
            max(abs(value) for value in gradients) if gradients else 0.0
        ),
        "reference_compared": reference_force_N is not None,
        "reference_checked_count": reference_checked_count,
        "comparison_stencils": (
            "all" if stencil_filter is None else sorted(stencil_filter)
        ),
        "force_abs_tolerance_N": abs_tol,
        "force_rel_tolerance": rel_tol,
        "reference_pass": reference_pass,
        "max_reference_force_abs_error_N": (
            max(selected_abs_errors) if selected_abs_errors else None
        ),
        "max_reference_force_rel_error": (
            max(selected_rel_errors) if selected_rel_errors else None
        ),
        "status": "ok" if ok else "needs_attention",
        "ok_for_force_sweep": ok,
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


def two_port_scattering_sweep_momentum_force_summary(
    frequency_Hz,
    s11_values,
    s21_values,
    power_incident_W=1.0,
    incident_direction=(1.0, 0.0, 0.0),
    transmitted_direction=None,
    speed=C0,
    passivity_tolerance=1.0e-12,
):
    """Summarize two-port scattering momentum force over a frequency sweep.

    The sweep assumes port-1 excitation with power-normalized S-parameters:
    ``R = |S11|^2`` and ``T = |S21|^2``.  Unlike the single-point helper, rows
    with ``R + T > 1`` are retained and flagged as passivity violations.  This
    is useful for measured, interpolated, or lightly noisy solver data where the
    diagnostic table is more useful than an immediate exception.
    """

    frequencies = [float(value) for value in frequency_Hz]
    s11 = [complex(value) for value in s11_values]
    s21 = [complex(value) for value in s21_values]
    if len(frequencies) != len(s11) or len(frequencies) != len(s21):
        raise ValueError("frequency_Hz, s11_values, and s21_values must have the same length")
    if not frequencies:
        raise ValueError("at least one frequency sample is required")
    power = float(power_incident_W)
    speed = float(speed)
    tolerance = float(passivity_tolerance)
    if power < 0.0:
        raise ValueError("power_incident_W must be >= 0")
    if speed <= 0.0:
        raise ValueError("speed must be > 0")
    if tolerance < 0.0:
        raise ValueError("passivity_tolerance must be >= 0")
    inc = _unit_vector(incident_direction, "incident_direction")
    if transmitted_direction is None:
        out = list(inc)
    else:
        out = _unit_vector(transmitted_direction, "transmitted_direction")
        if len(out) != len(inc):
            raise ValueError("incident_direction and transmitted_direction must have the same length")

    rows = []
    violation_rows = []
    momentum_scale = power / speed
    for idx, (frequency, gamma, tau) in enumerate(zip(frequencies, s11, s21)):
        if not math.isfinite(frequency) or frequency < 0.0:
            raise ValueError("frequency samples must be finite and >= 0")
        if not math.isfinite(gamma.real) or not math.isfinite(gamma.imag):
            raise ValueError("s11 values must be finite")
        if not math.isfinite(tau.real) or not math.isfinite(tau.imag):
            raise ValueError("s21 values must be finite")
        s11_mag = abs(gamma)
        s21_mag = abs(tau)
        reflectance = s11_mag * s11_mag
        transmittance = s21_mag * s21_mag
        outgoing_fraction = reflectance + transmittance
        absorptance = 1.0 - outgoing_fraction
        force = [
            momentum_scale * ((1.0 + reflectance) * inc[axis] - transmittance * out[axis])
            for axis in range(len(inc))
        ]
        axial_force = sum(value * axis for value, axis in zip(force, inc))
        row = {
            "index": idx,
            "frequency_Hz": frequency,
            "s11_real": gamma.real,
            "s11_imag": gamma.imag,
            "s11_magnitude": s11_mag,
            "s11_phase_rad": math.atan2(gamma.imag, gamma.real),
            "s11_phase_deg": math.degrees(math.atan2(gamma.imag, gamma.real)),
            "s21_real": tau.real,
            "s21_imag": tau.imag,
            "s21_magnitude": s21_mag,
            "s21_phase_rad": math.atan2(tau.imag, tau.real),
            "s21_phase_deg": math.degrees(math.atan2(tau.imag, tau.real)),
            "reflectance": reflectance,
            "transmittance": transmittance,
            "absorptance": absorptance,
            "outgoing_power_fraction": outgoing_fraction,
            "power_reflected_W": reflectance * power,
            "power_transmitted_W": transmittance * power,
            "power_absorbed_W": absorptance * power,
            "force_N": force,
            "force_magnitude_N": math.sqrt(sum(value * value for value in force)),
            "axial_force_along_incident_direction_N": axial_force,
            "passivity_excess_power_fraction": max(0.0, outgoing_fraction - 1.0),
            "passivity_ok": outgoing_fraction <= 1.0 + tolerance,
        }
        rows.append(row)
        if not row["passivity_ok"]:
            violation_rows.append(row)

    max_force_row = max(rows, key=lambda row: row["force_magnitude_N"])
    min_force_row = min(rows, key=lambda row: row["force_magnitude_N"])
    max_outgoing_row = max(rows, key=lambda row: row["outgoing_power_fraction"])
    mean_force = sum(row["force_magnitude_N"] for row in rows) / len(rows)
    mean_axial = sum(row["axial_force_along_incident_direction_N"] for row in rows) / len(rows)
    monotonic = all(
        frequencies[idx] < frequencies[idx + 1]
        for idx in range(len(frequencies) - 1)
    )

    return {
        "policy": "two_port_scattering_sweep_momentum_force_audit",
        "n_points": len(rows),
        "frequency_min_Hz": min(frequencies),
        "frequency_max_Hz": max(frequencies),
        "frequency_monotonic_increasing": monotonic,
        "power_incident_W": power,
        "speed_m_per_s": speed,
        "incident_direction": inc,
        "transmitted_direction": out,
        "passivity_tolerance": tolerance,
        "passivity_ok": not violation_rows,
        "passivity_violation_count": len(violation_rows),
        "max_outgoing_power_fraction": max_outgoing_row["outgoing_power_fraction"],
        "max_outgoing_power_frequency_Hz": max_outgoing_row["frequency_Hz"],
        "max_passivity_excess_power_fraction": max(
            row["passivity_excess_power_fraction"] for row in rows
        ),
        "mean_reflectance": sum(row["reflectance"] for row in rows) / len(rows),
        "mean_transmittance": sum(row["transmittance"] for row in rows) / len(rows),
        "mean_absorptance": sum(row["absorptance"] for row in rows) / len(rows),
        "mean_force_magnitude_N": mean_force,
        "mean_axial_force_N": mean_axial,
        "max_force_magnitude_N": max_force_row["force_magnitude_N"],
        "max_force_frequency_Hz": max_force_row["frequency_Hz"],
        "min_force_magnitude_N": min_force_row["force_magnitude_N"],
        "min_force_frequency_Hz": min_force_row["frequency_Hz"],
        "force_span_N": max_force_row["force_magnitude_N"] - min_force_row["force_magnitude_N"],
        "max_force_row": max_force_row,
        "min_force_row": min_force_row,
        "passivity_violation_rows": violation_rows,
        "status": "ok" if not violation_rows else "needs_attention",
        "rows": rows,
        "two_port_sweep_force_formula": "F=P/c*((1+|S11|^2)k_in-|S21|^2 k_out)",
    }


def two_port_sparameter_sweep_health_summary(
    frequency_Hz,
    s11_values,
    s21_values,
    s12_values=None,
    s22_values=None,
    power_incident_W=1.0,
    incident_direction=(1.0, 0.0, 0.0),
    transmitted_direction=None,
    speed=C0,
    passivity_tolerance=1.0e-12,
    reciprocity_tolerance=1.0e-6,
    return_symmetry_tolerance=None,
):
    """Audit a two-port S-parameter sweep for RF force post-processing.

    The base force/passivity table is computed from port-1 excitation using
    ``S11`` and ``S21``.  If ``S12`` is supplied, the summary also checks
    reciprocity as ``max |S21-S12|``.  If ``S22`` is supplied, return symmetry
    can be checked as ``max |S11-S22|``.
    """

    sweep = two_port_scattering_sweep_momentum_force_summary(
        frequency_Hz,
        s11_values,
        s21_values,
        power_incident_W=power_incident_W,
        incident_direction=incident_direction,
        transmitted_direction=transmitted_direction,
        speed=speed,
        passivity_tolerance=passivity_tolerance,
    )
    frequencies = [float(value) for value in frequency_Hz]
    s11 = [complex(value) for value in s11_values]
    s21 = [complex(value) for value in s21_values]
    reciprocity_tol = float(reciprocity_tolerance)
    if reciprocity_tol < 0.0:
        raise ValueError("reciprocity_tolerance must be >= 0")

    reciprocity_rows = []
    max_reciprocity_error = None
    max_reciprocity_frequency = None
    reciprocity_ok = None
    if s12_values is not None:
        s12 = [complex(value) for value in s12_values]
        if len(s12) != len(frequencies):
            raise ValueError("s12_values must have the same length as frequency_Hz")
        for frequency, forward, reverse in zip(frequencies, s21, s12):
            if not math.isfinite(reverse.real) or not math.isfinite(reverse.imag):
                raise ValueError("s12 values must be finite")
            error = abs(forward - reverse)
            reciprocity_rows.append({
                "frequency_Hz": frequency,
                "s21_real": forward.real,
                "s21_imag": forward.imag,
                "s12_real": reverse.real,
                "s12_imag": reverse.imag,
                "s21_s12_abs_error": error,
                "reciprocity_ok": error <= reciprocity_tol,
            })
        worst = max(reciprocity_rows, key=lambda row: row["s21_s12_abs_error"])
        max_reciprocity_error = worst["s21_s12_abs_error"]
        max_reciprocity_frequency = worst["frequency_Hz"]
        reciprocity_ok = max_reciprocity_error <= reciprocity_tol

    return_symmetry_tol = (
        None if return_symmetry_tolerance is None else float(return_symmetry_tolerance)
    )
    if return_symmetry_tol is not None and return_symmetry_tol < 0.0:
        raise ValueError("return_symmetry_tolerance must be >= 0")
    effective_return_symmetry_tol = return_symmetry_tol
    return_symmetry_rows = []
    max_return_symmetry_error = None
    max_return_symmetry_frequency = None
    return_symmetry_ok = None
    if s22_values is not None:
        s22 = [complex(value) for value in s22_values]
        if len(s22) != len(frequencies):
            raise ValueError("s22_values must have the same length as frequency_Hz")
        effective_return_symmetry_tol = (
            reciprocity_tol if return_symmetry_tol is None else return_symmetry_tol
        )
        for frequency, port1, port2 in zip(frequencies, s11, s22):
            if not math.isfinite(port2.real) or not math.isfinite(port2.imag):
                raise ValueError("s22 values must be finite")
            error = abs(port1 - port2)
            return_symmetry_rows.append({
                "frequency_Hz": frequency,
                "s11_real": port1.real,
                "s11_imag": port1.imag,
                "s22_real": port2.real,
                "s22_imag": port2.imag,
                "s11_s22_abs_error": error,
                "return_symmetry_ok": error <= effective_return_symmetry_tol,
            })
        worst = max(return_symmetry_rows, key=lambda row: row["s11_s22_abs_error"])
        max_return_symmetry_error = worst["s11_s22_abs_error"]
        max_return_symmetry_frequency = worst["frequency_Hz"]
        return_symmetry_ok = max_return_symmetry_error <= effective_return_symmetry_tol

    issues = []
    if not sweep["passivity_ok"]:
        issues.append("passivity violation in port-1 outgoing power")
    if reciprocity_ok is False:
        issues.append("S21/S12 reciprocity error exceeds tolerance")
    if return_symmetry_ok is False:
        issues.append("S11/S22 return symmetry error exceeds tolerance")

    return {
        "policy": "two_port_sparameter_sweep_health",
        "status": "ok" if not issues else "needs_attention",
        "issues": issues,
        "n_points": sweep["n_points"],
        "frequency_monotonic_increasing": sweep["frequency_monotonic_increasing"],
        "passivity_ok": sweep["passivity_ok"],
        "passivity_violation_count": sweep["passivity_violation_count"],
        "max_passivity_excess_power_fraction": sweep["max_passivity_excess_power_fraction"],
        "reciprocity_checked": s12_values is not None,
        "reciprocity_tolerance": reciprocity_tol,
        "reciprocity_ok": reciprocity_ok,
        "max_s21_s12_abs_error": max_reciprocity_error,
        "max_s21_s12_error_frequency_Hz": max_reciprocity_frequency,
        "return_symmetry_checked": s22_values is not None,
        "return_symmetry_tolerance": effective_return_symmetry_tol,
        "return_symmetry_ok": return_symmetry_ok,
        "max_s11_s22_abs_error": max_return_symmetry_error,
        "max_s11_s22_error_frequency_Hz": max_return_symmetry_frequency,
        "mean_reflectance": sweep["mean_reflectance"],
        "mean_transmittance": sweep["mean_transmittance"],
        "mean_absorptance": sweep["mean_absorptance"],
        "mean_force_magnitude_N": sweep["mean_force_magnitude_N"],
        "max_force_magnitude_N": sweep["max_force_magnitude_N"],
        "max_force_frequency_Hz": sweep["max_force_frequency_Hz"],
        "min_force_magnitude_N": sweep["min_force_magnitude_N"],
        "min_force_frequency_Hz": sweep["min_force_frequency_Hz"],
        "reciprocity_rows": reciprocity_rows,
        "return_symmetry_rows": return_symmetry_rows,
        "sweep": sweep,
    }


def one_port_reflection_momentum_force_summary(
    power_incident_W,
    s11,
    incident_direction=(1.0, 0.0, 0.0),
    speed=C0,
):
    """Vector force [N] from one-port reflection coefficient ``S11``.

    This is the common VNA/RF-solver one-port specialization of
    :func:`two_port_scattering_momentum_force_summary`: reflected power is
    ``|S11|^2 P_inc`` and there is no transmitted output port.  The phase of
    ``S11`` is still reported for bookkeeping, but the time-average momentum
    force depends on its magnitude only:

        F = (1 + |S11|^2) P_inc k_inc / c

    The matched-load and perfect-short limits are therefore ``P/c`` and
    ``2P/c`` along the incident propagation direction.
    """

    gamma = complex(s11)
    if not math.isfinite(gamma.real) or not math.isfinite(gamma.imag):
        raise ValueError("s11 must be finite")
    magnitude = abs(gamma)
    reflectance = magnitude * magnitude
    summary = two_port_scattering_momentum_force_summary(
        power_incident_W,
        incident_direction,
        reflectance=reflectance,
        transmittance=0.0,
        speed=speed,
    )
    return_loss = None if magnitude == 0.0 else max(0.0, -20.0 * math.log10(magnitude))
    power_delivered = summary["power_absorbed_W"]
    delivered_fraction = summary["absorptance"]
    mismatch_loss = None if delivered_fraction == 0.0 else max(
        0.0,
        -10.0 * math.log10(delivered_fraction),
    )
    summary.update({
        "s11_real": gamma.real,
        "s11_imag": gamma.imag,
        "s11_magnitude": magnitude,
        "s11_phase_rad": math.atan2(gamma.imag, gamma.real),
        "s11_phase_deg": math.degrees(math.atan2(gamma.imag, gamma.real)),
        "return_loss_dB": return_loss,
        "return_loss_is_infinite": magnitude == 0.0,
        "power_delivered_to_one_port_W": power_delivered,
        "mismatch_loss_dB": mismatch_loss,
        "mismatch_loss_is_infinite": delivered_fraction == 0.0,
        "one_port_force_formula": "F=(1+|S11|^2)P_inc k_inc/c",
    })
    return summary


def one_port_reflection_sweep_momentum_force_summary(
    frequency_Hz,
    s11_values,
    power_incident_W=1.0,
    incident_direction=(1.0, 0.0, 0.0),
    speed=C0,
    passivity_tolerance=1.0e-12,
):
    """Summarize one-port reflection momentum force over a frequency sweep.

    ``frequency_Hz`` and ``s11_values`` are paired samples from a one-port
    reflection sweep.  The force model is the same as
    :func:`one_port_reflection_momentum_force_summary`:

        F = (1 + |S11|^2) P_inc k_inc / c

    Unlike the single-point helper, this sweep audit records passivity
    violations instead of failing immediately, so measured/simulated data with a
    small ``|S11| > 1`` overshoot can still be diagnosed.
    """

    frequencies = [float(value) for value in frequency_Hz]
    gammas = [complex(value) for value in s11_values]
    if len(frequencies) != len(gammas):
        raise ValueError("frequency_Hz and s11_values must have the same length")
    if not frequencies:
        raise ValueError("at least one frequency sample is required")
    power = float(power_incident_W)
    speed = float(speed)
    tolerance = float(passivity_tolerance)
    if power < 0.0:
        raise ValueError("power_incident_W must be >= 0")
    if speed <= 0.0:
        raise ValueError("speed must be > 0")
    if tolerance < 0.0:
        raise ValueError("passivity_tolerance must be >= 0")
    inc = _unit_vector(incident_direction, "incident_direction")

    rows = []
    violation_rows = []
    momentum_scale = power / speed
    for idx, (frequency, gamma) in enumerate(zip(frequencies, gammas)):
        if not math.isfinite(frequency) or frequency < 0.0:
            raise ValueError("frequency samples must be finite and >= 0")
        if not math.isfinite(gamma.real) or not math.isfinite(gamma.imag):
            raise ValueError("s11 values must be finite")
        magnitude = abs(gamma)
        reflectance = magnitude * magnitude
        factor = 1.0 + reflectance
        axial_force = factor * momentum_scale
        row = {
            "index": idx,
            "frequency_Hz": frequency,
            "s11_real": gamma.real,
            "s11_imag": gamma.imag,
            "s11_magnitude": magnitude,
            "s11_phase_rad": math.atan2(gamma.imag, gamma.real),
            "s11_phase_deg": math.degrees(math.atan2(gamma.imag, gamma.real)),
            "reflectance": reflectance,
            "absorptance_one_port": 1.0 - reflectance,
            "return_loss_dB": None if magnitude == 0.0 else max(0.0, -20.0 * math.log10(magnitude)),
            "return_loss_is_infinite": magnitude == 0.0,
            "momentum_transfer_factor": factor,
            "axial_force_along_incident_direction_N": axial_force,
            "force_magnitude_N": abs(axial_force),
            "force_N": [axial_force * value for value in inc],
            "power_reflected_W": reflectance * power,
            "power_delivered_to_one_port_W": (1.0 - reflectance) * power,
            "passivity_excess_magnitude": max(0.0, magnitude - 1.0),
            "passivity_excess_reflectance": max(0.0, reflectance - 1.0),
            "passivity_ok": magnitude <= 1.0 + tolerance,
        }
        rows.append(row)
        if not row["passivity_ok"]:
            violation_rows.append(row)

    max_force_row = max(rows, key=lambda row: row["force_magnitude_N"])
    min_force_row = min(rows, key=lambda row: row["force_magnitude_N"])
    max_reflectance_row = max(rows, key=lambda row: row["reflectance"])
    mean_force = sum(row["force_magnitude_N"] for row in rows) / len(rows)
    mean_reflectance = sum(row["reflectance"] for row in rows) / len(rows)
    monotonic = all(
        frequencies[idx] < frequencies[idx + 1]
        for idx in range(len(frequencies) - 1)
    )

    return {
        "n_points": len(rows),
        "frequency_min_Hz": min(frequencies),
        "frequency_max_Hz": max(frequencies),
        "frequency_monotonic_increasing": monotonic,
        "power_incident_W": power,
        "speed_m_per_s": speed,
        "incident_direction": inc,
        "passivity_tolerance": tolerance,
        "passivity_ok": not violation_rows,
        "passivity_violation_count": len(violation_rows),
        "max_s11_magnitude": max_reflectance_row["s11_magnitude"],
        "max_reflectance": max_reflectance_row["reflectance"],
        "max_reflectance_frequency_Hz": max_reflectance_row["frequency_Hz"],
        "max_passivity_excess_magnitude": max(row["passivity_excess_magnitude"] for row in rows),
        "max_passivity_excess_reflectance": max(row["passivity_excess_reflectance"] for row in rows),
        "mean_reflectance": mean_reflectance,
        "mean_force_magnitude_N": mean_force,
        "max_force_magnitude_N": max_force_row["force_magnitude_N"],
        "max_force_frequency_Hz": max_force_row["frequency_Hz"],
        "min_force_magnitude_N": min_force_row["force_magnitude_N"],
        "min_force_frequency_Hz": min_force_row["frequency_Hz"],
        "force_span_N": max_force_row["force_magnitude_N"] - min_force_row["force_magnitude_N"],
        "max_force_row": max_force_row,
        "min_force_row": min_force_row,
        "passivity_violation_rows": violation_rows,
        "rows": rows,
        "one_port_sweep_force_formula": "F=(1+|S11|^2)P_inc k_inc/c",
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
    1.8 % -> 0.5 % under such refinement."""
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
    validated on a magnetised cylinder in a uniform field (tau = m x B0)."""
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
    benchmark F/L = mu0 I1 I2 / (2 pi d) to ~1 % (see validation_test/radia_mcp/test_planar_force.py).
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
    (validation_test/radia_mcp/test_axi_force.py)."""
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
    static :func:`maxwell_surface_force` reference (validation_test/radia_mcp/test_maxwell_surface_harmonic.py).
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
    Kelvin functions, 0.07 %; see validation_test/radia_mcp/test_planar_eddy.py)."""
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
    busbars: |F| = mu0 I1 I2/(2 pi d)."""
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
    independent, to 0.06 % (validation_test/radia_mcp/test_planar_inductance.py)."""
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
