"""Persistent C++ field evaluator for HDiv BDM1/BDM2 solutions.

The solved HDiv field needs neither a piecewise-constant magnetization collapse
nor internal-face sources:
M is linear per tet (NGSolve HDiv order=1 on tets = full (P1)^3, BDM1-type), so

    * internal faces carry NO charge (HDiv conformity: M.n continuous),
    * BDM1 tetrahedra use linear boundary charge and constant volume charge,
    * BDM2 tetrahedra use quadratic boundary charge and linear volume charge,

and both pieces have closed forms in the C++ TET production kernel.  Affine HEX
sources use an exact six-TET physical-polynomial decomposition; curved HEX and
WEDGE sources use their tensor-product or prism-polynomial charge bases and
mapped quadrature clouds.  Solve-time
materialization stores the immutable C++ source evaluator in the result.  Its
NumPy-buffer API performs no per-call source packing, evaluates all IMA terms in
one TaskManager region, and selects a quadrupole treecode for sufficiently large
target-source work.  Flat tet leaves retain the analytic kernel; curved tet
leaves retain P2 geometry and integrate the BDM1/BDM2 charge polynomial directly.
"""
import hashlib
import time

import numpy as np

import radia._radia_pybind as _rp

MU0 = 4.0e-7 * np.pi
_FIELD_TREE_LEAF = 32
_FIELD_TREE_THETA = 0.05
_FIELD_TREE_MIN_SOURCES = 256
_FIELD_TREE_AUTO_MIN_WORK = 500_000_000
_FIELD_TREE_RELATIVE_TOLERANCE = 1.0e-5
_FIELD_TREE_PROBE_COUNT = 16

_FIELD_TREE_DEFAULTS = {
    "leaf_size": _FIELD_TREE_LEAF,
    "theta": _FIELD_TREE_THETA,
    "tree_min_sources": _FIELD_TREE_MIN_SOURCES,
    "auto_min_work": _FIELD_TREE_AUTO_MIN_WORK,
    "tree_relative_tolerance": _FIELD_TREE_RELATIVE_TOLERANCE,
    "probe_count": _FIELD_TREE_PROBE_COUNT,
}
_FIELD_TREE_INTEGER_OPTIONS = {
    "leaf_size": 1,
    "tree_min_sources": 1,
    "auto_min_work": 0,
    "probe_count": 1,
}


def _normalized_tree_options(tree_options):
    if not isinstance(tree_options, dict):
        raise TypeError("field-evaluator tree_options must be a dict")
    unknown = set(tree_options) - set(_FIELD_TREE_DEFAULTS)
    if unknown:
        raise ValueError(
            f"unknown field-evaluator tree options: {sorted(unknown)}")

    merged = dict(_FIELD_TREE_DEFAULTS)
    merged.update(tree_options)
    for name, minimum in _FIELD_TREE_INTEGER_OPTIONS.items():
        value = merged[name]
        numeric = float(value)
        if (isinstance(value, (bool, np.bool_)) or not np.isfinite(numeric)
                or not numeric.is_integer() or numeric < minimum):
            raise ValueError(
                f"field-evaluator {name} must be an integer >= {minimum}")
        merged[name] = int(numeric)
    for name, minimum, inclusive in (
        ("theta", 0.0, False),
        ("tree_relative_tolerance", 0.0, True),
    ):
        value = float(merged[name])
        valid_range = value >= minimum if inclusive else value > minimum
        if not np.isfinite(value) or not valid_range:
            relation = ">=" if inclusive else ">"
            raise ValueError(
                f"field-evaluator {name} must be finite and {relation} {minimum}")
        merged[name] = value
    return merged


def _create_field_evaluator(gram, coefficients, order, *, tree_options=None):
    """Materialize one immutable C++ source evaluator from configured geometry."""
    started = time.perf_counter()
    options = (
        dict(_FIELD_TREE_DEFAULTS)
        if tree_options is None
        else _normalized_tree_options(tree_options)
    )
    evaluator = gram.create_field_evaluator(
        np.ascontiguousarray(coefficients, dtype=np.float64),
        int(options["leaf_size"]), float(options["theta"]),
        int(options["tree_min_sources"]), int(options["auto_min_work"]),
        float(options["tree_relative_tolerance"]), int(options["probe_count"]))
    stats = dict(evaluator.stats())
    stats["source_kind"] = f"{stats.pop('source_representation')}-bdm{int(order)}"
    stats["build_wall_s"] = time.perf_counter()-started
    return evaluator, stats


def _materialize_field_evaluator(res, tree_options=None):
    """Build and cache the immutable C++ source evaluator.

    ``tree_options`` optionally overrides the Barnes-Hut controls
    ``(leaf_size, theta, tree_min_sources, auto_min_work,
    tree_relative_tolerance, probe_count)`` as a dict; each distinct
    override set is cached separately.  The strict default
    ``theta=0.05`` effectively disables far-field grouping (measured:
    tree == direct to 4e-15 on the C-type), which is the safe default;
    batch consumers that VERIFY the tree against direct may relax it.
    """
    if not isinstance(res, dict):
        raise TypeError("vim field evaluator requires Solve's result dict")
    gfM = res.get("gfM")
    if gfM is None:
        raise ValueError(
            "vim.FieldFromSolution: res carries no 'gfM' GridFunction -- pass the dict "
            "returned by vim.Solve/rad.Solve unmodified.")
    if int(res.get("order", -1)) not in (1, 2):
        raise NotImplementedError(
            "vim.FieldFromSolution: production supports HDiv order in {1,2} "
            f"(got order={res.get('order')!r}).")
    if tree_options is not None:
        merged = _normalized_tree_options(tree_options)
        key = tuple(sorted(merged.items()))
        cache = res.setdefault("_field_evaluator_custom", {})
        if not isinstance(cache, dict):
            raise TypeError(
                "vim.FieldFromSolution: custom field-evaluator cache is invalid")
        cached = cache.get(key)
        if cached is not None:
            return cached
        gram = res.get("_charge_gram")
        coefficients = res.get("_m_coefficients")
        if gram is None or coefficients is None:
            raise ValueError(
                "vim.FieldFromSolution: result does not carry the configured "
                "C++ operator and solution vector.")
        evaluator = gram.create_field_evaluator(
            np.ascontiguousarray(coefficients, dtype=np.float64),
            int(merged["leaf_size"]), float(merged["theta"]),
            int(merged["tree_min_sources"]), int(merged["auto_min_work"]),
            float(merged["tree_relative_tolerance"]),
            int(merged["probe_count"]))
        cache[key] = evaluator
        return evaluator
    cached = res.get("_field_evaluator")
    if cached is not None:
        return cached
    gram = res.get("_charge_gram")
    coefficients = res.get("_m_coefficients")
    if gram is None or coefficients is None:
        raise ValueError(
            "vim.FieldFromSolution: result does not carry the configured C++ operator and solution vector; "
            "pass the current vim.Solve/rad.Solve result unmodified.")
    evaluator, stats = _create_field_evaluator(
        gram, coefficients, int(res.get("order", 1)))
    res["_field_evaluator"] = evaluator
    res["field_evaluator_stats"] = stats
    res["field_evaluator_build_wall_s"] = stats["build_wall_s"]
    return evaluator


def field_from_solution(res, points, algorithm="auto", *, tree_options=None):
    """Demagnetizing field H_demag (A/m) of a solved HDiv-VIM magnetization at
    ``points`` (N,3), evaluated from the BDM1/BDM2 solution directly -- no per-element
    constant-M collapse, hence none of the near-surface piecewise-constant ripple of
    ``rad.Fld`` on the write-back elements (the O(h) bumps measured at standoff ~
    element size disappear identically; see the module docstring).

    ``res`` is the dict returned by ``vim.Solve`` / ``rad.Solve`` on a MeshSoftIron
    (it must carry the ``gfM`` GridFunction; order 1 or 2 on the production topology).
    Iron contribution only:

        B outside the iron = MU0 * (H_ext + H_demag)
        B inside  the iron = MU0 * (H_ext + H_demag + M)

    ``algorithm="direct"`` is the exact discrete source sum.  The default
    ``"auto"`` uses direct evaluation for ordinary batches and considers the
    quadrupole tree only above the large-work threshold; representative points
    must satisfy the configured direct-reference tolerance and show a measured
    speed benefit before the full batch uses the tree."""
    pts = np.ascontiguousarray(np.asarray(points, float).reshape(-1, 3))
    evaluator = _materialize_field_evaluator(res, tree_options)
    return np.asarray(evaluator.field(pts, str(algorithm)), float)/(4.0*np.pi)


def field_coefficient_from_solution(
    res, algorithm="direct", *, reflection_normal=None
):
    """Return the persistent HDiv demagnetizing field as an NGSolve CF.

    This is the zero-copy coupling surface for another independently meshed
    HDiv body.  Independent spaces preserve normal-magnetization jumps at
    touching permanent-magnet/iron and segmented-magnet interfaces.

    ``reflection_normal`` symmetrizes the source about the plane through the
    origin with that normal.  The reflected magnetic field is transformed as
    an axial vector, ``det(R) R B = -R B``.  This is the B/H counterpart of the
    full-volume vector-potential source symmetrization used by EarlyTimes.
    """
    if algorithm not in ("direct", "tree"):
        raise ValueError(
            "vim.FieldCoefficientFromSolution: algorithm must be 'direct' or 'tree'")
    evaluator = _materialize_field_evaluator(res)
    normal = []
    if reflection_normal is not None:
        normal = np.asarray(reflection_normal, dtype=float)
        if normal.shape != (3,) or not np.all(np.isfinite(normal)):
            raise ValueError("reflection_normal must be a finite three-vector")
        norm = float(np.linalg.norm(normal))
        if norm <= 0.0:
            raise ValueError("reflection_normal must be nonzero")
        normal = (normal / norm).tolist()
    return _rp._HDivFieldCoefficient(evaluator, str(algorithm), normal)


def _mapped_vector_values(value, mapped_rule, count, label):
    sampled = np.asarray(value(mapped_rule), dtype=float)
    if sampled.shape == (3, count):
        sampled = sampled.T
    sampled = sampled.reshape(count, 3)
    if not np.all(np.isfinite(sampled)):
        raise RuntimeError(f"{label} contains non-finite values")
    return np.ascontiguousarray(sampled)


def _target_adaptive_integration_order(
    mesh,
    element,
    target_points,
    minimum_order,
    maximum_order,
):
    """Choose a conservative even-order increment from source/target clearance."""
    vertices = np.asarray(
        [mesh[vertex].point for vertex in element.vertices], dtype=float
    )
    center = np.mean(vertices, axis=0)
    radius = float(
        np.max(np.linalg.norm(vertices - center, axis=1), initial=0.0)
    )
    if radius <= 0.0:
        return int(maximum_order)
    clearance = max(
        0.0,
        float(
            np.min(np.linalg.norm(target_points - center, axis=1)) - radius
        ),
    )
    ratio = clearance / radius
    # One two-order increment for every factor-of-two approach from four
    # element radii.  A target intersecting the conservative element sphere
    # receives the requested maximum order.
    if ratio <= 0.0:
        levels = (int(maximum_order) - int(minimum_order) + 1) // 2
    else:
        levels = max(0, int(np.ceil(np.log2(4.0 / ratio))))
    selected = int(minimum_order) + 2 * levels
    return min(int(maximum_order), max(int(minimum_order), selected))


_TET_REFERENCE_SAMPLES = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (0.25, 0.25, 0.25),
)
_TET_FACE_VERTICES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
_AFFINE_RELATIVE_TOLERANCE = 1.0e-9


def _cross_matrix_right(normal):
    """Return ``N`` with ``N @ v == v x normal`` (note the operand order)."""
    n_x, n_y, n_z = normal
    return np.array(
        [[0.0, n_z, -n_y], [-n_z, 0.0, n_x], [n_y, -n_x, 0.0]], dtype=float
    )


def _affine_equivalent_current_sources(gfM):
    """Return exact equivalent-current sources of an affine (BDM1) tet field.

    ``A(r) = mu0/(4 pi) [ INT_V (curl M)/R dV' + INT_dV (M x n)/R dS' ]`` is an
    identity, not a quadrature.  On a straight tetrahedron carrying an affine
    ``M`` the volume current is constant and each face current is affine, so the
    integral closes analytically at any standoff.  Both sides of every interior
    face are retained: HDiv keeps ``M.n`` continuous but deliberately allows the
    tangential jump that carries the equivalent sheet current.

    The element geometry and the affine magnetization are both read through
    NGSolve's own element transformation, so no reference-vertex ordering
    convention is assumed.  A curved element or a non-affine magnetization is a
    hard error, not a silently degraded source.
    """
    import ngsolve as ng

    mesh = gfM.space.mesh
    rule = ng.IntegrationRule(
        [tuple(point) for point in _TET_REFERENCE_SAMPLES],
        [1.0 / 24.0] * len(_TET_REFERENCE_SAMPLES),
    )
    coordinate = ng.CoefficientFunction((ng.x, ng.y, ng.z))
    count = len(_TET_REFERENCE_SAMPLES)
    tetrahedron_vertices = []
    tetrahedron_current = []
    triangle_vertices = []
    triangle_current = []
    triangle_current_gradient = []
    maximum_geometry_defect = 0.0
    maximum_affine_defect = 0.0
    maximum_magnetization = 0.0
    for element in mesh.Elements(ng.VOL):
        if element.type != ng.ET.TET:
            raise NotImplementedError(
                "vim.VectorPotentialCoefficientFromSolution: the exact "
                "equivalent-current construction supports straight "
                f"tetrahedra only (got element type {element.type!r})"
            )
        mapped = mesh.GetTrafo(element)(rule)
        points = _mapped_vector_values(
            coordinate, mapped, count, "mapped tetrahedron sample points"
        )
        magnetization = _mapped_vector_values(
            gfM, mapped, count, "HDiv magnetization samples"
        )
        vertices = points[:4]
        scale = float(
            np.max(np.linalg.norm(vertices - vertices.mean(axis=0), axis=1))
        )
        if not scale > 0.0:
            raise RuntimeError("degenerate tetrahedron in the HDiv source mesh")
        geometry_defect = float(
            np.linalg.norm(points[4] - vertices.mean(axis=0))
        )
        maximum_geometry_defect = max(maximum_geometry_defect, geometry_defect / scale)

        design = np.column_stack((np.ones(4), vertices))
        coefficients = np.linalg.solve(design, magnetization[:4])
        predicted = np.concatenate(([1.0], points[4])) @ coefficients
        magnitude = float(np.max(np.abs(magnetization), initial=0.0))
        maximum_magnetization = max(maximum_magnetization, magnitude)
        maximum_affine_defect = max(
            maximum_affine_defect,
            float(np.max(np.abs(predicted - magnetization[4]), initial=0.0))
            / max(magnitude, 1.0),
        )

        constant = coefficients[0]
        # gradient[i][j] = d M_i / d x_j
        gradient = coefficients[1:].T
        tetrahedron_vertices.append(vertices)
        tetrahedron_current.append(
            [
                gradient[2][1] - gradient[1][2],
                gradient[0][2] - gradient[2][0],
                gradient[1][0] - gradient[0][1],
            ]
        )
        for face in _TET_FACE_VERTICES:
            corners = vertices[list(face)]
            opposite = vertices[
                next(index for index in range(4) if index not in face)
            ]
            normal = np.cross(corners[1] - corners[0], corners[2] - corners[0])
            norm = float(np.linalg.norm(normal))
            if not norm > 0.0:
                raise RuntimeError("degenerate face in the HDiv source mesh")
            normal = normal / norm
            if float(np.dot(normal, opposite - corners[0])) > 0.0:
                normal = -normal
            cross_right = _cross_matrix_right(normal)
            triangle_vertices.append(corners)
            triangle_current.append(cross_right @ constant)
            triangle_current_gradient.append(cross_right @ gradient)

    if not tetrahedron_vertices:
        raise ValueError("HDiv vector-potential source mesh contains no volume")
    if maximum_geometry_defect > _AFFINE_RELATIVE_TOLERANCE:
        raise NotImplementedError(
            "vim.VectorPotentialCoefficientFromSolution: the exact "
            "equivalent-current construction requires straight tetrahedra "
            f"(relative curvature defect {maximum_geometry_defect:.3e})"
        )
    if maximum_affine_defect > _AFFINE_RELATIVE_TOLERANCE:
        raise NotImplementedError(
            "vim.VectorPotentialCoefficientFromSolution: the exact "
            "equivalent-current construction requires an affine (BDM1) "
            "magnetization; the supplied space is not affine per element "
            f"(relative defect {maximum_affine_defect:.3e})"
        )
    return {
        "tetrahedron_vertices": np.ascontiguousarray(
            tetrahedron_vertices, dtype=np.float64
        ),
        "tetrahedron_current": np.ascontiguousarray(
            tetrahedron_current, dtype=np.float64
        ),
        "triangle_vertices": np.ascontiguousarray(
            triangle_vertices, dtype=np.float64
        ),
        "triangle_current": np.ascontiguousarray(
            triangle_current, dtype=np.float64
        ),
        "triangle_current_gradient": np.ascontiguousarray(
            triangle_current_gradient, dtype=np.float64
        ),
        "maximum_relative_curvature_defect": maximum_geometry_defect,
        "maximum_relative_affine_defect": maximum_affine_defect,
        "maximum_magnetization_A_m": maximum_magnetization,
    }


def _reflect_equivalent_current_sources(sources, normal):
    """Median-plane-symmetrize the exact equivalent-current source at half weight.

    Positions reflect as polar vectors.  ``M`` is axial, so both the volume
    current ``curl M`` and the sheet current ``M x n`` reflect as polar vectors;
    this is the same axial-source convention as the B route.
    """
    reflection = np.eye(3) - 2.0 * np.outer(normal, normal)
    reflected = {}
    for key in ("tetrahedron_vertices", "triangle_vertices"):
        original = sources[key]
        reflected[key] = np.ascontiguousarray(
            np.concatenate((original, original @ reflection.T))
        )
    for key in ("tetrahedron_current", "triangle_current"):
        original = sources[key]
        reflected[key] = np.ascontiguousarray(
            0.5 * np.concatenate((original, original @ reflection.T))
        )
    original = sources["triangle_current_gradient"]
    reflected["triangle_current_gradient"] = np.ascontiguousarray(
        0.5
        * np.concatenate(
            (original, np.einsum("ij,njk,kl->nil", reflection, original, reflection))
        )
    )
    return reflected


def vector_potential_coefficient_from_solution(
    res,
    integration_order=8,
    *,
    construction="auto",
    reflection_normal=None,
    target_points_m=None,
    maximum_integration_order=None,
):
    """Return exterior ``A`` integrated directly from the solved HDiv ``M``.

    The NGSolve ``GridFunction`` is evaluated on mapped element quadrature
    rules, including its native HDiv orientation and Piola transform.  The
    immutable native coefficient then evaluates

    ``A(r) = mu0/(4*pi) * integral M(r') x (r-r') / |r-r'|^3 dV'``.

    This is a source integral, not a reconstruction from ``curl(A)`` or from
    median-plane samples.  It is deliberately an *exterior* coefficient: all
    EarlyTimes beam-tube points must remain outside the iron.  The caller owns
    the surrounding :class:`ngsolve.TaskManager`.

    ``reflection_normal`` optionally applies full-volume median-plane
    symmetrisation.  Source points are reflected as polar positions, while
    magnetization is transformed as an axial vector.  The original and
    reflected source clouds receive one-half weight each; this retains all
    upper/lower volume information without averaging HCurl traces.

    ``construction`` selects the source representation explicitly.

    ``"exact"`` uses the analytic equivalent-current identity
    ``A = mu0/(4 pi) [INT (curl M)/R dV + INT (M x n)/R dS]`` with the closed-form
    tetrahedron/triangle Newtonian potentials, matching the analytic kernels the
    B route already uses.  It is the only construction that stays accurate when
    a target sits much closer to a source element than that element's own size,
    which is the EarlyTimes beam-tube geometry.  It requires a straight
    tetrahedral mesh carrying an affine (order-1/BDM1) magnetization and ignores
    ``integration_order``/``target_points_m``.

    ``"quadrature"`` keeps the point-dipole cloud below.  It is the available
    representation for element/space classes without an analytic kernel, and it
    remains useful as an independent cross-check of the exact construction.

    ``"auto"`` (the default) selects ``"exact"`` for an order-1 tetrahedral HDiv
    space and ``"quadrature"`` otherwise.  A curved tetrahedral mesh is a hard
    error under ``"exact"``; request ``"quadrature"`` explicitly for it.  The
    selected construction is always recorded in
    ``res["vector_potential_coefficient_stats"]``.

    ``integration_order`` is an explicit convergence control and is the
    minimum order in target-adaptive mode.  There is no fixed production
    order: accept a map only after consecutive orders change the independent
    A/B field and tracking observables by less than their declared tolerances.

    Passing ``target_points_m`` together with ``maximum_integration_order``
    activates conservative target-adaptive quadrature.  Each source element
    is assigned an even order increment from the clearance between its vertex
    bounding sphere and the supplied beam-volume points.  Elements within the
    conservative sphere receive the maximum order; the order drops by two for
    every factor-of-two increase in clearance and never falls below
    ``integration_order``.  Supply the full loft-chain vertex cloud, not only
    the design orbit, so the certified aperture controls the source accuracy.
    """
    if not isinstance(res, dict):
        raise TypeError("vim vector potential requires Solve's result dict")
    gfM = res.get("gfM")
    if gfM is None:
        raise ValueError(
            "vim.VectorPotentialCoefficientFromSolution: res carries no "
            "'gfM' GridFunction -- pass the dict returned by vim.Solve/"
            "rad.Solve unmodified."
        )
    if construction not in ("auto", "exact", "quadrature"):
        raise ValueError(
            "vim.VectorPotentialCoefficientFromSolution: construction must be "
            "'auto', 'exact', or 'quadrature'"
        )
    normal = None
    if reflection_normal is not None:
        normal = np.asarray(reflection_normal, dtype=float)
        if normal.shape != (3,) or not np.all(np.isfinite(normal)):
            raise ValueError("reflection_normal must be a finite three-vector")
        norm = float(np.linalg.norm(normal))
        if norm <= 0.0:
            raise ValueError("reflection_normal must be nonzero")
        normal = normal / norm
    normal_key = None if normal is None else tuple(np.round(normal, decimals=15))

    selected = construction
    if selected == "auto":
        import ngsolve as ng

        tetrahedral = all(
            element.type == ng.ET.TET for element in gfM.space.mesh.Elements(ng.VOL)
        )
        selected = (
            "exact"
            if tetrahedral and int(res.get("order", -1)) == 1
            else "quadrature"
        )
    if selected == "exact":
        cache_key = ("exact-equivalent-current", normal_key)
        cache = res.setdefault("_vector_potential_coefficients", {})
        if cache_key in cache:
            return cache[cache_key]
        sources = _affine_equivalent_current_sources(gfM)
        diagnostics = {
            key: sources.pop(key)
            for key in (
                "maximum_relative_curvature_defect",
                "maximum_relative_affine_defect",
                "maximum_magnetization_A_m",
            )
        }
        if normal is not None:
            sources = _reflect_equivalent_current_sources(sources, normal)
        coefficient = _rp._HDivExactVectorPotentialCoefficient(
            sources["tetrahedron_vertices"].ravel(),
            sources["tetrahedron_current"].ravel(),
            sources["triangle_vertices"].ravel(),
            sources["triangle_current"].ravel(),
            sources["triangle_current_gradient"].ravel(),
        )
        res.setdefault("vector_potential_coefficient_stats", {})[cache_key] = {
            "construction": "analytic-equivalent-current-tet-triangle",
            "requested_construction": construction,
            "tetrahedron_count": len(sources["tetrahedron_current"]),
            "triangle_count": len(sources["triangle_current"]),
            "element_count": int(gfM.space.mesh.ne),
            "reflection_normal": None if normal is None else normal.tolist(),
            "full_volume_reflection_symmetrized": normal is not None,
            **diagnostics,
        }
        cache[cache_key] = coefficient
        return coefficient

    order = int(integration_order)
    if isinstance(integration_order, bool) or order != integration_order or order < 1:
        raise ValueError("integration_order must be a positive integer")
    targets = None
    maximum_order = None
    if target_points_m is not None:
        targets = np.ascontiguousarray(
            np.asarray(target_points_m, dtype=float).reshape(-1, 3)
        )
        if len(targets) == 0 or not np.all(np.isfinite(targets)):
            raise ValueError("target_points_m must contain finite 3D points")
        if maximum_integration_order is None:
            raise ValueError(
                "maximum_integration_order is required with target_points_m"
            )
        maximum_order = int(maximum_integration_order)
        if (
            isinstance(maximum_integration_order, bool)
            or maximum_order != maximum_integration_order
            or maximum_order < order
        ):
            raise ValueError(
                "maximum_integration_order must be an integer not smaller "
                "than integration_order"
            )
    elif maximum_integration_order is not None:
        raise ValueError(
            "target_points_m is required with maximum_integration_order"
        )
    if targets is None:
        cache_key = order if normal_key is None else (order, normal_key)
    else:
        target_digest = hashlib.sha256(targets.view(np.uint8)).hexdigest()
        cache_key = (
            "target-adaptive",
            order,
            maximum_order,
            targets.shape,
            target_digest,
            normal_key,
        )
    cache = res.setdefault("_vector_potential_coefficients", {})
    if cache_key in cache:
        return cache[cache_key]

    import ngsolve as ng

    mesh = gfM.space.mesh
    coordinate = ng.CoefficientFunction((ng.x, ng.y, ng.z))
    rules = {}
    source_points = []
    integrated_magnetization = []
    integrated_volume = 0.0
    order_histogram = {}
    for element in mesh.Elements(ng.VOL):
        element_order = (
            order
            if targets is None
            else _target_adaptive_integration_order(
                mesh,
                element,
                targets,
                order,
                maximum_order,
            )
        )
        order_histogram[element_order] = order_histogram.get(element_order, 0) + 1
        rule_key = (element.type, element_order)
        rule = rules.get(rule_key)
        if rule is None:
            rule = ng.IntegrationRule(element.type, element_order)
            rules[rule_key] = rule
        transformation = mesh.GetTrafo(element)
        mapped = transformation(rule)
        count = len(rule)
        points = _mapped_vector_values(
            coordinate, mapped, count, "mapped quadrature points"
        )
        magnetization = _mapped_vector_values(
            gfM, mapped, count, "HDiv magnetization samples"
        )
        weights = np.fromiter(
            (
                float(point.weight) * float(transformation(point).measure)
                for point in rule
            ),
            dtype=np.float64,
            count=count,
        )
        if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
            raise RuntimeError(
                "HDiv vector-potential quadrature has a non-positive mapped weight"
            )
        source_points.append(points)
        integrated_magnetization.append(magnetization * weights[:, None])
        integrated_volume += float(np.sum(weights))
    if not source_points or integrated_volume <= 0.0:
        raise ValueError("HDiv vector-potential source mesh contains no volume")

    points = np.ascontiguousarray(np.vstack(source_points), dtype=np.float64)
    moments = np.ascontiguousarray(
        np.vstack(integrated_magnetization), dtype=np.float64
    )
    if normal is not None:
        reflected_points = points - 2.0 * (points @ normal)[:, None] * normal
        # An axial vector transforms under a reflection R as det(R) R M = -R M.
        reflected_moments = (
            -moments + 2.0 * (moments @ normal)[:, None] * normal
        )
        points = np.ascontiguousarray(
            np.vstack((points, reflected_points)), dtype=np.float64
        )
        moments = np.ascontiguousarray(
            0.5 * np.vstack((moments, reflected_moments)), dtype=np.float64
        )
    coefficient = _rp._HDivVectorPotentialCoefficient(points, moments)
    res.setdefault("vector_potential_coefficient_stats", {})[cache_key] = {
        "construction": "ngsolve-mapped-hdiv-magnetization-volume-integral",
        "requested_construction": construction,
        "integration_order": order,
        "maximum_integration_order": maximum_order,
        "target_adaptive": targets is not None,
        "target_point_count": 0 if targets is None else len(targets),
        "element_integration_order_histogram": dict(sorted(order_histogram.items())),
        "source_count": len(points),
        "element_count": int(mesh.ne),
        "integrated_volume_m3": integrated_volume,
        "integrated_magnetic_moment_A_m2": np.sum(moments, axis=0).tolist(),
        "reflection_normal": None if normal is None else normal.tolist(),
        "full_volume_reflection_symmetrized": normal is not None,
    }
    cache[cache_key] = coefficient
    return coefficient


def magnetization_from_solution(res, points):
    """Evaluate the BDM1 magnetization inside the solved mesh, zero outside."""
    import ngsolve as ng

    gfM = res.get("gfM") if isinstance(res, dict) else None
    if gfM is None:
        raise ValueError("vim magnetization evaluation requires Solve's result dict")
    pts = np.asarray(points, float).reshape(-1, 3)
    values = np.zeros((len(pts), 3), dtype=float)
    mesh = gfM.space.mesh
    if not len(pts):
        return values
    # NGSolve's ndarray MeshPoint representation exposes nr=-1 for points
    # outside the volume, allowing one vectorized lookup and one batched GF
    # evaluation without a Python point loop.
    with ng.TaskManager():
        mapped = mesh(pts[:, 0], pts[:, 1], pts[:, 2])
        if not (isinstance(mapped, np.ndarray) and mapped.dtype.names and "nr" in mapped.dtype.names):
            raise RuntimeError(
                "vim magnetization evaluation requires NGSolve's vectorized MeshPoint array API")
        valid = np.asarray(mapped["nr"] >= 0)
        if np.any(valid):
            values[valid] = np.asarray(gfM(mapped[valid]), float).reshape(-1, 3)
    return values
