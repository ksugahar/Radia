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


def _create_field_evaluator(gram, coefficients, order):
    """Materialize one immutable C++ source evaluator from configured geometry."""
    started = time.perf_counter()
    evaluator = gram.create_field_evaluator(
        np.ascontiguousarray(coefficients, dtype=np.float64),
        _FIELD_TREE_LEAF, _FIELD_TREE_THETA,
        _FIELD_TREE_MIN_SOURCES, _FIELD_TREE_AUTO_MIN_WORK,
        _FIELD_TREE_RELATIVE_TOLERANCE, _FIELD_TREE_PROBE_COUNT)
    stats = dict(evaluator.stats())
    stats["source_kind"] = f"{stats.pop('source_representation')}-bdm{int(order)}"
    stats["build_wall_s"] = time.perf_counter()-started
    return evaluator, stats


def _materialize_field_evaluator(res):
    """Build and cache the immutable C++ source evaluator exactly once."""
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


def field_from_solution(res, points, algorithm="auto"):
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
    evaluator = _materialize_field_evaluator(res)
    return np.asarray(evaluator.field(pts, str(algorithm)), float)/(4.0*np.pi)


def field_coefficient_from_solution(res, algorithm="direct"):
    """Return the persistent HDiv demagnetizing field as an NGSolve CF.

    This is the zero-copy coupling surface for another independently meshed
    HDiv body.  Independent spaces preserve normal-magnetization jumps at
    touching permanent-magnet/iron and segmented-magnet interfaces.
    """
    if algorithm not in ("direct", "tree"):
        raise ValueError(
            "vim.FieldCoefficientFromSolution: algorithm must be 'direct' or 'tree'")
    evaluator = _materialize_field_evaluator(res)
    return _rp._HDivFieldCoefficient(evaluator, str(algorithm))


def _mapped_vector_values(value, mapped_rule, count, label):
    sampled = np.asarray(value(mapped_rule), dtype=float)
    if sampled.shape == (3, count):
        sampled = sampled.T
    sampled = sampled.reshape(count, 3)
    if not np.all(np.isfinite(sampled)):
        raise RuntimeError(f"{label} contains non-finite values")
    return np.ascontiguousarray(sampled)


def vector_potential_coefficient_from_solution(
    res, integration_order=8, *, reflection_normal=None
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

    ``integration_order`` is an explicit convergence control.  Order 8 is the
    production baseline for BDM1/BDM2 on affine iron elements; a C-type map is
    accepted only after the order-10 result changes the HCurl/Lie observables
    by less than their declared tolerances.
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
    order = int(integration_order)
    if isinstance(integration_order, bool) or order != integration_order or order < 1:
        raise ValueError("integration_order must be a positive integer")
    normal = None
    if reflection_normal is not None:
        normal = np.asarray(reflection_normal, dtype=float)
        if normal.shape != (3,) or not np.all(np.isfinite(normal)):
            raise ValueError("reflection_normal must be a finite three-vector")
        norm = float(np.linalg.norm(normal))
        if norm <= 0.0:
            raise ValueError("reflection_normal must be nonzero")
        normal = normal / norm
    cache_key = (
        order
        if normal is None
        else (order, tuple(np.round(normal, decimals=15)))
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
    for element in mesh.Elements(ng.VOL):
        rule = rules.get(element.type)
        if rule is None:
            rule = ng.IntegrationRule(element.type, order)
            rules[element.type] = rule
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
        "integration_order": order,
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
