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
    stats["source_kind"] = "%s-bdm%d" % (
        stats.pop("source_representation"), int(order))
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
            "(got order=%r)." % (res.get("order"),))
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
