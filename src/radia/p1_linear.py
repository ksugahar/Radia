"""Linear total-A magnetostatics with one lowest-order beta-zero AMS-PCG solve.

Call outside ``ngsolve.TaskManager``, as for ``p1_newton``: AMS builds its
hierarchy outside a parallel region. This function owns its assembly/solve
regions; the caller selects the NGSolve thread count.
"""
from __future__ import annotations

from collections.abc import Mapping
import math
import operator
import re
import time

import ngsolve as ng
import numpy as np

from radia._residual_gate import RELATIVE_LIMIT

MU0 = 4.0e-7 * math.pi


def _native():
    import radia.sparsesolv_ngsolve as native

    for name in ("LowestOrderCurlSystem", "LowestOrderGradient", "NativePCG",
                 "HypreBasedAMSPreconditioner", "TaskManagerActive"):
        if not hasattr(native, name):
            raise RuntimeError(f"radia.sparsesolv_ngsolve lacks {name}; rebuild the native module")
    return native


def solve_p1_linear(mesh, *, current_cf, current_materials, mu_r_dict=None,
                    dirichlet="outer", cg_tolerance=1e-8, cg_max_iterations=2000,
                    mixed_precision=True, observation_points=None):
    """Solve ``int nu curl A.curl v = int J.v`` on straight 3D tetrahedra.

    ``mu_r_dict`` maps material labels to positive finite constant scalar
    relative permeabilities; omitted labels are vacuum. ``current_cf`` is a
    real vector CoefficientFunction in A/m^2 on ``current_materials``. Its
    assembled load must be orthogonal to the discrete gradients on free
    dofs. An incompatible load raises: no projection or mass gauge is added.
    The compatibility check includes all vertex gradients, as in p1_newton;
    it can reject a current reaching the Dirichlet boundary even when a
    weaker compatibility condition would suffice. Interior sources avoid
    this restriction.

    Uses HCurl(order=1, nograds=True), homogeneous tangential Dirichlet data,
    beta-zero AMS and exactly one native PCG call. The original masked
    equation residual must satisfy both the requested tolerance and the
    repository relative-residual limit. Nonconvergence raises.

    Returns ``A``, ``B_cf``, ``H_cf``, elementwise ``nu``, ``fes`` and ``stats``
    with phase timings, iterations and independently checked true residual.
    Optional observations are returned as ``observation_B_T``. This function
    must be entered outside TaskManager, even for a zero current.
    """
    native = _native()
    if native.TaskManagerActive():
        raise RuntimeError("solve_p1_linear must be called outside ngsolve.TaskManager: "
                           "AMS setup requires it; the solver owns its parallel regions")
    tolerance = float(cg_tolerance)
    if not math.isfinite(tolerance) or not 0.0 < tolerance < 1.0:
        raise ValueError("cg_tolerance must lie in (0, 1)")
    if isinstance(cg_max_iterations, (bool, np.bool_)):
        raise ValueError("cg_max_iterations must be a positive integer")
    try:
        max_iterations = operator.index(cg_max_iterations)
    except TypeError as exc:
        raise ValueError("cg_max_iterations must be a positive integer") from exc
    if not 0 < max_iterations <= np.iinfo(np.int32).max:
        raise ValueError("cg_max_iterations must be a positive int32")
    if not isinstance(mixed_precision, (bool, np.bool_)):
        raise ValueError("mixed_precision must be boolean")
    if mesh.dim != 3 or mesh.GetCurveOrder() != 1:
        raise ValueError("solve_p1_linear requires straight three-dimensional tetrahedra")
    if getattr(current_cf, "dim", None) != 3 or getattr(current_cf, "is_complex", False):
        raise ValueError("current_cf must be a real vector CoefficientFunction of dimension 3")
    materials = tuple(mesh.GetMaterials())
    names = (current_materials,) if isinstance(current_materials, str) else tuple(current_materials)
    if not names or not set(names) <= set(materials):
        raise ValueError("current_materials names must exist in the mesh")
    permeability = {} if mu_r_dict is None else mu_r_dict
    if not isinstance(permeability, Mapping) or not set(permeability) <= set(materials):
        raise ValueError("mu_r_dict must map existing material names to scalar relative permeabilities")
    relative_mu = np.asarray([float(permeability.get(name, 1.0)) for name in materials])
    if not np.all(np.isfinite(relative_mu)) or np.any(relative_mu <= 0.0):
        raise ValueError("relative permeabilities must be finite and positive")
    observation = None
    if observation_points is not None:
        observation = np.asarray(observation_points, dtype=float)
        if observation.ndim != 2 or observation.shape[1] != 3 or not np.all(np.isfinite(observation)):
            raise ValueError("observation_points must be a finite (n, 3) array")

    started = time.perf_counter()
    timings = {}
    t0 = time.perf_counter()
    elements = mesh.ngmesh.Elements3D().NumPy()
    if len(elements) != mesh.ne or np.any(elements["np"] != 4):
        raise ValueError("solve_p1_linear requires straight tetrahedra")
    material_index = np.asarray(elements["index"], dtype=np.int64) - 1
    if np.any(material_index < 0) or np.any(material_index >= len(materials)):
        raise ValueError("unexpected volume material numbering")
    nu_values = np.ascontiguousarray(1.0 / (MU0 * relative_mu[material_index]))
    if not np.all(np.isfinite(nu_values)) or np.any(nu_values <= 0.0):
        raise ValueError("permeabilities give non-finite or non-positive reluctivity")
    fes = ng.HCurl(mesh, order=1, dirichlet=dirichlet, nograds=True)
    free = np.fromiter(fes.FreeDofs(), dtype=bool, count=fes.ndof)
    solution = ng.GridFunction(fes, name="A_p1_linear")
    timings["space_material_s"] = time.perf_counter() - t0

    with ng.TaskManager():
        t0 = time.perf_counter()
        system = native.LowestOrderCurlSystem(fes, coefficient=nu_values)
        if not np.all(np.isfinite(system["volume"])) or np.any(system["volume"] <= 0.0):
            raise ValueError("non-positive or non-finite element volume")
        matrix = system["matrix"]
        del system  # Element/Newton geometry arrays are not needed by a linear solve.
        timings["element_system_s"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        gradient = native.LowestOrderGradient(fes)
        timings["gradient_s"] = time.perf_counter() - t0
        t0 = time.perf_counter()
        current_form = ng.LinearForm(fes)
        current_form += ng.InnerProduct(current_cf, fes.TestFunction()) * ng.dx(
            definedon=mesh.Materials("|".join(map(re.escape, names))))
        current_form.Assemble()
        rhs = current_form.vec
        load = rhs.FV().NumPy()
        if not np.all(np.isfinite(load)):
            raise ValueError("current load contains non-finite values")
        timings["current_assembly_s"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    from scipy.sparse import csr_matrix

    values, columns, offsets = gradient.CSR()
    g = csr_matrix((np.array(values), np.array(columns), np.array(offsets)),
                   shape=(gradient.height, gradient.width))
    free_load = np.where(free, load, 0.0)
    defect = float(np.linalg.norm(g.T @ free_load))
    scale = float(np.linalg.norm(abs(g).T @ np.abs(free_load)))
    if not math.isfinite(defect) or not math.isfinite(scale):
        raise ValueError("non-finite discrete-gradient load compatibility")
    if defect > 1e-10 * scale:
        raise ValueError("current load is not orthogonal to the discrete gradients; "
                         "use a weakly divergence-free current (radia.meshed_current)")
    timings["compatibility_s"] = time.perf_counter() - t0
    del g, values, columns, offsets

    t0 = time.perf_counter()
    coordinates = np.asarray(mesh.ngmesh.Coordinates())
    coord_x, coord_y, coord_z = (coordinates[:, k].tolist() for k in range(3))
    timings["coordinates_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    ams = native.HypreBasedAMSPreconditioner(
        mat=matrix, grad_mat=gradient, freedofs=fes.FreeDofs(),
        coord_x=coord_x, coord_y=coord_y, coord_z=coord_z,
        cycle_type=1, print_level=0, beta_zero=True, reuse_hierarchy=False,
        mixed_precision=bool(mixed_precision))
    timings["ams_setup_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    pcg = native.NativePCG(matrix, ams, fes.FreeDofs())
    timings["pcg_setup_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    with ng.TaskManager():
        iterations, reported, converged = pcg.Solve(rhs, solution.vec, tolerance, max_iterations)
    timings["solve_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    residual = rhs.CreateVector()
    with ng.TaskManager():
        residual.data = rhs - matrix * solution.vec
    residual_norm = float(np.linalg.norm(residual.FV().NumPy()[free]))
    rhs_norm = float(np.linalg.norm(free_load))
    relative = residual_norm / rhs_norm if rhs_norm > 0.0 else (0.0 if residual_norm == 0.0 else math.inf)
    limit = min(RELATIVE_LIMIT, tolerance * (1.0 + 1e-9))
    if (not converged or not math.isfinite(float(reported)) or not math.isfinite(relative)
            or relative > limit or not np.all(np.isfinite(solution.vec.FV().NumPy()))):
        raise RuntimeError(f"linear AMS-PCG did not satisfy the original-equation residual: "
                           f"{relative:.2e}, limit {limit:.2e}, iterations {iterations}")
    timings["residual_check_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    B_cf = ng.curl(solution)
    nu = ng.GridFunction(ng.L2(mesh, order=0), name="reluctivity")
    if nu.space.ndof != mesh.ne:
        raise RuntimeError("unexpected order-zero L2 layout")
    nu.vec.FV().NumPy()[:] = nu_values
    result = {"A": solution, "B_cf": B_cf, "H_cf": nu * B_cf, "nu": nu, "fes": fes}
    if observation is not None:
        result["observation_B_T"] = np.asarray([[float(c) for c in B_cf(mesh(*p))]
                                               for p in observation])
    timings["field_result_s"] = time.perf_counter() - t0
    result["stats"] = {
        "method": "linear lowest-order total-A beta-zero AMS-PCG", "converged": True,
        "cg_iterations": int(iterations), "linear_solves": 1, "cg_tolerance": tolerance,
        "final_residual_relative": relative, "reported_residual_relative": float(reported),
        "residual_norm": residual_norm, "rhs_norm": rhs_norm, "residual_limit": limit,
        "gradient_load_defect": defect / scale if scale > 0.0 else 0.0,
        "beta_zero": True, "gauge_epsilon": 0.0, "mixed_precision": bool(mixed_precision),
        "ndof": int(fes.ndof), "free_dofs": int(np.count_nonzero(free)),
        "phases_s": timings, "total_s": time.perf_counter() - started,
    }
    return result
