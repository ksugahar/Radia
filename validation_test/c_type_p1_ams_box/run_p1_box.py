"""First-order finite-box lane for the C-type dipole: reduced-A with AMS, and mixed Omega.

This is the accelerator-scale *solver* lane, deliberately separate from the
three-engine accuracy lane in ``validation_test/c_type_three_engine``: the
mesh is one finite air box of straight tetrahedra with the coil meshed, the
finite-element order is one, and the question is how fast and how accurately
the two H1/HCurl formulations reach the C-type gap field on such a mesh.

Engines on the same tetrahedra, the same CoilBuilder source and the same
monotone PCHIP B(H) law as the three-engine lane:

* ``reduced_a``: lowest-order Nedelec reduced vector potential
  (``B = B_s + curl A_r``, ``A_r x n = 0`` on the box).  The nonlinear loop is
  Newton on the element-constant flux density (tabulated approximation to the
  PCHIP inverse and its derivative, Armijo backtracking on the residual) or
  the damped/Anderson Picard update of the three-engine lane.  Each linear
  system is solved by Radia's compiled auxiliary-space Maxwell preconditioner
  (``radia.sparsesolv_ngsolve.HypreBasedAMSPreconditioner``, updated in place
  between iterations) with conjugate gradients to a true relative residual,
  warm-started from the current iterate; or by the METIS SPD PARDISO direct
  solve for a cross-check.
* ``mixed_omega``: first-order total/reduced scalar potential with the coil
  air as the reduced region, iron as the total region, the total-Hodge
  source split, natural ``B.n = 0`` on the box and the ``GND`` point gauge.

The reference for the accuracy column is a three-engine result of the Kelvin
lane (``--reference``): its consensus gap-core field is open-boundary, so the
deviation reported here bundles the box truncation with the discretisation.

Timing is wall time per phase inside one process with a fixed thread count;
the mesh is loaded once and the source is evaluated once.  Every run writes a
machine-readable JSON.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE_ROOT = HERE.parents[1] / "src"
if (SOURCE_ROOT / "radia" / "__init__.py").is_file():
    sys.path.insert(0, str(SOURCE_ROOT))

import ngsolve as ng  # noqa: E402
import numpy as np  # noqa: E402
import radia as rad  # noqa: E402
from ngsolve.krylovspace import CGSolver  # noqa: E402

from radia.kelvin_solver import (  # noqa: E402
    MixedOmegaPicardNotConverged,
    project_source_total_hodge,
    solve_magnetostatic_mixed_total_reduced_omega_kelvin,
    solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin,
)
from radia.picard_acceleration import ConstrainedAndersonAccelerator  # noqa: E402
from radia.mixed_omega_newton import (  # noqa: E402
    solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin,
    MixedOmegaNewtonNotConverged)
from radia.vector_potential_solver import (  # noqa: E402
    _build_nu_of_b_interpolator,
    direct_inverse_type,
)

MU0 = 4.0e-7 * math.pi
NU0 = 1.0 / MU0
GAUGE_EPSILON = 1.0e-6
THREE_ENGINE_RUNNER = HERE.parent / "c_type_three_engine" / "run_three_engine.py"
DEFAULT_BH = Path(rad.__file__).resolve().parent / "panels" / "samples" / "em_sample_bh.txt"
REQUIRED_MATERIALS = ("iron", "coil", "air")
CENTRE_INDEX = 40  # (x, y, z) = (0, 0, 0) in the shared 9 x 3 x 3 stencil


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_three_engine_adapters():
    """The C-type coil, observation stencil and comparison metrics, by path."""
    spec = importlib.util.spec_from_file_location("_radia_c_type_engines", THREE_ENGINE_RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {THREE_ENGINE_RUNNER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def progress(event: str, **values) -> None:
    print(json.dumps({"event": event, **values}, sort_keys=True, default=str), flush=True)


def check_mesh_contract(vol: Path, mesh: ng.Mesh) -> dict:
    contract_path = vol.with_suffix(".json")
    if not contract_path.is_file():
        raise FileNotFoundError(contract_path)
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if contract.get("vol_sha256") != sha256(vol):
        raise RuntimeError(f"{vol.name} differs from its contract {contract_path.name}")
    materials = set(mesh.GetMaterials())
    if not {"iron", "air"} <= materials <= set(REQUIRED_MATERIALS):
        raise RuntimeError(f"mesh materials {sorted(materials)} are not iron/air(/coil)")
    boundaries = set(mesh.GetBoundaries())
    for name in ("outer", "iron_air_interface"):
        if name not in boundaries:
            raise RuntimeError(f"boundary {name!r} missing from {sorted(boundaries)}")
    if "GND" not in set(mesh.GetBBBoundaries()):
        raise RuntimeError("GND vertex gauge missing")
    if mesh.GetCurveOrder() != 1:
        raise RuntimeError("this lane uses straight tetrahedra only")
    contract["contract_path"] = str(contract_path)
    return contract


def iron_elements_with_centroids(mesh: ng.Mesh) -> tuple[np.ndarray, np.ndarray]:
    """Element numbers and vertex-mean centroids of the iron tets, vectorised.

    Same element order and the same left-to-right vertex sum as a loop over
    ``mesh.Elements(VOL)``, hence bit-identical centroids.
    """
    elements = mesh.ngmesh.Elements3D().NumPy()
    iron_index = [index + 1 for index, name in enumerate(mesh.GetMaterials()) if name == "iron"]
    numbers = np.flatnonzero(np.isin(elements["index"], iron_index)).astype(np.int64)
    if numbers.size == 0:
        raise RuntimeError("mesh has no iron elements")
    if np.any(elements["np"][numbers] != 4):
        raise RuntimeError("this lane uses straight tetrahedra only")
    nodes = elements["nodes"][numbers][:, :4].astype(np.int64) - 1
    points = np.array([point.p for point in mesh.ngmesh.Points()], dtype=float)
    corners = points[nodes]
    # Python 3.12's float sum() is Neumaier-compensated; reproduce it so the
    # centroids stay bit-identical to the former per-element loop.
    total = np.zeros((len(numbers), 3))
    compensation = np.zeros_like(total)
    for k in range(4):
        term = corners[:, k]
        running = total + term
        compensation += np.where(np.abs(total) >= np.abs(term),
                                 (total - running) + term, (term - running) + total)
        total = running
    total = np.where(compensation != 0.0, total + compensation, total)
    return numbers, total / 4



class SoftIronLaw:
    """The shared monotone PCHIP B(H) law seen from the flux-density side.

    ``nu(|B|)`` approximates the production reduced-A Picard inverse
    (:func:`radia.vector_potential_solver._build_nu_of_b_interpolator`) by
    linear interpolation on a dense grid. This approximation also enters the
    residual and thus the converged field. ``dH/dB`` is a numerical derivative
    of the tabulated H values and steers the Newton direction.
    """

    def __init__(self, bh_table, *, samples: int = 4001):
        nu_of_b, b_saturation = _build_nu_of_b_interpolator(bh_table)
        self.b_saturation = float(b_saturation)
        self.nu_initial = float(nu_of_b(0.0))
        grid = np.concatenate(([0.0], np.geomspace(1.0e-7, 20.0 * self.b_saturation, samples)))
        nu = np.asarray([nu_of_b(value) for value in grid], dtype=float)
        h = nu * grid
        dhdb = np.gradient(h, grid)
        dhdb[0] = nu[0]
        if not np.all(np.isfinite(nu)) or np.any(nu <= 0.0) or np.any(dhdb <= 0.0):
            raise RuntimeError("the inverted B-H law is not positive and monotone")
        self._grid, self._nu, self._dhdb = grid, nu, dhdb

    def reluctivity(self, magnitude: np.ndarray) -> np.ndarray:
        return np.interp(magnitude, self._grid, self._nu)

    def differential_reluctivity(self, magnitude: np.ndarray) -> np.ndarray:
        return np.interp(magnitude, self._grid, self._dhdb)


class PiecewiseLinearIronLaw:
    """Explicit piecewise-linear H(B), with a vacuum-slope high-field tail.

    This is a different constitutive contract from SoftIronLaw's PCHIP B(H).
    Select it only when the input contract specifies this interpolation.
    """

    def __init__(self, bh_table):
        table = np.asarray(bh_table, dtype=float)
        if (table.ndim != 2 or table.shape[1] != 2 or len(table) < 2
                or not np.isfinite(table).all() or np.any(table[0] != 0)
                or np.any(np.diff(table, axis=0) <= 0)):
            raise ValueError("strictly increasing finite [H,B] table starting at [0,0] required")
        self._h, self._b = table[:, 0], table[:, 1]
        self._slopes = np.diff(self._h) / np.diff(self._b)
        self.nu_initial = float(self._slopes[0])
        self.b_saturation = float(self._b[-1])

    def reluctivity(self, magnitude):
        b = np.asarray(magnitude, dtype=float)
        if not np.isfinite(b).all() or np.any(b < 0):
            raise ValueError("finite nonnegative flux magnitude required")
        h = np.interp(b, self._b, self._h)
        h = np.where(b > self._b[-1], self._h[-1] + NU0*(b-self._b[-1]), h)
        return np.divide(h, b, out=np.full_like(b, self.nu_initial), where=b > 0)

    def differential_reluctivity(self, magnitude):
        b = np.asarray(magnitude, dtype=float)
        if not np.isfinite(b).all() or np.any(b < 0):
            raise ValueError("finite nonnegative flux magnitude required")
        index = np.clip(np.searchsorted(self._b, b, side="right")-1, 0, len(self._slopes)-1)
        return np.where(b >= self._b[-1], NU0, self._slopes[index])


class _TrueResidualCG(CGSolver):
    """Keep NGSolve's CG recurrence; stop on the original free-DOF equation."""

    def __init__(self, *, rhs, solution, free, tolerance, check_interval=1, **options):
        if isinstance(check_interval, bool) or int(check_interval) != check_interval or check_interval < 1:
            raise ValueError("check_interval must be a positive integer")
        self._check_interval = int(check_interval)
        super().__init__(tol=tolerance, **options)
        self._rhs = rhs
        self._solution = solution
        self._free = free
        self._relative_tolerance = tolerance
        self._rhs_norm = max(float(np.linalg.norm(rhs.FV().NumPy()[free])), 1e-30)
        self._true_residual = rhs.CreateVector()

    def CheckResidual(self, _preconditioned_residual):
        self.iterations += 1
        if (self.iterations != 1 and self.iterations < self.maxiter
                and (self.iterations - 1) % self._check_interval):
            return False
        self._true_residual.data = self._rhs - self.mat * self._solution
        relative = float(np.linalg.norm(self._true_residual.FV().NumPy()[self._free]) / self._rhs_norm)
        if not math.isfinite(relative):
            raise RuntimeError("non-finite true AMS-CG residual")
        self.residuals.append(relative)
        return relative <= self._relative_tolerance or self.iterations >= self.maxiter


class _GradientComplement(ng.BaseMatrix):
    """Diagnostic orthogonal P0 B P0; K and its residual remain unchanged."""

    def __init__(self, preconditioner, matrix, gradient, nodal_factor, free):
        super().__init__()
        self.preconditioner = preconditioner
        self.matrix = matrix
        self.gradient = gradient
        self.nodal_factor = nodal_factor
        self.free = free
        self.input = matrix.CreateRowVector()
        self.output = matrix.CreateColVector()

    def IsComplex(self): return False
    def Height(self): return self.matrix.height
    def Width(self): return self.matrix.width
    def CreateRowVector(self): return self.matrix.CreateRowVector()
    def CreateColVector(self): return self.matrix.CreateColVector()

    def project(self, vector):
        data = vector.copy()
        data[~self.free] = 0.0
        if self.nodal_factor is not None:
            data -= self.gradient @ self.nodal_factor.solve(self.gradient.T @ data)
        return data

    def Mult(self, x, y):
        self.input.FV().NumPy()[:] = self.project(x.FV().NumPy())
        self.preconditioner.Mult(self.input, self.output)
        y.FV().NumPy()[:] = self.project(self.output.FV().NumPy())


class ReducedAP1Box:
    """Reduced-A, order one, on a finite box; AMS or direct linear solves."""

    def __init__(self, mesh: ng.Mesh, source_cf, *, linear_solver: str,
                 cg_tolerance: float, cg_max_iterations: int, ams_num_smooth: int,
                 source_projection_order: int, ams_update_every: int = 1,
                 ic_shift: float = 1.05, gauge_epsilon: float = GAUGE_EPSILON,
                 ams_preconditioner_shift: float = 0.0,
                 ams_project_gradients: bool = False, ams_beta_zero: bool = False,
                 cg_check_interval: int = 1, ams_print_level: int = 0,
                 outer_boundary: str = "source_flux", zero_source: bool = False,
                 algebraic_residual: bool = True):
        """``source_cf`` is the vacuum source flux density B_s as a vector CF.

        ``zero_source=True`` declares ``source_cf`` identically zero (total-A):
        its iron projection and centroid values are then zero without being
        evaluated, which avoids a point search per iron element.

        ``linear_solver``: ``"ams"`` (compiled auxiliary-space Maxwell
        preconditioner + CG), ``"iccg"`` (Radia's compiled shifted incomplete
        Cholesky CG of the same sparsesolv module, factorised per system with
        ``ic_shift``), or ``"direct"`` (METIS SPD PARDISO).
        ``ams_update_every`` rebuilds the AMS hierarchy only every N-th linear
        system (1 = every system); a lagged preconditioner trades CG iterations
        for setup time and is recorded as such.
        """
        if linear_solver not in ("ams", "iccg", "direct"):
            raise ValueError("linear_solver must be 'ams', 'iccg' or 'direct'")
        if isinstance(cg_check_interval, bool) or int(cg_check_interval) != cg_check_interval or cg_check_interval < 1:
            raise ValueError("cg_check_interval must be a positive integer")
        if cg_check_interval != 1 and linear_solver != "ams":
            raise ValueError("cg_check_interval requires AMS")
        self.cg_check_interval = int(cg_check_interval)
        if ams_print_level not in (0, 1):
            raise ValueError("ams_print_level must be 0 or 1")
        self.ams_print_level = int(ams_print_level)
        if cg_tolerance <= 0.0 or cg_max_iterations < 1 or ams_num_smooth < 1:
            raise ValueError("cg_tolerance, cg_max_iterations and ams_num_smooth must be positive")
        if int(ams_update_every) < 1:
            raise ValueError("ams_update_every must be positive")
        if not ic_shift >= 1.0:
            raise ValueError("ic_shift must be at least 1")
        if not math.isfinite(gauge_epsilon) or gauge_epsilon < 0.0:
            raise ValueError("gauge_epsilon must be non-negative")
        if not math.isfinite(ams_preconditioner_shift) or ams_preconditioner_shift < 0.0:
            raise ValueError("ams_preconditioner_shift must be finite and non-negative")
        if ams_preconditioner_shift and (linear_solver != "ams" or gauge_epsilon != 0.0):
            raise ValueError("ams_preconditioner_shift requires AMS and gauge_epsilon=0")
        if ams_project_gradients and not ams_preconditioner_shift:
            raise ValueError("ams_project_gradients requires a positive AMS preconditioner shift")
        if ams_beta_zero and (linear_solver != "ams" or gauge_epsilon != 0.0
                              or ams_preconditioner_shift or ams_project_gradients):
            raise ValueError("ams_beta_zero requires unshifted AMS with gauge_epsilon=0 and no projection")
        if gauge_epsilon == 0.0 and not (linear_solver == "iccg" or
                (linear_solver == "ams" and (ams_preconditioner_shift > 0.0 or ams_beta_zero))):
            raise ValueError("gauge_epsilon=0 requires linear_solver='iccg', shifted AMS or beta-zero AMS")
        if outer_boundary not in ("source_flux", "natural_total"):
            raise ValueError("unknown outer boundary policy")
        self.outer_boundary = outer_boundary
        self.mesh = mesh
        self.linear_solver = linear_solver
        self.ic_shift = float(ic_shift)
        self.gauge_epsilon = float(gauge_epsilon)
        self.ams_preconditioner_shift = float(ams_preconditioner_shift)
        self._ams_shift_form = None
        self.ams_project_gradients = bool(ams_project_gradients)
        self.ams_beta_zero = bool(ams_beta_zero)
        self._gradient_projection = None
        self.cg_tolerance = float(cg_tolerance)
        self.cg_max_iterations = int(cg_max_iterations)
        self.ams_num_smooth = int(ams_num_smooth)
        self.ams_update_every = int(ams_update_every)
        self._ams_systems_seen = 0
        self.timing: dict[str, float] = {}
        started = time.perf_counter()
        self.fes = ng.HCurl(mesh, order=1,
                            dirichlet="outer" if outer_boundary == "source_flux" else "",
                            nograds=True)
        self.free = np.fromiter(self.fes.FreeDofs(), dtype=bool, count=self.fes.ndof)
        self.nu_space = ng.L2(mesh, order=0)
        self.nu_gf = ng.GridFunction(self.nu_space)
        self.iron_numbers, centroids = iron_elements_with_centroids(mesh)
        self.iron = self.iron_numbers  # one entry per iron element
        self.timing["spaces"] = time.perf_counter() - started
        # The exact Radia source, evaluated once on the iron and kept as a
        # discontinuous polynomial field; the reduced right-hand side only
        # ever integrates B_s over iron, where B_s is smooth.
        started = time.perf_counter()
        self.source_cf = source_cf
        self.zero_source = bool(zero_source)
        bs_space = ng.L2(mesh, order=int(source_projection_order), dim=3,
                         definedon=mesh.Materials("iron"))
        self.source_gf = ng.GridFunction(bs_space)
        if not self.zero_source:
            with ng.TaskManager():
                self.source_gf.Set(self.source_cf, definedon=mesh.Materials("iron"))
        self.timing["source_projection"] = time.perf_counter() - started
        self.source_projection_order = int(source_projection_order)
        # Element-constant helpers: the iron element numbers, the exact source at
        # their centroids (once), the element volumes, and order-zero fields for
        # the Newton rank-one coefficient and the element flux density.
        started = time.perf_counter()
        if self.zero_source:
            self.source_at_centroids = np.zeros((len(self.iron_numbers), 3))
        else:
            mips = mesh(centroids[:, 0], centroids[:, 1], centroids[:, 2])
            self.source_at_centroids = np.asarray(self.source_cf(mips), dtype=float).reshape(-1, 3)
        with ng.TaskManager():
            self.volumes = np.asarray(ng.Integrate(ng.CF(1.0), mesh, ng.VOL, element_wise=True),
                                      dtype=float)
        self._bvec_gfs = [ng.GridFunction(self.nu_space) for _ in range(3)]
        self.q_gf = ng.GridFunction(self.nu_space)
        self.timing["element_helpers"] = time.perf_counter() - started
        self._ams = None
        self._gradient = None
        # Residual and element flux by sparse products with the element curl
        # (built lazily, once); False keeps per-call form assembly.
        self.algebraic_residual = bool(algebraic_residual)
        self._algebra_cache = None
        self._residual_vector = None

    # ------------------------------------------------------------------ linear algebra
    def _ams_prepare(self, matrix):
        import radia.sparsesolv_ngsolve as ssn

        if self.ams_preconditioner_shift:
            if self._ams_shift_form is None:
                u, v = self.fes.TnT()
                self._ams_shift_form = ng.BilinearForm(self.fes, symmetric=True)
                self._ams_shift_form += ng.InnerProduct(u, v) * ng.dx
                with ng.TaskManager():
                    self._ams_shift_form.Assemble()
                self._ams_mass_values = self._ams_shift_form.mat.AsVector().FV().NumPy().copy()
            shifted = self._ams_shift_form.mat
            _, columns, offsets = matrix.CSR()
            _, mass_columns, mass_offsets = shifted.CSR()
            if not (np.array_equal(columns, mass_columns) and np.array_equal(offsets, mass_offsets)):
                raise RuntimeError("AMS mass shift requires matching sparse matrix ordering")
            shifted.AsVector().FV().NumPy()[:] = (matrix.AsVector().FV().NumPy()
                + self.ams_preconditioner_shift * NU0 * self._ams_mass_values)
            matrix = shifted  # Only hierarchy setup sees this matrix; CG keeps K.

        if self._gradient is None:
            gradient, h1 = self.fes.CreateGradient()
            if int(h1.ndof) != int(self.mesh.nv):
                raise RuntimeError("AMS needs one H1 DOF per vertex")
            xyz = np.asarray(self.mesh.ngmesh.Coordinates(), dtype=float)
            self._gradient = (gradient, xyz)
        gradient, xyz = self._gradient
        # Construction and Update must run outside TaskManager (sparsesolv contract).
        if self._ams is None:
            self._ams = ssn.HypreBasedAMSPreconditioner(
                mat=matrix, grad_mat=gradient, freedofs=self.fes.FreeDofs(),
                coord_x=xyz[:, 0].tolist(), coord_y=xyz[:, 1].tolist(),
                coord_z=xyz[:, 2].tolist(), cycle_type=1, print_level=self.ams_print_level,
                num_smooth=self.ams_num_smooth,
                **({"beta_zero": True} if self.ams_beta_zero else {}))
            self._ams_lagged = False
        elif self._ams_systems_seen % self.ams_update_every == 0:
            self._ams.Update(matrix)
            self._ams_lagged = False
        else:
            self._ams_lagged = True
        self._ams_systems_seen += 1
        if self.ams_project_gradients:
            if self._gradient_projection is None:
                from scipy.sparse import csr_matrix, diags
                from scipy.sparse.linalg import splu
                values, columns, offsets = gradient.CSR()
                discrete_gradient = csr_matrix((np.array(values), np.array(columns), np.array(offsets)),
                                               shape=(gradient.height, gradient.width))
                nodal_space = ng.H1(self.mesh, order=1, dirichlet="outer")
                nodal_free = np.array(list(nodal_space.FreeDofs()), dtype=bool)
                restricted = diags(self.free.astype(float)) @ discrete_gradient[:, nodal_free]
                factor = splu((restricted.T @ restricted).tocsc()) if restricted.shape[1] else None
                self._gradient_projection = restricted, factor
            return _GradientComplement(self._ams, matrix, *self._gradient_projection, self.free)
        return self._ams

    def _solve_linear(self, matrix, rhs, solution_vec, *, warm_start: bool) -> dict:
        """Solve ``matrix x = rhs`` for the free DOFs to a true relative residual."""
        record: dict[str, float | int | None] = {}
        rhs_norm = float(np.linalg.norm(rhs.FV().NumPy()[self.free]))
        if rhs_norm == 0.0:
            solution_vec[:] = 0.0
            return {"preconditioner_s": 0.0, "solve_s": 0.0, "cg_iterations": 0,
                    "relative_residual": 0.0}
        if self.linear_solver == "ams":
            t0 = time.perf_counter()
            pre = self._ams_prepare(matrix)
            record["preconditioner_s"] = time.perf_counter() - t0
            t0 = time.perf_counter()
            # Restarting from a preconditioned-norm stop with a tighter relative
            # tolerance can over-solve a singular system and amplify nullspace
            # roundoff. Check the original equation at each CG iteration instead.
            residual = rhs.CreateVector()
            solver = _TrueResidualCG(mat=matrix, pre=pre, maxiter=self.cg_max_iterations,
                rhs=rhs, solution=solution_vec, free=self.free,
                tolerance=self.cg_tolerance, check_interval=self.cg_check_interval, printrates=False)
            with ng.TaskManager():
                solver.Solve(rhs=rhs, sol=solution_vec, initialize=not warm_start)
                residual.data = rhs - matrix * solution_vec
            iterations = int(solver.iterations)
            true_relative = float(np.linalg.norm(residual.FV().NumPy()[self.free]) / rhs_norm)
            if not math.isfinite(true_relative) or true_relative > self.cg_tolerance:
                raise RuntimeError(
                    f"AMS-CG did not reach the true relative residual {self.cg_tolerance:.1e} "
                    f"(reached {true_relative:.2e} after {iterations} iterations)")
            record["solve_s"] = time.perf_counter() - t0
            record["cg_iterations"] = iterations
            record["cg_restarts"] = 0
            record["true_residual_checks"] = len(solver.residuals)
            record["preconditioner_lagged"] = bool(self._ams_lagged)
            record["relative_residual"] = true_relative
        elif self.linear_solver == "iccg":
            import radia.sparsesolv_ngsolve as ssn

            t0 = time.perf_counter()
            solver = ssn.SparseSolvSolver(
                matrix, method="ICCG", freedofs=self.fes.FreeDofs(), tol=self.cg_tolerance,
                maxiter=self.cg_max_iterations, shift=self.ic_shift, save_best_result=True,
                use_abmc=True)
            record["preconditioner_s"] = time.perf_counter() - t0
            t0 = time.perf_counter()
            if not warm_start:
                solution_vec[:] = 0.0
            iterations = 0
            residual = rhs.CreateVector()
            true_relative = float("inf")
            for attempt in range(6):
                with ng.TaskManager():
                    result = solver.Solve(rhs, solution_vec)
                    residual.data = rhs - matrix * solution_vec
                iterations += int(result.iterations)
                true_relative = float(np.linalg.norm(residual.FV().NumPy()[self.free]) / rhs_norm)
                if true_relative <= self.cg_tolerance:
                    break
                solver.tol = solver.tol * 0.1
            else:
                raise RuntimeError(
                    f"shifted ICCG did not reach the true relative residual {self.cg_tolerance:.1e} "
                    f"(reached {true_relative:.2e} after {iterations} iterations)")
            record["solve_s"] = time.perf_counter() - t0
            record["cg_iterations"] = iterations
            record["cg_restarts"] = attempt
            record["relative_residual"] = true_relative
        else:
            t0 = time.perf_counter()
            with ng.TaskManager():
                inverse = matrix.Inverse(self.fes.FreeDofs(), inverse=direct_inverse_type())
                solution_vec.data = inverse * rhs
                residual = rhs.CreateVector()
                residual.data = rhs - matrix * solution_vec
            record["preconditioner_s"] = 0.0
            record["solve_s"] = time.perf_counter() - t0
            record["cg_iterations"] = None
            record["relative_residual"] = float(
                np.linalg.norm(residual.FV().NumPy()[self.free]) / rhs_norm)
        if not math.isfinite(record["relative_residual"]):
            raise RuntimeError("non-finite linear residual")
        return record

    # ------------------------------------------------------------------ algebraic path
    def _algebra(self) -> dict:
        """Operators for residual and element flux without per-call form assembly.

        The curl of a lowest-order Nedelec function is constant per element, so
        the element curl is ``(B_k A) / w`` with ``B_k[e, i] = int_e curl(phi_i)_k
        psi_e`` (``psi_e`` the order-zero L2 function of element e, ``w_e = int_e
        psi_e``), and ``int_e nu curl A . curl v = nu_e vol_e curlA_e . curlv_e``.
        Built once; every later evaluation is sparse matrix-vector products.
        """
        if self._algebra_cache is not None:
            return self._algebra_cache
        from scipy.sparse import csr_matrix

        started = time.perf_counter()
        mesh = self.mesh
        scalar = ng.L2(mesh, order=0)
        if scalar.ndof != mesh.ne or scalar.GetDofNrs(ng.ElementId(ng.VOL, mesh.ne - 1))[0] != mesh.ne - 1:
            raise RuntimeError("unexpected order-0 L2 dof layout")
        u = self.fes.TrialFunction()
        w = scalar.TestFunction()
        weight_form = ng.LinearForm(scalar)
        weight_form += w * ng.dx
        curls, curls_t = [], []
        with ng.TaskManager():
            weight_form.Assemble()
            for k in range(3):
                form = ng.BilinearForm(trialspace=self.fes, testspace=scalar)
                form += ng.curl(u)[k] * w * ng.dx
                form.Assemble()
                values, columns, offsets = form.mat.CSR()
                matrix = csr_matrix((np.array(values), np.array(columns), np.array(offsets)),
                                    shape=(scalar.ndof, self.fes.ndof))
                curls.append(matrix)
                curls_t.append(matrix.T.tocsr())
        weights = weight_form.vec.FV().NumPy().copy()
        if np.any(weights <= 0.0):
            raise RuntimeError("non-positive order-zero L2 weights")
        algebra = {"curl": curls, "curl_t": curls_t, "weight": weights,
                   "scale": self.volumes / weights ** 2, "mass": None, "constant": None,
                   "source_mean": None}
        if self.gauge_epsilon > 0.0:
            a, b = self.fes.TnT()
            mass = ng.BilinearForm(self.fes, symmetric=False)
            mass += self.gauge_epsilon * NU0 * ng.InnerProduct(a, b) * ng.dx
            with ng.TaskManager():
                mass.Assemble()
            values, columns, offsets = mass.mat.CSR()
            algebra["mass"] = csr_matrix((np.array(values), np.array(columns), np.array(offsets)),
                                         shape=(self.fes.ndof, self.fes.ndof))
        if not self.zero_source:
            with ng.TaskManager():
                means = np.stack([np.asarray(ng.Integrate(self.source_gf[k], mesh, ng.VOL,
                                                          element_wise=True), dtype=float)
                                  for k in range(3)], axis=1)
            algebra["source_mean"] = means[self.iron_numbers] / self.volumes[self.iron_numbers, None]
        if self.outer_boundary == "natural_total":
            v = self.fes.TestFunction()
            boundary = ng.LinearForm(self.fes)
            boundary += -NU0 * ng.InnerProduct(ng.Cross(ng.specialcf.normal(3), self.source_cf),
                                                v.Trace()) * ng.ds("outer", bonus_intorder=4)
            with ng.TaskManager():
                boundary.Assemble()
            algebra["constant"] = boundary.vec.FV().NumPy().copy()
        algebra["build_s"] = time.perf_counter() - started
        self.timing["algebra_build"] = algebra["build_s"]
        self._algebra_cache = algebra
        return algebra

    def _element_curl(self, solution) -> np.ndarray:
        algebra = self._algebra()
        coefficients = solution.vec.FV().NumPy()
        return np.stack([matrix @ coefficients for matrix in algebra["curl"]], axis=1) \
            / algebra["weight"][:, None]

    # ------------------------------------------------------------------ material state
    def _element_flux(self, solution) -> tuple[np.ndarray, np.ndarray]:
        """Element-constant B in iron: exact source at the centroid + curl A."""
        if self.algebraic_residual:
            curl = self._element_curl(solution)
            b_iron = self.source_at_centroids + curl[self.iron_numbers]
            return b_iron, np.linalg.norm(b_iron, axis=1)
        with ng.TaskManager():
            curl = np.stack([
                np.asarray(ng.Integrate(ng.curl(solution)[k], self.mesh, ng.VOL,
                                        element_wise=True), dtype=float)
                for k in range(3)], axis=1) / self.volumes[:, None]
        b_iron = self.source_at_centroids + curl[self.iron_numbers]
        return b_iron, np.linalg.norm(b_iron, axis=1)

    def _set_material(self, nu_iron: np.ndarray, *, law: SoftIronLaw | None = None,
                      b_iron: np.ndarray | None = None, magnitude: np.ndarray | None = None) -> None:
        values = self.nu_gf.vec.FV().NumPy()
        values[:] = NU0
        values[self.iron_numbers] = nu_iron
        if law is None:
            return
        # Newton rank-one coefficient (dH/dB - nu)/|B|^2 and the element field.
        dhdb = law.differential_reluctivity(magnitude)
        safe = np.where(magnitude > 1.0e-9, magnitude, 1.0)
        q = np.where(magnitude > 1.0e-9, (dhdb - nu_iron) / (safe * safe), 0.0)
        q_values = self.q_gf.vec.FV().NumPy()
        q_values[:] = 0.0
        q_values[self.iron_numbers] = q
        for k in range(3):
            component = self._bvec_gfs[k].vec.FV().NumPy()
            component[:] = 0.0
            component[self.iron_numbers] = b_iron[:, k]

    def _picard_forms(self):
        u, v = self.fes.TnT()
        a = ng.BilinearForm(self.fes, symmetric=True)
        a += self.nu_gf * ng.InnerProduct(ng.curl(u), ng.curl(v)) * ng.dx
        if self.gauge_epsilon > 0.0:
            a += self.gauge_epsilon * NU0 * ng.InnerProduct(u, v) * ng.dx
        f = ng.LinearForm(self.fes)
        f += (NU0 - self.nu_gf) * ng.InnerProduct(self.source_gf, ng.curl(v)) * ng.dx("iron")
        if self.outer_boundary == "natural_total":
            f += NU0 * ng.InnerProduct(ng.Cross(ng.specialcf.normal(3), self.source_cf), v.Trace()) * ng.ds("outer", bonus_intorder=4)
        return a, f

    def _residual(self, solution) -> tuple:
        """Nonlinear residual R(A) = K(nu(A)) A - f(nu(A)) for the material set."""
        if self.algebraic_residual:
            algebra = self._algebra()
            coefficients = solution.vec.FV().NumPy()
            nu = self.nu_gf.vec.FV().NumPy()
            values = np.zeros(self.fes.ndof)
            iron = self.iron_numbers
            for k in range(3):
                load = nu * algebra["scale"] * (algebra["curl"][k] @ coefficients)
                if algebra["source_mean"] is not None:
                    load[iron] += ((nu[iron] - NU0) * (self.volumes[iron] / algebra["weight"][iron])
                                   * algebra["source_mean"][:, k])
                values += algebra["curl_t"][k] @ load
            if algebra["mass"] is not None:
                values += algebra["mass"] @ coefficients
            if algebra["constant"] is not None:
                values += algebra["constant"]
            if self._residual_vector is None:
                self._residual_vector = solution.vec.CreateVector()
            residual = self._residual_vector.CreateVector()
            residual.FV().NumPy()[:] = values
            return residual, float(np.linalg.norm(values[self.free]))
        v = self.fes.TestFunction()
        r = ng.LinearForm(self.fes)
        r += self.nu_gf * ng.InnerProduct(ng.curl(solution), ng.curl(v)) * ng.dx
        if not self.zero_source:  # identically zero otherwise; skipping it changes no bit
            r += (self.nu_gf - NU0) * ng.InnerProduct(self.source_gf, ng.curl(v)) * ng.dx("iron")
        if self.outer_boundary == "natural_total":
            r += -NU0 * ng.InnerProduct(ng.Cross(ng.specialcf.normal(3), self.source_cf), v.Trace()) * ng.ds("outer", bonus_intorder=4)
        if self.gauge_epsilon > 0.0:
            r += self.gauge_epsilon * NU0 * ng.InnerProduct(solution, v) * ng.dx
        with ng.TaskManager():
            r.Assemble()
        return r.vec, float(np.linalg.norm(r.vec.FV().NumPy()[self.free]))

    def _jacobian(self):
        u, v = self.fes.TnT()
        bvec = ng.CF(tuple(self._bvec_gfs))
        J = ng.BilinearForm(self.fes, symmetric=True)
        J += self.nu_gf * ng.InnerProduct(ng.curl(u), ng.curl(v)) * ng.dx
        J += self.q_gf * ng.InnerProduct(bvec, ng.curl(u)) * ng.InnerProduct(bvec, ng.curl(v)) * ng.dx("iron")
        if self.gauge_epsilon > 0.0:
            J += self.gauge_epsilon * NU0 * ng.InnerProduct(u, v) * ng.dx
        return J

    # ------------------------------------------------------------------ drivers
    def run_linear(self, mu_r: float, observation: np.ndarray) -> tuple:
        started = time.perf_counter()
        self._set_material(np.full(len(self.iron), NU0 / float(mu_r)))
        solution = ng.GridFunction(self.fes, name="A_reduced")
        solution.vec[:] = 0.0
        a, f = self._picard_forms()
        t0 = time.perf_counter()
        with ng.TaskManager():
            a.Assemble()
            f.Assemble()
        entry = {"iteration": 1, "assemble_s": time.perf_counter() - t0}
        entry.update(self._solve_linear(a.mat, f.vec, solution.vec, warm_start=False))
        field = self._observe(solution, observation)
        stats = {"method": "linear", "converged": True, "iterations": 1, "history": [entry]}
        return field, stats, time.perf_counter() - started

    def run_picard(self, law: SoftIronLaw, *, relax: float, anderson_depth: int,
                   tolerance: float, max_iterations: int, observation: np.ndarray) -> tuple:
        started = time.perf_counter()
        nu_current = np.full(len(self.iron), law.nu_initial)
        self._set_material(nu_current)
        solution = ng.GridFunction(self.fes, name="A_reduced")
        solution.vec[:] = 0.0
        accelerator = ConstrainedAndersonAccelerator(
            depth=int(anderson_depth), relaxation=float(relax),
            lower=min(NU0, law.nu_initial), upper=max(NU0, law.nu_initial), transform="log")
        history, b_previous, converged, final_change = [], None, False, None
        for iteration in range(1, int(max_iterations) + 1):
            entry = {"iteration": iteration}
            a, f = self._picard_forms()
            t0 = time.perf_counter()
            with ng.TaskManager():
                a.Assemble()
                f.Assemble()
            entry["assemble_s"] = time.perf_counter() - t0
            entry.update(self._solve_linear(a.mat, f.vec, solution.vec, warm_start=iteration > 1))
            t0 = time.perf_counter()
            b_iron, magnitude = self._element_flux(solution)
            nu_target = law.reluctivity(magnitude)
            accelerator.lower = min(accelerator.lower, float(np.min(nu_target)))
            accelerator.upper = max(accelerator.upper, float(np.max(nu_target)))
            nu_current = np.asarray(accelerator.step(nu_current, nu_target), dtype=float)
            self._set_material(nu_current)
            entry["material_update_s"] = time.perf_counter() - t0
            entry["nu_min"], entry["nu_max"] = float(np.min(nu_current)), float(np.max(nu_current))
            if b_previous is not None:
                final_change = float(np.max(np.abs(magnitude - b_previous)) / law.b_saturation)
                entry["relative_B_change"] = final_change
            history.append(entry)
            progress("picard", engine="reduced_a", iteration=iteration,
                     cg_iterations=entry.get("cg_iterations"), relative_B_change=final_change)
            if b_previous is not None and final_change < tolerance:
                converged = True
                break
            b_previous = magnitude.copy()
        field = self._observe(solution, observation)
        stats = {
            "method": "Picard", "converged": bool(converged), "iterations": len(history),
            "final_relative_change": final_change, "tolerance": float(tolerance),
            "maximum_iterations": int(max_iterations), "relaxation": float(relax),
            "anderson_depth": int(anderson_depth), "history": history,
            "maximum_linear_relative_residual": max(row["relative_residual"] for row in history),
            "anderson": accelerator.stats(),
        }
        return field, stats, time.perf_counter() - started

    def run_newton(self, law: SoftIronLaw, *, newton_tolerance: float, tolerance: float,
                   max_iterations: int, max_halvings: int, observation: np.ndarray,
                   inexact_linear: bool = False, linear_floor: bool = True) -> tuple:
        """Newton on the element-constant flux density with Armijo backtracking.

        ``linear_floor``: never ask the linear solve for an absolute residual
        below ``0.1 * newton_tolerance * |R_0|`` -- beyond the nonlinear target
        it buys nothing, and an ungauged singular system cannot deliver it
        (round-off leaves a gradient component near 1e-7 of |R_k|). The
        relative inner tolerance is then capped at 0.5; the nonlinear gates are
        unchanged.
        """
        if not 0.0 < float(newton_tolerance) < 1.0:
            raise ValueError("newton_tolerance must lie in (0, 1)")
        started = time.perf_counter()
        solution = ng.GridFunction(self.fes, name="A_reduced")
        solution.vec[:] = 0.0
        update = solution.vec.CreateVector()
        trial = ng.GridFunction(self.fes)
        b_iron, magnitude = self._element_flux(solution)
        self._set_material(law.reluctivity(magnitude), law=law, b_iron=b_iron, magnitude=magnitude)
        residual_vec, residual_norm = self._residual(solution)
        residual_0 = residual_norm
        history, converged, final_change = [], False, None
        if residual_0 == 0.0:
            return self._observe(solution, observation), {
                "method": "Newton", "converged": True, "iterations": 0,
                "inexact_linear": bool(inexact_linear),
                "final_relative_change": 0.0, "tolerance": float(tolerance),
                "newton_tolerance": float(newton_tolerance),
                "final_residual_relative": 0.0,
                "maximum_iterations": int(max_iterations), "history": [],
                "maximum_linear_relative_residual": 0.0,
            }, time.perf_counter() - started
        # Iterative solvers: one form for the whole run, its coefficients (nu, q,
        # element B) updated in place, reassembly reusing the sparsity graph.
        # The direct cross-check keeps a fresh form per step: the sparse direct
        # factorisation keeps state on the matrix and rejects a reassembled one.
        reuse_jacobian = self.linear_solver != "direct"
        J = self._jacobian() if reuse_jacobian else None
        for iteration in range(1, int(max_iterations) + 1):
            entry = {"iteration": iteration, "residual_relative": residual_norm / residual_0}
            t0 = time.perf_counter()
            if not reuse_jacobian:
                J = self._jacobian()
            with ng.TaskManager():
                J.Assemble()
            entry["assemble_s"] = time.perf_counter() - t0
            negative = residual_vec.CreateVector()
            negative.data = -1.0 * residual_vec
            update[:] = 0.0
            fixed_tolerance = self.cg_tolerance
            # Tighten inner solves with the nonlinear residual; outer gates stay fixed.
            if inexact_linear:
                self.cg_tolerance = max(fixed_tolerance, min(0.01, 0.1 * residual_norm / residual_0))
            if linear_floor:
                floor = 0.1 * float(newton_tolerance) * residual_0 / residual_norm
                entry["linear_floor"] = floor
                self.cg_tolerance = min(0.5, max(self.cg_tolerance, floor))
            entry["linear_tolerance"] = self.cg_tolerance
            try:
                entry.update(self._solve_linear(J.mat, negative, update, warm_start=False))
            finally:
                self.cg_tolerance = fixed_tolerance
            # Armijo backtracking on the true residual norm.
            t0 = time.perf_counter()
            alpha, accepted = 1.0, None
            for halving in range(int(max_halvings) + 1):
                trial.vec.data = solution.vec + alpha * update
                b_trial, magnitude_trial = self._element_flux(trial)
                self._set_material(law.reluctivity(magnitude_trial), law=law,
                                   b_iron=b_trial, magnitude=magnitude_trial)
                trial_vec, trial_norm = self._residual(trial)
                if trial_norm <= (1.0 - 1.0e-4 * alpha) * residual_norm:
                    accepted = (alpha, trial_vec, trial_norm, magnitude_trial)
                    break
                alpha *= 0.5
            if accepted is None:
                # Take the last (smallest) step anyway and let the residual
                # history show the stall rather than hide it.
                accepted = (alpha * 2.0, trial_vec, trial_norm, magnitude_trial)
                entry["line_search_failed"] = True
            alpha, residual_vec, new_norm, magnitude_new = accepted
            solution.vec.data = solution.vec + alpha * update
            entry["line_search_s"] = time.perf_counter() - t0
            entry["step_length"] = alpha
            entry["line_search_evaluations"] = halving + 1
            final_change = float(np.max(np.abs(magnitude_new - magnitude)) / law.b_saturation)
            entry["relative_B_change"] = final_change
            entry["residual_relative_after"] = new_norm / residual_0
            nu_now = self.nu_gf.vec.FV().NumPy()[self.iron_numbers]
            entry["nu_min"], entry["nu_max"] = float(np.min(nu_now)), float(np.max(nu_now))
            history.append(entry)
            magnitude, residual_norm = magnitude_new, new_norm
            progress("newton", engine="reduced_a", iteration=iteration, step=alpha,
                     cg_iterations=entry.get("cg_iterations"),
                     residual_relative=entry["residual_relative_after"],
                     relative_B_change=final_change)
            if residual_norm <= newton_tolerance * residual_0 and final_change <= tolerance:
                converged = True
                break
        field = self._observe(solution, observation)
        stats = {
            "method": "Newton", "converged": bool(converged), "iterations": len(history),
            "inexact_linear": bool(inexact_linear), "linear_floor": bool(linear_floor),
            "final_relative_change": final_change, "tolerance": float(tolerance),
            "newton_tolerance": float(newton_tolerance),
            "final_residual_relative": residual_norm / residual_0,
            "maximum_iterations": int(max_iterations), "history": history,
            "maximum_linear_relative_residual": max(row["relative_residual"] for row in history),
        }
        return field, stats, time.perf_counter() - started

    def _observe(self, solution, observation):
        total = self.source_cf + ng.curl(solution)
        return np.asarray([[float(c) for c in total(self.mesh(*map(float, p)))]
                           for p in observation], dtype=float)

    def describe(self) -> dict:
        return {
            "formulation": "HCurl reduced-A, order 1, nograds gauge",
            "boundary": ("n x H_total = 0 on the box; source boundary load included"
                         if self.outer_boundary == "natural_total" else
                         "A_r x n = 0 on the box (source flux passes through)"),
            "linear_solver": {
                "ams": "AMS(compiled, sparsesolv)+CG to a true relative residual, warm-started",
                "iccg": f"shifted IC(0) CG (compiled, sparsesolv, shift {self.ic_shift}, ABMC) "
                        "to a true relative residual, warm-started",
                "direct": direct_inverse_type(),
            }[self.linear_solver],
            "cg_relative_tolerance": self.cg_tolerance if self.linear_solver != "direct" else None,
            "cg_check_interval": self.cg_check_interval,
            "ams_print_level": self.ams_print_level,
            "ams_num_smooth": self.ams_num_smooth if self.linear_solver == "ams" else None,
            "ams_update_every": self.ams_update_every if self.linear_solver == "ams" else None,
            "ams_preconditioner_shift": self.ams_preconditioner_shift,
            "ams_project_gradients": self.ams_project_gradients,
            "ams_beta_zero": self.ams_beta_zero,
            "source": "exact Radia B_s projected once to L2 order "
                      f"{self.source_projection_order} on iron; exact at observation points",
            "ndof": int(self.fes.ndof),
            "gauge_epsilon": self.gauge_epsilon,
            "gauge": ("nograds space + mass regularisation" if self.gauge_epsilon > 0.0
                      else "nograds space only: singular, consistent right-hand side"),
            "setup_timing_s": self.timing,
        }


class TotalAP1Box(ReducedAP1Box):
    """Total-A with a meshed coil current and homogeneous outer tangential A.

    The caller supplies a charge-conserving current density in A/m^2 and
    verifies its cross-section current. No analytical source field is added.
    """

    def __init__(self, mesh, current_cf, **settings):
        if settings.get("outer_boundary", "source_flux") != "source_flux":
            raise ValueError("total-A requires homogeneous tangential A on outer")
        if "coil" not in mesh.GetMaterials():
            raise ValueError("total-A requires a meshed coil region")
        # zero_source=False keeps the generic evaluation path (diagnostic).
        zero_source = settings.pop("zero_source", True)
        super().__init__(mesh, ng.CF((0, 0, 0)), zero_source=zero_source, **settings)
        self.current_load = ng.LinearForm(self.fes)
        self.current_load += ng.InnerProduct(current_cf, self.fes.TestFunction()) * ng.dx("coil")
        with ng.TaskManager():
            self.current_load.Assemble()

    def _picard_forms(self):
        matrix, _ = super()._picard_forms()
        return matrix, self.current_load

    def _residual(self, solution):
        residual, _ = super()._residual(solution)
        residual.data -= self.current_load.vec
        return residual, float(np.linalg.norm(residual.FV().NumPy()[self.free]))

    def describe(self):
        result = super().describe()
        result.update(formulation="HCurl total-A, order 1, meshed coil J",
                      boundary="A_total x n = 0 on outer", source="volume integral J dot v on coil")
        return result


def coil_current_phi(mesh: ng.Mesh, coil_manifest: dict, cut_radius: float | None) -> dict:
    """A-phi current of the meshed racetrack: one cut across the first straight leg.

    The cut plane passes through the middle of the leg at +x with normal +y,
    the leg the CoilBuilder path starts on; the net current through it is the
    CoilBuilder current. The default in-plane radius lies halfway between the
    section's half diagonal and the in-plane distance to the opposite leg, so
    the plane's second crossing of the loop is excluded.
    """
    from radia.meshed_current import solve_closed_coil_current_phi

    centre = np.asarray(coil_manifest["centre_m"], dtype=float)
    leg_x = 0.5 * float(coil_manifest["straight_x_m"]) + float(coil_manifest["radius_m"])
    width, height = map(float, coil_manifest["cross_section_m"])
    half_diagonal = math.hypot(0.5 * width, 0.5 * height)
    opposite = 2.0 * leg_x - 0.5 * width
    if not half_diagonal < opposite:
        raise RuntimeError("coil legs too close for a single-leg cut")
    radius = 0.5 * (half_diagonal + opposite) if cut_radius is None else float(cut_radius)
    if not half_diagonal < radius < opposite:
        raise ValueError(f"cut radius must lie in ({half_diagonal:.4g}, {opposite:.4g}) m")
    origin = centre + np.array([leg_x, 0.0, 0.0])
    with ng.TaskManager():
        result = solve_closed_coil_current_phi(
            mesh, current_A=float(coil_manifest["current_A"]), cut_origin=origin,
            cut_normal=(0.0, 1.0, 0.0), cut_radius_m=radius, materials=("coil",))
    result["cut"] = {"origin_m": origin.tolist(), "normal": [0.0, 1.0, 0.0], "radius_m": radius}
    return result


def solve_mixed_omega_box(mesh: ng.Mesh, coil: int, material, *, nonlinear: bool,
                          relaxation: float, anderson_depth: int, tolerance: float,
                          max_iterations: int, observation: np.ndarray,
                          bonus_intorder: int, nonlinear_method: str = "picard",
                          order: int = 1) -> tuple:
    started = time.perf_counter()
    source_h = rad.RadiaField(coil, "h")
    timing = {}
    t0 = time.perf_counter()
    with ng.TaskManager():
        hodge = project_source_total_hodge(mesh, source_h, ("iron",), order=order,
                                           bonus_intorder=bonus_intorder)
    timing["source_hodge_projection_s"] = time.perf_counter() - t0
    reduced = tuple(name for name in ("air", "coil") if name in set(mesh.GetMaterials()))
    common = dict(
        reduced_materials=reduced, total_materials=("iron",),
        interface_boundary="iron_air_interface", order=order, dirichlet_bbbnd="GND",
        bonus_intorder=bonus_intorder, kelvin_mats=(), kelvin_match_exact=True,
        total_source_h=hodge["harmonic_field"], total_source_materials=("iron",),
    )
    t0 = time.perf_counter()
    stats = {}
    with ng.TaskManager():
        if nonlinear:
            solver = solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin
            iteration_options = dict(relaxation=float(relaxation), anderson_depth=int(anderson_depth),
                                     material_update_order=order - 1 if order > 1 else None)
            if nonlinear_method == "newton":
                if anderson_depth != 0:
                    raise ValueError("Newton does not use Anderson mixing")
                solver = solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin
                iteration_options = dict(progress_callback=lambda row: progress("newton", engine="mixed_omega", **row))
            elif nonlinear_method != "picard":
                raise ValueError("nonlinear_method must be picard or newton")
            try:
                result = solver(
                    mesh, source_h, hodge["potential"], 1.0, (10.0, 0.0, 0.0),
                    bh_table=material, nonlinear_materials=("iron",),
                    tolerance=float(tolerance), max_iterations=int(max_iterations),
                    **iteration_options,
                    observation_points=observation, **common)
            except (MixedOmegaPicardNotConverged, MixedOmegaNewtonNotConverged) as exc:
                # The carried state has the iteration history but no field:
                # report the failure with its history and a NaN field.
                stats = dict(exc.state["nonlinear_stats"])
                stats["converged"] = False
                stats["error"] = str(exc)[:500]
                timing["solve_s"] = time.perf_counter() - t0
                field = np.full((len(observation), 3), float("nan"))
                description = {
                    "formulation": f"H1 mixed total/reduced Omega, order {order}",
                    "boundary": "natural B.n = 0 on the box; GND vertex gauge",
                    "iron_relative_harmonic_norm": float(hodge["relative_harmonic_norm"]),
                    "linear_solver": "PARDISO (symmetric indefinite saddle point)",
                    "timing_s": timing, "ndof": None,
                }
                return field, stats, description, time.perf_counter() - started
            stats = dict(result["nonlinear_stats"])
        else:
            result = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
                mesh, source_h, hodge["potential"], 1.0, (10.0, 0.0, 0.0),
                mu_r_by_material={"iron": float(material)}, **common)
            stats = {"converged": True, "iterations": 1}
    timing["solve_s"] = time.perf_counter() - t0
    field = np.asarray([[float(c) for c in result["B_cf"](mesh(*map(float, p)))]
                        for p in observation], dtype=float)
    description = {
        "formulation": f"H1 mixed total/reduced Omega, order {order}",
        "boundary": "natural B.n = 0 on the box; GND vertex gauge",
        "source": "exact Radia H_s in the reduced air/coil; total-Hodge split in iron",
        "iron_relative_harmonic_norm": float(hodge["relative_harmonic_norm"]),
        "ndof": int(result["fes"].ndof),
        "linear_solver": "PARDISO (symmetric indefinite saddle point)",
        "phase_timings_seconds": result.get("phase_timings_seconds"),
        "timing_s": timing,
    }
    return field, stats, description, time.perf_counter() - started


def compare(points, fields: dict, reference: dict | None, core_half_length: float, adapters):
    core = np.abs(points[:, 0]) <= core_half_length + 1e-14
    projected = {name: adapters.median_plane_projection(points, value)
                 for name, value in fields.items()}
    out = {
        "pairwise_median_projected_gap_core": adapters.pairwise_metrics(projected, core),
        "pairwise_raw_full_tube": adapters.pairwise_metrics(fields, np.ones(len(points), bool)),
        "centre_field_T": {name: value[CENTRE_INDEX].tolist() for name, value in fields.items()},
    }
    if reference is not None:
        ref_points = np.asarray(reference["observation_points_m"], dtype=float)
        if not np.allclose(ref_points, points, atol=1e-12):
            raise RuntimeError("reference observation stencil differs")
        ref_fields = {f"reference_{name}": np.asarray(value, dtype=float)
                      for name, value in reference["fields_T"].items()}
        ref_projected = {name: adapters.median_plane_projection(points, value)
                         for name, value in ref_fields.items()}
        versus = {}
        for ours, value in projected.items():
            for theirs, ref in ref_projected.items():
                versus[f"{ours}__vs__{theirs}"] = {
                    "gap_core_relative_rms": adapters.relative_rms(ref[core], value[core]),
                    "centre_relative_difference": float(
                        np.linalg.norm(value[CENTRE_INDEX] - ref[CENTRE_INDEX])
                        / np.linalg.norm(ref[CENTRE_INDEX])),
                }
        out["versus_reference"] = versus
        out["reference_centre_field_T"] = {
            name: v[CENTRE_INDEX].tolist() for name, v in ref_fields.items()}
        out["reference_file"] = reference["_path"]
        out["reference_sha256"] = reference["_sha256"]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ams-print-level", type=int, choices=(0, 1), default=0,
                        help="Native AMS setup and sampled application timing diagnostics")
    parser.add_argument("--cg-check-interval", type=int, default=1,
                        help="AMS true residual check interval; final check is mandatory")
    parser.add_argument("--inexact-linear", action="store_true",
                        help="Adapt Newton inner tolerance without relaxing final convergence gates")
    parser.add_argument("--legacy-residual", action="store_true",
                        help="assemble the Newton residual and element flux as forms on every call "
                             "(default: sparse products with the element curl, built once)")
    parser.add_argument("--no-linear-floor", action="store_true",
                        help="Newton: drop the absolute inner-solve floor 0.1*newton_tolerance*|R_0| "
                             "(needed by ungauged systems at tight rules; on by default)")
    parser.add_argument("--vol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bh-table", type=Path, default=DEFAULT_BH)
    parser.add_argument("--mode", choices=("linear", "nonlinear"), default="nonlinear")
    parser.add_argument("--mu-r", type=float, default=1000.0)
    parser.add_argument("--engines", default="reduced_a,mixed_omega",
                        help="comma list of reduced_a, total_a (meshed coil, A-phi current), mixed_omega")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--reduced-a-solver", choices=("ams", "iccg", "direct"), default="ams")
    parser.add_argument("--ic-shift", type=float, default=1.05,
                        help="shift of the incomplete Cholesky factorisation (iccg)")
    parser.add_argument("--coil-cut-radius", type=float, default=None,
                        help="total_a: in-plane radius of the A-phi cut (default between the "
                             "section half diagonal and the opposite leg)")
    parser.add_argument("--arc-rel-tol", type=float, default=None,
                        help="relative tolerance of the arc-current source quadrature "
                             "(Radia 'PrcArc', [1e-12, 1e-3]); default keeps Radia's 1e-9")
    parser.add_argument("--gauge-epsilon", type=float, default=GAUGE_EPSILON,
                        help="operator mass regularisation; 0 with iccg, beta-zero or shifted AMS")
    parser.add_argument("--ams-preconditioner-shift", type=float, default=0.0,
                        help="AMS-only mass shift sigma*nu_0*M; requires --gauge-epsilon 0")
    parser.add_argument("--ams-beta-zero", action="store_true",
                        help="Pure curl-curl AMS without gradient correction; requires gauge-epsilon 0")
    parser.add_argument("--ams-project-gradients", action="store_true",
                        help="diagnostic P0 B_shift P0 with a cached nodal direct solve")
    parser.add_argument("--nonlinear-method", choices=("newton", "picard"), default="newton")
    parser.add_argument("--cg-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--cg-max-iterations", type=int, default=2000)
    parser.add_argument("--ams-num-smooth", type=int, default=1)
    parser.add_argument("--ams-update-every", type=int, default=1,
                        help="rebuild the AMS hierarchy every N-th linear system")
    parser.add_argument("--source-projection-order", type=int, default=2)
    parser.add_argument("--outer-boundary", choices=("source_flux", "natural_total"), default="source_flux")
    parser.add_argument("--relax", type=float, default=0.3)
    parser.add_argument("--anderson-depth", type=int, default=0)
    parser.add_argument("--newton-tolerance", type=float, default=1.0e-6,
                        help="relative nonlinear residual at which Newton stops")
    parser.add_argument("--line-search-max-halvings", type=int, default=6)
    parser.add_argument("--omega-relax", type=float, default=0.3)
    parser.add_argument("--omega-nonlinear-method", choices=("picard", "newton"), default="picard")
    parser.add_argument("--omega-order", type=int, choices=(1, 2), default=1)
    parser.add_argument("--omega-anderson-depth", type=int, default=0)
    parser.add_argument("--omega-bonus-intorder", type=int, default=4)
    parser.add_argument("--tolerance", type=float, default=2.0e-5,
                        help="max |dB| / B_sat between iterations at which a loop is converged")
    parser.add_argument("--max-iterations", type=int, default=80)
    parser.add_argument("--gap-core-half-length", type=float, default=0.010)
    parser.add_argument("--reference", type=Path, default=None,
                        help="three-engine result JSON of the Kelvin lane")
    options = parser.parse_args()

    try:
        output = run(options)
    except Exception as exc:
        output = {
            "schema": "radia.validation.c-type-p1-box-bench.v1",
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "machine": platform.node(), "python": sys.version,
            "ngsolve": ng.__version__, "completed": False, "passed": False,
            "error": {"type": type(exc).__name__, "message": str(exc)},
            "options": {key: str(value) if isinstance(value, Path) else value
                        for key, value in vars(options).items()},
        }
        options.output.parent.mkdir(parents=True, exist_ok=True)
        options.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
        raise
    if not output["passed"]:
        raise SystemExit(1)


def run(options) -> dict:

    if options.threads > 0:
        ng.SetNumThreads(options.threads)
    engines = [name.strip() for name in options.engines.split(",") if name.strip()]
    if not engines or len(set(engines)) != len(engines):
        raise ValueError("engines must be a nonempty list without duplicates")
    for name in engines:
        if name not in ("reduced_a", "total_a", "mixed_omega"):
            raise ValueError(f"unknown engine {name!r}")
    adapters = load_three_engine_adapters()
    nonlinear = options.mode == "nonlinear"
    material = (np.loadtxt(options.bh_table, dtype=float)[:, :2].tolist()
                if nonlinear else options.mu_r)

    t0 = time.perf_counter()
    mesh = ng.Mesh(str(options.vol.resolve()))
    mesh_load_s = time.perf_counter() - t0
    contract = check_mesh_contract(options.vol.resolve(), mesh)
    rad.UtiDelAll()
    coil, coil_manifest = adapters.build_coil()
    if options.arc_rel_tol is not None:
        # Global Radia setting; applies to every later source evaluation of this run.
        rad.FldCmpPrc(f"PrcArc->{float(options.arc_rel_tol):.17g}")
    coil_manifest = dict(coil_manifest, arc_rel_tol=options.arc_rel_tol)
    points = adapters.observation_points()
    for point in points:
        if not mesh(*map(float, point)):
            raise RuntimeError(f"observation point {point} outside the mesh")

    reference = None
    if options.reference is not None:
        reference = json.loads(options.reference.read_text(encoding="utf-8"))
        reference["_path"] = str(options.reference.resolve())
        reference["_sha256"] = sha256(options.reference)

    fields: dict[str, np.ndarray] = {}
    diagnostics: dict[str, dict] = {}
    potential_settings = dict(
        linear_solver=options.reduced_a_solver,
        cg_tolerance=options.cg_tolerance, cg_max_iterations=options.cg_max_iterations,
        cg_check_interval=options.cg_check_interval,
        ams_print_level=options.ams_print_level,
        ams_num_smooth=options.ams_num_smooth,
        source_projection_order=options.source_projection_order,
        outer_boundary=options.outer_boundary,
        ams_update_every=options.ams_update_every, ic_shift=options.ic_shift,
        gauge_epsilon=options.gauge_epsilon,
        ams_preconditioner_shift=options.ams_preconditioner_shift,
        ams_project_gradients=options.ams_project_gradients,
        ams_beta_zero=options.ams_beta_zero,
        algebraic_residual=not options.legacy_residual)

    def run_potential_engine(name, engine, extra=None):
        if not nonlinear:
            field, stats, runtime = engine.run_linear(options.mu_r, points)
        elif options.nonlinear_method == "newton":
            field, stats, runtime = engine.run_newton(
                SoftIronLaw(material), newton_tolerance=options.newton_tolerance,
                tolerance=options.tolerance, max_iterations=options.max_iterations,
                max_halvings=options.line_search_max_halvings, observation=points,
                inexact_linear=options.inexact_linear, linear_floor=not options.no_linear_floor)
        else:
            field, stats, runtime = engine.run_picard(
                SoftIronLaw(material), relax=options.relax,
                anderson_depth=options.anderson_depth, tolerance=options.tolerance,
                max_iterations=options.max_iterations, observation=points)
        fields[name] = field
        totals: dict[str, float] = {}
        for row in stats["history"]:
            for key in ("assemble_s", "preconditioner_s", "solve_s", "material_update_s",
                        "line_search_s"):
                if key in row:
                    totals[key] = totals.get(key, 0.0) + float(row[key])
        diagnostics[name] = {
            **engine.describe(), **(extra or {}), "nonlinear": nonlinear, "nonlinear_stats": stats,
            "runtime_s": runtime, "mesh_elements": int(mesh.ne), "mesh_vertices": int(mesh.nv),
            "phase_totals_s": totals,
            "total_cg_iterations": sum(int(row["cg_iterations"] or 0) for row in stats["history"]),
        }
        progress("engine_complete", engine=name, runtime_s=runtime,
                 converged=stats["converged"], iterations=stats["iterations"])

    if "reduced_a" in engines:
        progress("engine_start", engine="reduced_a", mesh_elements=int(mesh.ne))
        engine = ReducedAP1Box(mesh, rad.RadiaField(coil, "b"), **potential_settings)
        run_potential_engine("reduced_a", engine)
    if "total_a" in engines:
        progress("engine_start", engine="total_a", mesh_elements=int(mesh.ne))
        current = coil_current_phi(mesh, coil_manifest, options.coil_cut_radius)
        t0 = time.perf_counter()
        engine = TotalAP1Box(mesh, current["current"], **potential_settings)
        setup_s = time.perf_counter() - t0
        run_potential_engine("total_a", engine, {
            "coil_current": current["stats"], "coil_current_cut": current["cut"],
            "engine_setup_s": setup_s})
    if "mixed_omega" in engines:
        progress("engine_start", engine="mixed_omega", mesh_elements=int(mesh.ne))
        field, stats, description, runtime = solve_mixed_omega_box(
            mesh, coil, material, nonlinear=nonlinear, relaxation=options.omega_relax,
            anderson_depth=options.omega_anderson_depth, tolerance=options.tolerance,
            max_iterations=options.max_iterations, observation=points,
            bonus_intorder=options.omega_bonus_intorder,
            nonlinear_method=options.omega_nonlinear_method, order=options.omega_order)
        fields["mixed_omega"] = field
        diagnostics["mixed_omega"] = {**description, "nonlinear": nonlinear,
                                      "nonlinear_stats": stats, "runtime_s": runtime,
                                      "mesh_elements": int(mesh.ne), "mesh_vertices": int(mesh.nv)}
        progress("engine_complete", engine="mixed_omega", runtime_s=runtime,
                 converged=stats.get("converged"), iterations=stats.get("iterations"))

    comparison = compare(points, fields, reference, options.gap_core_half_length, adapters)
    output = {
        "schema": "radia.validation.c-type-p1-box-bench.v1",
        "completed": True,
        "passed": all(row["nonlinear_stats"].get("converged") is True
                      for row in diagnostics.values()),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "machine": platform.node(),
        "python": sys.version,
        "ngsolve": ng.__version__,
        "threads": int(options.threads),
        "mode": options.mode,
        "mesh": {"vol": str(options.vol.resolve()), "vol_sha256": contract["vol_sha256"],
                 "elements": int(mesh.ne), "vertices": int(mesh.nv),
                 "load_seconds": mesh_load_s, "contract": contract},
        "bh_table": None if not nonlinear else str(options.bh_table.resolve()),
        "bh_table_sha256": None if not nonlinear else sha256(options.bh_table),
        "mu_r": None if nonlinear else options.mu_r,
        "coil": coil_manifest,
        "observation_points_m": points.tolist(),
        "engines": diagnostics,
        "fields_T": {name: value.tolist() for name, value in fields.items()},
        "comparison": comparison,
        "gap_core_half_length_m": options.gap_core_half_length,
        "options": {key: (str(value) if isinstance(value, Path) else value)
                    for key, value in vars(options).items()},
    }
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(output, indent=2, default=str) + "\n", encoding="utf-8")
    progress("complete", output=str(options.output),
             runtimes={name: row["runtime_s"] for name, row in diagnostics.items()})
    return output


if __name__ == "__main__":
    main()
