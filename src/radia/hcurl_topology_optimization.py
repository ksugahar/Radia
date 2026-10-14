"""Analytic HCurl Eddy-Bubble geometry adjoints and LP updates."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.optimize import linprog

from .vim._hcurl_tet_interaction import (
    HCurlCellVolumeInteraction,
    HCurlTetVolumeInteraction,
    HCurlHMatrixOperator,
    SampleNgsolveHCurlCellSubtetVelocities,
)
from .vim._vim import _csr


@dataclass(frozen=True)
class HCurlResistanceLinearization:
    matrix: np.ndarray
    jacobian: np.ndarray


@dataclass(frozen=True)
class HCurlJouleAdjointLinearization:
    state: np.ndarray
    adjoint: np.ndarray
    objective: float
    gradient: np.ndarray
    resistance: np.ndarray
    resistance_jacobian: np.ndarray
    inductance: object


@dataclass(frozen=True)
class HCurlJouleLPUpdate:
    design: np.ndarray
    delta: np.ndarray
    objective_gradient: np.ndarray
    status: str


@dataclass(frozen=True)
class HCurlJouleLoadCase:
    frequency_hz: float
    rhs: np.ndarray
    weight: float = 1.0
    rhs_jacobian: np.ndarray | None = None


@dataclass(frozen=True)
class HCurlMultiFrequencyJouleLinearization:
    cases: tuple[HCurlJouleAdjointLinearization, ...]
    objective: float
    gradient: np.ndarray
    weights: np.ndarray


@dataclass(frozen=True)
class HCurlSheetLPStep:
    linearization: HCurlMultiFrequencyJouleLinearization
    update: object


@dataclass(frozen=True)
class HCurlConductivityInterpolation:
    solid: float
    void: float
    power: float = 3.0

    def evaluate(self, activation):
        rho = np.asarray(activation, dtype=float)
        if np.any(~np.isfinite(rho)) or np.any((rho < 0) | (rho > 1)):
            raise ValueError("activation must be finite and lie in [0,1]")
        solid, void, power = float(self.solid), float(self.void), float(self.power)
        if not (np.isfinite(solid) and np.isfinite(void) and solid > void > 0):
            raise ValueError("conductivities must satisfy solid > void > 0")
        if not np.isfinite(power) or power < 1:
            raise ValueError("conductivity power must be at least one")
        sigma = void + (solid-void)*rho**power
        dsigma = (solid-void)*power*rho**(power-1)
        return sigma, dsigma


def linearize_hcurl_multifrequency_joule_loss(
    *, inductance, resistance, resistance_jacobian, cell_vertex_velocities,
    load_cases,
) -> HCurlMultiFrequencyJouleLinearization:
    """Linearize preassembled HCurl Joule systems for multiple load cases.

    This is the UI-neutral numerical contract shared with the native MATLAB
    MEX entry point.  NGSolve-specific assembly belongs in the ``_from_ngsolve``
    adapter below; this function owns the complex state/adjoint convention.
    """
    if not isinstance(inductance, HCurlHMatrixOperator):
        raise TypeError("inductance must be HCurlHMatrixOperator")
    n = inductance.mode_count
    r_matrix = np.asarray(resistance, dtype=float)
    if r_matrix.shape != (n, n) or np.any(~np.isfinite(r_matrix)):
        raise ValueError(f"resistance must be finite with shape ({n},{n})")
    dr = np.asarray(resistance_jacobian, dtype=float)
    if dr.ndim != 3 or dr.shape[1:] != (n, n) or np.any(~np.isfinite(dr)):
        raise ValueError(f"resistance_jacobian must be finite with shape (q,{n},{n})")
    velocity = np.asarray(cell_vertex_velocities, dtype=float)
    if velocity.ndim == 3:
        velocity = velocity[None]
    if velocity.ndim != 4 or velocity.shape[0] != dr.shape[0]:
        raise ValueError("cell_vertex_velocities must contain one entry per derivative")
    cases = tuple(load_cases)
    if not cases:
        raise ValueError("load_cases must not be empty")

    l_dense = inductance.to_dense()
    results = []
    weights = np.empty(len(cases))
    total_objective = 0.0
    total_gradient = np.zeros(dr.shape[0])
    for index, case in enumerate(cases):
        if not isinstance(case, HCurlJouleLoadCase):
            raise TypeError("load_cases entries must be HCurlJouleLoadCase")
        weight = float(case.weight)
        omega = 2.0*np.pi*float(case.frequency_hz)
        if not np.isfinite(weight) or weight < 0:
            raise ValueError("load-case weights must be finite and non-negative")
        if not np.isfinite(omega) or omega <= 0:
            raise ValueError("frequency_hz must be positive")
        load = np.asarray(case.rhs, dtype=complex).reshape(-1)
        if load.shape != (n,):
            raise ValueError(f"rhs must have shape ({n},)")
        if case.rhs_jacobian is None:
            db = np.zeros((dr.shape[0], n), dtype=complex)
        else:
            db = np.asarray(case.rhs_jacobian, dtype=complex)
            if db.shape != (dr.shape[0], n):
                raise ValueError(
                    f"rhs_jacobian must have shape ({dr.shape[0]},{n})"
                )

        system = r_matrix + 1j*omega*l_dense
        state = np.linalg.solve(system, load)
        r_state = r_matrix@state
        objective = 0.5*float(np.real(np.vdot(state, r_state)))
        adjoint = np.linalg.solve(system.conj().T, r_state)
        dl = np.asarray(inductance.directional_contractions(
            velocity, adjoint, state
        ))
        gradient = np.empty(dr.shape[0])
        for direction, dr_matrix in enumerate(dr):
            dr_state = dr_matrix@state
            gradient[direction] = 0.5*np.real(np.vdot(state, dr_state)) + np.real(
                np.vdot(adjoint, db[direction]-dr_state)
                - 1j*omega*dl[direction]
            )
        results.append(HCurlJouleAdjointLinearization(
            state, adjoint, objective, gradient, r_matrix, dr, inductance,
        ))
        weights[index] = weight
        total_objective += weight*objective
        total_gradient += weight*gradient
    return HCurlMultiFrequencyJouleLinearization(
        tuple(results), float(total_objective), total_gradient, weights
    )


def linearize_hcurl_multifrequency_activation_joule_loss(
    *, inductance, cell_curl_grams, activation, load_cases, conductivity,
    inductance_power=1.0,
) -> HCurlMultiFrequencyJouleLinearization:
    """Linearize a preassembled cellwise conductivity/inductance topology."""
    if not isinstance(inductance, HCurlHMatrixOperator):
        raise TypeError("inductance must be HCurlHMatrixOperator")
    if not isinstance(conductivity, HCurlConductivityInterpolation):
        raise TypeError("conductivity must be HCurlConductivityInterpolation")
    n = inductance.mode_count
    rho = np.asarray(activation, dtype=float).reshape(-1)
    grams = np.asarray(cell_curl_grams, dtype=float)
    if grams.shape != (len(rho), n, n) or np.any(~np.isfinite(grams)):
        raise ValueError(
            f"cell_curl_grams must be finite with shape ({len(rho)},{n},{n})"
        )
    sigma, dsigma = conductivity.evaluate(rho)
    inverse_sigma = 1.0/sigma
    inverse_sigma_derivative = -dsigma/sigma**2
    r_matrix = np.einsum("e,eij->ij", inverse_sigma, grams)
    l_dense = inductance.activation_to_dense(rho, power=inductance_power)
    cases = tuple(load_cases)
    if not cases:
        raise ValueError("load_cases must not be empty")

    results = []
    weights = np.empty(len(cases))
    total_objective = 0.0
    total_gradient = np.zeros(len(rho))
    for index, case in enumerate(cases):
        if not isinstance(case, HCurlJouleLoadCase):
            raise TypeError("load_cases entries must be HCurlJouleLoadCase")
        weight = float(case.weight)
        omega = 2.0*np.pi*float(case.frequency_hz)
        if not np.isfinite(weight) or weight < 0:
            raise ValueError("load-case weights must be finite and non-negative")
        if not np.isfinite(omega) or omega <= 0:
            raise ValueError("frequency_hz must be positive")
        load = np.asarray(case.rhs, dtype=complex).reshape(-1)
        if load.shape != (n,):
            raise ValueError(f"rhs must have shape ({n},)")
        if case.rhs_jacobian is None:
            db = np.zeros((len(rho), n), dtype=complex)
        else:
            db = np.asarray(case.rhs_jacobian, dtype=complex)
            if db.shape != (len(rho), n):
                raise ValueError(f"rhs_jacobian must have shape ({len(rho)},{n})")

        system = r_matrix + 1j*omega*l_dense
        state = np.linalg.solve(system, load)
        r_state = r_matrix@state
        objective = 0.5*float(np.real(np.vdot(state, r_state)))
        adjoint = np.linalg.solve(system.conj().T, r_state)
        gram_state = np.einsum("eij,j->ei", grams, state)
        direct = 0.5*np.real(np.einsum("i,ei->e", np.conj(state), gram_state))
        implicit_r = -np.real(
            np.einsum("i,ei->e", np.conj(adjoint), gram_state)
        )
        implicit_b = np.real(db@np.conj(adjoint))
        dl = np.asarray(inductance.activation_contractions(
            rho, adjoint, state, power=inductance_power
        ))
        gradient = (
            (direct+implicit_r)*inverse_sigma_derivative + implicit_b
            - np.real(1j*omega*dl)
        )
        results.append(HCurlJouleAdjointLinearization(
            state, adjoint, objective, gradient, r_matrix,
            np.empty((0,)+r_matrix.shape), inductance,
        ))
        weights[index] = weight
        total_objective += weight*objective
        total_gradient += weight*gradient
    return HCurlMultiFrequencyJouleLinearization(
        tuple(results), float(total_objective), total_gradient, weights
    )


def _ngsolve_cellwise_coefficient(mesh, values, materials=None):
    import ngsolve as ng
    selected = None if materials is None else {str(x) for x in (
        (materials,) if isinstance(materials, str) else materials
    )}
    elements = [element for element in mesh.Elements(ng.VOL)
        if selected is None or str(element.mat) in selected]
    data = np.asarray(values, dtype=float).reshape(-1)
    if data.shape != (len(elements),):
        raise ValueError(f"cell values must have shape ({len(elements)},)")
    space = ng.L2(mesh, order=0)
    field = ng.GridFunction(space)
    field.vec[:] = 0.0
    for value, element in zip(data, elements):
        dofs = space.GetDofNrs(element)
        if len(dofs) != 1:
            raise RuntimeError("NGSolve L2 order-0 cell contract changed")
        field.vec[dofs[0]] = float(value)
    return field, tuple(elements)


def assemble_ngsolve_hcurl_resistance_shape_tangents(
    fes, vectors, deformation_modes, *, conductivity=1.0, definedon=None
) -> HCurlResistanceLinearization:
    """Assemble reduced ``R`` and analytic HCurl-Piola geometry tangents."""
    import ngsolve as ng

    basis = np.asarray(vectors)
    if basis.ndim == 1:
        basis = basis[:, None]
    if basis.ndim != 2 or basis.shape[0] != fes.ndof:
        raise ValueError(f"vectors must have shape ({fes.ndof}, n_mode)")
    sigma = conductivity
    if np.isscalar(sigma):
        sigma = float(sigma)
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError("conductivity must be positive")
    u, v = fes.TnT()
    measure = ng.dx if definedon is None else ng.dx(definedon=definedon)
    curl_u, curl_v = ng.curl(u), ng.curl(v)
    form = ng.BilinearForm(fes)
    form += (curl_u*curl_v/sigma)*measure
    form.Assemble()
    parent = _csr(form)
    reduced = np.asarray(basis.conj().T @ (parent @ basis))
    derivatives = []
    for mode in tuple(deformation_modes):
        gradient = ng.Grad(mode)
        divergence = ng.div(mode)
        tangent = ng.BilinearForm(fes)
        tangent += (
            ((gradient*curl_u)*curl_v + curl_u*(gradient*curl_v)
             - divergence*curl_u*curl_v)/sigma
        )*measure
        tangent.Assemble()
        matrix = _csr(tangent)
        derivatives.append(np.asarray(basis.conj().T @ (matrix @ basis)))
    jacobian = (
        np.stack(derivatives)
        if derivatives else np.empty((0, reduced.shape[0], reduced.shape[1]))
    )
    return HCurlResistanceLinearization(reduced, jacobian)


def linearize_hcurl_joule_loss_from_ngsolve(
    *, mesh, fes, vectors, interaction, deformation_modes, frequency_hz,
    rhs, conductivity=1.0, rhs_jacobian=None, materials=None,
) -> HCurlJouleAdjointLinearization:
    """Close ``GetTrafo -> dR,dL -> complex adjoint -> Joule gradient``."""
    if not isinstance(interaction, (HCurlTetVolumeInteraction, HCurlCellVolumeInteraction)):
        raise TypeError("interaction must be an HCurl volume interaction")
    inductance = interaction.matrix
    if not isinstance(inductance, HCurlHMatrixOperator):
        raise NotImplementedError("analytic shape adjoints require matrix_free=True")
    modes = tuple(deformation_modes)
    resistance = assemble_ngsolve_hcurl_resistance_shape_tangents(
        fes, vectors, modes, conductivity=conductivity, definedon=materials
    )
    subtet_velocity = SampleNgsolveHCurlCellSubtetVelocities(
        mesh, modes, interaction, materials=materials
    )
    return linearize_hcurl_multifrequency_joule_loss(
        inductance=inductance, resistance=resistance.matrix,
        resistance_jacobian=resistance.jacobian,
        cell_vertex_velocities=subtet_velocity,
        load_cases=(HCurlJouleLoadCase(
            frequency_hz, np.asarray(rhs), rhs_jacobian=rhs_jacobian
        ),),
    ).cases[0]


def linearize_hcurl_multifrequency_joule_loss_from_ngsolve(
    *, mesh, fes, vectors, interaction, deformation_modes, load_cases,
    conductivity=1.0, materials=None,
) -> HCurlMultiFrequencyJouleLinearization:
    """Aggregate analytic Joule adjoints over frequencies and excitations.

    Geometry-dependent ``R``, ``dR``, ``L`` and GetTrafo velocities are built
    once.  Each load case then contributes its weighted complex adjoint without
    materialising any directional ``dL`` matrix.
    """
    if not isinstance(interaction, (HCurlTetVolumeInteraction, HCurlCellVolumeInteraction)):
        raise TypeError("interaction must be an HCurl volume interaction")
    inductance = interaction.matrix
    if not isinstance(inductance, HCurlHMatrixOperator):
        raise NotImplementedError("analytic shape adjoints require matrix_free=True")
    cases = tuple(load_cases)
    if not cases:
        raise ValueError("load_cases must not be empty")
    modes = tuple(deformation_modes)
    resistance = assemble_ngsolve_hcurl_resistance_shape_tangents(
        fes, vectors, modes, conductivity=conductivity, definedon=materials
    )
    subtet_velocity = SampleNgsolveHCurlCellSubtetVelocities(
        mesh, modes, interaction, materials=materials
    )
    return linearize_hcurl_multifrequency_joule_loss(
        inductance=inductance, resistance=resistance.matrix,
        resistance_jacobian=resistance.jacobian,
        cell_vertex_velocities=subtet_velocity, load_cases=cases,
    )


def linearize_hcurl_multifrequency_activation_joule_loss_from_ngsolve(
    *, mesh, fes, vectors, interaction, activation, load_cases,
    conductivity, inductance_power=1.0, materials=None,
) -> HCurlMultiFrequencyJouleLinearization:
    """Analytic eddy-bubble material-topology adjoint for every parent cell."""
    import ngsolve as ng

    if not isinstance(conductivity, HCurlConductivityInterpolation):
        raise TypeError("conductivity must be HCurlConductivityInterpolation")
    inductance = interaction.matrix
    if not isinstance(inductance, HCurlHMatrixOperator):
        raise NotImplementedError("activation adjoints require matrix_free=True")
    rho = np.asarray(activation, dtype=float).reshape(-1)
    sigma, dsigma = conductivity.evaluate(rho)
    inverse_sigma = 1.0/sigma
    inverse_sigma_derivative = -dsigma/sigma**2
    inverse_cf, elements = _ngsolve_cellwise_coefficient(
        mesh, inverse_sigma, materials
    )
    basis = np.asarray(vectors)
    if basis.ndim == 1:
        basis = basis[:, None]
    if basis.shape[0] != fes.ndof:
        raise ValueError(f"vectors must have shape ({fes.ndof}, n_mode)")
    u, v = fes.TnT()
    measure = ng.dx if materials is None else ng.dx(definedon=materials)
    form = ng.BilinearForm(fes)
    form += inverse_cf*ng.curl(u)*ng.curl(v)*measure
    form.Assemble()
    parent = _csr(form)
    resistance = np.asarray(basis.conj().T@(parent@basis))
    l_dense = inductance.activation_to_dense(rho, power=inductance_power)
    modes = []
    for column in range(basis.shape[1]):
        field = ng.GridFunction(fes)
        field.vec.FV().NumPy()[:] = np.asarray(basis[:, column], dtype=float)
        modes.append(ng.curl(field))
    cases = tuple(load_cases)
    if not cases:
        raise ValueError("load_cases must not be empty")
    results = []
    weights = np.empty(len(cases))
    total_objective = 0.0
    total_gradient = np.zeros(len(rho))
    element_numbers = np.asarray([element.nr for element in elements], dtype=int)
    for index, case in enumerate(cases):
        if not isinstance(case, HCurlJouleLoadCase):
            raise TypeError("load_cases entries must be HCurlJouleLoadCase")
        weight = float(case.weight)
        if not np.isfinite(weight) or weight < 0:
            raise ValueError("load-case weights must be finite and non-negative")
        omega = 2*np.pi*float(case.frequency_hz)
        if not np.isfinite(omega) or omega <= 0:
            raise ValueError("frequency_hz must be positive")
        operator = resistance + 1j*omega*l_dense
        load = np.asarray(case.rhs, dtype=complex).reshape(-1)
        if load.shape != (operator.shape[0],):
            raise ValueError(f"rhs must have shape ({operator.shape[0]},)")
        state = np.linalg.solve(operator, load)
        objective = 0.5*float(np.real(np.vdot(state, resistance@state)))
        adjoint = np.linalg.solve(operator.conj().T, resistance@state)
        curl_state = sum((state[j]*modes[j] for j in range(len(modes))), ng.CF((0,0,0)))
        curl_adjoint = sum((adjoint[j]*modes[j] for j in range(len(modes))), ng.CF((0,0,0)))
        direct_all = np.asarray(ng.Integrate(
            ng.InnerProduct(curl_state, curl_state), mesh, element_wise=True
        ))
        adjoint_all = np.asarray(ng.Integrate(
            ng.InnerProduct(curl_adjoint, curl_state), mesh, element_wise=True
        ))
        dL = np.asarray(inductance.activation_contractions(
            rho, adjoint, state, power=inductance_power
        ))
        if case.rhs_jacobian is None:
            db = np.zeros((len(rho), len(load)), dtype=complex)
        else:
            db = np.asarray(case.rhs_jacobian, dtype=complex)
            if db.shape != (len(rho), len(load)):
                raise ValueError(
                    f"rhs_jacobian must have shape ({len(rho)},{len(load)})"
                )
        direct = 0.5*np.real(direct_all[element_numbers])*inverse_sigma_derivative
        implicit_r = -np.real(adjoint_all[element_numbers])*inverse_sigma_derivative
        implicit_b = np.real(db@np.conj(adjoint))
        gradient = direct + implicit_r + implicit_b - np.real(1j*omega*dL)
        results.append(HCurlJouleAdjointLinearization(
            state, adjoint, objective, gradient, resistance,
            np.empty((0,)+resistance.shape), inductance,
        ))
        weights[index] = weight
        total_objective += weight*objective
        total_gradient += weight*gradient
    return HCurlMultiFrequencyJouleLinearization(
        tuple(results), float(total_objective), total_gradient, weights
    )


def linearize_and_solve_hcurl_sheet_joule_lp(
    *, state, mesh, fes, vectors, interaction, deformation_modes, load_cases,
    design_mode_jacobian, area, volume_max, displacement_move,
    thickness_move, activation_move, thickness_bounds,
    conductivity=1.0, materials=None, additional_sheet_gradient=None,
    laplacian=None, curvature_limit=None,
) -> HCurlSheetLPStep:
    """Close multi-frequency HCurl adjoints into the HEX-sheet LP contract.

    ``design_mode_jacobian`` maps packed sheet increments
    ``[normal displacement, thickness, activation]`` to the analytic
    deformation coordinates.  Material/topology sensitivities can be supplied
    explicitly through ``additional_sheet_gradient``; they are never silently
    replaced by zero physics.
    """
    from .sheet_metal_optimization import solve_sheet_metal_lp

    linearization = linearize_hcurl_multifrequency_joule_loss_from_ngsolve(
        mesh=mesh, fes=fes, vectors=vectors, interaction=interaction,
        deformation_modes=deformation_modes, load_cases=load_cases,
        conductivity=conductivity, materials=materials,
    )
    n = np.asarray(state.normal_displacement).size
    mapping = np.asarray(design_mode_jacobian, dtype=float)
    if mapping.shape != (linearization.gradient.size, 3*n):
        raise ValueError(
            f"design_mode_jacobian must have shape ({linearization.gradient.size},{3*n})"
        )
    sheet_gradient = mapping.T@linearization.gradient
    if additional_sheet_gradient is not None:
        extra = np.asarray(additional_sheet_gradient, dtype=float).reshape(-1)
        if extra.shape != (3*n,):
            raise ValueError(f"additional_sheet_gradient must have shape ({3*n},)")
        sheet_gradient = sheet_gradient + extra
    update = solve_sheet_metal_lp(
        state.normal_displacement, state.thickness, state.activation,
        sheet_gradient, area, volume_max=volume_max,
        displacement_move=displacement_move, thickness_move=thickness_move,
        activation_move=activation_move, thickness_bounds=thickness_bounds,
        laplacian=laplacian, curvature_limit=curvature_limit,
    )
    return HCurlSheetLPStep(linearization, update)


def linearize_and_solve_hcurl_activation_sheet_joule_lp(
    *, state, mesh, fes, vectors, interaction, load_cases, conductivity,
    area, volume_max, activation_move, inductance_power=1.0, materials=None,
) -> HCurlSheetLPStep:
    """Close pure material-topology sensitivities into activation-only LP."""
    from .sheet_metal_optimization import solve_sheet_metal_lp

    linearization = linearize_hcurl_multifrequency_activation_joule_loss_from_ngsolve(
        mesh=mesh, fes=fes, vectors=vectors, interaction=interaction,
        activation=state.activation, load_cases=load_cases,
        conductivity=conductivity, inductance_power=inductance_power,
        materials=materials,
    )
    n = np.asarray(state.activation).size
    if linearization.gradient.shape != (n,):
        raise RuntimeError("activation gradient and sheet-cell count differ")
    gradient = np.r_[np.zeros(2*n), linearization.gradient]
    frozen_move = np.finfo(float).eps
    update = solve_sheet_metal_lp(
        state.normal_displacement, state.thickness, state.activation,
        gradient, area, volume_max=volume_max, displacement_move=frozen_move,
        thickness_move=frozen_move, activation_move=activation_move,
        thickness_bounds=(float(np.min(state.thickness)),
                          float(np.max(state.thickness))),
    )
    return HCurlSheetLPStep(linearization, update)


def optimize_hcurl_eddy_bubble_hex_sheet(
    initial_state, *, build_step_inputs, **hex_driver_options,
):
    """Run the deform-first/Cubit HEX driver with HCurl Joule adjoints.

    ``build_step_inputs(state)`` rebuilds the NGSolve/HCurl objects owned by
    the application after each deformation or remesh and returns keyword
    arguments for :func:`linearize_and_solve_hcurl_sheet_joule_lp`.  The
    returned step already exposes the ``.update`` contract consumed by the
    shared HEX-sheet driver.
    """
    from .sheet_metal_optimization import optimize_hex_sheet_topology

    if not callable(build_step_inputs):
        raise TypeError("build_step_inputs must be callable")
    if "linearize_step" in hex_driver_options:
        raise ValueError("linearize_step is supplied by the HCurl driver")

    def linearize_step(state):
        inputs = dict(build_step_inputs(state))
        inputs.setdefault("state", state)
        inputs.setdefault("mesh", state.mesh)
        return linearize_and_solve_hcurl_sheet_joule_lp(**inputs)

    return optimize_hex_sheet_topology(
        initial_state, linearize_step=linearize_step, **hex_driver_options
    )


def optimize_hcurl_eddy_bubble_activation_hex_sheet(
    initial_state, *, build_step_inputs, **hex_driver_options,
):
    """Run pure activation topology with deform-first/Cubit remeshing."""
    from .sheet_metal_optimization import optimize_hex_sheet_topology

    if not callable(build_step_inputs):
        raise TypeError("build_step_inputs must be callable")
    if "linearize_step" in hex_driver_options:
        raise ValueError("linearize_step is supplied by the HCurl driver")

    def linearize_step(state):
        inputs = dict(build_step_inputs(state))
        inputs.setdefault("state", state)
        inputs.setdefault("mesh", state.mesh)
        return linearize_and_solve_hcurl_activation_sheet_joule_lp(**inputs)

    return optimize_hex_sheet_topology(
        initial_state, linearize_step=linearize_step, **hex_driver_options
    )


def solve_hcurl_joule_lp(
    design, objective_gradient, *, move_limit, lower_bounds=None,
    upper_bounds=None, volume_gradient=None, volume_limit=None,
) -> HCurlJouleLPUpdate:
    """Solve one trust-region LP step in the supplied analytic design modes."""
    x = np.asarray(design, dtype=float).reshape(-1)
    gradient = np.asarray(objective_gradient, dtype=float).reshape(-1)
    if gradient.shape != x.shape:
        raise ValueError("objective_gradient must match design")
    move = np.broadcast_to(np.asarray(move_limit, dtype=float), x.shape)
    if np.any(move < 0):
        raise ValueError("move_limit must be non-negative")
    lo = np.full_like(x, -np.inf) if lower_bounds is None else np.broadcast_to(lower_bounds, x.shape)
    hi = np.full_like(x, np.inf) if upper_bounds is None else np.broadcast_to(upper_bounds, x.shape)
    bounds = list(zip(np.maximum(-move, lo-x), np.minimum(move, hi-x)))
    A_ub = b_ub = None
    if volume_gradient is not None or volume_limit is not None:
        if volume_gradient is None or volume_limit is None:
            raise ValueError("volume_gradient and volume_limit must be supplied together")
        vg = np.asarray(volume_gradient, dtype=float).reshape(-1)
        if vg.shape != x.shape:
            raise ValueError("volume_gradient must match design")
        A_ub = vg[None, :]
        b_ub = np.asarray([float(volume_limit)])
    solved = linprog(gradient, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    if not solved.success:
        raise RuntimeError(f"HCurl Joule LP failed: {solved.message}")
    delta = np.asarray(solved.x)
    return HCurlJouleLPUpdate(x+delta, delta, gradient, str(solved.message))


__all__ = [
    "HCurlResistanceLinearization", "HCurlJouleAdjointLinearization",
    "HCurlJouleLPUpdate", "HCurlJouleLoadCase",
    "HCurlMultiFrequencyJouleLinearization", "HCurlSheetLPStep",
    "HCurlConductivityInterpolation",
    "assemble_ngsolve_hcurl_resistance_shape_tangents",
    "linearize_hcurl_multifrequency_joule_loss",
    "linearize_hcurl_multifrequency_activation_joule_loss",
    "linearize_hcurl_joule_loss_from_ngsolve", "solve_hcurl_joule_lp",
    "linearize_hcurl_multifrequency_joule_loss_from_ngsolve",
    "linearize_hcurl_multifrequency_activation_joule_loss_from_ngsolve",
    "linearize_and_solve_hcurl_sheet_joule_lp",
    "linearize_and_solve_hcurl_activation_sheet_joule_lp",
    "optimize_hcurl_eddy_bubble_hex_sheet",
    "optimize_hcurl_eddy_bubble_activation_hex_sheet",
]
