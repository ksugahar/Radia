import numpy as np
import pytest

pytest.importorskip("radia._radia_pybind")
from radia.beam import (
    Boris2,
    CartesianState,
    ClassicalRK4,
    LorentzEquation,
    ParticleSpecies,
    ReferenceParticle,
    Tracker,
    TrackPlan,
    UniformField,
    ZeroField,
)


def _state(momentum):
    return CartesianState([0.0, 0.0, 0.0], momentum)


def test_native_reference_particle_uses_si_and_signed_rigidity():
    reference = ReferenceParticle.from_kinetic_energy_ev(
        ParticleSpecies.proton(), 220e6
    )
    assert reference.kinetic_energy_j == pytest.approx(
        220e6 * 1.602176634e-19
    )
    assert reference.magnetic_rigidity_t_m == pytest.approx(
        reference.momentum_kg_m_s / reference.species.charge_c
    )
    electron = ReferenceParticle.from_kinetic_energy_ev(
        ParticleSpecies.electron(), 1e6
    )
    assert electron.magnetic_rigidity_t_m < 0.0


def test_native_field_rhs_and_one_step_are_independently_inspectable():
    species = ParticleSpecies.proton()
    field = UniformField([0.0, 0.0, 0.7])
    equation = LorentzEquation(species, field, independent="time")
    state = _state([2e-19, 0.0, 0.0])

    sample = field.evaluate([1.0, 2.0, 3.0])
    rhs = equation.rhs(0.0, state)
    step = ClassicalRK4().step(equation, 0.0, state, 1e-12)

    np.testing.assert_array_equal(sample.magnetic_t, [0.0, 0.0, 0.7])
    assert rhs.dkinetic_momentum_kg_m_s[1] < 0.0
    assert step.independent_after == pytest.approx(1e-12)
    assert step.state_after.time_s == pytest.approx(1e-12)


def test_native_path_length_rhs_has_unit_tangent_and_unit_path_rate():
    equation = LorentzEquation(
        ParticleSpecies.proton(), ZeroField(), independent="path_length"
    )
    rhs = equation.rhs(0.0, _state([0.0, 3e-19, 4e-19]))
    np.testing.assert_allclose(rhs.dposition_m, [0.0, 0.6, 0.8], atol=1e-15)
    assert rhs.dpath_length_m == pytest.approx(1.0)


def test_native_azimuth_rhs_uses_explicit_cylindrical_rate():
    equation = LorentzEquation(
        ParticleSpecies.proton(), ZeroField(), independent="azimuth"
    )
    state = CartesianState([2.0, 0.0, 0.0], [0.0, 3e-19, 0.0])
    rhs = equation.rhs(0.0, state)
    invariant = equation.invariants(state)
    np.testing.assert_allclose(rhs.dposition_m, [0.0, 2.0, 0.0], atol=2e-15)
    assert rhs.dtime_s == pytest.approx(2.0 / invariant.speed_m_s)
    assert rhs.dpath_length_m == pytest.approx(2.0)


def test_native_boris_closes_uniform_field_orbit_and_preserves_momentum():
    species = ParticleSpecies.proton()
    magnetic_t = 0.7
    momentum = 3e-19
    equation = LorentzEquation(
        species, UniformField([0.0, 0.0, magnetic_t]), independent="time"
    )
    state = _state([momentum, 0.0, 0.0])
    invariant = equation.invariants(state)
    omega = species.charge_c * magnetic_t / (
        invariant.relativistic_gamma * species.rest_mass_kg
    )
    period = 2 * np.pi / omega
    plan = TrackPlan()
    plan.start = 0.0
    plan.stop = period
    plan.maximum_step = period / 8000
    plan.maximum_steps = 8100

    trajectory = Tracker(equation, Boris2()).track(state, plan)

    final = trajectory.samples[-1]
    assert np.linalg.norm(final.kinetic_momentum_kg_m_s) == pytest.approx(
        momentum, rel=2e-13
    )
    np.testing.assert_allclose(final.position_m[:2], 0.0, atol=3e-6)
    assert trajectory.summary.accepted_steps == 8000
    assert trajectory.summary.momentum_conservation_applicable
    assert trajectory.summary.maximum_relative_momentum_error < 2e-13


def test_native_electric_tracking_does_not_claim_momentum_conservation():
    equation = LorentzEquation(
        ParticleSpecies.proton(),
        UniformField([0.0, 0.0, 0.0], electric_v_m=[1e4, 0.0, 0.0]),
    )
    plan = TrackPlan()
    plan.start = 0.0
    plan.stop = 1e-9
    plan.maximum_step = 1e-10
    trajectory = Tracker(equation, ClassicalRK4()).track(
        _state([1e-19, 0.0, 0.0]), plan
    )
    assert not trajectory.summary.momentum_conservation_applicable
    assert np.isnan(trajectory.summary.maximum_relative_momentum_error)


def test_native_tracking_rejects_zero_step_and_step_limit():
    equation = LorentzEquation(ParticleSpecies.proton(), ZeroField())
    state = _state([1e-19, 0.0, 0.0])
    with pytest.raises(ValueError, match="nonzero"):
        ClassicalRK4().step(equation, 0.0, state, 0.0)

    plan = TrackPlan()
    plan.start = 0.0
    plan.stop = 1.0
    plan.maximum_step = 0.1
    plan.maximum_steps = 2
    with pytest.raises(RuntimeError, match="maximum_steps"):
        Tracker(equation, ClassicalRK4()).track(state, plan)


def test_native_tracking_uses_interval_scale_for_terminal_tolerance():
    equation = LorentzEquation(ParticleSpecies.proton(), ZeroField())
    state = _state([1e-19, 0.0, 0.0])
    start = 1e6
    stop = start + 1e-8
    span = stop - start
    plan = TrackPlan()
    plan.start = start
    plan.stop = stop
    plan.maximum_step = span / 10
    plan.maximum_steps = 10

    trajectory = Tracker(equation, ClassicalRK4()).track(state, plan)

    assert trajectory.summary.accepted_steps == 10
    assert trajectory.summary.independent_stop == stop
