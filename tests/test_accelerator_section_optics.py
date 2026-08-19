"""Section-optics inverse design against closed forms and structural identities.

Everything here runs without a field solve, so the design model can be
checked before it is pointed at a magnet where a wrong sign or a missing
factor would be invisible among the physics.

The straight slab of constant gradient is the one section whose map is
known exactly,

    K = g / (B rho),  R = [[cos wL, sin wL / w], [-w sin wL, cos wL]]

with ``w = sqrt(K)``, so ``R[1,0] -> -K L`` as ``L -> 0`` and therefore
``d R[1,0] / d b_1(j) -> -ds_j / (B rho)``.  Those statements exercise the
specification, the analytic Jacobian, the weighted inversion and the
linearization gate in turn.

The composition identity ``M = M_after . M_S . M_before`` is checked
separately because the whole design model rests on it being exact rather
than approximate.

Map builds through the public multipole route carry the design adjoint
unconditionally and cost seconds each, so the slab is deliberately short
and its Jacobian is built once for the module.  The full analytic-versus
-difference comparison over every column lives in the validation lane;
here it is a spot check on a few columns.
"""
import numpy as np
import pytest

from radia.accelerator_section_optics import (
    COMPONENT_COUNT, NORMAL_BLOCKS, aperture_field_metric,
    constrained_minimum_norm_step, focusing_coupling,
    minimum_norm_field_difference, multipole_section_spec,
    multipole_section_spec_jacobian,
    multipole_section_spec_jacobian_by_differences,
    multipole_section_transfer, scatter_multipole_difference,
    section_composition_defect, snap_breaks_to_section,
    verify_multipole_linearization,
)

RIGIDITY = 0.25
GRADIENT = 15.0
SEGMENTS = 4
LENGTH = 0.040
HALF_APERTURE = 0.012


@pytest.fixture(scope="module")
def slab():
    """A straight slab of constant normal quadrupole, and its exact map."""
    lengths = np.full(SEGMENTS, LENGTH / SEGMENTS)
    curvature = np.zeros(SEGMENTS)
    response = np.zeros((COMPONENT_COUNT, SEGMENTS))
    response[1, :] = GRADIENT
    strength = GRADIENT / RIGIDITY
    root = np.sqrt(strength)
    flat = response.reshape(-1)
    jacobian, index = multipole_section_spec_jacobian(
        flat, lengths, curvature, RIGIDITY, 0, SEGMENTS)
    return {
        "response": flat, "lengths": lengths, "curvature": curvature,
        "spec": multipole_section_spec(flat, lengths, curvature, RIGIDITY,
                                       0, SEGMENTS),
        "jacobian": jacobian, "index": index,
        "weight": aperture_field_metric(index, lengths, HALF_APERTURE),
        "exact_horizontal": -root * np.sin(root * LENGTH),
        "exact_vertical": root * np.sinh(root * LENGTH),
    }


def test_slab_specification_matches_the_closed_form(slab):
    spec = slab["spec"]
    assert spec[0] == pytest.approx(slab["exact_horizontal"], rel=2.0e-3)
    assert spec[1] == pytest.approx(slab["exact_vertical"], rel=2.0e-3)
    # a pure quadrupole carries no bend
    assert spec[2] == pytest.approx(0.0, abs=1.0e-12)


def test_analytic_jacobian_agrees_with_central_differences():
    """Two segments only: each map build carries the design adjoint and costs
    seconds, and the agreement is a per-column statement that does not need a
    long slab.  The full comparison over a production chain is a validation run.
    """
    short = 2
    lengths = np.full(short, LENGTH / short)
    curvature = np.zeros(short)
    response = np.zeros((COMPONENT_COUNT, short))
    response[1, :] = GRADIENT
    flat = response.reshape(-1)
    analytic, index = multipole_section_spec_jacobian(
        flat, lengths, curvature, RIGIDITY, 0, short, blocks_used=(0, 1))
    differenced, other = multipole_section_spec_jacobian_by_differences(
        flat, lengths, curvature, RIGIDITY, 0, short, blocks_used=(0, 1))
    assert np.array_equal(index, other)
    assert analytic.shape == differenced.shape == (3, 2 * short)
    for column in range(analytic.shape[1]):
        scale = max(float(np.max(np.abs(differenced[:, column]))), 1.0e-30)
        assert np.max(np.abs(analytic[:, column]
                             - differenced[:, column])) / scale < 1.0e-4


def test_jacobian_rows_match_the_thin_lens_limit(slab):
    jacobian, index = slab["jacobian"], slab["index"]
    assert jacobian.shape == (3, len(NORMAL_BLOCKS) * SEGMENTS)
    quadrupole = [column for column, (block, _s) in enumerate(index)
                  if int(block) == 1]
    dipole = [column for column, (block, _s) in enumerate(index)
              if int(block) == 0]
    slice_length = slab["lengths"][0]
    # Thin lens is the L -> 0 limit; over a finite slab each slice sits at
    # a different betatron phase, so the MEAN column is what converges.
    assert float(np.mean(jacobian[0, quadrupole])) == pytest.approx(
        -slice_length / RIGIDITY, rel=0.05)
    assert float(np.mean(jacobian[1, quadrupole])) == pytest.approx(
        slice_length / RIGIDITY, rel=0.05)
    # the bend row is exact, not differenced
    assert jacobian[2, dipole[0]] == pytest.approx(slice_length, rel=1.0e-12)
    assert np.max(np.abs(jacobian[2, quadrupole])) == pytest.approx(
        0.0, abs=1.0e-12)


def test_analytic_jacobian_zeroes_the_inert_blocks(slab):
    """Only the dipole and quadrupole move a linear-map specification."""
    jacobian, index = slab["jacobian"], slab["index"]
    for block in (3, 5, 7):
        columns = [column for column, (b, _s) in enumerate(index)
                   if int(b) == block]
        assert np.max(np.abs(jacobian[:, columns])) == 0.0


def test_straight_slab_locks_the_two_focusing_planes(slab):
    cosine, multiplier = focusing_coupling(slab["jacobian"], slab["weight"])
    # One gradient sets both planes with opposite sign; only the spread of
    # betatron phase along the slab breaks the tie, so the price of
    # holding one plane while changing the other is large but finite.
    assert abs(cosine) == pytest.approx(1.0, abs=5.0e-3)
    assert 10.0 < multiplier < 2000.0


@pytest.mark.parametrize("fraction", [0.05, 0.20])
def test_recovered_difference_delivers_the_request(slab, fraction):
    # The vertical is left free: this section cannot decouple the planes,
    # and demanding it produces a huge step for no reason.
    requested = np.array([fraction * slab["spec"][0], 0.0, 0.0])
    difference, kept = minimum_norm_field_difference(
        slab["jacobian"][[0, 2]], requested[[0, 2]], slab["weight"])
    full = scatter_multipole_difference(difference, slab["index"],
                                        (COMPONENT_COUNT, SEGMENTS))
    _before, _after, delivered, _ratio = verify_multipole_linearization(
        slab["response"], slab["lengths"], slab["curvature"], RIGIDITY,
        0, SEGMENTS, full, requested)
    assert delivered[0] / requested[0] == pytest.approx(1.0, abs=0.02)
    assert delivered[2] == pytest.approx(0.0, abs=1.0e-12)
    # Only the dipole and quadrupole are design variables here.
    assert int(kept.sum()) == 2 * SEGMENTS
    for block in (3, 5, 7):
        assert np.max(np.abs(full[block])) == 0.0


def test_filter_holds_uninformative_columns_at_zero():
    """A column that is only round-off must not become a design variable.

    The analytic multipole Jacobian returns exact zeros for the blocks
    that cannot move a linear-map specification.  A DIFFERENCED Jacobian
    -- which is what the chain representation must use -- returns
    round-off instead, and a minimum-norm solution will happily spend it,
    the more so because the aperture metric makes high-order coefficients
    cheap.  The contract locked here is that the filter makes the result
    independent of that noise: the uninformative columns come back exactly
    zero, and the informative ones are unchanged by their presence.
    """
    columns = 6
    jacobian = np.zeros((2, columns))
    jacobian[0, :2] = [-2.0e-2, -2.0e-2]          # quadrupole authority
    jacobian[1, 2:4] = [5.0e-3, 5.0e-3]           # bend row
    jacobian[0, 4:] = [3.0e-16, -1.0e-16]         # difference round-off
    weight = np.array([1.4e-4, 1.4e-4, 5.0e-3, 5.0e-3, 4.3e-18, 4.3e-18])
    requested = np.array([-1.0e-2, 0.0])
    filtered, kept = minimum_norm_field_difference(jacobian, requested,
                                                   weight)
    unfiltered, kept_all = minimum_norm_field_difference(
        jacobian, requested, weight, controllability_tolerance=0.0)
    assert int(kept.sum()) == 4 and int(kept_all.sum()) == columns
    # the noise columns are spent when they are kept, and silent when not
    assert np.max(np.abs(unfiltered[4:])) > 0.0
    assert np.max(np.abs(filtered[4:])) == 0.0
    # and dropping them does not disturb the columns that carry the design
    assert np.allclose(filtered[:4], unfiltered[:4], rtol=1.0e-6)


def test_zero_jacobian_row_is_refused(slab):
    broken = slab["jacobian"].copy()
    broken[0] = 0.0
    with pytest.raises(ValueError, match="identically zero Jacobian"):
        minimum_norm_field_difference(broken, np.zeros(3), slab["weight"])


def test_section_composition_identity_is_exact(slab):
    """``M = M_after . M_S . M_before`` is the design model's foundation."""
    split = 2
    lengths, curvature = slab["lengths"], slab["curvature"]
    before = multipole_section_transfer(slab["response"], lengths, curvature,
                                        RIGIDITY, 0, split)
    section = multipole_section_transfer(slab["response"], lengths, curvature,
                                         RIGIDITY, split, SEGMENTS)
    whole = multipole_section_transfer(slab["response"], lengths, curvature,
                                       RIGIDITY, 0, SEGMENTS)
    assert section_composition_defect(before, section, whole) < 1.0e-12
    # and each block is symplectic in its own right
    for matrix in (before, section, whole):
        assert np.linalg.det(matrix[:2, :2]) == pytest.approx(1.0, abs=1.0e-9)
        assert np.linalg.det(matrix[2:4, 2:4]) == pytest.approx(1.0,
                                                                abs=1.0e-9)


def test_snap_breaks_pins_a_geometric_section():
    """A section is a fixed arc-length interval, not a field threshold."""
    breaks = np.linspace(0.0, 0.190, 11)
    target = 0.1548
    snapped, index, moved = snap_breaks_to_section(breaks, target)
    assert snapped[index] == target
    assert moved == pytest.approx(abs(breaks[index] - target))
    assert np.all(np.diff(snapped) > 0.0)
    assert np.array_equal(np.delete(snapped, index), np.delete(breaks, index))
    with pytest.raises(ValueError, match="chain end"):
        snap_breaks_to_section(breaks, 0.0)


def test_dense_and_diagonal_metrics_agree(slab):
    """The dense solver must reduce exactly to the diagonal one."""
    requested = np.array([-1.0e-2, 0.0])
    diagonal, kept_a = minimum_norm_field_difference(
        slab["jacobian"][[0, 2]], requested, slab["weight"])
    dense, kept_b = constrained_minimum_norm_step(
        slab["jacobian"][[0, 2]], requested, np.diag(slab["weight"]))
    assert np.array_equal(kept_a, kept_b)
    assert np.allclose(diagonal, dense, rtol=0.0, atol=1.0e-18)
