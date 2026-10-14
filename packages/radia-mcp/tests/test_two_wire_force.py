r"""Ampere's force law: force per unit length between two parallel wires -- regression test.

Two parallel wires a distance d apart carrying currents I1, I2 attract/repel with
F = mu0 I1 I2/(2 pi d) (two_wire_force_per_length); like currents attract. The analytic primitive of
the image-wire force (image_force_wire_iron = this with the image at 2h). The FE leg gets the force
as the Lorentz body force int J x B over one wire (only the OTHER wire's field nets a force; the
symmetric self-field integrates to zero). Tool-independent (Ampere's force law).
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (two_wire_force_per_length, image_force_wire_iron,
                                           MU0)
from radia_mcp.radia_ngsolve.force import (
    parallel_wire_lorentz_force_summary,
    parallel_wire_virtual_work_force_summary,
    planar_lorentz_force_summary,
)

a, d, I = 0.002, 0.02, 1.0


def test_two_wire_force_closed_form():
    # Ampere's force law + the sign/relationship conventions
    F = two_wire_force_per_length(I, I, d)
    assert F == pytest.approx(MU0 * I * I / (2 * math.pi * d))      # mu0 I1 I2/(2 pi d)
    assert F == pytest.approx(1.0e-5)                               # I=1, d=20mm -> 1e-5 N/m
    assert F > 0                                                    # like currents (I1 I2>0) attract
    assert two_wire_force_per_length(I, -I, d) < 0                  # opposite currents repel (sign flip)
    assert two_wire_force_per_length(2 * I, I, d) == pytest.approx(2 * F)   # linear in each current
    assert two_wire_force_per_length(I, I, 2 * d) == pytest.approx(F / 2)   # 1/d falloff
    # image-wire force is the special case I1=I2=I at separation 2h
    h = d / 2
    assert image_force_wire_iron(I, h) == pytest.approx(two_wire_force_per_length(I, I, 2 * h))


def test_planar_lorentz_block_force_matches_two_wire_field():
    area = math.pi * a * a
    jz_right = I / area
    by_from_left = MU0 * I / (2.0 * math.pi * d)
    row = planar_lorentz_force_summary(jz_right, (0.0, by_from_left), area_m2=area)
    pair = parallel_wire_lorentz_force_summary(I, I, d)

    assert row["current_A"] == pytest.approx(I)
    assert row["force_per_depth_N_per_m"] == pytest.approx([-two_wire_force_per_length(I, I, d), 0.0])
    assert row["force_magnitude_per_depth_N_per_m"] == pytest.approx(two_wire_force_per_length(I, I, d))
    assert row["force_per_depth_N_per_m"] == pytest.approx(pair["force_on_wire2_N_per_m"])


def test_parallel_wire_lorentz_force_summary_tracks_direction_and_sign():
    like = parallel_wire_lorentz_force_summary(I, I, d)
    assert like["interaction"] == "attraction"
    assert like["field_from_wire1_at_wire2_T"] == pytest.approx([0.0, MU0 * I / (2.0 * math.pi * d)])
    assert like["signed_ampere_force_per_length_N_per_m"] == pytest.approx(two_wire_force_per_length(I, I, d))
    assert like["force_on_wire2_N_per_m"] == pytest.approx([-1.0e-5, 0.0])
    assert like["force_on_wire1_N_per_m"] == pytest.approx([1.0e-5, 0.0])

    vertical = parallel_wire_lorentz_force_summary(I, I, (0.0, d))
    assert vertical["field_from_wire1_at_wire2_T"] == pytest.approx([-MU0 * I / (2.0 * math.pi * d), 0.0])
    assert vertical["force_on_wire2_N_per_m"] == pytest.approx([0.0, -1.0e-5])

    opposite = parallel_wire_lorentz_force_summary(I, -I, d)
    assert opposite["interaction"] == "repulsion"
    assert opposite["force_on_wire2_N_per_m"] == pytest.approx([1.0e-5, 0.0])

    with pytest.raises(ValueError):
        parallel_wire_lorentz_force_summary(I, I, 0.0)
    with pytest.raises(ValueError):
        parallel_wire_lorentz_force_summary(I, I, (0.0, 0.0))


def test_parallel_wire_virtual_work_force_matches_lorentz_direction():
    like = parallel_wire_virtual_work_force_summary(I, I, d, displacement_step_m=d * 1.0e-4)
    expected = two_wire_force_per_length(I, I, d)

    assert like["interaction"] == "attraction"
    assert like["analytic_radial_force_per_length_N_per_m"] == pytest.approx(-expected)
    assert like["virtual_work_radial_force_per_length_N_per_m"] == pytest.approx(-expected, rel=1.0e-8)
    assert like["virtual_work_force_on_wire2_N_per_m"] == pytest.approx(like["lorentz_force_on_wire2_N_per_m"], rel=1.0e-8)
    assert like["force_rel_error"] < 1.0e-8
    assert like["coenergy_minus_J_per_m"] > like["coenergy_plus_J_per_m"]

    vertical_repulsion = parallel_wire_virtual_work_force_summary(I, -I, (0.0, d))
    assert vertical_repulsion["interaction"] == "repulsion"
    assert vertical_repulsion["analytic_radial_force_per_length_N_per_m"] > 0.0
    assert vertical_repulsion["virtual_work_force_on_wire2_N_per_m"][1] > 0.0

    with pytest.raises(ValueError):
        parallel_wire_virtual_work_force_summary(I, I, d, displacement_step_m=0.0)
    with pytest.raises(ValueError):
        parallel_wire_virtual_work_force_summary(I, I, d, displacement_step_m=d)
    with pytest.raises(ValueError):
        parallel_wire_virtual_work_force_summary(I, I, d, reference_separation_m=0.0)


def test_planar_lorentz_force_rejects_bad_inputs():
    with pytest.raises(ValueError):
        planar_lorentz_force_summary(1.0, (0.0, 1.0, 0.0))
    with pytest.raises(ValueError):
        planar_lorentz_force_summary(1.0, (0.0, 1.0), area_m2=-1.0)
