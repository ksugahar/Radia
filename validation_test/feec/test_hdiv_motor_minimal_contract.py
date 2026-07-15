"""Minimal radia-motor contract for planar HDiv-VIM.

This is the production reduced reluctance-motor gate: a saliency rotor, one
cached body operator, and torque read independently from Maxwell stress,
magnetization volume coupling, and fixed-current coenergy.
"""
from __future__ import annotations

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import OCCGeometry, WorkPlane  # noqa: E402

from radia.motor_hdiv import HDivReducedMotor  # noqa: E402

MU0 = 4.0e-7 * np.pi


def test_planar_hdiv_motor_saliency_angle_sweep_contract():
    a_el, b_el = 0.2, 0.1
    chi = 1000.0
    H0 = 8.0e4
    theta = np.radians(25.0)
    mesh = ng.Mesh(OCCGeometry(WorkPlane().Ellipse(a_el, b_el).Face(), dim=2).GenerateMesh(maxh=b_el / 2.5))
    area = np.pi * a_el * b_el
    Na, Nb = b_el / (a_el + b_el), a_el / (a_el + b_el)

    def closed_form_torque(th):
        Ha, Hb = H0 * np.cos(th), H0 * np.sin(th)
        Ma = chi * Ha / (1.0 + chi * Na)
        Mb = chi * Hb / (1.0 + chi * Nb)
        return MU0 * area * (Ma * Hb - Mb * Ha)

    torque_sets = []
    with ng.TaskManager():
        motor = HDivReducedMotor(mesh, chi + 1.0)
        # A fixed global field and rotor angles -/+theta correspond to local
        # field angles +/-theta.
        for rotor_angle in (-theta, theta):
            state = motor.solve_angle(rotor_angle, (H0, 0.0))
            torque_sets.append((
                motor.maxwell_torque(state, 0.28, circle_points=1440),
                state.torque_volume_Nm,
                motor.virtual_work_torque(
                    rotor_angle, (H0, 0.0), delta_angle=np.radians(0.1)),
            ))

    assert motor.gram_build_count == 1
    positive, negative = torque_sets
    assert positive[0] * negative[0] < 0.0, \
        f"saliency torque must reverse with angle: {torque_sets}"
    odd_rel = abs(positive[0] + negative[0]) / max(abs(positive[0]), 1e-30)
    assert odd_rel < 2e-3, f"torque is not odd in rotor angle: {torque_sets}"
    for routes in torque_sets:
        spread = max(routes)-min(routes)
        assert spread / max(abs(value) for value in routes) < 5e-5, routes
    ref = closed_form_torque(theta)
    assert abs(positive[0] - ref) / abs(ref) < 2e-2, \
        f"HDiv motor torque {positive[0]:.3g} vs closed form {ref:.3g}"
