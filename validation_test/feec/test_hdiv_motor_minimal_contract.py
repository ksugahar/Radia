"""Minimal radia-motor contract for planar HDiv-VIM.

This is not a full motor model.  It is the smallest rotating-machine-style
gate: a saliency body, a rotating applied field in the rotor frame, one cached
body operator, and torque that is odd in electrical angle and matches the
closed-form reluctance torque.
"""
from __future__ import annotations

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import OCCGeometry, WorkPlane  # noqa: E402

from radia.vim import Solve  # noqa: E402

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

    body = None
    torques = []
    with ng.TaskManager():
        for th in (theta, -theta):
            Ha, Hb = H0 * np.cos(th), H0 * np.sin(th)
            if body is None:
                res = Solve(mesh, chi + 1.0, ng.CoefficientFunction((Ha, Hb)))
                body = res["body"]
                m = res["m"]
            else:
                m = body.solve_linear(chi, body.project(ng.CoefficientFunction((Ha, Hb))))
            Mx, My = body.M_avg(m)
            torques.append(MU0 * area * (Mx * Hb - My * Ha))

    assert torques[0] * torques[1] < 0.0, f"saliency torque must reverse with angle: {torques}"
    odd_rel = abs(torques[0] + torques[1]) / max(abs(torques[0]), 1e-30)
    assert odd_rel < 2e-3, f"torque is not odd in rotor angle: {torques}"
    ref = closed_form_torque(theta)
    assert abs(torques[0] - ref) / abs(ref) < 2e-2, \
        f"HDiv motor torque {torques[0]:.3g} vs closed form {ref:.3g}"
