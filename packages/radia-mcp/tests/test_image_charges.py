r"""Method-of-images electrostatic forces -- analytic regression test.

Grounded conducting sphere (single image q'=-q a/d at b=a^2/d) and an infinite line
charge above a grounded plane (image -lambda at depth h). Closed-form helpers only; no
field solver. (Jackson, *Classical Electrodynamics*, Secs. 2.1-2.3; Griffiths, Sec. 3.2.)
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.image_charges import (
    EPS0,
    grounded_sphere_image_force,
    line_charge_grounded_plane_image,
)


def test_grounded_sphere_image_force():
    q, a, d = 1.0e-9, 0.05, 0.20
    out = grounded_sphere_image_force(q, a, d)

    # exact compact closed form: F = -(1/4pi eps0) q^2 a d / (d^2-a^2)^2
    eps = EPS0
    expected = -(1.0 / (4.0 * math.pi * eps)) * q * q * a * d / (d * d - a * a) ** 2
    assert math.isclose(out["force"], expected, rel_tol=1e-12)
    assert out["force"] < 0.0                                   # attractive

    # image charge and position match q'=-q a/d, b=a^2/d
    assert math.isclose(out["image_charge"], -q * a / d, rel_tol=1e-12)
    assert math.isclose(out["image_position"], a * a / d, rel_tol=1e-12)

    # GATE 1 (independent): the image-Coulomb expression F=(1/4pi eps) q q'/(d-b)^2
    # must equal the compact closed form to machine precision.
    q_img = -q * a / d
    b = a * a / d
    f_coulomb = (1.0 / (4.0 * math.pi * eps)) * q * q_img / (d - b) ** 2
    assert math.isclose(out["force"], f_coulomb, rel_tol=1e-12)

    # GATE 2 (independent): far field d>>a, F = -(1/4pi eps) q^2 a d/(d^2-a^2)^2
    # -> -(1/4pi eps) q^2 a / d^3 (induced-dipole 1/d^3 law). So |F| d^3 / a
    # approaches the constant (1/4pi eps) q^2; the ratio tightens as d grows.
    def far(dd):
        return abs(grounded_sphere_image_force(q, a, dd)["force"]) * dd ** 3 / a
    const = (1.0 / (4.0 * math.pi * eps)) * q * q
    assert math.isclose(far(2.0), const, rel_tol=4e-3)
    assert math.isclose(far(4.0), const, rel_tol=1e-3)          # tighter as d grows

    # GATE 3 (independent): dielectric scaling F(eps_r=k) = F(1)/k
    k = 2.5
    out_k = grounded_sphere_image_force(q, a, d, eps_r=k)
    assert math.isclose(out_k["force"], out["force"] / k, rel_tol=1e-12)

    # ValueError path: charge inside the sphere (d <= a)
    with pytest.raises(ValueError):
        grounded_sphere_image_force(q, a, a)


def test_line_charge_grounded_plane_image():
    lam, h = 3.0e-9, 0.01
    out = line_charge_grounded_plane_image(lam, h)

    # exact: F' = -lambda^2/(4 pi eps0 h)
    expected = -lam * lam / (4.0 * math.pi * EPS0 * h)
    assert math.isclose(out["force_per_length"], expected, rel_tol=1e-12)
    assert out["force_per_length"] < 0.0                        # attractive

    # GATE 1 (independent): 1/h scaling -- doubling the height halves the force.
    out_2h = line_charge_grounded_plane_image(lam, 2.0 * h)
    assert math.isclose(out_2h["force_per_length"], out["force_per_length"] / 2.0,
                        rel_tol=1e-12)

    # GATE 2 (independent): quadratic in lambda -- tripling lambda gives 9x the force.
    out_3lam = line_charge_grounded_plane_image(3.0 * lam, h)
    assert math.isclose(out_3lam["force_per_length"], 9.0 * out["force_per_length"],
                        rel_tol=1e-12)

    # GATE 3 (independent): dielectric scaling F'(eps_r=k) = F'(1)/k
    k = 4.0
    out_k = line_charge_grounded_plane_image(lam, h, eps_r=k)
    assert math.isclose(out_k["force_per_length"], out["force_per_length"] / k,
                        rel_tol=1e-12)

    # ValueError path: non-positive height
    with pytest.raises(ValueError):
        line_charge_grounded_plane_image(lam, 0.0)


if __name__ == "__main__":
    test_grounded_sphere_image_force()
    test_line_charge_grounded_plane_image()
    print("[OK] method-of-images forces: sphere & line-over-plane closed forms verified.")
