"""Phase F-2 test: HDiv div-free hex basis matches Mathematica reference.

The reference values come from `scripts/export_basis_test_values.wls`, which
runs the canonical `ShapesHDivHexDivFree` from
`hex_bem_hdiv_order3_overnight.wls` at 5 test points for orders 1..3, and
dumps the 3-component vector value of every basis function at every point.

This test asserts that `radia_vim.HDivDivFreeHexBasis(order).evaluate(idx, x, y, z)`
agrees with that Mathematica reference to a relative tolerance of ~1e-12
(double precision floor; absolute tolerance for components that are zero).
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import pytest

from radia_vim import HDivDivFreeHexBasis


REF_PATH = Path(__file__).parent / "basis_test_values.json"


def load_reference():
    if not REF_PATH.exists():
        pytest.skip(f"Mathematica reference not generated: {REF_PATH}\n"
                    f"Run: wolframscript -file scripts/export_basis_test_values.wls")
    return json.loads(REF_PATH.read_text(encoding="utf-8"))


def assert_vec_close(got, ref, *, rtol=1e-12, atol=1e-12, label=""):
    """numpy.allclose-style: diff <= atol + rtol * |ref|."""
    for c in range(3):
        diff = abs(got[c] - ref[c])
        threshold = atol + rtol * abs(ref[c])
        assert diff <= threshold, (
            f"{label}: comp {c} got={got[c]:+.16e} ref={ref[c]:+.16e} "
            f"diff={diff:.3e}  threshold={threshold:.3e}  "
            f"(rtol={rtol}, atol={atol})")


def test_basis_orders_1_2_3():
    data = load_reference()
    test_points = data["test_points"]
    for entry in data["by_order"]:
        order = entry["order"]
        n_dofs_ref = entry["n_dofs"]
        values_ref = entry["values"]   # shape: [n_pt][n_dof][3]

        basis = HDivDivFreeHexBasis(order=order)
        assert basis.order == order
        # Expected n_dofs = 2*p^3 + 3*p^2
        expected_n = 2 * order**3 + 3 * order**2
        assert basis.n_dofs == expected_n == n_dofs_ref

        for p_idx, pt in enumerate(test_points):
            x, y, z = pt
            for idx in range(n_dofs_ref):
                got = basis.evaluate(idx, x, y, z)
                ref = values_ref[p_idx][idx]
                assert_vec_close(got, ref,
                                 label=f"order={order} pt={pt} dof={idx}")


def test_n_dofs_formula():
    """n_dofs = 2*p^3 + 3*p^2 for orders 1..6."""
    expected = {1: 2 + 3, 2: 16 + 12, 3: 54 + 27, 4: 128 + 48,
                5: 250 + 75, 6: 432 + 108}
    for p, n in expected.items():
        b = HDivDivFreeHexBasis(order=p)
        assert b.n_dofs == n, f"order={p}: got {b.n_dofs}, expected {n}"


def test_dof_key_blocks():
    """Block layout: A (p^3) | B (p^3) | C (p^2) | D (p^2) | E (p^2)."""
    p = 3
    b = HDivDivFreeHexBasis(order=p)
    from radia_vim import HexBasisBlock
    counts = {block: 0 for block in (HexBasisBlock.A, HexBasisBlock.B,
                                     HexBasisBlock.C, HexBasisBlock.D,
                                     HexBasisBlock.E)}
    for idx in range(b.n_dofs):
        counts[b.dof_key(idx).block] += 1
    assert counts[HexBasisBlock.A] == p**3
    assert counts[HexBasisBlock.B] == p**3
    assert counts[HexBasisBlock.C] == p**2
    assert counts[HexBasisBlock.D] == p**2
    assert counts[HexBasisBlock.E] == p**2


def test_invalid_order_raises():
    with pytest.raises(Exception):
        HDivDivFreeHexBasis(order=0)


def test_jn_zero_on_boundary_face_x_eq_0():
    """At x=0 (ξ=-1) the x-component of every basis function must vanish.

    This verifies the J·n=0 BC built into the integrated-Legendre basis.
    """
    b = HDivDivFreeHexBasis(order=3)
    # Try a few (y, z) values
    for y in (0.1, 0.5, 0.9):
        for z in (0.1, 0.5, 0.9):
            for idx in range(b.n_dofs):
                v = b.evaluate(idx, 0.0, y, z)
                assert abs(v[0]) < 1e-13, (
                    f"J.n nonzero at x=0: dof={idx} (y,z)=({y},{z}) vx={v[0]:.3e}")


def test_jn_zero_on_boundary_face_x_eq_1():
    b = HDivDivFreeHexBasis(order=3)
    for y in (0.1, 0.5, 0.9):
        for z in (0.1, 0.5, 0.9):
            for idx in range(b.n_dofs):
                v = b.evaluate(idx, 1.0, y, z)
                assert abs(v[0]) < 1e-13, (
                    f"J.n nonzero at x=1: dof={idx} (y,z)=({y},{z}) vx={v[0]:.3e}")


def test_jn_zero_all_six_faces():
    b = HDivDivFreeHexBasis(order=3)
    test_pts = [
        (0.0, 0.5, 0.5),  (1.0, 0.5, 0.5),       # x faces
        (0.5, 0.0, 0.5),  (0.5, 1.0, 0.5),       # y faces
        (0.5, 0.5, 0.0),  (0.5, 0.5, 1.0),       # z faces
    ]
    normals = [(1,0,0), (1,0,0), (0,1,0), (0,1,0), (0,0,1), (0,0,1)]
    for (x, y, z), n in zip(test_pts, normals):
        for idx in range(b.n_dofs):
            v = b.evaluate(idx, x, y, z)
            jn = v[0]*n[0] + v[1]*n[1] + v[2]*n[2]
            assert abs(jn) < 1e-13, (
                f"J.n!=0: dof={idx} pt=({x},{y},{z}) n={n} jn={jn:.3e}")


if __name__ == "__main__":
    test_basis_orders_1_2_3()
    test_n_dofs_formula()
    test_dof_key_blocks()
    test_invalid_order_raises()
    test_jn_zero_on_boundary_face_x_eq_0()
    test_jn_zero_on_boundary_face_x_eq_1()
    test_jn_zero_all_six_faces()
    print("All Phase F-2 basis tests passed.")
