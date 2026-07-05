"""Fast unit tests for sym_<bc>_<axis> boundary-condition label plumbing.

Locks the formulation-agnostic symmetry naming introduced 2026-04-25:

    sym_bn=0_<axis>  -- B . n = 0 on axis=0 plane (field parallel)
                        Omega -> natural (no-op)
                        A     -> Dirichlet (A x n = 0)
                        Radia -> image sign '+'
    sym_ht=0_<axis>  -- H x n = 0 on axis=0 plane (field perpendicular)
                        Omega -> Dirichlet (Omega=const)
                        A     -> natural
                        Radia -> image sign '-'

These tests are fast (<1s total), Cubit-free, and cover the three
Python touch-points:

  1. `add_kelvin.sym_sideset_name` / `parse_sym_label`
  2. `add_kelvin.add_kelvin_cubit(reduction=...)` validation branches
     (1/8 NotImplementedError, invalid bc, invalid axis, offset_dir
     conflict).
  3. `calc_accel_hdiv.ima_from_mesh_labels` auto-assembly.

The Cubit-headless end-to-end run (reduction mode actually producing
`sym_*_*` sidesets in a .vol) lives in a separate slow test.
"""
from __future__ import annotations

import os
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PANELS = os.path.join(ROOT, "src", "radia", "panels")
if PANELS not in sys.path:
    sys.path.insert(0, PANELS)


# ==========================================================
# add_kelvin.sym_sideset_name / parse_sym_label
# ==========================================================

def test_sym_sideset_name_canonical():
    from add_kelvin import sym_sideset_name
    assert sym_sideset_name("x", "bn=0") == "sym_bn=0_x"
    assert sym_sideset_name("y", "ht=0") == "sym_ht=0_y"
    assert sym_sideset_name("Z", "bn=0") == "sym_bn=0_z"   # lowercased


def test_sym_sideset_name_rejects_bad_axis():
    from add_kelvin import sym_sideset_name
    with pytest.raises(ValueError, match="axis"):
        sym_sideset_name("w", "bn=0")


def test_sym_sideset_name_rejects_bad_bc():
    from add_kelvin import sym_sideset_name
    with pytest.raises(ValueError, match="bc"):
        sym_sideset_name("x", "sym_tangential")
    with pytest.raises(ValueError, match="bc"):
        sym_sideset_name("x", "dirichlet")


def test_parse_sym_label_canonical_roundtrip():
    from add_kelvin import sym_sideset_name, parse_sym_label
    for axis in ("x", "y", "z"):
        for bc in ("bn=0", "ht=0"):
            label = sym_sideset_name(axis, bc)
            assert parse_sym_label(label) == (axis, bc), label


def test_parse_sym_label_legacy():
    """Legacy axis-less labels survive for pre-2026-04-25 samples."""
    from add_kelvin import parse_sym_label
    assert parse_sym_label("sym_tangential") == ("", "bn=0")
    assert parse_sym_label("sym_normal") == ("", "ht=0")


def test_parse_sym_label_non_sym():
    from add_kelvin import parse_sym_label
    for label in ("", "air", "kelvin_int", "GND", "iron",
                  "sym_bogus_x", "sym_bn=0_w", None, 42):
        assert parse_sym_label(label) is None, label


# ==========================================================
# add_kelvin.add_kelvin_cubit(reduction=...) validation
# ==========================================================

def test_reduction_and_symmetry_are_mutually_exclusive():
    """Passing both symmetry= and reduction= must raise early."""
    from add_kelvin import add_kelvin_cubit
    with pytest.raises(ValueError, match="mutually exclusive"):
        add_kelvin_cubit(R=0.1, symmetry=["z"],
                         reduction={"x": "bn=0"})


def test_reduction_eighth_all_bn_rejected():
    """1/8 reduction with all 3 planes set to bn=0 means B parallel to
    three mutually perpendicular planes, which forces B = 0 everywhere.
    Must fail loud with a physics-aware message."""
    from add_kelvin import add_kelvin_cubit
    with pytest.raises(ValueError, match="physically impossible"):
        add_kelvin_cubit(
            R=0.1,
            reduction={"x": "bn=0", "y": "bn=0", "z": "bn=0"})


def test_reduction_eighth_with_ht_axis_validates():
    """1/8 with at least one ht=0 axis must pass argument validation
    (no more NotImplementedError -- supported as of 2026-04-25).

    Cubit is importable on the LAB machine so the helper proceeds
    past validation into the actual Cubit work, which then fails
    because no air block exists in the empty Cubit session.  We
    accept any RuntimeError / NotImplementedError NOT containing
    the obsolete 'reduction supports 1 or 2 axes' / '1/8 reduction
    (x+y+z) is not yet supported' messages.
    """
    from add_kelvin import add_kelvin_cubit
    # Catch broadly -- the goal is just to assert that the
    # *validation* block doesn't raise.  If validation rejected this
    # case, we'd get a ValueError saying "reduction supports 1 or 2"
    # or NotImplementedError saying "1/8 ... not yet supported"; both
    # would be caught here, so we additionally check the message.
    try:
        add_kelvin_cubit(
            R=0.1,
            reduction={"x": "ht=0", "y": "ht=0", "z": "bn=0"})
    except (RuntimeError, ImportError, ModuleNotFoundError):
        # Expected: cubit not in a valid state (no air block).
        return
    except (ValueError, NotImplementedError) as e:
        msg = str(e).lower()
        assert "1 or 2 axes" not in msg, \
            f"validation still rejecting 1/8: {e}"
        assert "not yet supported" not in msg, \
            f"validation still has obsolete 'not yet supported': {e}"
        # Any OTHER ValueError is also unexpected for this valid input.
        raise


def test_reduction_empty_rejected():
    from add_kelvin import add_kelvin_cubit
    with pytest.raises(ValueError, match="reduction dict is empty"):
        add_kelvin_cubit(R=0.1, reduction={})


def test_reduction_invalid_bc_rejected():
    from add_kelvin import add_kelvin_cubit
    with pytest.raises(ValueError, match="bn=0|ht=0"):
        add_kelvin_cubit(R=0.1, reduction={"x": "sym_normal"})


def test_reduction_invalid_axis_rejected():
    from add_kelvin import add_kelvin_cubit
    with pytest.raises(ValueError, match="axis"):
        add_kelvin_cubit(R=0.1, reduction={"w": "bn=0"})


def test_reduction_offset_dir_conflict_rejected():
    """If the user sets offset_dir along a reduction axis, the Kelvin
    sphere would poke through the symmetry plane -- must fail loud."""
    from add_kelvin import add_kelvin_cubit
    with pytest.raises(ValueError, match="conflicts with reduction"):
        add_kelvin_cubit(
            R=0.1,
            reduction={"x": "bn=0", "z": "ht=0"},
            offset_dir="x")        # x is a reduction axis


def test_reduction_type_rejected():
    from add_kelvin import add_kelvin_cubit
    with pytest.raises(TypeError, match="must be a dict"):
        add_kelvin_cubit(R=0.1, reduction=["x"])


# ==========================================================
# calc_accel_hdiv.ima_from_mesh_labels
# ==========================================================

def test_ima_auto_quarter_xz_ht_bn():
    """1/4 xz model: sym_ht=0_x (antisym about x=0) + sym_bn=0_z
    (symmetric about z=0) -> Radia IMA '-x+z'."""
    from calc_accel_hdiv import ima_from_mesh_labels
    bnds = ["yoke", "air", "sym_ht=0_x", "sym_bn=0_z"]
    assert ima_from_mesh_labels(bnds) == "-x+z"


def test_ima_auto_full_symmetric_quadrant():
    """1/4 xz model with BOTH faces bn=0 (flux parallel) -> '+x+z'."""
    from calc_accel_hdiv import ima_from_mesh_labels
    assert ima_from_mesh_labels(
        ["sym_bn=0_x", "sym_bn=0_z"]) == "+x+z"


def test_ima_auto_half_x_only():
    """1/2 x model with sym_ht=0_x -> '-x'."""
    from calc_accel_hdiv import ima_from_mesh_labels
    assert ima_from_mesh_labels(["air", "sym_ht=0_x"]) == "-x"


def test_ima_auto_xyz_order_canonical():
    """Output axes appear in x/y/z order regardless of input ordering."""
    from calc_accel_hdiv import ima_from_mesh_labels
    assert ima_from_mesh_labels(
        ["sym_bn=0_z", "sym_ht=0_x", "sym_bn=0_y"]) == "-x+y+z"


def test_ima_auto_empty_when_no_sym():
    from calc_accel_hdiv import ima_from_mesh_labels
    assert ima_from_mesh_labels(["yoke", "air", "kelvin_int"]) == ""
    assert ima_from_mesh_labels([]) == ""


def test_ima_auto_ignores_legacy_axisless():
    """'sym_tangential' / 'sym_normal' have no axis so they cannot
    contribute to an IMA string; they are intentionally ignored here
    (the user should pass --ima explicitly for such legacy samples)."""
    from calc_accel_hdiv import ima_from_mesh_labels
    assert ima_from_mesh_labels(["sym_tangential", "sym_normal"]) == ""


def test_ima_auto_accepts_string_input():
    """Single string is treated as a 1-element list."""
    from calc_accel_hdiv import ima_from_mesh_labels
    assert ima_from_mesh_labels("sym_ht=0_x") == "-x"


def test_ima_auto_unknown_axis_ignored():
    """sym_*_w etc. must not leak into the output."""
    from calc_accel_hdiv import ima_from_mesh_labels
    assert ima_from_mesh_labels(
        ["sym_ht=0_w", "sym_bn=0_x"]) == "+x"
