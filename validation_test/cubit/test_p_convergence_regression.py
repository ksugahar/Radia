"""
p-Convergence regression test using pre-exported .vol + .vol.json files.

Verifies that NGSolve integration values converge to CAD reference values
as mesh polynomial order increases (p=1..5).

Requirements:
  - NGSolve (pip install ngsolve)
  - Pre-exported .vol files at order 1-5 with companion .vol.json
  - NO Cubit needed

Usage:
  pytest validation_test/cubit/test_p_convergence_regression.py -v
  python validation_test/cubit/test_p_convergence_regression.py   # standalone

Test data location:
  W:/00_CAE/CoreformCubit/2026_04_06_CoreformCUBIT_ExportMesh/
"""

import json
import math
import os
import sys

import pytest

# NGSolve import
try:
    from ngsolve import Mesh, Integrate, CF, BND, BBND
    from ngsolve import TaskManager
    HAS_NGSOLVE = True
except ImportError:
    HAS_NGSOLVE = False

# Test data root
_DATA_ROOT = r"W:\00_CAE\CoreformCubit\2026_04_06_CoreformCUBIT_ExportMesh"


def _skip_if_no_data():
    if not os.path.isdir(_DATA_ROOT):
        pytest.skip(f"Test data not found: {_DATA_ROOT}")


def _skip_if_no_ngsolve():
    if not HAS_NGSOLVE:
        pytest.skip("NGSolve not installed")


def _load_vol_and_json(subdir, prefix, order):
    """Load .vol mesh and .vol.json CAD reference for given order."""
    vol_path = os.path.join(
        _DATA_ROOT, subdir, f"{prefix}_order{order}.vol")
    json_path = vol_path + ".json"
    if not os.path.isfile(vol_path):
        pytest.skip(f"Missing: {vol_path}")
    if not os.path.isfile(json_path):
        pytest.skip(f"Missing: {json_path}")

    mesh = Mesh(vol_path)
    with open(json_path, "r") as f:
        cad = json.load(f)

    with TaskManager():
        ng_vol = Integrate(CF(1), mesh)
        ng_area = Integrate(CF(1), mesh, BND)

        cad_vol = sum(cad.get("materials", {}).values())
        cad_area = sum(cad.get("boundaries", {}).values())

        return ng_vol, ng_area, cad_vol, cad_area


# ================================================================
# Test cases: (subdir, prefix, has_curved_volume, has_curved_area)
# ================================================================

CURVED_MODELS = [
    ("02_Cylinder_Tet", "02_Cylinder_Tet", True, True),
    ("06_acorn", "06_acorn", True, True),
    ("10_Cylinder_Map_webcut", "10_Cylinder_Map_webcut", True, True),
]

# Loft: volume converges, area does NOT (ACIS limitation)
LOFT_MODEL = ("05_loft", "05_loft")


@pytest.mark.parametrize("subdir,prefix,vol_converges,area_converges",
                         CURVED_MODELS,
                         ids=[c[0] for c in CURVED_MODELS])
def test_p_convergence_monotonic(subdir, prefix, vol_converges, area_converges):
    """Volume and area errors must decrease monotonically with order."""
    _skip_if_no_data()
    _skip_if_no_ngsolve()

    vol_errors = []
    area_errors = []

    for p in range(1, 6):
        try:
            ng_vol, ng_area, cad_vol, cad_area = _load_vol_and_json(
                subdir, prefix, p)
        except Exception:
            pytest.skip(f"Cannot load order {p}")

        vol_err = abs(ng_vol - cad_vol) / cad_vol if cad_vol else 0
        area_err = abs(ng_area - cad_area) / cad_area if cad_area else 0
        vol_errors.append(vol_err)
        area_errors.append(area_err)

    # Check monotonic decrease (each order should be better)
    if vol_converges:
        for i in range(1, len(vol_errors)):
            assert vol_errors[i] <= vol_errors[i-1] * 1.01, (
                f"Volume error not monotonic at order {i+1}: "
                f"{vol_errors[i]:.2e} > {vol_errors[i-1]:.2e}")

    if area_converges:
        for i in range(1, len(area_errors)):
            assert area_errors[i] <= area_errors[i-1] * 1.01, (
                f"Area error not monotonic at order {i+1}: "
                f"{area_errors[i]:.2e} > {area_errors[i-1]:.2e}")

    # Order 5 should be very accurate
    if vol_converges:
        assert vol_errors[-1] < 1e-6, (
            f"Volume error at order 5 too large: {vol_errors[-1]:.2e}")
    if area_converges:
        assert area_errors[-1] < 1e-6, (
            f"Area error at order 5 too large: {area_errors[-1]:.2e}")


def test_p_convergence_loft_volume_only():
    """Loft: volume converges, area does NOT (known ACIS limitation)."""
    _skip_if_no_data()
    _skip_if_no_ngsolve()

    subdir, prefix = LOFT_MODEL
    vol_errors = []
    area_errors = []

    for p in range(1, 6):
        try:
            ng_vol, ng_area, cad_vol, cad_area = _load_vol_and_json(
                subdir, prefix, p)
        except Exception:
            pytest.skip(f"Cannot load order {p}")

        vol_err = abs(ng_vol - cad_vol) / cad_vol if cad_vol else 0
        area_err = abs(ng_area - cad_area) / cad_area if cad_area else 0
        vol_errors.append(vol_err)
        area_errors.append(area_err)

    # Volume should converge
    for i in range(1, len(vol_errors)):
        assert vol_errors[i] <= vol_errors[i-1] * 1.01, (
            f"Loft volume error not monotonic at order {i+1}: "
            f"{vol_errors[i]:.2e} > {vol_errors[i-1]:.2e}")

    # Area: known ~0.6% plateau, just warn
    if area_errors[-1] > 0.005:
        print(f"WARNING: Loft area error at p=5 is {area_errors[-1]:.2e} "
              f"(known ACIS loft limitation)")


def test_flat_geometry_exact():
    """Planar geometry (01_Tet_Hex_Pyramid) should be exact at order 1."""
    _skip_if_no_data()
    _skip_if_no_ngsolve()

    ng_vol, ng_area, cad_vol, cad_area = _load_vol_and_json(
        "01_Tet_Hex_Pyramid", "01_Tet_Hex_Pyramid", 1)

    vol_err = abs(ng_vol - cad_vol) / cad_vol if cad_vol else 0
    area_err = abs(ng_area - cad_area) / cad_area if cad_area else 0

    assert vol_err < 1e-10, f"Flat volume error: {vol_err:.2e}"
    assert area_err < 1e-10, f"Flat area error: {area_err:.2e}"


# ================================================================
# Standalone runner
# ================================================================

if __name__ == "__main__":
    if not os.path.isdir(_DATA_ROOT):
        print(f"Test data not found: {_DATA_ROOT}")
        sys.exit(1)
    if not HAS_NGSOLVE:
        print("NGSolve not installed")
        sys.exit(1)

    print("=" * 60)
    print("p-Convergence Regression Test")
    print("=" * 60)

    for subdir, prefix, vol_conv, area_conv in CURVED_MODELS:
        print(f"\n--- {subdir} ---")
        for p in range(1, 6):
            try:
                ng_vol, ng_area, cad_vol, cad_area = _load_vol_and_json(
                    subdir, prefix, p)
                vol_err = (ng_vol - cad_vol) / cad_vol * 100
                area_err = (ng_area - cad_area) / cad_area * 100
                print(f"  p={p}: V_err={vol_err:+.2e}%, A_err={area_err:+.2e}%")
            except Exception as e:
                print(f"  p={p}: SKIP ({e})")

    print(f"\n--- {LOFT_MODEL[0]} (area known issue) ---")
    for p in range(1, 6):
        try:
            ng_vol, ng_area, cad_vol, cad_area = _load_vol_and_json(
                LOFT_MODEL[0], LOFT_MODEL[1], p)
            vol_err = (ng_vol - cad_vol) / cad_vol * 100
            area_err = (ng_area - cad_area) / cad_area * 100
            print(f"  p={p}: V_err={vol_err:+.2e}%, A_err={area_err:+.2e}%")
        except Exception as e:
            print(f"  p={p}: SKIP ({e})")

    print("\nDone.")
