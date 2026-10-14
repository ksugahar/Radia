"""
Simple Radia Test - Basic functionality verification

Tests:
  - Module import
  - Version query
  - Rectangular magnet creation
  - Magnetization setting
  - Field computation
"""

import pytest


@pytest.mark.basic
def test_import_radia():
    """Verify radia module can be imported."""
    import radia as rad

    assert hasattr(rad, "__version__")
    assert hasattr(rad, "ObjRecMag")


@pytest.mark.basic
def test_version(radia_module):
    """Verify version is returned."""
    version = radia_module.UtiVer()
    assert version is not None


@pytest.mark.basic
def test_rectangular_magnet(radia_clean):
    """Create a rectangular magnet and verify field computation."""
    rad = radia_clean

    # Create 10x10x10 mm magnet with magnetization (1000 A/m in Z)
    magnet = rad.magnet_box([0, 0, 0], [10, 10, 10], [0, 0, 1000])
    assert isinstance(magnet, int)
    assert magnet > 0

    # Compute field at (0, 0, 20) mm
    field = rad.Fld(magnet, "b", [0, 0, 20])
    assert len(field) == 3

    # Bz should be positive (above magnet magnetized in +z)
    assert field[2] > 0

    # Bx, By should be ~0 by symmetry
    assert abs(field[0]) < 1e-10
    assert abs(field[1]) < 1e-10
