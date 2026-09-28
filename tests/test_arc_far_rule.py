"""Fixed azimuthal Gauss rule for thick arcs away from the conductor.

The reference is the independent tensor Gauss volume integral of
test_arc_section_regression, so the fixed rule is checked against the
physics and not only against the adaptive rule it replaces.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
import radia as rad

_spec = importlib.util.spec_from_file_location(
    "arc_section_regression", Path(__file__).with_name("test_arc_section_regression.py"))
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
volume_reference = _module.volume_reference

RI, RO, HEIGHT, DENSITY = .005, .040, .105, 5.44e5


@pytest.fixture(autouse=True)
def exact_coordinates():
    rad.FldLenRndSw('off')
    yield
    rad.FldLenRndSw('on')


@pytest.mark.parametrize('angles', [(0., np.pi / 2), (0.3, 1.2)])
@pytest.mark.parametrize('point', [(.12, .05, .02), (-.05, .09, .0), (.02, -.08, .15),
                                   (.30, -.20, -.10), (.06, .07, .06)])
def test_fixed_far_rule_matches_the_volume_integral(angles, point):
    expected = volume_reference(point, RI, RO, HEIGHT, angles, DENSITY, order=64)
    obj = rad.ObjArcCur([0, 0, 0], [RI, RO], list(angles), HEIGHT, 40, 'man', 'z', DENSITY)
    actual = np.array(rad.Fld(obj, 'b', list(point)))
    assert np.linalg.norm(actual - expected) <= 1e-9 * np.linalg.norm(expected)


def test_arc_quadrature_accuracy_is_not_adjustable():
    """The removed PrcArc option is an unknown precision key and fails loudly."""
    with pytest.raises(Exception):
        rad.FldCmpPrc('PrcArc->1e-6')
