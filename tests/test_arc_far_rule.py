"""Fixed azimuthal Gauss rule for thick arcs away from the conductor.

The reference is the independent tensor Gauss volume integral of
test_arc_section_regression, so the fixed rule is checked against the
physics and not only against the adaptive rule it replaces.  A second check
compares both rules in separate processes (the switch is read once).
"""
import importlib.util
import json
import os
import subprocess
import sys
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


_CHILD = r"""
import json, sys
import numpy as np
import radia as rad
rad.FldLenRndSw('off')
rng = np.random.default_rng(3)
points = rng.uniform([-.3, -.3, -.3], [.3, .3, .3], size=(400, 3))
obj = rad.ObjArcCur([0, 0, 0], [.005, .040], [0.0, 1.5707963267948966], .105, 40, 'man', 'z', 5.44e5)
print(json.dumps(np.asarray(rad.Fld(obj, 'b', points.tolist())).tolist()))
"""


def _child_field(far_rule):
    env = dict(os.environ)
    env.pop("RADIA_ARC_FAR_RULE", None)
    if not far_rule:
        env["RADIA_ARC_FAR_RULE"] = "0"
    out = subprocess.run([sys.executable, "-c", _CHILD], env=env, check=True,
                         capture_output=True, text=True).stdout
    return np.asarray(json.loads(out.strip().splitlines()[-1]))


def test_far_rule_agrees_with_the_adaptive_rule_everywhere():
    fast, adaptive = _child_field(True), _child_field(False)
    scale = np.linalg.norm(adaptive, axis=1).max()
    assert np.linalg.norm(fast - adaptive, axis=1).max() <= 1e-8 * scale
