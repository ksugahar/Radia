"""RadiaField(precision=...) configures Radia's B/A precision instead of failing."""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402

DEFAULT = "PrcB->0.0001,PrcA->0.001"  # radTCompCriterium defaults


@pytest.fixture(autouse=True)
def restore_precision():
    rad.FldLenRndSw("off")
    yield
    rad.FldCmpPrc(DEFAULT)
    rad.FldLenRndSw("on")


def arc():
    rad.UtiDelAll()
    return rad.ObjArcCur([0, 0, 0], [0.035, 0.070], [0.2, 1.8], 0.105, 4, "man", "z", 1e6)


@pytest.mark.parametrize("precision", [1e-4, 1e-6, 1e-9, 2.5e-13])
def test_precision_is_accepted_at_any_positive_magnitude(precision):
    field = rad.RadiaField(arc(), "b", precision=precision)
    assert field.precision == precision
    with ng.TaskManager():
        mesh = ng.Mesh(ng.unit_cube.GenerateMesh(maxh=0.5))
        value = np.asarray(field(mesh(0.5, 0.5, 0.5)))
    assert np.all(np.isfinite(value)) and np.linalg.norm(value) > 0


@pytest.mark.parametrize("field_type", ["b", "h", "a", "m", "phi"])
def test_precision_works_for_every_field_type(field_type):
    rad.RadiaField(arc(), field_type, precision=1e-7)


def test_precision_does_not_change_the_arc_field():
    with ng.TaskManager():
        mesh = ng.Mesh(ng.unit_cube.GenerateMesh(maxh=0.5))
        point = mesh(0.2, 0.3, 0.4)
        plain = np.asarray(rad.RadiaField(arc(), "b")(point))
        rad.FldCmpPrc(DEFAULT)
        tight = np.asarray(rad.RadiaField(arc(), "b", precision=1e-9)(point))
    np.testing.assert_array_equal(plain, tight)


@pytest.mark.parametrize("precision", [0.0, -1e-6, float("nan"), float("inf")])
def test_invalid_precision_fails_before_touching_radia(precision):
    with pytest.raises(ValueError, match="precision"):
        rad.RadiaField(arc(), "b", precision=precision)
