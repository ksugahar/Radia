"""Guards: explicit volume topology, contact isolation and material loss accounting.

Assumes: exact rectilinear metal cells, peak phasors, linear nonmagnetic MQS.
These contracts establish no physical accuracy or adopted surrogate baseline.
"""
import importlib.util
from pathlib import Path
import sys
import numpy as np
import pytest

_spec = importlib.util.spec_from_file_location(
    "_volume_peec_contract", Path(__file__).parents[1]/"src/radia/volume_peec.py")
_module = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _module
_spec.loader.exec_module(_module)
build = _module.build_volume_peec


def bar(sigma=2.):
    return build(([0.,1.,3.],[0.,2.],[0.,4.]),
                 np.ones((2,1,1), dtype=bool), sigma,
                 source_faces=[(0,0,0,0,-1)], sink_faces=[(1,0,0,0,1)])


def test_contacts_do_not_merge_interior_nodes_and_kcl_is_exact():
    model = bar()
    assert len(model.nodes_m) == 4
    assert len(model.branches) == 3
    current = np.array([1.,1.,1.])
    divergence = model.incidence() @ current
    expected = np.array([0.,0.,-1.,1.])
    np.testing.assert_array_equal(divergence, expected)
    # Integrated series resistance: length/(sigma*area), independent of split.
    assert model.resistance_ohm.sum() == 3/(2*8)


def test_material_jump_and_cell_losses_preserve_series_halves():
    model = bar(np.array([2.,4.]).reshape(2,1,1))
    current = np.array([1.+2j]*3)
    expected_r = (1/2+2/4)/8
    assert model.resistance_ohm.sum() == expected_r
    np.testing.assert_allclose(model.mean_cell_power_W(current),
                               .5*5*np.array([1/2,2/4])/8, rtol=2e-15)
    assert model.mean_cell_power_W(current).sum() == pytest.approx(
        model.mean_branch_power_W(current).sum(), rel=2e-15)


def test_transverse_branches_exist_without_equipotential_internal_planes():
    model = build(([0.,1.,2.],[0.,1.,2.],[0.,1.,2.]),
                  np.ones((2,2,2), dtype=bool), 1.,
                  source_faces=[(0,0,0,0,-1)], sink_faces=[(1,1,1,0,1)])
    assert len(model.nodes_m) == 10
    assert len(model.branches) == 14
    assert set(model.axes) == {0,1,2}
    assert np.count_nonzero(model.support_cells[:,1] < 0) == 2
    with pytest.raises(ValueError):
        model.nodes_m[0,0] = 7


def test_invalid_geometry_material_contacts_fail_loud():
    mask = np.ones((2,1,1), dtype=bool)
    args = (([0.,1.,3.],[0.,2.],[0.,4.]), mask, 1.)
    with pytest.raises(ValueError, match="exterior"):
        build(*args, source_faces=[(0,0,0,0,1)], sink_faces=[(1,0,0,0,1)])
    with pytest.raises(ValueError, match="exactly one"):
        build(*args, source_faces=[(0,0,0,0,-1)], sink_faces=[(0,0,0,0,-1)])
    with pytest.raises(ValueError, match="conductivity"):
        bar(0.)
    with pytest.raises(ValueError, match="connected"):
        build(([0.,1.,2.,3.],[0.,1.],[0.,1.]),
              np.array([True,False,True]).reshape(3,1,1), 1.,
              source_faces=[(0,0,0,0,-1)], sink_faces=[(2,0,0,0,1)])
    with pytest.raises(ValueError, match="current"):
        bar().mean_cell_power_W([float("nan")]*3)


def test_native_geometry_preserves_support_and_series_resistance():
    model = bar(np.array([2.,4.]).reshape(2,1,1))
    geometry = model.native_geometry()
    np.testing.assert_array_equal(geometry['widths'], [4.,4.,4.])
    np.testing.assert_array_equal(geometry['heights'], [2.,2.,2.])
    np.testing.assert_array_equal(geometry['directions'], [[1.,0.,0.]]*3)
    np.testing.assert_array_equal(geometry['centers'], model.support_boxes_m.mean(axis=1))
    np.testing.assert_allclose(geometry['lengths']/(
        geometry['sigmas']*geometry['widths']*geometry['heights']),
        model.resistance_ohm, rtol=2e-15)
