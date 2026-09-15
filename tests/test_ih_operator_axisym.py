from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
import ngsolve as ng
from radia.simulink.ih_operator_assembly import IHOperatorAssemblyOptions, _assemble_thermal_operators


def fixture():
    path = Path(__file__).parents[1]/'validation_test/panels/fixtures/heat_workpiece_cylinder_R25_H25_axisym.vol'
    mesh=ng.Mesh(str(path))
    for i,name in enumerate(mesh.GetBoundaries()):
        if name != 'axis': mesh.ngmesh.SetBCName(i,'sibc')
    return mesh


def test_axisym_operators_revolved_mass_and_surface_power():
    mesh=fixture()
    options=IHOperatorAssemblyOptions(axisymmetric_thermal_vol='selected.vol')
    with patch.object(ng,'Mesh',return_value=mesh), ng.TaskManager():
        result=_assemble_thermal_operators(Path('selected.vol'),np.ones(mesh.nv)*125,options)
    volume=np.pi*.025**2*.025
    area=2*np.pi*.025*.025+2*np.pi*.025**2
    assert sum(result.mass_value)==pytest.approx(7800*467*volume,rel=1e-12)
    assert result.heat_power_W==pytest.approx(125*area,rel=1e-12)
    assert sum(result.convection_value)==pytest.approx(area,rel=1e-12)
    assert abs(sum(result.stiffness_value))<1e-10


def test_axis_role_is_rejected_even_when_label_is_shared():
    mesh=fixture();mesh.ngmesh.SetBCName(3,'sibc')
    with patch.object(ng,'Mesh',return_value=mesh), ng.TaskManager():
        with pytest.raises(ValueError,match='r=0'):
            _assemble_thermal_operators(Path('selected.vol'),np.ones(mesh.nv),
                IHOperatorAssemblyOptions(axisymmetric_thermal_vol='selected.vol'))


def test_p2_constant_and_mixed_surface_load():
    mesh = fixture()
    options = IHOperatorAssemblyOptions(axisymmetric_thermal_vol='selected.vol', thermal_order=2).checked()
    with patch.object(ng, 'Mesh', return_value=mesh), ng.TaskManager():
        result = _assemble_thermal_operators(Path('selected.vol'), np.ones(mesh.nv)*125, options)
    c = result.constant_coefficients
    assert len(c) > mesh.nv
    assert np.linalg.norm(c[mesh.nv:]) < 1e-10
    assert np.dot(c, result.temperature_cell_weights_J_per_K) == pytest.approx(7800*467*np.pi*.025**2*.025)
    load = np.array(result.heat_to_temperature_projection).reshape(len(c), -1) @ result.unit_heat_density_W_per_m3
    assert np.dot(c, load) == pytest.approx(result.heat_power_W, rel=1e-12)
    evaluation = result.temperature_evaluation
    values = np.zeros(evaluation['n_samples'])
    np.add.at(values, evaluation['rows'], np.array(evaluation['values'])*c[evaluation['cols']])
    assert np.max(np.abs(values-1)) < 1e-11
    from scipy.sparse import csr_matrix
    stiffness = csr_matrix((result.stiffness_value,result.stiffness_col,result.stiffness_row_ptr))
    convection = csr_matrix((result.convection_value,result.convection_col,result.convection_row_ptr))
    steady = c*(293.15+125/10)
    residual = (stiffness+10*convection)@steady-load-10*293.15*(convection@c)
    assert np.linalg.norm(residual) < 1e-9
