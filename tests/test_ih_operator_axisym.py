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
