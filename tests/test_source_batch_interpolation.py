"""An explicitly selected native batch evaluator preserves source coefficients."""
import numpy as np
import pytest
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
from radia.kelvin_solver import interpolate_source_field


def test_explicit_batch_source_matches_coefficient_evaluation():
    mesh=MakeStructured3DMesh(nx=1,ny=1,nz=1)
    source=ng.CF((ng.x+2*ng.y,ng.y-ng.z,ng.z+.3))
    calls=[]
    def evaluate(points):
        calls.append(points.copy())
        x,y,z=points.T
        return np.column_stack((x+2*y,y-z,z+.3))
    with ng.TaskManager():
        direct=interpolate_source_field(mesh,source,mesh.Materials('.*'),order=2)
        batch=interpolate_source_field(mesh,source,mesh.Materials('.*'),order=2,point_evaluator=evaluate)
    assert len(calls)==1 and len(calls[0])==batch['nodes']
    for point in ((.12,.17,.23),(.71,.42,.35)):
        np.testing.assert_allclose(batch['field'](mesh(*point)),direct['field'](mesh(*point)),atol=1e-13,rtol=1e-13)


@pytest.mark.parametrize('invalid',[lambda p:np.zeros(len(p)),lambda p:np.full((len(p),3),np.nan),lambda p:np.ones((len(p),3),dtype=complex)])
def test_invalid_batch_source_fails(invalid):
    mesh=MakeStructured3DMesh(nx=1,ny=1,nz=1)
    with ng.TaskManager(),pytest.raises(ValueError,match='point_evaluator'):
        interpolate_source_field(mesh,ng.CF((0,0,1)),mesh.Materials('.*'),order=1,point_evaluator=invalid)
