import sys,json,hashlib,argparse
from pathlib import Path
from dataclasses import asdict
from unittest.mock import patch
import numpy as np
import ngsolve as ng
import radia
root=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser()
parser.add_argument('--thermal-order',type=int,choices=(1,2),default=1)
parser.add_argument('--output',type=Path,default=Path('C:/temp/ih-axisym-native/axisymmetric'))
args=parser.parse_args()
radia.__path__.insert(0,str(root/'src/radia'))
import radia.simulink
radia.simulink.__path__.insert(0,str(root/'src/radia/simulink'))
from radia.simulink.ih_operator_assembly import IHOperatorAssemblyOptions,_assemble_thermal_operators
mesh=ng.Mesh(str(root/'validation_test/panels/fixtures/heat_workpiece_cylinder_R25_H25_axisym.vol'))
for i,name in enumerate(mesh.GetBoundaries()):
    if name!='axis':mesh.ngmesh.SetBCName(i,'sibc')
opts=IHOperatorAssemblyOptions(axisymmetric_thermal_vol='fixture.vol',sample_time_s=.1,thermal_order=args.thermal_order)
with ng.TaskManager(),patch.object(ng,'Mesh',return_value=mesh):
    op=_assemble_thermal_operators(Path('fixture.vol'),np.ones(mesh.nv)*125,opts)
cfg=dict(schema='radia.ih.simulink.native_sfunction.v1',backend='matlab-level2+radia-mex-handles',
    python_fallback=False,eddy_solver='fem',thermal_solver='fem',bh_mode='linear',
    n_eddy_unknown=1,eddy_matrix_real=[1.],eddy_matrix_imag=[0.],eddy_rhs_real=[1.],eddy_rhs_imag=[0.],
    sample_time_s=.1,thermal_tolerance=1e-13,thermal_max_iterations=2000,
    convection_W_per_m2K=10.,rotation_mode='none',angle_origin_rad=0.)
fields=asdict(op)
for key in ('n_temperature','mass_row_ptr','mass_col','mass_value','stiffness_row_ptr','stiffness_col','stiffness_value','convection_row_ptr','convection_col','convection_value','heat_to_temperature_projection'):cfg[key]=fields[key]
cfg.update(n_heat=len(op.heat_dofs),heat_projection=op.unit_heat_density_W_per_m3.tolist(),heat_cell_weights=op.heat_cell_weights_m3.tolist(),temperature_cell_weights=op.temperature_cell_weights_J_per_K.tolist(),initial_temperature_K=[293.15]*op.n_temperature)
if op.constant_coefficients is not None:
    cfg.update(temperature_representation='ngsolve-h1-coefficients',
        temperature_constant_coefficients=op.constant_coefficients.tolist(),
        initial_temperature_K=(293.15*op.constant_coefficients).tolist(),
        temperature_evaluation=op.temperature_evaluation)
with ng.TaskManager():
    fes=ng.H1(mesh,order=args.thermal_order);u,v=fes.TnT();w=2*np.pi*ng.x
    m=ng.BilinearForm(fes);m+=w*7800*467*u*v*ng.dx;m.Assemble()
    a=ng.BilinearForm(fes);a+=w*(7800*467*u*v+.1*46.6*ng.grad(u)*ng.grad(v))*ng.dx+.1*w*10*u*v*ng.ds('sibc');a.Assemble()
    f=ng.LinearForm(fes);f+=w*(125+10*293.15)*v*ng.ds('sibc');f.Assemble()
    inverse=a.mat.Inverse(fes.FreeDofs(),inverse='sparsecholesky');temp=ng.GridFunction(fes);temp.Set(293.15);history=[]
    for i in range(100):
        b=temp.vec.CreateVector();b.data=m.mat*temp.vec+.1*f.vec;temp.vec.data=inverse*b;history.append(temp.vec.FV().NumPy().tolist())
out=args.output;out.mkdir(parents=True,exist_ok=True)
(out/'chain.json').write_text(json.dumps(dict(config=cfg,reference_K=history,current_A=1.,unit_power_W=op.heat_power_W,assembled_power_W=op.heat_power_W,nodes=mesh.nv,source_sha256={
    'ih_operator_assembly':hashlib.sha256((root/'src/radia/simulink/ih_operator_assembly.py').read_bytes()).hexdigest(),
    'driver':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},scope='axisymmetric production thermal operators; uniform input; in-memory boundary relabeling, not production VOL gate')))
