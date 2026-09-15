"""Real BEM -> production thermal assembly -> independent NGSolve transient.

In-memory coarse volume fixture; deliberately not a CAD/VOL export gate.
The axisymmetric geometry is represented in 3D because the current native
operator assembler accepts only 3D volumes. No claim of mesh convergence.
"""
import argparse
import json
import hashlib
import sys
from pathlib import Path
from unittest.mock import patch
from dataclasses import asdict
import numpy as np
import ngsolve as ng
from netgen.occ import Cylinder, Pnt, Vec, OCCGeometry
import radia

ROOT = Path(__file__).resolve().parents[2]
radia.__path__.insert(0, str(ROOT/'src/radia'))
import radia.simulink
radia.simulink.__path__.insert(0, str(ROOT/'src/radia/simulink'))
from radia.simulink.ih_operator_assembly import IHOperatorAssemblyOptions, _assemble_thermal_operators
sys.path.insert(0, str(ROOT/'src/radia/panels'))
import calc_inductance as calc


def run(out, bore):
    out.mkdir(parents=True, exist_ok=True)
    body = Cylinder(Pnt(0,0,-.0125),Vec(0,0,1),r=.025,h=.025)
    if bore:
        body = body-Cylinder(Pnt(0,0,-.013),Vec(0,0,1),r=bore,h=.026)
    body.mat('workpiece')
    for face in body.faces: face.name='sibc'
    mesh = ng.Mesh(OCCGeometry(body).GenerateMesh(maxh=.0045))
    label = str(out/'in_memory.vol')
    original = ng.Mesh
    def load(value,*args,**kw):
        return mesh if isinstance(value,str) and value == label else original(value,*args,**kw)
    args = calc.build_argparser().parse_args(['--coil-solver','peec','--coil-step','unused.step',
        '--vol',label,'--wp-label','sibc','--frequency','1000','--current','1',
        '--sigma','5.8e7','--mu-r','1','--h1-order','1','--wp-bem-backend','intree-dense',
        '--msh-output',str(out/'heat.msh')])
    theta=np.linspace(0,2*np.pi,721)
    points=np.column_stack([.03*np.cos(theta),.03*np.sin(theta),theta*0])
    source=dict(source_type='filament',paths=[np.stack([points[:-1],points[1:]],axis=1)],I_fil=np.array([1+0j]))
    opts=IHOperatorAssemblyOptions(sample_time_s=.1,workpiece_relative_permeability=1,
        workpiece_conductivity_S_per_m=5.8e7,frequency_hz=1000)
    with patch.object(ng,'Mesh',side_effect=load):
        result=calc._solve_workpiece_weak_coupled(args,source)
        fes=ng.H1(mesh,order=1)
        flux=ng.GridFunction(fes); flux.Load(result['qsurf_sol'])
        thermal=_assemble_thermal_operators(Path(label),flux.vec.FV().NumPy().copy(),opts)
    assert abs(thermal.heat_power_W/result['P_wp']-1)<1e-8
    fields=asdict(thermal)
    cfg={key:fields[key] for key in ('n_temperature','mass_row_ptr','mass_col','mass_value',
        'stiffness_row_ptr','stiffness_col','stiffness_value','convection_row_ptr','convection_col',
        'convection_value','heat_to_temperature_projection')}
    cfg.update(schema='radia.ih.simulink.native_sfunction.v1',backend='matlab-level2+radia-mex-handles',
        python_fallback=False,eddy_solver='peec',thermal_solver='fem',bh_mode='linear',
        n_eddy_unknown=1,n_heat=len(thermal.heat_dofs),eddy_matrix_real=[1.],eddy_matrix_imag=[0.],
        eddy_rhs_real=[1.],eddy_rhs_imag=[0.],heat_projection=thermal.unit_heat_density_W_per_m3.tolist(),
        heat_cell_weights=thermal.heat_cell_weights_m3.tolist(),
        temperature_cell_weights=thermal.temperature_cell_weights_J_per_K.tolist(),
        initial_temperature_K=[293.15]*fes.ndof,sample_time_s=.1,thermal_tolerance=1e-13,
        thermal_max_iterations=2000,convection_W_per_m2K=10.,rotation_mode='none',angle_origin_rad=0.)
    # Independent forms/load: do not rebuild reference matrices from exported CSR.
    u,v=fes.TnT(); mass=ng.BilinearForm(fes); mass+=7800*467*u*v*ng.dx; mass.Assemble()
    system=ng.BilinearForm(fes)
    system+=7800*467*u*v*ng.dx+.1*46.6*ng.grad(u)*ng.grad(v)*ng.dx+.1*10*u*v*ng.ds
    system.Assemble(); inverse=system.mat.Inverse(fes.FreeDofs(),inverse='sparsecholesky')
    current=100.; loadform=ng.LinearForm(fes)
    loadform+=(current**2*flux+10*293.15)*v*ng.ds; loadform.Assemble()
    temp=ng.GridFunction(fes);temp.Set(293.15)
    history=[]
    for step in range(100):
        rhs=temp.vec.CreateVector();rhs.data=mass.mat*temp.vec+.1*loadform.vec
        temp.vec.data=inverse*rhs;history.append(temp.vec.FV().NumPy().tolist())
    payload=dict(config=cfg,reference_K=history,current_A=current,unit_power_W=result['P_wp'],
        assembled_power_W=thermal.heat_power_W,bore_m=bore,nodes=fes.ndof,
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__), ROOT/'src/radia/panels/calc_inductance.py',
                      ROOT/'src/radia/simulink/ih_operator_assembly.py',
                      ROOT/'src/radia/bem_sibc_solver.py', ROOT/'src/radia/bem_loop_extension.py',
                      ROOT/'src/radia/workpiece_surface.py']},
        scope='coarse in-memory volume; real weak BEM and production thermal assembly; no CAD/VOL gate')
    (out/'chain.json').write_text(json.dumps(payload),encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();ng.SetNumThreads(4)
    with ng.TaskManager():
        for bore,name in [(0.,'solid'),(.0075,'bored')]: run(args.output/name,bore)
