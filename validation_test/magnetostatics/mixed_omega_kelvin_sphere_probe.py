"""Uniform applied field / permeable sphere: Kelvin response-only experiment.

An independently authored analytic case. The unbounded uniform applied field
is kept analytic; only the decaying reaction is represented in the Kelvin ball.
This experiment must pass before claiming this source split is supported.
"""
import argparse
import json
from pathlib import Path
import sys
import time
import types


def main():
    parser=argparse.ArgumentParser(__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--maxh',type=float,default=.35)
    parser.add_argument('--order',type=int,choices=(1,2,3),default=1)
    parser.add_argument('--geometry-order',type=int,choices=(1,2,3),default=1)
    args=parser.parse_args()
    import numpy as np
    import ngsolve as ng
    from netgen.occ import Sphere,Pnt,Glue,OCCGeometry
    from netgen.meshing import IdentificationType
    root=Path(__file__).resolve().parents[2]/'src'/'radia'
    namespace=types.ModuleType('radia');namespace.__path__=[str(root)]
    sys.modules['radia']=namespace
    from radia.kelvin_material import make_kelvin_mu_cf, MU_0
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin as solve
    ng.SetNumThreads(1)
    start=time.perf_counter()
    record=dict(schema='mixed-omega.kelvin-response-sphere.v1',passed=False,
                mu_r=100.,sphere_radius=.3,kelvin_radius=1.,ngsolve=ng.__version__,
                order=args.order,geometry_order=args.geometry_order)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    try:
        iron=Sphere(Pnt(0,0,0),.3);iron.mat('iron');iron.faces.name='interface';iron.maxh=args.maxh/2
        outer=Sphere(Pnt(0,0,0),1.);outer.faces.name='kelvin_int'
        air=outer-iron;air.mat('air')
        exterior=Sphere(Pnt(3,0,0),1.);exterior.mat('kelvin');exterior.faces.name='kelvin_ext'
        for vertex in exterior.vertices:
            if vertex.p[2]>0:vertex.name='GND'
        geo=Glue([iron,air,exterior])
        faces={face.name:face for solid in geo.solids for face in solid.faces if face.name}
        faces['kelvin_int'].Identify(faces['kelvin_ext'],'periodic',IdentificationType.PERIODIC)
        with ng.TaskManager():
            mesh=ng.Mesh(OCCGeometry(geo).GenerateMesh(maxh=args.maxh,grading=.6))
            mesh.Curve(args.geometry_order)
            space=ng.H1(mesh,order=args.order)
            source=ng.GridFunction(space);source.Set(-ng.z)
            zero=ng.GridFunction(space)
            mu=make_kelvin_mu_cf(mesh,1.,(3,0,0),mu_r_by_material={'iron':100.})
            result=solve(mesh,ng.CF((0,0,1)),source,1.,(3,0,0),mu_cf=mu,
                reduced_materials=('air',),total_materials=('iron','kelvin'),
                interface_boundary='interface',kelvin_interface_boundary='kelvin_int',
                reduced_normal_flux=MU_0*ng.specialcf.normal(3)[2],
                reduced_flux_boundary='kelvin_int',
                kelvin_source_potential=zero,order=args.order)
        points=[(0,0,0),(.1,0,0),(0,.1,.1),(.55,.1,.1),(.1,.6,.15)]
        actual=[];exact=[]
        alpha=99/102
        for point in points:
            actual.append(np.array(result['H_cf'](mesh(*point))).tolist())
            r=np.asarray(point,dtype=float);rho=np.linalg.norm(r)
            expected=np.array([0.,0.,3/102]) if rho<.3 else np.array([0.,0.,1.])+alpha*.3**3*(3*r*r[2]/rho**5-np.array([0.,0.,1.])/rho**3)
            exact.append(expected.tolist())
        error=float(np.linalg.norm(np.array(actual)-exact)/np.linalg.norm(exact))
        record.update(elements=mesh.ne,actual_H=actual,exact_H=exact,points=points,
            relative_rms=error,residual=result.get('linear_residual'),elapsed_seconds=time.perf_counter()-start,
            acceptance_relative_rms=.01,passed=error<.01)
        assert record['passed'],f'Kelvin sphere error {error}'
    except Exception as error:
        record['failure']=repr(error)
        raise
    finally:
        args.out.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
