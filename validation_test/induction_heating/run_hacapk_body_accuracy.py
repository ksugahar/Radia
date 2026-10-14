"""Isolate ACA map error from iterative stopping error on a self-authored ring.

This small validation temporarily retains the exact construction matrices.
Every compressed/direct product uses the SAME assembled entries. It does
not measure production memory or provide an induced-norm error theorem.
"""
from __future__ import annotations
import argparse,hashlib,json,platform
from pathlib import Path
import ngsolve as ng
import numpy as np
from radia.bem_sibc_solver import ScalarBIESIBCSolver
from radia.bem_loop_work import _body_layer
from run_loop_work_ring import ring_mesh


def run(level,epsilons,threads):
    mesh,points,triangles=ring_mesh(level)
    if len(points)>7000 or len(triangles)>14000:
        raise ValueError('The retained-entry validation keeps the construction guard')
    theta=np.arctan2(points[:,1],points[:,0])
    phi=np.arctan2(points[:,2],np.linalg.norm(points[:,:2],axis=1)-.03)
    rng=np.random.default_rng(271828)
    probes=dict(constant=np.ones(len(points)),x=points[:,0],z=points[:,2],
        theta7_phi3=np.cos(7*theta)*np.cos(3*phi),
        random_real=rng.normal(size=len(points)),
        random_complex=rng.normal(size=len(points))+1j*rng.normal(size=len(points)))
    probes={key:vector/np.linalg.norm(vector) for key,vector in probes.items()}
    records=[]
    with ng.TaskManager():
        for epsilon in epsilons:
            solver=ScalarBIESIBCSolver(mesh,order=1,assemble_dense=True,
                use_intree_bem=True,use_intree_hacapk=True,hacapk_aca_eps=epsilon,
                intree_geom_order=1,intree_singular_n_q=6,intree_regular_quad_degree=7)
            solver._loop_body_backend='hacapk'
            if solver.ndof!=len(points):raise RuntimeError('Unexpected nodal ordering')
            layers={}
            for name in ('SL','DL'):
                exact=getattr(solver,name)
                handle=getattr(solver,'_'+name+'_hacapk')
                handle.ReleaseDenseEntries()
                compressed=_body_layer(solver,name)
                errors={}
                for label,vector in probes.items():
                    for transpose in (False,True):
                        reference=(exact.T if transpose else exact)@vector
                        actual=(compressed.T if transpose else compressed)@vector
                        relative=np.linalg.norm(actual-reference)/max(np.linalg.norm(reference),np.finfo(float).tiny)
                        if not np.isfinite(relative):raise RuntimeError('Nonfinite map error')
                        errors[label+('_transpose' if transpose else '_forward')]=float(relative)
                layers[name]=dict(errors=errors,max_relative=max(errors.values()),stats=handle.GetStats())
            records.append(dict(controls=solver.hacapk_controls.copy(),layers=layers))
    return dict(schema='radia.hacapk-body-accuracy.v1',faces=len(triangles),vertices=len(points),
        geometry_sha256=hashlib.sha256(points.tobytes()+triangles.tobytes()).hexdigest(),
        os=platform.system(),threads=threads,ngsolve_version=ng.__version__,
        regular_quad_degree=7,singular_n_q=6,probe_seed=271828,
        norm='Euclidean relative action norm, normalized deterministic vectors; not an induced norm bound',
        probes=list(probes),records=records)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--level',type=int,default=64)
    parser.add_argument('--eps',type=float,nargs='+',default=[1e-10,1e-12])
    parser.add_argument('--threads',type=int,default=1)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.threads<1 or any(not np.isfinite(e) or e<=0 for e in args.eps):
        parser.error('Threads and finite ACA tolerances must be positive')
    ng.SetNumThreads(args.threads)
    result=run(args.level,args.eps,args.threads)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({r['controls']['aca_eps']:{k:v['max_relative'] for k,v in r['layers'].items()}
                      for r in result['records']}),flush=True)


if __name__=='__main__':main()
