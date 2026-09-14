import importlib.util
import json
import hashlib
import platform
from pathlib import Path
import numpy as np
import ngsolve as ng

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
ROOT = args.input
OUT = args.out
OUT.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location('avg', ROOT/'05_qsurf_chain/axisymmetrize.py')
avg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(avg)
ng.SetNumThreads(4)
result = {'averaging': {}, 'initialization': {}, 'em': {},
          'status': 'diagnostic; EM independent local-field validation incomplete',
          'runtime': {'host': platform.node(), 'ngsolve': ng.__version__},
          'input_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [ROOT/'05_qsurf_chain/axisymmetrize.py',
                                    *sorted((ROOT/'data').glob('*.sol')),
                                    ROOT/'data/work_occ_h0.0018.vol']}}

with ng.TaskManager():
    mesh = ng.Mesh(str(ROOT/'data/work_occ_h0.0018.vol'))
    face = np.array([e.index for e in mesh.Elements(ng.BND)])
    centre, area, verts, region = avg.element_data(mesh)
    fes = ng.H1(mesh, order=1)
    def load(path):
        g=ng.GridFunction(fes);g.Load(str(path));return g
    def face_power(g):
        e=np.array(ng.Integrate(g,mesh,ng.BND,element_wise=True))
        return np.bincount(face,weights=e,minlength=27)
    def rel(a,b):
        return float(np.sqrt(ng.Integrate((a-b)**2,mesh,ng.BND)/ng.Integrate(b*b,mesh,ng.BND)))
    for case,rawfile in [('A','qsurf_caseA_raw_10645A.sol'),('B','qsurf_caseB_raw_9757A.sol')]:
        path=ROOT/'data'/rawfile
        raw=load(path);p0=face_power(raw);rows=[];prev=None
        for width in [2.,1.,.5,.25,.125]:
            dest=OUT/f'{case}_{width}.sol'
            info=avg.axisymmetrize(str(path),str(dest),bin_mm=width)
            g=load(dest);p=face_power(g)
            again=OUT/f'{case}_{width}_twice.sol'
            avg.axisymmetrize(str(dest),str(again),bin_mm=width)
            twice=load(again)
            row={'bin_mm':width,'cells':info['n_bins'],'total_relative_error':float(p.sum()/p0.sum()-1),
                 'inner_percent':float(100*p[19]/p.sum()),
                 'max_face_power_change_percentage_points':float(100*np.max(abs(p-p0))/p0.sum()),
                 'repeat_relative_surface_L2':rel(twice,g),
                 'change_from_previous_width_relative_surface_L2':None if prev is None else rel(g,prev)}
            rows.append(row);prev=g
            print(case,row,flush=True)
        result['averaging'][case]={'raw_inner_percent':float(100*p0[19]/p0.sum()),'rows':rows}
        em=json.loads((ROOT/f'reference/em_case{case}.json').read_text())
        port=.5*em['current_A']**2*em['delta_R_mOhm']*1e-3
        result['em'][case]={'reported_workpiece_power_W':em['P_wp_W'],
                           'port_increment_power_W':port,'ratio':em['P_wp_W']/port,
                           'esim_converged':em['esim_converged'],
                           'status':'inconsistency requiring EM rerun; not independent local-field validation'}
    # Homogeneous difference equation exactly isolates initial-condition error
    # for the received linear constant-property thermal model, independent of q.
    for order in [1,2,3]:
        f=ng.H1(mesh,order=order);u,v=f.TnT()
        good=ng.GridFunction(f);good.Set(ng.CF(25))
        bad=ng.GridFunction(f);bad.vec[:]=25
        err=ng.GridFunction(f);err.vec.data=bad.vec-good.vec
        volume=ng.Integrate(ng.CF(1),mesh)
        def stats():
            return {'mean_error_C':float(ng.Integrate(err,mesh)/volume),
                    'rms_error_C':float(np.sqrt(ng.Integrate(err*err,mesh)/volume)),
                    'vertex_max_abs_error_C':float(np.max(abs(err.vec.FV().NumPy()[:mesh.nv])))}
        row={'initial':stats()}
        if order>1:
            a=ng.BilinearForm(f,symmetric=True)
            a+=46.6*ng.grad(u)*ng.grad(v)*ng.dx+10*u*v*ng.ds
            a.Assemble()
            m=ng.BilinearForm(f,symmetric=True);m+=7800*467*u*v*ng.dx;m.Assemble()
            mat=m.mat.CreateMatrix();mat.AsVector().data=m.mat.AsVector()+.1*a.mat.AsVector()
            inv=mat.Inverse(f.FreeDofs(),inverse='sparsecholesky')
            rhs=err.vec.CreateVector()
            for _ in range(50):
                rhs.data=m.mat*err.vec
                err.vec.data=inv*rhs
        row['after_5s']=stats()
        result['initialization'][str(order)]=row
        print('INITIALIZATION',order,row,flush=True)
        (OUT/'result.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result['em'],indent=2),flush=True)
