"""Independent ideal-unit-sphere mean/core magnetization check, CPU FE evaluation only."""
from pathlib import Path
import argparse, hashlib, json
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=Path,required=True)
    parser.add_argument('--cpu-h5',type=Path,required=True)
    parser.add_argument('--gpu-result',nargs=2,action='append',default=[],metavar=('ROLE','JSON'))
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False,parents=True)
    import ngsolve as ng
    from scipy.optimize import brentq
    from radia.vim._nonlinear import _bh_inverse_funcs
    from hdiv_hdf5 import read_hdiv_h5
    ng.SetNumThreads(1)
    result=dict(accepted=False,scope='Ideal unit-sphere mean and declared core M stencil; no whole-field or operator-residual certification',units='M and H A/m')
    try:
        # Independently authored public unit-sphere fixture; reject unrelated inputs.
        mesh_sha=hashlib.sha256(args.mesh.read_bytes()).hexdigest()
        if mesh_sha not in {
            '74da487538d92c370393f149cbff5fed0b3e89a28af8382e7fbb42e01cfacf77',
            '681e3d16060a9ecb34d968923c960104f3eac6986848b57275beda80b150b6cb',
        }:
            raise ValueError('This physical gate is scoped to the declared public P2 unit sphere')
        ref,_=read_hdiv_h5(args.cpu_h5);field=ref['applied_field'].ravel();strength=np.linalg.norm(field)
        if strength<=0:raise ValueError('Nonzero uniform applied field required')
        curves,_,_=_bh_inverse_funcs(ref['table'][:,0],ref['table'][:,1])
        mag=brentq(lambda value:float(curves(np.array([value]))[0][0])*value+value/3-strength,0,3*strength,xtol=1e-8)
        expected=field/strength*mag
        points=np.array([[x+.013,y-.009,z+.017] for x in [-.5,-.25,0,.25,.5] for y in [-.5,-.25,0,.25,.5] for z in [-.5,-.25,0,.25,.5]])
        points=points[np.linalg.norm(points,axis=1)<=.65]
        cases=[('cpu_reference',ref['coefficients'].ravel(),ref['mean'].ravel())]
        for role,path in args.gpu_result:
            gpu=json.loads(Path(path).read_text())
            if not gpu['accepted']:raise ValueError('Unaccepted GPU result input')
            s=gpu['solution'];cases.append((role,np.asarray(s['coefficients']),np.asarray(s['mean_magnetization'])))
        result.update(mesh_sha256=mesh_sha,expected_uniform_magnetization=expected.tolist(),core_points=points.tolist(),relative_limit=.01,cases={})
        with ng.TaskManager():
            mesh=ng.Mesh(str(args.mesh));fes=ng.HDiv(mesh,order=2);gf=ng.GridFunction(fes)
            for role,m,mean in cases:
                if m.size!=fes.ndof:raise ValueError('Coefficient/FE dimension mismatch')
                gf.vec.FV().NumPy()[:]=m
                values=np.array([tuple(gf(mesh(*point))) for point in points])
                core=float(np.sqrt(np.mean(np.sum((values-expected)**2,axis=1)))/np.linalg.norm(expected))
                avg=float(np.linalg.norm(mean-expected)/np.linalg.norm(expected))
                result['cases'][role]=dict(core_relative_rms=core,mean_relative_error=avg)
                if not np.isfinite(core+avg) or max(core,avg)>.01:raise ValueError('Ideal sphere physical mean/core gate failed')
        result['accepted']=True
    except BaseException as ex:
        result['error']=repr(ex);raise
    finally:
        (args.out/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
