"""Independent NumPy direct Laplace sum for CUDA scalar/dipole FMM."""
from pathlib import Path
import argparse,hashlib,json,socket,time
import numpy as np
from scipy.io import savemat,loadmat

def direct(sources,targets,charge,dipole):
    nc=charge.shape[1];p=np.zeros((len(targets),nc));field=np.zeros((len(targets),nc,3))
    for start in range(0,len(targets),64):
        r=targets[start:start+64,None,:]-sources[None,:,:];r2=np.sum(r*r,axis=2)
        if np.any(r2==0):raise RuntimeError('Coincident reference point')
        inv=1/np.sqrt(r2);inv3=inv/r2
        for c in range(nc):
            d=dipole[:,c*3:c*3+3];dr=np.einsum('tsj,sj->ts',r,d)
            p[start:start+64,c]=np.sum(charge[None,:,c]*inv+dr*inv3,axis=1)
            field[start:start+64,c]=np.sum(charge[None,:,c,None]*r*inv3[:,:,None]+3*r*dr[:,:,None]*(inv3/r2)[:,:,None]-d[None,:,:]*inv3[:,:,None],axis=1)
    return p,field.reshape(len(targets),3*nc)

def main(args):
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[2]
    record=dict(schema='radia.laplace-fmm-numeric.v1',accepted=False,host=socket.gethostname(),timing_accepted=False,
        kernel='1/r and source dipoles; field=-grad potential; no physical scaling',
        scope='numeric point densities only; not singular Galerkin or full FEM-BEM validation',samples=[])
    engine=None
    try:
        receipt=json.loads((args.binary/'build.json').read_text());record['build']=receipt
        if hashlib.sha256((args.binary/'radia_cuda_fmm.mexw64').read_bytes()).hexdigest()!=receipt['binary_sha256']:raise RuntimeError('Binary hash mismatch')
        import matlab.engine
        engine=matlab.engine.start_matlab('-nodesktop -nosplash');engine.addpath(str(Path(__file__).parent),nargout=0)
        rng=np.random.default_rng(19870909);sources=rng.uniform(-1,1,(args.sources,3));targets=rng.uniform(-1.2,1.2,(args.targets,3))
        charge=rng.normal(size=(args.sources,3));dipole=rng.normal(size=(args.sources,9))*.03
        t=time.perf_counter();expectedP,expectedE=direct(sources,targets,charge,dipole);record['cpu_direct_seconds']=time.perf_counter()-t
        for order in args.orders:
            path=out/f'input-{order}.mat';savemat(path,dict(sources=sources,targets=targets,charge=charge,dipole=dipole,order=float(order)))
            output=out/f'output-{order}.mat';engine.fmm_numeric(str(path),str(args.binary.resolve()),str(root/'packages/radia-gpu/matlab'),str(output),nargout=0)
            result=loadmat(output,simplify_cells=True)
            p=float(np.linalg.norm(result['potential']-expectedP)/np.linalg.norm(expectedP));e=float(np.linalg.norm(result['field']-expectedE)/np.linalg.norm(expectedE))
            sample=dict(order=order,potential_relative=p,field_relative=e,report=result['report']);record['samples'].append(sample);print('order',order,'potential',p,'field',e,flush=True)
        if max(record['samples'][-1][k] for k in ('potential_relative','field_relative'))>1e-5:raise RuntimeError('Highest-order direct-sum gate failed')
        record['accepted']=True
    except BaseException as e:record['error']=repr(e);raise
    finally:
        if engine is not None:engine.quit()
        (out/'result.json').write_text(json.dumps(record,indent=2,default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v.item())+'\n',encoding='utf-8')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--sources',type=int,default=512);p.add_argument('--targets',type=int,default=128);p.add_argument('--orders',type=int,nargs='+',default=[4,6,8]);main(p.parse_args())
