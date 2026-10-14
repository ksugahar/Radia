"""Independent native NGSolve Galerkin SL/DL oracle for corrected CUDA FMM."""
from pathlib import Path
import argparse,json,hashlib,os,sys,socket,time
import numpy as np
import ngsolve as ng
from ngsolve.bem import LaplaceSL,LaplaceDL
from scipy.io import loadmat
from scipy.sparse import coo_matrix

def main(args):
    root=Path(__file__).resolve().parents[2];out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    record=dict(schema='radia.native-boundary-fmm.v1',accepted=False,host=socket.gethostname(),timing_accepted=False,
        scope='Flat exterior P1 Galerkin SL/DL with native exact near corrections; independent NGSolve quadrature')
    engine=None
    try:
        import netgen,matlab.engine
        net=Path(netgen.__file__).parent;os.environ['PATH']=os.pathsep.join([str(Path(sys.prefix)/'bin'),str(net),str(net.parent/'ngsolve_openblas'),str(Path(sys.base_prefix)),os.environ['PATH']])
        sys.path.insert(0,str(root/'packages/cubit-mesh-export/src'));from cubit_mesh_export.check import check_consistency
        checked=check_consistency(args.vol,contract=Path(__file__).with_name('native_hcurl_labels.json'),strict_labels=True)
        (out/'vol_check.json').write_text(json.dumps(checked,indent=2)+'\n')
        if not checked['passed']:raise RuntimeError('Checked .vol gate failed')
        mesh=ng.Mesh(str(args.vol));fes=ng.H1(mesh,order=1);u,v=fes.TnT();ng.SetNumThreads(4)
        start=time.perf_counter()
        with ng.TaskManager():
            sl=(LaplaceSL(u.Trace()*ng.ds(bonus_intorder=16),use_fmm=False)*v.Trace()*ng.ds(bonus_intorder=16)).mat
            dl=(LaplaceDL(u.Trace()*ng.ds(bonus_intorder=16),use_fmm=False)*v.Trace()*ng.ds(bonus_intorder=16)).mat
        record['cpu_assembly_seconds']=time.perf_counter()-start
        def dense(a):
            i,j,x=a.COO();return coo_matrix((x,(i,j)),shape=(a.height,a.width)).toarray()
        SL=dense(sl);DL=dense(dl);record['mesh_sha256']=hashlib.sha256(args.vol.read_bytes()).hexdigest()
        receipt=json.loads((args.binary/'build.json').read_text());record['build']=receipt
        if hashlib.sha256((args.binary/'radia_cuda_fmm.mexw64').read_bytes()).hexdigest()!=receipt['binary_sha256']:raise RuntimeError('Binary hash mismatch')
        engine=matlab.engine.start_matlab('-nodesktop -nosplash');engine.addpath(str(Path(__file__).parent),nargout=0)
        engine.fmm_boundary(str(args.vol.resolve()),str(args.binary.resolve()),str(root/'packages/radia-gpu/matlab'),str(out/'native.mat'),nargout=0)
        result=loadmat(out/'native.mat',simplify_cells=True);ids=np.asarray(result['global_dofs'],dtype=int);x=result['x'];expectedSL=SL[np.ix_(ids,ids)]@x;expectedDL=DL[np.ix_(ids,ids)]@x
        metrics=dict(SL_relative=float(np.linalg.norm(result['slResult']-expectedSL)/np.linalg.norm(expectedSL)),DL_relative=float(np.linalg.norm(result['dlResult']-expectedDL)/np.linalg.norm(expectedDL)))
        record.update(metrics=metrics,report=result['report']);print(json.dumps(metrics),flush=True)
        if max(metrics.values())>1e-6:raise RuntimeError('Independent Galerkin boundary operator gate failed')
        record['accepted']=True
    except BaseException as e:record['error']=repr(e);raise
    finally:
        if engine is not None:engine.quit()
        (out/'result.json').write_text(json.dumps(record,indent=2,default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v.item())+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--vol',type=Path,required=True);main(p.parse_args())
