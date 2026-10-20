"""Analytic sphere layer potentials and circular-current Biot-Savart checks."""
from pathlib import Path
import argparse,hashlib,json,socket,os,sys
import numpy as np
from scipy.io import savemat,loadmat

def main(args):
    root=Path(__file__).resolve().parents[2];out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    record=dict(schema='radia.fmm-physical-adapters.v1',accepted=False,host=socket.gethostname(),timing_accepted=False,
        scope='smooth sphere potentials and axial circular current; not singular Galerkin traces or full FEM-BEM solve')
    engine=None
    try:
        mu,wm=np.polynomial.legendre.leggauss(24);phi=(np.arange(48)+.5)*2*np.pi/48
        z=np.repeat(mu,48);theta=np.tile(phi,24);sphere=np.column_stack([np.sqrt(1-z*z)*np.cos(theta),np.sqrt(1-z*z)*np.sin(theta),z]);weights=np.repeat(wm,48)*2*np.pi/48
        rng=np.random.default_rng(1247);directions=rng.normal(size=(64,3));directions/=np.linalg.norm(directions,axis=1)[:,None]
        outside=directions*np.linspace(1.8,3,64)[:,None];inside=directions[:32]*.3;targets=np.vstack([outside,inside])
        count=2048;angle=(np.arange(count)+.5)*2*np.pi/count;coil=np.column_stack([np.cos(angle),np.sin(angle),np.zeros(count)])
        current=np.column_stack([-np.sin(angle),np.cos(angle),np.zeros(count)])*2*np.pi/count
        height=np.linspace(.2,2,64);coilTargets=np.column_stack([np.zeros(64),np.zeros(64),height])
        savemat(out/'input.mat',dict(sphere=sphere,weights=weights[:,None],targets=targets,coil=coil,current=current,coilTargets=coilTargets))
        receipt=json.loads((args.binary/'build.json').read_text());record['build']=receipt
        if hashlib.sha256((args.binary/'radia_cuda_fmm.mexw64').read_bytes()).hexdigest()!=receipt['binary_sha256']:raise RuntimeError('Binary hash mismatch')
        import netgen,matlab.engine
        net=Path(netgen.__file__).parent;os.environ['PATH']=os.pathsep.join([str(Path(sys.prefix)/'bin'),str(net),str(net.parent/'ngsolve_openblas'),str(Path(sys.base_prefix)),os.environ['PATH']])
        engine=matlab.engine.start_matlab('-nodesktop -nosplash');engine.addpath(str(Path(__file__).parent),nargout=0)
        engine.fmm_physics(str(out/'input.mat'),str(args.binary.resolve()),str(root/'packages/radia-simulink/matlab'),str(out/'output.mat'),nargout=0)
        result=loadmat(out/'output.mat',simplify_cells=True);r=np.linalg.norm(targets,axis=1)
        expectedSL=np.where(r>1,1/r,1.);expectedDL=np.where(r>1,0.,-1.)
        expectedB=np.column_stack([np.zeros(64),np.zeros(64),1.25663706212e-6/(2*(1+height**2)**1.5)])
        expectedCloud=np.zeros_like(targets);expectedCloud[:64]=4*np.pi*outside/np.linalg.norm(outside,axis=1)[:,None]**3
        metrics=dict(single_layer_max_absolute=float(np.max(np.abs(result['slPotential']-expectedSL))),
            double_layer_max_absolute=float(np.max(np.abs(result['dlPotential']-expectedDL))),
            B_relative=float(np.linalg.norm(result['B']-expectedB)/np.linalg.norm(expectedB)),
            A_axis_max_absolute=float(np.max(np.abs(result['A']))),
            cloud_relative=float(np.linalg.norm(result['cloudField']-expectedCloud)/np.linalg.norm(expectedCloud)))
        record.update(metrics=metrics,report=result['report'])
        if max(metrics[k] for k in ('single_layer_max_absolute','double_layer_max_absolute','B_relative','cloud_relative'))>1e-5:raise RuntimeError('Independent analytic physical adapter gate failed')
        if metrics['A_axis_max_absolute']>1e-10:raise RuntimeError('Vector potential symmetry gate failed')
        record['accepted']=True;print(json.dumps(metrics,indent=2),flush=True)
    except BaseException as e:record['error']=repr(e);raise
    finally:
        if engine is not None:engine.quit()
        (out/'result.json').write_text(json.dumps(record,indent=2,default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v.item())+'\n',encoding='utf-8')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);main(p.parse_args())
