"""Resident nonlinear-transfer cache acceptance with changed SPD values."""
import argparse,json,os,shutil,tempfile
from pathlib import Path
import numpy as np
from scipy.sparse import diags
import ngsolve as ng
from mixed_omega_amg_resident import ResidentAMG
from mixed_omega_amg_cache import load_nonlinear_transfers


def main():
    p=argparse.ArgumentParser(__doc__);p.add_argument('--native-directory',type=Path,required=True);p.add_argument('--gpu-dll-directory',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    handle=os.add_dll_directory(str(a.gpu_dll_directory));record=dict(passed=False);temporary=tempfile.TemporaryDirectory();cache=Path(temporary.name)/'transfers.h5'
    try:
        matrix=diags([-np.ones(299),2.1*np.ones(300),-np.ones(299)],[-1,0,1]).tocsr();candidate=np.ones(300);runs=[]
        for mode in ['build','reuse']:
            current=matrix if mode=='build' else (matrix+diags(np.linspace(.1,1.,300))).tocsr()
            solver=ResidentAMG(current.indptr.astype(np.int32),current.indices.astype(np.int32),current.data,directory=a.native_directory,candidate=candidate,damping=1.5,nonlinear_cache_path=cache,nonlinear_cache_mode=mode)
            try:
                x,info=solver.solve(np.ones(300));res=float(np.linalg.norm(current@x-np.ones(300))/np.sqrt(300));ref=np.linalg.solve(current.toarray(),np.ones(300));error=float(np.linalg.norm(x-ref)/np.linalg.norm(ref));assert res<1e-6 and error<1e-8
                runs.append(dict(mode=mode,residual=res,error=error,info=info,receipt=solver.receipt))
            finally:solver.close()
        assert runs[-1]['receipt']['cache']['updated_operator'] and runs[-1]['receipt']['cache']['galerkin_operators_rebuilt']
        import h5py
        with h5py.File(cache,'r') as f:header=json.loads(f.attrs['header'])
        rejected=[]
        for name,path,b,parameters,layout in [('layout',cache,candidate[:,None],header['parameters'],'0'*64),('candidate',cache,(candidate*2)[:,None],header['parameters'],header['layout']),('settings',cache,candidate[:,None],{},header['layout'])]:
            try:load_nonlinear_transfers(path,b,parameters=parameters,layout=layout)
            except ValueError:rejected.append(name)
            else:raise AssertionError('Invalid nonlinear cache accepted: '+name)
        corrupt=cache.with_name(cache.stem+'-corrupt.h5');shutil.copyfile(cache,corrupt)
        with h5py.File(corrupt,'r+') as f:f['p0/val'][0]*=1.01
        try:load_nonlinear_transfers(corrupt,candidate[:,None],parameters=header['parameters'],layout=header['layout'])
        except ValueError:rejected.append('payload')
        else:raise AssertionError('Corrupt payload accepted')
        record.update(passed=True,runs=runs,expected_rejections=rejected)
    except Exception as error:record['failure']=repr(error);raise
    finally:
        temporary.cleanup();handle.close()
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(record,indent=2)+'\n')

if __name__=='__main__':main()
