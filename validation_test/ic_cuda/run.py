"""Independent NumPy IC oracle and native lifecycle acceptance, owned Engine."""
from pathlib import Path
import argparse,hashlib,json,subprocess

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def check(engine,binary_dir,matlab_dir,out):
    import numpy as np
    from scipy.io import savemat,loadmat
    out=Path(out);out.mkdir(parents=True,exist_ok=False);binary=Path(binary_dir)/'radia_cuda_ic.mexw64'
    A=np.array([[4.,1.,.5],[1.,3.,.2],[.5,.2,2.]]);b=np.array([1.,2.,3.]);fixture=out/'python-reference.mat'
    # Independent square-root IC(0) oracle, distinct from CUDA weighted LDLT.
    sparse_A=5*np.eye(25)
    for i in range(25):
        for j in (i+1 if i%5<4 else 25,i+5):
            if j<25:sparse_A[i,j]=sparse_A[j,i]=-1
    L=np.zeros_like(sparse_A)
    for i in range(25):
        for j in range(i):
            if sparse_A[i,j]!=0:L[i,j]=(sparse_A[i,j]-np.dot(L[i,:j],L[j,:j]))/L[j,j]
        L[i,i]=np.sqrt(1.05*sparse_A[i,i]-np.dot(L[i,:i],L[i,:i]))
    sparse_b=np.linspace(.2,1.7,25)
    sparse_pre=np.linalg.solve(L.T,np.linalg.solve(L,sparse_b))
    savemat(fixture,dict(A=A,b=b[:,None],expected=np.linalg.solve(A,b)[:,None],expected_pre=np.linalg.solve(A+.05*np.diag(np.diag(A)),b)[:,None],sparse_A=sparse_A,sparse_b=sparse_b[:,None],sparse_expected=np.linalg.solve(sparse_A,sparse_b)[:,None],sparse_expected_pre=sparse_pre[:,None]))
    engine.addpath(str(Path(__file__).parent),nargout=0)
    engine.validate_ic(str(fixture),str(Path(binary_dir).resolve()),str(Path(matlab_dir).resolve()),str(out/'numerical.mat'),nargout=0)
    native=loadmat(out/'numerical.mat',simplify_cells=True)['result'];assert native['passed']
    module=out/'loaded-module.json'
    subprocess.run(['powershell','-NoProfile','-ExecutionPolicy','Bypass','-File',str(Path(__file__).with_name('loaded_module.ps1')),'-ProcessId',str(int(engine.feature('getpid'))),'-ModuleName','radia_cuda_ic.mexw64','-OutputPath',str(module)],check=True,capture_output=True)
    loaded=json.loads(module.read_text(encoding='utf-8-sig'))
    if Path(loaded['path']).resolve()!=binary.resolve() or loaded['sha256']!=sha(binary):raise RuntimeError('Actual loaded IC binary mismatch')
    libraries={}
    for name in ('cusparse64_12.dll','cublas64_12.dll'):
        path=out/(name+'.json')
        subprocess.run(['powershell','-NoProfile','-ExecutionPolicy','Bypass','-File',str(Path(__file__).with_name('loaded_module.ps1')),'-ProcessId',str(int(engine.feature('getpid'))),'-ModuleName',name,'-OutputPath',str(path)],check=True,capture_output=True)
        actual=json.loads(path.read_text(encoding='utf-8-sig'))
        if not Path(actual['path']).resolve().is_relative_to(Path(engine.matlabroot()).resolve()):raise RuntimeError('CUDA library resolved outside the selected MATLAB runtime')
        libraries[name]=actual['sha256']
    result=dict(schema='radia.gpu.ic-accuracy.v1',passed=True,binary_sha256=sha(binary),fixture_sha256=sha(fixture),numpy=np.__version__,independent_cpu='NumPy original-A solve, shifted full-pattern solve and independent square-root sparse IC(0) oracle',apply_relative=float(native['apply_relative']),sparse_apply_relative=float(native['sparse_apply_relative']),solution_relative_max=float(np.max(np.asarray(native['solve_errors'])[:,0])),true_relative_residual_max=float(np.max(np.asarray(native['solve_errors'])[:,1])),rejection_count=int(native['rejection_count']),live_count_restored=bool(native['live_count_restored']),frontend_checked=True,device=native['device'])
    result['lifecycle']=dict(engine.validate_ic_lifecycle(nargout=1))
    result['runtime_library_sha256']=libraries
    (out/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return result

if __name__=='__main__':
    import matlab.engine
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary-dir',type=Path,required=True);p.add_argument('--matlab-dir',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();engine=None
    try:
        engine=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');engine.restoredefaultpath(nargout=0);print(json.dumps(check(engine,a.binary_dir,a.matlab_dir,a.out),indent=2))
    finally:
        if engine is not None:engine.quit()
