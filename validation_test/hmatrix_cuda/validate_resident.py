"""Independent small H-matrix algebra and CPU pybind sampled-Gram checks."""
from pathlib import Path
import argparse,ctypes,hashlib,importlib,json,subprocess,sys
import numpy as np
import psutil
from scipy.io import savemat
from scipy.sparse import csc_matrix
ROOT=Path(__file__).resolve().parents[2]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def structs(rows):
    out=np.empty((1,len(rows)),dtype=[(key,'O') for key in rows[0]])
    for i,row in enumerate(rows):
        for key,value in row.items():out[key][0,i]=value
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--artifact',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);a=parser.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    record=dict(schema='radia.hmatrix-resident-numerical.v1',accepted=False,a100_used=False,timing_acceptance=False,scope='Resident symmetric leaf algebra and public sampled-Gram CPU parity; not Matrix Build or HDiv-MMM physics acceptance')
    eng=None
    source_files=(Path(__file__),Path(__file__).with_name('checkResidentHMatrix.m'))
    source_hashes={str(p.relative_to(ROOT)):sha(p) for p in source_files}
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        record['cpu_mean_percent']=psutil.cpu_percent(interval=2)
        if record['cpu_mean_percent']>=30 or psutil.virtual_memory().available<8*2**30:raise RuntimeError('Host budget unavailable')
        receipt=json.loads((a.artifact/'build.json').read_text(encoding='utf-8'));assert receipt['accepted'];assert sha(a.artifact/'radia_cuda_hmatrix.mexw64')==receipt['binary_sha256'];record['build']=receipt
        rng=np.random.default_rng(20261012);cases=[]
        for mode in ('dense','lowrank','mixed'):
            n=9;matrix=np.zeros((n,n));leaves=[]
            for i in range(3):
                for j in range(i,3):
                    rank=0 if i==j or mode=='dense' or (mode=='mixed' and i==0 and j==2) else 2
                    if rank:
                        a1=rng.normal(size=(3,rank))*.03;a2=rng.normal(size=(3,rank))*.03;block=a2@a1.T
                    else:
                        block=rng.normal(size=(3,3))*.01
                        if i==j:block=(block+block.T)*.5+np.eye(3)*.2
                        a1=block.ravel(order='C').reshape(-1,1);a2=np.empty((0,1))
                    matrix[3*i:3*i+3,3*j:3*j+3]=block
                    if i!=j:matrix[3*j:3*j+3,3*i:3*i+3]=block.T
                    leaves.append(dict(row=float(3*i+1),column=float(3*j+1),rows=3.,columns=3.,rank=float(rank),a1=a1,a2=a2))
            permutation=rng.permutation(n);scale=rng.uniform(.5,1.5,size=n);original=np.zeros_like(matrix);original[np.ix_(permutation,permutation)]=matrix;original=scale[:,None]*original*scale[None,:]
            initial=rng.normal(size=(n,1));first=original@initial;third=original@(original@first)
            cases.append(dict(name=mode,config=dict(schema_version=1.,operator_kind='gram',dimension=float(n),permutation=(permutation+1).astype(float).reshape(-1,1),scale=scale[:,None],leaves=structs(leaves)),initial=initial,expected_first=first,expected_third=third))
        # Public CPU pybind route: regularized sampled-Laplace Gram is an explicit
        # numerical scope, never a replacement for integrated physical FE entries.
        import ngsolve as ng
        ng.SetNumThreads(1);native=importlib.import_module('radia._radia_pybind')
        points=rng.uniform(-1,1,size=(8,3));weights=rng.uniform(.1,.2,size=8)
        with ng.TaskManager():
            gram=native._ChargeGramHMatrix.from_sampled_laplace(points.ravel(),weights,.1,1e-12,64,2.,True)
            matrix=np.array([[gram.entry(i,j) for j in range(8)] for i in range(8)])
            initial=rng.normal(size=8);first=np.asarray(gram.matvec_sym(initial));third=np.asarray(gram.matvec_sym(gram.matvec_sym(first)))
        cases.append(dict(name='public-pybind-sampled-Gram',config=dict(schema_version=1.,operator_kind='gram',dimension=8.,permutation=np.arange(1,9,dtype=float)[:,None],scale=np.ones((8,1)),leaves=structs([dict(row=1.,column=1.,rows=8.,columns=8.,rank=0.,a1=matrix.ravel(order='C')[:,None],a2=np.empty((0,1)))])),initial=initial[:,None],expected_first=first[:,None],expected_third=third[:,None]))
        # Preserve the physical charge-map composition, not just the Gram apply.
        gram_case=cases[-1]
        for name,charge in [('demag-dense-map',rng.normal(size=(8,5))*.2),('demag-sparse-map',rng.normal(size=(8,5))*.2)]:
            if 'sparse' in name:charge[abs(charge)<.15]=0.
            config=dict(gram_case['config']);config.update(operator_kind='demag',charge_map=csc_matrix(charge) if 'sparse' in name else charge)
            x=rng.normal(size=5)
            with ng.TaskManager():
                def apply(v):return charge.T@np.asarray(gram.matvec_sym(charge@v))
                y=apply(x);y3=apply(apply(y))
            cases.append(dict(name=name,config=config,initial=x[:,None],expected_first=y[:,None],expected_third=y3[:,None]))
        record['CPU_oracle']=dict(module='radia._radia_pybind',binary_sha256=sha(native.__file__),ngsolve=ng.__version__,note='Installed independent CPU reference; sampled Gram only')
        inputfile=a.out/'inputs.mat';savemat(inputfile,dict(cases=structs(cases)));record['input_sha256']=sha(inputfile)
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');pid=int(eng.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        eng.addpath(str(a.artifact),str(Path(__file__).parent),nargout=0)
        eng.checkResidentHMatrix(str(inputfile),str(a.out/'matlab.json'),nargout=0)
        result=json.loads((a.out/'matlab.json').read_text(encoding='utf-8'));assert result['accepted']
        if source_hashes!={str(p.relative_to(ROOT)):sha(p) for p in source_files}:raise RuntimeError('Validation source changed while running')
        record.update(accepted=True,result=result,matlab=eng.version(),source_sha256=source_hashes,memory_limitations='Exact owned device buffers; no process/global peak or idle timing measured')
    except BaseException as error:record['error']=f'{type(error).__name__}: {error}';raise
    finally:
        if eng is not None:eng.quit()
        (a.out/'result.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
