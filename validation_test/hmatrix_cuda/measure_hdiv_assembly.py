"""Guarded actual RT0 whole-solve ABBA diagnostic on idle 100, not development/build host acceptance."""
from pathlib import Path
import argparse,ctypes,hashlib,importlib.util,json,os,statistics,subprocess,sys

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser()
    for key in ['before','after','projection','mesh','reference','out']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[2];eng=None
    result=dict(accepted=False,classification='Idle-100 diagnostic, not formal gpu evaluation worker acceptance',scope='Eight affine tetra RT0 FE setup through CUDA HACApK/material solve/result recovery/release; not general HDiv',cold_limitations='First payload per module in an initialized owned MATLAB process, not cold process',memory_limitations='Exact owned buffers plus available-memory samples; excludes cuBLAS/runtime and host/global peaks')
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        spec=importlib.util.spec_from_file_location('guard',root/'validation_test/ams_cuda/graph_validation.py');guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
        guard.assert_idle_host()
        result['hashes']={k:sha(v) for k,v in [('mesh',a.mesh),('cpu_reference',a.reference),('projection',a.projection/'radia_cuda_projection.mexw64'),('wrapper',root/'packages/radia-gpu/matlab/+radia/+gpu/solveHDivRT0.m'),('matlab_runner',Path(__file__).with_suffix('.m').with_name('compareHDivAssembly.m')),('python_runner',Path(__file__))]}
        for key in ['before','after']:
            artifact=getattr(a,key);receipt=json.loads((artifact/'build.json').read_text());assert receipt['accepted'] and sha(artifact/'radia_cuda_hmatrix.mexw64')==receipt['binary_sha256'];result[key+'_build']=receipt
        prefix=Path(sys.base_prefix);os.environ['PATH']=';'.join(map(str,[prefix,prefix/'bin',prefix/'Library/bin',prefix/'Lib/site-packages/netgen',prefix/'Lib/site-packages/netgen_mesher.libs',prefix/'Lib/site-packages/ngsolve.libs',prefix/'Lib/site-packages/ngsolve_openblas']))+';'+os.environ['PATH']
        os.environ['MKL_NUM_THREADS']='1';os.environ['MKL_THREADING_LAYER']='SEQUENTIAL'
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');eng.addpath(str(a.projection),str(Path(__file__).parent),str(root/'packages/radia-gpu/matlab'),nargout=0)
        guard.assert_idle_host();eng.compareHDivAssembly(str(a.mesh),str(a.reference),str(a.before),str(a.after),str(a.out/'matlab.json'),nargout=0);guard.assert_idle_host()
        r=json.loads((a.out/'matlab.json').read_text());assert r['accepted'];result['result']=r;result['matlab']=eng.version()
        result['median_whole_seconds']={str(k):statistics.median(s['whole_seconds'] for s in r['samples'] if s['mode']==k) for k in [1,2]};result['accepted']=True
    except BaseException as e:result['error']=f'{type(e).__name__}: {e}';raise
    finally:
        if eng is not None:eng.quit()
        (a.out/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
