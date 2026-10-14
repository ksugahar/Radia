"""Idle-100 operator residency diagnostic; never formal development/build host performance acceptance."""
from pathlib import Path
import argparse,hashlib,importlib.util,json,statistics,subprocess
import psutil
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--artifact',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    record=dict(schema='radia.hmatrix-residency-diagnostic.v1',accepted=False,timing_acceptance=False,a100_used=False,classification='Idle-100 diagnostic, not gpu evaluation worker performance acceptance',scope='Preassembled numeric leaf setup and three repeated operator applications, not FE Matrix Build or HDiv-MMM solve')
    files=[Path(__file__),Path(__file__).with_name('profileResidentHMatrix.m'),ROOT/'packages/radia-gpu/matlab/+radia/+gpu/HMatrixOperator.m'];hashes={str(p.relative_to(ROOT)):sha(p) for p in files};eng=None
    try:
        spec=importlib.util.spec_from_file_location('guard',ROOT/'validation_test/ams_cuda/graph_validation.py');guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
        guard.assert_idle_host();record['cpu_mean_percent']=psutil.cpu_percent(interval=2)
        if record['cpu_mean_percent']>=10 or psutil.virtual_memory().available<8*2**30:raise RuntimeError('Idle CPU/memory budget unavailable')
        build=json.loads((a.artifact/'build.json').read_text(encoding='utf-8'));assert build['accepted'] and sha(a.artifact/'radia_cuda_hmatrix.mexw64')==build['binary_sha256'];record['build']=build;record['input_sha256']=sha(a.inputs)
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');eng.addpath(str(a.artifact),str(Path(__file__).parent),str(ROOT/'packages/radia-gpu/matlab'),nargout=0)
        guard.assert_idle_host();eng.profileResidentHMatrix(str(a.inputs),str(a.out/'matlab.json'),nargout=0);guard.assert_idle_host()
        result=json.loads((a.out/'matlab.json').read_text(encoding='utf-8'));assert result['accepted']
        if hashes!={str(p.relative_to(ROOT)):sha(p) for p in files}:raise RuntimeError('Timing source changed during run')
        record.update(accepted=True,result=result,matlab=eng.version(),source_sha256=hashes,median_whole_seconds={mode:statistics.median(row['whole_seconds'] for row in result['samples'] if row['mode']==mode) for mode in ('rebuild','resident')})
    except BaseException as error:record['error']=f'{type(error).__name__}: {error}';raise
    finally:
        if eng is not None:eng.quit()
        (a.out/'result.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
