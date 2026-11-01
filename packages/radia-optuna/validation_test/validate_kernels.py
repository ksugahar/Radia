"""Formal stored Optuna5 density/selection comparison; never regenerates truth."""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import sys
import ctypes

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def current_context():
    # Independent non-creating driver observation, not the kernel's counters.
    driver=ctypes.WinDLL('nvcuda.dll',winmode=0x800)
    get=driver.cuCtxGetCurrent
    get.argtypes=[ctypes.POINTER(ctypes.c_void_p)];get.restype=ctypes.c_int
    context=ctypes.c_void_p();status=get(ctypes.byref(context))
    if status not in (0,3):raise RuntimeError(f'Current-context probe failed: {status}')
    return {'driver_status':status,'current_context':context.value or 0}
def descriptor(d):
    categorical=d['categorical']
    return dict(weights=d['weights'],categorical=categorical,
                probabilities=[x for row in d.get('probabilities',[]) for x in row],
                choices=len(d['probabilities'][0]) if categorical else 0,
                mu=d.get('mu',[]),sigma=d.get('sigma',[]),
                low=d.get('internal_low',0),high=d.get('internal_high',1),
                step=d.get('step') if d.get('step') is not None else float('nan'),
                log_scale=d.get('log',False))
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--backend',choices=['cpu','gpu'],required=True)
    p.add_argument('--build-receipt',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    a=p.parse_args()
    before=current_context() if a.backend=='cpu' else None
    from radia_optuna.kernels import load
    k=load(a.backend)
    build=json.loads(a.build_receipt.read_text())
    binary=Path(k.__file__).resolve()
    assert binary.parent==(ROOT/'src/radia_optuna').resolve(),'Wrong native import origin'
    assert sha(binary)==build['python_extension_sha256s'].get(binary.name),'Native binding hash mismatch'
    oracle=ROOT/'oracles/optuna50_cuda_oracle.json'
    frozen=sha(oracle);data=json.loads(oracle.read_text());assert data['optuna_version']=='5.0.0'
    cases=data['cases'];records=[]
    for entry in cases:
        samples=[x for row in entry['samples'] for x in row]
        lower=[descriptor(d) for d in entry['below']]
        upper=[descriptor(d) for d in entry['above']]
        below=k.log_pdf(a.backend,samples,lower)
        above=k.log_pdf(a.backend,samples,upper)
        assert len(below)==len(entry['below_pdf']) and len(above)==len(entry['above_pdf']),entry['name']
        error=max(abs(x-y) for got,ref in [(below,entry['below_pdf']),(above,entry['above_pdf'])] for x,y in zip(got,ref))
        # Existing formal MATLAB validate_optuna_cuda density limit, unchanged.
        assert error<=1e-10,(entry['name'],error,1e-10)
        selected=k.best(a.backend,samples,lower,upper)
        assert selected+1==entry['best_index'],(entry['name'],selected+1,entry['best_index'])
        assert k.best(a.backend,samples,lower,lower)==0,'First-maximum tie order'
        records.append({'name':entry['name'],'density_max_abs_error':error,
            'below_bits':[struct.pack('>d',x).hex() for x in below],
            'above_bits':[struct.pack('>d',x).hex() for x in above],'best_index_one_based':selected+1})
    assert sha(oracle)==frozen,'Stored oracle changed'
    stats=k.statistics()
    if a.backend=='cpu':
        assert stats['cuda_device']==-1 and stats['density_launches']==0 and stats['host_to_device_bytes']==0
        assert before['current_context']==0 and current_context()['current_context']==0,'CPU created a context'
    else:assert stats['density_launches']>0 and stats['selection_launches']>0
    lifecycle=None
    if a.backend=='gpu':
        from runtime_lifecycle import validate
        lifecycle=validate(k)
    k.shutdown() # release owned allocations before interpreter/DLL teardown
    source={f.relative_to(ROOT).as_posix():sha(f) for f in sorted((ROOT/'native').glob('*')) if f.is_file()}
    assert source==build['source_hashes'],'Build receipt source mismatch'
    report={'schema':'radia.optuna.python-kernel-verification.v1','ok':True,'backend':a.backend,
            'source_sha':build['source_sha'],'recipe_hashes':build['recipe_hashes'],
            'lifecycle':lifecycle,
            'CPU_driver_observation':{'before':before,'after':current_context()} if a.backend=='cpu' else None,
            'source_hashes':source,'archive_sha256':build['archive_sha256'],
            'oracle_sha256':frozen,'cases':records,'statistics':stats,
            'scope':'Stored density/selection cases only; model-build, wrapper bitwise parity and MATLAB RNG remain separate gates'}
    a.report.parent.mkdir(parents=True,exist_ok=True)
    a.report.write_text(json.dumps(report,indent=2),encoding='utf-8',newline='\n')
if __name__=='__main__':main()
