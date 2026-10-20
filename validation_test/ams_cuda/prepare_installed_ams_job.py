"""Stage a private installed-wheel acceptance job; never alter a live runtime."""
from pathlib import Path
import argparse,hashlib,json,shutil

def main(a):
    root=Path(__file__).resolve().parents[2];out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    for n in ('ams_tuning.py','ams_tuning_sample.m','ams_loaded_module.ps1','nonlinear_resident.py','measure_native_cpu_gpu.py'):
        shutil.copy2(Path(__file__).parent/n,out/n)
    shutil.copy2(root/'packages/radia-simulink/verify_wheel.py',out/'verify_wheel.py')
    shutil.copy2(a.wheel,out/a.wheel.name)
    (out/'baseline').mkdir();shutil.copy2(a.baseline_binary,out/'baseline/radia_cuda_ams.mexw64')
    (out/'cpu').mkdir();shutil.copy2(a.cpu_reference,out/'cpu/_radia_ams_cpu_reference.pyd')
    shutil.copy2(a.vol,out/'model.vol');shutil.copy2(a.vol_check,out/'vol_check.json')
    check=json.loads(a.vol_check.read_text(encoding='utf-8-sig'))
    if not check.get('passed'):raise RuntimeError('Required versioned check-vol receipt did not pass')
    shutil.copytree(a.runtime,out/'runtime')
    owner=dict(owner='ams-wheel-acceptance',job=a.job_id,cleanup_trigger='Recover every input/result/log and verify hashes before guarded job-root removal',
        vol_sha256=hashlib.sha256(a.vol.read_bytes()).hexdigest(),vol_check_sha256=hashlib.sha256(a.vol_check.read_bytes()).hexdigest(),
        input_provenance='Independently authored public NGSolve unit cube; source file and checker receipt retained in private job')
    (out/'owner.json').write_text(json.dumps(owner,indent=2)+'\n',encoding='utf-8')
    manifest={p.relative_to(out).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*') if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(staged=True,files=len(manifest),wheel_sha256=manifest[a.wheel.name],vol_sha256=owner['vol_sha256'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('wheel','baseline-binary','cpu-reference','vol','vol-check','runtime','out'):p.add_argument('--'+key,required=True,type=Path)
    p.add_argument('--job-id',required=True);main(p.parse_args())
