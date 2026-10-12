"""Compile both C dispatch contracts on 100; no GPU/physical acceptance claim."""
from pathlib import Path
import argparse,ctypes,hashlib,json,socket,subprocess

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);parser.add_argument('--compression',action='store_true');a=parser.parse_args()
    if socket.gethostname().lower()!='t7920':raise RuntimeError('Native build belongs to 100')
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
    root=Path(__file__).resolve().parents[2];a.out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).with_name('check_aca_regression.c') if a.compression else Path(__file__).with_suffix('.c');header=root/'src/ext/HACApK/cHACApK_calc_entry_ij.h'
    extra=[root/'packages/radia-gpu/native/hmatrix_aca_host.c',root/'src/ext/HACApK/cHACApK_lib.c'] if a.compression else []
    dependencies=extra+[root/'src/ext/HACApK/cHACApK_aca.inc'] if a.compression else []
    result=dict(schema='radia.hacapk-entry-batch-contract.v1',accepted=False,a100_used=False,
                scope='Shared C ACA search regression only; no CUDA/whole solver acceptance' if a.compression else 'C active entry extraction only; no CUDA callback/ACA compression/physics acceptance',cases={},
                source_sha256={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,header,Path(__file__)]+dependencies})
    try:
        for mode in ['scalar','batch']:
            define='/DRAD_HACAPK_BATCH_ENTRIES' if mode=='batch' else ''
            exe=a.out/(mode+'.exe');batch=a.out/(mode+'.cmd')
            inputs=' '.join('"'+str(p)+'"' for p in [source]+extra)
            batch.write_text('@echo off\ncall "C:/BuildTools/VS2022-CUDA/VC/Auxiliary/Build/vcvars64.bat" >nul\n'
                f'cl /nologo /TC /MD /O2 {define} /I"{header.parent}" {inputs} /Fe:"{exe}" /Fo:"{a.out}/"\n',encoding='utf-8')
            build=subprocess.run(['cmd','/c',str(batch)],capture_output=True)
            (a.out/(mode+'-build.log')).write_bytes(build.stdout+build.stderr);build.check_returncode()
            run=subprocess.run([str(exe)],capture_output=True,text=True,check=True)
            result['cases'][mode]=json.loads(run.stdout)
        result['accepted']=True
    finally:(a.out/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
