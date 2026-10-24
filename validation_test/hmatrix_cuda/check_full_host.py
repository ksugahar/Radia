"""Build full existing C HACApK with supported isolated MKL, no TaskManager."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'packages/radia-simulink/tools'))
from vs2022_toolchain import discover_vs, compiler_bin
SELECTOR = Path(__file__).resolve().parents[2] / "packages/radia-simulink/tools/vs2022_toolchain.py"
import argparse, ctypes, hashlib, json, os, socket, subprocess

def main():
    p=argparse.ArgumentParser();p.add_argument('--mkl-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--deferred-dense',action='store_true',help='Opt-in C ownership/leaf-order contract only; no GPU implementation');p.add_argument('--frame-admission',action='store_true',help='Pre-fill rejection and retry C ABI contract, no CUDA');a=p.parse_args()
    if a.frame_admission and a.deferred_dense:raise ValueError('Select one optional contract')
    if socket.gethostname().lower()!='t7920':raise RuntimeError('Build belongs to 100')
    root=Path(__file__).resolve().parents[2];a.out.mkdir(parents=True,exist_ok=False)
    sources=[Path(__file__).with_suffix('.c'),root/'packages/radia-simulink/native/hmatrix_serial_host.c']+[root/'src/ext/HACApK'/s for s in ['cHACApK_base.c','cHACApK_cpp_impl.c','cHACApK_lib.c','cHACApK_pca_cluster.c']]
    tracked=sources+[SELECTOR,Path(__file__),root/'src/ext/HACApK/cHACApK_aca.inc']+list((root/'src/ext/HACApK').glob('*.h'))+[root/'src/core/rad_hacapk_parallel.h']
    result=dict(schema='radia.hacapk-full-host.v1',accepted=False,scope='Synthetic full native cluster/leaf/ACA infrastructure; not CUDA or HDiv physical acceptance',source_commit=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),source_sha256={s.relative_to(root).as_posix():hashlib.sha256(s.read_bytes()).hexdigest() for s in tracked},cases={},timing_measured=False)
    vs=discover_vs()
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        if not (a.mkl_root/'include/mkl.h').is_file() or not (a.mkl_root/'lib/mkl_rt.lib').is_file():raise RuntimeError('Supported mkl-devel headers and import library required')
        result['mkl_import_library_sha256']=hashlib.sha256((a.mkl_root/'lib/mkl_rt.lib').read_bytes()).hexdigest()
        result['mkl_header_sha256']=hashlib.sha256((a.mkl_root/'include/mkl.h').read_bytes()).hexdigest()
        env=dict(os.environ,MKL_NUM_THREADS='1',MKL_THREADING_LAYER='SEQUENTIAL',OMP_NUM_THREADS='1')
        env['PATH']=str(a.mkl_root/'bin')+';C:/Program Files/Python312/Library/bin;'+env['PATH']
        for mode in (['frame'] if a.frame_admission else ['deferred'] if a.deferred_dense else ['scalar','batch']):
            exe=a.out/(mode+'.exe');batch=a.out/(mode+'.cmd')
            inputs=' '.join('"'+str(s)+'"' for s in sources)
            define='/DRAD_HACAPK_BATCH_ENTRIES /DRAD_HACAPK_FRAME_ADMISSION' if mode=='frame' else '/DRAD_HACAPK_BATCH_ENTRIES /DRAD_HACAPK_DENSE_LEAF_SINK' if mode=='deferred' else '/DRAD_HACAPK_BATCH_ENTRIES' if mode=='batch' else ''
            batch.write_text('@echo off\ncall "'+str(vs/'VC/Auxiliary/Build/vcvars64.bat')+'" >nul\n'+f'cl /nologo /TC /MD /O2 /Gy /DRADIA_BLAS_MKL=1 {define} /I"{root/"src/ext/HACApK"}" /I"{root/"src/core"}" /I"{a.mkl_root/"include"}" {inputs} /Fo:"{a.out}/" /Fe:"{exe}" /link /OPT:REF /LIBPATH:"{a.mkl_root/"lib"}" mkl_rt.lib\n',encoding='utf-8')
            built=subprocess.run(['cmd','/c',str(batch)],cwd=a.out,env=env,capture_output=True)
            (a.out/(mode+'-build.log')).write_bytes(built.stdout+built.stderr);built.check_returncode()
            run=subprocess.run([str(exe)],cwd=a.out,env=env,capture_output=True)
            (a.out/(mode+'-run.log')).write_bytes(run.stdout+run.stderr);run.check_returncode()
            result['cases'][mode]=json.loads(run.stdout.decode().splitlines()[-1])
            result['cases'][mode]['binary_sha256']=hashlib.sha256(exe.read_bytes()).hexdigest()
        if result['source_sha256']!={s.relative_to(root).as_posix():hashlib.sha256(s.read_bytes()).hexdigest() for s in tracked}:raise RuntimeError('Source changed during validation')
        result['accepted']=True
    except BaseException as error:
        result['error']=f'{type(error).__name__}: {error}';raise
    finally:(a.out/'result.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
