"""Isolated, low-priority persistent IH mexcuda build on 100."""
from pathlib import Path
import argparse,ctypes,hashlib,json,re,subprocess,sys,tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ih_payload import REQUIRED_SOURCES, NATIVE_SOURCES

def build_mode(transaction_test):
    return dict(transaction_test_build=transaction_test,distribution_eligible=not transaction_test,
                output_directory_name='build-ih-transaction-test' if transaction_test else 'build-ih-cuda')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--require-clean-source',action='store_true');parser.add_argument('--transaction-test',action='store_true');parser.add_argument('--describe-build-mode',action='store_true');args=parser.parse_args()
    mode=build_mode(args.transaction_test)
    if args.describe_build_mode:print(json.dumps(mode));return
    root=Path(__file__).resolve().parents[3];package=root/'packages/radia-gpu';base=root/mode['output_directory_name'];base.mkdir(exist_ok=True);out=Path(tempfile.mkdtemp(prefix='job-',dir=base));print('IH CUDA build directory: '+str(out),flush=True)
    record=dict(schema='radia.ih-cuda-build.v1',accepted=False,scope='Research complete flat-P1 reaction and element-temperature amplitude B-H coupling',features=['element-temperature-bh','complete-flat-p1-reaction','accepted-element-impedance','transactional-temperature-range'],a100_used=False)
    record.update(mode)
    eng=None
    try:
        record['source_commit']=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
        record['source_dirty']=bool(subprocess.check_output(['git','-C',str(root),'status','--porcelain'],text=True).strip())
        if args.require_clean_source and record['source_dirty']:raise RuntimeError('Clean source required')
        sources=[root/name for name in sorted(REQUIRED_SOURCES)]
        gateway=package/'native/radia_cuda_ih.cu'
        native_text='\n'.join((package/'native'/name).read_text(encoding='utf-8') for name in NATIVE_SOURCES)
        forbidden=r'\b(TaskManager|RegionTaskManager|Py_Initialize|PyImport_ImportModule|pybind11)\b'
        if re.search(forbidden,native_text):raise RuntimeError('IH native-only backend contract violated')
        record['backend_contract_verification']=dict(method='declared build graph and native source token check',native_sources_checked=len(NATIVE_SOURCES),host_define='RAD_MATLAB_HOST',link_libraries=['cublas'],limitations='Static build contract; not numerical or whole-wheel qualification')
        hashes=lambda:{p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
        record['source_sha256']=hashes()
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash');pid=int(eng.feature('getpid'))
        subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        config=Path(eng.matlabroot())/'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
        if not config.is_file():raise RuntimeError('Supported MATLAB CUDA MSVC2022 configuration missing')
        flags='NVCCFLAGS=--std=c++17 --fmad=false -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_80,code=compute_80'
        definitions=['-DRADIA_IH_TRANSACTION_TEST'] if args.transaction_test else []
        eng.mexcuda('-v','-R2018a','-f',str(config),flags,'-I'+str(package/'native'),'-I'+str(root/'src/core'),'-DRAD_MATLAB_HOST',*definitions,'-lcublas','-outdir',str(out),str(gateway),str(root/'src/core/rad_bem_galerkin.cpp'),str(root/'src/core/rad_biot_savart_filaments.cpp'),nargout=0)
        if hashes()!=record['source_sha256'] or subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()!=record['source_commit']:raise RuntimeError('Source changed during build')
        record.update(accepted=True,matlab=eng.version(),binary_sha256=hashlib.sha256((out/'radia_cuda_ih.mexw64').read_bytes()).hexdigest(),cuda_targets=['sm_75','sm_80','compute_80 PTX'],cpu_fallback=False,python_runtime_required=False,taskmanager_regions=False)
    except BaseException as error:record['error']=str(error);raise
    finally:
        if eng is not None:eng.quit()
        (out/'build.json').write_text(json.dumps(record,indent=2)+'\n')
if __name__=='__main__':main()
