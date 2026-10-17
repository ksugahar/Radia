"""Element-specific ESIM material plumbing; no independent physics certification."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'packages/radia-gpu/tools'))
from vs2022_toolchain import discover_vs, compiler_bin
SELECTOR = Path(__file__).resolve().parents[2] / "packages/radia-gpu/tools/vs2022_toolchain.py"
import argparse,ctypes,hashlib,json,os,socket,subprocess

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=False);root=Path(__file__).resolve().parents[2]
    r=dict(accepted=False,scope='Three generic constitutive block cells; no complete physical reaction assembly',research_preview=True,
           full_peec_bem_connected=False,simulink_accepted=False,release_accepted=False,
           independent_physics_validation=False,timing_accepted=False)
    try:
        if socket.gethostname().lower()!='t7920':raise RuntimeError('CUDA builds belong to 100')
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        os.environ.update(OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
        header=root/'packages/radia-gpu/native';src=Path(__file__).with_name('esim_coupling_smoke.cu')
        sources=[src,SELECTOR,Path(__file__)]+[header/n for n in ['ih_temperature_material.cuh','ih_esim_cell.cuh','ih_thermal_kernels.cuh','ih_resident_kernels.cuh','ih_esim_coupling.cuh','ih_coupled_state.cuh','ih_complex_solver.cuh','ih_bicgstab_kernels.cuh','ih_preconditioner.cuh']]
        hashes={str(s.relative_to(root)):sha(s) for s in sources};r['source_sha256']=hashes
        r['source_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
        r['source_dirty']=bool(subprocess.check_output(['git','status','--porcelain'],cwd=root,text=True).strip())
        compiler=str(compiler_bin(discover_vs()))
        quote=lambda path: str(path).replace('\\','/').replace("'","''")
        script=a.out/'check.m'
        script.write_text(f"setenv('NVCC_CCBIN','{compiler}');setenv('PATH',['{compiler};' getenv('PATH')]);\n"
          "config=fullfile(matlabroot,'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml');\n"
          f"mexcuda('-v','-R2018a','-f',config,'NVCCFLAGS=--fmad=false -gencode=arch=compute_75,code=sm_75 -gencode=arch=compute_80,code=sm_80','-I{quote(header)}','-lcublas','-outdir','{quote(a.out)}','{quote(src)}');\n"
          f"addpath('{quote(a.out)}');\n"
          "rows=esim_coupling_smoke([300;400;500]);assert(rows(1)<=1e-8 && rows(2)<=1e-8 && rows(3)>=2 && rows(5)<=1e-10 && rows(6)<=1e-10 && rows(7)==5);rejected=false;try,unused=esim_coupling_smoke([299;299;500]);catch e,rejected=contains(e.message,'outside declared table');end;assert(rejected);\n"
          f"f=fopen('{quote(a.out/'native.json')}','w');assert(f>=0);closer=onCleanup(@()fclose(f));fprintf(f,'%s',jsonencode(struct('rows',rows,'matlab',version)));clear closer;\n",encoding='utf-8')
        with (a.out/'matlab.log').open('wb') as log:
            subprocess.run(['C:/Program Files/MATLAB/R2026a/bin/matlab.exe','-singleCompThread','-batch',f"run('{quote(script)}')"],stdout=log,stderr=subprocess.STDOUT,check=True,creationflags=subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        r.update(json.loads((a.out/'native.json').read_text()))
        r.update(binary_sha256=sha(a.out/'esim_coupling_smoke.mexw64'),execution_route='Owned MATLAB batch',
                 repeated_calls=5,out_of_table_rejected=True,accepted=True)
        if hashes!={str(s.relative_to(root)):sha(s) for s in sources}:raise RuntimeError('Source changed during test')
    except BaseException as e:r['accepted']=False;r['error']=str(e);raise
    finally:
        (a.out/'result.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
