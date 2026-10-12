"""Small non-timing IH CUDA algebra check; always retains a JSON result."""
from pathlib import Path
import ctypes,hashlib,json,os,subprocess,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OUT=Path('C:/temp/ih-cuda-resident-smoke-20261012')

def main():
    OUT.mkdir(exist_ok=True)
    record=dict(schema='radia.ih-cuda-algebra.v1',accepted=False,
        scope='Resident algebra only; not a coupled IH/Simulink solver acceptance',
        timing_acceptance=False,a100_used=False)
    eng=None
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        import matlab.engine,matlab
        eng=matlab.engine.start_matlab('-nodesktop -nosplash')
        pid=int(eng.feature('getpid'))
        subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        compiler=Path('C:/BuildTools/VS2022-CUDA/VC/Tools/MSVC/14.44.35207/bin/Hostx64/x64')
        os.environ['NVCC_CCBIN']=str(compiler)
        eng.setenv('NVCC_CCBIN',str(compiler),nargout=0)
        eng.setenv('PATH',str(compiler)+';'+eng.getenv('PATH'),nargout=0)
        config=Path(eng.matlabroot())/'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
        src=Path(__file__).with_name('resident_smoke.cu')
        header=ROOT/'packages/radia-gpu/native/ih_resident_kernels.cuh'
        record['source_dirty']=bool(subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True).strip())
        record['source_commit']=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()
        record['source_sha256']={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in (src,header,header.with_name("ih_bicgstab_kernels.cuh"),header.with_name("ih_complex_solver.cuh"),header.with_name("ih_thermal_kernels.cuh"),header.with_name("ih_coupled_state.cuh"),Path(__file__))}
        eng.mexcuda('-v','-R2018a','-f',str(config),'NVCCFLAGS=--fmad=false -gencode=arch=compute_75,code=sm_75',
            '-I'+str(header.parent),'-lcublas','-outdir',str(OUT),str(src),nargout=0)
        eng.addpath(str(OUT),nargout=0)
        rng=np.random.default_rng(20261012)
        a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        p=rng.normal(size=(6,3))+1j*rng.normal(size=(6,3))
        x=rng.normal(size=(3,1))+1j*rng.normal(size=(3,1))
        cases=[]
        for current in (0.,1.,-2.5):
            y,q,info,temperature=eng.resident_smoke(matlab.double(a.tolist(),is_complex=True),
                matlab.double(p.tolist(),is_complex=True),matlab.double(x.tolist(),is_complex=True),current,nargout=4)
            expected=np.linalg.solve(a,a@(current*x));loss=.5*np.sum(np.abs(p@expected).reshape(2,3)**2,axis=1)
            k=np.array([[2.,-1.,0.],[-1.,2.,-1.],[0.,-1.,2.]])
            initial=np.array([300.,301.,302.]);mapping=np.array([[1.,0.],[.5,.5],[0.,1.]])
            load=mapping@loss
            thermal_expected=initial+np.linalg.solve(np.eye(3)+.01*k+.002*np.eye(3),.01*(load-k@initial+.2*(293.15-initial)))
            resident_expected=initial.copy()
            for step in range(10):
                resident_expected+=np.linalg.solve(np.eye(3)+.01*k+.002*np.eye(3),.01*(load-k@resident_expected+.2*(293.15-resident_expected)))
            re=float(np.linalg.norm(np.asarray(info["resident_temperature"]).ravel()-resident_expected)/np.linalg.norm(resident_expected))
            if re>1e-10 or int(info["resident_steps"])!=10 or not all(info[name] for name in ("reset_checked","invalid_update_checked","missing_material_rejected")):
                raise RuntimeError("Persistent coupled state mismatch")
            te=float(np.linalg.norm(np.asarray(temperature).ravel()-thermal_expected)/np.linalg.norm(thermal_expected))
            actual=np.asarray(y);actual_q=np.asarray(q).ravel()
            ye=float(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected),1.))
            qe=float(np.linalg.norm(actual_q-loss)/max(np.linalg.norm(loss),1.))
            if not np.isfinite(actual).all() or not np.isfinite(actual_q).all() or max(ye,qe,te)>1e-10:raise RuntimeError('CPU/CUDA algebra mismatch')
            cases.append(dict(current_A=current,complex_operator_relative=ye,quadratic_loss_relative=qe,resident_repeats=int(info["resident_repeats"]),bicgstab_identity_first_correction_checked=True,thermal_relative=te,persistent_temperature_relative=re,persistent_steps=10,reset_checked=True,invalid_update_checked=True,missing_material_rejected=True,eddy_true_residual=float(info["eddy_true_residual"]),thermal_true_residual=float(info["thermal_true_residual"])))
        record.update(accepted=True,cases=cases,matlab=eng.version(),binary_sha256=hashlib.sha256((OUT/'resident_smoke.mexw64').read_bytes()).hexdigest())
    except BaseException as error:
        record['error']=str(error)
        raise
    finally:
        if eng is not None:eng.quit()
        (OUT/'result.json').write_text(json.dumps(record,indent=2)+'\n')
if __name__=='__main__':main()
