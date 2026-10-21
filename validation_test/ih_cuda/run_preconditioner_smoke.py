"""Check fixed FP32 preconditioning with FP64 outer arithmetic (algebra only)."""
from pathlib import Path
import argparse,ctypes,hashlib,json,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[2]

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('build_directory',type=Path);args=parser.parse_args();out=args.build_directory.resolve()
    record=dict(schema='radia.ih-fp32-preconditioner.v1',accepted=False,scope='Synthetic coupled operators only; actual PEEC-BEM qualification pending',timing_acceptance=False,transfer_inclusive_time_s=None,peak_memory_bytes=None,a100_used=False)
    eng=None
    try:
        receipt=json.loads((out/'build.json').read_text(encoding='utf-8'));assert receipt['accepted']
        record['build']=receipt
        binary=out/'radia_cuda_ih.mexw64';assert hashlib.sha256(binary.read_bytes()).hexdigest()==receipt['binary_sha256']
        sources=[Path(__file__),Path(__file__).with_name('makeResidentAlgebraConfig.m')]
        record['validation_source_sha256']={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash');pid=int(eng.feature('getpid'))
        subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        eng.addpath(str(out),str(Path(__file__).parent),nargout=0)
        a=np.array([[4+1j,.2-.1j,0],[.1+.2j,5-.5j,.3],[.2,0,6+.3j]])
        x0=np.array([1+.2j,-.3+.4j,.7-.1j]);p=np.array([[1,.2,0],[0,1,.1],[.3,0,1],[1,0,.4],[.2,1,0],[0,.3,1]])
        k=np.array([[2.,-1.,0.],[-1.,2.,-1.],[0.,-1.,2.]])
        projection=np.array([[1.,0.],[.5,.5],[0.,1.]])
        cases=[];baseline={}
        for precision in ('none','fp64','fp32'):
            eng.eval("cfg=makeResidentAlgebraConfig(); cfg.bh_mode='linear'; cfg.preconditioner_inverse=diag(1./diag(cfg.eddy)); cfg.preconditioner_inverse(1,2)=.002+.001i; cfg.preconditioner_inverse(2,3)=-.001i; cfg.preconditioner_inverse(3,1)=.001;",nargout=0)
            eng.workspace['precision_mode']=precision;eng.eval('cfg.preconditioner_precision=precision_mode;',nargout=0)
            for current in (0.,1.,-2.5):
                eng.workspace['drive']=current;eng.eval("h=radia_cuda_ih('create',cfg);",nargout=0)
                try:
                    eng.eval("for step=1:10,radia_cuda_ih('update',h,drive,293.15);end;[temperature,heat,info]=radia_cuda_ih('outputs',h); state=radia_cuda_ih('electromagnetic_state',h);",nargout=0)
                    actual=np.asarray(eng.workspace['state']).ravel();temperature=np.asarray(eng.workspace['temperature']).ravel();heat=np.asarray(eng.workspace['heat']).ravel();info=eng.workspace['info']
                    rhs=a@x0*current;expected=np.linalg.solve(a,rhs);loss=.5*np.sum(np.abs(p@expected).reshape(2,3)**2,axis=1)
                    expected_temperature=np.array([300.,301.,302.])
                    for _ in range(10):expected_temperature+=np.linalg.solve(np.eye(3)+.01*k+.002*np.eye(3),.01*(projection@loss-k@expected_temperature+.2*(293.15-expected_temperature)))
                    error=lambda x,y:float(np.linalg.norm(x-y)/max(np.linalg.norm(y),1.))
                    xe,qe,te=error(actual,expected),error(heat,loss),error(temperature,expected_temperature)
                    residual=float(np.linalg.norm(a@actual-rhs)/max(np.linalg.norm(rhs),1e-30))
                    if not all(np.isfinite(v) for v in (xe,qe,te,residual)) or max(xe,qe,te)>1e-9 or residual>1e-10:raise RuntimeError('FP64 CPU reference or true residual mismatch')
                    if precision=='none':baseline[current]=(actual.copy(),temperature.copy(),heat.copy())
                    ref=baseline[current]
                    cases.append(dict(preconditioner_precision=precision,outer_precision='fp64',operator_precision='fp64',true_residual_precision='fp64',current_A=current,state_cpu_relative=xe,heat_cpu_relative=qe,temperature_cpu_relative=te,state_fp64_baseline_relative=error(actual,ref[0]),temperature_fp64_baseline_relative=error(temperature,ref[1]),heat_fp64_baseline_relative=error(heat,ref[2]),independent_true_residual=residual,native_true_residual=float(info['eddy_true_residual']),eddy_iterations=int(info['eddy_iterations']),thermal_iterations=int(info['thermal_iterations']),owned_device_buffer_bytes=int(info['owned_device_bytes']),memory_limitations='Owned numerical buffers only; excludes CUDA graph/cuBLAS/allocator overhead and is not a sampled peak',steps=int(info['steps'])))
                    eng.eval("radia_cuda_ih('reset',h);[t0,q0,resetInfo]=radia_cuda_ih('outputs',h);assert(isequal(t0,cfg.initial_temperature_K)&&all(q0==0)&&resetInfo.steps==0);",nargout=0)
                finally:eng.eval("radia_cuda_ih('destroy',h);",nargout=0)
        if eng.radia_cuda_ih('live_count')!=0:raise RuntimeError('IH context leak')
        record.update(accepted=True,cases=cases,matlab=eng.version(),contexts_released=True)
    except BaseException as error:record['error']=str(error);raise
    finally:
        if eng is not None:eng.quit()
        (out/'preconditioner-result.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
if __name__=='__main__':main()
