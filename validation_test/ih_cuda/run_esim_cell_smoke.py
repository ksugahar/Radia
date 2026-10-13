"""Differential ESIM finite-cell physics check, not full PEEC-BEM acceptance."""
from pathlib import Path
import ctypes,hashlib,importlib.util,json,subprocess,sys,types
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OUT=Path('C:/temp/ih-esim-cell-cuda-20261012')

def main():
    OUT.mkdir(exist_ok=True)
    record=dict(schema='radia.ih-esim-cell-cuda.v1',accepted=False,scope='Finite slab/cylinder cells only; full PEEC-BEM, pybind11 physical and Simulink acceptance pending',timing_acceptance=False,a100_used=False)
    eng=None
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        # Isolate the unmodified pure-Python finite-difference oracle from the
        # installed Radia native runtime. This is explicitly not pybind11 parity.
        package=types.ModuleType('radia');package.__path__=[str(ROOT/'src/radia')];sys.modules['radia']=package
        from radia.esim_cell_problem import ESIMFiniteSlabSolver
        import matlab,matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash');pid=int(eng.feature('getpid'))
        subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        compiler='C:/BuildTools/VS2022-CUDA/VC/Tools/MSVC/14.44.35207/bin/Hostx64/x64'
        eng.setenv('NVCC_CCBIN',compiler,nargout=0);eng.setenv('PATH',compiler+';'+eng.getenv('PATH'),nargout=0)
        config=Path(eng.matlabroot())/'toolbox/parallel/gpu/extern/src/mex/win64/nvcc_msvcpp2022.xml'
        src=Path(__file__).with_name('esim_cell_smoke.cu');header=ROOT/'packages/radia-gpu/native/ih_esim_cell.cuh'
        sources=[src,header,Path(__file__),ROOT/'src/radia/esim_cell_problem.py',ROOT/'src/radia/bh_law.py']
        record.update(source_commit=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),source_dirty=bool(subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True).strip()),source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
        eng.mexcuda('-v','-R2018a','-f',str(config),'NVCCFLAGS=--fmad=false -gencode=arch=compute_75,code=sm_75','-I'+str(header.parent),'-outdir',str(OUT),str(src),nargout=0)
        eng.addpath(str(OUT),nargout=0)
        fields=np.array([1e-3,1.,100.,10000.]);cases=[]
        bh=np.array([[0.,0.],[10.,.02],[100.,.18],[1000.,1.],[10000.,1.5],[100000.,1.65]])
        for geometry in ('slab','cylinder'):
            for material in ('linear','complex','nonlinear'):
                options=dict(half_thickness=.001,sigma=5e6,frequency=7000.,n_nodes=100,geometry=geometry)
                if material=='linear':options['mu_r']=100.
                elif material=='complex':options['complex_mu']=(100.,5.)
                else:options['bh_curve']=bh
                solver=ESIMFiniteSlabSolver(**options)
                knots=[];coeff=[];saturation=0.
                if solver.bh_interp is not None:
                    knots=solver.bh_interp.H_data.tolist();coeff=solver.bh_interp._B_interp.c.T.tolist();saturation=solver.bh_interp._saturation_magnetization
                params=[100.,50.,float(geometry=='cylinder'),.001,solver.rho,solver.omega,1e-6,.5,float(np.real(solver.mu_initial)),float(np.imag(solver.mu_initial)),saturation]
                z,info,h=eng.esim_cell_smoke(matlab.double(fields.reshape(-1,1).tolist()),matlab.double([params]),matlab.double([knots]) if knots else matlab.double([]),matlab.double(coeff),nargout=3)
                z=np.asarray(z).ravel();info=np.asarray(info);h=np.asarray(h)
                for i,field in enumerate(fields):
                    cpu=solver.solve(float(field));assert cpu['converged']
                    ze=float(abs(z[i]-cpu['Z'])/abs(cpu['Z']));he=float(np.linalg.norm(h[:,i]-cpu['H_solution'])/np.linalg.norm(cpu['H_solution']))
                    power=np.array([cpu['P_prime'],cpu['Q_prime'],cpu['P_magnetic']]);pe=float(np.linalg.norm(info[i,:3]-power)/max(np.linalg.norm(power),1e-30))
                    if max(ze,he,pe)>1e-9 or info[i,4]>1e-6 or info[i,5]!=cpu['iterations']:raise RuntimeError(f'ESIM CPU/CUDA mismatch: {geometry}/{material}/{field}: {ze},{he},{pe}')
                    cases.append(dict(geometry=geometry,material=material,surface_field_A_per_m=float(field),impedance_relative=ze,field_relative=he,power_relative=pe,true_residual=float(info[i,4]),iterations=int(info[i,5])))
        # A deliberately insufficient nonlinear iteration budget must fail loudly.
        params[1]=1.
        rejected=False
        try:eng.esim_cell_smoke(matlab.double([[1000.]]),matlab.double([params]),matlab.double([knots]),matlab.double(coeff),nargout=3)
        except matlab.engine.MatlabExecutionError:rejected=True
        if not rejected:raise RuntimeError('Nonconverged ESIM cell accepted')
        record.update(accepted=True,cases=cases,nonconvergence_rejected=True,matlab=eng.version(),binary_sha256=hashlib.sha256((OUT/'esim_cell_smoke.mexw64').read_bytes()).hexdigest())
    except BaseException as error:
        record['error']=str(error);raise
    finally:
        if eng is not None:eng.quit()
        (OUT/'result.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
