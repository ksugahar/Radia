from pathlib import Path
import os,sys,json,time,hashlib,importlib,traceback,cmath,math
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent; SOURCE=ROOT/'source'
SHA='cf77cd6e02ebc1d002af8adb37a0b2cedd58d7d9'
A=[25,49,97,193,385]; B=[1e-3,3e-4,1e-4,3e-5,1e-5]
CONDITIONS=dict(half_depth_m=.005,sigma_S_per_m=5e6,frequency_Hz=222e3,mu_r=100.,H_peak_A_per_m=1000.)
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,data):(ROOT/name).write_text(json.dumps(data,indent=2,allow_nan=False),encoding='utf-8')
def load():
 import ngsolve,netgen,radia,numpy as np,scipy,threadpoolctl
 assert Path(radia.__file__).resolve().is_relative_to(ROOT/'venv')
 ngsolve.SetNumThreads(1)
 radia.__path__.insert(0,str(SOURCE/'src/radia'))
 for name in ('radia.esim_cell_problem','radia.esim_depth'):sys.modules.pop(name,None)
 module=importlib.import_module('radia.esim_cell_problem');depth=importlib.import_module('radia.esim_depth')
 assert Path(module.__file__).resolve()==(SOURCE/'src/radia/esim_cell_problem.py').resolve()
 assert Path(depth.__file__).resolve()==(SOURCE/'src/radia/esim_depth.py').resolve()
 global controller
 controller=threadpoolctl.threadpool_limits(limits=1)
 pools=threadpoolctl.threadpool_info();assert all(p['num_threads']<=1 for p in pools)
 identity=dict(python=sys.executable,prefix=sys.prefix,radia=str(radia.__file__),radia_version=radia.__version__,ngsolve=str(ngsolve.__file__),ngsolve_version=ngsolve.__version__,netgen=str(netgen.__file__),numpy=np.__version__,scipy=scipy.__version__,threadpools=pools,source={str(p.relative_to(SOURCE)):digest(p) for p in (SOURCE/'src/radia/esim_cell_problem.py',SOURCE/'src/radia/esim_depth.py')})
 return module,depth,np,identity

def measure(label,dry=False):
 t=time.perf_counter();module,depth,np,identity=load();cls=module.ESIMFiniteSlabSolver
 from scipy.constants import mu_0
 attrs={n:callable(getattr(cls,n,None)) for n in ('solve','_solve_on_mesh','_solve_linear_system','_compute_power_losses')}
 assert all(attrs.values()) and callable(depth.depth_mesh)
 sys.path.insert(0,str(SOURCE));from tools.regression_lock import load_baseline,assert_regression
 baseline=load_baseline(SOURCE/'tests/fixtures/esim_finite_slab_b49251dac.json')
 def new(n,tol):return cls(.005,sigma=5e6,frequency=222e3,mu_r=100,geometry='slab',n_nodes=n,mesh_tol=tol,max_nodes=32769,max_refinements=16,cell_max_iter=200)
 probe=new(24,1e-3);assert all(hasattr(probe,n) for n in ('mesh_points','rho','delta','n_nodes','mu_initial'))
 if dry:
  save('library-dry-run.json',dict(identity=identity,attributes=attrs,solve_executed=False,baseline_sha256=digest(SOURCE/'tests/fixtures/esim_finite_slab_b49251dac.json'),passed=True));return
 mu=float(mu_0)*100.;rho=1/5e6;gamma=cmath.sqrt(1j*2*math.pi*222e3*mu/rho);zref=rho*gamma*cmath.tanh(gamma*.005)
 out=dict(label=label,source_sha=SHA,identity=identity,conditions=CONDITIONS,mu_SI=mu,reference=dict(real=zref.real,imag=zref.imag,formula='rho*sqrt(i*omega*mu/rho)*tanh(a*sqrt(i*omega*mu/rho))'),A=[],B=[],surrogate=None)
 save(label+'.json',out)
 def one(kind,n,tol):
  solver=new(n,tol);begin=time.perf_counter()
  result=solver._solve_on_mesh(1000.,tol=1e-6,max_iter=200,relaxation=.5) if kind=='A' else solver.solve(1000.)
  z=complex(result['Z']);x=np.asarray(result['mesh_points']);h=float(np.diff(x).max()/.005)
  item=dict(kind=kind,initial_nodes=n,mesh_tol=tol,final_nodes=len(x),h_max_normalized=h,h_min_m=float(np.diff(x).min()),Z_real=z.real,Z_imag=z.imag,error_complex=float(abs(z-zref)/abs(zref)),P_prime=float(result['P_prime']),Q_prime=float(result['Q_prime']),residual=float(getattr(solver,'_last_linear_residual')),converged=bool(result['converged']),mesh_converged=result.get('mesh_converged'),mesh_history=result.get('mesh_history'),runtime_s=time.perf_counter()-begin,mesh_points_m=x.tolist(),H_real=np.real(result['H_solution']).tolist(),H_imag=np.imag(result['H_solution']).tolist())
  assert math.isfinite(item['error_complex']) and math.isfinite(h)
  return item
 try:
  for n in A:
   out['A'].append(one('A',n,1e-3));save(label+'.json',out)
  for tol in B:
   out['B'].append(one('B',24,tol));save(label+'.json',out)
  value=out['B'][0]['Z_real']
  try:assert_regression({'Re_Z_ohm':value},baseline,conditions=baseline['conditions'])
  except AssertionError as exc:out['surrogate']=dict(status='FAIL',value=value,error=str(exc))
  else:out['surrogate']=dict(status='PASS',value=value)
 except Exception:
  out['exception']=traceback.format_exc()
 finally:
  out['elapsed_s']=time.perf_counter()-t;save(label+'.json',out)

if __name__=='__main__':measure(sys.argv[1],len(sys.argv)>2 and sys.argv[2]=='--dry-run')
