from pathlib import Path
import sys,os,json,time,hashlib,cmath,math,traceback
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT/'source';SHA='cf77cd6e02ebc1d002af8adb37a0b2cedd58d7d9'
def save(n,v):(ROOT/n).write_text(json.dumps(v,indent=2,allow_nan=False),encoding='utf-8')
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load():
 import radia,ngsolve,netgen,numpy as np,scipy,threadpoolctl,importlib
 assert Path(radia.__file__).resolve().is_relative_to(ROOT/'venv')
 ngsolve.SetNumThreads(1);radia.__path__.insert(0,str(SOURCE/'src/radia'))
 for name in ('radia.esim_cell_problem','radia.esim_depth'):sys.modules.pop(name,None)
 module=importlib.import_module('radia.esim_cell_problem');depth=importlib.import_module('radia.esim_depth')
 assert Path(module.__file__).resolve()==(SOURCE/'src/radia/esim_cell_problem.py').resolve()
 assert Path(depth.__file__).resolve()==(SOURCE/'src/radia/esim_depth.py').resolve()
 global controller
 controller=threadpoolctl.threadpool_limits(limits=1);pools=threadpoolctl.threadpool_info();assert all(p['num_threads']<=1 for p in pools)
 identity=dict(python=sys.executable,prefix=sys.prefix,radia=str(radia.__file__),radia_version=radia.__version__,ngsolve=str(ngsolve.__file__),ngsolve_version=ngsolve.__version__,netgen=str(netgen.__file__),numpy=np.__version__,scipy=scipy.__version__,threadpools=pools,source={str(p.relative_to(SOURCE)):digest(p) for p in (SOURCE/'src/radia/esim_cell_problem.py',SOURCE/'src/radia/esim_depth.py')})
 return module,depth,np,identity

def profile_error(np,z,H,mu,order):
 # Nodal finite-volume solution is reconstructed as continuous piecewise-linear
 # interpolant. This is not a newly assembled FE solve. Integrate each interval.
 nodes,weights=np.polynomial.legendre.leggauss(order);fraction=(nodes+1)*.5
 dz=np.diff(z);zz=z[:-1,None]+dz[:,None]*fraction[None,:]
 Hh=H[:-1,None]*(1-fraction[None,:])+H[1:,None]*fraction[None,:]
 k=cmath.sqrt(1j*2*math.pi*222e3*mu*5e6);x=.005-zz
 reference=1000*np.cosh(k*x)/cmath.cosh(k*.005)
 w=dz[:,None]*weights[None,:]*.5
 numerator=float(np.sum(w*np.abs(Hh-reference)**2));denominator=float(np.sum(w*np.abs(reference)**2))
 return dict(relative_complex_L2=math.sqrt(numerator/denominator),numerator_integral=numerator,denominator_integral=denominator,quadrature_order=order,k_real=k.real,k_imag=k.imag,x_convention='centre x=a-z; z starts at surface',interpolant='continuous piecewise linear on actual production CV nodes')

def main(label,dry=False):
 start=time.perf_counter();module,depth,np,identity=load();cls=module.ESIMFiniteSlabSolver
 attrs={k:callable(getattr(cls,k,None)) for k in ('_solve_on_mesh','_solve_linear_system','_compute_power_losses')};assert all(attrs.values())
 solver=cls(.005,sigma=5e6,frequency=222e3,mu_r=100,geometry='slab',n_nodes=25)
 assert all(hasattr(solver,n) for n in ('mesh_points','mu_initial','rho','n_nodes'))
 if dry:
  q,w=np.polynomial.legendre.leggauss(32);assert len(q)==32 and abs(float(w.sum())-2)<1e-12
  save('library-dry-run.json',dict(passed=True,solve_executed=False,attributes=attrs,identity=identity));return
 from scipy.constants import mu_0
 mu=float(mu_0)*100;out=dict(label=label,source_sha=SHA,identity=identity,mu_SI=mu,conditions=dict(half_depth_m=.005,sigma=5e6,freq_Hz=222e3,mu_r=100,H0=1000,time_convention='exp(+j omega t)'),levels=[],phase_observation=None)
 save(label+'.json',out)
 try:
  for n in (25,49,97,193,385):
   solver=cls(.005,sigma=5e6,frequency=222e3,mu_r=100,geometry='slab',n_nodes=n)
   t=time.perf_counter();result=solver._solve_on_mesh(1000.,tol=1e-6,max_iter=200,relaxation=.5)
   z=np.asarray(result['mesh_points']);H=np.asarray(result['H_solution']);Z=complex(result['Z'])
   item=dict(initial_nodes=n,final_nodes=len(z),mesh_points_m=z.tolist(),H_real=H.real.tolist(),H_imag=H.imag.tolist(),h_max_normalized=float(np.diff(z).max()/.005),converged=bool(result['converged']),linear_residual=float(solver._last_linear_residual),Z_real=Z.real,Z_imag=Z.imag,first_interior_ImH_over_H0=float(H[1].imag/1000.),runtime_s=time.perf_counter()-t)
   out['levels'].append(item);save(label+'.json',out) # raw result before gates/postprocessing
   item['profile32']=profile_error(np,z,H,mu,32);item['profile16']=profile_error(np,z,H,mu,16);save(label+'.json',out)
  historic=json.loads((ROOT/'historical-A25.json').read_text());old=float(historic['H_imag'][1]/1000.);new=out['levels'][0]['first_interior_ImH_over_H0']
  out['phase_observation']=dict(quantity='Im(H at first interior A25 node)/H0',historic_value=old,actual=new,difference=new-old,relative_difference=abs((new-old)/old),historic_source_sha=SHA,historical_result_sha256=digest(ROOT/'historical-A25.json'),candidate_only=True,baseline_status='proposed; no approved new phase baseline/tolerance; no permanent lock changed',proposed_rtol=1e-12,proposed_atol=0.)
 except Exception:out['exception']=traceback.format_exc()
 finally:out['elapsed_s']=time.perf_counter()-start;save(label+'.json',out)
if __name__=='__main__':main(sys.argv[1],len(sys.argv)>2 and sys.argv[2]=='--dry-run')
