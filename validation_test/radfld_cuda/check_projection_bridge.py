"""Validate actual native target replay against independent NGSolve Set."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,sys,subprocess
import numpy as np
from element_cases import REPRESENTATIONS, element_case

def source_case(representation, count, placement):
 if representation=='cloud':
  rng=np.random.default_rng(20261012);xyz=np.full(3,-20. if placement=='far' else -.5)+rng.uniform(-.2,.2,(count,3));q=rng.uniform(.5,1.5,count)
  return None,xyz,q
 spec,_=element_case(representation,count,0)
 if placement=='far':
  # Translate geometry and polynomial origins together; coefficients unchanged.
  curved=representation=='curved_tet'
  for key,columns in [('volume',30 if curved else 12),('surface',18 if curved else 9)]:
   spec[key][:,:columns]-=20.
  for key in ('volume_origins','surface_origins'):
   if key in spec:spec[key]-=20.
 return spec,None,None

def main():
 p=argparse.ArgumentParser();p.add_argument('build_dir',type=Path);p.add_argument('--field-build',type=Path);p.add_argument('--reference-dir',type=Path);p.add_argument('--algorithm',choices=['direct','tree','tree_mixed32'],default='direct');p.add_argument('--representation',choices=('cloud',)+REPRESENTATIONS,default='cloud');p.add_argument('--sources',type=int,default=256);p.add_argument('--placement',choices=['far','native'],default='far');p.add_argument('--result',type=Path);a=p.parse_args();root=Path(__file__).resolve().parents[2];out=a.build_dir
 if not 1<=a.sources<=4096:p.error('Small numerical checks require 1..4096 sources')
 if bool(a.field_build)!=bool(a.reference_dir):p.error('--field-build and --reference-dir must be supplied together')
 ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
 import ngsolve as ng
 ng.SetNumThreads(1)
 from netgen.csg import unit_cube
 job=out/'validation';job.mkdir(exist_ok=True);meshng=unit_cube.GenerateMesh(maxh=.7);meshng.Save(str(job/'unit.vol'));mesh=ng.Mesh(meshng)
 refs={};spaces={}
 with ng.TaskManager():
  for name,space in [('hdiv',ng.HDiv(mesh,order=1)),('hcurl',ng.HCurl(mesh,order=1)),('vectorh1',ng.VectorH1(mesh,order=1))]:
   spaces[name]=space;gf=ng.GridFunction(space);gf.Set(ng.CF((ng.x,ng.y,ng.z)));refs[name]=gf.vec.FV().NumPy().copy()
 ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
 prefix=Path(sys.base_prefix);paths=[prefix,prefix/'bin',prefix/'Lib/site-packages/netgen',prefix/'Lib/site-packages/netgen_mesher.libs',prefix/'Lib/site-packages/ngsolve.libs',prefix/'Lib/site-packages/ngsolve_openblas']
 os.environ['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+os.environ['PATH']
 cases=[]
 record=dict(accepted=False,scope='Actual SetValues projection with independently evaluated RadFld sources; not whole application acceptance',representation=a.representation,placement=a.placement,cases=cases,timing_acceptance=False,a100_used=False,ngsolve=ng.__version__,binary_sha256=hashlib.sha256((out/'radia_cuda_projection.mexw64').read_bytes()).hexdigest());eng=None
 record.update(validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),relative_limit=1e-4,whole_transfer_inclusive_seconds=None,memory_limitations='Exact owned device buffers and cudaMemGetInfo samples only; host/process/device global peak not measured. Incidental per-call times are diagnostic, not a continuous full-operation timing.')
 try:
  import matlab.engine
  eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread')
  pid=int(eng.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
  eng.addpath(str(out),str(Path(__file__).parent),nargout=0);eng.workspace['dependency_path']=os.pathsep.join(map(str,paths));eng.eval("setenv('PATH',[dependency_path pathsep getenv('PATH')]);",nargout=0);eng.workspace['vol_path']=str(job/'unit.vol')
  for name,ref in refs.items():
   eng.workspace['kind']=name;eng.eval("h=radia_cuda_projection('create',vol_path,kind,1);",nargout=0)
   try:
    eng.eval("[v,cache,info]=projectBatchedField(h,@(p)p,\"analytic-v1\");[v2,cache,reused]=projectBatchedField(h,@(p)error('cache miss'),\"analytic-v1\",cache);assert(isequal(v,v2)&&reused.evaluation.reused);",nargout=0)
    actual=np.asarray(eng.workspace['v']).ravel();error=float(np.linalg.norm(actual-ref)/max(np.linalg.norm(ref),1e-30));assert error<1e-12
    cases.append(dict(space=name,relative_coefficient_error=error,repeat_equal=True,targets=int(eng.workspace['info']['target_count'])))
    if a.field_build and a.reference_dir:
     sys.path.insert(0,str(a.reference_dir));import _radia_field_reference as oracle
     cpu_receipt=json.loads((a.reference_dir/'reference.json').read_text(encoding='utf-8'));assert hashlib.sha256(Path(oracle.__file__).read_bytes()).hexdigest()==cpu_receipt['binary_sha256']
     field_receipt=json.loads((a.field_build/'build.json').read_text(encoding='utf-8'));assert hashlib.sha256((a.field_build/'radia_cuda_field.mexw64').read_bytes()).hexdigest()==field_receipt['binary_sha256']
     record['field_build']=field_receipt;record['cpu_field_reference']=cpu_receipt
     spec,xyz,q=source_case(a.representation,a.sources,a.placement);empty=np.zeros(0)
     record['input_sha256']={key:hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest() for key,value in (spec.items() if spec else [('xyz',xyz),('strength',q)]) if isinstance(value,np.ndarray)}
     record['input_sha256']['vol']=hashlib.sha256((job/'unit.vol').read_bytes()).hexdigest()
     with ng.TaskManager():
      f=oracle.Field.from_elements(spec,empty,empty,empty,8,.05) if spec else oracle.Field(xyz,q,empty,empty,empty,8,.05)
      eng.eval("mapped_targets=radia_cuda_projection('targets',h);",nargout=0);points=np.asarray(eng.workspace['mapped_targets']);expected=f.field(points,'direct')
     import matlab
     eng.workspace['cpu_values']=matlab.double(expected.tolist())
     if spec:
      eng.workspace['element_source']={key:matlab.double(value.tolist()) if isinstance(value,np.ndarray) else value for key,value in spec.items()}
      create="gpu_field=radia.gpu.FieldElements(element_source,LeafSize=8);"
     else:
      eng.workspace['source_xyz']=matlab.double(xyz.tolist());eng.workspace['source_q']=matlab.double(q[:,None].tolist());create="gpu_field=radia.gpu.FieldCloud(source_xyz,source_q,LeafSize=8);"
     eng.workspace['field_algorithm']=a.algorithm;eng.addpath(str(a.field_build),str(root/'packages/radia-gpu/matlab'),nargout=0)
     eng.eval("[cpu_vector,~]=radia_cuda_projection('project',h,cpu_values);"+create,nargout=0)
     try:
      eng.eval('[probe_values,device_memory]=gpu_field.field(mapped_targets,field_algorithm);',nargout=0)
      probe=np.asarray(eng.workspace['probe_values']);field_error=float(np.linalg.norm(probe-expected)/max(np.linalg.norm(expected),1e-30))
      cases[-1]['relative_field_error']=field_error;cases[-1]['device_memory']={key:float(value) for key,value in eng.workspace['device_memory'].items() if isinstance(value,(int,float))}
      assert field_error<=(1e-12 if a.algorithm=='direct' else 1e-4)
      if a.algorithm=='tree_mixed32':
       with ng.TaskManager():cpu_tree=f.field(points,'tree')
       eng.eval("[fp64_values,fp64_memory]=gpu_field.field(mapped_targets,'tree');[fp64_vector,~]=radia_cuda_projection('project',h,fp64_values);",nargout=0)
       fp64=np.asarray(eng.workspace['fp64_values']);parity=float(np.linalg.norm(fp64-cpu_tree)/max(np.linalg.norm(cpu_tree),1e-30));assert parity<=1e-12
       added=float(np.linalg.norm(probe-fp64)/max(np.linalg.norm(fp64),1e-30));assert added<=1e-4
       cases[-1].update(fp64_tree_CPU_parity=parity,relative_field_difference_from_FP64_tree=added,fp64_tree_relative_to_CPU_direct=float(np.linalg.norm(fp64-expected)/max(np.linalg.norm(expected),1e-30)),fp64_tree_device_memory={key:float(value) for key,value in eng.workspace['fp64_memory'].items() if isinstance(value,(int,float))})
      eng.eval('[gpu_vector,field_cache,field_info]=projectBatchedField(h,@(p)gpu_field.field(p,field_algorithm),\"cloud-v1\");[again,field_cache,reuse_info]=projectBatchedField(h,@(p)error(\"unexpected evaluation\"),\"cloud-v1\",field_cache);assert(isequal(again,gpu_vector));',nargout=0)
      gv=np.asarray(eng.workspace['gpu_vector']).ravel();cv=np.asarray(eng.workspace['cpu_vector']).ravel();e=float(np.linalg.norm(gv-cv)/max(np.linalg.norm(cv),1e-30));assert e<=(1e-12 if a.algorithm=='direct' else 1e-4)
      with ng.TaskManager():
       space=spaces[name];u,v=space.TnT();mass=ng.BilinearForm(space);mass+=ng.SymbolicBFI(ng.InnerProduct(u,v));mass.Assemble();d=ng.GridFunction(space);c=ng.GridFunction(space);d.vec.FV().NumPy()[:]=gv-cv;c.vec.FV().NumPy()[:]=cv;l2=float(np.sqrt(ng.InnerProduct(d.vec,mass.mat*d.vec)/max(ng.InnerProduct(c.vec,mass.mat*c.vec),1e-300)));assert l2<=1e-4
      cases[-1].update(algorithm=a.algorithm,source_count=a.sources*(2 if spec else 1),relative_FE_L2_error=l2)
      cases[-1]['gpu_projection_relative_error']=e;cases[-1]['gpu_reuse_passed']=True
     finally:eng.eval('delete(gpu_field)',nargout=0)

   finally:eng.eval("radia_cuda_projection('destroy',h)",nargout=0)
  assert eng.radia_cuda_projection('live_count')==0
  record.update(accepted=True,cases=cases,contexts_released=True)
 except BaseException as e:record['error']=str(e);raise
 finally:
  if eng is not None:eng.quit()
  (a.result or out/'validation.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
