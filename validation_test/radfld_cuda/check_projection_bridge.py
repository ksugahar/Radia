"""Validate actual native target replay against independent NGSolve Set."""
from pathlib import Path
import argparse,ctypes,hashlib,json,os,sys
import numpy as np

def main():
 p=argparse.ArgumentParser();p.add_argument('build_dir',type=Path);p.add_argument('--field-build',type=Path);p.add_argument('--reference-dir',type=Path);p.add_argument('--algorithm',choices=['direct','tree','tree_mixed32'],default='direct');a=p.parse_args();root=Path(__file__).resolve().parents[2];out=a.build_dir
 import ngsolve as ng
 from netgen.csg import unit_cube
 job=out/'validation';job.mkdir(exist_ok=True);meshng=unit_cube.GenerateMesh(maxh=.7);meshng.Save(str(job/'unit.vol'));mesh=ng.Mesh(meshng)
 refs={};spaces={}
 with ng.TaskManager():
  for name,space in [('hdiv',ng.HDiv(mesh,order=1)),('hcurl',ng.HCurl(mesh,order=1)),('vectorh1',ng.VectorH1(mesh,order=1))]:
   spaces[name]=space;gf=ng.GridFunction(space);gf.Set(ng.CF((ng.x,ng.y,ng.z)));refs[name]=gf.vec.FV().NumPy().copy()
 ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
 prefix=Path(sys.base_prefix);paths=[prefix,prefix/'bin',prefix/'Lib/site-packages/netgen',prefix/'Lib/site-packages/netgen_mesher.libs',prefix/'Lib/site-packages/ngsolve.libs',prefix/'Lib/site-packages/ngsolve_openblas']
 os.environ['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+os.environ['PATH']
 record=dict(accepted=False,scope='Analytic vector coefficient target replay, not RadFld physical/GPU acceptance',timing_acceptance=False,a100_used=False,ngsolve=ng.__version__,binary_sha256=hashlib.sha256((out/'radia_cuda_projection.mexw64').read_bytes()).hexdigest());eng=None
 try:
  import matlab.engine
  eng=matlab.engine.start_matlab('-nodesktop -nosplash');eng.addpath(str(out),str(Path(__file__).parent),nargout=0);eng.workspace['dependency_path']=os.pathsep.join(map(str,paths));eng.eval("setenv('PATH',[dependency_path pathsep getenv('PATH')]);",nargout=0);eng.workspace['vol_path']=str(job/'unit.vol');cases=[]
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
     rng=np.random.default_rng(20261012);xyz=np.array([-20.,-20.,-20.])+rng.uniform(-.2,.2,(256,3));q=rng.uniform(.5,1.5,256);empty=np.zeros(0)
     with ng.TaskManager():
      f=oracle.Field(xyz,q,empty,empty,empty,8,.05)
      eng.eval("mapped_targets=radia_cuda_projection('targets',h);",nargout=0);points=np.asarray(eng.workspace['mapped_targets']);expected=f.field(points,'direct')
     import matlab
     eng.workspace['cpu_values']=matlab.double(expected.tolist());eng.workspace['source_xyz']=matlab.double(xyz.tolist());eng.workspace['source_q']=matlab.double(q[:,None].tolist())
     eng.workspace['field_algorithm']=a.algorithm;eng.addpath(str(a.field_build),str(root/'packages/radia-gpu/matlab'),nargout=0)
     eng.eval("[cpu_vector,~]=radia_cuda_projection('project',h,cpu_values);gpu_field=radia.gpu.FieldCloud(source_xyz,source_q,LeafSize=8);",nargout=0)
     try:
      eng.eval('[gpu_vector,field_cache,field_info]=projectBatchedField(h,@(p)gpu_field.field(p,field_algorithm),\"cloud-v1\");[again,field_cache,reuse_info]=projectBatchedField(h,@(p)error(\"unexpected evaluation\"),\"cloud-v1\",field_cache);assert(isequal(again,gpu_vector));',nargout=0)
      gv=np.asarray(eng.workspace['gpu_vector']).ravel();cv=np.asarray(eng.workspace['cpu_vector']).ravel();e=float(np.linalg.norm(gv-cv)/max(np.linalg.norm(cv),1e-30));assert e<=(1e-12 if a.algorithm=='direct' else 1e-4)
      with ng.TaskManager():
       space=spaces[name];u,v=space.TnT();mass=ng.BilinearForm(space);mass+=ng.SymbolicBFI(ng.InnerProduct(u,v));mass.Assemble();d=ng.GridFunction(space);c=ng.GridFunction(space);d.vec.FV().NumPy()[:]=gv-cv;c.vec.FV().NumPy()[:]=cv;l2=float(np.sqrt(ng.InnerProduct(d.vec,mass.mat*d.vec)/max(ng.InnerProduct(c.vec,mass.mat*c.vec),1e-300)));assert l2<=1e-4
      cases[-1].update(algorithm=a.algorithm,source_count=len(q),relative_FE_L2_error=l2)
      cases[-1]['retained_cloud_gpu_projection_relative_error']=e;cases[-1]['gpu_reuse_passed']=True
     finally:eng.eval('delete(gpu_field)',nargout=0)

   finally:eng.eval("radia_cuda_projection('destroy',h)",nargout=0)
  assert eng.radia_cuda_projection('live_count')==0
  record.update(accepted=True,cases=cases,contexts_released=True)
 except BaseException as e:record['error']=str(e);raise
 finally:
  if eng is not None:eng.quit()
  (out/'validation.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
