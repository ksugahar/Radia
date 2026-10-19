"""Progressive actual projection size, precision, batching and reuse diagnostics."""
from pathlib import Path
import tempfile
import argparse,ctypes,hashlib,importlib.util,json,os,statistics,subprocess,sys,time
import numpy as np
import psutil
from check_projection_bridge import source_case

ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--sources',type=int,required=True);p.add_argument('--maxh',type=float,required=True);p.add_argument('--representation',choices=['cloud','tet','curved_tet'],default='cloud');p.add_argument('--numerical-only',action='store_true');a=p.parse_args()
 if not 1<=a.sources<=4096 or not .08<=a.maxh<=.7:p.error('Use bounded incremental cases (sources<=4096, maxh>=.08)')
 if a.numerical_only and (a.sources>1024 or a.maxh<.25):p.error('Non-timing checks are limited to 1024 sources and maxh>=.25')
 if a.numerical_only:ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
 a.out.mkdir(exist_ok=False,parents=True)
 field=ROOT/'build-field-cuda/job-45q9rhs_';projection=Path(tempfile.gettempdir())/'projection-cuda-l8o_j38p';reference=Path(tempfile.gettempdir())/'radfld-elements-reference-20261012-1910'
 receipts={label:json.loads((folder/name).read_text(encoding='utf-8')) for label,folder,name in [('field',field,'build.json'),('projection',projection,'build.json'),('CPU',reference,'reference.json')]}
 for label,folder,binary in [('field',field,'radia_cuda_field.mexw64'),('projection',projection,'radia_cuda_projection.mexw64'),('CPU',reference,'_radia_field_reference.pyd')]:assert sha(folder/binary)==receipts[label]['binary_sha256']
 spec=importlib.util.spec_from_file_location('idle_guard',ROOT/'validation_test/ams_cuda/graph_validation.py');guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
 record=dict(schema='radia.projection-scaling-reuse.v1',accepted=False,timing_acceptance=False,classification='Idle-100 RTX diagnostic, not gpu evaluation worker performance acceptance',a100_used=False,sources=a.sources,representation=a.representation,maxh=a.maxh,receipts=receipts,source_sha256={str(q.relative_to(ROOT)):sha(q) for q in [Path(__file__),Path(__file__).with_name('profileProjectionReuse.m'),Path(__file__).with_name('check_projection_bridge.py'),Path(__file__).with_name('projectBatchedField.m'),Path(__file__).with_name('evaluateProjectionTargets.m')]})
 if a.numerical_only:record['classification']='Nonidle small numerical check; timings invalid and excluded'
 def idle():
  if not a.numerical_only:guard.assert_idle_host()
 eng=None
 try:
  idle();cpu=psutil.cpu_percent(interval=2,percpu=True);record['idle_cpu_mean_percent']=sum(cpu)/len(cpu)
  if record['idle_cpu_mean_percent']>=(50 if a.numerical_only else 10):raise RuntimeError('Host CPU load exceeds the selected diagnostic budget')
  assert psutil.virtual_memory().available>8*1024**3
  import ngsolve as ng
  from netgen.csg import unit_cube
  ng.SetNumThreads(1);meshng=unit_cube.GenerateMesh(maxh=a.maxh);vol=a.out/'unit.vol';meshng.Save(str(vol));record['mesh_sha256']=sha(vol);mesh=ng.Mesh(meshng)
  prefix=Path(sys.base_prefix);paths=[prefix,prefix/'bin',prefix/'Lib/site-packages/netgen',prefix/'Lib/site-packages/netgen_mesher.libs',prefix/'Lib/site-packages/ngsolve.libs',prefix/'Lib/site-packages/ngsolve_openblas'];os.environ['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+os.environ['PATH']
  sys.path.insert(0,str(reference));import _radia_field_reference as oracle
  import matlab.engine
  import matlab
  eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');eng.addpath(str(field),str(projection),str(ROOT/'packages/radia-simulink/matlab'),str(Path(__file__).parent),nargout=0)
  if a.numerical_only:
   pid=int(eng.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
  eng.workspace['dependency_path']=os.pathsep.join(map(str,paths));eng.eval("setenv('PATH',[dependency_path pathsep getenv('PATH')]);",nargout=0)
  record['matlab']=eng.version();record['gpu']=json.loads(eng.eval("jsonencode(struct('name',gpuDevice().Name,'total_memory',gpuDevice().TotalMemory,'available_memory',gpuDevice().AvailableMemory))"));assert record['gpu']['available_memory']>max(2*1024**3,.2*record['gpu']['total_memory'])
  eng.workspace['check_path']=str(Path(__file__).with_name('profileProjectionReuse.m'));record['code_analysis']=json.loads(eng.eval("jsonencode(checkcode(check_path,'-id'))"))
  assert not [x for x in (record['code_analysis'] if isinstance(record['code_analysis'],list) else [record['code_analysis']]) if x and x.get('id') not in ('MOCUP','NASGU')]
  eng.workspace['vol_path']=str(vol);eng.eval("h=radia_cuda_projection('create',vol_path,'hdiv',1);targets=radia_cuda_projection('targets',h);radia_cuda_projection('destroy',h);",nargout=0);targets=np.asarray(eng.workspace['targets']);record['targets']=len(targets)
  source,xyz,q=source_case(a.representation,a.sources,'far');empty=np.zeros(0)
  arrays={key:value for key,value in source.items() if isinstance(value,np.ndarray)} if source else dict(xyz=xyz,strength=q);record['input_sha256']={key:hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest() for key,value in arrays.items()};record['targets_sha256']=hashlib.sha256(targets.tobytes()).hexdigest()
  with ng.TaskManager():f=oracle.Field.from_elements(source,empty,empty,empty,8,.05) if source else oracle.Field(xyz,q,empty,empty,empty,8,.05);expected=f.field(targets,'direct')
  if source:eng.workspace['source']={key:matlab.double(value.tolist()) if isinstance(value,np.ndarray) else value for key,value in source.items()}
  else:eng.workspace['source']={'representation':'cloud','xyz':matlab.double(xyz.tolist()),'strength':matlab.double(q[:,None].tolist())}
  eng.workspace['cpu_values']=matlab.double(expected.tolist());eng.workspace['json_path']=str(a.out/'matlab.json')
  idle();eng.eval('[measurements,vectors]=profileProjectionReuse(vol_path,source,cpu_values,json_path);',nargout=0);idle()
  measurements=json.loads((a.out/'matlab.json').read_text(encoding='utf-8'));assert measurements['accepted'];record['measurements']=measurements
  eng.eval("h=radia_cuda_projection('create',vol_path,'hdiv',1);[cpu_vector,~]=radia_cuda_projection('project',h,cpu_values);radia_cuda_projection('destroy',h);assert(radia_cuda_projection('live_count')==0);",nargout=0);cv=np.asarray(eng.workspace['cpu_vector']).ravel();vectors=np.asarray(eng.workspace['vectors'])
  with ng.TaskManager():
   space=ng.HDiv(mesh,order=1);u,v=space.TnT();mass=ng.BilinearForm(space);mass+=ng.SymbolicBFI(ng.InnerProduct(u,v));mass.Assemble();base=ng.GridFunction(space);diff=ng.GridFunction(space);base.vec.FV().NumPy()[:]=cv;den=max(ng.InnerProduct(base.vec,mass.mat*base.vec),1e-300);errors=[]
   for column in range(vectors.shape[1]):diff.vec.FV().NumPy()[:]=vectors[:,column]-cv;errors.append(float(np.sqrt(ng.InnerProduct(diff.vec,mass.mat*diff.vec)/den)))
  record['projected_FE_L2_relative_errors']=errors;assert max(errors)<=1e-4
  record['median_amortized_step_seconds']={mode:statistics.median(s['amortized_step_seconds'] for s in measurements['samples'] if s['mode']==mode) for mode in ('rebuild64','chunk64','batch64','batch32','cache64','cache32')}
  if a.numerical_only:
   record['median_amortized_step_seconds']=None
   measurements.pop('precision_validation_seconds',None);measurements.pop('warmup_sequences',None)
   for sample in measurements['samples']:
    sample.pop('whole_seconds',None);sample.pop('amortized_step_seconds',None)
    sample.get('device_memory',{}).pop('total_seconds',None)
  record.update(accepted=True,contexts_released=True,limitations=['Far synthetic source placement, actual order-1 HDiv projection; not general physical/scaling acceptance.','Cache comparisons reuse unchanged fields only; rebuild is required after source changes.','Host/native peak unmeasured; device owned buffers and sampled free memory are distinguished.'])
 except BaseException as error:record['error']=f'{type(error).__name__}: {error}';raise
 finally:
  if eng is not None:eng.quit()
  (a.out/'result.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
