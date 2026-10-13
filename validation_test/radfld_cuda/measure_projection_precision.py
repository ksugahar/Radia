"""Idle-100 diagnostic ABBA through MATLAB; Python is validation-only."""
from pathlib import Path
import argparse, hashlib, importlib.util, json, os, statistics, subprocess, sys, time
import numpy as np
import psutil
from check_projection_bridge import source_case
from element_cases import REPRESENTATIONS

ROOT=Path(__file__).resolve().parents[2]
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 p=argparse.ArgumentParser();p.add_argument('--field-build',type=Path,required=True);p.add_argument('--projection-build',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 a.out.mkdir(parents=True,exist_ok=False)
 spec=importlib.util.spec_from_file_location('idle_guard',ROOT/'validation_test/ams_cuda/graph_validation.py');guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
 evidence=json.loads((Path(__file__).parent/'results/elements-projection-precision-20261012.json').read_text(encoding='utf-8'))
 field=json.loads((a.field_build/'build.json').read_text(encoding='utf-8'));projection=json.loads((a.projection_build/'build.json').read_text(encoding='utf-8'))
 assert digest(a.field_build/'radia_cuda_field.mexw64')==field['binary_sha256']==evidence['common_provenance']['field_build']['binary_sha256']
 assert digest(a.projection_build/'radia_cuda_projection.mexw64')==projection['binary_sha256']==evidence['common_provenance']['binary_sha256']
 prefix=Path(sys.base_prefix);paths=[prefix,prefix/'bin',prefix/'Lib/site-packages/netgen',prefix/'Lib/site-packages/netgen_mesher.libs',prefix/'Lib/site-packages/ngsolve.libs',prefix/'Lib/site-packages/ngsolve_openblas']
 os.environ['PATH']=os.pathsep.join(map(str,paths))+os.pathsep+os.environ['PATH']
 record=dict(schema='radia.projection-precision-timing.v1',accepted=False,timing_acceptance=False,classification='100 RTX 5000 idle diagnostic; not mdx performance acceptance',a100_used=False,field_build=field,projection_build=projection,cpu_precision_evidence_sha256=digest(Path(__file__).parent/'results/elements-projection-precision-20261012.json'),source_sha256={str(q.relative_to(ROOT)):digest(q) for q in [Path(__file__),Path(__file__).with_name('compareProjectionPrecision.m'),Path(__file__).with_name('check_projection_bridge.py'),Path(__file__).with_name('element_cases.py')]},cases=[])
 eng=None
 record['validation_commit']=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip();record['validation_dirty']=bool(subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True).strip())
 try:
  guard.assert_idle_host()
  cpus=psutil.cpu_percent(interval=2,percpu=True);record['idle_probe']=dict(logical_processors=len(cpus),mean_cpu_percent=sum(cpus)/len(cpus),max_cpu_percent=max(cpus),available_host_bytes=psutil.virtual_memory().available,processor_load_caveat='Socket CIM reports conflicted with per-logical-CPU samples; do not use the socket average alone.')
  if record['idle_probe']['mean_cpu_percent']>10:raise RuntimeError('Host CPU load is too high for timing')
  import matlab.engine
  import matlab
  started=time.perf_counter();eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');record['matlab_startup_seconds']=time.perf_counter()-started
  eng.addpath(str(a.field_build),str(a.projection_build),str(ROOT/'packages/radia-gpu/matlab'),str(Path(__file__).parent),nargout=0)
  eng.workspace['dependency_path']=os.pathsep.join(map(str,paths));eng.eval("setenv('PATH',[dependency_path pathsep getenv('PATH')]);",nargout=0)
  eng.workspace['check_path']=str(Path(__file__).with_name('compareProjectionPrecision.m'));record['code_analysis']=json.loads(eng.eval("jsonencode(checkcode(check_path,'-id'))"));assert not record['code_analysis']
  record['matlab']=eng.version();record['gpu']=json.loads(eng.eval("jsonencode(struct('name',gpuDevice().Name,'total_memory',gpuDevice().TotalMemory,'available_memory',gpuDevice().AvailableMemory))"))
  if record['gpu']['available_memory']<max(2*1024**3,.2*record['gpu']['total_memory']):raise RuntimeError('GPU reserve not available')
  vol=a.projection_build/'validation/unit.vol';record['mesh_sha256']=digest(vol);eng.workspace['vol_path']=str(vol)
  assert all(record['mesh_sha256']==check['input_sha256']['vol'] for check in evidence['checks'])
  for kind in ('cloud',)+REPRESENTATIONS:
   guard.assert_idle_host()
   source,xyz,q=source_case(kind,256 if kind=='cloud' else 16,'far')
   raw=a.out/(kind+'.json');eng.workspace['json_path']=str(raw)
   if source:
    eng.workspace['element_source']={key:matlab.double(value.tolist()) if isinstance(value,np.ndarray) else value for key,value in source.items()}
    command="compareProjectionPrecision(vol_path,zeros(1,3),1,json_path,ElementSource=element_source);"
    arrays={key:value for key,value in source.items() if isinstance(value,np.ndarray)}
   else:
    eng.workspace['source_xyz']=matlab.double(xyz.tolist());eng.workspace['source_q']=matlab.double(q[:,None].tolist());command='compareProjectionPrecision(vol_path,source_xyz,source_q,json_path);';arrays=dict(xyz=xyz,strength=q)
   eng.eval(command,nargout=0)
   case=json.loads(raw.read_text(encoding='utf-8'));assert case['accepted'];case['representation']=kind;case['input_sha256']={key:hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest() for key,value in arrays.items()}
   case['median_whole_seconds']={alg:statistics.median(s['whole_seconds'] for s in case['samples'] if s['algorithm']==alg) for alg in ('tree','tree_mixed32')}
   case['FP64_over_mixed32_ratio']=case['median_whole_seconds']['tree']/case['median_whole_seconds']['tree_mixed32'];record['cases'].append(case)
   guard.assert_idle_host()
  eng.eval("assert(radia_cuda_projection('live_count')==0);",nargout=0)
  record.update(accepted=True,contexts_released=True,limitations=['Small fixed source and projection cases; not scaling or broad performance acceptance.','Memory: exact owned device buffers plus CUDA free-memory samples, not process/global peak; host peak unmeasured.','Cold samples are first calls per algorithm in one MATLAB session; second algorithm shares initialized CUDA runtime.'])
 except BaseException as error:record['error']=str(error);raise
 finally:
  if eng is not None:eng.quit()
  (a.out/'result.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
