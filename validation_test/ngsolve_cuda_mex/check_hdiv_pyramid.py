"""Bounded single-pyramid HDiv construction and mass/div-div assembly probe.

Not mixed-mesh conformity, interpolation, MMM physics or GPU acceptance.
Pass the output JSON path as the sole argument. No runtime kernel compilation.
"""
import json,time,traceback,sys,importlib
from pathlib import Path
out=Path(sys.argv[1]);result={'schema':'radia.ngsolve2608-probe.v1','accepted':False,'cases':[]}
try:
 import ngsolve as ng,netgen
 from netgen.meshing import Mesh,MeshPoint,Pnt,Element3D,Element2D,FaceDescriptor
 ng.SetNumThreads(2)
 result.update(ngsolve=ng.__version__,ngsolve_path=ng.__file__,netgen_path=netgen.__file__)
 for order in [0,1,2,3]:
  rec={'order':order};start=time.perf_counter()
  try:
   nm=Mesh(dim=3);nm.SetMaterial(1,'domain')
   ids=[nm.Add(MeshPoint(Pnt(*p))) for p in [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1)]]
   nm.Add(Element3D(1,ids))
   for i,face in enumerate([[0,3,2,1],[0,1,4],[1,2,4],[2,3,4],[3,0,4]],1):
    fd=nm.Add(FaceDescriptor(surfnr=i,domin=1,domout=0,bc=i));nm.SetBCName(i-1,'boundary');nm.Add(Element2D(fd,[ids[j] for j in face]))
   mesh=ng.Mesh(nm);fes=ng.HDiv(mesh,order=order);u,v=fes.TnT()
   with ng.TaskManager():
    a=ng.BilinearForm(fes,symmetric=True);a+=ng.InnerProduct(u,v)*ng.dx+ng.div(u)*ng.div(v)*ng.dx;a.Assemble()
   vals,cols,rows=a.mat.CSR()
   import numpy as np
   values=np.asarray(vals);assert np.isfinite(values).all()
   rec.update(accepted=True,ndof=fes.ndof,nnz=len(values),volume=float(ng.Integrate(1,mesh)),element_types=[str(e.type) for e in mesh.Elements()])
  except Exception as e:rec.update(accepted=False,error=repr(e),traceback=traceback.format_exc())
  rec['seconds']=time.perf_counter()-start;result['cases'].append(rec)
 try:
  gpu=importlib.import_module('ngsolve.gpu');dev=gpu.GetGPUDevice()
  result['commonGPU']={'module':gpu.__file__,'backend':str(gpu.backend),'device':str(dev),'has_float64':bool(dev.has_float64),'runtime_compile_executed':False}
 except Exception as e:result['commonGPU']={'import_error':repr(e)}
 result['accepted']=all(x['accepted'] for x in result['cases'])
finally:
 out.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
