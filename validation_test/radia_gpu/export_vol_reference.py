"""Export an NGSolve CPU reference from a hash-verified .vol input bundle.

Usage: python export_vol_reference.py BUNDLE_DIRECTORY
The bundle requires poisson.vol, labels.json, mesh-check.json and identity.json.
The producer runs check-vol before transfer; no Cubit runtime is required here.
This Poisson case assumes the explicitly declared unit-cube geometry.
"""
from pathlib import Path
import json,hashlib,math,sys,time,statistics
import numpy as np
import ngsolve as ng
from scipy.sparse import coo_matrix
from scipy.io import savemat
base=Path(sys.argv[1]).resolve() if len(sys.argv)==2 else Path(__file__).resolve().parent
identity=json.loads((base/'identity.json').read_text())
for name,digest in identity['inputs'].items():
 assert hashlib.sha256((base/name).read_bytes()).hexdigest()==digest,name
check=json.loads((base/'mesh-check.json').read_text());assert check['passed']
ng.SetNumThreads(1)
with ng.TaskManager():
 mesh=ng.Mesh(str(base/'poisson.vol'))
 coordinates=np.asarray([tuple(vertex.point) for vertex in mesh.vertices])
 if not (np.allclose(coordinates.min(axis=0),0,atol=1e-12,rtol=0) and np.allclose(coordinates.max(axis=0),1,atol=1e-12,rtol=0) and abs(ng.Integrate(1,mesh)-1)<1e-12):
  raise ValueError('This declared Poisson reference requires the unit cube')
 contract=json.loads((base/'labels.json').read_text())
 assert set(mesh.GetMaterials())==set(contract['required']['materials'])
 assert set(mesh.GetBoundaries())==set(contract['required']['boundaries'])
 fes=ng.H1(mesh,order=2,dirichlet='.*');u,v=fes.TnT()
 exact=ng.sin(math.pi*ng.x)*ng.sin(math.pi*ng.y)*ng.sin(math.pi*ng.z)
 a=ng.BilinearForm(ng.grad(u)*ng.grad(v)*ng.dx).Assemble()
 f=ng.LinearForm(3*math.pi**2*exact*v*ng.dx).Assemble()
 reference=f.vec.CreateVector()
 started=time.perf_counter();inverse=a.mat.Inverse(fes.FreeDofs(),inverse='sparsecholesky');factor_seconds=time.perf_counter()-started
 solve_seconds=[]
 for _ in range(6):
  started=time.perf_counter();reference.data=inverse*f.vec;solve_seconds.append(time.perf_counter()-started)
 free=np.asarray(list(fes.FreeDofs()),dtype=bool)
 rows,columns,values=a.mat.COO();matrix=coo_matrix((np.asarray(values),(np.asarray(rows),np.asarray(columns))),shape=(fes.ndof,fes.ndof)).tocsr()[free][:,free].tocsr();matrix.sort_indices()
 rhs=f.vec.FV().NumPy()[free].copy();x=reference.FV().NumPy()[free].copy()
 residual=float(np.linalg.norm(matrix@x-rhs)/np.linalg.norm(rhs));assert residual<1e-6
 savemat(base/'system.mat',{'rowptr':matrix.indptr.astype(np.int32),'columns':matrix.indices.astype(np.int32),'values':matrix.data,'rhs':rhs,'diagonal_inverse':1/matrix.diagonal(),'reference':x},oned_as='column')
report={'schema':'radia.gpu.ngsolve-vol-reference/v1','passed':True,'ngsolve':ng.__version__,'python':sys.version.split()[0],'source_commit':identity['source_commit'],'mesh_sha256':identity['inputs']['poisson.vol'],'mesh_preflight_passed':True,'label_contract_schema':contract['schema'],'cpu_threads':1,'ndof':fes.ndof,'free_dofs':int(free.sum()),'nnz':matrix.nnz,'cpu_true_relative_residual':residual,'system_sha256':hashlib.sha256((base/'system.mat').read_bytes()).hexdigest(),'gpu_checked':False,'timing_checked':True,'formal_performance_checked':False,'cpu_solver':'sparsecholesky','cpu_factor_seconds':factor_seconds,'cpu_cold_solve_seconds':solve_seconds[0],'cpu_warm_solve_seconds':solve_seconds[1:],'cpu_warm_median_seconds':statistics.median(solve_seconds[1:]),'cubit_execution':False,'execution_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(base/'cpu-result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
