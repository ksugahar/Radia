import os
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
import argparse,ast,hashlib,json,platform,sys,time
from pathlib import Path
root=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description='Compare an explicitly supplied earlier trace helper with the current solve-local cache.')
parser.add_argument('--baseline-source', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args=parser.parse_args()
sys.path.insert(0,str(root/'tests'))
import radia,ngsolve as ng,numpy as np
from scipy.sparse import coo_matrix,diags
from radia.kelvin_solver import _matching_trace_direct_inverse as updated
from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin
import test_kelvin_mixed_omega as fixture
tree=ast.parse(args.baseline_source.read_text())
node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='_matching_trace_direct_inverse')
ns={'np':np}
exec(compile(ast.Module(body=[node],type_ignores=[]),str(args.baseline_source),'exec'),ns)
baseline=ns['_matching_trace_direct_inverse']
ng.SetNumThreads(8);ng.SetHeapSize(10000000)
record=dict(host=platform.node(),ngsolve=ng.__version__,threads=8,cases=[],
    method='Identical matrix graphs, three different positive primal coefficients; fresh baseline vs cache updates',
    runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    old_source_sha256=hashlib.sha256(args.baseline_source.read_bytes()).hexdigest(),
    new_source_sha256=hashlib.sha256(Path(sys.modules['radia.kelvin_solver'].__file__).read_bytes()).hexdigest())
for maxh in (.2,.12,.09):
    mesh,source,potential,_=fixture._picard_case(maxh=maxh)
    with ng.TaskManager():
        result=solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh,source,potential,1.,(3.,0.,0.),reduced_materials=('reduced',),
            total_materials=('total',),interface_boundary='source_total_interface',order=2,
            dirichlet_bbbnd='outer',mu_r_by_material={'total':1000.},return_system=True)
        matrix=result['system']['bilinear_form'].mat;fes=result['fes']
        rows,cols,values=matrix.COO()
        assembled=coo_matrix((np.asarray(values),(np.asarray(rows),np.asarray(cols))),shape=(fes.ndof,fes.ndof)).tocsr()
        if isinstance(matrix,ng.la.SparseMatrixSymmetricdouble):
            assembled=assembled+assembled.T-diags(assembled.diagonal())
        entries=assembled.tocoo();end=fes.Range(2).start
        primal=(entries.row<end)&(entries.col<end)
        free=np.asarray(list(fes.FreeDofs()),dtype=bool)
        exact=matrix.CreateColVector();exact.FV().NumPy()[:]=np.random.default_rng(831).normal(size=fes.ndof)
        exact.FV().NumPy()[~free]=0
        rhs=matrix.CreateColVector();res=matrix.CreateColVector()
        case=dict(maxh=maxh,ndof=fes.ndof,baseline=[],updated=[])
        for name,solver in (('baseline',baseline),('updated',updated)):
            cache={}
            for scale in (.9,1.1,.7):
                mat=ng.la.SparseMatrixdouble.CreateFromCOO(entries.row,entries.col,
                    entries.data*np.where(primal,scale,1.),*assembled.shape)
                rhs.data=mat*exact
                start=time.perf_counter()
                inverse=solver(mat,fes,order=2,**({'cache':cache} if name=='updated' else {}))
                seconds=time.perf_counter()-start
                answer=inverse*rhs;res.data=mat*answer-rhs
                residual=float(np.linalg.norm(res.FV().NumPy()[free])/np.linalg.norm(rhs.FV().NumPy()[free]))
                assert residual<1e-10
                case[name].append(dict(scale=scale,seconds=seconds,residual=residual))
            if name=='updated':case['cached_primal_factor']='primal_factor' in cache
    record['cases'].append(case)
    args.output.write_text(json.dumps(record,indent=2))
    print(json.dumps(case),flush=True)
