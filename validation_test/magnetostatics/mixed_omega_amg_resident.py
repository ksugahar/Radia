"""Validation-only CPU hierarchy setup and device-resident AMG execution."""
import hashlib
import importlib.util
import json
from pathlib import Path
import time
import numpy as np

class ResidentAMG:
    def __init__(self,row,col,values,*,directory,candidate,strength=0.,damping=.8,cache_path=None,cache_mode='build',updated_transfers=None,nonlinear_cache_path=None,nonlinear_cache_mode='build'):
        from radia.ngsolve_sparse import real_sparse_matrix_type
        if nonlinear_cache_mode=='reuse' and nonlinear_cache_path is None:raise ValueError('Nonlinear transfer cache path required')
        if nonlinear_cache_mode not in ('build','reuse'):raise ValueError('Explicit nonlinear AMG cache mode required')
        if nonlinear_cache_path is not None and (cache_path is not None or updated_transfers is not None):raise ValueError('Distinct nonlinear transfer cache required')
        if nonlinear_cache_path is not None and nonlinear_cache_mode=='build' and Path(nonlinear_cache_path).exists():raise FileExistsError('Refuse to overwrite nonlinear AMG cache')
        if cache_mode not in {'build','reuse','reuse-hierarchy'} or (cache_mode!='build' and cache_path is None):raise ValueError('Explicit AMG cache mode/path required')
        if updated_transfers is not None and (cache_path is not None or cache_mode!='build'):
            raise ValueError('Updated-operator transfer reuse is distinct from unchanged-operator cache reuse')
        if cache_mode=='build' and updated_transfers is None and not (nonlinear_cache_path is not None and nonlinear_cache_mode=='reuse'):
            import pyamg
            pyamg_version=pyamg.__version__
        else:
            from importlib.metadata import version
            pyamg_version=version('pyamg')
        import ngsolve as ng
        from scipy.sparse import csr_matrix
        if not 0<=strength<=1 or not 0<damping<2:raise ValueError('Invalid AMG strength/damping')
        if cache_path is not None and cache_mode=='build' and Path(cache_path).exists():raise FileExistsError('Refuse to overwrite an AMG cache')
        started=time.perf_counter()
        binaries=list(Path(directory).glob('_resident_amg*.pyd'))
        if len(binaries)!=1:raise RuntimeError('Exactly one resident AMG module is required')
        binary=binaries[0];binary_hash=hashlib.sha256(binary.read_bytes()).hexdigest()
        build_path=Path(directory)/'resident-amg-build.json'
        build_bytes=build_path.read_bytes();build=json.loads(build_bytes)
        if build.get('binary_sha256')!=binary_hash:raise RuntimeError('Resident AMG build receipt hash mismatch')
        spec=importlib.util.spec_from_file_location('_resident_amg',binary)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        from mixed_omega_amg_cache import nonlinear_layout
        layout=nonlinear_layout(row,col)
        a=csr_matrix((values,col,row)).copy();a.sum_duplicates();a.eliminate_zeros();a.sort_indices()
        difference=a-a.T
        if difference.nnz and np.max(np.abs(difference.data))>1e-12*np.max(np.abs(a.data)):raise RuntimeError('AMG requires symmetric operator')
        b=np.ascontiguousarray(candidate,dtype=np.float64).reshape(-1,1)
        if b.shape[0]!=a.shape[0] or not np.isfinite(b).all() or not np.linalg.norm(b):raise ValueError('Invalid near-nullspace candidate')
        from mixed_omega_amg_cache import load_transfers,save_transfers,symmetric_row_bound
        parameters=dict(ngsolve=ng.__version__,pyamg=pyamg_version,strength=strength,damping=damping,
                        adapter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        cache_receipt=dict(reused=False)
        if nonlinear_cache_path is not None and nonlinear_cache_mode=='reuse':
            from mixed_omega_amg_cache import load_nonlinear_transfers
            updated_transfers,cache_receipt=load_nonlinear_transfers(nonlinear_cache_path,b,parameters=parameters,layout=layout)
        if updated_transfers is not None:
            from types import SimpleNamespace
            hierarchy=[];current=a
            for transfer in updated_transfers:
                transfer=transfer.tocsr(copy=False)
                if transfer.shape[0]!=current.shape[0] or not np.isfinite(transfer.data).all():raise ValueError('Updated AMG transfer mismatch')
                restriction=transfer.T.tocsr()
                hierarchy.append(SimpleNamespace(A=current,P=transfer,R=restriction))
                current=(restriction@current@transfer).tocsr();current.sum_duplicates();current.sort_indices()
            hierarchy.append(SimpleNamespace(A=current))
            cache_receipt.update(updated_operator=True,transfers_reused=True,galerkin_operators_rebuilt=True,fine_operator_is_current=True)
        elif cache_mode in {'reuse','reuse-hierarchy'}:
            hierarchy,cache_receipt=load_transfers(cache_path,a,b,parameters=parameters,rebuild_galerkin=cache_mode=='reuse')
        else:
            rng=np.random.get_state()
            try:
                np.random.seed(1729)
                hierarchy=pyamg.smoothed_aggregation_solver(a,B=b,symmetry='symmetric',strength=('symmetric',{'theta':strength}),max_coarse=32,max_levels=12).levels
            finally:np.random.set_state(rng)
        self.transfers=tuple(level.P for level in hierarchy[:-1])
        setup_seconds=time.perf_counter()-started
        levels=[]
        def export(m):
            m=m.tocsr();m.sum_duplicates();m.sort_indices()
            if m.nnz>np.iinfo(np.int32).max or max(m.shape)>np.iinfo(np.int32).max:raise OverflowError('AMG CSR int32 limit')
            return dict(n=m.shape[0],m=m.shape[1],row=np.ascontiguousarray(m.indptr,dtype=np.int32),col=np.ascontiguousarray(m.indices,dtype=np.int32),val=np.ascontiguousarray(m.data,dtype=np.float64))
        sizes=[]
        for i,level in enumerate(hierarchy):
            diagonal=level.A.diagonal()
            if np.any(diagonal<=0) or not np.isfinite(diagonal).all():raise RuntimeError('AMG diagonal must be positive')
            # Certified upper bound on normalized operator spectral radius.
            bound=symmetric_row_bound(level.A)
            item=dict(a=export(level.A),smoother=np.ascontiguousarray(damping/bound/diagonal))
            if i+1<len(hierarchy):
                trans=level.R-level.P.T
                if trans.nnz and np.max(abs(trans.data))>1e-13:raise RuntimeError('AMG restriction must be transpose prolongation')
                item.update(p=export(level.P),r=export(level.R))
            levels.append(item);sizes.append(dict(n=level.A.shape[0],nnz=level.A.nnz,smoother_bound=bound))
        coarse=hierarchy[-1].A.tocoo();n=coarse.shape[0]
        if n>128:raise RuntimeError('AMG coarse action exceeds bounded dense export')
        # This validation adapter calls NGSolve only for coarse preparation;
        # it owns that CPU region, separately from the caller's FE assembly.
        with ng.TaskManager():
            matrix=real_sparse_matrix_type().CreateFromCOO(coarse.row.tolist(),coarse.col.tolist(),coarse.data.tolist(),n,n)
            inverse=matrix.Inverse(inverse='sparsecholesky');rhs=matrix.CreateColVector();answer=matrix.CreateColVector();dense=np.empty((n,n))
            for i in range(n):
                rhs[:]=0.;rhs[i]=1.;answer.data=inverse*rhs;dense[:,i]=answer.FV().NumPy()
        defect=float(np.linalg.norm(coarse.tocsr()@dense-np.eye(n),ord=np.inf))
        if not np.isfinite(defect) or defect>1e-8:raise RuntimeError('AMG coarse inverse action failed')
        dense=(dense+dense.T)/2
        if np.linalg.eigvalsh(dense)[0]<=0:raise RuntimeError('AMG coarse action is not SPD')
        digest=hashlib.sha256()
        for level in levels:
            for name in ('a','p','r'):
                if name in level:
                    for component in ('row','col','val'):digest.update(level[name][component].view(np.uint8))
            digest.update(level['smoother'].view(np.uint8))
        digest.update(np.ascontiguousarray(dense).view(np.uint8))
        before=time.perf_counter();self.native=module.Resident(levels,np.ascontiguousarray(dense.ravel()))
        device_upload_seconds=time.perf_counter()-before
        if cache_path is not None and cache_mode=='build':
            saved=time.perf_counter();cache_receipt.update(saved=True,payload_sha256=save_transfers(cache_path,hierarchy,b,parameters=parameters),save_seconds=time.perf_counter()-saved)
        if nonlinear_cache_path is not None and nonlinear_cache_mode=='build':
            from mixed_omega_amg_cache import save_nonlinear_transfers
            saved=time.perf_counter();cache_receipt.update(saved=True,transfers_only=True,payload_sha256=save_nonlinear_transfers(nonlinear_cache_path,self.transfers,b,parameters=parameters,layout=layout),save_seconds=time.perf_counter()-saved)
        self.operator = a
        self.receipt=dict(hierarchy_seed=1729,hierarchy_sha256=digest.hexdigest(),adapter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),strength_threshold=strength,smoother_damping=damping,matlab_engine_loaded='matlab.engine' in __import__('sys').modules,pyamg_version=pyamg_version,pyamg_imported='pyamg' in __import__('sys').modules,binary_sha256=binary_hash,build_receipt_sha256=hashlib.sha256(build_bytes).hexdigest(),hierarchy_levels=sizes,hierarchy_setup_seconds=setup_seconds,coarse_inverse_defect=defect,coarse_route='NGSolve sparsecholesky action exported once',device_upload_seconds=device_upload_seconds,prepare_seconds=time.perf_counter()-started,near_nullspace='NGSolve H1 constant coefficients restricted to free DOFs',preconditioner='symmetric Galerkin V-cycle, two pre/post damped-Jacobi sweeps, fixed SPD coarse action',cpu_hierarchy_setup=True,cpu_preconditioner_execution=False,host_orchestration='graph launch and convergence status per eight iterations',engine_runtime_required=False,cache=cache_receipt)
    def solve(self, rhs):
        rhs = np.ascontiguousarray(rhs, dtype=np.float64)
        solution, info = self.native.solve(rhs, 1e-12, 10000)
        started = time.perf_counter()
        defect = self.operator @ solution - rhs
        denominator = float(np.linalg.norm(rhs))
        residual = float(np.linalg.norm(defect)) / denominator if denominator else float(np.linalg.norm(defect))
        if not np.isfinite(residual) or residual > 1e-6:
            raise RuntimeError(f'Independent host original residual gate: {residual}')
        info = dict(info)
        info.update(host_true_relative_residual=residual, host_certificate_gate=1e-6,
                    host_residual_seconds=time.perf_counter()-started)
        return solution, info
    def close(self):self.native.close();return dict(released=True)
