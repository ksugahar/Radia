"""Independent small H-matrix algebra and CPU pybind sampled-Gram checks."""
from pathlib import Path
import argparse,ctypes,hashlib,importlib,json,subprocess,sys
import numpy as np
import psutil
from scipy.io import savemat
from scipy.sparse import csc_matrix
ROOT=Path(__file__).resolve().parents[2]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def structs(rows):
    out=np.empty((1,len(rows)),dtype=[(key,'O') for key in rows[0]])
    for i,row in enumerate(rows):
        for key,value in row.items():out[key][0,i]=value
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--artifact',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--aca-block-size',type=int,default=0);parser.add_argument('--full-assembly',action='store_true');parser.add_argument('--full-max-rank',type=int,default=32);parser.add_argument('--full-epsilon',type=float,default=1e-10);parser.add_argument('--physical-refinements',type=int,choices=[0,1,2],default=0);a=parser.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    record=dict(schema='radia.hmatrix-resident-numerical.v1',accepted=False,a100_used=False,timing_acceptance=False,scope='Resident symmetric leaf algebra and public sampled/integrated-tet Gram CPU parity; bounded native constant-charge block assembly only; not general HACApK Matrix Build or full HDiv-MMM physics acceptance')
    if not 1<=a.full_max_rank<=4096:raise ValueError('Full hierarchy rank budget must be 1..4096')
    if not 0<a.full_epsilon<1:raise ValueError('Full hierarchy epsilon must be in (0,1)')
    record['full_epsilon']=a.full_epsilon
    parser_rank = a.full_max_rank
    record['full_max_rank']=parser_rank
    eng=None
    source_files=(Path(__file__),Path(__file__).with_name('checkResidentHMatrix.m'),ROOT/'packages/radia-gpu/matlab/+radia/+gpu/HMatrixOperator.m')
    source_hashes={str(p.relative_to(ROOT)):sha(p) for p in source_files}
    try:
        ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(),0x4000)
        record['cpu_mean_percent']=psutil.cpu_percent(interval=2)
        if record['cpu_mean_percent']>=30 or psutil.virtual_memory().available<8*2**30:raise RuntimeError('Host budget unavailable')
        receipt=json.loads((a.artifact/'build.json').read_text(encoding='utf-8'));assert receipt['accepted'];assert sha(a.artifact/'radia_cuda_hmatrix.mexw64')==receipt['binary_sha256'];record['build']=receipt
        rng=np.random.default_rng(20261012);cases=[]
        for mode,widths in [('dense',[3,3,3]),('lowrank',[3,3,3]),('mixed',[3,3,3]),('mixed-uneven',[1,17,257])]:
            n=sum(widths);starts=np.cumsum([0]+widths);matrix=np.zeros((n,n));leaves=[]
            for i in range(3):
                for j in range(i,3):
                    ni,nj=widths[i],widths[j];ri,rj=starts[i],starts[j]
                    rank=0 if i==j or mode=='dense' or (mode.startswith('mixed') and i==0 and j==2) else min(7,ni,nj)
                    if rank:
                        a1=rng.normal(size=(nj,rank))*.03;a2=rng.normal(size=(ni,rank))*.03;block=a2@a1.T
                    else:
                        block=rng.normal(size=(ni,nj))*.01
                        if i==j:block=(block+block.T)*.5+np.eye(ni)*.2
                        a1=block.ravel(order='C').reshape(-1,1);a2=np.empty((0,1))
                    matrix[ri:ri+ni,rj:rj+nj]=block
                    if i!=j:matrix[rj:rj+nj,ri:ri+ni]=block.T
                    leaves.append(dict(row=float(ri+1),column=float(rj+1),rows=float(ni),columns=float(nj),rank=float(rank),a1=a1,a2=a2))
            permutation=rng.permutation(n);scale=rng.uniform(.5,1.5,size=n);original=np.zeros_like(matrix);original[np.ix_(permutation,permutation)]=matrix;original=scale[:,None]*original*scale[None,:]
            initial=rng.normal(size=(n,1));first=original@initial;third=original@(original@first)
            cases.append(dict(name=mode,config=dict(schema_version=1.,operator_kind='gram',dimension=float(n),permutation=(permutation+1).astype(float).reshape(-1,1),scale=scale[:,None],leaves=structs(leaves)),initial=initial,expected_first=first,expected_third=third))
        # Public CPU pybind route: regularized sampled-Laplace Gram is an explicit
        # numerical scope, never a replacement for integrated physical FE entries.
        import ngsolve as ng
        ng.SetNumThreads(1);native=importlib.import_module('radia._radia_pybind')
        points=rng.uniform(-1,1,size=(8,3));weights=rng.uniform(.1,.2,size=8)
        with ng.TaskManager():
            gram=native._ChargeGramHMatrix.from_sampled_laplace(points.ravel(),weights,.1,1e-12,64,2.,True)
            matrix=np.array([[gram.entry(i,j) for j in range(8)] for i in range(8)])
            initial=rng.normal(size=8);first=np.asarray(gram.matvec_sym(initial));third=np.asarray(gram.matvec_sym(gram.matvec_sym(first)))
        cases.append(dict(name='public-pybind-sampled-Gram',config=dict(schema_version=1.,operator_kind='gram',dimension=8.,permutation=np.arange(1,9,dtype=float)[:,None],scale=np.ones((8,1)),leaves=structs([dict(row=1.,column=1.,rows=8.,columns=8.,rank=0.,a1=matrix.ravel(order='C')[:,None],a2=np.empty((0,1)))])),initial=initial[:,None],expected_first=first[:,None],expected_third=third[:,None]))
        # Preserve the physical charge-map composition, not just the Gram apply.
        gram_case=cases[-1]
        for name,charge in [('demag-dense-map',rng.normal(size=(8,5))*.2),('demag-sparse-map',rng.normal(size=(8,5))*.2)]:
            if 'sparse' in name:charge[abs(charge)<.15]=0.
            config=dict(gram_case['config']);config.update(operator_kind='demag',charge_map=csc_matrix(charge) if 'sparse' in name else charge)
            x=rng.normal(size=5)
            with ng.TaskManager():
                def apply(v):return charge.T@np.asarray(gram.matvec_sym(charge@v))
                y=apply(x);y3=apply(apply(y))
            cases.append(dict(name=name,config=config,initial=x[:,None],expected_first=y[:,None],expected_third=y3[:,None]))
        # Integrated analytic tetrahedron operator: keep its production normalization.
        tetra=np.array([[0.,0.,0.],[.4,0.,0.],[0.,.4,0.],[0.,0.,.4]])
        tetra2=.5*tetra+np.array([.6,.1,.1]);cells=np.concatenate([tetra.ravel(),tetra2.ravel()])
        with ng.TaskManager():
            integrated=native._ChargeGramHMatrix(cells,np.empty(0),2,eps=1e-12,leaf=64,eta=2.,near_factor=1e30,build=True,far_quad=0)
            raw=np.array([[integrated.entry(i,j) for j in range(2)] for i in range(2)])
            sigma=np.asarray(integrated.charge_sigma());x=rng.normal(size=2);y=np.asarray(integrated.matvec_sym(x));y3=np.asarray(integrated.matvec_sym(integrated.matvec_sym(y)))
        normalized=raw/sigma[:,None]/sigma[None,:]
        config=dict(schema_version=1.,operator_kind='gram',dimension=2.,permutation=np.array([[1.],[2.]]),scale=sigma[:,None],leaves=structs([dict(row=1.,column=1.,rows=2.,columns=2.,rank=0.,a1=normalized.ravel(order='C')[:,None],a2=np.empty((0,1)))]))
        cases.append(dict(name='public-pybind-analytic-tet-Gram',config=config,initial=x[:,None],expected_first=y[:,None],expected_third=y3[:,None]))
        # GPU builds these entries from geometry, not CPU-generated coefficients.
        cases.append(dict(name='native-CUDA-analytic-tet-block',config=dict(schema_version=1.,operator_kind='gram',assembly='analytic_constant_block',cells=cells[:,None],faces=np.empty((0,1))),initial=x[:,None],expected_first=y[:,None],expected_third=y3[:,None]))
        triangles=np.concatenate([tetra[[0,1,2]].ravel(),tetra2[[0,2,3]].ravel()])
        with ng.TaskManager():
            mixedgram=native._ChargeGramHMatrix(cells,triangles,2,eps=1e-12,leaf=64,eta=2.,near_factor=1e30,build=True,far_quad=0)
            x4=rng.normal(size=4);first4=np.asarray(mixedgram.matvec_sym(x4));third4=np.asarray(mixedgram.matvec_sym(mixedgram.matvec_sym(first4)))
        cases.append(dict(name='native-CUDA-analytic-tet-triangle-block',config=dict(schema_version=1.,operator_kind='gram',assembly='analytic_constant_block',cells=cells[:,None],faces=triangles[:,None]),initial=x4[:,None],expected_first=first4[:,None],expected_third=third4[:,None]))
        record['CPU_oracle']=dict(module='radia._radia_pybind',binary_sha256=sha(native.__file__),ngsolve=ng.__version__,note='Independent installed public CPU reference: sampled and integrated constant-tet/triangle Gram')
        with ng.TaskManager():
            entry_rows=np.tile(np.arange(4),150);entry_columns=np.tile(np.array([3,2,1,0]),150)
            entry_expected=np.array([mixedgram.entry(int(i),int(j)) for i,j in zip(entry_rows,entry_columns)])
            aca_matrix=np.array([[mixedgram.entry(i,j) for j in range(4)] for i in range(4)])
        entry_probe=dict(cells=cells[:,None],faces=triangles[:,None],rows=(entry_rows+1.).reshape(-1,1),columns=(entry_columns+1.).reshape(-1,1),expected=entry_expected[:,None],aca_matrix=aca_matrix)
        payload=dict(cases=structs(cases),entry_probe=entry_probe)
        if a.aca_block_size:
            n=a.aca_block_size
            if not 65<=n<=256:raise ValueError('Admissible-block validation requires 65..256 charges per cluster')
            small=tetra*.25
            centers=np.array([[.3*(i%11),.3*(i//11),0.] for i in range(n)])
            larger=np.concatenate([(small+c).ravel() for c in centers]+[(small+c+np.array([30.,0.,0.])).ravel() for c in centers])
            with ng.TaskManager():
                largegram=native._ChargeGramHMatrix(larger,np.empty(0),2*n,eps=1e-10,leaf=64,eta=2.,near_factor=1e30,build=True,far_quad=0)
                block_reference=np.array([[largegram.entry(i,n+j) for j in range(n)] for i in range(n)])
                stats=largegram.stats();leaf_stats={key:stats[key] for key in ['n_dof','n_leaves','n_lowrank','n_dense','max_rank']}
                if leaf_stats['n_lowrank']<1:raise RuntimeError('CPU leaf builder did not select an admissible low-rank leaf')
            payload['aca_large_probe']=dict(cells=larger[:,None],faces=np.empty((0,1)),split=float(n),expected=block_reference,cpu_leaf_stats=leaf_stats)
            if a.full_assembly:
                with ng.TaskManager():
                    direct=np.array([[largegram.entry(i,j) for j in range(2*n)] for i in range(2*n)])
                x=rng.normal(size=(2*n,1))
                payload['full_assembly_probe']=dict(config=dict(schema_version=1.,cells=larger[:,None],faces=np.empty((0,1)),operator_kind='gram',epsilon=a.full_epsilon,leaf_size=64.,eta=2.,max_rank=float(parser_rank)),initial=x,expected=direct@x,expected_matrix=direct)
        elif a.full_assembly:raise ValueError('--full-assembly requires --aca-block-size')
        solver_cases=[]
        for test in cases:
            if not test['name'].startswith('demag-'):continue
            charge=test['config']['charge_map'];charge=charge.toarray() if hasattr(charge,'toarray') else charge
            with ng.TaskManager():
                demag=np.column_stack([charge.T@np.asarray(gram.matvec_sym(charge[:,k])) for k in range(charge.shape[1])])
            trial=rng.normal(size=(charge.shape[1],charge.shape[1]));mass=trial.T@trial+np.eye(charge.shape[1]);rhs=rng.normal(size=charge.shape[1]);inv_chi=2.5
            expected=np.linalg.solve(inv_chi*mass+demag,rhs)
            solver_cases.append(dict(config=test['config'],mass=csc_matrix(mass),rhs=rhs[:,None],inv_chi=inv_chi,expected=expected[:,None],operator=inv_chi*mass+demag))
        # Actual RT0 FE material case, using a public checked tetra fixture.
        # Python FE assembly is a validation fixture only; native MATLAB FE
        # setup remains required before this can become a GPU user workflow.
        import radia.vim._vim as cpu_vim
        meshfile=ROOT/'tests/fixtures/beam/affine_field_tetra.vol'
        mesh=ng.Mesh(str(meshfile))
        if a.physical_refinements:
            for _ in range(a.physical_refinements):mesh.Refine()
            meshfile=a.out/'physical.vol';mesh.ngmesh.Save(str(meshfile));mesh=ng.Mesh(str(meshfile))
        if mesh.dim!=3 or tuple(mesh.GetMaterials())!=('air',) or tuple(mesh.GetBoundaries())!=('boundary',):raise RuntimeError('Public tetra label contract mismatch')
        with ng.TaskManager():
            fes=ng.HDiv(mesh,order=0);cb=cpu_vim._charge_basis(fes,4)
            physical=native._ChargeGramHMatrix(cb['cell_verts'],cb['face_verts'],cb['n_el'],eps=1e-12,leaf=64,eta=2.,near_factor=1e30,build=True,far_quad=0)
            B=cb['B'];mass=cb['M_mass']
            physical.configure_charge_map(np.asarray(B.indptr,dtype=np.int32),np.asarray(B.indices,dtype=np.int32),np.asarray(B.data,dtype=np.float64),fes.ndof)
            physical.configure_mass_matrix_ngsolve(cb['M_mass_ngsolve']);physical.configure_geometry_mass_matrix_ngsolve(cb['M_mass_ngsolve'])
            load=ng.LinearForm(fes);load+=ng.CoefficientFunction((1.,.2,.1))*fes.TestFunction()*ng.dx;load.Assemble();rhs=np.array(load.vec.FV().NumPy());inv_chi=2.5
            reference=physical.solve_configured_linear_material_mass_riesz(inv_chi,rhs,1e-11,100)
            demag=np.column_stack([B.T@np.asarray(physical.matvec_sym(B[:,k].toarray().ravel())) for k in range(fes.ndof)])
            system=inv_chi*mass.toarray()+demag
            expected=np.asarray(reference['m'])
            if np.linalg.norm(rhs-system@expected)/np.linalg.norm(rhs)>1e-10:raise RuntimeError('CPU physical true residual rejected')
            means=[]
            for c in range(3):
                form=ng.LinearForm(fes);form+=fes.TestFunction()[c]*ng.dx;form.Assemble();means.append(np.array(form.vec.FV().NumPy())/ng.Integrate(1,mesh))
        config=dict(schema_version=1.,operator_kind='demag',assembly='analytic_constant_block',cells=np.asarray(cb['cell_verts'])[:,None],faces=np.asarray(cb['face_verts'])[:,None],charge_map=B.tocsc())
        for row in solver_cases:row.update(mean_operator=np.zeros((0,row['rhs'].size)),scope='Sampled-Gram material algebra')
        solver_cases.append(dict(config=config,mass=mass.tocsc(),rhs=rhs[:,None],inv_chi=inv_chi,expected=expected[:,None],operator=system,mean_operator=np.asarray(means),scope='Public affine tetra RT0 FE material solve; Python-generated FE validation fixture, not native MATLAB setup'))
        record['physical_fixture']=dict(mesh_sha256=sha(meshfile),cpu_vim_sha256=sha(cpu_vim.__file__),cpu_preconditioner='NGSolve mass-Riesz sparsecholesky',gpu_preconditioner='Explicit inv_chi mass diagonal',cpu_true_residual=float(np.linalg.norm(rhs-system@expected)/np.linalg.norm(rhs)),cpu_iterations=int(reference['iters']),native_matlab_fe_setup_accepted=False)
        payload['solver_cases']=structs(solver_cases)
        inputfile=a.out/'inputs.mat';savemat(inputfile,payload);record['input_sha256']=sha(inputfile)
        import matlab.engine
        eng=matlab.engine.start_matlab('-nodesktop -nosplash -singleCompThread');pid=int(eng.feature('getpid'));subprocess.run(['powershell','-NoProfile','-Command',f'(Get-Process -Id {pid}).PriorityClass="BelowNormal"'],check=True,capture_output=True)
        eng.addpath(str(a.artifact),str(Path(__file__).parent),str(ROOT/'packages/radia-gpu/matlab'),nargout=0)
        eng.checkResidentHMatrix(str(inputfile),str(a.out/'matlab.json'),nargout=0)
        result=json.loads((a.out/'matlab.json').read_text(encoding='utf-8'));assert result['accepted']
        if source_hashes!={str(p.relative_to(ROOT)):sha(p) for p in source_files}:raise RuntimeError('Validation source changed while running')
        record.update(accepted=True,result=result,matlab=eng.version(),source_sha256=source_hashes,memory_limitations='Exact owned device buffers; no process/global peak or idle timing measured')
    except BaseException as error:record['error']=f'{type(error).__name__}: {error}';raise
    finally:
        if eng is not None:eng.quit()
        (a.out/'result.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
