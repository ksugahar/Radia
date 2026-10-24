"""Separated-body Kelvin Omega / HDiv-MMM coupling investigation.

The MMM body is strictly outside the physical Kelvin ball. Its field is
therefore harmonic throughout the FEM volume. The Kelvin ball represents
only the FEM body's decaying reaction; the applied source's outer flux is
removed explicitly. No source charge is silently placed inside source-free
FEM air. This prototype uses private trace/native internals for investigation;
it is not a production coupling API.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time


def main():
    from radia.ngsolve_sparse import real_sparse_matrix_type
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('--runtime-source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--mode',choices=('parity','coupled','repeated'),default='parity')
    p.add_argument('--backend',choices=('cpu','gpu'),default='cpu')
    p.add_argument('--cpu-solver',choices=('direct','cg-local','cg-h1amg','cg-bddc-h1amg'),default='direct')
    p.add_argument('--threads',type=int,choices=range(1,17),default=1)
    p.add_argument('--radius',type=float,default=.18)
    p.add_argument('--maxh',type=float,default=.18)
    p.add_argument('--order',type=int,choices=(1,2,3),default=2)
    p.add_argument('--voxel-n',type=int,default=7)
    p.add_argument('--reuse',action='store_true',help='Reuse the CPU factor or prepared host CSR arrays; GPU transfers remain explicit')
    p.add_argument('--reuse-surface-load',action='store_true',help='Reuse the fixed connected-interface mass/constant data in this sphere probe')
    p.add_argument('--interpolate-source',action='store_true',help='Explicit P2 source interpolation with direct native-point accuracy gate')
    p.add_argument('--crosscheck-mmm',action='store_true',help='Independent two-body all-MMM solution')
    p.add_argument('--harmonic-source-degree',type=int,choices=(0,4,6,8),default=0,
                   help='Explicit local harmonic source representation, independently gated against native samples')
    p.add_argument('--lift-bonus-offset',type=int,choices=range(9),default=0)
    p.add_argument('--lift-load-route',choices=('integrated','assembled'),default='assembled')
    p.add_argument('--quadrature-audit',action='store_true')
    p.add_argument('--operator-audit',action='store_true')
    p.add_argument('--reference-source-load',choices=('volume','surface_flux'),default='volume')
    p.add_argument('--refinement-steps',type=int,choices=(0,1,2,3,4,5),default=0)
    p.add_argument('--condition-estimate',action='store_true')
    p.add_argument('--parity-trace',choices=('projected','copied'),default='projected')
    p.add_argument('--parity-flux',choices=('mean-corrected','raw'),default='mean-corrected')
    p.add_argument('--bonus-intorder',type=int,choices=(10,20,30,40),default=20)
    p.add_argument('--gpu-dll-directory',type=Path,action='append',default=[])
    p.add_argument('--native-dll-directory',type=Path,action='append',default=[])
    args=p.parse_args()
    if args.quadrature_audit and not args.operator_audit:raise ValueError('Quadrature audit requires operator audit')
    if args.lift_bonus_offset and args.lift_load_route!='integrated':raise ValueError('Lift quadrature offset requires integrated diagnostic route')
    if args.reference_source_load!='volume' and args.mode!='parity':raise ValueError('Reference source route is a parity diagnostic')
    if (args.operator_audit or args.condition_estimate or args.refinement_steps) and args.mode!='parity':raise ValueError('Operator diagnostics require parity mode')
    if args.operator_audit and not args.reuse_surface_load:raise ValueError('Operator audit requires explicit surface-load route')
    import numpy as np
    import ngsolve as ng
    if ng.__version__!='6.2.2608':raise RuntimeError('This probe requires the qualified 6.2.2608 ABI')
    from scipy.sparse import coo_matrix
    from netgen.occ import Sphere,Pnt,OCCGeometry,Glue
    from netgen.meshing import IdentificationType
    sys.path.insert(0,str(args.runtime_source))
    dlls=[os.add_dll_directory(str(x.resolve())) for x in args.native_dll_directory]
    if args.backend=='gpu':dlls.extend(os.add_dll_directory(str(x.resolve())) for x in args.gpu_dll_directory)
    elif args.gpu_dll_directory:raise ValueError('GPU DLL directories require --backend gpu')
    import radia
    from radia.vim import HDivSolver,FieldCoefficientFromSolution,FieldFromSolution,CoupledBody,SolveCoupled
    if args.backend=='gpu':
        from radia_simulink.csr import solve_spd
        import radia_simulink
    from radia.kelvin_material import MU_0,make_kelvin_mu_cf
    candidate=Path(__file__).resolve().parents[2]/'src/radia/kelvin_solver.py'
    spec=importlib.util.spec_from_file_location('radia._mixed_candidate',candidate)
    omega=importlib.util.module_from_spec(spec);spec.loader.exec_module(omega)
    ng.SetNumThreads(args.threads)
    record=dict(schema='mixed-omega.separated-coupling.v1',passed=False,backend=args.backend,
                ngsolve=ng.__version__,mode=args.mode,radius=args.radius,history=[])
    record.update(fem_backend=args.backend,mmm_backend='cpu',threads=args.threads,formulation_parity_status='OPEN')
    if args.backend=='gpu':record['csr_manifest']=radia_simulink.verify()
    record.update(bonus_intorder=args.bonus_intorder,parity_trace=args.parity_trace,parity_flux=args.parity_flux)
    if args.mode!='parity' and (args.parity_trace!='projected' or args.parity_flux!='mean-corrected'):raise ValueError('Parity route changes require parity mode')
    if args.parity_flux=='raw' and not args.reuse_surface_load:raise ValueError('Raw parity flux requires explicit surface-load route')
    record['cpu_solver']=args.cpu_solver
    record['reference_source_load']=args.reference_source_load
    record['lift_load_route']=args.lift_load_route
    record['lift_bonus_offset']=args.lift_bonus_offset
    if args.backend=='gpu' and args.cpu_solver!='direct':raise ValueError('CPU preconditioning is not a GPU fallback')
    record['source_hashes']={str(path):hashlib.sha256(path.read_bytes()).hexdigest()
                            for path in (Path(__file__),candidate,Path(radia._radia_pybind.__file__))}
    started=time.perf_counter()
    record.update(order=args.order,maxh=args.maxh,reuse=args.reuse,reuse_surface_load=args.reuse_surface_load,
                  interpolate_source=args.interpolate_source,harmonic_source_degree=args.harmonic_source_degree)
    if args.interpolate_source and args.harmonic_source_degree:raise ValueError('Select one explicit source approximation')
    assert .85-args.radius>.5 and .05<args.radius<.3
    args.out.parent.mkdir(parents=True,exist_ok=True)
    try:
        iron=Sphere(Pnt(0,0,0),.2);iron.mat('iron');iron.faces.name='interface';iron.maxh=args.maxh/2
        outer=Sphere(Pnt(0,0,0),.5);outer.faces.name='kelvin_int'
        air=outer-iron;air.mat('air')
        exterior=Sphere(Pnt(1.5,0,0),.5);exterior.mat('kelvin');exterior.faces.name='kelvin_ext'
        for vertex in exterior.vertices:
            if vertex.p[2]>0:vertex.name='GND'
        geometry=Glue([iron,air,exterior])
        faces={f.name:f for s in geometry.solids for f in s.faces if f.name}
        faces['kelvin_int'].Identify(faces['kelvin_ext'],'periodic',IdentificationType.PERIODIC)
        with ng.TaskManager():
            mesh=ng.Mesh(OCCGeometry(geometry).GenerateMesh(maxh=args.maxh,grading=.6));mesh.Curve(args.order)
            fes=ng.Periodic(ng.H1(mesh,order=args.order,dirichlet_bbbnd='GND',
                               **({'wb_withedges':False} if args.cpu_solver=='cg-bddc-h1amg' else {})))
            if mesh.ne>50000 or fes.ndof>200000:raise RuntimeError('Declared mesh/DOF validation envelope exceeded')
            mu=make_kelvin_mu_cf(mesh,.5,(1.5,0,0),mu_r_by_material={'iron':100.})
            u,v=fes.TnT();a=ng.BilinearForm(fes,symmetric=True);a+=mu*ng.grad(u)*ng.grad(v)*ng.dx(bonus_intorder=args.bonus_intorder)
            preconditioner=None
            if args.cpu_solver.startswith('cg-'):
                if args.cpu_solver=='cg-bddc-h1amg':
                    preconditioner=ng.Preconditioner(a,'bddc',coarsetype='h1amg',inverse='sparsecholesky')
                else:
                    preconditioner=ng.Preconditioner(a,args.cpu_solver[3:],inverse='sparsecholesky')
            a.Assemble()
            iron_matrix=None
            if args.lift_load_route=='assembled' or args.operator_audit:
                iron_matrix=ng.BilinearForm(fes,symmetric=True)
                iron_matrix+=mu*ng.grad(u)*ng.grad(v)*ng.dx(definedon=mesh.Materials('iron'),bonus_intorder=args.bonus_intorder)
                iron_matrix.Assemble()
        mask=np.array(list(fes.FreeDofs()),dtype=bool)
        n=int(mask.sum())
        matrix=None
        # CPU Krylov uses NGSolve's assembled operator directly. Exporting and
        # rebuilding a second sparse matrix wastes setup time and host memory.
        if preconditioner is None:
            ii,jj,vv=a.mat.COO();full=coo_matrix((vv,(ii,jj)),shape=(fes.ndof,fes.ndof)).tocsr()
            compact=full[mask][:,mask].tocoo()
            if n>np.iinfo(np.int32).max or compact.nnz>np.iinfo(np.int32).max:raise OverflowError('Native CSR int32 limit exceeded')
            matrix=real_sparse_matrix_type().CreateFromCOO(compact.row.tolist(),compact.col.tolist(),compact.data.tolist(),n,n)
        record['operator_route']='native-ngsolve' if matrix is None else 'compact-export'
        record['mesh_and_matrix_seconds']=time.perf_counter()-started
        prepare_started=time.perf_counter()
        cached_inverse=None
        prepared=None
        if args.reuse:
            if args.backend=='cpu' and args.cpu_solver=='direct':
                with ng.TaskManager():cached_inverse=matrix.Inverse(inverse='sparsecholesky')
            elif args.backend=='gpu':
                from radia_simulink import _csr_cuda
                record['gpu_binding_sha256']=hashlib.sha256(Path(_csr_cuda.__file__).read_bytes()).hexdigest()
                vals,cols,rows=matrix.CSR()
                prepared=(np.ascontiguousarray(rows,dtype=np.int32),np.ascontiguousarray(cols,dtype=np.int32),
                          np.ascontiguousarray(vals),np.ascontiguousarray(1/full[mask][:,mask].diagonal()))
        record['solver_prepare_seconds']=time.perf_counter()-prepare_started
        normal=ng.specialcf.normal(3)
        source_sign=omega.material_set_boundary_signs(mesh,('air',))['interface']
        interface_measure=ng.ds(definedon=mesh.Boundaries('interface'),bonus_intorder=args.bonus_intorder)
        if args.reuse_surface_load:
            with ng.TaskManager():
                boundary_mass=ng.LinearForm(fes);boundary_mass+=v.Trace()*interface_measure;boundary_mass.Assemble()
                constant=ng.GridFunction(fes);constant.Set(1.)
                constant_values=np.array(constant.vec.FV().NumPy(),copy=True)
                mass_values=np.array(boundary_mass.vec.FV().NumPy(),copy=True)
                interface_area=float(np.dot(constant_values,mass_values))
                if not np.isfinite(interface_area) or interface_area<=0:raise RuntimeError('Invalid connected interface area')
        timings=[];record['solves']=timings
        def fem(source):
            with ng.TaskManager():
                phase=time.perf_counter()
                projection=omega.project_source_interface_potential(mesh,source,'interface',order=args.order,relative_tolerance=.01)
                lift=ng.GridFunction(fes)
                omega._copy_matching_h1_trace_to_global(mesh,projection['potential'],lift,'interface',order=args.order)
                projection_seconds=time.perf_counter()-phase
                phase=time.perf_counter()
                load=ng.LinearForm(fes)
                source_values=None
                if args.reuse_surface_load:
                    flux=MU_0*source_sign*ng.InnerProduct(source,normal)
                    raw=ng.LinearForm(fes);raw+=flux*v.Trace()*interface_measure;raw.Assemble()
                    absolute=ng.LinearForm(fes);absolute+=ng.IfPos(flux,flux,-flux)*v.Trace()*interface_measure;absolute.Assemble()
                    source_values=np.array(raw.vec.FV().NumPy(),copy=True)
                    net=float(np.dot(constant_values,source_values))
                    norm=float(np.dot(constant_values,absolute.vec.FV().NumPy()))
                    ratio=abs(net)/max(norm,1e-300)
                    if not np.isfinite(ratio) or ratio>1e-6:raise RuntimeError('Source flux compatibility gate failed')
                    if args.parity_flux=='mean-corrected':source_values-=net/interface_area*mass_values
                    balance=dict(relative_flux_imbalance=ratio)
                else:
                    boundary_load,balance=omega.outward_source_flux_form(
                        mesh,source,{'interface':source_sign},v.Trace(),materials=('air',),
                        order=args.order,bonus_intorder=args.bonus_intorder,scale=MU_0)
                    load+=boundary_load
                if args.lift_load_route=='integrated':load+=-mu*ng.grad(lift)*ng.grad(v)*ng.dx(definedon=mesh.Materials('iron'),bonus_intorder=args.bonus_intorder+args.lift_bonus_offset)
                load.Assemble()
                if args.lift_load_route=='assembled':load.vec.data-=iron_matrix.mat*lift.vec
                if source_values is not None:load.vec.FV().NumPy()[:]+=source_values
                rhs=load.vec.CreateVector() if matrix is None else matrix.CreateColVector()
                rhs.FV().NumPy()[:]=load.vec.FV().NumPy() if matrix is None else load.vec.FV().NumPy()[mask]
                rhs_seconds=time.perf_counter()-phase
                t=time.perf_counter()
                if preconditioner is not None:
                    iteration=ng.CGSolver(a.mat,preconditioner.mat,tol=1e-12,maxiter=10000,printrates=False)
                    solved=load.vec.CreateVector();solved.data=iteration*load.vec
                    values=np.array(solved.FV().NumPy()[mask],copy=True)
                    info=dict(backend='cpu',solver=args.cpu_solver,iterations=iteration.iterations,preconditioner_reused=True)
                elif cached_inverse is not None:
                    solved=rhs.CreateVector();solved.data=cached_inverse*rhs
                    values=np.array(solved.FV().NumPy(),copy=True)
                    info=dict(backend='cpu',solver='NGSolve sparsecholesky',factor_reused=True)
                elif prepared is not None:
                    row,col,vals,invd=prepared
                    values,info=_csr_cuda.solve(row,col,vals,np.ascontiguousarray(rhs.FV().NumPy()),invd,1e-12,10000)
                elif args.backend=='cpu':
                    inverse=matrix.Inverse(inverse='sparsecholesky')
                    direct=rhs.CreateVector();direct.data=inverse*rhs
                    values=np.array(direct.FV().NumPy(),copy=True)
                    info=dict(backend='cpu',solver='NGSolve sparsecholesky',factor_reused=False)
                else:
                    values,info=solve_spd(matrix,rhs,backend='gpu',tolerance=1e-12,max_iterations=10000)
                # Check the original assembled free equations for every route,
                # independently of any compact CSR representation.
                original_solution=load.vec.CreateVector();original_solution[:]=0.
                original_solution.FV().NumPy()[mask]=values
                residue=load.vec.CreateVector();residue.data=load.vec-a.mat*original_solution
                true_r=float(np.linalg.norm(residue.FV().NumPy()[mask])/np.linalg.norm(load.vec.FV().NumPy()[mask]))
                if not np.isfinite(true_r) or true_r>1e-6:raise RuntimeError('Original FEM residual failed')
                info=dict(info,true_relative_residual=true_r)
                timings.append(dict(seconds=time.perf_counter()-t,info=info,
                                    projection_seconds=projection_seconds,rhs_assembly_seconds=rhs_seconds,
                                    source_relative_flux_imbalance=balance['relative_flux_imbalance']))
                state=ng.GridFunction(fes);state.vec.FV().NumPy()[mask]=values
                ht=-ng.grad(state)-ng.grad(lift)
                mean=np.array([ng.Integrate(ht[k],mesh,definedon=mesh.Materials('iron')) for k in range(3)])/ng.Integrate(1,mesh,definedon=mesh.Materials('iron'))
            return state,lift,projection,mean
        applied=ng.CF((0,0,1))
        state,lift,projection,mean=fem(applied)
        record['uncoupled_analytic_mean_relative_error']=float(np.linalg.norm(mean-[0,0,3/102])/(3/102))
        if args.mode=='repeated':
            reference=mean.copy()
            for _ in range(3):state,lift,projection,mean=fem(applied)
            record.update(mean_H=mean.tolist(),repeat_mean_difference=float(np.linalg.norm(reference-mean)),
                          analytic_mean_relative_error=float(np.linalg.norm(mean-[0,0,3/102])/(3/102)),
                          passed=bool(np.linalg.norm(reference-mean)<1e-8 and np.linalg.norm(mean-[0,0,3/102])/(3/102)<.01))
            points=[(0,0,0),(.1,0,0),(0,.1,.1),(.3,.1,.1),(-.3,.05,.1),(1.6,.1,.2),(1.5,0,.4)]
            actual=[];exact=[]
            for xyz in points:
                point=mesh(*xyz);position=np.array(xyz)
                if xyz[0]>1:
                    actual.append(np.array(ng.grad(state)(point)))
                    local=position-[1.5,0,0];rho=np.linalg.norm(local);unit=local/rho
                    physical=.25*local/rho**2;distance=np.linalg.norm(physical)
                    reaction=99/102*.2**3*(3*physical*physical[2]/distance**5-np.array([0.,0.,1.])/distance**3)
                    exact.append(-(.5/rho)**2*(reaction-2*unit*np.dot(unit,reaction)))
                elif np.linalg.norm(position)<.2:
                    actual.append(np.array((-ng.grad(state)-ng.grad(lift))(point)))
                    exact.append(np.array([0,0,3/102]))
                else:
                    actual.append(np.array((applied-ng.grad(state))(point)))
                    distance=np.linalg.norm(position)
                    exact.append(np.array([0.,0.,1.])+99/102*.2**3*(3*position*position[2]/distance**5-np.array([0.,0.,1.])/distance**3))
            error=float(np.linalg.norm(np.array(actual)-exact)/np.linalg.norm(exact))
            record.update(probe_points=points,probe_H=np.array(actual).tolist(),probe_exact_H=np.array(exact).tolist(),probe_relative_rms=error)
            record['passed']=bool(record['passed'] and error<.01)
        elif args.mode=='parity':
            with ng.TaskManager():
                trace_difference=projection['potential']-lift
                record['interface_trace_l2_difference']=float(np.sqrt(ng.Integrate(trace_difference**2*interface_measure,mesh)))
                reference=omega.solve_magnetostatic_mixed_total_reduced_omega_kelvin(
                    mesh,applied,lift if args.parity_trace=='copied' else projection['potential'],.5,(1.5,0,0),mu_cf=mu,
                    reduced_materials=('air',),total_materials=('iron','kelvin'),
                    interface_boundary='interface',kelvin_interface_boundary='kelvin_int',
                    kelvin_source_potential=ng.CF(0),order=args.order,bonus_intorder=args.bonus_intorder,return_system=True,
                    reduced_normal_flux=MU_0*normal[2],reduced_flux_boundary='kelvin_int',
                    reduced_source_load=args.reference_source_load)
            with ng.TaskManager():
                system=reference['system'];saddle=system['bilinear_form'];saddle_load=system['linear_form']
                saddle_state=reference['solution'];saddle_space=saddle_state.space
                saddle_free=np.array(list(saddle_space.FreeDofs()),dtype=bool)
                record['reference_factorization']='matching-trace constraint elimination; sparsecholesky only on trace mass and constrained primal SPD matrices; no indefinite Cholesky'
                record['reference_original_residual']=reference['linear_residual']
                record['refinement_steps']=args.refinement_steps
                inverse=omega._matching_trace_direct_inverse(saddle.mat,saddle_space,order=args.order)
                def saddle_residual():
                    residual=saddle_load.vec.CreateVector();residual.data=saddle_load.vec-saddle.mat*saddle_state.vec
                    return residual
                record['saddle_refinement_residuals']=[]
                for step in range(args.refinement_steps+1):
                    residual=saddle_residual()
                    relative=float(np.linalg.norm(residual.FV().NumPy()[saddle_free])/np.linalg.norm(saddle_load.vec.FV().NumPy()[saddle_free]))
                    record['saddle_refinement_residuals'].append(relative)
                    if step<args.refinement_steps:saddle_state.vec.data+=inverse*residual
                record['reference_refined_block_residual']=omega._report_direct_solve_residual(saddle,saddle_load,saddle_state,saddle_space,block_names=('phi_reduced','phi_total','interface_constraint'))
                if args.operator_audit:
                    from scipy.sparse import coo_matrix,diags
                    from scipy.sparse.linalg import norm as sparse_norm
                    def csr(original):
                        from radia.ngsolve_sparse import real_sparse_matrix_type
                        rows,cols,vals=original.COO()
                        result=coo_matrix((np.asarray(vals),(np.asarray(rows),np.asarray(cols))),shape=(original.height,original.width)).tocsr()
                        if isinstance(original,real_sparse_matrix_type(symmetric=True)):result=result+result.T-diags(result.diagonal())
                        return result
                    locations={};ranges=[saddle_space.Range(i).start for i in (0,1)]
                    entities=((ng.VERTEX,mesh.vertices),(ng.EDGE,mesh.edges),(ng.FACE,mesh.faces))
                    for kind,group in entities:
                        for entity in group:
                            node=ng.NodeId(kind,entity.nr)
                            global_ids=[x for x in fes.GetDofNrs(node) if x>=0]
                            for component in (0,1):
                                local_ids=[x for x in saddle_space.components[component].GetDofNrs(node) if x>=0]
                                if not local_ids:continue
                                if len(local_ids)!=len(global_ids):raise RuntimeError('Nonmatching audit entity modes')
                                for local_id,global_id in zip(local_ids,global_ids):
                                    row=ranges[component]+local_id
                                    if not saddle_free[row]:continue
                                    if row in locations and locations[row]!=global_id:raise RuntimeError('Ambiguous audit mapping')
                                    locations[row]=global_id
                    primal_end=saddle_space.Range(2).start
                    needed=set(np.flatnonzero(saddle_free[:primal_end]))
                    if set(locations)!=needed:raise RuntimeError('Unmapped active audit primal modes')
                    rr=np.array(list(locations));cc=np.array([locations[x] for x in rr])
                    transform=coo_matrix((np.ones(len(rr)),(rr,cc)),shape=(saddle_space.ndof,fes.ndof)).tocsr()
                    k=csr(saddle.mat);ag=csr(a.mat)
                    effective=(transform.T@k@transform).tocsr()
                    lifted=np.zeros(saddle_space.ndof)
                    total_rows=rr[rr>=ranges[1]];total_cols=np.array([locations[x] for x in total_rows])
                    lifted[total_rows]=lift.vec.FV().NumPy()[total_cols]
                    target=np.asarray(saddle_load.vec.FV().NumPy())
                    condensed_rhs=np.asarray(transform.T@(target-k@lifted))
                    # Reassemble the declared SPD load on this very same mesh.
                    own=ng.LinearForm(fes)
                    own+=MU_0*source_sign*ng.InnerProduct(applied,normal)*v.Trace()*interface_measure
                    own+=-mu*ng.grad(lift)*ng.grad(v)*ng.dx(definedon=mesh.Materials('iron'),bonus_intorder=args.bonus_intorder+args.lift_bonus_offset)
                    own.Assemble();own_rhs=np.array(own.vec.FV().NumPy(),copy=True)
                    if args.parity_flux=='mean-corrected':
                        flux_vector=ng.LinearForm(fes);flux_vector+=MU_0*source_sign*ng.InnerProduct(applied,normal)*v.Trace()*interface_measure;flux_vector.Assemble()
                        net=float(np.dot(constant_values,flux_vector.vec.FV().NumPy()))
                        own_rhs-=net/interface_area*mass_values
                    consistent=own.vec.CreateVector();consistent.data=-iron_matrix.mat*lift.vec
                    flux_only=ng.LinearForm(fes);flux_only+=MU_0*source_sign*ng.InnerProduct(applied,normal)*v.Trace()*interface_measure;flux_only.Assemble()
                    consistent.data+=flux_only.vec
                    consistent_rhs=np.array(consistent.FV().NumPy(),copy=True)
                    if args.parity_flux=='mean-corrected':consistent_rhs-=net/interface_area*mass_values
                    record['lift_load_audit']=dict(integrated_vs_assembled_l2=float(np.linalg.norm((own_rhs-consistent_rhs)[mask])),assembled_vs_reference_l2=float(np.linalg.norm((condensed_rhs-consistent_rhs)[mask])))
                    if args.quadrature_audit:
                        record['shared_quadrature_lift_audit']=[]
                        for rule_order in (12,20,30):
                            measure=ng.dx(definedon=mesh.Materials('iron'),intrules={ng.ET.TET:ng.IntegrationRule(ng.ET.TET,rule_order)})
                            operator=ng.BilinearForm(fes,symmetric=True);operator+=mu*ng.grad(u)*ng.grad(v)*measure;operator.Assemble()
                            form=ng.LinearForm(fes);form+=mu*ng.grad(lift)*ng.grad(v)*measure;form.Assemble()
                            action=form.vec.CreateVector();action.data=operator.mat*lift.vec
                            delta=form.vec.FV().NumPy()-action.FV().NumPy()
                            record['shared_quadrature_lift_audit'].append(dict(rule_order=rule_order,absolute_difference=float(np.linalg.norm(delta[mask])),relative_difference=float(np.linalg.norm(delta[mask])/np.linalg.norm(action.FV().NumPy()[mask]))))

                    paired_state=ng.GridFunction(fes)
                    paired_rhs=own.vec.CreateVector();paired_rhs.FV().NumPy()[:]=consistent_rhs
                    paired_inverse=a.mat.Inverse(fes.FreeDofs(),inverse='sparsecholesky')
                    paired_state.vec.data=paired_inverse*paired_rhs
                    paired_points=[(0,0,0),(.1,0,0),(.3,.1,.1),(1.6,.1,.2)]
                    def sample(field_state):
                        result=[]
                        for xyz in paired_points:
                            point=mesh(*xyz)
                            cf=ng.grad(field_state) if xyz[0]>1 else (-ng.grad(field_state)-ng.grad(lift) if np.linalg.norm(xyz)<.2 else applied-ng.grad(field_state))
                            result.append(np.array(cf(point)).tolist())
                        return result
                    paired_reference=[list(reference['H_cf'](mesh(*xyz))) for xyz in paired_points]
                    paired_actual=sample(paired_state);paired_existing=sample(state)
                    paired_max=float(max(np.linalg.norm(np.array(x)-y) for x,y in zip(paired_actual,paired_reference)))
                    record['same_mesh_lift_comparison']=dict(mesh_scope='single in-memory mesh and identical source lift for every compared state',points=paired_points,existing_H=paired_existing,assembled_lift_H=paired_actual,reference_H=paired_reference,assembled_max_absolute_difference=paired_max,passed=paired_max<1e-8,gate=1e-8)
                    if args.lift_load_route=='assembled':own_rhs=consistent_rhs
                    delta_matrix=(effective-ag)[mask][:,mask]
                    defect=k@lifted-target
                    lm_rows=np.arange(primal_end,saddle_space.ndof)[saddle_free[primal_end:]]
                    record['operator_audit']=dict(
                        matrix_relative_frobenius=float(sparse_norm(delta_matrix)/sparse_norm(ag[mask][:,mask])),
                        rhs_relative_difference=float(np.linalg.norm((condensed_rhs-own_rhs)[mask])/np.linalg.norm(own_rhs[mask])),
                        rhs_absolute_difference=float(np.linalg.norm((condensed_rhs-own_rhs)[mask])),
                        copied_lift_constraint_relative_defect=float(np.linalg.norm(defect[lm_rows])/np.linalg.norm(target[lm_rows])),
                        projected_state_reduced_residual=float(np.linalg.norm((condensed_rhs-effective@state.vec.FV().NumPy())[mask])/np.linalg.norm(condensed_rhs[mask])),
                        reference_trace_jump_l2=float(np.sqrt(abs(ng.Integrate((reference['phi_total']-reference['phi_reduced']-projection['potential'])**2*interface_measure,mesh)))))

                if args.condition_estimate:
                    from scipy.sparse.linalg import LinearOperator,onenormest
                    def estimate(original,free,apply_inverse):
                        count=int(free.sum())
                        def action(values,invert=False):
                            rhs=original.CreateColVector();rhs[:]=0.;rhs.FV().NumPy()[free]=values
                            answer=rhs.CreateVector();answer.data=apply_inverse*rhs if invert else original*rhs
                            return np.array(answer.FV().NumPy()[free],copy=True)
                        op=LinearOperator((count,count),matvec=action,rmatvec=action,dtype=float)
                        inv=LinearOperator((count,count),matvec=lambda x:action(x,True),rmatvec=lambda x:action(x,True),dtype=float)
                        return float(onenormest(op)*onenormest(inv))
                    record['saddle_condition_1_estimate']=estimate(saddle.mat,saddle_free,inverse)
                    primal_inverse=a.mat.Inverse(fes.FreeDofs(),inverse='sparsecholesky')
                    record['spd_condition_1_estimate']=estimate(a.mat,mask,primal_inverse)
                    record['condition_estimate_scope']='lower estimates from scipy onenormest and implemented inverse; not certified condition bounds'
            points=[(0,0,0),(.1,0,0),(.3,.1,.1),(1.6,.1,.2)]
            errors=[]
            for xyz in points:
                point=mesh(*xyz)
                if xyz[0]>1:
                    field=np.array(ng.grad(state)(point))
                elif np.linalg.norm(xyz)<.2:
                    field=np.array((-ng.grad(state)-ng.grad(lift))(point))
                else:field=np.array((applied-ng.grad(state))(point))
                errors.append(float(np.linalg.norm(field-reference['H_cf'](point))))
            record['evaluation_points']=[list(x) for x in points]
            record['evaluation_metric']='max Euclidean absolute H vector difference in A/m'
            record.update(max_absolute_field_difference=max(errors),field_errors=errors,mean_H=mean.tolist(),
                          analytic_mean_H=[0,0,3/102],passed=max(errors)<1e-8)
        else:
            if record['uncoupled_analytic_mean_relative_error']>.01:
                raise RuntimeError('Uncoupled Kelvin discretization does not meet the frozen 1% accuracy gate')
            coupled_started=time.perf_counter()
            body=Sphere(Pnt(0,0,.85),args.radius);body.mat('iron')
            with ng.TaskManager():
                body_mesh=ng.Mesh(OCCGeometry(body).GenerateMesh(maxh=args.radius))
                body_solver=HDivSolver(body_mesh,order=1)
                result=body_solver.Solve(mu_r=100.,H_ext=applied,linear_solver='cpp-cg',preconditioner='mass-riesz',tol=1e-11)
            coordinates=np.asarray([vertex.point for vertex in body_mesh.vertices],dtype=np.float64)
            cells=np.asarray([[vertex.nr for vertex in element.vertices] for element in body_mesh.Elements(ng.VOL)],dtype=np.int64)
            record['mmm_mesh_sha256']=hashlib.sha256(coordinates.tobytes()+cells.tobytes()).hexdigest()
            record['initial_mmm_M_avg']=np.asarray(result['M_avg']).tolist()
            def reaction_at(state,physical):
                rho=np.linalg.norm(physical);unit=physical/rho
                mapped=np.array([1.5,0,0])+.25*physical/rho**2
                hc=np.array(ng.grad(state)(mesh(*mapped)))
                return -(.5/rho)**2*(hc-2*unit*np.dot(unit,hc))
            def reaction_voxel(state):
                lo=np.array([-args.radius,-args.radius,.85-args.radius])*1.
                hi=np.array([args.radius,args.radius,.85+args.radius])*1.
                grids=[np.linspace(lo[k],hi[k],args.voxel_n) for k in range(3)]
                values=np.zeros((args.voxel_n,args.voxel_n,args.voxel_n,3))
                for iz,z in enumerate(grids[2]):
                    for iy,y in enumerate(grids[1]):
                        for ix,x in enumerate(grids[0]):
                            values[iz,iy,ix]=reaction_at(state,np.array([x,y,z]))
                coefficient=ng.CF(tuple(ng.VoxelCoefficient(tuple(lo),tuple(hi),np.ascontiguousarray(values[:,:,:,k]),linear=True) for k in range(3)))
                directions=np.array([[x,y,z] for x in (-1,0,1) for y in (-1,0,1) for z in (-1,0,1) if (x,y,z)!=(0,0,0)],dtype=float)
                points=np.array([0,0,.85])+.7*args.radius*directions/np.linalg.norm(directions,axis=1)[:,None]
                exact=np.array([reaction_at(state,point) for point in points])
                approximate=np.array([coefficient(body_mesh(*point)) for point in points])
                error=float(np.linalg.norm(approximate-exact)/np.sqrt(len(points)))
                if error>1e-4:raise RuntimeError(f'Kelvin reaction transfer RMS {error} exceeds 1e-4 A/m')
                record.setdefault('reaction_interpolation_absolute_rms',[]).append(error)
                return coefficient
            def native_source(result):
                with ng.TaskManager():field=FieldCoefficientFromSolution(result,algorithm='direct')
                if args.harmonic_source_degree:
                    from scipy.linalg import null_space
                    degree=args.harmonic_source_degree;scale=.24
                    exponents=[(i,j,k) for i in range(degree+1) for j in range(degree+1-i) for k in range(degree+1-i-j)]
                    low=[e for e in exponents if sum(e)<=degree-2]
                    laplace=np.zeros((len(low),len(exponents)));lookup={e:i for i,e in enumerate(low)}
                    for col,e in enumerate(exponents):
                        for axis in range(3):
                            if e[axis]>=2:
                                reduced=list(e);reduced[axis]-=2
                                laplace[lookup[tuple(reduced)],col]+=e[axis]*(e[axis]-1)
                    harmonic=null_space(laplace)
                    count=2*(degree+1)**2+12
                    index=np.arange(count);zz=1-2*(index+.5)/count;angle=index*np.pi*(3-np.sqrt(5))
                    points=scale*np.column_stack((np.sqrt(1-zz*zz)*np.cos(angle),np.sqrt(1-zz*zz)*np.sin(angle),zz))
                    def design(points):
                        xyz=points/scale;blocks=[]
                        for axis in range(3):
                            block=np.zeros((len(points),len(exponents)))
                            for col,e in enumerate(exponents):
                                if e[axis]:
                                    reduced=list(e);reduced[axis]-=1
                                    block[:,col]=-e[axis]/scale*np.prod(xyz**np.array(reduced),axis=1)
                            blocks.append(block)
                        return np.stack(blocks,axis=1).reshape(3*len(points),-1)
                    phase=time.perf_counter()
                    exact=FieldFromSolution(result,points,algorithm='direct')
                    weights=np.linalg.lstsq(design(points)@harmonic,exact.reshape(-1),rcond=None)[0]
                    coefficients=harmonic@weights
                    xyz=(ng.x/scale,ng.y/scale,ng.z/scale);components=[]
                    for axis in range(3):
                        terms=[]
                        for coefficient,e in zip(coefficients,exponents):
                            if e[axis]:
                                reduced=list(e);reduced[axis]-=1
                                term=-float(coefficient)*e[axis]/scale
                                for k in range(3):term=term*xyz[k]**reduced[k]
                                terms.append(term)
                        components.append(sum(terms))
                    source=ng.CF(tuple(components)).Compile()
                    directions=np.array([[x,y,z] for x in (-1,0,1) for y in (-1,0,1) for z in (-1,0,1) if (x,y,z)!=(0,0,0)],dtype=float)
                    held=.22*directions/np.linalg.norm(directions,axis=1)[:,None]
                    held_exact=FieldFromSolution(result,held,algorithm='direct')
                    approximation=(design(held)@coefficients).reshape(-1,3)
                    error=float(np.linalg.norm(approximation-held_exact)/np.sqrt(len(held)))
                    if error>1e-5:raise RuntimeError(f'Harmonic source holdout RMS {error} exceeds 1e-5 A/m')
                    record.setdefault('harmonic_source',[]).append(dict(degree=degree,sample_count=count,
                        holdout_count=len(held),holdout_absolute_rms=error,laplace_coefficient_residual=float(np.linalg.norm(laplace@coefficients)),seconds=time.perf_counter()-phase))
                    return source
                if args.interpolate_source:
                    with ng.TaskManager():
                        data=omega.interpolate_source_field(mesh,field,mesh.Materials('air'),order=2,
                            point_evaluator=lambda points:FieldFromSolution(result,points,algorithm='direct'))
                    directions=np.array([[x,y,z] for x in (-1,0,1) for y in (-1,0,1) for z in (-1,0,1) if (x,y,z)!=(0,0,0)],dtype=float)
                    points=.22*directions/np.linalg.norm(directions,axis=1)[:,None]
                    exact=FieldFromSolution(result,points,algorithm='direct')
                    approximate=np.array([data['field'](mesh(*point)) for point in points])
                    error=float(np.linalg.norm(approximate-exact)/np.sqrt(len(points)))
                    if error>1e-4:raise RuntimeError(f'Interpolated source absolute RMS {error} exceeds 1e-4 A/m')
                    record.setdefault('source_interpolation',[]).append(dict(nodes=data['nodes'],absolute_rms=error,timings=data['timings_seconds']))
                    return data['field']
                return field
            def original_mmm_residual(result,external):
                space=result['gfM'].space
                with ng.TaskManager():
                    uu,vv=space.TnT();mass=ng.BilinearForm(space);mass+=uu*vv*ng.dx;mass.Assemble()
                    load=ng.LinearForm(space);load+=external*vv*ng.dx;load.Assemble()
                    vector=mass.mat.CreateColVector();vector.FV().NumPy()[:]=result['_m_coefficients']
                    product=mass.mat.CreateColVector();product.data=mass.mat*vector
                    ax=np.array(product.FV().NumPy())/99+np.array(result['_charge_gram'].apply_configured_demag(np.ascontiguousarray(result['_m_coefficients']),True))
                    rhs=np.array(load.vec.FV().NumPy())
                return float(np.linalg.norm(ax-rhs)/np.linalg.norm(rhs))
            for iteration in range(1,31):
                iteration_started=time.perf_counter()
                old=np.array(result['_m_coefficients'],copy=True)
                field=native_source(result)
                state,lift,projection,mean=fem(applied+field)
                external=applied+reaction_voxel(state)
                with ng.TaskManager():
                    result=body_solver.Solve(mu_r=100.,H_ext=external,linear_solver='cpp-cg',preconditioner='mass-riesz',tol=1e-11)
                step=float(np.linalg.norm(result['_m_coefficients']-old)/np.linalg.norm(result['_m_coefficients']))
                record['history'].append(dict(iteration=iteration,relative_step=step,mean_H=mean.tolist(),seconds=time.perf_counter()-iteration_started))
                print(json.dumps(record['history'][-1]),flush=True)
                if step<1e-8:
                    field=native_source(result)
                    state,lift,projection,mean=fem(applied+field)
                    true_residual=original_mmm_residual(result,applied+reaction_voxel(state))
                    if true_residual<1e-6:break
            else:raise RuntimeError('Coupled iteration failed to converge')
            record.update(passed=True,coupled_original_mmm_residual=true_residual,
                          mean_H=mean.tolist(),mmm_M_avg=np.array(result['M_avg']).tolist(),
                          voxel_nodes=args.voxel_n,iterations=iteration,mmm_operator_builds=body_solver.operator_build_count)
            record['coupled_seconds']=time.perf_counter()-coupled_started
            record['coupled_only_passed']=True
            if args.crosscheck_mmm:
                record.update(passed=False,stage='independent_all_mmm')
                args.out.write_text(json.dumps(record,indent=2)+'\n')
                reference_started=time.perf_counter()
                reference_body=Sphere(Pnt(0,0,0),.2);reference_body.mat('iron')
                with ng.TaskManager():
                    reference_mesh=ng.Mesh(OCCGeometry(reference_body).GenerateMesh(maxh=.14))
                    options=dict(linear_solver='cpp-cg',preconditioner='mass-riesz',tol=1e-11)
                    independent=SolveCoupled([
                        CoupledBody(reference_mesh,'reference-a',mu_r=100.,solve_options=options),
                        CoupledBody(body_mesh,'reference-b',mu_r=100.,solve_options=options)],H_ext=applied,tol=1e-8,maxit=30)
                    a_ref,b_ref=independent['bodies']
                    a_field=np.array(a_ref['M_avg'])/99
                    err_a=float(np.linalg.norm(mean-a_field)/np.linalg.norm(a_field))
                    err_b=float(np.linalg.norm(np.array(result['M_avg'])-b_ref['M_avg'])/np.linalg.norm(b_ref['M_avg']))
                    residual_a=original_mmm_residual(a_ref,applied+FieldCoefficientFromSolution(b_ref,algorithm='direct'))
                    residual_b=original_mmm_residual(b_ref,applied+FieldCoefficientFromSolution(a_ref,algorithm='direct'))
                record['independent_all_mmm']=dict(relative_field_a=err_a,relative_magnetization_b=err_b,
                    iterations=independent['iterations'],original_residual_a=residual_a,original_residual_b=residual_b)
                record['independent_reference_seconds']=time.perf_counter()-reference_started
                if max(err_a,err_b)>.01 or max(residual_a,residual_b)>1e-6:
                    record['passed']=False
                    raise RuntimeError('Independent coupled MMM accuracy/residual gate failed')
                record.update(passed=True,stage='complete')
        record.update(fem_elements=mesh.ne,free_dofs=n,solves=timings,elapsed_seconds=time.perf_counter()-started)
        assert record['passed'],record
    except Exception as error:
        record['failure']=repr(error);raise
    finally:
        record['elapsed_seconds']=time.perf_counter()-started
        args.out.write_text(json.dumps(record,indent=2)+'\n')
        for dll in dlls:dll.close()
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
