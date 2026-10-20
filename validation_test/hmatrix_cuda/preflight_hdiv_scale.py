"""Public-mesh FE capacity probe, without material maps, Gram matrices or GPU work."""
from pathlib import Path
import argparse, contextlib, hashlib, json, runpy


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--refinements',type=int,choices=[0,1],default=0)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    repo=Path(__file__).resolve().parents[2];result=dict(accepted=False,scope='CPU FE metadata and current software capacity only; no solve, GPU allocation or timing prediction')
    def check(mesh,label):
        report=a.out/(label+'-check.json')
        # Reuse the already initialized one-thread NGSolve runtime. A nested
        # interpreter would initialize a second default-size local heap.
        checker=runpy.run_path(str(repo/'packages/cubit-mesh-export/src/cubit_mesh_export/check.py'))
        args=[str(mesh),'--strict-labels','--contract',str(repo/'validation_test/hmatrix_cuda/curved_projection_validation.contract.json'),'--tet-only','--report-json',str(report)]
        with (a.out/(label+'-check.log')).open('w') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
            status=checker['main'](args)
        if status:raise RuntimeError('Canonical mesh check failed: '+label)
        data=json.loads(report.read_text());assert data['passed'];return data
    try:
        import ngsolve as ng
        from ngsolve.comp import IntegrationRuleSpace
        ng.SetNumThreads(1);check(a.mesh,'source');source_hash=hashlib.sha256(a.mesh.read_bytes()).hexdigest()
        mesh=ng.Mesh(str(a.mesh));target=a.mesh
        if a.refinements:
            with ng.TaskManager():mesh.Refine()
            target=a.out/'refined.vol';mesh.ngmesh.Save(str(target))
        checked=check(target,'target');mesh=ng.Mesh(str(target));before=hashlib.sha256(target.read_bytes()).hexdigest()
        with ng.TaskManager():
            fes=ng.HDiv(mesh,order=2);samples=IntegrationRuleSpace(mesh,order=6)
            ne=mesh.ne;nf=sum(1 for _ in mesh.Elements(ng.BND));nd=fes.ndof
            assert 0<ne<=16384 and 0<nd<=300000
            shape_entries=0;tangent=0;chunk_points=[];chunk=0
            for index,element in enumerate(mesh.Elements(ng.VOL)):
                local=len(fes.GetDofNrs(element));points=len(samples.GetDofNrs(element))
                shape_entries+=3*local*points;tangent+=local*local;chunk+=points
                if (index+1)%128==0:chunk_points.append(chunk);chunk=0
            if chunk:chunk_points.append(chunk)
            tangent=min(tangent,nd*nd);charges=4*ne+6*nf
            material=24*shape_entries+288*tangent+384*max(chunk_points)+64*2**20
            resident=16*charges*charges+64*2**20;solver=128*tangent+512*nd+64*2**20
        assert hashlib.sha256(target.read_bytes()).hexdigest()==before
        native=repo/'packages/radia-simulink/native';full=(native/'hmatrix_full_assembly.cuh').read_text();elements=(native/'hdiv_element_setup.cuh').read_text()
        assert all(x in full for x in ['n>4608','(std::size_t(384)<<20)']) and all(x in elements for x in ['ne>4096','nd>16384','chargeCount>8192'])
        gates={'native_elements':ne<=4096,'native_dofs':nd<=16384,'fe_charge_count':charges<=8192,'resident_charge_count':charges<=4608,'hacapk_signed_product':charges*charges<=2147483647//2,'hacapk_conservative_allocation':16*charges*charges<=384*2**20,'legacy_dense_mass_guard':8*nd*nd<=2*2**30,'legacy_dense_charge_guard':8*nd*charges<=2*2**30,'material_points':samples.ndof<=2000000,'material_array_budget':material<=2**30,'material_solve_combined_array_budget':material+resident+solver<=2*2**30}
        result.update(accepted=True,ngsolve_version=ng.__version__,threads=1,refinements=a.refinements,source_mesh_sha256=source_hash,mesh_sha256=before,curve_order=checked['mesh']['curve_order'],elements=ne,boundary_elements=nf,dofs=nd,charge_count=charges,material_points=samples.ndof,potential_material_shape_entries=shape_entries,material_tangent_nonzero_bound=tangent,largest_chunk_points=max(chunk_points),dense_upper_triangle_coefficient_bytes=8*charges*(charges+1)//2,current_material_array_bound_bytes=material,current_resident_worst_case_bound_bytes=resident,current_solver_bound_bytes=solver,current_gates=gates,current_pipeline_eligible=all(gates.values()),failed_current_gates=[k for k,v in gates.items() if not v],native_source_sha256={str(q.relative_to(repo)).replace('\\','/'):hashlib.sha256(q.read_bytes()).hexdigest() for q in [native/'hmatrix_full_assembly.cuh',native/'hdiv_element_setup.cuh']},limitations=['Bounds do not predict actual compressed storage or CUDA/host/global peak','Geometry cache capacity and work remain unestimated','Refined-mesh metadata does not certify physical FE error','Above-2GiB payloads would require HDF5 if actually exchanged; none are allocated here','No numerical acceptance or 100k-DoF performance claim'])
    except BaseException as ex:result['error']=repr(ex);raise
    finally:(a.out/'result.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
