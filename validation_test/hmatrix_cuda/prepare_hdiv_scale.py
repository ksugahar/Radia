"""Independent public sphere fixtures and bounded memory estimates only."""
from pathlib import Path
import argparse,hashlib,json,re


def estimate(dofs,cells,faces,order,tet_points,tri_points):
    if min(dofs,cells,faces,tet_points,tri_points)<1 or order not in (1,2):raise ValueError('Invalid scale dimensions')
    charges=cells*(1 if order==1 else 4)+faces*(3 if order==1 else 6)
    mass=8*dofs*dofs;mapping=8*charges*dofs;gram=8*charges*charges
    points=cells*tet_points+faces*tri_points
    cache=8*charges*points
    tile_sources=min(charges,max(1,(256*2**20)//(8*points)));resident_cache=8*tile_sources*points
    return dict(dofs=dofs,cells=cells,faces=faces,order=order,charges=charges,mass_dense_bytes=mass,charge_map_dense_bytes=mapping,
                gram_dense_bytes=gram,potential_cache_bytes=cache,estimated_resident_potential_cache_bytes=resident_cache,cache_source_tile_count=tile_sources,host_staging_two_copies_bytes=2*(mass+mapping),
                factor_bytes=mass,device_lower_bound_bytes=gram+resident_cache+mass,
                limitations='Lower bounds exclude cuSOLVER workspace, sparse maps, allocator/runtime, geometry and host Python/MATLAB temporaries',
                native_dimension_gate=dofs<=16384 and cells<=4096 and charges<=8192,
                mass_cholesky_array_gate=mass<=2*2**30,
                dense_gram_array_gate=gram<=2*2**30,
                dense_array_gate=max(mass,mapping)<=2*2**30 and 2*(mass+mapping)<=6*2**30,
                cache_workspace_lower_bound_gate=resident_cache<=2*2**30)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--maxh',type=float,nargs='+',default=[.6,.45,.35]);p.add_argument('--quadrature-file',type=Path);p.add_argument('--sphere-separation',type=float,help='Two independently authored unit spheres; centre distance in metres, strictly greater than 2');a=p.parse_args()
    import math
    if any(not math.isfinite(h) or h<=0 for h in a.maxh):p.error('maxh must be finite and positive')
    if a.sphere_separation is not None and (not math.isfinite(a.sphere_separation) or a.sphere_separation<=2):p.error('sphere-separation must be finite and greater than 2 metres')
    a.out.mkdir(exist_ok=False,parents=True)
    import ngsolve as ng
    from netgen.csg import CSGeometry,Sphere,Pnt
    ng.SetNumThreads(1)
    source=a.quadrature_file or Path(__file__).resolve().parents[2]/'packages/radia-simulink/matlab/+radia/+gpu/curvedHDivQuadrature.m'
    section=source.read_text().split('elseif order == 2')[1]
    def count(name):return len(re.findall(r';',section.split('q.'+name+' = [',1)[1].split('];',1)[0]))
    nt,ns=count('tet_rule'),count('tri_rule')
    results=[]
    for maxh in a.maxh:
        geometry=CSGeometry()
        centres=[0.] if a.sphere_separation is None else [-a.sphere_separation/2,a.sphere_separation/2]
        for x in centres:geometry.Add(Sphere(Pnt(x,0,0),1).mat('magnetic_body').bc('exterior'))
        mesh=ng.Mesh(geometry.GenerateMesh(maxh=maxh));mesh.Curve(2)
        prefix='sphere-' if a.sphere_separation is None else 'two-spheres-'
        path=a.out/(prefix+str(maxh).replace('.','p')+'.vol');mesh.ngmesh.Save(str(path));loaded=ng.Mesh(str(path))
        with ng.TaskManager():fes=ng.HDiv(loaded,order=2)
        record=estimate(fes.ndof,loaded.ne,len(list(loaded.Elements(ng.BND))),2,nt,ns)
        record.update(mesh=path.name,mesh_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),maxh=maxh,curve_order=loaded.GetCurveOrder(),materials=loaded.GetMaterials(),boundaries=loaded.GetBoundaries(),ngsolve=ng.__version__,units='SI; unit radius in metres',source='Independently authored elementary CSG sphere geometry',sphere_centres_m=[[x,0,0] for x in centres],sphere_radius_m=1.,hacapk_curved_charge_gate=record['charges']<=4608,full_potential_cache_2gib_lower_bound_gate=record['potential_cache_bytes']<=2*2**30)
        results.append(record);(a.out/'scale-estimates.json').write_text(json.dumps(dict(accepted=True,scope='Input generation and memory lower bounds only; no numerical acceptance',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),quadrature_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cases=results),indent=2))
if __name__=='__main__':main()
