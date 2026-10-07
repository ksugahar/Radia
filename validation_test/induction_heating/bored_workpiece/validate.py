"""Synthetic stepped axial bore: independent FEM and scalar BIE with one loop DOF.
Run from the repository root; scratch meshes are never part of the evidence.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import datetime
import importlib.metadata
from pathlib import Path
import platform
import subprocess
import sys
import time
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'src/radia/panels'))
import ngsolve as ng
import numpy as np
import scipy
import radia
import radia.bem_loop_extension as loop
from radia.bem_sibc_solver import ScalarBIESIBCSolver, compute_phi_inc_surface_poisson
from radia.biot_savart import h_segments_batch

FREQ, SIGMA, MUR, CURRENT = 100000., 5e6, 10., 100.
OMEGA = 2*math.pi*FREQ
DELTA = math.sqrt(2/(OMEGA*4e-7*math.pi*MUR*SIGMA))
ZS = (1+1j)/(SIGMA*DELTA)
FILLET = .0018
LIMITS = dict(true_residual=1e-6, faraday=1e-6, unit_jump=.005,
              bem_power_balance=.02, fem_power_balance=1e-5,
              bem_fem_power=.02, mesh_power_change=.05,
              carrier_width=.02, frozen_error_min=.20,
              delta_over_fillet=.15, elapsed_s=1200, total_elapsed_s=3600, bem_nodes=7000)


def part():
    from netgen.occ import Cylinder, Pnt, Z
    p = Cylinder(Pnt(0,0,-.025),Z,.030,.025)+Cylinder(Pnt(0,0,0),Z,.024,.025)
    p = p-Cylinder(Pnt(0,0,-.03),Z,.010,.06)-Cylinder(Pnt(0,0,0),Z,.014,.03)
    circles = [e for e in p.edges if abs(e.center.x)<1e-8 and abs(e.center.y)<1e-8]
    return p.MakeFillet(circles,FILLET)


def coils():
    theta = np.arange(129)*2*np.pi/128
    paths=[]
    for z in (-.018,0.,.018):
        p=np.c_[.040*np.cos(theta),.040*np.sin(theta),np.full(len(theta),z)]
        p[-1]=p[0]
        paths.append(list(zip(p[:-1],p[1:])))
    return paths, [CURRENT]*3


def old_point_trace(pts,tris,areas,normals,H,M_inv):
    vn=np.zeros_like(pts)
    for k in range(3): np.add.at(vn,tris[:,k],areas[:,None]*normals)
    vn/=np.linalg.norm(vn,axis=1)[:,None]
    return -np.einsum('ij,ij->i',H(pts),vn)


def compact(o):
    keys=('P_total','P_frozen','P_reaction','theta_jump','linear_residual_rel',
          'faraday_residual_rel','section_anchor','automatic_carrier')
    d={k:o[k] for k in keys}
    d.update(alpha_A=[o['alpha'].real,o['alpha'].imag],
             power_balance_rel=abs(o['P_reaction']/o['P_total']-1),
             unit_jump_error=abs(abs(o['theta_jump'])-1))
    if o['carrier_diagnostics']:
        d['carrier_diagnostics']={k:o['carrier_diagnostics'][k] for k in ('candidate_count','minimum_wall_distance','clearance_tolerance','candidates')}
    d['loop_work_diagnostics']=o['loop_work_diagnostics']
    return d


def bem(h,paths,currents,compare):
    from netgen.occ import OCCGeometry, Glue
    import netgen.meshing as nm
    from surface_mesh_extract import orient_surface_triangles
    t=time.perf_counter()
    cad=part()
    original=ng.Mesh(OCCGeometry(Glue(list(cad.faces))).GenerateMesh(
        maxh=h,curvaturesafety=1))
    pts=np.array([list(v.point) for v in original.vertices])
    tris=np.array([[v.nr for v in e.vertices] for e in original.Elements(ng.BND)])
    tris,_=orient_surface_triangles(pts,tris)
    if len(pts)>LIMITS['bem_nodes']: raise RuntimeError('dense node limit exceeded')
    surface=nm.Mesh(dim=3)
    fd=surface.Add(nm.FaceDescriptor(surfnr=1,domin=1,bc=1))
    ids=[surface.Add(nm.MeshPoint(nm.Pnt(*p))) for p in pts]
    for tri in tris: surface.Add(nm.Element2D(fd,[ids[i] for i in tri]))
    solver=ScalarBIESIBCSolver(ng.Mesh(surface),order=1,assemble_dense=True,
        use_intree_bem=True,intree_geom_order=1,intree_singular_n_q=6,
        intree_regular_quad_degree=7)
    H=sum(I*h_segments_batch(p,pts) for p,I in zip(paths,currents))
    phi,graderr=compute_phi_inc_surface_poisson(pts,tris,H,max_grad_residual=.10)
    A=lambda x: loop.A_from_filaments(x,paths,currents)
    auto=loop.solve_loop_extended(solver,phi,ZS,OMEGA,A)
    area=float(np.linalg.norm(np.cross(pts[tris[:,1]]-pts[tris[:,0]],pts[tris[:,2]]-pts[tris[:,0]]),axis=1).sum()/2)
    out=dict(nodes=len(pts),triangles=len(tris),geometry_order=1,area_m2=area,
             phi_gradient_residual=graderr,automatic=compact(auto),carriers=[])
    if compare:
        anchors=[None,(.020,-.013),(.025,-.010),(.019,.013)]
        for anchor in anchors:
            kw={}
            clearance=auto['carrier_diagnostics']['minimum_wall_distance'] if anchor is None else None
            if anchor is not None:
                r,z=anchor; th=np.arange(256)*2*np.pi/256
                ring=np.c_[r*np.cos(th),r*np.sin(th),np.full(256,z)]
                if abs(loop._surface_winding(ring[0],pts[tris]))<.99:
                    raise RuntimeError('manual carrier is outside material')
                clearance=min(loop._segment_surface_distance(x,y,pts[tris]) for x,y in zip(ring,np.roll(ring,-1,axis=0)))
                if clearance<=1e-8: raise RuntimeError('manual carrier touches surface')
                kw=dict(section_anchor=anchor,carrier_ring=ring)
            new=auto if anchor is None else loop.solve_loop_extended(solver,phi,ZS,OMEGA,A,**kw)
            with patch.object(loop,'_project_ring_neumann',old_point_trace):
                old=loop.solve_loop_extended(solver,phi,ZS,OMEGA,A,**kw)
            out['carriers'].append(dict(anchor=anchor,clearance_m=clearance,l2=compact(new),point=compact(old)))
        out['widths']={}
        for mode in ('l2','point'):
            p=np.array([r[mode]['P_total'] for r in out['carriers']])
            out['widths'][mode]=dict(power_width_W=float(np.ptp(p)),
                                    power_width_rel=float(np.ptp(p)/p.mean()))
    out['elapsed_s']=time.perf_counter()-t
    return out


def fem(h,paths,currents,work,air_radius=.18):
    from netgen.occ import OCCGeometry, Sphere, Pnt, Cylinder, Z, Glue
    from radia.em_material import EMMaterial
    import calc_fem_kelvin as K
    t=time.perf_counter()
    with ng.TaskManager():
        p=part(); p.faces.name='sibc'; p.faces.maxh=h
        ring=Cylinder(Pnt(0,0,-.06),Z,.065,.12)
        near=ring-p; near.mat('air'); near.maxh=.012
        sphere=Sphere(Pnt(0,0,0),air_radius); sphere.faces.name='outer'
        far=sphere-ring; far.mat('air')
        mesh=ng.Mesh(OCCGeometry(Glue([near,far])).GenerateMesh(maxh=.045,curvaturesafety=1.5))
        volume0=float(ng.Integrate(1,mesh))
        area0=float(ng.Integrate(1,mesh,ng.BND,definedon=mesh.Boundaries('sibc')))
        mesh.Curve(2)
        volume2=float(ng.Integrate(1,mesh))
        area2=float(ng.Integrate(1,mesh,ng.BND,definedon=mesh.Boundaries('sibc')))
        vol=work/f'air_{h:g}_{air_radius:g}.vol';mesh.ngmesh.Save(str(vol))
    from cubit_mesh_export.check import check_consistency
    contract={'schema':'radia.vol-label-contract.v1','application':'bored_workpiece',
              'required':{'materials':['air'],'boundaries':['sibc','outer']}}
    checked=check_consistency(str(vol),contract=contract)
    if not checked['passed']: raise RuntimeError('air mesh failed check-vol: '+str(checked['warnings']))
    audit={}
    # Observe the independent FEM's assembled source work, without changing it.
    # Im(conj(a)^T f)*omega/2 equals the SIBC dissipated power.
    def profile(frame,event,arg):
        if event=='return' and frame.f_code is K.solve_fem.__code__ and isinstance(arg,dict) and 'P_total' in arg:
            v=frame.f_locals
            audit['source_work_W']=float(.5*OMEGA*np.vdot(v['gfu'].vec.FV().NumPy(),v['rhs_vec'].FV().NumPy()).imag)
    def trace(frame,event,arg):
        if frame.f_code is K.solve_fem.__code__:
            frame.f_trace_lines=False
            profile(frame,event,arg)
            return trace
        return None
    previous=sys.gettrace(); sys.settrace(trace)
    try:
        with ng.TaskManager():
            r=K.solve_fem(vol_file=str(vol),fes_order=2,frequency=FREQ,
                mat=EMMaterial(name='linear alloy',sigma=SIGMA,mu_r=MUR),half_thickness=.005,
                solver='sparsecholesky',nthreads=0,filaments=(paths,currents))
    finally: sys.settrace(previous)
    if 'error' in r: raise RuntimeError(r['error'])
    out={k:r[k] for k in ('P_total','ndof','ne','linear_true_relative_residual','Z_s')}
    out.update(audit,mesh_check={k:checked[k] for k in ('schema','passed','labels','quality','boundary_domain_ownership')},
               volume_flat_m3=volume0,volume_curved_m3=volume2,geometry_order=2,area_flat_m2=area0,area_curved_m2=area2,
               air_radius_m=air_radius,elapsed_s=time.perf_counter()-t)
    out['power_balance_rel']=abs(out['source_work_W']/out['P_total']-1)
    return out


def provenance():
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    pkg=ROOT/'src/radia'
    if Path(radia.__file__).resolve().parent!=pkg.resolve(): raise RuntimeError('runtime must use checkout')
    files=[Path(__file__),pkg/'bem_loop_extension.py',pkg/'bem_sibc_solver.py',
           pkg/'cohomology.py',pkg/'bem_loop_work.py',pkg/'panels/calc_fem_kelvin.py',pkg/'biot_savart.py']
    return dict(commit=subprocess.check_output(['git','-c','safe.directory=*','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        native_sha256={p.name:sha(p) for p in sorted(pkg.iterdir()) if p.suffix in ('.pyd','.dll','.so')},
        source_dirty=bool(subprocess.check_output(['git','-c','safe.directory=*','status','--porcelain','--','src'],cwd=ROOT,text=True).strip()),
        utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        installed_native_distribution=importlib.metadata.version('radia'),
        mesh_checker_version=importlib.metadata.version('cubit-mesh-export'),
        python=platform.python_version(),ngsolve=ng.__version__,numpy=np.__version__,
        scipy=scipy.__version__,host=platform.node(),threads=2,
        radia_import='src/radia/__init__.py')


def evaluate(out):
    gates={}
    finest_b=max(out['meshes'],key=lambda r:r['bem']['nodes'])
    finest_f=max(out['meshes'],key=lambda r:r['fem']['ndof'])['fem']
    out['acceptance_reference']=dict(bem_nodes=finest_b['bem']['nodes'],fem_ndof=finest_f['ndof'])
    for i,row in enumerate(out['meshes']):
        b,f=row['bem']['automatic'],row['fem']
        gates[f'mesh_{i}']=dict(mesh_quality=f['mesh_check']['passed'],
            positive_power=min(b['P_total'],b['P_reaction'],f['P_total'],f['source_work_W'])>0,
            bem_residual=b['linear_residual_rel']<=LIMITS['true_residual'],
            faraday=b['faraday_residual_rel']<=LIMITS['faraday'],unit_jump=b['unit_jump_error']<=LIMITS['unit_jump'],
            fem_residual=f['linear_true_relative_residual']<=LIMITS['true_residual'],
            fem_balance=f['power_balance_rel']<=LIMITS['fem_power_balance'],
            frozen_error=row['frozen_fem_error_rel']>=LIMITS['frozen_error_min'],automatic=b['automatic_carrier'])
    best=finest_b['bem']['automatic']
    gates['finest_solution']=dict(
        bem_balance=best['power_balance_rel']<=LIMITS['bem_power_balance'],
        agreement=abs(best['P_total']/finest_f['P_total']-1)<=LIMITS['bem_fem_power'])
    balances=[r['bem']['automatic']['power_balance_rel'] for r in sorted(out['meshes'],key=lambda r:r['bem']['nodes'])]
    gates['balance_monotonic']=len(balances)>=3 and all(b<a for a,b in zip(balances,balances[1:]))
    if len(out['meshes'])>1:
        out['mesh_power_change']={}
        for mode in ('bem','fem'):
            last=[row[mode]['automatic']['P_total'] if mode=='bem' else row[mode]['P_total'] for row in out['meshes'][-2:]]
            out['mesh_power_change'][mode]=abs(last[0]/last[1]-1)
        gates['mesh_sensitivity']=all(v<=LIMITS['mesh_power_change'] for v in out['mesh_power_change'].values())
    widths=finest_b['bem'].get('widths')
    if widths:
        gates['carrier_solution_checks']=all(
            r[mode]['linear_residual_rel']<=LIMITS['true_residual'] and
            r[mode]['faraday_residual_rel']<=LIMITS['faraday'] and
            r[mode]['unit_jump_error']<=LIMITS['unit_jump'] and
            r[mode]['power_balance_rel']<=LIMITS['bem_power_balance']
            for r in finest_b['bem']['carriers'] for mode in ('l2','point'))
        gates['carrier_width']=widths['l2']['power_width_rel']<=LIMITS['carrier_width']
        gates['projection_reduces_width']=widths['l2']['power_width_rel']<widths['point']['power_width_rel']
    gates['runtime']=all(r['elapsed_s']<=LIMITS['elapsed_s'] for r in out['runs'])
    gates['total_runtime']=out['elapsed_s']<=LIMITS['total_elapsed_s']
    gates['skin_depth']=DELTA/FILLET<=LIMITS['delta_over_fillet']
    nodes=[r['bem']['nodes'] for r in out['meshes']]
    dofs=[r['fem']['ndof'] for r in out['meshes']]
    gates['refinement']=(len(nodes)>=3 and all(b/a>=1.25 for a,b in zip(nodes,nodes[1:]))
                         and max(dofs)/min(dofs)>=1.10 and max(nodes)<=LIMITS['bem_nodes'])
    out['gates']=gates
    def passed(v):return all(passed(x) for x in v.values()) if isinstance(v,dict) else bool(v)
    out['all_recorded_gates_pass']=passed(gates)
    out['carrier_comparison_executed']=bool(widths)
    out['complete_validation']=bool(widths) and len(out['meshes'])>=3 and gates['refinement']
    out['validation_pass']=out['all_recorded_gates_pass'] and out['complete_validation']


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--work',type=Path)
    ap.add_argument('--merge',type=Path,nargs='+',help='Combine matching split-run evidence')
    ap.add_argument('--fem-sizes',help='Independent FEM maxh values, one per BEM size')
    ap.add_argument('--out',type=Path,default=Path(__file__).with_name('results.json'))
    ap.add_argument('--sizes',default='.006,.0025,.0017')
    ap.add_argument('--no-carriers',action='store_true')
    a=ap.parse_args()
    if a.merge:
        records=[json.loads(p.read_text(encoding='utf-8')) for p in a.merge]
        first=records[0]
        for r in records:
            if (r['conditions']!=first['conditions'] or r['limits']!=LIMITS or
                any(r['provenance'][k]!=first['provenance'][k] for k in
                    ('source_sha256','native_sha256','commit','python','ngsolve','numpy','scipy'))):
                raise ValueError('split runs have incompatible conditions, limits or runtime sources')
        out=dict(first)
        out['meshes']=sorted([row for r in records for row in r['meshes']],key=lambda row:row['bem']['nodes'])
        if len({row['maxh_m'] for row in out['meshes']})!=len(out['meshes']):
            raise ValueError('duplicate BEM mesh sizes')
        out['runs']=[run for r in records for run in r['runs']]
        out['elapsed_s']=sum(run['elapsed_s'] for run in out['runs'])
        evaluate(out)
        a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        print(json.dumps(out['gates'],indent=2))
        return 0 if out['validation_pass'] else 1
    if a.work is None: ap.error('--work is required for computation')
    a.work.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter();ng.SetNumThreads(2)
    paths,currents=coils()
    out=dict(provenance=provenance(),limits=LIMITS,
        conditions=dict(frequency_Hz=FREQ,sigma_S_m=SIGMA,mu_r=MUR,current_peak_A=CURRENT,
          outer_radii_m=[.030,.024],bore_radii_m=[.010,.014],height_m=.050,step_z_m=0,
          fillet_radius_m=FILLET,skin_depth_m=DELTA,delta_over_fillet=DELTA/FILLET,
          coil_radius_m=.040,coil_z_m=[-.018,0,.018],segments_per_loop=128,
          coil_endpoint_gap_m=0,phasor='peak, exp(+i omega t)',fem_order=2,bem_order=1),meshes=[])
    def save():
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    sizes=[float(x) for x in a.sizes.split(',')]
    fem_sizes=[float(x) for x in a.fem_sizes.split(',')] if a.fem_sizes else sizes
    if len(fem_sizes)!=len(sizes): ap.error('--fem-sizes must match --sizes')
    for index,h in enumerate(sizes):
        print('MESH',h,flush=True)
        row=dict(maxh_m=h,fem_maxh_m=fem_sizes[index]);out['meshes'].append(row)
        with ng.TaskManager(): row['bem']=bem(h,paths,currents,index==len(sizes)-1 and not a.no_carriers)
        save();print('BEM',row['bem']['nodes'],row['bem']['automatic']['P_total'],flush=True)
        row['fem']=fem(fem_sizes[index],paths,currents,a.work)
        row['bem_fem_power_rel']=abs(row['bem']['automatic']['P_total']/row['fem']['P_total']-1)
        row['frozen_fem_error_rel']=abs(row['bem']['automatic']['P_frozen']/row['fem']['P_total']-1)
        save();print('FEM',row['fem']['P_total'],flush=True)
    out['elapsed_s']=time.perf_counter()-started
    out['runs']=[dict(provenance=out['provenance'],elapsed_s=out['elapsed_s'])]
    evaluate(out)
    save();print(json.dumps(out['gates'],indent=2),flush=True)
    return 0 if out['all_recorded_gates_pass'] else 1

if __name__=='__main__':
    raise SystemExit(main())
