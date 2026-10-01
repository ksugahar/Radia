"""Focused numerical cases behind the result-reading application gallery."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
import numpy as np
import ngsolve as ng
import radia

ROOT = Path(__file__).resolve().parents[2]


def run(case):
    ng.SetNumThreads(2)
    if case == 'winding':
        sys.path.insert(0, str(ROOT/'docs/stream_function'))
        import demo_sphere_fe_direct as s
        mesh, fes, fi = s.build_surface(.15, .03, 3)
        obs, heldout = s.make_dsv(.08, 7), s.make_dsv(.08, 9)
        target, reference = s.target_field(obs, 'z2'), s.target_field(heldout, 'z2')
        A = s.assemble_A(fes, obs)
        psi = s.h1_min_seminorm(fes, fi, A, target)
        residual = float(np.linalg.norm(A@psi-target)/np.linalg.norm(target))
        grid, phg, thg = s.sample_to_grid(psi, mesh, .15, 96, 60)
        path, ph, th, cid, turns = s.single_stroke(grid, phg, thg, .15, 14)
        unit = s.bz_fast(path, 1., heldout)
        current = float(unit@reference/(unit@unit))
        error = float(np.linalg.norm(current*unit-reference)/np.linalg.norm(reference))
        checks = {'continuous_residual': residual < 1e-6, 'wire_relative_rms': error < .06}
        data = dict(path_m=path.tolist(), current_A=current, turns=turns,
                    continuous_relative_residual=residual, wire_relative_rms=error,
                    target='Z2', sphere_radius_m=.15, dsv_radius_m=.08,
                    checks=checks, scope='Undistorted spherical winding; no global manufacturability qualification.')
    elif case == 'lift':
        sys.path.insert(0, str(ROOT/'docs/maglev/demos/sphere'))
        import maglev_sphere_force as m
        m.main()  # Includes low/high-frequency and modal-reduction assertions.
        data=json.loads((ROOT/'validation_test/maglev/demos/sphere/maglev_sphere_force_results.json').read_text())
        data['checks']={'analytic_and_modal_assertions': True}
        data['scope']='Induced-dipole conducting sphere in a prescribed field gradient; not a closed-loop levitation simulation.'
    elif case == 'particles':
        from radia.coil_builder import saddle_coil
        from radia.particle_tracking import ParticleSpecies, TrackingBox, track_lorentz_ivp, velocity_from_kinetic_voltage
        objs=[]
        for phi,sign in ((0.,1.),(180.,-1.)):
            c=saddle_coil(sign*5000.,phi_center_deg=phi,radius=.04,length=.12,span_deg=70.,bend_radius=.006,width=.01,height=.008)
            objs.extend(c.to_radia(arc_max_segment_length=.004))
        magnet=radia.ObjCnt(objs)
        def field(x,y,z):return tuple(radia.Fld(magnet,'b',[float(x),float(y),float(z)]))
        species=ParticleSpecies(charge_c=-1.602176634e-19,mass_kg=9.1093837015e-31)
        tracks=[]
        for energy in (15.,20.,25.):
            v=velocity_from_kinetic_voltage(species,energy*1e6,[0,1,0],relativistic=True)
            tr=track_lorentz_ivp(species,[0,-.15,0],v,np.linspace(0,1.6e-9,240),magnetic_flux_density_t=field,relativistic=True,stop_box=TrackingBox(minimum_m=(-.08,-.16,-.08),maximum_m=(.08,.16,.08)),rtol=1e-9,atol=1e-12)
            tracks.append(dict(energy_MeV=energy,track=tr))
        checks={'magnetic_energy_conservation':all(t['track']['maximum_relative_kinetic_energy_drift']<1e-6 for t in tracks)}
        data=dict(tracks=tracks,checks=checks,scope='Saddle-coil electron dispersion only; quadrupole and edge-focusing studies are separate.')
    else: raise ValueError(case)
    assert all(data['checks'].values()), data['checks']
    return dict(schema='radia.selected_showcase.v1',case=case,generated_at_utc=datetime.now(timezone.utc).isoformat(),host=platform.node(),versions=dict(radia=radia.__version__,ngsolve=ng.__version__,numpy=np.__version__),driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),threads=2,result=data)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',choices=['winding','lift','particles']);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    d=run(a.case);a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item()),encoding='utf-8')
    print(a.output)
