"""Exercise the real workpiece stage and saved heat, using in-memory fixtures.

Only the file-load boundary is substituted; production extraction, incident
projection, BIE/loop solve, heat projection, reciprocity and writers execute.
This is not a certification of coil CAD, VOL export, or native Simulink MEX.
"""
import sys
import json
import hashlib
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
import ngsolve as ng
from netgen.occ import Cylinder, Pnt, Vec, OCCGeometry, Glue

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'/'radia'/'panels'))
import calc_inductance as calc


@pytest.mark.parametrize('bore,backend', [(0., 'intree-dense'), (0., 'hacapk'), (.0075, 'intree-dense')])
def test_production_workpiece_heat_and_reaction(tmp_path, bore, backend):
    ng.SetNumThreads(4)
    with ng.TaskManager():
        body = Cylinder(Pnt(0,0,-.0125), Vec(0,0,1), r=.025, h=.025)
        if bore:
            body = body - Cylinder(Pnt(0,0,-.013), Vec(0,0,1), r=bore, h=.026)
        for face in body.faces:
            face.name = 'sibc'
        fixture = ng.Mesh(OCCGeometry(Glue(list(body.faces))).GenerateMesh(maxh=.003))
        label = str(tmp_path/'in_memory.vol')
        args = calc.build_argparser().parse_args([
            '--coil-solver','peec','--coil-step','not-used.step','--vol',label,
            '--wp-label','sibc','--frequency','1000','--current','1',
            '--sigma','5.8e7','--mu-r','1','--h1-order','1',
            '--wp-bem-backend',backend, '--msh-output',str(tmp_path/'heat.msh')])
        angle = np.linspace(0,2*np.pi,721)
        path = np.stack([.03*np.cos(angle), .03*np.sin(angle), np.zeros_like(angle)], axis=1)
        paths = [np.stack([path[:-1],path[1:]],axis=1)]
        source = dict(source_type='filament',paths=paths,I_fil=np.array([1.+0j]))
        original_mesh = ng.Mesh
        def load(value, *a, **kw):
            return fixture if isinstance(value,str) and value == label else original_mesh(value,*a,**kw)
        with patch.object(ng, 'Mesh', side_effect=load):
            result = calc._solve_workpiece_weak_coupled(args,source)
        assert result['wp_genus'] == int(bool(bore))
        assert result['wp_power_balance_relative_error'] < .02
        assert result['qsurf_method'] == 'solved-total-field-lumped-P1'
        saved = ng.GridFunction(ng.H1(fixture,order=1))
        saved.Load(result['qsurf_sol'])
        assert np.min(saved.vec.FV().NumPy()) >= 0
        power = ng.Integrate(saved,fixture,ng.BND)
        assert power == pytest.approx(result['P_wp'],rel=1e-10)
        mesh_text = (tmp_path/'heat.msh').read_text()
        assert '$MeshFormat\n4.1' in mesh_text.replace('\r\n','\n')
        assert 'q_surf_W_per_m2' in mesh_text
        assert result['delta_R_mOhm'] > 0
        result['calc_source_sha256'] = hashlib.sha256(Path(calc.__file__).read_bytes()).hexdigest()
        result['saved_heat_power_W'] = float(power)
        (tmp_path/'result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
        # Unsupported higher-order EM postprocessing must not reach assembly
        # or silently restore the old incident-field heat approximation.
        args.h1_order = 2
        with patch.object(ng, 'Mesh', side_effect=load):
            with pytest.raises(ValueError, match='P1'):
                calc._solve_workpiece_weak_coupled(args,source)
