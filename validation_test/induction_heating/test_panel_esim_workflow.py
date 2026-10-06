"""Weak-coil production handoff with a prescribed nonlinear passive cell law.

This checks orchestration, not accuracy of a ferromagnetic B-H material model.
The file loader and cell constitutive law are explicit test substitutes;
surface extraction, excitation, weighted solve, reaction and SOL I/O are real.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest
import ngsolve as ng
from netgen.occ import Cylinder, Pnt, Vec, OCCGeometry, Glue

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src/radia/panels'))
from radia.panels import calc_inductance as calc
import em_material


@pytest.mark.parametrize('backend,converges', [('intree-dense', True), ('hacapk', True), ('intree-dense', False)])
def test_panel_esim_workflow_preserves_heat_and_panel_identity(tmp_path, monkeypatch, backend, converges):
    ng.SetNumThreads(4)
    with ng.TaskManager():
        body = Cylinder(Pnt(0, 0, -.0125), Vec(0, 0, 1), r=.025, h=.025)
        for face in body.faces:
            face.name = 'sibc'
        fixture = ng.Mesh(OCCGeometry(Glue(list(body.faces))).GenerateMesh(maxh=.005))
        label = str(tmp_path/'case.vol')
        fixture.ngmesh.Save(label)
        bh = tmp_path/'linear-bh.txt'
        bh.write_text('0 0\n100 0.0001256637061435917\n1000 0.001256637061435917\n', encoding='utf-8')
        z0 = (1+1j)*np.sqrt(2*np.pi*1000*4e-7*np.pi/(2*5.8e7))

        class PrescribedCell:
            def solve(self, h):
                return {'Z': z0*(1+.2*h/(1+h)), 'converged': True}

        # The panel module currently imports this dependency by its standalone
        # script name; patch that exact module, not a second package instance.
        material = em_material
        monkeypatch.setattr(material.EMMaterial, 'create_esim_solver',
                            lambda *a, **kw: PrescribedCell())
        args = calc.build_argparser().parse_args([
            '--coil-solver', 'peec', '--coil-step', 'not-used.step', '--vol', label,
            '--wp-label', 'sibc', '--frequency', '1000', '--current', '1',
            '--sigma', '5.8e7', '--mu-r', '1', '--h1-order', '1',
            '--wp-bem-backend', backend, '--impedance-model', 'esim',
            '--bh-file', str(bh), '--half-thickness', '.003', '--esim-per-panel',
            '--esim-max-iter', '25', '--esim-tol', '.001'])
        angle = np.linspace(0, 2*np.pi, 361)
        path = np.stack([.03*np.cos(angle), .03*np.sin(angle), np.zeros_like(angle)], axis=1)
        paths = [np.stack([path[:-1], path[1:]], axis=1)]
        source = dict(source_type='filament', paths=paths, I_fil=np.array([1+0j]))
        original_mesh = ng.Mesh

        def load(value, *a, **kw):
            return fixture if isinstance(value, str) and value == label else original_mesh(value, *a, **kw)

        monkeypatch.setattr(ng, 'Mesh', load)
        if not converges:
            args.esim_max_iter = 1
            args.esim_tol = 1e-10
            with pytest.raises(RuntimeError, match='ESIM did not converge'):
                calc._solve_workpiece_weak_coupled(args, source)
            assert not list(tmp_path.glob('*qsurf.sol'))
            return
        result = calc._solve_workpiece_weak_coupled(args, source)
        assert result['esim_converged']
        assert result['esim_fixed_point_relative_error'] <= args.esim_tol
        assert result['esim_impedance_layout'] == 'BND-element-order'
        values = np.array(result['esim_per_panel_Z_s_real'])
        assert np.ptp(values) > .001*values.mean()
        assert len(values) == len(result['esim_per_panel_centroids'])
        heat = ng.GridFunction(ng.H1(fixture, order=1))
        heat.Load(result['qsurf_sol'])
        assert ng.Integrate(heat, fixture, ng.BND) == pytest.approx(result['P_wp'], rel=1e-10)
        assert np.min(heat.vec.FV().NumPy()) >= 0
        from radia import ih_thermal
        ih_thermal.verify_field_pair(result['qsurf_sol'], label, fixture, 1,
                                     quantity=ih_thermal.QSURF_QUANTITY)
        h = ng.GridFunction(ng.H1(fixture, order=1))
        h.Load(result['ht_sol'])
        assert np.all(np.isfinite(h.vec.FV().NumPy()))
        (tmp_path/'summary.json').write_text(json.dumps({
            'backend': backend, 'power_W': result['P_wp'],
            'reaction_balance': result['wp_power_balance_relative_error'],
            'panels': len(values), 'iterations': result['esim_iterations'],
            'fixed_point_error': result['esim_fixed_point_relative_error'],
        }, indent=2), encoding='utf-8')
