"""Carrier sensitivity of the weak trace on a smooth torus and a creased tube.

The old point trace is injected only into this controlled historical
comparison. Production calls below always use the unmodified solver.
Run with QT_EVIDENCE_DIR to retain machine-readable numerical evidence.
"""
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
import pytest
import ngsolve as ng
import radia.bem_loop_extension as loop
from test_loop_extension_ring import ring_setup, Z_S, OMEGA, MU_0


def surface_solver(pts, tris):
    import netgen.meshing as nm
    from radia.bem_sibc_solver import ScalarBIESIBCSolver
    mesh = nm.Mesh(dim=3)
    fd = mesh.Add(nm.FaceDescriptor(surfnr=1, domin=1, bc=1))
    ids = [mesh.Add(nm.MeshPoint(nm.Pnt(*p))) for p in pts]
    for t in tris:
        mesh.Add(nm.Element2D(fd, [ids[i] for i in t]))
    return ScalarBIESIBCSolver(ng.Mesh(mesh), order=1, assemble_dense=True,
        use_intree_bem=True, intree_geom_order=1, intree_singular_n_q=6,
        intree_regular_quad_degree=7)


@pytest.fixture(scope="module")
def stepped_setup():
    # Piecewise-linear meridional profile, genuinely sharp circular creases.
    corners = np.array([[.024,-.004],[.036,-.004],[.036,0.],
                        [.033,0.],[.033,.004],[.024,.004]])
    profile = []
    for a,b in zip(corners, np.roll(corners,-1,axis=0)):
        n = int(np.ceil(np.linalg.norm(b-a)/.0015))
        profile.extend(a+(b-a)*j/n for j in range(n))
    count, m = 128, len(profile)
    pts = np.array([[r*np.cos(t),r*np.sin(t),z]
        for t in np.arange(count)*2*np.pi/count for r,z in profile])
    tris = []
    for i in range(count):
        for j in range(m):
            a=i*m+j; b=((i+1)%count)*m+j
            c=((i+1)%count)*m+(j+1)%m; d=i*m+(j+1)%m
            tris.extend([(a,b,c),(a,c,d)])
    with ng.TaskManager():
        solver = surface_solver(pts, np.array(tris))
    def A_inc(x):
        return .5*MU_0*np.c_[-x[:,1],x[:,0],np.zeros(len(x))]
    return solver, -pts[:,2].astype(complex), A_inc


def old_point_trace(pts, tris, areas, normals, H, M_inv, **kwargs):
    vnorm = np.zeros_like(pts)
    for k in range(3):
        np.add.at(vnorm, tris[:,k], areas[:,None]*normals)
    vnorm /= np.linalg.norm(vnorm,axis=1)[:,None]
    return -np.einsum('ij,ij->i', H(pts), vnorm)


def compact(out):
    return {k: out[k] for k in ('P_total','P_reaction','theta_jump',
        'linear_residual_rel','faraday_residual_rel')} | {
        'alpha_real':out['alpha'].real, 'alpha_imag':out['alpha'].imag,
        'power_balance_rel':abs(out['P_reaction']/out['P_total']-1)}


@pytest.mark.parametrize('geometry', ['torus','stepped'])
def test_carrier_sensitivity(geometry, request, monkeypatch):
    solver, phi, A_inc = request.getfixturevalue(
        'ring_setup' if geometry=='torus' else 'stepped_setup')
    anchors = ([(.030,0.),(.029,0.),(.031,0.),(.030,.001)] if geometry=='torus'
        else [(.029,-.002),(.027,-.002),(.033,-.002),(.029,.002)])
    project = loop._project_ring_neumann
    rows = []
    started = time.perf_counter()
    for anchor in anchors:
        theta = np.arange(256)*2*np.pi/256
        r,z = anchor
        ring = np.c_[r*np.cos(theta),r*np.sin(theta),np.full(256,z)]
        row = {'anchor':anchor}
        for mode in ('old','qt6','qt12'):
            with monkeypatch.context() as patch:
                if mode=='old':
                    patch.setattr(loop, '_project_ring_neumann', old_point_trace)
                elif mode=='qt12':
                    patch.setattr(loop, '_project_ring_neumann',
                        lambda *a: project(*a, quadrature_n=12))
                with ng.TaskManager():
                    out = loop.solve_loop_extended(solver,phi,Z_S,OMEGA,A_inc,
                        section_anchor=anchor,carrier_ring=ring)
            row[mode] = compact(out)
        rows.append(row)
    widths = {}
    for mode in ('old','qt6','qt12'):
        p = np.array([r[mode]['P_total'] for r in rows])
        alpha = np.array([complex(r[mode]['alpha_real'],r[mode]['alpha_imag']) for r in rows])
        widths[mode] = {'P_width_W':float(np.ptp(p)),
            'P_width_percent':float(100*np.ptp(p)/p.mean()),
            'alpha_width_A':float(np.max(np.abs(alpha[:,None]-alpha)))}
    quadrature_error = max(abs(r['qt6']['P_total']/r['qt12']['P_total']-1) for r in rows)
    evidence = dict(geometry=geometry, nv=solver.mesh.nv, rows=rows,
        widths=widths, quadrature_power_relative_error=quadrature_error,
        elapsed_s=time.perf_counter()-started,host=platform.node(),
        ngsolve=ng.__version__, production_module=Path(loop.__file__).resolve().relative_to(
            Path(__file__).resolve().parents[2]).as_posix())
    if os.getenv('QT_EVIDENCE_DIR'):
        path = Path(os.environ['QT_EVIDENCE_DIR']); path.mkdir(parents=True,exist_ok=True)
        (path/f'qt_{geometry}.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
    for row in rows:
        for mode in ('old','qt6','qt12'):
            out = row[mode]
            assert out['linear_residual_rel'] < 1e-6
            assert out['faraday_residual_rel'] < 1e-6
            assert abs(abs(out['theta_jump'])-1) < 5e-3
            # The existing coarse analytic-ring validation uses 10%.
            assert out['P_total']>0 and out['P_reaction']>0
            assert out['power_balance_rel'] < .1
    assert quadrature_error < 1e-4
    # A weak trace need not reduce total discretization error on a smooth
    # mesh: Theta interpolation and the unchanged Faraday row also contribute.
    # The crease-dominated case, not the torus, tests the reduction hypothesis.
    if geometry == 'stepped':
        assert widths['qt6']['P_width_percent'] < widths['old']['P_width_percent']
    # Residual carrier dependence is finite; no claim of exact invariance.
    assert widths['qt6']['P_width_percent'] < 2.
