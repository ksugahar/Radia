"""Frozen ESIM assembly contracts; mandatory native dependencies, no fixture skips."""
from pathlib import Path
import json
import numpy as np
import ngsolve as ng
import pytest
from netgen.occ import Box, OCCGeometry, Pnt
from radia.simulink import ih_operator_assembly as assembly
from radia.panels.calc_inductance import _extract_bnd_only_inline
from radia.surface_impedance import read_panel_impedance


def options(tmp_path, **kwargs):
    bh = tmp_path/'bh.txt'
    bh.write_text('0 0\n100 0.01\n')
    return assembly.IHOperatorAssemblyOptions(zs_mode='per-panel-esim', esim_bh_file=str(bh), **kwargs).checked()


@pytest.mark.parametrize('coupling', ['weak', 'strong'])
@pytest.mark.parametrize('evaluator', ['table', 'direct'])
def test_existing_esim_solver_options_and_reference_are_forwarded(tmp_path, coupling, evaluator):
    opts=options(tmp_path, coupling_mode=coupling, esim_panel_evaluator=evaluator, esim_reference_current_A=3)
    argv=assembly._unit_current_argv(tmp_path/'wp.vol', tmp_path/'coil.vol', 'bem-a', opts, tmp_path/'field.msh')
    assert '--esim-per-panel' in argv
    for flag, value in [('--current','3'),('--impedance-model','esim'),('--bh-file',opts.esim_bh_file),('--esim-panel-evaluator',evaluator)]:
        assert argv[argv.index(flag)+1] == value


@pytest.mark.parametrize('kwargs,match', [({'esim_reference_current_A':0},'reference'),({'esim_half_thickness_m':-1},'thickness'),({'esim_tolerance':float('nan')},'tolerance'),({'esim_max_iter':0},'budget|positive'),({'esim_drive_relative_band':.2},'band'),({'frequency_hz':-1},'frequency'),({'esim_reference_phase_rad':.2},'zero-phase')])
def test_invalid_esim_contract_fails_before_solve(tmp_path, kwargs, match):
    with pytest.raises(ValueError,match=match): options(tmp_path,**kwargs)


def test_missing_material_and_conflicting_modes_fail(tmp_path):
    with pytest.raises(ValueError,match='B-H'):
        assembly.IHOperatorAssemblyOptions(zs_mode='per-panel-esim').checked()
    bh=tmp_path/'bh'; bh.write_text('x')
    zs=tmp_path/'zs'; zs.write_text('{}')
    with pytest.raises(ValueError,match='also'):
        assembly.IHOperatorAssemblyOptions(zs_mode='per-panel-esim',esim_bh_file=str(bh),panel_zs_file=str(zs)).checked()


def test_certified_artifact_is_mesh_frequency_and_material_bound(tmp_path):
    shape=Box(Pnt(0,0,0),Pnt(.01,.01,.01)); shape.faces.name='sibc'
    vol=tmp_path/'wp.vol'; OCCGeometry(shape).GenerateMesh(maxh=.008).Save(str(vol))
    opts=options(tmp_path)
    surface=_extract_bnd_only_inline(ng.Mesh(str(vol)),'sibc')
    z=np.full(surface.GetNE(ng.BND),.01+.02j)
    payload=dict(current_A=opts.esim_reference_current_A,esim_converged=True,esim_per_panel=True,esim_impedance_layout='BND-element-order',esim_fixed_point_relative_error=1e-5,esim_iterations=3,esim_panel_evaluation=dict(mode='table',final_certification='direct-all-panels'),esim_per_panel_Z_s_real=z.real.tolist(),esim_per_panel_Z_s_imag=z.imag.tolist())
    provenance=assembly._write_esim_panel_artifact(vol,opts,payload,tmp_path)
    artifact=Path(provenance['artifact_file'])
    assert provenance['artifact_sha256']==assembly._sha256(artifact)
    assert provenance['bh_sha256']==assembly._sha256(Path(opts.esim_bh_file))
    assert provenance['zs_semantics']=='frozen_at_reference_current'
    assert np.array_equal(read_panel_impedance(artifact,surface,frequency_hz=opts.frequency_hz).values,z)
    with pytest.raises(ValueError,match='frequency'):
        read_panel_impedance(artifact,surface,frequency_hz=2*opts.frequency_hz)
    payload['wp_genus']=1
    with pytest.raises(RuntimeError,match='active genus-1'):
        assembly._write_esim_panel_artifact(vol,opts,payload,tmp_path)
    payload['wp_genus']=0
    payload['current_A']=2
    with pytest.raises(RuntimeError,match='operating current'):
        assembly._write_esim_panel_artifact(vol,opts,payload,tmp_path)
    payload['current_A']=opts.esim_reference_current_A
    for invalid in (-1, float('nan')):
        payload['esim_fixed_point_relative_error']=invalid
        with pytest.raises(RuntimeError,match='certification'):
            assembly._write_esim_panel_artifact(vol,opts,payload,tmp_path)
    payload['esim_fixed_point_relative_error']=1e-5
    for invalid in (1.5, True, 16):
        payload['esim_iterations']=invalid
        with pytest.raises(RuntimeError,match='certification'):
            assembly._write_esim_panel_artifact(vol,opts,payload,tmp_path)
    payload['esim_iterations']=3
    payload['esim_fixed_point_relative_error']=.1
    with pytest.raises(RuntimeError,match='certification'):
        assembly._write_esim_panel_artifact(vol,opts,payload,tmp_path)
