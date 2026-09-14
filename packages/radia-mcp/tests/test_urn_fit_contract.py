"""MCP boundary checks only: no optimizer or private research data required."""

import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest
from radia_mcp.radia_ngsolve.knowledge import urn


@pytest.mark.parametrize("freqs,z,options", [
    ([], [], {}), ([1], [1j], {}), ([[1, 2]], [[1j, 2j]], {}),
    ([1, 2], [1j], {}), ([0, 2], [1j, 2j], {}),
    ([-1, 2], [1j, 2j], {}), ([1, float('nan')], [1j, 2j], {}),
    ([1, 2], [1j, complex('inf')], {}), ([1, 2], [0, 0], {}),
    ([1, 2], [1j, 2j], {'n_epochs': 0}),
    ([1, 2], [1j, 2j], {'n_restarts': 1.5}),
    ([1, 2], [1j, 2j], {'n_debye': -1}),
    ([1, 2], [1j, 2j], {'n_debye': True}),
    ([1, 2], [1j, 2j], {'n_debye': 0, 'n_cole_cole': 0, 'n_warburg': 0}),
    ([1, 2], [1j, 2j], {'sparsity_weight': float('nan')}),
])
def test_invalid_fit_inputs_fail_before_optional_solver_import(monkeypatch, freqs, z, options):
    monkeypatch.setitem(sys.modules, 'torch', None)
    monkeypatch.setitem(sys.modules, 'radia.urn', None)
    with pytest.raises(ValueError):
        urn.run_urn_fit(freqs, z, **options)


@pytest.mark.parametrize('content,options', [
    ('1,2,3\n', {}), ('1,2\n2,3\n', {}),
    ('1,2,3\n2,3,4\n', {'freq_col': -1}),
    ('1,2,3\n2,3,4\n', {'real_col': 0}),
    ('1,2,3\n2,3,4\n', {'skip_rows': 0.5}),
])
def test_bad_csv_contract_does_not_start_training(tmp_path, monkeypatch, content, options):
    source = tmp_path/'input.csv'
    source.write_text(content, encoding='utf-8')
    monkeypatch.setattr(urn, 'run_urn_fit', lambda *a, **k: pytest.fail('training invoked'))
    with pytest.raises(ValueError):
        urn.urn_fit_from_csv(source, **options)


@pytest.mark.parametrize('export_fails', [False, True])
def test_csv_routes_to_legacy_model_and_only_writes_real_netlist(tmp_path, monkeypatch, export_fails):
    calls = []
    backend = ModuleType('radia.urn')
    backend.URNConfig = lambda **kwargs: kwargs

    def train(freqs, values, cfg, verbose):
        calls.append((freqs.copy(), values.copy(), cfg))
        model = type('Model', (), {
            '__call__': lambda self, omega: SimpleNamespace(detach=lambda: SimpleNamespace(numpy=lambda: values)),
            'get_active_components': lambda self: {},
        })()
        return model

    def export(*args):
        if export_fails:
            raise RuntimeError('unsupported realization')
        return '* synthetic test circuit\n.end\n'

    backend.train_urn = train
    backend.generate_spice_netlist = export
    monkeypatch.setitem(sys.modules, 'radia', ModuleType('radia'))
    monkeypatch.setitem(sys.modules, 'radia.urn', backend)
    monkeypatch.setitem(sys.modules, 'torch', SimpleNamespace(tensor=np.asarray, float64=np.float64))
    source, output = tmp_path/'input.csv', tmp_path/'out.cir'
    source.write_text('f,re,im\n2,4,6\n1,3,5\n1,3.1,5.1\n', encoding='utf-8')
    if export_fails:
        output.write_text('existing circuit', encoding='utf-8')
    report = urn.urn_fit_from_csv(source, skip_rows=1, spice_out=output)
    assert calls[0][2]['n_debye'] == 3
    np.testing.assert_array_equal(calls[0][0], [2, 1, 1])  # Keep order and repeated observations.
    np.testing.assert_array_equal(calls[0][1], [4+6j, 3+5j, 3.1+5.1j])
    assert 'not maximum error' in report
    assert 'not uniquely identified physical parts' in report
    if export_fails:
        assert 'SPICE synthesis failed' in report
        assert output.read_text() == 'existing circuit'
    else:
        assert output.read_text() == '* synthetic test circuit\n.end\n'


def test_research_settings_are_not_mcp_defaults():
    overview = urn.get_urn_documentation('overview')
    assert 'not the urn_fit default' in overview
    assert 'does not invoke Y-domain fitting' in overview
    assert 'not universal defaults' in overview


def test_registered_tool_keeps_the_same_csv_boundary(tmp_path):
    from radia_mcp.radia_ngsolve.server import urn_fit

    source = tmp_path/'one_row.csv'
    source.write_text('1,2,3\n', encoding='utf-8')
    with pytest.raises(ValueError, match='at least two rows'):
        urn_fit(str(source))
