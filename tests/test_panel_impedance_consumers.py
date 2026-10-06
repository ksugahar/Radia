"""New triangle fields must fail before a legacy consumer writes nodal SOLs."""
import importlib.util
import json
from pathlib import Path

import pytest


def test_legacy_sol_export_refuses_triangle_layout_before_writing(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[1]
    path = root/'validation_test/ih_esim_benchmark/save_per_panel_sols.py'
    spec = importlib.util.spec_from_file_location('legacy_sol_export', path)
    consumer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(consumer)
    monkeypatch.setattr(consumer, 'SWEEP', tmp_path)
    volume = tmp_path/'not-read.vol'
    volume.write_text('mesh access must not occur')
    monkeypatch.setattr(consumer, 'VOL', volume)
    (tmp_path/'new_per_panel.json').write_text(json.dumps({
        'esim_impedance_layout':'BND-element-order',
        'esim_per_panel_Z_s_real':[1]*8, 'esim_per_panel_Z_s_imag':[1]*8,
        'esim_per_panel_H_t':[1]*8}))
    with pytest.raises(ValueError, match='legacy vertex consumer'):
        consumer.save_sols_for_case('new')
    assert not list(tmp_path.glob('*.sol'))
