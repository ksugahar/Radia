"""Pure audit arithmetic distinguishes throughput closure from storage accuracy."""
import importlib.util
from pathlib import Path
import sys

import pytest

SOURCE = Path(__file__).resolve().parents[1] / 'src/radia/ih_heat_transient.py'
_spec = importlib.util.spec_from_file_location('_ih_audit_contract', SOURCE)
iht = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = iht
_spec.loader.exec_module(iht)


@pytest.mark.parametrize('incoming,loss,stored,expected', [
    (1000., 1000., 0., 0.),
    (1000., 1000., .001, 1e-6),
    (0., 0., 1., 1.),
    (0., 0., 0., 0.),
    (0., 10., -10., 0.),
])
def test_energy_audit_handles_steady_state_and_zero_input(incoming, loss, stored, expected):
    audit = iht.TransientAudit(energy_in_J=incoming, energy_loss_J=loss,
                              energy_stored_J=stored)
    assert audit.as_dict()['energy_balance_relative_error'] == pytest.approx(expected)


@pytest.mark.parametrize('incoming,loss,stored,defect,closure,storage_error', [
    (1e6, .99e6, 10050., 50., 5e-5, 50./10050.),
    (1000., 999., .9, -.1, -1e-4, -.1),
    (1000., 1000., .001, .001, 1e-6, 1.),
    (0., 0., 0., 0., 0., 0.),
])
def test_storage_diagnostic_does_not_silently_change_throughput_gate(
        incoming, loss, stored, defect, closure, storage_error):
    audit = iht.TransientAudit(energy_in_J=incoming, energy_loss_J=loss,
                              energy_stored_J=stored)
    result = audit.as_dict()
    assert result['energy_balance_residual_J'] == pytest.approx(defect)
    assert result['energy_balance_relative_error'] == pytest.approx(closure)
    assert result['energy_storage_relative_error'] == pytest.approx(storage_error)
    assert result['energy_balance_normalization'] == 'max(abs(input), abs(loss), abs(stored))'
    assert result['energy_storage_normalization'] == 'max(abs(input-loss), abs(stored))'
    # Exercise the actual gate without constructing a numerical stepper.
    stepper = object.__new__(iht.NonlinearHeatStepper)
    stepper.audit = audit
    assert stepper.check_energy() == pytest.approx(closure)
