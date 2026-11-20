import importlib.util
from pathlib import Path

import pytest

TOOL = Path(__file__).resolve().parents[1] / 'tools' / 'release_quad.py'
SPEC = importlib.util.spec_from_file_location('release_quad_state', TOOL)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


@pytest.mark.parametrize('matched', [True, False])
def test_preflight_checks_only_solver_version_pair(monkeypatch, matched):
    versions = dict.fromkeys([
        'radia', 'radia.__version__', 'radia-mcp'], '1.0.0')
    if not matched:
        versions['radia.__version__'] = '1.0.1'
    monkeypatch.setattr(module, '_read_repo_versions', lambda *_args: versions)
    checked_main = []
    monkeypatch.setattr(module, '_check_main_synced',
                        lambda **_: checked_main.append(True))
    assert module.cmd_preflight(None) == (0 if matched else 2)
    assert bool(checked_main) == matched


def test_compute_host_installer_passes_metadata_pins_to_native_checks(monkeypatch):
    """Guards: dependency identity precedes installation. Assumes: synthetic pins, no remote action."""
    import ast
    monkeypatch.setattr(module, "_solver_dependency_pins",
                        lambda: {"ngsolve": "sentinel-ng", "netgen-mesher": "sentinel-netgen"})
    tree = ast.parse(module._editable_install_program(
        {"source": "C:/Sources/test", "commit": "a" * 40,
         "python": "C:/Runtime/python.exe", "binding": {"hostname": "test"}}))
    args = next(ast.literal_eval(node.value) for node in tree.body
                if isinstance(node, ast.Assign) and node.targets[0].id == "native_args")
    assert args == ["native-check", "C:/Sources/test", "a" * 40,
                    "sentinel-ng", "sentinel-netgen"]


def test_solver_dependency_versions_come_from_exact_metadata_pins():
    assert module._solver_dependency_pins() == {
        "ngsolve": "6.2.2608",
        "netgen-mesher": "6.2.2608",
    }


def test_retired_simulink_candidate_route_is_absent():
    """Guards: the retired CPU MEX candidate route stays deleted; radia-simulink owns its release."""
    for name in ("cmd_simulink_candidate", "_run_simulink_candidate_target",
                 "_verify_simulink_candidate_state", "_write_simulink_state", "SIMULINK_TARGETS"):
        assert not hasattr(module, name), name
