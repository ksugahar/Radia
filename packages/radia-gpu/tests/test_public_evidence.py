import importlib.util
from pathlib import Path


def test_public_poisson_keeps_numerics_and_omits_private_provenance():
    path = Path(__file__).resolve().parents[1] / 'prepare_ci_evidence.py'
    spec = importlib.util.spec_from_file_location('public_gpu_evidence', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = dict(schema='radia.ngsolve-cuda-poisson.v1', passed=True, taskmanager_used=False,
                  hostname='private-host', python='C:/private/python.exe',
                  error='C:/private/failure', cases=[{'gpu_true_relative_residual': 1e-12}],
                  build=dict(hostname='private-host', python='C:/private/python.exe',
                             compiler_configuration='C:/private/compiler.xml',
                             source_commit='abc', binary_sha256='def', release='2026a'))
    result = module.poisson_public(source)
    assert result['taskmanager_used'] is False
    assert module.ams_public(dict(taskmanager_used=False))['taskmanager_used'] is False
    assert result['cases'] == source['cases']
    assert result['build'] == dict(source_commit='abc', binary_sha256='def', release='2026a')
    assert not {'hostname', 'python', 'error'} & result.keys()
