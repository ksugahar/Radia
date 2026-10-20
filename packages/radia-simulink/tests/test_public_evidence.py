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
    setup = module.setup_public(dict(source, schema='radia.ams-setup-comparison.v1',
        accepted=True, timing_classification='diagnostic',formal_performance_acceptance='pending'))
    assert setup['cases'] == source['cases']
    assert setup['taskmanager_used'] is False
    assert setup['formal_performance_acceptance'] == 'pending'
    assert not {'hostname','host','python','error'} & setup.keys()
    fmm = module.fmm_public(dict(source, accepted=True,
        metrics={'field_relative':1e-9}, report={'reload_stale_handle_rejected':True}))
    assert fmm['metrics']['field_relative'] == 1e-9
    assert fmm['report']['reload_stale_handle_rejected'] is True
    assert not {'hostname','host','python','error'} & fmm.keys()
