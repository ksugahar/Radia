"""Routing mutations test actual workflows, entrypoints and selector closure."""
import copy
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('routing_contract', ROOT/'tools/ci_runner_routing.py')
routing = importlib.util.module_from_spec(spec)
spec.loader.exec_module(routing)


@pytest.fixture
def tree(tmp_path):
    manifest = json.loads((ROOT/'tools/ci_runner_routing.json').read_text(encoding="utf-8-sig"))
    paths, _ = routing.light_sources(ROOT, manifest['light']['radia.yml:fast-contracts'])
    tier=json.loads((ROOT/"tests/test_tier_manifest.json").read_text(encoding="utf-8"))
    paths.update(p for values in tier["impact_rules"].values() for p in values)
    paths |= {'tests/test_tier_manifest.json','tools/ci_runner_routing.json'}
    paths |= {p.relative_to(ROOT).as_posix() for p in (ROOT/'.github/workflows').glob('*.y*ml')}
    for path in paths:
        dst = tmp_path/path
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/path, dst)
    return tmp_path, manifest


def alter_job(root, callback):
    path = root/'.github/workflows/radia.yml'
    data = yaml.safe_load(path.read_text(encoding="utf-8-sig"))
    callback(data['jobs'])
    path.write_text(yaml.safe_dump(data))


def test_stage3_real_inventory(tree):
    root, m = tree
    assert routing.validate(root,m)['jobs'] == len(m['jobs'])


@pytest.mark.parametrize('labels', [
    ['self-hosted','Windows','X64','100'],
    ['self-hosted','Windows','X64','100','ci-light'],
    ['self-hosted','Windows','X64','ci-light','ci-native'],
    ['self-hosted','Windows','X64','ci-matlab-gpu'],
])
def test_wrong_light_route(tree,labels):
    root,m=tree
    alter_job(root,lambda jobs:jobs['fast-contracts'].update({'runs-on':labels}))
    with pytest.raises(ValueError,match='capability'):routing.validate(root,m)


def test_unknown_actual_job_rejected(tree):
    root,m=tree
    alter_job(root,lambda jobs:jobs.update({'unclassified':{'runs-on':'ubuntu-latest','steps':[]}}))
    with pytest.raises(ValueError,match='unclassified'):routing.validate(root,m)


@pytest.mark.parametrize('command',['matlab -batch version','python -c "import matlab.engine"','nvcc kernel.cu','cmake --build build','pwsh ./Build.ps1','python hidden.py'])
def test_heavy_or_unaudited_direct_command(tree,command):
    root,m=tree
    alter_job(root,lambda jobs:jobs['fast-contracts']['steps'].append({'run':command}))
    with pytest.raises(ValueError):routing.validate(root,m)


def test_transitive_entrypoint_matlab_mutation(tree):
    root,m=tree
    path=root/'validation_test/esrf_three_engine/check_case6_inputs.py'
    path.write_text(path.read_text(encoding="utf-8-sig")+'\nimport matlab.engine\n',encoding='utf-8')
    with pytest.raises(ValueError,match='heavy light dependency'):routing.validate(root,m)


def test_transitive_native_compile_mutation(tree):
    root,m=tree
    path=root/'tools/audit_taskmanager.py'
    path.write_text(path.read_text(encoding="utf-8-sig")+"\nimport subprocess\nsubprocess.run(['nvcc','kernel.cu'])\n",encoding='utf-8')
    with pytest.raises(ValueError,match='closure changed'):routing.validate(root,m)


def test_heavy_requirement_rejected(tree):
    root,m=tree
    m['light']['radia.yml:fast-contracts']['requirements'].append('matlabengine')
    with pytest.raises(ValueError,match='requirement'):routing.validate(root,m)


def test_final_requires_host_and_exact_capability(tree):
    root,m=tree
    with pytest.raises(ValueError,match='capability'):routing.validate(root,m,'final')
    for file in (root/'.github/workflows').glob('*.y*ml'):
        data=yaml.safe_load(file.read_text(encoding="utf-8-sig"))
        for name,job in data['jobs'].items():
            rule=m['jobs'][f'{file.name}:{name}']
            if rule['target']!='github-hosted':
                job['runs-on']=['self-hosted','Windows','X64',rule['target']]
                if rule['target'] != 'gpu-a100': job['runs-on'].insert(3,'100')
        file.write_text(yaml.safe_dump(data))
    assert routing.validate(root,m,'final')['mode']=='final'


def test_source_only_conftest_does_not_probe_and_rejects_native():
    # Fresh process, intercept optional probes before executing actual conftest.
    code = '''import builtins,runpy,sys
original=builtins.__import__
def checked(name,*a,**kw):
 if name.split('.')[0] in {'ngsolve','netgen','cubit','cubit_mesh_export','matlab'}:
  raise AssertionError('native probe attempted:'+name)
 return original(name,*a,**kw)
builtins.__import__=checked
runpy.run_path('tests/conftest.py')
builtins.__import__=original
try: __import__('ngsolve')
except RuntimeError as e: assert 'forbidden' in str(e)
else: raise AssertionError('native import accepted')
'''
    env={**os.environ,'RADIA_SOURCE_ONLY_CONTRACTS':'1'}
    result=subprocess.run([sys.executable,'-c',code],cwd=ROOT,env=env,capture_output=True,text=True)
    assert result.returncode==0,result.stderr

@pytest.mark.parametrize('changed',['src/radia/lamination.py','tests/test_lamination_api.py'])
def test_lamination_native_coverage_is_retained(changed):
    spec=importlib.util.spec_from_file_location('routing_impact',ROOT/'tools/sparsesolv_ci_impact.py')
    impact=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(impact)
    manifest=impact.manifest()
    workflow=yaml.safe_load((ROOT/'.github/workflows/sparsesolv.yml').read_text(encoding='utf-8'))
    for event in ('push','pull_request'):
        assert changed in workflow[True][event]['paths']
    assert impact.lane_required('sparsesolv',{changed},lambda p:(None,None),manifest)[0]
    assert 'tests/test_lamination_api.py' in json.loads((ROOT/'tests/test_tier_manifest.json').read_text(encoding='utf-8'))['profiles']['solver-numerics']['paths']
    # Removing the trigger from the projection actually loses native admission.
    mutant=copy.deepcopy(manifest)
    for event in ('push','pull_request'):
        mutant['lanes']['sparsesolv']['events'][event]['paths'].remove(changed)
        mutant['lanes']['sparsesolv']['events'][event]['paths'].append('!' + changed)
    assert not impact.lane_required('sparsesolv',{changed},lambda p:(None,None),mutant)[0]

@pytest.mark.parametrize('key', ['radia-simulink.yml:cuda-mex', 'radia-optuna.yml:build-test', 'radia-optuna.yml:installed-wheel-matlab-e2e'])
def test_gpu_jobs_cannot_lose_rtx_capability(tree, key):
    root, manifest = tree
    filename, name = key.split(':')
    path = root / '.github/workflows' / filename
    document = yaml.safe_load(path.read_text(encoding='utf-8'))
    document['jobs'][name]['runs-on'].remove('gpu-rtx')
    path.write_text(yaml.safe_dump(document), encoding='utf-8')
    with pytest.raises(ValueError, match='capability'):
        routing.validate(root, manifest)


def test_a100_is_compute_only_and_never_a_pr_gate():
    workflow = yaml.safe_load((ROOT / '.github/workflows/gpu-validation.yml').read_text(encoding='utf-8'))
    assert set(workflow[True]) == {'schedule', 'workflow_dispatch'}
    job = workflow['jobs']['validate']
    assert job['runs-on'] == ['self-hosted', 'Windows', 'X64', 'gpu-a100']
    commands = '\n'.join(step.get('run', '') for step in job['steps'])
    assert 'ci_gpu_artifact.py fetch' in commands
    assert '--sha ${{ steps.source.outputs.sha }}' in commands
    assert "$sha = '${{ github.sha }}'" in commands
    assert "github.event_name }}' -eq 'schedule'" in commands
    assert 'ci_gpu_artifact.py select' in commands
    assert workflow['concurrency']['group'].endswith('${{ github.event_name }}')
    assert "needs.validate.result == 'failure'" in workflow['jobs']['nightly-failure']['if']
    assert not any(token in commands.lower() for token in ['build.ps1', 'mexcuda', 'cmake', 'pip install build'])
    assert 'issues: write' in (ROOT / '.github/workflows/gpu-validation.yml').read_text(encoding='utf-8')
