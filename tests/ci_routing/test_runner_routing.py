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
    paths.update(p for profile in tier["profiles"].values() for p in profile["paths"])
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
    # A deliberate mode/runner promotion re-audits the changed job definition.
    declaration = m['light']['radia.yml:fast-contracts']
    paths, _ = routing.light_sources(root, declaration)
    declaration['audit_sha256'] = routing.audit_digest(root, routing.workflows(root)['radia.yml:fast-contracts'], paths, declaration)
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


@pytest.mark.parametrize("mutation", ["tier-format", "description", "non-light-budget", "ownership-rule", "existing-alternative"])
def test_validated_data_only_changes_preserve_executable_seal(tree, mutation):
    root, manifest = tree
    tier_path = root / "tests/test_tier_manifest.json"
    tier = json.loads(tier_path.read_text(encoding="utf-8"))
    if mutation == "description":
        tier["profiles"]["fast-contracts"]["description"] = "DATA description only"
    elif mutation == "non-light-budget":
        tier["profiles"]["sparsesolv"]["max_elapsed_seconds"] += 1
    elif mutation == "existing-alternative":
        tier["impact_rules"]["src/synthetic_C2.py"] = [tier["profiles"]["fast-contracts"]["paths"][0]]
    elif mutation == "ownership-rule":
        path = root / "tools/ci_distribution_manifest.json"
        owners = json.loads(path.read_text(encoding="utf-8"))
        owners["rules"].append({"pattern": "validation_records/C2_synthetic/**", "owners": ["radia"]})
        path.write_text(json.dumps(owners), encoding="utf-8")
    tier_path.write_text(json.dumps(tier, sort_keys=True, separators=(",", ":")), encoding="utf-8")
    assert routing.validate(root, manifest)["jobs"] == len(manifest["jobs"])


@pytest.mark.parametrize("mutation", ["new-test", "existing-fast", "repo-to-fast", "new-import", "conftest", "requirements", "guard-exception", "active-budget", "job-timeout"])
def test_genuine_execution_changes_still_require_f2_review(tree, mutation):
    root, manifest = tree
    declaration = manifest["light"]["radia.yml:fast-contracts"]
    path = root / "tests/test_tier_manifest.json"
    tier = json.loads(path.read_text(encoding="utf-8"))
    if mutation == "new-test":
        (root / "tests/test_C2_new.py").write_text("def test_new(): pass\n", encoding="utf-8")
        tier["impact_rules"]["src/synthetic_C2.py"] = ["tests/test_C2_new.py"]
    elif mutation in {"existing-fast", "repo-to-fast"}:
        candidate = "tests/test_ci_execution_policy.py"
        assert (root / candidate).is_file()
        if mutation == "repo-to-fast":
            candidate = "tests/test_validation_records.py"
            tier["profiles"]["repo-contracts"]["paths"].remove(candidate)
        tier["profiles"]["fast-contracts"]["paths"].append(candidate)
    elif mutation == "new-import":
        selected = root / tier["profiles"]["fast-contracts"]["paths"][0]
        (selected.parent / "C2_helper.py").write_text("VALUE = 1\n", encoding="utf-8")
        selected.write_text(selected.read_text(encoding="utf-8") + "\nimport C2_helper\n", encoding="utf-8")
    elif mutation == "conftest":
        path2 = root / "tests/conftest.py"
        path2.write_text(path2.read_text(encoding="utf-8") + "\n# changed execution context\n", encoding="utf-8")
    elif mutation == "requirements":
        declaration["requirements"].append("synthetic-light-dependency")
    elif mutation == "guard-exception":
        declaration["deferred_or_mocked_imports"]["src/radia/ih_heat_transient.py"].append("synthetic_module")
    elif mutation == "active-budget":
        tier["profiles"]["fast-contracts"]["max_elapsed_seconds"] += 1
    else:
        alter_job(root, lambda jobs: jobs["fast-contracts"].update({"timeout-minutes": 99}))
    path.write_text(json.dumps(tier), encoding="utf-8")
    with pytest.raises(ValueError, match="closure changed"):
        routing.validate(root, manifest)


@pytest.mark.parametrize("mutation", ["unknown-field", "unsafe-path", "cycle", "missing-file", "unknown-owner", "missing-coverage-contract"])
def test_data_contracts_fail_independently_of_raw_byte_sealing(tree, mutation):
    root, manifest = tree
    path = root / "tests/test_tier_manifest.json"
    tier = json.loads(path.read_text(encoding="utf-8"))
    if mutation == "unknown-field":
        tier["profiles"]["fast-contracts"]["run"] = "arbitrary command"
    elif mutation == "unsafe-path":
        tier["impact_rules"]["../outside.py"] = []
    elif mutation == "cycle":
        tier["profiles"]["fast-contracts"]["extends"] = "fast-contracts"
    elif mutation == "missing-file":
        tier["profiles"]["solver-heavy"]["paths"].append("tests/nonexistent_C2.py")
    elif mutation == "unknown-owner":
        op = root / "tools/ci_distribution_manifest.json"
        owners = json.loads(op.read_text(encoding="utf-8"))
        owners["rules"].append({"pattern": "validation_records/C2_synthetic/**", "owners": ["unknown-owner"]})
        op.write_text(json.dumps(owners), encoding="utf-8")
    else:
        tier["profiles"]["repo-contracts"]["paths"].remove("tests/test_test_tier_policy.py")
    path.write_text(json.dumps(tier), encoding="utf-8")
    with pytest.raises(ValueError):
        routing.validate(root, manifest)


def test_actual_light_selection_cannot_bypass_audited_union():
    spec = importlib.util.spec_from_file_location("C2_runner", ROOT / "tools/run_test_tier.py")
    runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
    with pytest.raises(ValueError, match="outside audited alternatives"):
        runner.check_light_selection(["tests/test_ci_execution_policy.py"])


def test_other_json_inputs_are_not_blanket_exempted(tree):
    root, manifest = tree
    declaration = manifest["light"]["radia.yml:fast-contracts"]
    path = root / "tools/C2_reviewed_input.json"
    path.write_text('{"value":1}', encoding="utf-8")
    declaration["audited_dynamic_sources"].append("tools/C2_reviewed_input.json")
    paths, _ = routing.light_sources(root, declaration)
    declaration["audit_sha256"] = routing.audit_digest(root, routing.workflows(root)["radia.yml:fast-contracts"], paths, declaration)
    path.write_text('{"value":2}', encoding="utf-8")
    with pytest.raises(ValueError, match="closure changed"):
        routing.validate(root, manifest)
