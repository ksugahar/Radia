"""Guards: numerical workflows retain native provenance and required execution paths.

Assumes: source-only workflow and selector contracts; no native import or numerical
accuracy is certified by these checks.
"""

import json
import importlib.util
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SUITE = "src/ext/sparsesolv/tests/test_sparsesolv.py"
SPEC = importlib.util.spec_from_file_location("ams_impact", ROOT / "tools/sparsesolv_ci_impact.py")
IMPACT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(IMPACT)
PROJECT = '''[build-system]
requires = ["setuptools"]
[project]
dependencies = ["ngsolve==6.2.2606"]
[project.optional-dependencies]
test = ["pytest"]
dev = ["ruff"]
beam = ["scipy"]
'''


def test_sparsesolv_is_in_native_release_and_focused_tiers_only():
    profiles = json.loads((ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8"))["profiles"]
    assert SUITE in profiles["native-smoke"]["paths"]
    parallel = "tests/test_ams_internal_setup_parallel.py"
    additional = ["src/ext/sparsesolv/tests/test_hermitian_iccg.py",
                  "src/ext/sparsesolv/tests/test_ams_repeatability.py"]
    assert profiles["sparsesolv"]["paths"] == [SUITE, *additional, parallel]
    assert all(path in profiles["native-smoke"]["paths"] for path in additional)
    assert parallel in profiles["native-smoke"]["paths"]
    assert SUITE not in profiles["fast-contracts"]["paths"]


def test_solver_workflow_builds_before_testing_on_mdx():
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text(encoding="utf-8"))
    triggers = workflow.get("on", workflow.get(True))
    for event in ("push", "pull_request"):
        assert "src/ext/sparsesolv/**" in triggers[event]["paths"]
        assert "pyproject.toml" in triggers[event]["paths"]
    job = workflow["jobs"]["ams-regression"]
    assert "100" in job["runs-on"]
    for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        assert job["env"][variable] == "1", "Bound BLAS before subprocess imports"
    assert "head.repo.full_name == github.repository" in job["if"]
    script = "\n".join(step.get("run", "") for step in job["steps"])
    assert script.index("--target sparsesolv_ngsolve") < script.index("--profile sparsesolv")
    assert "s.has_compact_ams()" in script
    assert "p.parent ==" in script


def test_ams_unavailable_is_not_a_silent_skip():
    source = (ROOT / SUITE).read_text(encoding="utf-8")
    assert "pytest.skip(" not in source
    assert source.count("assert has_compact_ams()") == 3


@pytest.mark.parametrize("old,new,skip", [
    ('test = ["pytest"]', 'test = ["pytest", "setuptools"]', True),
    ('dev = ["ruff"]', 'dev = []', True),
    ('beam = ["scipy"]', 'beam = []', False),
    ('6.2.2606', '6.2.2607', False),
    ('requires = ["setuptools"]', 'requires = ["setuptools>=80"]', False),
    ('test = ["pytest"]', 'test = "invalid"', False),
    ('test = ["pytest"]', 'test = ["pytest"] # comment', False),
])
def test_native_impact_only_exempts_test_and_dev_extras(old, new, skip):
    assert IMPACT.test_dependencies_only(PROJECT, PROJECT.replace(old, new)) is skip


def test_native_impact_preserves_mixed_configuration_changes():
    after = PROJECT.replace('["pytest"]', '[]').replace('6.2.2606', '6.2.2607')
    assert not IMPACT.test_dependencies_only(PROJECT, after)


@pytest.mark.parametrize("event,base", [("push", "a" * 40), ("pull_request", "HEAD^1")])
def test_native_impact_uses_complete_event_range(event, base):
    calls = []
    def git(*args):
        calls.append(args)
        return {("diff", "--name-only", "-z", base, "HEAD", "--"): "pyproject.toml\0",
                ("show", f"{base}:pyproject.toml"): PROJECT,
                ("show", "HEAD:pyproject.toml"): PROJECT.replace('["pytest"]', '[]')}[args]
    assert IMPACT.native_required(event, {"before": base}, git)[0] is False
    assert calls[0][3] == base


def test_markdown_only_change_skips_native_ams():
    docs = "README.md\0src/ext/sparsesolv/docs/algorithms.md\0"
    assert IMPACT.native_required("pull_request", {}, lambda *args: docs) == (
        False, "documentation-only change")
    mixed = "README.md\0src/ext/sparsesolv/test.cpp\0"
    assert IMPACT.native_required("pull_request", {}, lambda *args: mixed)[0]


@pytest.mark.parametrize("paths", ["", "pyproject.toml\0Build.ps1\0", "src/ext/sparsesolv/test.cpp\0"])
def test_native_impact_does_not_exempt_other_files(paths):
    assert IMPACT.native_required("push", {"before": "a" * 40}, lambda *args: paths)[0]


def test_manual_native_impact_always_runs():
    assert IMPACT.native_required("workflow_dispatch", {}, None)[0]
    assert IMPACT.native_required("push", {"before": "0" * 40}, None)[0]


@pytest.mark.parametrize("error", [ValueError("invalid TOML"), OSError("missing git base")])
def test_unknown_native_impact_fails_closed(tmp_path, monkeypatch, error):
    event = tmp_path / "event.json"
    event.write_text('{"before":"base"}')
    output = tmp_path / "output"
    monkeypatch.setenv("GITHUB_EVENT_PATH", str(event))
    monkeypatch.setenv("GITHUB_EVENT_NAME", "push")
    monkeypatch.setenv("GITHUB_OUTPUT", str(output))
    def fail(*args):
        raise error
    monkeypatch.setattr(IMPACT, "required_for_event", fail)
    IMPACT.main([])
    assert output.read_text(encoding="utf-8").startswith("required=true\nreason=")
    assert "fail closed" in output.read_text(encoding="utf-8")


def test_native_impact_gates_expensive_steps_and_watches_itself():
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text(encoding="utf-8"))
    triggers = workflow.get("on", workflow.get(True))
    for event in ("push", "pull_request"):
        assert "tools/sparsesolv_ci_impact.py" in triggers[event]["paths"]
        assert set(triggers[event]["paths"]) == set(IMPACT.NATIVE_PATHS) | {
            "pyproject.toml", "tests/test_tier_manifest.json", "!**/*.md"}
        assert triggers[event]["paths"][-1] == "!**/*.md"
    steps = workflow["jobs"]["ams-regression"]["steps"]
    selector = next(s for s in steps if s.get("id") == "native-impact")
    assert "python tools/ci_numerical_plan.py --github-event" in selector["run"]
    for step in steps:
        if step.get("name", "").startswith(("Create isolated", "Build ", "Verify source", "Upload solver")):
            assert "steps.native-impact.outputs.required == 'true'" in step["if"]


def manifest():
    return {"schema": "radia.test-tier-manifest.v1", "impact_rules": {}, "profiles": {
        "sparsesolv": {"extends": "base", "paths": [SUITE], "max_elapsed_seconds": 30},
        "solver-numerics": {"paths": []}, "solver-heavy": {"paths": []},
        "native-build-contracts": {"paths": []},
        "base": {"paths": ["tests/base.py"]}, "fast": {"paths": ["tests/old.py"]}}}


@pytest.mark.parametrize("change,expected", [
    (lambda d: d["impact_rules"].update({"tools/new.py": ["tests/new.py"]}), False),
    (lambda d: d["profiles"]["fast"]["paths"].append("tests/new.py"), True),
    (lambda d: d["profiles"]["sparsesolv"]["paths"].append("tests/new.py"), False),
    (lambda d: d["profiles"]["base"]["paths"].append("tests/new.py"), False),
    (lambda d: d["profiles"]["sparsesolv"].update(max_elapsed_seconds=60), False),
    (lambda d: d["profiles"]["sparsesolv"].update(extends="fast"), False),
])
def test_manifest_compares_ams_and_all_ancestors(change, expected):
    before, after = manifest(), manifest()
    change(after)
    assert IMPACT.ams_tier_projection_unchanged(json.dumps(before), json.dumps(after)) is expected


@pytest.mark.parametrize("change", [
    lambda d: d.update(schema="unknown"),
    lambda d: d.update(future_global_option=True),
    lambda d: d["profiles"].pop("sparsesolv"),
    lambda d: d["profiles"].pop("base"),
    lambda d: d["profiles"]["base"].update(extends="sparsesolv"),
    lambda d: d["profiles"]["base"].update(paths="bad"),
    lambda d: d["profiles"]["base"].update(paths=[SUITE, SUITE]),
    lambda d: d["profiles"]["sparsesolv"].update(max_elapsed_seconds=-1),
])
def test_invalid_manifest_fails_closed(change):
    before, after = manifest(), manifest()
    change(after)
    with pytest.raises((KeyError, ValueError, TypeError)):
        IMPACT.ams_tier_projection_unchanged(json.dumps(before), json.dumps(after))


@pytest.mark.parametrize("extra,required", [
    ("tests/new.py\0tools/new.py\0", False),
    ("src/ext/sparsesolv/test.cpp\0", True),
    ("tests/matlab/test_beam_transfer_mex.m\0", False),
    ("src/core/rad_beam_transfer.cpp\0", True),
    ("packages/radia-simulink/matlab/+radia/+internal/callMex.m\0", False),
    ("tools/run_test_tier.py\0", True),
    (".github/workflows/sparsesolv.yml\0", True),
])
def test_manifest_edit_with_other_files(extra, required):
    path = "tests/test_tier_manifest.json"
    before, after = manifest(), manifest()
    # Exercise unrelated profile data; impact rules are read by native contracts.
    after["profiles"]["fast"]["description"] = "unrelated metadata"
    def git(*args):
        if args[0] == "diff":
            return path + "\0" + extra
        return json.dumps(before if args[1].startswith("base:") else after)
    assert IMPACT.native_required("push", {"before": "base"}, git)[0] is required


def test_mixed_manifest_and_runtime_dependency_change_still_runs():
    def git(*args):
        if args[0] == "diff":
            return "tests/test_tier_manifest.json\0pyproject.toml\0tests/new.py\0"
        if args[1].endswith("pyproject.toml"):
            return PROJECT if args[1].startswith("base:") else PROJECT.replace("6.2.2606", "6.2.2607")
        return json.dumps(manifest())
    assert IMPACT.native_required("push", {"before": "base"}, git)[0]


@pytest.mark.parametrize("source", [
    "src/radia/mixed_omega_newton.py", "src/radia/eddy_aphi.py",
    "src/radia/eddy_axisym_ring.py", "src/radia/sibc_corner_patch.py",
    "src/radia/vector_potential_solver.py", "src/radia/topology_optimization.py",
    "src/radia/ih_thermal.py", "src/radia/ih_axisym_coupled.py",
    "src/radia/panels/calc_fem_kelvin.py",
    "src/radia/coil_from_cad.py",
])
def test_numerical_source_changes_select_built_ngsolve_lane(source):
    assert IMPACT.native_required("push", {"before": "a" * 40},
                                  lambda *args: source + "\0")[0]
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text(encoding="utf-8"))
    steps = workflow["jobs"]["ams-regression"]["steps"]
    numerical = next(s for s in steps if s.get("name") == "Run required NGSolve numerical contracts")
    assert "--profile solver-numerics" in numerical["run"]
    assert numerical["env"]["RADIA_TESTS_ALLOW_PARTIAL"] == "0"
    assert "find('skipped') is None" in numerical["run"]



def test_ih_thermal_and_complex_ams_contracts_run_with_ngsolve():
    profiles = json.loads((ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8"))["profiles"]
    required = {"tests/test_ih_axisym_coupled.py", "tests/test_ih_thermal_transfer.py",
                "tests/test_ih_thermal_exposure.py",
                "validation_test/induction_heating/test_fem_ams_parity.py"}
    assert required <= set(profiles["solver-numerics"]["paths"])
    assert required.isdisjoint(profiles["fast-contracts"]["paths"])


def test_numerical_tier_changes_cannot_skip_native_lane():
    before, after = manifest(), manifest()
    after["profiles"]["solver-numerics"] = {"paths": ["tests/test_mixed_omega_newton.py"]}
    assert not IMPACT.ams_tier_projection_unchanged(json.dumps(before), json.dumps(after))


def test_heavy_derivative_validation_remains_required_after_tier_split():
    path = "validation_test/topology_optimization/test_hex_cluster_derivative.py"
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text(encoding="utf-8"))
    steps = workflow["jobs"]["ams-regression"]["steps"]
    heavy = next(s for s in steps if s.get("id") == "heavy-validation")
    assert "tools/validation_source_gate.py" in heavy["run"]
    gate = (ROOT / "tools/validation_source_gate.py").read_text(encoding="utf-8")
    assert "PROFILE = 'solver-heavy'" in gate
    assert "'tools/run_test_tier.py','--profile',PROFILE" in gate
    assert heavy["env"]["RADIA_TESTS_ALLOW_PARTIAL"] == "0"
    assert "('failure','error','skipped')" in gate
    assert "if not cases or any" in gate
    assert "steps.native-impact.outputs.required == 'true'" in heavy["if"]
    profiles = json.loads((ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8"))["profiles"]
    assert path in profiles["solver-heavy"]["paths"]
    assert path not in profiles["solver-numerics"]["paths"]
    assert IMPACT.native_required("push", {"before": "a" * 40},
                                  lambda *args: path + "\0")[0]
    before, after = manifest(), manifest()
    after["profiles"]["solver-heavy"] = {"paths": [path]}
    assert not IMPACT.ams_tier_projection_unchanged(json.dumps(before), json.dumps(after))


def test_cpu_lane_excludes_foreign_runtime_and_integration_gates():
    workflow = yaml.load((ROOT / ".github/workflows/sparsesolv.yml").read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    excluded = [
        "src/runtime/ih/radia_ih_thermal.cpp",
        "tests/cpp/test_radia_mex_sparse_symmetry.cpp",
        "tests/test_sparsesolv_ci_contract.py",
    ]
    for event in ("push", "pull_request"):
        for path in excluded:
            assert not IMPACT.CONTRACT.selected(path, workflow["on"][event]), path
    profiles = json.loads((ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8"))["profiles"]
    for name in ("solver-numerics", "solver-heavy", "native-build-contracts"):
        assert set(excluded).isdisjoint(profiles[name]["paths"])
    steps = workflow["jobs"]["ams-regression"]["steps"]
    scripts = "\n".join(step.get("run", "") for step in steps)
    assert "test_radia_mex_sparse_symmetry" not in scripts
    assert "test_radia_ih_thermal_acceptance" not in scripts
    numerical = next(step for step in steps if step.get("name") == "Run required NGSolve numerical contracts")
    assert "packages\\radia-mcp\\src" in numerical["env"]["PYTHONPATH"]


def test_native_solver_lane_builds_required_axisymmetric_extension():
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text(encoding="utf-8"))
    steps = workflow["jobs"]["ams-regression"]["steps"]
    build = next(s for s in steps if s.get("name") == "Build current package and SparseSolv")
    assert ".\\Build.ps1 -AxiFemOnly" in build["run"]
    assert "--target sparsesolv_ngsolve peec_matrices" in build["run"]
    verify = next(s for s in steps if s.get("name") == "Verify source provenance and run required AMS tests")
    assert "radia.peec_matrices" in verify["run"]
    assert "ExtensionFileLoader" in verify["run"]


def test_axifem_build_links_python_from_the_base_interpreter():
    script = (ROOT / "Build.ps1").read_text(encoding="utf-8")
    assert "sys.base_prefix" in script
    assert 'Join-Path $pyBasePrefix "libs\\$pyLib"' in script
    assert '"$PythonImportLibrary"' in script
    assert '"$pyPrefix\\libs\\$pyLib"' not in script


@pytest.mark.parametrize("changed", ["src/radia/p1_linear.py", "tests/test_p1_linear.py"])
def test_p1_linear_routes_to_solver_numerics_and_native_ci(changed):
    manifest = json.loads((ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8"))
    test = "tests/test_p1_linear.py"
    assert test in manifest["impact_rules"]["src/radia/p1_linear.py"]
    assert test in manifest["profiles"]["solver-numerics"]["paths"]
    assert changed in IMPACT.NATIVE_PATHS
    assert IMPACT.native_required(
        "push", {"before": "a" * 40}, lambda *args: changed + "\0"
    )[0]
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text(encoding="utf-8"))
    triggers = workflow.get("on", workflow.get(True))
    for event in ("push", "pull_request"):
        assert changed in triggers[event]["paths"]


ROOT_CASES = [
    ('pyproject.toml', '[project]\nversion="1"\n', '[project]\nversion="2"\n', {'radia'}),
    ('pyproject.toml', '[tool.setuptools.package-data]\nradia_simulink=["*.m"]\n',
     '[tool.setuptools.package-data]\nradia_simulink=["*.cu"]\n', {'radia-simulink'}),
    ('Build.ps1', "param()\nif ($OptunaMexOnly) { matlab -batch 'old' }\n",
     "param()\nif ($OptunaMexOnly) { matlab -batch 'new' }\n", {'radia-optuna'}),
    ('Build.ps1', 'if ($RadiaOnly) { python old }', 'if ($RadiaOnly) { python new }', {'radia'}),
    ('Build.ps1', 'if ($AmsCudaOnly) { python old }', 'if ($AmsCudaOnly) { python new }', {'radia-simulink'}),
    ('CMakeLists.txt', 'option(RADIA_BUILD_MATLAB_MEX "enable" OFF)',
     'option(RADIA_BUILD_MATLAB_MEX "enable" ON)', {'radia-simulink'}),
    ('CMakeLists.txt', 'target_compile_definitions(_radia_pybind PRIVATE OLD)',
     'target_compile_definitions(_radia_pybind PRIVATE NEW)', {'radia'}),
    ('MANIFEST.in', 'recursive-include packages/radia-simulink *.m\n',
     'recursive-include packages/radia-simulink *.cu\n', {'radia-simulink'}),
    ('tests/test_tier_manifest.json',
     json.dumps({'schema':'radia.test-tier-manifest.v1','profiles':{'gpu-unit':{
         'paths':['packages/radia-simulink/tests/test_ic_payload.py'],'max_elapsed_seconds':1}}}),
     json.dumps({'schema':'radia.test-tier-manifest.v1','profiles':{'gpu-unit':{
         'paths':['packages/radia-simulink/tests/test_ic_payload.py'],'max_elapsed_seconds':2}}}), {'radia-simulink'}),
]


CODE_LIKE = ('Build.ps1', 'CMakeLists.txt', 'MANIFEST.in')


@pytest.mark.parametrize('path,before,after,owners', [c for c in ROOT_CASES if c[0] not in CODE_LIKE])
def test_root_content_diff_selects_only_evidenced_owners(path,before,after,owners):
    assert IMPACT.root_owners(path,before,after,IMPACT.manifest())[0] == owners


@pytest.mark.parametrize('path,before,after,owners', [c for c in ROOT_CASES if c[0] in CODE_LIKE])
def test_code_like_root_changes_select_all_file_owners(path,before,after,owners):
    data=IMPACT.manifest()
    assert IMPACT.root_owners(path,before,after,data)[0] == set(IMPACT.CONTRACT.owner_rule(data,path)['owners'])


@pytest.mark.parametrize('path', ['pyproject.toml', 'tests/test_tier_manifest.json'])
def test_real_root_comments_and_whitespace_select_no_build(path):
    before=(ROOT/path).read_text(encoding='utf-8')
    after=before+'\n' if path.endswith('.json') else before+'\n# comment only\n'
    assert IMPACT.root_owners(path,before,after,IMPACT.manifest())[0] == set()


@pytest.mark.parametrize('path', ['Build.ps1', 'CMakeLists.txt', 'MANIFEST.in'])
def test_code_like_root_formatting_is_never_exempt(path):
    """PowerShell whitespace can change meaning, so even comment-only edits select all owners."""
    before=(ROOT/'Build.ps1').read_text(encoding='utf-8')
    after=before+'\n# comment only\n'
    assert IMPACT.root_owners('Build.ps1',before,after,IMPACT.manifest())[0] == set(IMPACT.CONTRACT.DISTRIBUTIONS)


@pytest.mark.parametrize('path,before,after', [
    ('pyproject.toml','[project]\nfuture=1','[project]\nfuture=2'),
    ('Build.ps1','$shared=1','$shared=2'),
    ('CMakeLists.txt','set(UNKNOWN 1)','set(UNKNOWN 2)'),
    ('MANIFEST.in','global-exclude *.so','global-exclude *.dll'),
    ('tests/test_tier_manifest.json','{}','{"unknown":true}'),
])
def test_unattributed_root_edits_select_all_file_owners(path,before,after):
    assert IMPACT.root_owners(path,before,after,IMPACT.manifest())[0] == set(IMPACT.CONTRACT.DISTRIBUTIONS)


@pytest.mark.parametrize('path,broken', [
    ('pyproject.toml','[broken'),('Build.ps1','if ($AmsCudaOnly) {'),
    ('CMakeLists.txt','add_library('),('MANIFEST.in','unknown directive'),
    ('tests/test_tier_manifest.json','{'),
])
def test_unparsable_deleted_or_unavailable_roots_fail_closed(path,broken):
    before=(ROOT/path).read_text(encoding='utf-8');data=IMPACT.manifest()
    for old,new in [(before,broken),(before,None),(None,before)]:
        assert IMPACT.root_owners(path,old,new,data)[0] == set(IMPACT.CONTRACT.DISTRIBUTIONS)


def test_rename_status_preserves_both_ownership_boundaries():
    raw='R100\0packages/radia-simulink/old.py\0packages/radia-optuna/new.py\0D\0Build.ps1\0'
    assert IMPACT.changed_paths(raw)=={'packages/radia-simulink/old.py','packages/radia-optuna/new.py','Build.ps1'}
    with pytest.raises(ValueError): IMPACT.changed_paths('R100\0old\0')
    data=IMPACT.manifest()
    for lane in ('radia-simulink','radia-optuna'):
        assert IMPACT.lane_required(lane,IMPACT.changed_paths(raw),lambda p:(None,None),data)[0]


def test_directive_order_and_quoted_hash_are_not_comment_only_edits():
    data=IMPACT.manifest()
    a='include packages/radia-simulink/source.m\nexclude packages/radia-simulink/source.m\n'
    b='exclude packages/radia-simulink/source.m\ninclude packages/radia-simulink/source.m\n'
    for path, before, after in (('MANIFEST.in', a, b),
                                ('CMakeLists.txt', 'option(RADIA_BUILD_MATLAB_MEX "a#b" OFF)',
                                 'option(RADIA_BUILD_MATLAB_MEX "a#c" OFF)')):
        assert IMPACT.root_owners(path,before,after,data)[0] == set(IMPACT.CONTRACT.owner_rule(data,path)['owners'])


def test_ams_legacy_exemption_cannot_override_unparsable_or_unknown_content():
    before=manifest();after=manifest();after['impact_rules']='invalid'
    assert IMPACT.lane_required('sparsesolv',{'tests/test_tier_manifest.json'},
        lambda p:(json.dumps(before),json.dumps(after)),IMPACT.manifest())[0]


@pytest.mark.parametrize('lane', ['radia-simulink','radia-optuna','sparsesolv','radia-mcp-matrix'])
def test_missing_base_or_blob_lookup_cannot_skip_a_lane(lane):
    data=IMPACT.manifest()
    assert IMPACT.required_for_event(lane,'push',{'before':'0'*40},None,data)[0]
    assert IMPACT.required_for_event(lane,'pull_request',{},None,data)[0]
    def missing(path): raise OSError('blob unavailable')
    assert IMPACT.lane_required(lane,{'Build.ps1'},missing,data)[0]


@pytest.mark.parametrize('event', ['push','pull_request'])
def test_router_uses_full_push_range_or_pr_merge_base(event):
    calls=[];base='a'*40
    payload={'before':base,'pull_request':{'base':{'sha':base}}}
    def git(*args):
        calls.append(args)
        if args[0]=='merge-base': return 'b'*40
        if args[0]=='diff': return 'M\0pyproject.toml\0'
        return '[project]\nversion="'+('1' if args[1].startswith(('a','b')) else '2')+'"\n'
    assert not IMPACT.required_for_event('radia-simulink',event,payload,git,IMPACT.manifest())[0]
    diff=next(c for c in calls if c[0]=='diff')
    assert diff[4]==('b'*40 if event=='pull_request' else base)


def _router_lane(command):
    """Exact --lane argument of the router command (not a substring match)."""
    tokens=command.split()
    values=[tokens[i+1] for i,t in enumerate(tokens[:-1]) if t=='--lane']
    return values[0] if len(values)==1 else None


def test_router_lane_is_matched_exactly():
    assert _router_lane('python tools/sparsesolv_ci_impact.py --lane radia-simulink')=='radia-simulink'
    assert _router_lane('python tools/sparsesolv_ci_impact.py --lane radia-simulink-release')!='radia-simulink'
    assert _router_lane('python tools/sparsesolv_ci_impact.py --lane a --lane b') is None


def test_all_distribution_gates_are_hosted_and_preserve_full_manual_release_runs():
    data=IMPACT.manifest()
    for lane,definition in data['lanes'].items():
        if definition['scope']=='repository_contract' or not definition['content_diff_active']: continue
        workflow=yaml.load((ROOT/definition['workflow']).read_text(encoding='utf-8'),Loader=yaml.BaseLoader)
        if 'ci-impact' not in workflow['jobs']: continue  # tag/manual-only workflows
        gate=workflow['jobs']['ci-impact']
        assert gate['runs-on']=='ubuntu-latest'
        assert gate['steps'][0]['with']['fetch-depth']=='0'
        # The router step must run; lane-local steps (e.g. eqnedit64's prose gate) may follow it.
        routers=[s.get('run','') for s in gate['steps'] if 'tools/sparsesolv_ci_impact.py' in s.get('run','')]
        assert len(routers)==1 and _router_lane(routers[0])==lane, (lane, routers)
        for name,job in workflow['jobs'].items():
            if name in ('ci-impact', 'nightly-failure'):continue
            assert 'ci-impact' in job['needs']
            assert "needs.ci-impact.outputs.required == 'true'" in job['if']
        assert IMPACT.required_for_event(lane,'workflow_dispatch',{},None,data)[0]
        assert IMPACT.required_for_event(lane,'push',{'ref':'refs/tags/package-v1'},None,data)[0]


@pytest.mark.parametrize("bad", [True, False, float("nan"), float("inf"), 0, -1, "30"])
def test_invalid_tier_budget_fails_closed(bad):
    """JSON true equals 1 in Python; a budget edited to true must not read as unchanged."""
    before = (ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8")
    data = json.loads(before)
    name = next(n for n, p in data["profiles"].items() if "max_elapsed_seconds" in p)
    data["profiles"][name]["max_elapsed_seconds"] = bad
    after = json.dumps(data)
    with pytest.raises(ValueError):
        IMPACT._tier_owners(before, after, IMPACT.CONTRACT.load(ROOT),
                            IMPACT.CONTRACT.load(ROOT)["rules"][0].get("content_diff", {}))
    if name == "sparsesolv":
        with pytest.raises(ValueError):
            IMPACT.ams_tier_projection_unchanged(before, after)


def test_numeric_one_and_true_are_distinct_changes():
    assert IMPACT._canonical({"max_elapsed_seconds": 1}) != IMPACT._canonical({"max_elapsed_seconds": True})


def test_toml_number_to_bool_is_a_change():
    before = "[tool.pytest.ini_options]\nverbosity_assertions = 1\n"
    after = "[tool.pytest.ini_options]\nverbosity_assertions = true\n"
    assert not IMPACT.test_dependencies_only(before, after)
    assert IMPACT._canonical(1) != IMPACT._canonical(True)


def test_huge_integer_budget_is_invalid_not_a_crash():
    assert IMPACT._valid_budget(10 ** 400) is False


def test_typed_canonical_separates_dates_from_look_alike_strings():
    import datetime
    when = datetime.date(2026, 10, 15)
    for look_alike in (f"date:{when!r}", "2026-10-15", repr(when)):
        assert IMPACT._canonical({"k": when}) != IMPACT._canonical({"k": look_alike})
    before = "[tool.pytest.ini_options]\ncustom_date = 2026-10-15\n"
    after = '[tool.pytest.ini_options]\ncustom_date = "date:datetime.date(2026, 10, 15)"\n'
    assert not IMPACT.test_dependencies_only(before, after)


def test_malformed_event_payload_fails_closed():
    assert IMPACT.required_for_event("radia-simulink", "push", [], None, {})[0] is True


@pytest.mark.parametrize("value", [None, "tests/x.py", [1], {"a": []}])
def test_malformed_impact_rule_fails_closed(value):
    before = (ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8")
    data = json.loads(before)
    data.setdefault("impact_rules", {})["src/radia/new_router_probe.py"] = value
    with pytest.raises(ValueError):
        IMPACT._tier_owners(before, json.dumps(data), IMPACT.CONTRACT.load(ROOT), {})


def test_added_and_removed_impact_rules_select_owners():
    before = (ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8")
    data = json.loads(before)
    data.setdefault("impact_rules", {})["src/radia/new_router_probe.py"] = ["tests/test_sparsesolv_ci_contract.py"]
    added = json.dumps(data)
    manifest = IMPACT.CONTRACT.load(ROOT)
    assert IMPACT._tier_owners(before, added, manifest, {})
    assert IMPACT._tier_owners(added, before, manifest, {})


def test_powershell_block_comment_is_not_a_cmake_comment():
    import json as _json
    manifest = IMPACT.CONTRACT.load(ROOT)
    old = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    new = old + "\n<# INVALID CMAKE #>\n"
    required = IMPACT.lane_required("radia-simulink", {"CMakeLists.txt"}, lambda p: (old, new), manifest)
    assert required[0] is True
    assert IMPACT.code_tokens(old, "cmake") != IMPACT.code_tokens(new, "cmake")


@pytest.mark.parametrize("profile", [{}, {"description": "only prose"}])
def test_added_or_removed_unknown_profile_fails_closed(profile):
    before = (ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8")
    data = json.loads(before)
    data["profiles"]["unknown_new_profile"] = profile
    after = json.dumps(data)
    manifest = IMPACT.CONTRACT.load(ROOT)
    for old, new in ((before, after), (after, before)):
        with pytest.raises(ValueError):
            IMPACT._tier_owners(old, new, manifest, {})


def test_joining_powershell_statements_is_a_change():
    manifest = IMPACT.CONTRACT.load(ROOT)
    old = (ROOT / "Build.ps1").read_text(encoding="utf-8")
    marker = '$ErrorActionPreference = "Stop"'
    assert marker in old
    head, tail = old.split(marker, 1)
    joined = head + marker + " " + tail.lstrip("\r\n")
    assert IMPACT.lane_required("radia-simulink", {"Build.ps1"}, lambda p: (old, joined), manifest)[0] is True
    # Indentation and comment-only edits stay exempt.
    reindented = old.replace(marker, "    " + marker + "  # comment", 1)
    assert IMPACT.code_tokens(old, "powershell") == IMPACT.code_tokens(reindented, "powershell")


def test_powershell_formatting_change_selects_all_owners():
    """Trailing space after a continuation backtick breaks PowerShell; never exempt it."""
    manifest = IMPACT.CONTRACT.load(ROOT)
    old = (ROOT / "Build.ps1").read_text(encoding="utf-8")
    line = next(l for l in old.splitlines() if l.rstrip().endswith("`"))
    new = old.replace(line, line.rstrip() + " ", 1)
    assert new != old
    for lane in ("radia-simulink", "sparsesolv"):
        assert IMPACT.lane_required(lane, {"Build.ps1"}, lambda p: (old, new), manifest)[0] is True


def test_cmake_variable_reference_split_selects_owners():
    manifest = IMPACT.CONTRACT.load(ROOT)
    old = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    assert "${Python3_EXECUTABLE}" in old
    new = old.replace("${Python3_EXECUTABLE}", "$ {Python3_EXECUTABLE}", 1)
    for lane in ("radia-simulink", "sparsesolv"):
        assert IMPACT.lane_required(lane, {"CMakeLists.txt"}, lambda p: (old, new), manifest)[0] is True


def test_attributed_edit_cannot_mask_a_global_cmake_edit():
    manifest = IMPACT.CONTRACT.load(ROOT)
    old = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    new = old.replace("${Python3_EXECUTABLE}", "$ {Python3_EXECUTABLE}", 1).replace("retired", "unsupported", 1)
    assert IMPACT.lane_required("sparsesolv", {"CMakeLists.txt"}, lambda p: (old, new), manifest)[0] is True


def test_separated_acceptance_has_an_executing_lane_and_input_triggers():
    """Moving a numerical gate must preserve execution, not just ownership labels."""
    data = IMPACT.manifest()
    tier = json.loads((ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8"))["profiles"]
    homes = {
        "sparsesolv": [
            "tests/test_ih_operator_assembly.py", "tests/test_ih_operator_axisym.py",
            "tests/test_ih_fin_surface_boundary.py", "tests/test_ih_panel_esim_assembly.py",
            "validation_test/induction_heating/test_ih_operator_physics_golden.py",
            "validation_test/radia_mcp/test_aphi_gauge_point.py",
            "validation_test/radia_mcp/test_axi_eddy_current_drive.py",
            "validation_test/radia_mcp/test_eddy_driven_region.py",
            "validation_test/radia_mcp/test_hoibc_weingarten_curvature.py"],
        "radia-simulink": ["tests/cpp/test_radia_ih_thermal_acceptance.cpp",
                      "tests/cpp/test_radia_mex_sparse_symmetry.cpp"],
    }
    for lane, paths in homes.items():
        definition = data["lanes"][lane]
        workflow = yaml.load((ROOT / definition["workflow"]).read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
        scripts = "\n".join(step.get("run", "") for job in workflow["jobs"].values()
                            for step in job.get("steps", []))
        for path in paths:
            assert all(IMPACT.CONTRACT.selected(path, definition["events"][e])
                       for e in ("push", "pull_request")), (lane, path)
            assert set(IMPACT.CONTRACT.owner_rule(data, path)["owners"]) & set(definition["owners"])
            if lane == "sparsesolv":
                profiles = [name for name in ("solver-numerics", "solver-heavy")
                            if path in tier[name]["paths"]]
                assert profiles
                for name in profiles:
                    if name == "solver-heavy":
                        assert "tools/validation_source_gate.py" in scripts
                        gate = (ROOT / "tools/validation_source_gate.py").read_text(encoding="utf-8")
                        assert "PROFILE = 'solver-heavy'" in gate
                        assert "'tools/run_test_tier.py','--profile',PROFILE" in gate
                    else:
                        assert "--profile " + name in scripts
            elif lane == "radia-simulink":
                project = (ROOT / "packages/radia-simulink/tests/native/CMakeLists.txt").read_text(encoding="utf-8")
                assert path in project
                assert Path(path).stem in scripts
                assert "ctest --test-dir" in scripts and "--no-tests=error" in scripts
                assert Path(path).stem not in (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
                assert "target_compile_options(test_radia_mex_sparse_symmetry PRIVATE /W4 /permissive-)" in project
                assert "target_compile_options(test_radia_mex_sparse_symmetry PRIVATE -Wall -Wextra -Wpedantic)" in project
                assert '-G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=cl' in scripts
                assert 'vcvars64.bat' in scripts and 'call `"$vcvars`" && cmake' in scripts
                assert '-G "Visual Studio' not in scripts
    inputs = {
        "sparsesolv": ["packages/cubit-mesh-export/src/cubit_mesh_export/check.py",
                       "packages/radia-simulink/src/radia_simulink/simulink/ih_operator_assembly.py",
                       "packages/radia-mcp/src/radia_mcp/radia_ngsolve/solve.py",
                       "packages/radia-mcp/src/radia_mcp/radia_ngsolve/_direct.py",
                       "packages/radia-mcp/src/radia_mcp/radia_ngsolve/_vendor/residual_gate.py",
                       "packages/radia-mcp/src/radia_mcp/peec/hoibc_knowledge.py"],
        "radia-simulink": ["src/runtime/ih/radia_ih_thermal.cpp", "src/runtime/ih/radia_ih_runtime.cpp",
                      "src/runtime/ih/radia_ih_transport.cpp", "src/runtime/ih/radia_ih_thermal.h",
                      "src/core/radia_mex_sparse_symmetry.h"],
    }
    for lane, paths in inputs.items():
        for path in paths:
            assert all(IMPACT.CONTRACT.selected(path, data["lanes"][lane]["events"][e])
                       for e in ("push", "pull_request")), (lane, path)
    mcp_workflow = yaml.load((ROOT / ".github/workflows/radia-mcp-matrix.yml").read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    assert "numerical-validation" not in mcp_workflow["jobs"]
    # #421 runs the router contract once in the always-on Ubuntu policy job.
    spec = importlib.util.spec_from_file_location("tier_coverage", ROOT / "tools/run_test_tier.py")
    tier_router = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tier_router)
    router_test = "tests/test_sparsesolv_ci_contract.py"
    selected = tier_router.select_impact_tests(tier["fast-contracts"]["paths"], [router_test],
                                               profile_name="fast-contracts")
    assert router_test not in selected
    assert router_test not in tier["repo-contracts"]["paths"]
    assert all(IMPACT.CONTRACT.selected(router_test, data["lanes"]["policy-lint"]["events"][event])
               for event in ("push", "pull_request"))


def test_data_separation_keeps_existing_physics_route_coverage_check(monkeypatch, tmp_path):
    """Reuse the existing route contract; do not invent a second gate table."""
    data = json.loads((ROOT / "tests/test_tier_manifest.json").read_text(encoding="utf-8"))
    data["profiles"]["solver-numerics"]["paths"].remove("tests/test_p1_linear.py")
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests/test_tier_manifest.json").write_text(json.dumps(data), encoding="utf-8")
    monkeypatch.setitem(globals(), "ROOT", tmp_path)
    with pytest.raises(AssertionError):
        test_p1_linear_routes_to_solver_numerics_and_native_ci("src/radia/p1_linear.py")
