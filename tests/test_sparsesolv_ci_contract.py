"""Keep native solver coverage connected without importing native dependencies."""

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
    profiles = json.loads((ROOT / "tests/test_tier_manifest.json").read_text())["profiles"]
    assert SUITE in profiles["native-smoke"]["paths"]
    parallel = "tests/test_ams_internal_setup_parallel.py"
    additional = ["src/ext/sparsesolv/tests/test_hermitian_iccg.py",
                  "src/ext/sparsesolv/tests/test_ams_repeatability.py"]
    assert profiles["sparsesolv"]["paths"] == [SUITE, *additional, parallel]
    assert all(path in profiles["native-smoke"]["paths"] for path in additional)
    assert parallel in profiles["native-smoke"]["paths"]
    assert SUITE not in profiles["fast-contracts"]["paths"]


def test_solver_workflow_builds_before_testing_on_mdx():
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
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
    source = (ROOT / SUITE).read_text()
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
    assert output.read_text().startswith("required=true\nreason=")
    assert "fail closed" in output.read_text()


def test_native_impact_gates_expensive_steps_and_watches_itself():
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
    triggers = workflow.get("on", workflow.get(True))
    for event in ("push", "pull_request"):
        assert "tools/sparsesolv_ci_impact.py" in triggers[event]["paths"]
        assert set(triggers[event]["paths"]) == set(IMPACT.NATIVE_PATHS) | {
            "pyproject.toml", "tests/test_tier_manifest.json", "!**/*.md"}
        assert triggers[event]["paths"][-1] == "!**/*.md"
    steps = workflow["jobs"]["ams-regression"]["steps"]
    selector = next(s for s in steps if s.get("id") == "native-impact")
    assert "python tools/sparsesolv_ci_impact.py" in selector["run"]
    for step in steps:
        if step.get("name", "").startswith(("Create isolated", "Build ", "Verify source", "Upload solver")):
            assert "steps.native-impact.outputs.required == 'true'" in step["if"]


def manifest():
    return {"schema": "radia.test-tier-manifest.v1", "impact_rules": {}, "profiles": {
        "sparsesolv": {"extends": "base", "paths": [SUITE], "max_elapsed_seconds": 30},
        "base": {"paths": ["tests/base.py"]}, "fast": {"paths": ["tests/old.py"]}}}


@pytest.mark.parametrize("change,expected", [
    (lambda d: d["impact_rules"].update({"tools/new.py": ["tests/new.py"]}), True),
    (lambda d: d["profiles"]["fast"]["paths"].append("tests/new.py"), True),
    (lambda d: d["profiles"]["sparsesolv"]["paths"].append("tests/new.py"), False),
    (lambda d: d["profiles"]["base"]["paths"].append("tests/new.py"), False),
    (lambda d: d["profiles"]["sparsesolv"].update(max_elapsed_seconds=60), False),
    (lambda d: d["profiles"]["sparsesolv"].update(extends="fast"), False),
    (lambda d: d.update(future_global_option=True), False),
])
def test_manifest_compares_ams_and_all_ancestors(change, expected):
    before, after = manifest(), manifest()
    change(after)
    assert IMPACT.ams_manifest_unchanged(json.dumps(before), json.dumps(after)) is expected


@pytest.mark.parametrize("change", [
    lambda d: d.update(schema="unknown"),
    lambda d: d["profiles"].pop("sparsesolv"),
    lambda d: d["profiles"].pop("base"),
    lambda d: d["profiles"]["base"].update(extends="sparsesolv"),
    lambda d: d["profiles"]["base"].update(paths="bad"),
    lambda d: d["profiles"]["base"].update(paths=[SUITE]),
    lambda d: d["profiles"]["sparsesolv"].update(max_elapsed_seconds=-1),
])
def test_invalid_manifest_fails_closed(change):
    before, after = manifest(), manifest()
    change(after)
    with pytest.raises((KeyError, ValueError, TypeError)):
        IMPACT.ams_manifest_unchanged(json.dumps(before), json.dumps(after))


@pytest.mark.parametrize("extra,required", [
    ("tests/new.py\0tools/new.py\0", False),
    ("src/ext/sparsesolv/test.cpp\0", True),
    ("tests/matlab/test_beam_transfer_mex.m\0", False),
    ("src/core/rad_beam_transfer.cpp\0", True),
    ("packages/radia-gpu/matlab/+radia/+internal/callMex.m\0", False),
    ("tools/run_test_tier.py\0", True),
    (".github/workflows/sparsesolv.yml\0", True),
])
def test_manifest_edit_with_other_files(extra, required):
    path = "tests/test_tier_manifest.json"
    before, after = manifest(), manifest()
    after["impact_rules"]["tools/new.py"] = ["tests/new.py"]
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
    "packages/radia-gpu/src/radia_gpu/simulink/ih_operator_assembly.py",
    "src/radia/coil_from_cad.py",
    "packages/cubit-mesh-export/src/cubit_mesh_export/check.py",
])
def test_numerical_source_changes_select_built_ngsolve_lane(source):
    assert IMPACT.native_required("push", {"before": "a" * 40},
                                  lambda *args: source + "\0")[0]
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
    steps = workflow["jobs"]["ams-regression"]["steps"]
    numerical = next(s for s in steps if s.get("name") == "Run required NGSolve numerical contracts")
    assert "--profile solver-numerics" in numerical["run"]
    assert numerical["env"]["RADIA_TESTS_ALLOW_PARTIAL"] == "0"
    assert "find('skipped') is None" in numerical["run"]


@pytest.mark.parametrize("source", [
    "src/core/radia_mex_sparse_symmetry.h",
    "tests/cpp/test_radia_mex_sparse_symmetry.cpp",
])
def test_mex_sparse_symmetry_sources_require_the_native_lane(source):
    assert IMPACT.native_required("push", {"before": "a" * 40},
                                  lambda *args: source + "\0")[0]


def test_ih_thermal_and_complex_ams_contracts_run_with_ngsolve():
    profiles = json.loads((ROOT / "tests/test_tier_manifest.json").read_text())["profiles"]
    required = {"tests/test_ih_axisym_coupled.py", "tests/test_ih_thermal_transfer.py",
                "tests/test_ih_thermal_exposure.py",
                "tests/test_ih_operator_assembly.py",
                "tests/test_ih_operator_axisym.py",
                "tests/test_ih_fin_surface_boundary.py",
                "validation_test/induction_heating/test_fem_ams_parity.py"}
    assert required <= set(profiles["solver-numerics"]["paths"])
    assert required.isdisjoint(profiles["fast-contracts"]["paths"])


def test_numerical_tier_changes_cannot_skip_native_lane():
    before, after = manifest(), manifest()
    after["profiles"]["solver-numerics"] = {"paths": ["tests/test_mixed_omega_newton.py"]}
    assert not IMPACT.ams_manifest_unchanged(json.dumps(before), json.dumps(after))


def test_heavy_derivative_validation_remains_required_after_tier_split():
    path = "validation_test/topology_optimization/test_hex_cluster_derivative.py"
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
    steps = workflow["jobs"]["ams-regression"]["steps"]
    heavy = next(s for s in steps if s.get("name") == "Run required heavier numerical validation")
    assert "--profile solver-heavy" in heavy["run"]
    assert heavy["env"]["RADIA_TESTS_ALLOW_PARTIAL"] == "0"
    assert "find('skipped') is None" in heavy["run"]
    assert "steps.native-impact.outputs.required == 'true'" in heavy["if"]
    profiles = json.loads((ROOT / "tests/test_tier_manifest.json").read_text())["profiles"]
    assert path in profiles["solver-heavy"]["paths"]
    assert path not in profiles["solver-numerics"]["paths"]
    assert IMPACT.native_required("push", {"before": "a" * 40},
                                  lambda *args: path + "\0")[0]
    before, after = manifest(), manifest()
    after["profiles"]["solver-heavy"] = {"paths": [path]}
    assert not IMPACT.ams_manifest_unchanged(json.dumps(before), json.dumps(after))


def test_mex_symmetry_cpp_runs_through_generator_independent_ctest():
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
    steps = workflow["jobs"]["ams-regression"]["steps"]
    build = next(s for s in steps if s.get("name") == "Build current package and SparseSolv")
    assert "ctest --test-dir build-msvc -C Release" in build["run"]
    assert "--no-tests=error" in build["run"]
    assert "radia\\.mex_sparse_symmetry\\.cpp$" in build["run"]
    assert "build-msvc\\Release\\test_radia_mex_sparse_symmetry.exe" not in build["run"]


def test_native_solver_lane_builds_required_axisymmetric_extension():
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
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
    manifest = json.loads((ROOT / "tests/test_tier_manifest.json").read_text())
    test = "tests/test_p1_linear.py"
    assert test in manifest["impact_rules"]["src/radia/p1_linear.py"]
    assert test in manifest["profiles"]["solver-numerics"]["paths"]
    assert changed in IMPACT.NATIVE_PATHS
    assert IMPACT.native_required(
        "push", {"before": "a" * 40}, lambda *args: changed + "\0"
    )[0]
    workflow = yaml.safe_load((ROOT / ".github/workflows/sparsesolv.yml").read_text())
    triggers = workflow.get("on", workflow.get(True))
    for event in ("push", "pull_request"):
        assert changed in triggers[event]["paths"]


ROOT_CASES = [
    ('pyproject.toml', '[project]\nversion="1"\n', '[project]\nversion="2"\n', {'radia'}),
    ('pyproject.toml', '[tool.setuptools.package-data]\nradia_gpu=["*.m"]\n',
     '[tool.setuptools.package-data]\nradia_gpu=["*.cu"]\n', {'radia-gpu'}),
    ('Build.ps1', "param()\nif ($OptunaMexOnly) { matlab -batch 'old' }\n",
     "param()\nif ($OptunaMexOnly) { matlab -batch 'new' }\n", {'radia-optuna'}),
    ('Build.ps1', 'if ($RadiaOnly) { python old }', 'if ($RadiaOnly) { python new }', {'radia'}),
    ('Build.ps1', 'if ($AmsCudaOnly) { python old }', 'if ($AmsCudaOnly) { python new }', {'radia-gpu'}),
    ('CMakeLists.txt', 'option(RADIA_BUILD_MATLAB_MEX "enable" OFF)',
     'option(RADIA_BUILD_MATLAB_MEX "enable" ON)', {'radia-gpu'}),
    ('CMakeLists.txt', 'target_compile_definitions(_radia_pybind PRIVATE OLD)',
     'target_compile_definitions(_radia_pybind PRIVATE NEW)', {'radia'}),
    ('MANIFEST.in', 'recursive-include packages/radia-gpu *.m\n',
     'recursive-include packages/radia-gpu *.cu\n', {'radia-gpu'}),
    ('tests/test_tier_manifest.json',
     json.dumps({'schema':'radia.test-tier-manifest.v1','profiles':{'gpu-unit':{
         'paths':['packages/radia-gpu/tests/test_ic_payload.py'],'max_elapsed_seconds':1}}}),
     json.dumps({'schema':'radia.test-tier-manifest.v1','profiles':{'gpu-unit':{
         'paths':['packages/radia-gpu/tests/test_ic_payload.py'],'max_elapsed_seconds':2}}}), {'radia-gpu'}),
]


@pytest.mark.parametrize('path,before,after,owners', ROOT_CASES)
def test_root_content_diff_selects_only_evidenced_owners(path,before,after,owners):
    assert IMPACT.root_owners(path,before,after,IMPACT.manifest())[0] == owners


@pytest.mark.parametrize('path', IMPACT.CONTRACT.ROOT_CONFIGS)
def test_real_root_comments_and_whitespace_select_no_build(path):
    before=(ROOT/path).read_text(encoding='utf-8')
    after=before+'\n' if path.endswith('.json') else before+'\n# comment only\n'
    assert IMPACT.root_owners(path,before,after,IMPACT.manifest())[0] == set()


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
    raw='R100\0packages/radia-gpu/old.py\0packages/radia-optuna/new.py\0D\0Build.ps1\0'
    assert IMPACT.changed_paths(raw)=={'packages/radia-gpu/old.py','packages/radia-optuna/new.py','Build.ps1'}
    with pytest.raises(ValueError): IMPACT.changed_paths('R100\0old\0')
    data=IMPACT.manifest()
    for lane in ('radia-gpu','radia-optuna'):
        assert IMPACT.lane_required(lane,IMPACT.changed_paths(raw),lambda p:(None,None),data)[0]


def test_directive_order_and_quoted_hash_are_not_comment_only_edits():
    data=IMPACT.manifest()
    a='include packages/radia-gpu/source.m\nexclude packages/radia-gpu/source.m\n'
    b='exclude packages/radia-gpu/source.m\ninclude packages/radia-gpu/source.m\n'
    assert IMPACT.root_owners('MANIFEST.in',a,b,data)[0]=={'radia-gpu'}
    assert IMPACT.root_owners('CMakeLists.txt',
        'option(RADIA_BUILD_MATLAB_MEX "a#b" OFF)',
        'option(RADIA_BUILD_MATLAB_MEX "a#c" OFF)',data)[0]=={'radia-gpu'}


def test_ams_legacy_exemption_cannot_override_unparsable_or_unknown_content():
    before=manifest();after=manifest();after['impact_rules']='invalid'
    assert IMPACT.lane_required('sparsesolv',{'tests/test_tier_manifest.json'},
        lambda p:(json.dumps(before),json.dumps(after)),IMPACT.manifest())[0]


@pytest.mark.parametrize('lane', ['radia-gpu','radia-optuna','sparsesolv','radia-mcp-matrix'])
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
    assert not IMPACT.required_for_event('radia-gpu',event,payload,git,IMPACT.manifest())[0]
    diff=next(c for c in calls if c[0]=='diff')
    assert diff[4]==('b'*40 if event=='pull_request' else base)


def test_all_distribution_gates_are_hosted_and_preserve_full_manual_release_runs():
    data=IMPACT.manifest()
    for lane,definition in data['lanes'].items():
        if definition['scope']=='repository_contract' or not definition['content_diff_active']: continue
        workflow=yaml.load((ROOT/definition['workflow']).read_text(encoding='utf-8'),Loader=yaml.BaseLoader)
        if 'ci-impact' not in workflow['jobs']: continue  # tag/manual-only workflows
        gate=workflow['jobs']['ci-impact']
        assert gate['runs-on']=='ubuntu-latest'
        assert gate['steps'][0]['with']['fetch-depth']=='0'
        assert f'--lane {lane}' in gate['steps'][-1]['run']
        for name,job in workflow['jobs'].items():
            if name=='ci-impact':continue
            assert 'ci-impact' in job['needs']
            assert "needs.ci-impact.outputs.required == 'true'" in job['if']
        assert IMPACT.required_for_event(lane,'workflow_dispatch',{},None,data)[0]
        assert IMPACT.required_for_event(lane,'push',{'ref':'refs/tags/package-v1'},None,data)[0]
