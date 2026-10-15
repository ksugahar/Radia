"""A new workflow must explicitly declare whether it belongs to solver CI."""
from fnmatch import fnmatchcase
import copy
import importlib.util
from pathlib import Path
import re

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("ci_distribution_contract", ROOT / "tools/ci_distribution_contract.py")
CONTRACT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONTRACT)
MANIFEST = CONTRACT.load(ROOT)
SOLVER = {"radia-fast.yml", "build-test.yml", "release.yml"}
INDEPENDENT = {Path(lane["workflow"]).name for lane in MANIFEST["lanes"].values()} - SOLVER


def test_every_workflow_has_an_explicit_owner():
    CONTRACT.validate_schema(MANIFEST)
    names = {p.name for p in (ROOT / ".github/workflows").glob("*.yml")}
    assert names == {Path(lane["workflow"]).name for lane in MANIFEST["lanes"].values()}
    for name, lane in MANIFEST["lanes"].items():
        document = yaml.load((ROOT / lane["workflow"]).read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
        assert CONTRACT.workflow_events(document) == lane["events"], name
        assert lane["scope"] == ("repository_contract" if name == "policy-lint" else "distribution")
        assert CONTRACT.owner_rule(MANIFEST, lane["workflow"])["owners"] == lane["owners"]
        for reference in lane["selector_refs"]:
            assert (ROOT / reference.split(":", 1)[0]).is_file(), reference


def _files():
    # These two authored files are not in git ls-files until the reviewer stages
    # PR1. Existing tracked paths are always read from Git, never an inventory.
    return sorted(set(CONTRACT.tracked(ROOT)) | {
        "tools/ci_distribution_contract.py", "tools/ci_distribution_manifest.json"})


def test_tracked_coverage_and_target_contract_with_explicit_debt():
    files = _files()
    for path in files:
        CONTRACT.owner_rule(MANIFEST, path)
    violations = CONTRACT.issues(MANIFEST, files)
    CONTRACT.exception_coverage(MANIFEST, violations)
    CONTRACT.check_shrink_only(MANIFEST["migration_exceptions"], EXCEPTION_CEILING)


def test_exception_removal_is_allowed_but_growth_and_stale_entries_fail():
    exceptions = MANIFEST["migration_exceptions"]
    CONTRACT.check_shrink_only(exceptions[1:], EXCEPTION_CEILING)
    expanded = copy.deepcopy(exceptions)
    expanded[0]["pattern"] = "**"
    with pytest.raises(ValueError, match="new or expanded"):
        CONTRACT.check_shrink_only(expanded, EXCEPTION_CEILING)
    new = copy.deepcopy(exceptions[0]); new["id"] = "UNAPPROVED"
    with pytest.raises(ValueError, match="new or expanded"):
        CONTRACT.check_shrink_only([new], EXCEPTION_CEILING)
    with pytest.raises(ValueError, match="resolved/unused"):
        CONTRACT.exception_coverage(MANIFEST, set())


@pytest.mark.parametrize("change", ["owner", "reason", "selection", "root_catchall", "duplicate", "empty_attribution", "unknown_attribution"])
def test_schema_rejects_ambiguous_or_unjustified_ownership(change):
    candidate = copy.deepcopy(MANIFEST)
    shared = next(r for r in candidate["rules"] if len(r["owners"]) > 1)
    if change == "owner":
        candidate["rules"][0]["owners"] = ["unknown"]
    elif change == "reason":
        del shared["reason"]
    elif change == "selection":
        shared["selection"] = "implicit_fallback"
    elif change == "root_catchall":
        candidate["rules"].append({"pattern": "**", "owners": ["repo"]})
    elif change in ("empty_attribution", "unknown_attribution"):
        candidate["root_content_diff"]["attribution"]["Build.ps1"]["OptunaMexOnly"] = [] if change == "empty_attribution" else ["unknown"]
    else:
        candidate["rules"].append(copy.deepcopy(candidate["rules"][0]))
    with pytest.raises(ValueError):
        CONTRACT.validate_schema(candidate)


def test_ordered_matching_preserves_actions_markdown_reinclusion():
    projection = {"paths": ["packages/radia-gpu/**", "!**/*.md", "packages/radia-gpu/README.md"]}
    assert CONTRACT.selected("packages/radia-gpu/native/solver.cu", projection)
    assert CONTRACT.selected("packages/radia-gpu/README.md", projection)
    assert not CONTRACT.selected("packages/radia-gpu/handover.md", projection)
    assert not CONTRACT.matches("packages/a/b.py", "packages/*.py")
    assert CONTRACT.matches("README.md", "**/*.md")
    assert CONTRACT.matches("docs/ltspice/README.md", "docs/**/*.md")
    with pytest.raises(ValueError, match="unowned"):
        CONTRACT.owner_rule(MANIFEST, "unclassified-root-file")
    assert CONTRACT.owner_rule(MANIFEST, "packages/radia-mcp/src/radia_mcp/presentation/tools.py")["owners"] == ["radia-mcp"]


def test_shared_ams_closure_and_native_consumers_are_real():
    closure, = MANIFEST["shared_include_closures"]
    actual = CONTRACT.include_closure(ROOT, closure["seeds"], closure["prefix"])
    assert actual == closure["paths"]
    build = (ROOT / "packages/radia-gpu/tools/build_ams_cuda.py").read_text(encoding="utf-8")
    for seed in closure["seeds"]:
        assert seed in build
    shared = {r["pattern"] for r in MANIFEST["rules"] if len(r["owners"]) > 1}
    assert {p for p in shared if p.startswith("src/ext/sparsesolv/")} == set(actual)
    assert CONTRACT.owner_rule(MANIFEST, "src/ext/sparsesolv/include/sparsesolv/preconditioners/ic_preconditioner.hpp")["owners"] == ["radia"]
    helper = "tools/native_build_provenance.ps1"
    assert helper in (ROOT / "packages/radia-gpu/matlab/verify_radia_mex_provenance.py").read_text(encoding="utf-8")
    assert "$NativeBuildProvenanceScript" in (ROOT / "Build.ps1").read_text(encoding="utf-8")
    assert CONTRACT.owner_rule(MANIFEST, helper)["owners"] == ["radia", "radia-gpu"]
    # Derive the shared core closure from real GPU quote includes and the build
    # receipt, rather than marking every CPU kernel as GPU-owned.
    seeds = set()
    for source in (ROOT / "packages/radia-gpu/native").iterdir():
        if source.suffix not in (".cu", ".cuh", ".h", ".hpp", ".cpp", ".c"):
            continue
        for include in re.findall(r'^\s*#\s*include\s*"([^"]+)"', source.read_text(encoding="utf-8"), re.M):
            if (ROOT / "src/core" / include).is_file():
                seeds.add("src/core/" + include)
    hmatrix_build = (ROOT / "packages/radia-gpu/tools/build_hmatrix_cuda.py").read_text(encoding="utf-8")
    # Every GPU build script that names or compiles a src/core file seeds the closure.
    for script in sorted((ROOT / "packages/radia-gpu/tools").glob("build_*_cuda.py")):
        seeds.update("src/core/" + name for name in
                     re.findall(r"src/core/([A-Za-z0-9_.]+)", script.read_text(encoding="utf-8")))
    core = set(CONTRACT.include_closure(ROOT, sorted(seeds), "src/core"))
    assert {p for p in shared if p.startswith("src/core/")} == core
    assert len(core) == 14
    # The radia-gpu lane must start on every file its native build reads.
    lane = MANIFEST["lanes"]["radia-gpu"]["events"]
    for event in ("push", "pull_request"):
        missing = [p for p in sorted(core) if not CONTRACT.selected(p, lane[event])]
        assert not missing, (event, missing)
    assert "src/ext/HACApK" in hmatrix_build
    expected = set(CONTRACT.ROOT_CONFIGS) | set(actual) | core | {
        "src/ext/HACApK/**", helper, "tests/test_native_build_provenance.py",
        "tests/matlab/fixtures/ltspice_rc.cir", "tools/sparsesolv_ci_impact.py",
        "tools/ci_distribution_contract.py", "tools/ci_distribution_manifest.json",
        "packages/radia-mcp/src/radia_mcp/ltspice/**",
        "packages/radia-mcp/src/radia_mcp/common/__init__.py",
        "packages/radia-mcp/src/radia_mcp/common/lazy_call.py",
        "packages/radia-mcp/src/radia_mcp/common/status.py"}
    assert shared == expected
    # Integration evidence: radia-gpu LTspice tests cross-check the radia-mcp LTspice gates.
    gpu_ltspice_tests = list((ROOT / "packages/radia-gpu/tests/ltspice").glob("test_*.py"))
    assert any("from radia_mcp.ltspice" in path.read_text(encoding="utf-8") for path in gpu_ltspice_tests)
    # The radia_mcp modules the LTspice server imports are part of that integration surface.
    server = (ROOT / "packages/radia-mcp/src/radia_mcp/ltspice/server.py").read_text(encoding="utf-8")
    imported = set(re.findall(r"^\s*from radia_mcp\.common(?:\.(\w+))? import", server, re.M))
    modules = {"packages/radia-mcp/src/radia_mcp/common/__init__.py"} | {
        f"packages/radia-mcp/src/radia_mcp/common/{name}.py" for name in imported if name}
    if "" in imported:  # package-level names resolve lazily through common/__init__.py
        modules.add("packages/radia-mcp/src/radia_mcp/common/status.py")
    lane = MANIFEST["lanes"]["radia-gpu"]["events"]
    for module in sorted(modules):
        assert module in shared, module
        assert all(CONTRACT.selected(module, lane[event]) for event in ("push", "pull_request")), module
    objective = ROOT / "tests/matlab/fixtures/radia_optuna_ltspice_objective.m"
    assert "ltspice_rc.cir" in objective.read_text(encoding="utf-8")
    assert "ltspice_rc.cir" in (ROOT / "tests/matlab/test_ltspice_workflow.m").read_text(encoding="utf-8")


def test_root_prose_and_release_events_do_not_become_native_path_gates():
    for lane in MANIFEST["lanes"].values():
        for event, value in lane["events"].items():
            if event not in ("push", "pull_request"):
                assert "paths" not in value and "paths-ignore" not in value
            elif lane["native"] and ("branches" in value or event == "pull_request"):
                assert not CONTRACT.selected("README.md", value)
                assert not CONTRACT.selected("AGENTS.md", value)
    native = MANIFEST["lanes"]["build-test"]["events"]
    assert "pull_request" not in native and native["push"]["tags"] == ["v*"]
    assert MANIFEST["root_content_diff"]["unknown"] == "all_owners"
    assert MANIFEST["root_content_diff"]["prose"] == "no_native_build"


def test_solver_triggers_exclude_only_independent_workflows_and_policy_tests():
    workflow = yaml.load((ROOT / ".github/workflows/radia-fast.yml").read_text(), Loader=yaml.BaseLoader)
    for event in ("push", "pull_request"):
        ignored = workflow["on"][event]["paths-ignore"]
        excluded = lambda path: any(fnmatchcase(path, pattern) for pattern in ignored)
        assert all(excluded(".github/workflows/" + name) for name in INDEPENDENT)
        assert all(not excluded(".github/workflows/" + name) for name in SOLVER)
        assert excluded("tests/policy/test_release_workflow_ref_gate.py")
        assert excluded("packages/eqnedit64/src/eqnedit64/api.py")
        assert excluded("tools/eqnedit64/src/office_html.cpp")
        assert excluded("packages/radia-mcp/src/radia_mcp/presentation/tools.py")
        for path in ("src/radia/__init__.py", "Build.ps1", "tests/test_wheel_package_policy.py", "tools/ci_preflight.py"):
            assert not excluded(path), path


def test_policy_contracts_run_in_policy_lane():
    workflow = (ROOT / ".github/workflows/policy-lint.yml").read_text()
    assert "python -m pytest tests/policy/ -q" in workflow
    manifest = (ROOT / "tests/test_tier_manifest.json").read_text()
    assert "test_release_workflow_ref_gate.py" not in manifest


def test_policy_lane_runs_for_every_distribution():
    workflow = yaml.load((ROOT / ".github/workflows/policy-lint.yml").read_text(), Loader=yaml.BaseLoader)
    for event in ("push", "pull_request"):
        assert "paths-ignore" not in workflow["on"][event]
        assert "paths" not in workflow["on"][event]


# Approved PR1 exception ceiling. Never regenerate this to make a failure pass.
# Deleting manifest entries is allowed; new scopes require an explicit policy review.
EXCEPTION_CEILING = {
    'CI-MIG-001': 'a7aa2794057c02349a7c6632956d8bf139c06832c85578bef01d06699ec083e0',
    'CI-MIG-002': 'e8d32b61efc3bb76b82fb88dd4d506863bcf2817833a95c397cd511c538ba791',
    'CI-MIG-003': 'b7f2e03d0aa4d997da616ce7137b309e3086982f08dc2ac30070fe2db4442cbd',
    'CI-MIG-004': '9076fe16621f02ff4e2c4d41e70ddf54103fd51272f3b39d0eab5f5e6979028b',
    'CI-MIG-005': '04dea927dfef1f3ab3c64ad7a1d66779fa934c99e388801d9893b91cf300a7e8',
    'CI-MIG-006': '50a7f0b75a0b16c2ee9a276c28d160d3e82660b6276a7a7bd140d6746e026bed',
    'CI-MIG-007': 'ddf3efa61eca7ef9a096d4b6a1cabcd537baa9002a4986c6ca2a81ed5166f6d7',
    'CI-MIG-008': '4d1662d7d4628b43292029de276630009c9c4d6c174c44fd4f49cf6f59854980',
    'CI-MIG-009': 'e4e457d648638f7bb056bb54b4d282841d7b49d066f9d667b442b658536e3dbb',
    'CI-MIG-010': '6be02f534bd117a55ede67d4041fd3639d891226905228e25964ec7237b1cfc9',
    'CI-MIG-011': 'aa049d39dc6cf44ebc3bce510b0ef93a122b6dc9663c6b68a5d1b6f2bd52dcd1',
    'CI-MIG-012': '8df2de88496418c6a21091c85f25f4581a459a1deac2e4a054c2e96577ebb878',
    'CI-MIG-013': '3b60bb5aabcfe6e4e09d2e526eb51028b30ead36be8a4ab5218d3e62a903a55e',
    'CI-MIG-014': 'a784ec4a716a9558c1d1ead52613e9062d81362bf47ab9f68b045258bb31cd6c',
    'CI-MIG-015': '1a49afc997ab5ce46a3366c2c5b2424789f1f564ee43e23706f8fe6456635fe4',
    'CI-MIG-016': 'f1c00ed294100bb51adcf785ff1bc9935003b1d858e0287b7e42f2d33ab8486f',
    'CI-MIG-017': 'bb2756904c2f52fc7cc7b40e0cf3c694a902149870369f8d9ffe0e89d253c137',
    'CI-MIG-018': 'eab6dc5f578128a541d0cc70c2bc88ffa61b3410c44b2b03697b2dc57f9efcd7',
    'CI-MIG-019': '3073f83cdbd6ba39abe36e44b08a1950a3e2aa60256c64f84b4e065611ef45a0',
    'CI-MIG-020': '806f98679f8ada4df64ac903593d881c4010ea6ae43e8a51975d8e2f6a54fc76',
    'CI-MIG-021': 'a4e804f774594e5ea98b6a219b629dc9d3ce1040aea630df19b801817694c646',
    'CI-MIG-022': 'ffd4f478062a027aff741766dcd07e90edd2652a3717801ff5d0f18349778bcd',
    'CI-MIG-023': 'aa5650cecfbad950493460d176bf18eef77f0a9ee0485cf794881d4ffee1626f',
    'CI-MIG-024': '71bd9424fc9c02bec174eb028b158dc14b9bcb5fece5b9c1bae6ea99ef442c53',
    'CI-MIG-025': '3e97cbc4e5223480a51b0583c2ac777657594ac482723d20c61b0009e05ddf37',
    'CI-MIG-026': 'a4cc77d3b5da4d0218aab1272129cb8156eb0f218dc917b7dc69e5d7863b9075',
    'CI-MIG-027': 'ff4a9bbac8a94fa8d85672013f956f5816367b3136a309543d41c1102fb8216f',
    'CI-MIG-028': 'a6ea1cb87ae0ffd766cbd0eda947ce03c547ef896d096649b3b6567a78be3315',
    'CI-MIG-029': '05052b41d45e13dbac4d9b2e616b7596e497dbd79ca215a9910f3c4a82771678',
    'CI-MIG-030': '581601619fd8a694edfc34ed83583fd050f5b340a32640ea8a130b43ff5d82b2',
    'CI-MIG-031': '9b4f9bb7703e133d3346af7c2351164b55a185aef54fbf55574362d6ea3ae9eb',
    'CI-MIG-032': '87fd8d503e562dff3b4fbe20b832efaed686ed22d299ea5b71e0e4514b12efd6',
    'CI-MIG-033': '4dcda1dfc920e3fac051bff219e6bf69182498b23e061b695889ddc4a1e21e56',
    'CI-MIG-034': 'dc05bbe2426de3f62500e80254e0d97cf617c2687dde7e9150f2bdb207ef4c56',
    'CI-MIG-035': '8d2d9ff97dbb84de8eb55f28298ef271064475f35a4257464c65022f3033fcee',
    'CI-MIG-036': 'dcac224cea2a7980adf6eecb92591e8e472933e25d1402623d09744dc4821fb8',
    'CI-MIG-037': 'bfb15c7cfe4d8f357f2241a81f882ecfe0eec3b40243229ca24a4d5e0aeedafc',
    'CI-MIG-038': '7ae4b091120032a02f925c93c9ba906d77be57e57596bab2ff47fb43865f5004',
    'CI-MIG-039': 'b9b1128989e6c73fe9367dd092bbc96909e4226c1421815d91b2f4653bea2054',
    'CI-MIG-040': '6a1fed0a83df72d3f2d420c8eeefeaaba992578f8bb4cd37e77bb6ed36dd4b05',
    'CI-MIG-041': 'a3f0a55494f31bd8fd0ed610603f6b379336b087808545b1ff150815bdc3883d',
    'CI-MIG-042': 'a2454b440bd700e8d0470f6f223fb61ffde7b568f7f15cc9b2c004fa773847cf',
    'CI-MIG-043': '402c788ce573d79e076ad8c615f2a5c26c20cf89c0ffc76b1aee726704c9e5b1',
    'CI-MIG-044': 'd5c19ba93ab4c6d1c62a0bcd01b194d7de9b6ec8a48f073d80a5164cb08be9d0',
    'CI-MIG-045': 'c4802f4f7b4f4dd36f20090974e87d6c5649278a321bcda94ae94836198a6d12',
    'CI-MIG-046': 'c8076b1f98283c035b534bc4d3e0562d00c4997193ef921872bed7717e84db0e',
    'CI-MIG-047': '311762c3c3cdad9ea6d018103c0e8a5924683e8a7248d9a52efab8bdb89ed8fc',
    'CI-MIG-048': '2b4364bb05ff1d65cafde8146e9356ea9ca25abdd020472767dd08380b90a24a',
    'CI-MIG-049': 'c78c10f223cdb1cd2b54be8740607b9676a65d79f131646c5f346949df4401fb',
    'CI-MIG-050': '5c99e56ff54fc1fc288149e2eca94aa5785e530f7fddd2cb914f6ec06eabeeb9',
    'CI-MIG-051': '50a4f3edfe06b3666e8ca6337a03f7c2bc5f1f0bf0e3e7f5afe213f84b83be80',
    'CI-MIG-052': 'c6172e3626a6bb3e25270f11b572c82e68e4d96937d69e36d9b400cb2eb59bd1',
    'CI-MIG-053': '1ea282436c4fa4d882df925dd0ae7e8675863ffe558a835235eaba0a4c2cd838',
    'CI-MIG-054': 'b0ccd858a76015732b1274ff77e8912a0d60c7934dcb7916be96bd50f6bf4976',
    'CI-MIG-055': 'baa682c1bb43b94d9e3ee46c71b34fcfa5c8307784641eb64e7044535ad6faea',
    'CI-MIG-056': 'b2677020f60c8110e0376bd8dafbc5d4d9f0adddd1541cfac2994ac1950124e6',
    'CI-MIG-057': '9e5ee4ba3b06fff7f8b0359f6caf11978b30bc1ce39ca59a510897ad928bf69f',
    'CI-MIG-058': '2eeac18b32bb81edf5a6b2ea70ef580dfe590b43702b2251bbee2399439559f3',
    'CI-MIG-059': '43ef0696e529011eb4b722d87b800ee6b4c7543f22bd8bcc4503f2412d3f5e52',
    'CI-MIG-060': '82ee1b24098a27337b66125e781726ed518bb43876937e613c0fbcf1816c8b46',
    'CI-MIG-061': '693dbcd3064b632a06b6e47d476f8958e30661000527c117fb70b1116909e3e5',
    'CI-MIG-062': 'fd8ee10e73875d1b82dd525f1188152e83fed6e4d17bfbc9bba259b07fc75e3e',
    'CI-MIG-063': 'a3adfbf8a65247bacc1460e18e715b80036cfc6dfddf3eafc020f4f3ba0d7558',
    'CI-MIG-064': 'f2730e569642439f646909d42b6eadf13849597964e260cc0c8992c0932c544d',
    'CI-MIG-065': 'ee78278b20de27d0e938eb15fc6f220f1b2328ea3d9c281cdc607b4973286841',
    'CI-MIG-066': 'f8f31e192253a56a2b4846dac61b7919696775b0ef921bfda62da33d97f05814',
    'CI-MIG-067': 'dd74e1b1f8f6981f7914aa60db74c7ed7a3e93ffd8a5c956c71c9c4f6176085e',
    'CI-MIG-068': 'e89a73b7035593591e2cc2cfdd4fe56c50e5ce1f5d20ab4e9546d8c49deaf5dc',
    'CI-MIG-069': 'c564ba609eda79f36157a153e8110f9d96f4973a5fe5d90ef9816808665bb7f1',
    'CI-MIG-070': '34d4aae0659ec821d910d0ac579ccbe94ba616a5272ba75268fb5b67c1957387',
    'CI-MIG-071': 'b540a05c314c2604f6436d07e0f6d90ed6817a9e0a700b983efb1e2d50167855',
    'CI-MIG-072': '61fcb68d6cb7492d11bef6fbb6f4b9604d644dbe7d751108822e4e2e0d6b9622',
    'CI-MIG-073': '18992b1813e0f2840407021920799a68cb1d1f2d1ffa37a536c231ef3c7c7fc7',
    'CI-MIG-074': '8ce7b4952a4f53c1bb81b652f476219f989d0d2480cb771ed953193c99b1a348',
    'CI-MIG-075': 'f5fe8f53c10cd4ee8661a7a74fe859361b30cd664148878012848ab36dd5c25f',
    'CI-MIG-076': 'a85ab6c7424e1fe439820649d22a672683325139e2d715c83a8f9c49dcc8a71f',
    'CI-MIG-077': '22384f8e9d9d6ce55da27a1e3c6e0e59a54033bdd939b57691e50423c09b40f7',
    'CI-MIG-078': 'a6ad24f7fefc5644a3ecf37339047ad03c5c8cae67a897afe9dc2eef791c8cd7',
    'CI-MIG-079': '801bee4f0df68628bc2e14f4d49661b5b5a1d8548caf69034ae985b87206406e',
    'CI-MIG-080': '82d4ac62f1bb56e71b03dd17d016f77fc24e051ddb00079429349bf61681490a',
    'CI-MIG-081': 'a769c8eabc8b8c5d21ba75376d45918029b776c9325a09c4723301625a38141c',
    'CI-MIG-082': '0056a562d19d5a3ff48170c20e9883836e6a0ed2d50da42dbeb214e465644ee6',
    'CI-MIG-083': 'd80f2d57ed6a6d40fd25b3b44e8a1de02330c4d0a64ec6b8252214ea9c5c272d',
    'CI-MIG-084': '74ab1b3ccf1adea48ff5aa53fb2eccaac98c0b098b03d72ff12640ce997e6b80',
}


def test_job_conditions_are_single_expressions():
    """A `${{ }}` condition followed by stray text becomes a truthy string.

    That silently ran tag-only release jobs on pull requests, so every job
    condition that opens an expression must also end with it.
    """
    for workflow in sorted((ROOT / ".github/workflows").glob("*.yml")):
        document = yaml.load(workflow.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
        for name, job in (document.get("jobs") or {}).items():
            condition = str(job.get("if", "")).strip()
            if condition.startswith("${{"):
                assert condition.endswith("}}") and condition.count("${{") == 1, (workflow.name, name, condition)


def test_cubit_lane_reincludes_every_packaged_markdown_pattern():
    import tomllib
    package = ROOT / "packages/cubit-mesh-export"
    project = tomllib.loads((package / "pyproject.toml").read_text(encoding="utf-8"))
    expected = {f"packages/cubit-mesh-export/{project['project']['readme']}"}
    for module, patterns in project["tool"]["setuptools"]["package-data"].items():
        folder = "packages/cubit-mesh-export/src/" + module.replace(".", "/")
        expected |= {f"{folder}/{p}" for p in patterns if p.endswith(".md")}
    for event in ("push", "pull_request"):
        paths = MANIFEST["lanes"]["cubit-mesh-export"]["events"][event]["paths"]
        assert "!**/*.md" in paths
        assert expected <= set(paths[paths.index("!**/*.md") + 1:]), event
