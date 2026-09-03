"""Keep the mdx-only CI split and remote pre-push contract explicit."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_fast_ci_runs_only_on_mdx_and_native_is_a_named_release_lane():
    fast = (ROOT / ".github" / "workflows" / "radia-fast.yml").read_text(encoding="utf-8")
    native = (ROOT / ".github" / "workflows" / "build-test.yml").read_text(encoding="utf-8")
    optuna = (ROOT / ".github" / "workflows" / "radia-optuna.yml").read_text(encoding="utf-8")
    release = (ROOT / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")

    assert "name: Radia\n" in fast
    assert "runs-on: [self-hosted, Windows, X64, mdx]" in fast
    assert "Build with MSVC" not in fast
    assert "validation_test/" not in fast
    assert "name: Radia Native Release" in native
    assert "runs-on: [self-hosted, Windows, X64, mdx]" in native
    assert "mkl-devel" in native
    assert "pybind11==3.0.2" in native
    assert "ninja" in native
    assert '"trimesh>=4.0"' in native
    assert "twine==6.2.0" in native
    assert "pytest-rerunfailures" not in native
    assert "--reruns" not in native
    assert "tools/run_test_tier.py --profile native-smoke" in native
    assert "'tests/'" not in native
    assert '"MKLROOT=$mklRoot"' in native
    assert "MKLROOT=C:\\Program Files\\Python312\\Library" not in native
    assert "runs-on: [self-hosted, Windows, X64, mdx]" in optuna
    assert "windows-radia" not in optuna
    optuna_self_hosted = optuna.split("  installed-wheel-matlab-e2e:", 1)[0]
    assert "actions/setup-python" not in optuna_self_hosted
    assert "RADIA_OPTUNA_CI_PYTHON" in optuna_self_hosted
    assert "-m venv" in optuna_self_hosted
    assert "pip install pybind11==3.0.2 ninja" in optuna_self_hosted
    assert optuna.count("'.github/workflows/radia-optuna.yml'") == 2
    assert optuna.count("'tests/test_release_quad_optuna_candidate.py'") == 2
    assert 'workflows: ["Radia Native Release"]' in release


def test_regular_ci_never_selects_a_lab_runner():
    workflows = ROOT / ".github" / "workflows"
    signing_release = "release-eqnedit64-pypi.yml"

    for path in sorted(workflows.glob("*.yml")):
        if path.name == signing_release:
            continue
        runner_lines = [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip().startswith("runs-on:")
        ]
        for line in runner_lines:
            normalized = line.lower()
            assert "windows-radia" not in normalized, f"{path.name}: {line}"
            assert "lab" not in normalized, f"{path.name}: {line}"


def test_distribution_ci_is_change_scoped_and_mcp_full_suite_is_explicit():
    fast = (ROOT / ".github" / "workflows" / "radia-fast.yml").read_text(
        encoding="utf-8"
    )
    policy = (ROOT / ".github" / "workflows" / "policy-lint.yml").read_text(
        encoding="utf-8"
    )
    mcp = (ROOT / ".github" / "workflows" / "radia-mcp-matrix.yml").read_text(
        encoding="utf-8"
    )

    for independent_path in (
        "packages/radia-mcp/**",
        "packages/cubit-mesh-export/**",
        "packages/radia-optuna/**",
        "packages/eqnedit64/**",
    ):
        assert independent_path in fast

    assert 'paths-ignore:' in policy
    assert '"packages/radia-mcp/**"' in policy
    assert "Classify radia-mcp change scope" in mcp
    assert 'git diff --name-only "$BASE_SHA...$GITHUB_SHA"' in mcp
    assert "mode=targeted" in mcp
    assert "python_versions='[\"3.10\",\"3.11\",\"3.12\"]'" in mcp
    assert "python_versions='[\"3.12\"]'" in mcp
    assert "python-version: ${{ fromJSON(needs.scope.outputs.python_versions) }}" in mcp
    assert "Prepare impact-scoped test plan" in mcp
    assert "packages/radia-mcp/tools/select_ci_tests.py" in mcp
    assert '--changed-files-json "$CHANGED_FILES_JSON"' in mcp
    assert "Fetch impact comparison base" in mcp
    assert 'base_args=(--base "$CI_BASE_SHA")' in mcp
    assert "--full --output radia-mcp-ci-selection.json" in mcp
    assert 'plan["package_tests"]' in mcp
    assert 'env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"' in mcp
    assert 'else plan["package_tests"]' in mcp
    assert 'plan["server_selftests"]' in mcp
    assert 'plan["run_mcp_response_tests"]' in mcp
    assert "Meta health (all cataloged subpackages must import)" not in mcp
    assert "radia_mcp_health" not in mcp
    assert "pytest-rerunfailures" not in mcp
    assert "--reruns" not in mcp
    assert "mode=release" in mcp
    assert "needs.scope.outputs.build_wheel == 'true'" in mcp


def test_pre_push_runs_the_unpushed_candidate_on_mdx():
    hook = (ROOT / "tools" / "git-hooks" / "pre-push").read_text(encoding="utf-8")
    helper = (ROOT / "tools" / "ci_preflight_mdx.py").read_text(encoding="utf-8")

    assert "ci_preflight_mdx.py --base \"$remote_sha\" --head \"$local_sha\"" in hook
    assert '"bundle", "create"' in helper
    assert 'git_output("merge-base", head, main_ref)' in helper
    assert "$base = '{effective_base}'" in helper
    assert "scp" in helper
    assert "tools/ci_preflight.py --since $base" in helper
    assert "'mcp>=1.0,<2'" in helper
    assert "tools/run_test_tier.py --profile fast-contracts" in helper
    assert "upload_release_asset.py" not in hook
    assert "developer push must never upload mutable" in hook


def test_policy_twins_define_the_same_mdx_notebook_contract():
    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    claude = (ROOT / "CLAUDE.md").read_text(encoding="utf-8")
    heading = "### CI Execution, Validation Evidence, and Notebook Policy (2026-09-02)"

    def section(text: str) -> str:
        start = text.index(heading)
        return text[start:text.index("### ", start + len(heading))]

    policy = section(agents)
    assert policy == section(claude)
    normalized = " ".join(policy.split())
    assert "mdx gives CI and preflight work priority" in normalized
    assert "use hibino first when it is available" in normalized
    assert "the mdx CI queue is idle" in normalized
    assert "CI scope begins at the independently released distribution boundary" in normalized
    assert "generated-inventory checks run when the inventory or its generator changes" in normalized
    assert "normal pull-request and main-push CI runs a stable compact contract set" in normalized
    assert "tests selected from the changed source/test paths" in normalized
    assert "only the affected server selftests" in normalized
    assert "complete package pytest suite, all-server selftests" in normalized
    assert "explicit full-audit workflow" in normalized
    assert "normal CI optimizes for fast, high-signal feedback" in normalized
    assert "does not automatically rerun a failed deterministic test" in normalized
    assert "run only the relevant validation lane and retain its result JSON" in normalized
    assert "Do not keep two tests whose purpose and failure signal are the same" in normalized
    assert "an adjacent JSON and a runtime gate are not required" in normalized
    assert "A docs-only contract lane parses changed notebooks" in normalized
    assert "Developer pre-push hooks run only the impact-scoped mdx preflight" in normalized
