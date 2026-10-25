from pathlib import Path


ROOT=Path(__file__).resolve().parents[2]
RELEASE_WORKFLOWS=(
    ROOT/".github"/"workflows"/"release.yml",
)
OPTUNA_RELEASE_WORKFLOW=(
    ROOT/".github"/"workflows"/"release-radia-optuna.yml"
)
PACKAGE_CI_CONTEXTS=(
    (
        ROOT/".github"/"workflows"/"build-test.yml",
        "radia.ci-release-context.v1",
        "ci-release-context",
    ),
)


def test_pypi_distributions_have_independent_ci_boundaries():
    workflows={
        "Radia": ROOT/".github"/"workflows"/"radia.yml",
        "Radia Native Release": ROOT/".github"/"workflows"/"build-test.yml",
        "radia-mcp": ROOT/".github"/"workflows"/"radia-mcp-matrix.yml",
        "cubit-mesh-export": ROOT/".github"/"workflows"/"cubit-mesh-export.yml",
        "radia-optuna": ROOT/".github"/"workflows"/"radia-optuna.yml",
    }
    texts={name: path.read_text(encoding="utf-8")
           for name, path in workflows.items()}
    for name, text in texts.items():
        assert f"name: {name}" in text

    fast=texts["Radia"]
    native=texts["Radia Native Release"]
    assert "workflow_dispatch:" in fast
    import yaml
    labels = yaml.safe_load(fast)['jobs']['fast-contracts']['runs-on']
    assert labels == ['self-hosted', 'Windows', 'X64', 'ci-light']
    assert '100' not in labels
    assert "Build with MSVC" not in fast
    assert "compileall -q src/radia" in fast
    assert "compileall -q src/radia packages/radia-mcp/src" not in fast
    assert "workflow_dispatch:" in native
    assert "tags: ['v*']" in native
    trigger = native[:native.index("workflow_dispatch:")]
    assert "paths-ignore:" not in trigger
    assert "radia-mcp-wheel" not in native
    assert "cubit-mesh-export-wheel" not in native
    assert "radia-optuna-wheel" not in native
    assert "tools/run_test_tier.py --profile native-smoke" in native
    assert "-m pip install build wheel pytest mcp" in native
    assert "'tests/'" not in native
    assert "pytest-rerunfailures" not in native
    assert "--reruns" not in native

    policy=(ROOT/".github"/"workflows"/"policy-lint.yml").read_text(
        encoding="utf-8")
    assert "packages/radia-mcp/tools/policy_lint.py" not in policy
    assert "packages/radia-mcp/tools/policy_lint.py" in texts["radia-mcp"]

    assert "tags: ['radia-mcp-v*']" in texts["radia-mcp"]
    assert "tags: ['cubit-mesh-export-v*']" in texts["cubit-mesh-export"]
    assert "tags: ['radia-optuna-v*']" in texts["radia-optuna"]


def test_ci_records_exact_ref_context_before_release():
    for path, schema, artifact in PACKAGE_CI_CONTEXTS:
        workflow=path.read_text(encoding="utf-8")
        assert schema in workflow
        assert "GITHUB_REF_TYPE" in workflow
        assert "GITHUB_SHA" in workflow
        assert "GITHUB_RUN_ID" in workflow
        assert f"name: {artifact}" in workflow
        assert "GIT_TERMINAL_PROMPT = '0'" in workflow
        assert "git -c credential.helper= ls-remote" in workflow


def test_pypi_release_jobs_require_the_triggering_ci_to_be_a_tag_run():
    for path in RELEASE_WORKFLOWS:
        workflow=path.read_text(encoding="utf-8")
        assert "publication-hold:" in workflow
        assert "needs: verify-promotion" in workflow
        publish=workflow[workflow.index("  publish-pypi:"):]
        assert "github.event_name == 'workflow_dispatch'" in publish
        assert "github.ref == 'refs/heads/main'" in publish
        assert "needs.verify-promotion.result == 'success'" in publish
        assert "artifact-ids: ${{ needs.verify-promotion.outputs.artifact_id }}" in publish
        assert "run-id: ${{ inputs.ci_run_id }}" in workflow
        gate=(ROOT/"tools/verify_radia_promotion.py").read_text(encoding="utf-8")
        assert 'context.get("ref_type") == "tag"' in gate
        assert 'run.get("event") == "push"' in gate
        assert 'git/ref/tags/' in gate


def test_radia_mcp_publishes_only_from_its_successful_tag_ci():
    workflow = (ROOT / ".github" / "workflows" / "radia-mcp-matrix.yml").read_text(
        encoding="utf-8"
    )
    assert "workflow_run:" not in workflow
    assert "release-radia-mcp.yml" not in workflow
    assert "radia-mcp-ci-release-context" not in workflow
    assert "needs: [ci-impact, scope, selftest, minimum-sdk, release-evidence]" in workflow
    assert "Require a successful full audit on the tagged commit" in workflow
    assert "head_sha=${GITHUB_SHA}&event=workflow_dispatch&status=success" in workflow
    assert "full audit evidence found for ${GITHUB_SHA}" in workflow
    assert "github.ref_type == 'tag'" in workflow
    assert "startsWith(github.ref_name, 'radia-mcp-v')" in workflow
    assert "Download wheel built by this tag CI" in workflow
    assert "${GITHUB_REF_NAME#radia-mcp-v}" in workflow
    assert "pypa/gh-action-pypi-publish@release/v1" in workflow


def test_cubit_mesh_export_publishes_only_from_its_successful_tag_ci():
    workflow = (ROOT / ".github" / "workflows" / "cubit-mesh-export.yml").read_text(
        encoding="utf-8"
    )
    assert "workflow_run:" not in workflow
    assert "release-cubit-mesh-export.yml" not in workflow
    assert "cubit-mesh-export-ci-release-context" not in workflow
    assert "needs: [ci-impact, build-test]" in workflow
    assert "github.ref_type == 'tag'" in workflow
    assert "startsWith(github.ref_name, 'cubit-mesh-export-v')" in workflow
    assert "Download wheel built by this tag CI" in workflow
    assert "${GITHUB_REF_NAME#cubit-mesh-export-v}" in workflow
    assert "cubit_mesh_export.ccm" in workflow
    assert "cubit_mesh_curver.pyd" in workflow
    assert "pypa/gh-action-pypi-publish@release/v1" in workflow


def test_optuna_manual_release_selects_one_fully_successful_ci_run():
    workflow=OPTUNA_RELEASE_WORKFLOW.read_text(encoding="utf-8")
    assert "workflow_dispatch:" in workflow
    assert "workflow_run:" not in workflow
    assert "ci_run_id:" in workflow
    assert "candidate_sha256:" in workflow
    assert 'STATUS" != "completed"' in workflow
    assert 'CONCLUSION" != "success"' in workflow
    assert 'EVENT" != "push"' in workflow
    assert 'HEAD_BRANCH" != "main"' in workflow
    assert (
        'WORKFLOW_PATH" != ".github/workflows/radia-optuna.yml"'
        in workflow
    )
    assert 'HEAD_REPOSITORY" != "$GITHUB_REPOSITORY"' in workflow
    assert '"build-test"' in workflow
    assert '"radia-optuna installed-wheel MATLAB/Simulink E2E"' in workflow
    assert "run-id: ${{ steps.source.outputs.run_id }}" in workflow
    assert "ref: ${{ steps.source.outputs.sha }}" in workflow
    assert (
        'python3 packages/radia-optuna/verify_wheel.py "$WHL" '
        "--release-candidate --json"
        in workflow
    )
    assert '"packages/radia-optuna/matlab/optuna_mex.mexw64"' in workflow
    assert (
        '"radia_optuna/matlab/optuna_mex.mexw64": Path('
        in workflow
    )
    assert '"radia_optuna/matlab/LICENSE": Path("LICENSE")' not in workflow
    assert "target.write_bytes(archive.read(member))" in workflow
    assert 'EXPECTED_TAG="refs/tags/radia-optuna-v${VERSION}"' in workflow
    assert 'ACTUAL_SHA256=$(sha256sum "$WHL"' in workflow
    assert '"${ACTUAL_SHA256,,}" != "${CANDIDATE_SHA256,,}"' in workflow
    assert "contents: write" in workflow
    assert 'gh release create "$TAG_NAME" "$WHEEL"' in workflow


def test_optuna_manual_release_fails_loudly_when_ci_sha_has_no_tag():
    workflow=OPTUNA_RELEASE_WORKFLOW.read_text(encoding="utf-8")
    tag_check=workflow[workflow.index("- name: Confirm tag points"):]
    tag_check=tag_check[:tag_check.index("- name: Check out the exact")]
    assert 'echo "::error::No radia-optuna-v* tag points at $SHA"' in tag_check
    assert "exit 1" in tag_check








def test_native_assets_use_accepted_tag_artifact_after_promotion():
    workflow = RELEASE_WORKFLOWS[0].read_text(encoding="utf-8")
    upload = workflow.split("  upload-binaries:", 1)[1].split("  publication-hold:", 1)[0]
    assert "needs: verify-promotion" in upload
    assert "needs.verify-promotion.result == 'success'" in upload
    assert "run-id: ${{ inputs.ci_run_id }}" in upload
    assert "head_branch" not in upload
