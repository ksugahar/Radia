from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
RELEASE_WORKFLOWS=(
    ROOT/".github"/"workflows"/"release.yml",
    ROOT/".github"/"workflows"/"release-radia-mcp.yml",
    ROOT/".github"/"workflows"/"release-cubit-mesh-export.yml",
)
OPTUNA_RELEASE_WORKFLOW=(
    ROOT/".github"/"workflows"/"release-radia-optuna.yml"
)


def test_ci_records_exact_ref_context_before_release():
    workflow=(ROOT/".github"/"workflows"/"build-test.yml").read_text(
        encoding="utf-8")
    assert "schema = 'radia.ci-release-context.v1'" in workflow
    assert "if ($env:GITHUB_REF_TYPE -eq 'tag')" in workflow
    assert 'git ls-remote --tags $repoUrl @patterns' in workflow
    context=workflow[workflow.index("- name: Record exact CI ref context"):]
    context=context[:context.index("- name: Upload exact CI ref context")]
    assert "ls-remote --tags origin" not in context
    assert "ref_type = $env:GITHUB_REF_TYPE" in workflow
    assert "sha = $env:GITHUB_SHA" in workflow
    assert "run_id = [string]$env:GITHUB_RUN_ID" in workflow
    assert "release_tags = @($releaseTags | Sort-Object -Unique)" in workflow
    assert "name: ci-release-context" in workflow


def test_pypi_release_jobs_require_the_triggering_ci_to_be_a_tag_run():
    for path in RELEASE_WORKFLOWS:
        workflow=path.read_text(encoding="utf-8")
        context=workflow.index("Download exact CI ref context")
        verification=workflow.index("Verify triggering CI ran on a tag ref")
        tag_lookup=workflow.index("id: tag_check")
        assert context < verification < tag_lookup
        assert "github.event.workflow_run.event == 'push'" in workflow
        assert 'context.get("ref_type") == "tag"' in workflow
        assert 'context.get("release_tags")' in workflow
        assert "if: steps.ref_check.outputs.eligible == 'true'" in workflow
        assert "run-id: ${{ github.event.workflow_run.id }}" in workflow


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
        'WORKFLOW_PATH" != ".github/workflows/build-test.yml"'
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
    assert 'target = Path("matlab/optuna_mex.mexw64")' in workflow
    assert (
        'archive.read("radia_optuna/matlab/optuna_mex.mexw64")'
        in workflow
    )
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
