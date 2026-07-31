from pathlib import Path


ROOT=Path(__file__).resolve().parents[1]
RELEASE_WORKFLOWS=(
    ROOT/".github"/"workflows"/"release.yml",
    ROOT/".github"/"workflows"/"release-radia-mcp.yml",
    ROOT/".github"/"workflows"/"release-cubit-mesh-export.yml",
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
