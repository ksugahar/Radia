"""Guards: release jobs cannot publish unapproved or substituted artifacts.

Assumes: GitHub enforces the configured environment approval; remote tag and
PyPI availability are independent services, replaced only in offline contracts.
"""
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[2]
RELEASE_WORKFLOWS=(
    ROOT/".github"/"workflows"/"release.yml",
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
        assert "publication-hold:" not in workflow
        assert "workflow_run:" not in workflow
        assert "needs: verify-promotion" in workflow
        publish=workflow[workflow.index("  publish-pypi:"):]
        assert "github.event_name == 'workflow_dispatch'" in publish
        assert "github.ref_type == 'tag'" in publish
        assert "startsWith(github.ref, 'refs/tags/v')" in publish
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


def test_native_assets_use_accepted_tag_artifact_after_promotion():
    workflow = RELEASE_WORKFLOWS[0].read_text(encoding="utf-8")
    upload = workflow.split("  upload-binaries:", 1)[1].split("  verify-promotion:", 1)[0]
    assert "needs: [verify-promotion, publish-pypi]" in upload
    assert "needs.publish-pypi.result == 'success'" in upload
    assert "environment: pypi" not in upload
    assert "needs.verify-promotion.result == 'success'" in upload
    assert "artifact-ids: ${{ needs.verify-promotion.outputs.native_artifact_id }}" in upload
    assert "--after-approval" in upload
    assert "head_branch" not in upload


def test_radia_publication_graph_has_one_approval_and_no_failure_bypass():
    import yaml
    workflow = yaml.safe_load(RELEASE_WORKFLOWS[0].read_text(encoding="utf-8"))
    assert set(workflow.get("on", workflow.get(True))) == {"workflow_dispatch"}
    jobs = workflow["jobs"]
    assert set(jobs) == {"verify-promotion", "publish-pypi", "upload-binaries"}
    assert jobs["publish-pypi"]["environment"] == "pypi"
    assert "environment" not in jobs["verify-promotion"] and "environment" not in jobs["upload-binaries"]
    assert jobs["upload-binaries"]["needs"] == ["verify-promotion", "publish-pypi"]
    for name, job in jobs.items():
        condition = job["if"]
        assert "github.ref_type == 'tag'" in condition
        assert "always(" not in condition and "||" not in condition
        assert not job.get("continue-on-error")
        for step in job["steps"]:
            assert not step.get("continue-on-error")
            if step.get("uses", "").startswith("actions/checkout"):
                assert step["with"] == {"ref": "${{ github.sha }}", "persist-credentials": False}
        if name != "verify-promotion":
            audit = next(i for i,s in enumerate(job["steps"]) if "--after-approval" in s.get("run", ""))
            write = next(i for i,s in enumerate(job["steps"]) if s.get("uses", "").startswith("pypa/") or "upload_release_asset.py" in s.get("run", ""))
            assert audit < write


def test_radia_publisher_fails_existing_and_uses_env_for_shell_inputs():
    import re
    import yaml
    jobs = yaml.safe_load(RELEASE_WORKFLOWS[0].read_text(encoding="utf-8"))["jobs"]
    publish = jobs["publish-pypi"]
    assert publish["permissions"] == {"contents": "read", "id-token": "write"}
    upload = next(s for s in publish["steps"] if s.get("uses", "").startswith("pypa/gh-action-pypi-publish"))
    assert upload.get("with", {}).get("skip-existing", False) is False
    for job in jobs.values():
        for step in job["steps"]:
            assert not re.search(r"\$\{\{\s*inputs(?:\.|\[)", step.get("run", "")), step.get("name", step)
    native_upload = next(s for s in jobs["upload-binaries"]["steps"] if "upload_release_asset.py" in s.get("run", ""))["run"]
    # This shell step must itself fail on no files, even if the upstream
    # validator is accidentally weakened; never turn empty publication green.
    assert 'test -n "$FILES"' in native_upload
    assert not re.search(r"exit\s+0", native_upload)
