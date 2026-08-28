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
EQNEDIT64_WORKFLOW=(
    ROOT/".github"/"workflows"/"eqnedit64.yml"
)
EQNEDIT64_RELEASE_WORKFLOW=(
    ROOT/".github"/"workflows"/"release-eqnedit64-pypi.yml"
)
PACKAGE_CI_CONTEXTS=(
    (
        ROOT/".github"/"workflows"/"build-test.yml",
        "radia.ci-release-context.v1",
        "ci-release-context",
    ),
    (
        ROOT/".github"/"workflows"/"radia-mcp-matrix.yml",
        "radia-mcp.ci-release-context.v1",
        "radia-mcp-ci-release-context",
    ),
    (
        ROOT/".github"/"workflows"/"cubit-mesh-export.yml",
        "cubit-mesh-export.ci-release-context.v1",
        "cubit-mesh-export-ci-release-context",
    ),
)


def test_pypi_distributions_have_independent_ci_boundaries():
    workflows={
        "Radia": ROOT/".github"/"workflows"/"build-test.yml",
        "radia-mcp": ROOT/".github"/"workflows"/"radia-mcp-matrix.yml",
        "cubit-mesh-export": ROOT/".github"/"workflows"/"cubit-mesh-export.yml",
        "radia-optuna": ROOT/".github"/"workflows"/"radia-optuna.yml",
        "Eqnedit64": ROOT/".github"/"workflows"/"eqnedit64.yml",
    }
    texts={name: path.read_text(encoding="utf-8")
           for name, path in workflows.items()}
    for name, text in texts.items():
        assert f"name: {name}" in text

    radia=texts["Radia"]
    for package in ("eqnedit64", "radia-mcp", "cubit-mesh-export",
                    "radia-optuna"):
        assert f"packages/{package}/**" in radia
    assert "tests/test_release_quad_optuna_candidate.py" in radia
    assert "tools/release_quad.py" in radia
    assert ".agents/skills/release-eqnedit64/**" in radia
    assert "radia-mcp-wheel" not in radia
    assert "cubit-mesh-export-wheel" not in radia
    assert "radia-optuna-wheel" not in radia
    assert "--ignore=tests/test_cubit_installers.py" in radia
    assert "--ignore=tests/test_release_quad_optuna_candidate.py" in radia

    policy=(ROOT/".github"/"workflows"/"policy-lint.yml").read_text(
        encoding="utf-8")
    assert "packages/radia-mcp/tools/policy_lint.py" not in policy
    assert "packages/radia-mcp/tools/policy_lint.py" in texts["radia-mcp"]

    assert "tags: ['radia-mcp-v*']" in texts["radia-mcp"]
    assert "tags: ['cubit-mesh-export-v*']" in texts["cubit-mesh-export"]
    assert "tags: ['radia-optuna-v*']" in texts["radia-optuna"]
    assert "tags: ['eqnedit64-v*']" in texts["Eqnedit64"]


def test_ci_records_exact_ref_context_before_release():
    for path, schema, artifact in PACKAGE_CI_CONTEXTS:
        workflow=path.read_text(encoding="utf-8")
        assert schema in workflow
        assert "GITHUB_REF_TYPE" in workflow
        assert "GITHUB_SHA" in workflow
        assert "GITHUB_RUN_ID" in workflow
        assert f"name: {artifact}" in workflow


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
    assert '"matlab/optuna_mex.mexw64"' in workflow
    assert (
        '"radia_optuna/matlab/optuna_mex.mexw64": Path('
        in workflow
    )
    assert '"radia_optuna/matlab/LICENSE": Path("LICENSE")' in workflow
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


def test_eqnedit64_release_requires_exact_successful_tag_ci():
    ci=EQNEDIT64_WORKFLOW.read_text(encoding="utf-8")
    release=EQNEDIT64_RELEASE_WORKFLOW.read_text(encoding="utf-8")

    assert "Record exact Eqnedit64 CI ref context" in ci
    assert "eqnedit64.ci-release-context.v1" in ci
    assert "ref_type = $env:GITHUB_REF_TYPE" in ci
    assert "sha = $env:GITHUB_SHA" in ci
    assert "run_id = [string]$env:GITHUB_RUN_ID" in ci
    assert "name: eqnedit64-ci-release-context" in ci

    assert "workflow_run:" in release
    assert 'workflows: ["Eqnedit64"]' in release
    assert "github.event.workflow_run.conclusion == 'success'" in release
    assert "github.event.workflow_run.event == 'push'" in release
    assert "Download exact Eqnedit64 CI ref context" in release
    assert "run-id: ${{ github.event.workflow_run.id }}" in release
    assert "context.get(\"ref_type\") == \"tag\"" in release
    assert r"eqnedit64-v\d+\.\d+\.\d+" in release
    assert "ref: ${{ needs.qualify.outputs.sha }}" in release
    assert "runs-on: [self-hosted, windows-radia]" in release
    assert "$releaseRoots = @('C:\\Users\\Administrator\\OneDrive')" in release
    assert "Get-PSDrive -Name O -PSProvider FileSystem" in release
    assert "$releaseRoots = @('O:\\') + $releaseRoots" in release
    assert "$exe = Join-Path $releaseRoot 'Eqnedit64.exe'" in release
    assert "$manifestPath = Join-Path $releaseRoot 'Eqnedit64.release.json'" in release
    assert "The runner service uses LocalSystem" in release
    assert "eqnedit64.o-release.v1" in release
    assert "$manifest.source_sha -cne $env:EQNEDIT64_SHA" in release
    assert "name: eqnedit64-signed-standalone" in release
    assert "needs: [qualify, signed-standalone]" in release
    assert "Get-AuthenticodeSignature" in release
    assert "CN=ksugahar" in release
    assert "id-token: write" in release
    assert "pypa/gh-action-pypi-publish@release/v1" in release
    assert 'gh release create "$EQNEDIT64_TAG"' in release
    assert release.index("pypa/gh-action-pypi-publish@release/v1") < release.index(
        'gh release create "$EQNEDIT64_TAG"')


def test_eqnedit64_release_order_is_push_o_drive_then_tag():
    root = Path(__file__).resolve().parents[1]
    skill = (root / ".agents/skills/release-eqnedit64/SKILL.md").read_text(
        encoding="utf-8")
    sync = (
        root / ".agents/skills/release-eqnedit64/scripts/sync_to_o.ps1"
    ).read_text(encoding="utf-8")

    push_main = skill.index("Commit and push the versioned release source to `main`.")
    sync_o = skill.index("sync_to_o.ps1")
    push_tag = skill.index("Create the annotated tag")
    assert push_main < sync_o < push_tag

    assert "Release tag already exists; O: must be prepared before tag push" in sync
    assert "HEAD=$headSha origin/main=$originMainSha" in sync
    assert "eqnedit64.o-release.v1" in sync
    assert "source_sha = $normalizedSourceSha" in sync
