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
    assert "expected = tag.removeprefix('cubit-mesh-export-v')" in workflow
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


def _cubit_workflow():
    import yaml
    return yaml.safe_load((ROOT / ".github/workflows/cubit-mesh-export.yml").read_text(encoding="utf-8"))


def _cubit_audit_program(phase):
    jobs = _cubit_workflow()["jobs"]
    job = jobs["build-test" if phase == "before-approval" else "publish"]
    step = next(s for s in job["steps"] if s.get("name") == "Verify exact release identity " + phase)
    return step["run"].split("<<'PY'\n", 1)[1].rsplit("\nPY", 1)[0]


def test_cubit_approval_separates_audit_from_publication():
    jobs = _cubit_workflow()["jobs"]
    audit, publish = jobs["build-test"], jobs["publish"]
    assert "environment" not in audit
    assert "id-token" not in audit.get("permissions", {})
    before = next(step for step in audit["steps"] if step.get("name") == "Verify exact release identity before-approval")
    assert publish["environment"] == "pypi"
    for unsafe in ("always(", "cancelled(", "failure("):
        assert unsafe not in publish["if"]
    assert not any(step.get("continue-on-error", False) for job in jobs.values() for step in job["steps"])
    checkout = next(step for step in publish["steps"] if step.get("uses", "").startswith("actions/checkout"))
    assert checkout["with"]["persist-credentials"] is False
    build_steps = audit["steps"]
    names = [step.get("name") for step in build_steps]
    assert names.index("Build and verify cubit-mesh-export wheel") + 1 == names.index("Record freshly built wheel bytes")
    assert names.index("Record freshly built wheel bytes") + 1 == names.index("Upload cubit-mesh-export wheel")
    assert names.index("Upload freshly built wheel identity") < names.index("Check standalone Python sources")
    assert names.index("Upload freshly built wheel identity") < names.index("Verify installed standalone GUI assets without Radia")
    assert before["env"]["WHEEL_DIR"] == "packages/cubit-mesh-export/dist"
    after = next(step for step in publish["steps"] if step.get("name") == "Verify exact release identity after-approval")
    assert after["env"]["WHEEL_DIR"] == "dist"
    assert publish["needs"] == ["ci-impact", "build-test"]
    for job in (before, publish):
        assert "github.ref_type == 'tag'" in job["if"]
        assert "github.event_name == 'push'" in job["if"]
        assert "github.event_name == 'workflow_dispatch'" in job["if"]
        assert "cubit-mesh-export-v" in job["if"]
    downloads = [s for s in publish["steps"] if s.get("uses", "").startswith("actions/download-artifact")]
    assert {s["with"]["artifact-ids"] for s in downloads} == {
        "${{ needs.build-test.outputs.wheel_artifact_id }}",
        "${{ needs.build-test.outputs.build_identity_artifact_id }}",
        "${{ needs.build-test.outputs.audit_artifact_id }}",
    }
    steps = publish["steps"]
    assert steps[-1]["uses"].startswith("pypa/gh-action-pypi-publish")
    assert "skip-existing" not in steps[-1].get("with", {})
    assert _cubit_audit_program("before-approval") == _cubit_audit_program("after-approval")
    assert jobs["build-test"]["outputs"]["wheel_artifact_id"] == "${{ steps.wheel-upload.outputs.artifact-id }}"


def _cubit_candidate(tmp_path, monkeypatch):
    import hashlib
    import importlib.util
    import json
    import shutil
    import subprocess
    import urllib.error
    import urllib.request
    import zipfile

    monkeypatch.chdir(tmp_path)
    package = tmp_path / "packages/cubit-mesh-export"
    source = package / "src/cubit_mesh_export"
    source.mkdir(parents=True)
    provenance_path = source / "_native_provenance.py"
    shutil.copyfile(ROOT / "packages/cubit-mesh-export/src/cubit_mesh_export/_native_provenance.py", provenance_path)
    spec = importlib.util.spec_from_file_location("fixture_provenance", provenance_path)
    provenance = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(provenance)
    native = tmp_path / "src/cubit_plugin"
    native.mkdir(parents=True)
    (native / "fixture.cpp").write_bytes(b"// release identity fixture\n")
    (tmp_path / "tools").mkdir()
    (tmp_path / "tools/find_cubit.ps1").write_bytes(b"# source identity fixture\n")
    digest, count = provenance.native_source_digest(tmp_path)
    payloads = {"cubit_mesh_export.ccm": b"MZ signed CCM", "cubit_mesh_curver.pyd": b"MZ signed PYD"}
    manifest = {"schema": provenance.SCHEMA, "source": {"commit": "b" * 40, "tree_sha256": digest, "file_count": count, "hash_format": provenance.HASH_FORMAT}, "payloads": {name: {"size": len(data), "sha256": hashlib.sha256(data).hexdigest()} for name, data in payloads.items()}}
    manifest_data = json.dumps(manifest).encode()
    (source / "native_payloads.json").write_bytes(manifest_data)
    (source / "cubit_mesh_export.ccm").write_bytes(payloads["cubit_mesh_export.ccm"])
    (source / "cubit_mesh_curver.pyd").write_bytes(payloads["cubit_mesh_curver.pyd"])
    (package / "pyproject.toml").write_text('[project]\nversion="2.2.7"\n', encoding="utf-8")
    wheel = package / "dist/cubit_mesh_export-2.2.7-cp312-cp312-win_amd64.whl"
    wheel.parent.mkdir()
    (source / "__init__.py").write_bytes(b'__version__ = "2.2.7"\n')
    with zipfile.ZipFile(wheel, "w") as z:
        for name, data in payloads.items():
            z.writestr("cubit_mesh_export/" + name, data)
        z.writestr("cubit_mesh_export/native_payloads.json", manifest_data)
        z.writestr("cubit_mesh_export/__init__.py", (source / "__init__.py").read_bytes())
        z.writestr("cubit_mesh_export/_native_provenance.py", provenance_path.read_bytes())
        z.writestr("cubit_mesh_export-2.2.7.dist-info/METADATA", "Name: cubit-mesh-export\nVersion: 2.2.7\n")
        for name in ("solver_ready_sample.jou", "cubit_export_menu.py", "cubit_toolbar/toolbars/cubit_mesh_export_toolbar.ttb.tmpl"):
            local = source / "cubit_gui" / name
            local.parent.mkdir(parents=True, exist_ok=True)
            local.write_bytes(b"fixture")
            z.writestr("cubit_mesh_export/cubit_gui/" + name, b"fixture")
    values = {"GITHUB_EVENT_NAME": "push", "GITHUB_REF_TYPE": "tag", "GITHUB_REF": "refs/tags/cubit-mesh-export-v2.2.7", "GITHUB_REF_NAME": "cubit-mesh-export-v2.2.7", "GITHUB_SHA": "a" * 40, "GITHUB_REPOSITORY": "fixture/repo", "GITHUB_RUN_ID": "12", "GITHUB_RUN_ATTEMPT": "1", "WHEEL_ARTIFACT_ID": "123"}
    for key, value in values.items():
        monkeypatch.setenv(key, value)
    remote = {"sha": "a" * 40}
    git_blobs = {p.relative_to(tmp_path).as_posix(): p.read_bytes() for p in tmp_path.rglob("*") if p.is_file() and "dist" not in p.parts and p.suffix != ".pyd"}
    def git(command, **kwargs):
        if command == ["git", "rev-parse", "HEAD"]:
            return "a" * 40
        if command[:2] == ["git", "ls-tree"]:
            prefixes = command[command.index("--") + 1:]
            return b"\0".join(name.encode() for name in sorted(git_blobs) if any(name == prefix or name.startswith(prefix + "/") for prefix in prefixes)) + b"\0"
        if command[:2] == ["git", "show"]:
            key = command[2].split(":", 1)[1]
            if key not in git_blobs:
                raise subprocess.CalledProcessError(128, command)
            return git_blobs[key]
        assert "ls-remote" in command, command
        return "c" * 40 + "\trefs/tags/cubit-mesh-export-v2.2.7\n" + remote["sha"] + "\trefs/tags/cubit-mesh-export-v2.2.7^{}\n"
    monkeypatch.setattr(subprocess, "check_output", git)
    def absent(url, **kwargs):
        assert url == "https://pypi.org/pypi/cubit-mesh-export/2.2.7/json"
        raise urllib.error.HTTPError(url, 404, "not found", {}, None)
    monkeypatch.setattr(urllib.request, "urlopen", absent)
    monkeypatch.setenv("WHEEL_DIR", "packages/cubit-mesh-export/dist")
    exec(compile(_cubit_build_identity_program(), "workflow-built-identity", "exec"), {"__name__": "__main__"})
    return wheel, source, native, remote


def _run_cubit_audit(phase, monkeypatch):
    import shutil
    monkeypatch.setenv("AUDIT_PHASE", phase)
    job = _cubit_workflow()["jobs"]["build-test" if phase == "before-approval" else "publish"]
    step = next(s for s in job["steps"] if s.get("name") == "Verify exact release identity " + phase)
    monkeypatch.setenv("WHEEL_DIR", step["env"]["WHEEL_DIR"])
    if phase == "after-approval":
        # Model the artifact download to publish's distinct directory.
        Path("dist").mkdir(exist_ok=True)
        for wheel in Path("packages/cubit-mesh-export/dist").glob("*.whl"):
            shutil.copyfile(wheel, Path("dist") / wheel.name)
    exec(compile(_cubit_audit_program(phase), "workflow-cubit-audit", "exec"), {"__name__": "__main__"})


def test_cubit_actual_audit_accepts_exact_candidate_before_and_after_approval(tmp_path, monkeypatch):
    import json
    _cubit_candidate(tmp_path, monkeypatch)
    _run_cubit_audit("before-approval", monkeypatch)
    _run_cubit_audit("after-approval", monkeypatch)
    receipt = json.loads((tmp_path / "audit/release-audit.json").read_text())
    assert receipt["artifact_id"] == "123"
    assert receipt["source_sha"] == "a" * 40


@pytest.mark.parametrize("change, message", [
    ("main", "exact Cubit tag ref"),
    ("event", "release event"),
    ("checkout", "checkout does not match"),
    ("remote", "remote tag moved"),
    ("version", "source package version"),
    ("native", "native source digest differs"),
    ("ccm", "CCM differs"),
    ("artifact", "artifact/ref identity changed"),
    ("wheel", "immediate post-build identity"),
    ("existing", "already exists"),
    ("attempt", "artifact/ref identity changed"),
])
def test_cubit_actual_audit_rejects_substitution_and_existing_version(tmp_path, monkeypatch, change, message):
    import json
    import subprocess
    import urllib.request
    import contextlib
    wheel, source, native, remote = _cubit_candidate(tmp_path, monkeypatch)
    _run_cubit_audit("before-approval", monkeypatch)
    if change == "main":
        monkeypatch.setenv("GITHUB_REF_TYPE", "branch")
    elif change == "event":
        monkeypatch.setenv("GITHUB_EVENT_NAME", "workflow_run")
    elif change == "checkout":
        original = subprocess.check_output
        monkeypatch.setattr(subprocess, "check_output", lambda cmd, **kw: "d" * 40 if cmd == ["git", "rev-parse", "HEAD"] else original(cmd, **kw))
    elif change == "remote":
        remote["sha"] = "d" * 40
    elif change == "version":
        (source.parents[1] / "pyproject.toml").write_text('[project]\nversion="2.2.8"\n', encoding="utf-8")
    elif change == "native":
        (native / "fixture.cpp").write_text("// altered\n", encoding="utf-8")
    elif change == "ccm":
        (source / "cubit_mesh_export.ccm").write_bytes(b"altered")
    elif change == "artifact":
        monkeypatch.setenv("WHEEL_ARTIFACT_ID", "456")
    elif change == "wheel":
        # Valid ZIP with identical entries but different whole-file bytes.
        wheel.write_bytes(wheel.read_bytes() + b"substitution")
    elif change == "attempt":
        monkeypatch.setenv("GITHUB_RUN_ATTEMPT", "2")
    elif change == "existing":
        monkeypatch.setattr(urllib.request, "urlopen", lambda *a, **kw: contextlib.nullcontext())
    with pytest.raises(SystemExit, match=message):
        _run_cubit_audit("after-approval", monkeypatch)


def _cubit_build_identity_program():
    step = next(s for s in _cubit_workflow()["jobs"]["build-test"]["steps"] if s.get("name") == "Record freshly built wheel bytes")
    return step["run"].split("<<'PY'\n", 1)[1].rsplit("\nPY", 1)[0]


@pytest.mark.parametrize("error", ["url", "timeout", "403", "500"])
def test_cubit_pypi_network_errors_never_create_an_audit(tmp_path, monkeypatch, error):
    import socket
    import urllib.request
    import urllib.error
    _cubit_candidate(tmp_path, monkeypatch)
    def broken(url, **kwargs):
        if error == "url":
            raise urllib.error.URLError("offline")
        if error == "timeout":
            raise socket.timeout("timeout")
        raise urllib.error.HTTPError(url, int(error), "unavailable", {}, None)
    monkeypatch.setattr(urllib.request, "urlopen", broken)
    with pytest.raises((urllib.error.URLError, TimeoutError)):
        _run_cubit_audit("before-approval", monkeypatch)
    assert not (tmp_path / "audit/release-audit.json").exists()


@pytest.mark.parametrize("phase", ["before-approval", "after-approval"])
def test_cubit_injected_python_fails_even_with_rewritten_build_identity(tmp_path, monkeypatch, phase):
    import zipfile
    wheel, source, native, remote = _cubit_candidate(tmp_path, monkeypatch)
    if phase == "after-approval":
        _run_cubit_audit("before-approval", monkeypatch)
    with zipfile.ZipFile(wheel) as z:
        contents = {name: z.read(name) for name in z.namelist()}
    contents["cubit_mesh_export/__init__.py"] += b"\nprint('injected')\n"
    with zipfile.ZipFile(wheel, "w") as z:
        for name, content in contents.items():
            z.writestr(name, content)
    # Deliberately model a forged matching same-job hash record. The independent
    # tagged Git blob must still reject the Python injection.
    monkeypatch.setenv("WHEEL_DIR", "packages/cubit-mesh-export/dist")
    exec(compile(_cubit_build_identity_program(), "forged-built-identity", "exec"), {"__name__": "__main__"})
    with pytest.raises(SystemExit, match="wheel member differs from tagged source"):
        _run_cubit_audit(phase, monkeypatch)


def test_cubit_real_build_layout_does_not_require_a_repository_root_dist(tmp_path, monkeypatch):
    _cubit_candidate(tmp_path, monkeypatch)
    assert not (tmp_path / "dist").exists()
    _run_cubit_audit("before-approval", monkeypatch)
    _run_cubit_audit("after-approval", monkeypatch)
    assert (tmp_path / "dist").is_dir()


def test_cubit_materializes_only_job_stage_from_tagged_blobs(tmp_path, monkeypatch):
    _cubit_candidate(tmp_path, monkeypatch)
    source = tmp_path / "packages/cubit-mesh-export/src/cubit_mesh_export"
    init = source / "__init__.py"
    tagged = init.read_bytes()
    # Model Windows checkout conversion; materialization must not rewrite it.
    init.write_bytes(tagged.replace(b"\n", b"\r\n"))
    existing = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    runner = tmp_path / "runner-temp"
    runner.mkdir()
    monkeypatch.setenv("RUNNER_TEMP", str(runner))
    export = tmp_path / "github-env"
    monkeypatch.setenv("GITHUB_ENV", str(export))
    step = next(s for s in _cubit_workflow()["jobs"]["build-test"]["steps"] if s.get("name") == "Materialize exact tagged package bytes")
    program = step["run"].split("<<'PY'\n", 1)[1].rsplit("\nPY", 1)[0]
    exec(compile(program, "workflow-stage-source", "exec"), {"__name__": "__main__"})
    stage = runner / "cubit-wheel-source"
    assert (stage / "packages/cubit-mesh-export/src/cubit_mesh_export/__init__.py").read_bytes() == tagged
    assert (stage / "src/cubit_plugin/fixture.cpp").read_bytes() == b"// release identity fixture\n"
    assert (stage / "tools/find_cubit.ps1").read_bytes() == b"# source identity fixture\n"
    assert (stage / "packages/cubit-mesh-export/src/cubit_mesh_export/cubit_mesh_curver.pyd").read_bytes() == (source / "cubit_mesh_curver.pyd").read_bytes()
    assert export.read_text(encoding="utf-8").strip() == "CUBIT_BUILD_ROOT=" + str(stage)
    assert all(path.read_bytes() == content for path, content in existing.items())
