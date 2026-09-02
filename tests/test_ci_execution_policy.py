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
    assert "--ignore=tests/equation" in native
    assert '"MKLROOT=$mklRoot"' in native
    assert "MKLROOT=C:\\Program Files\\Python312\\Library" not in native
    assert "runs-on: [self-hosted, Windows, X64, mdx]" in optuna
    assert "windows-radia" not in optuna
    assert 'workflows: ["Radia Native Release"]' in release


def test_pre_push_runs_the_unpushed_candidate_on_mdx():
    hook = (ROOT / "tools" / "git-hooks" / "pre-push").read_text(encoding="utf-8")
    helper = (ROOT / "tools" / "ci_preflight_mdx.py").read_text(encoding="utf-8")

    assert "ci_preflight_mdx.py --base \"$remote_sha\" --head \"$local_sha\"" in hook
    assert '"bundle", "create"' in helper
    assert "scp" in helper
    assert "tools/ci_preflight.py --only policy,version" in helper
    assert "tests/test_docs_notebook_contract.py" in helper


def test_policy_twins_define_the_same_mdx_notebook_contract():
    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    claude = (ROOT / "CLAUDE.md").read_text(encoding="utf-8")
    heading = "### CI Execution, Validation Evidence, and Notebook Policy (2026-09-01)"

    def section(text: str) -> str:
        start = text.index(heading)
        return text[start:text.index("### ", start + len(heading))]

    policy = section(agents)
    assert policy == section(claude)
    normalized = " ".join(policy.split())
    assert "mdx gives CI and preflight work priority" in normalized
    assert "use hibino first when it is available" in normalized
    assert "the mdx CI queue is idle" in normalized
