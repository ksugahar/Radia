"""Guards: CI host routing, native release and pre-push boundaries remain explicit.

Assumes: publication is dispatched on an approved release tag, independently of
native build completion; runner capabilities are declared in workflow sources.
"""

import importlib.util
import re
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]


def _preflight_module():
    source = ROOT / "tools" / "ci_preflight.py"
    spec = importlib.util.spec_from_file_location("ci_preflight", source)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fast_ci_runs_only_on_100_and_native_is_a_named_release_lane():
    fast = (ROOT / ".github" / "workflows" / "radia.yml").read_text(encoding="utf-8")
    native = (ROOT / ".github" / "workflows" / "build-test.yml").read_text(encoding="utf-8")
    optuna = (ROOT / ".github" / "workflows" / "radia-optuna.yml").read_text(encoding="utf-8")
    release = (ROOT / ".github" / "workflows" / "release.yml").read_text(encoding="utf-8")

    assert "name: Radia\n" in fast
    import yaml
    labels = yaml.safe_load(fast)['jobs']['fast-contracts']['runs-on']
    assert labels == ['self-hosted', 'Windows', 'X64', 'ci-light']
    assert '100' not in labels
    assert "Build with MSVC" not in fast
    assert "validation_test/" not in fast.split("jobs:", 1)[1]
    assert "name: Radia Native Release" in native
    assert 'runs-on: [self-hosted, Windows, X64, "100"]' in native
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
    assert 'runs-on: [self-hosted, Windows, X64, "100", gpu-rtx]' in optuna
    assert "windows-radia" not in optuna
    optuna_self_hosted = optuna.split("  build-test:", 1)[1].split("  installed-wheel-matlab-e2e:", 1)[0]
    assert "actions/setup-python" not in optuna_self_hosted
    assert "RADIA_OPTUNA_CI_PYTHON" in optuna_self_hosted
    assert "-m venv" in optuna_self_hosted
    assert "pip install pybind11" not in optuna_self_hosted
    assert "Build and verify Python kernels before linking plain MEX" in optuna_self_hosted
    assert optuna.count("'.github/workflows/radia-optuna.yml'") == 2
    release_document = yaml.safe_load(release)
    assert set(release_document.get('on', release_document.get(True))) == {'workflow_dispatch'}
    assert "github.ref_type == 'tag'" in release_document['jobs']['verify-promotion']['if']
    assert 'tools/verify_radia_promotion.py' in release


def test_fast_ci_declares_bh_and_wheel_contract_dependencies():
    import yaml

    workflow = yaml.safe_load((ROOT / '.github/workflows/radia.yml').read_text(encoding='utf-8'))
    steps = workflow['jobs']['fast-contracts']['steps']
    setup = next(step['run'] for step in steps
                 if step.get('name') == 'Create isolated fast-CI environment')
    installs = [line.split() for line in setup.splitlines() if '-m pip install ' in line]
    assert any({'numpy', 'scipy', 'pytest', 'wheel'}.issubset(tokens) for tokens in installs)
    assert 'python -m venv' in setup
    assert 'ngsolve' not in setup.lower()


def test_no_workflow_selects_a_retired_runner():
    """No lane may ask for a retired or user-desktop runner, signing included.

    Only active self-hosted labels may appear.
    """
    workflows = ROOT / ".github" / "workflows"
    retired_labels = ("windows-radia", "lab")

    for path in sorted(workflows.glob("*.yml")):
        runner_lines = [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip().startswith("runs-on:")
        ]
        for line in runner_lines:
            normalized = line.lower()
            for label in retired_labels:
                assert label not in normalized, f"{path.name}: {line}"

    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    claude = (ROOT / "CLAUDE.md").read_text(encoding="utf-8")
    for policy in (agents, claude):
        assert "No workflow selects a retired or user-desktop runner." in policy


def test_policy_lint_enforces_the_retired_examples_boundary():
    policy = (ROOT / "tools" / "policy_lint.py").read_text(encoding="utf-8")

    assert '"src/radia/*.py", "src/radia/**/*.py"' in policy
    # Policy 7 lists the retired tier through the failing-loud ls-files
    # helper, so a broken git fails the check instead of passing it.
    assert '_git_ls_files("examples")' in policy
    assert "README.md in every example dir" not in policy


def test_distribution_ci_is_change_scoped_and_mcp_full_suite_is_explicit():
    fast = (ROOT / ".github" / "workflows" / "radia.yml").read_text(
        encoding="utf-8"
    )
    mcp = (ROOT / ".github" / "workflows" / "radia-mcp-matrix.yml").read_text(
        encoding="utf-8"
    )

    for independent_path in (
        "packages/radia-mcp/**",
        "packages/cubit-mesh-export/**",
        "packages/radia-optuna/**",
    ):
        assert independent_path not in fast

    # Universal Policy Lint triggers are owned by the workflow-owner contract
    # in tests/policy/test_workflow_owner_scope.py.
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
    assert 'env["RADIA_MCP_CI_SELECTION_JSON"]' in mcp
    # Directory discovery must run optional-dependency collection hooks;
    # RADIA_MCP_CI_SELECTION_JSON still restricts the collected test items.
    # The package's workflow-execution test checks exact file/node selection.
    assert 'targets = ["tests"]' in mcp
    assert '"-m", "not slow", *targets,' in mcp
    assert '"--confcutdir=tests/mcp_integration", *targets' in mcp
    assert 'metadata["project"]["optional-dependencies"]["md2html"]' in mcp
    assert 'check_package_test_boundary.py' in mcp
    assert 'plan["server_selftests"]' in mcp
    assert 'plan["run_mcp_response_tests"]' in mcp
    assert "Meta health (all cataloged subpackages must import)" not in mcp
    assert "radia_mcp_health" not in mcp
    assert "pytest-rerunfailures" not in mcp
    assert "--reruns" not in mcp
    assert "check_tools_doc" not in mcp
    assert "TOOLS.md drift gate" not in mcp
    assert "mode=release" in mcp
    assert "needs.scope.outputs.build_wheel == 'true'" in mcp


def test_preflight_uses_live_mcp_contracts_without_generated_inventory_gate():
    preflight = _preflight_module()

    source_keys = {
        key
        for key, _, _ in preflight._gates_for_changes(
            ["packages/radia-mcp/src/radia_mcp/matlab/runtime.py"]
        )
    }
    assert source_keys == {"policy", "publish-boundary", "version", "radia-mcp"}

    docs_keys = {
        key
        for key, _, _ in preflight._gates_for_changes(
            ["packages/radia-mcp/docs/architecture.md"]
        )
    }
    assert docs_keys == {"policy", "version"}

    generator_keys = {
        key
        for key, _, _ in preflight._gates_for_changes(
            ["packages/radia-mcp/scripts/gen_tools_doc.py"]
        )
    }
    assert generator_keys == {
        "policy",
        "publish-boundary",
        "version",
        "radia-mcp",
    }

    assert not hasattr(preflight, "gate_tools_md")
    gitignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
    assert "/packages/radia-mcp/docs/TOOLS.md" in gitignore
    assert not (ROOT / "packages" / "radia-mcp" / ".gitattributes").exists()


def test_preflight_catches_taskmanager_failures_before_fast_tests(monkeypatch):
    preflight = _preflight_module()
    for path in (
        "tests/test_panel_impedance_file.py",
        "validation_test/induction_heating/test_panel_esim_workflow.py",
        "src/radia/surface_impedance.py",
        "docs/solver_probe.py",
        "tools/audit_taskmanager.py",
    ):
        keys = {key for key, _, _ in preflight._gates_for_changes([path])}
        assert "taskmanager" in keys, path

    calls = []

    def failed_audit(command):
        calls.append(command)
        return 1, "tests/test_panel_impedance_file.py: caller lacks TaskManager"

    monkeypatch.setattr(preflight, "_sh", failed_audit)
    ok, detail = preflight.gate_taskmanager_audit()
    assert not ok
    assert "caller lacks TaskManager" in detail
    assert Path(calls[0][-1]).name == "audit_taskmanager.py"
    monkeypatch.setattr(preflight, "_changed_since", lambda ref: ["tests/test_panel_impedance_file.py"])
    assert preflight.main(["--only", "taskmanager"]) == 1


def test_active_docs_do_not_restore_generated_catalog_or_old_release_name():
    active_docs = [
        ROOT / "README.md",
        ROOT / "AGENTS.md",
        ROOT / "CLAUDE.md",
        ROOT / ".agents" / "skills" / "release-quad" / "SKILL.md",
        ROOT / ".agents" / "skills" / "verify-deploy" / "SKILL.md",
    ]

    for path in active_docs:
        text = path.read_text(encoding="utf-8")
        assert "[Generated MCP tool catalog]" not in text, path
        assert "`release-qud`" not in text, path

    release_skill = active_docs[3].read_text(encoding="utf-8")
    assert "do not restore a tracked generated `docs/TOOLS.md` gate" in release_skill


def test_solver_release_quad_excludes_cubit_generation():
    release_quad = (ROOT / "tools" / "release_quad.py").read_text(
        encoding="utf-8"
    )
    release_skill = (ROOT / ".agents" / "skills" / "release-quad" / "SKILL.md").read_text(
        encoding="utf-8"
    )
    assert "cubit-dual" not in release_skill
    assert "cubit-dual" not in release_quad
    assert "phase0" not in release_quad
    assert "cubit-plugin-install" not in release_quad
    assert "Stop-Process" not in release_quad
    assert "Never change radia-mcp or" in release_quad


def test_preflight_git_scope_is_safe_in_an_isolated_worktree(monkeypatch):
    preflight = _preflight_module()
    commands = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0, stdout="README.md\n", stderr="")

    monkeypatch.setattr(preflight.subprocess, "run", fake_run)

    assert preflight._changed_since("HEAD^") == ["README.md"]
    assert len(commands) == 4
    for command in commands:
        assert command[:3] == [
            "git",
            "-c",
            f"safe.directory={preflight.REPO}",
        ]


def test_pre_push_runs_the_unpushed_candidate_on_lab():
    hook = (ROOT / "tools" / "git-hooks" / "pre-push").read_text(encoding="utf-8")
    helper = (ROOT / "tools" / "ci_preflight_remote.py").read_text(encoding="utf-8")

    dispatcher = (ROOT / "tools" / "pre_push_check.py").read_text(encoding="utf-8")
    assert "python tools/pre_push_check.py" in hook
    assert '"--base", base, "--head", head' in dispatcher
    assert '"bundle", "create"' in helper
    assert 'git_output("merge-base", head, main_ref)' in helper
    assert "$base = '{effective_base}'" in helper
    assert "scp" in helper
    assert "tools/ci_preflight.py --since $base" in helper
    assert "'mcp>=1.0,<2'" in helper
    assert "numpy scipy pytest pyyaml setuptools wheel 'mcp>=1.0,<2'" in helper
    assert 'import mcp, numpy, pytest, scipy, setuptools, wheel, yaml' in helper
    assert "tools/run_test_tier.py --profile fast-contracts" in helper
    assert "candidate preflight helper is missing" in dispatcher
    assert '&& [ -f tools/ci_preflight_remote.py ]' not in hook
    assert "upload_release_asset.py" not in hook
    assert "developer push must never upload mutable" in hook


def test_policy_twins_define_the_same_ci_notebook_contract():
    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    claude = (ROOT / "CLAUDE.md").read_text(encoding="utf-8")
    heading = "### CI Execution, Validation Evidence, and Notebook Policy (2026-09-03)"

    def section(text: str) -> str:
        start = text.index(heading)
        return text[start:text.index("### ", start + len(heading))]

    policy = section(agents)
    assert policy == section(claude)
    normalized = " ".join(policy.split())
    assert "with the `100` runner label" in normalized
    # Hosts appear by role only; the role-to-host map is private configuration.
    assert "radia uses bound editable sources on build and validation roles" in normalized
    assert "isolated exact-wheel QA" in normalized
    assert "radia-mcp stays build-editable/validation-wheel" in normalized
    assert "cubit-mesh-export is build-only" in normalized
    assert "Use GPU workers for declared GPU builds and qualification" in normalized
    assert "retired targets stay excluded" in normalized
    assert "validation has no Cubit/cme" in normalized
    assert "NGSolve starts from checked exported .vol files" in normalized
    assert "CI scope begins at the independently released distribution boundary" in normalized
    assert "each server's live `tools/list` response are the tool-discovery source of truth" in normalized
    assert "Generated tool inventory snapshots are local diagnostics" in normalized
    assert "not committed artifacts or CI oracles" in normalized
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
    assert "Pre-push runs only impact-scoped build-target preflight" in normalized


def test_policy_twins_define_the_same_compute_host_routing():
    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    claude = (ROOT / "CLAUDE.md").read_text(encoding="utf-8")
    # Anchor on the heading, not on a dated POLICY line: the routing rule is
    # revised whenever the host fleet changes and the marker should survive it.
    start_marker = "### Compute Host Routing"
    end_marker = "**POLICY**: Every benchmark script saves a machine-readable JSON result."

    def routing_policy(text: str) -> str:
        start = text.index(start_marker)
        return text[start:text.index(end_marker, start)]

    policy = routing_policy(agents)
    assert policy == routing_policy(claude)
    normalized = " ".join(policy.split())
    # Roles, not hosts: probe commands, aliases and core ranges stay private.
    assert "Build owns CPU/secondary GPU development, CI/formal artifacts" in normalized
    assert "Validation is primary GPU development in isolated checkouts/venvs, separate from exact-wheel QA/runtime; no Cubit/cme" in normalized
    assert "JMAG evaluation runs no Radia" in normalized
    assert "Compare CPU/GPU on the same host" in normalized
    assert "Timings are diagnostic" in normalized
    assert "one heavy job at a time" in normalized.lower()
    assert "running jobs and available memory" in normalized
    assert "never stop another session's jobs" in normalized
    assert "build-host CPU <=8 physical cores/job and <=6 host-wide is unconditional" in normalized
    assert "no admission, owner-pool, idle-core or load gates" in normalized
    assert "Count all CI/build/validation/research trees, not compiler children; <=6 TOTAL" in normalized
    assert "GPU jobs use a separate slot, not a CPU slot" in normalized
    assert "record host CPU use, avoid oversubscription" in normalized
    assert "overflow to a free build GPU slot is pre-authorized" in normalized
    assert "<=1 GPU job/host, <=2 across both" in normalized
    assert "Keep artifact/architecture/ABI/runtime/owner scope" in normalized
    assert "no install/config or unverified substitution" in normalized
    assert "Never stop, repin or overwrite others' jobs/runtimes" in normalized
    assert "Outside-frame CPU and other remote, environment, deployment/publication actions retain authorization" in normalized
    assert "Use an isolated venv; preserve other/shared runtimes" in normalized
    assert "BelowNormal, owned no-breakaway containment and absolute deadlines" in normalized
    assert "Memory caps need specific justification" in normalized
    assert "`owner.json` records owner/source/runtime/start/deadline/job count/resources" in normalized
    assert "never PR gates" in normalized
    assert "Measure speed" in normalized
    assert "Historical timings retain provenance" in normalized
    assert "Keep host maps/originals/hashes/IPs/paths private; never relabel hosts" in normalized


def test_agent_policies_stay_compact_and_share_one_policy_body():
    agents_path = ROOT / "AGENTS.md"
    claude_path = ROOT / "CLAUDE.md"
    agents = agents_path.read_text(encoding="utf-8").splitlines()
    claude = claude_path.read_text(encoding="utf-8").splitlines()

    assert agents[0] == "# Codex - Radia Project Policy"
    assert claude[0] == "# Claude Code - Radia Project Policy"
    assert agents[1:] == claude[1:]
    assert len(agents) < 260
    # The public policy must stay small enough for an agent's context; private
    # session coordination lives in CLAUDE.local.md and user-level Codex rules.
    assert len("\n".join(agents).encode("utf-8")) < 24_000
    assert len("\n".join(claude).encode("utf-8")) < 24_000

    policy = "\n".join(agents)
    # Public text names roles, never physical hosts or private paths.
    for identifier in ("号機", "A100 host", "`packages/radia-gpu`"):
        assert identifier not in policy, identifier
    assert not re.search(r"\b\d{1,3}(?:\.\d{1,3}){3}\b", policy), "IPv4 address in public policy"
    assert "### Tests For AI-Authored Code" in policy
    assert "Investigations remain in Git" in policy
    assert "Never include, link, or use NGSolve TaskManager in MEX builds" in policy
    assert "mexcuda is retired; no fallback" in policy
    assert "link the same objects/library from plain mex/matlab_add_mex" in policy
    assert all("mexcuda" not in line.lower() or "retired" in line.lower() for line in agents)
    assert "NGSolve-inherited DLLs" in policy
    assert "CPU numerical default: 1 thread for independent optimization runs" in policy
    assert "multithreading supported, explicit 8-physical-core large-solve profile" in policy
    assert "GPU-first; measure speed per case" in policy
    assert "Log TaskManager/OpenMP/MKL limits" in policy
    assert "GPU uses NGSolve device/NVRTC, no TaskManager extension/fallback" in policy
    assert "Text: UTF-8/LF" in policy
    assert "Keep management, sequential RNG" in policy
    assert "Do not expand this file into a second manual" in policy
    assert "packages/radia-mcp/src/radia_mcp/**/knowledge/" in policy

    # Runtime ownership must survive both agents' policy updates. This checks
    # the documentation contract, not an unimplemented deployment mutex.
    assert "### Shared MCP Runtime Ownership" in policy
    normalized = " ".join(policy.split())
    assert "Developers may edit live MCP source and change its editable source" in normalized
    assert "no dedicated branch, frozen snapshot or separate deployment approval" in normalized
    assert "live loaded source separately" in normalized
    assert "Unknown or mixed evidence must stay unverified" in normalized
    assert "Manuscript folders contain generated `.bbl` only, never local `.bib` copies" in normalized
    assert "bibliography/data/references.bib" in normalized
    assert "commit-pinned source snapshot" not in normalized
    contract_path = "packages/radia-mcp/docs/operations/mcp-runtime-policy.md"
    assert contract_path in policy
    contract = (ROOT / contract_path).read_text(encoding="utf-8")
    for required in (
        "change its editable installation source",
        "no dedicated branch, immutable snapshot or separate",
        "mixed-generation",
        "registration-time",
        "Never mass-kill",
    ):
        assert required in contract


def test_ngsolve_mcp_install_guidance_matches_native_dependency_pin():
    knowledge = (
        ROOT
        / "packages"
        / "radia-mcp"
        / "src"
        / "radia_mcp"
        / "radia_ngsolve"
        / "knowledge"
        / "ngsolve.py"
    ).read_text(encoding="utf-8")

    pins = re.findall(r'"ngsolve==([0-9.]+)"', (ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert len(pins) == 1, pins
    pin = pins[0]
    assert f"ngsolve=={pin} netgen-mesher=={pin}" in knowledge, f"knowledge install line does not match pyproject pin {pin}"
    assert f"NGSolve {pin} uses ngsolve-openblas" in knowledge
    assert "Version 6.2.2604 includes: MKL" not in knowledge


def test_workflow_discovers_directory_and_preserves_exact_selectors(monkeypatch, tmp_path):
    import json
    import subprocess
    import textwrap
    from types import SimpleNamespace
    import pytest

    workflow = (ROOT / '.github/workflows/radia-mcp-matrix.yml').read_text(encoding='utf-8')
    step = workflow.split('      - name: Pytest (compact contracts plus affected package tests)', 1)[1]
    script = textwrap.dedent(step.split("          python - <<'PY'\n", 1)[1].split('          PY\n', 1)[0])
    selectors = ['tests/test_ac_internal_inductance.py', 'tests/test_meta_health.py::test_meta_catalog_has_at_least_30_servers']
    (tmp_path / 'radia-mcp-ci-selection.json').write_text(json.dumps({'package_tests': selectors}), encoding='utf-8')
    monkeypatch.chdir(tmp_path)
    calls = []
    def run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(subprocess, 'run', run)
    with pytest.raises(SystemExit) as stopped:
        exec(compile(script, '<workflow pytest step>', 'exec'), {})
    assert stopped.value.code == 0
    command, kwargs = calls[0]
    assert command[-1] == 'tests'
    assert not any(selector in command for selector in selectors)
    assert json.loads(kwargs['env']['RADIA_MCP_CI_SELECTION_JSON']) == selectors


def test_host_roles_agree_across_deployment_skills():
    root = ROOT
    cubit = " ".join((root / ".agents/skills/release-cubit-mesh-export/SKILL.md").read_text(encoding="utf-8").split())
    # The public role table documents CLI identifiers, not physical host names.
    # Bind its tiers to the actual deployment source instead of historical prose.
    guide = root / ".agents/skills/release-cubit-mesh-export/SKILL.md"
    documented_tiers = {}
    for line in guide.read_text(encoding="utf-8").splitlines():
        cells = [cell.strip().strip("`") for cell in line.strip("|").split("|")]
        if len(cells) == 3 and cells[2] in {"editable", "wheel"}:
            assert cells[1] not in documented_tiers, "duplicate deployment target"
            documented_tiers[cells[1]] = cells[2]
    source = root / "tools/release_cubit_dual.py"
    spec = importlib.util.spec_from_file_location("cubit_dual_policy_contract", source)
    assert spec and spec.loader
    dual = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dual)
    assert documented_tiers == dual.TIERS == {"100": "editable"}
    assert set(documented_tiers) == set(dual.TARGETS)
    assert "`tools/release_cubit_dual.py`" in cubit
    # Both release agents use the same portable acceptance/cleanup contract.
    import re
    links = re.findall(r"\]\((\.\./[^)]+)\)", guide.read_text(encoding="utf-8"))
    assert any((guide.parent / link).resolve() ==
               (root / ".agents/skills/cubit-toolbar-acceptance/SKILL.md").resolve()
               and (guide.parent / link).is_file() for link in links)


def _text_eol_module():
    spec = importlib.util.spec_from_file_location(
        "tracked_text_eol_contract", ROOT / "tools/check_text_line_endings.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_tracked_text_index_is_lf_and_attributes_preserve_evidence():
    """Guards: CRLF/mixed source blobs fail CI; binary/evidence bytes survive.

    Assumes: Git index is the publication input; working-copy autocrlf is not
    a reason to rewrite sealed evidence or another owner's files.
    """
    attributes = (ROOT / '.gitattributes').read_text(encoding='utf8')
    assert '* text=auto eol=lf' in attributes
    assert 'validation_records/** -text' in attributes
    result = _text_eol_module().check_index(ROOT)
    assert result['passed'], result['violations']


def test_text_eol_guard_rejects_crlf_and_mixed_including_cmd(tmp_path):
    """Guards: the index guard actually fails and respects binary/-text overrides.

    Assumes: synthetic Git blobs bypass checkout normalization deliberately;
    fixtures here are not solver baselines or frozen historical evidence.
    """
    import subprocess
    subprocess.run(['git', 'init', '-q', str(tmp_path)], check=True)
    def index(name, payload):
        oid = subprocess.run(['git', '-C', str(tmp_path), 'hash-object', '-w', '--stdin'],
                             input=payload, capture_output=True, check=True).stdout.decode().strip()
        subprocess.run(['git', '-C', str(tmp_path), 'update-index', '--add',
                        '--cacheinfo', f'100644,{oid},{name}'], check=True)
    index('.gitattributes', b'* text=auto eol=lf\n*.bin binary\nsealed.json -text\n')
    # Attributes are also present in the checkout, as for a normal CI checkout.
    (tmp_path / '.gitattributes').write_bytes(
        b'* text=auto eol=lf\n*.bin binary\nsealed.json -text\n')
    index('good.py', b'print("ok")\n')
    index('binary.bin', b'\0\r\n')
    index('sealed.json', b'{"historical": true}\r\n')
    guard = _text_eol_module()
    assert guard.check_index(tmp_path)['passed']
    for name, content in [('bad.py', b'x=1\r\n'), ('bad.cmd', b'echo bad\r\n'),
                          ('mixed.md', b'one\ntwo\r\n')]:
        index(name, content)
    result = guard.check_index(tmp_path)
    assert not result['passed']
    assert {v['path'] for v in result['violations']} == {'bad.py', 'bad.cmd', 'mixed.md'}
