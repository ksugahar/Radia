import re

from radia_mcp.radia_ngsolve.knowledge.install_deploy import (
    INSTALL_DEPLOY,
    get_install_deploy_documentation,
)
from radia_mcp.radia_ngsolve.knowledge.radia import RADIA_BUILD_AND_RELEASE
from radia_mcp.radia_ngsolve.knowledge.release_workflow import (
    RELEASE_WORKFLOW,
    get_release_workflow_documentation,
)


def test_install_deploy_topics_match_current_machine_roles():
    overview = get_install_deploy_documentation("overview")
    ci_compute = get_install_deploy_documentation("ci_compute")

    assert "CPU-only: no CUDA/nvcc/mexcuda, no NVIDIA/cuDSS packages, no GPU execution" in overview
    assert "Windows GPU validation host" in overview and "`gpu-a100`" in overview
    assert "GPU work never overflows to the build host" in " ".join(overview.split())
    assert "the only editable development host" in overview
    assert "verified artifacts" in overview
    for retired in (r"\bLAB\b", "md" + "x", "hi" + "bino", "on-demand GPU worker"):
        assert not re.search(retired, INSTALL_DEPLOY), retired
    assert "editable `radia`, `cubit-mesh-export`, `radia-mcp`" not in overview
    assert "impact-selected" in ci_compute
    assert "validation_test/" in ci_compute


def test_install_deploy_does_not_publish_retired_binary_copy_routes():
    assert "push_pyds_to_mdx" not in INSTALL_DEPLOY
    assert "download_binaries" not in INSTALL_DEPLOY
    assert "manually dropping" in INSTALL_DEPLOY
    assert "cubit-plugin-install" in get_install_deploy_documentation("cubit")


def test_build_release_manual_matches_the_isolated_mdx_pipeline():
    text = RADIA_BUILD_AND_RELEASE

    assert "run-local virtual environment" in text
    assert "tools/ci_preflight.py" in text
    assert "tools/release_quad.py" in text
    assert "Build_Wheel.ps1 -DryRun" in text
    assert "mkl>=2026,<2027" in text
    assert "NGSolve must be copied" not in text
    assert "robocopy" not in text
    assert "pytest -m basic" not in text
    assert "No GitHub CLI policy" not in text
    assert "binaries` release" not in text


def test_mcp_release_route_is_distinct_from_solver_release():
    mcp = get_install_deploy_documentation("mcp_release")
    solver = get_install_deploy_documentation("release")
    assert mcp.startswith("## mcp_release\n")
    assert "radia-mcp-v<VERSION>" in mcp
    assert "validation-host live source" in mcp
    assert "next-launch-pending" in mcp
    assert "a verified\nwheel and fresh import on the validation host" in mcp
    assert "editable development install\nand fresh import on the build host" in mcp
    assert "do not block release completion" in mcp
    assert "Failed\ninstallation/import still blocks" in mcp
    assert "tools/release_quad.py" not in mcp
    assert "tools/release_quad.py" in solver
    assert "topic `mcp_release`" in solver


def test_release_workflow_does_not_restore_coupled_mcp_policy():
    mcp = get_release_workflow_documentation("mcp_release")
    quality = get_release_workflow_documentation("mcp_quality_review")
    overview = get_release_workflow_documentation("overview")
    assert "## mcp_release" in mcp
    assert "next-launch-pending" in mcp
    assert "does not block\nrelease completion" in mcp
    assert "solver\nrelease installer" in mcp
    assert "release-quad Phase 8/9 checks" not in quality
    assert "validation-host live-source" in quality
    assert "radia-ngsolve, cubit, build123d" not in RELEASE_WORKFLOW
    assert "Cubit MCP is distributed by cubit-mesh-export" in overview
    assert "CPU-only editable development install" in overview
    assert "fixed release checkout" in overview
    assert "release_native_payloads.json" in overview
    assert "installed NGSolve/Netgen ABI" in overview
    assert "Bump 4 version files" not in RELEASE_WORKFLOW
    assert "Composite commit (HEREDOC, all packages" not in RELEASE_WORKFLOW
    assert "the radia-mcp Scripts/mcp-server-*.exe" not in RELEASE_WORKFLOW
    assert "Never kill every MCP" in RELEASE_WORKFLOW
    assert "Rows\ndo not have to share the same version" in RELEASE_WORKFLOW
    assert "never restore an older tree" in RELEASE_WORKFLOW
    assert "three-package reinstall" not in RELEASE_WORKFLOW
    assert "RADIA_RELEASE_PRESERVE_MCP_SOURCE" not in RELEASE_WORKFLOW
    assert "solver candidate lane never installs" in RELEASE_WORKFLOW
    assert "monorepo_lockstep" not in RELEASE_WORKFLOW
    assert "are the validation host (verified wheel) and the build host (editable) only" in RELEASE_WORKFLOW
    for lane in ("`cpu`", "`linux-python-gpu`", "`windows-python-gpu`", "`matlab`"):
        assert lane in RELEASE_WORKFLOW, lane
    for retired in (r"\bLAB\b", "md" + "x", "(?i)hi" + "bino", "QUAD", "four-machine", "simulink-candidate --"):
        assert not re.search(retired, RELEASE_WORKFLOW), retired
    assert "## version_pairs" in get_release_workflow_documentation("version_pairs")


def test_host_roles_agree_across_package_guidance_and_recovery():
    from radia_mcp.meta.bug_patterns import PATTERNS

    install = " ".join(INSTALL_DEPLOY.split())
    completion = " ".join(get_release_workflow_documentation("mcp_quality_review").split())
    assert "Only the build host uses editable installs" in install
    assert "verified wheel to the validation host" in install
    assert "verified non-editable wheel and fresh import on the validation host" in completion
    assert "editable development registration and fresh import on the build host" in completion
    incident = next(p for p in PATTERNS if p["id"] == "lab-editable-drift-after-pip-force-reinstall")
    assert "LAB-editable" not in incident["detection"]
    assert "On LAB, repair with the verified current wheel, never an editable pointer" in incident["prevention"]
    assert "LAB wheel / 100 editable release-dual" in incident["prevention"]
    assert "Source edits no longer flow to runtime" in incident["what"]
