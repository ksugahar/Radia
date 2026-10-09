#!/usr/bin/env python
"""The Radia source policies as ONE reusable module.

Single source of truth for BOTH:
  - .github/workflows/policy-lint.yml   (the CI gate)
  - tools/ci_preflight.py               (the local pre-push gate)

so the two can never drift (previously ci_preflight re-implemented the
policies inline -- if the workflow gained a policy, the local gate would
silently miss it).  Mirrors the historical inline bash greps; uses
`git grep` / `git ls-files` so it sees TRACKED working-tree content (what
CI checks out on a fresh runner).

    python tools/policy_lint.py        # run them all; exit 1 on any violation
    python tools/policy_lint.py --quiet

Policies (see CLAUDE.md):
  1 no FldUnits() in Radia source        5 no generated files at repo root
  2 no tracked binaries                  6 no legacy src/python import path
  3 no Helmholtz WAVE kernel in core     7 no tracked retired examples/ tier
  4 no CblasColMajor in core (allowlist)
  8 HDiv geometry/field pairs use the central capability table
  9 no scrubbed commercial-solver brand in code or stored results
"""
from __future__ import annotations

import os
import json
import re
import subprocess
import sys

# UNC-safe repo root (see tools/ci_preflight.py for the rationale: never
# resolve() on the LAB S: drive -- it canonicalises to a UNC form Python
# file reads reject).
_THIS = os.path.abspath(__file__)
REPO = os.path.dirname(os.path.dirname(_THIS))
_norm = REPO.replace("\\", "/")
if "192.168.11.100" in _norm and "/Radia/01_GitHub" in _norm:
    REPO = "S:" + _norm[_norm.index("/Radia/01_GitHub"):]

# Policy 9: brands whose scrub is complete, so a reappearance is a regression
# rather than unfinished work.  COMSOL, FEMM, JMAG, CST and ELF are NOT here
# yet -- they still appear widely in code and join this tuple as each one's
# scrub lands.  The formulation those names once carried is Simkin and
# Trowbridge's mixed total/reduced scalar potential (``simkin1980three``).
SCRUBBED_BRANDS = ("TOSCA",)

# Policy 10: solvers whose MCP servers must never be published, whatever the
# repository they live in.  This is a stricter rule than Policy 9 and does not
# wait for a brand's prose scrub: a catalog entry carries an install path and a
# repository URL, so it publishes the server rather than merely naming a tool.
COMMERCIAL_SOLVERS = ("COMSOL", "FEMM", "JMAG", "CST", "TOSCA", "Opera")

# Policy 9 allowlist: surfaces that name a brand for a legitimate reason --
# a citation of published work, or a declaration that something is absent.
BRAND_CITATION_ALLOW = (
    # the cited paper's own title, and the verbatim abstract of the paper the
    # formulation comes from
    "packages/radia-mcp/src/radia_mcp/accelerator/bibliography_index_knowledge.py",
    "packages/radia-mcp/src/radia_mcp/radia_ngsolve/bibliography_index_knowledge.py",
    "packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib",
    # conference-positioning advice, where naming the incumbent tools is the
    # substance of the advice rather than a validation provenance claim
    "packages/radia-mcp/src/radia_mcp/presentation/talk_feedback.py",
    # this file, which has to spell the brand out to forbid it
    "tools/policy_lint.py",
    "packages/radia-mcp/tools/policy_lint.py",
)

# Policy 4 allowlist: genuine LAPACK / HACApK column-major interop.
CBLAS_COLMAJOR_ALLOW = (
    "rad_mmm_matrices.cpp", "rad_relaxation_methods.cpp",
    "rad_hacapk.cpp", "rad_stream_function.cpp",
)


def _sh(cmd):
    if cmd[0] == "git":
        # CI service and interactive checkout owners may differ. Trust only
        # this invocation's repository, as ci_preflight's diff helper does.
        cmd = [cmd[0], "-c", f"safe.directory={REPO}", *cmd[1:]]
    p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    return p.returncode, (p.stdout or "") + (p.stderr or "" if p.returncode not in (0, 1) else "")


class LintExecutionError(RuntimeError):
    """The lint could not run a check, which is not the same as passing it."""


def _git_grep(pattern, pathspecs, extra=()):
    # git grep exits 0 with matches and 1 without; anything else means the
    # search itself failed (not a repository, bad pathspec, git missing).
    # That used to return None, which every caller read as "no violation" --
    # a broken git made all policies pass.  A check that could not run has
    # to fail the run, so this raises and main() turns it into exit 1.
    rc, out = _sh(["git", "grep", "-n", *extra, pattern, "--", *pathspecs])
    if rc not in (0, 1):
        raise LintExecutionError(
            f"git grep exited {rc} for pattern {pattern!r}; the policy "
            f"could not be checked and is therefore not passed: {out.strip()}")
    return [ln for ln in out.splitlines() if ln.strip()]


def _git_ls_files(*patterns):
    # Same contract as _git_grep: a listing that could not be produced is
    # not an empty listing.  Three policies read "no tracked files" as
    # "no violation", so a failing git used to pass all three.
    rc, out = _sh(["git", "ls-files", *patterns])
    if rc != 0:
        raise LintExecutionError(
            f"git ls-files exited {rc} for {patterns!r}; the policy could "
            f"not be checked and is therefore not passed")
    return [ln for ln in out.splitlines() if ln.strip()]


# Phase 1 exception: exact operational route files; no directory wildcard.
# Move these to private alias configuration in phase 2. Bibliography is an
# author/citation exception, not a machine provenance exception.
PUBLIC_RUNTIME_OPERATIONAL_REASONS = {
    "tests/test_privacy_source_contract.py": "Deliberate privacy metadata and numerical-drift detector fixtures.",
    '.agents/skills/api-inventory/inventory_workflow.js': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/build/SKILL.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/cubit-run/skill.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/debug-remote-user/SKILL.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/deploy/SKILL.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/pybind-mex-bridge/SKILL.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/pyside6-health/SKILL.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/release-cubit-mesh-export/SKILL.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/release-eqnedit64/SKILL.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/release-eqnedit64/scripts/assert_release_signature.ps1': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/release-eqnedit64/scripts/release_eqnedit64.ps1': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/release-quad/SKILL.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/simulink-app-health/SKILL.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.agents/skills/verify-deploy/SKILL.md': 'Existing operational deployment/build instructions and scripts; phase 2.',
    '.github/workflows/build-test.yml': 'Existing CI runner routing and runner-local scratch commands; phase 2.',
    '.github/workflows/matlab-engine-diagnostic.yml': 'Existing CI runner routing and runner-local scratch commands; phase 2.',
    '.github/workflows/radia-fast.yml': 'Existing CI runner routing and runner-local scratch commands; phase 2.',
    '.github/workflows/release.yml': 'Existing CI runner routing and runner-local scratch commands; phase 2.',
    'AGENTS.md': 'Routing policy retained verbatim until phase 2; result publication is not exempt.',
    'CLAUDE.md': 'Routing policy retained verbatim until phase 2; result publication is not exempt.',
    'CMakeLists.txt': 'Existing operational build/runtime scratch routing; phase 2.',
    'docs/induction_heating/axisymmetric_thermal_example.py': 'Executable demonstration scratch defaults; phase 2, prose already neutralized.',
    'matlab/+radia/+kicad/buildLTspiceBlock.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+ltspice/SimCommander.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+ltspice/SimRunner.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+ltspice/rlcResonanceObjective.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+ltspice/runRLCResonanceOptimizationDemo.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+optuna/LTspiceRunner.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+optuna/SheetMetalRunner.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+optuna/Study.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/addElectromagnetOptimizationSubsystem.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/addStreamFunctionOptimizationSubsystem.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/buildElectromagnetOptimizationModel.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/buildHystereticLTspiceBlock.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/buildLTspiceBlock.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/buildLibrary.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/buildOptunaTeachingModel.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/buildStreamFunctionOptimizationModel.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/runApplication.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/runLTspiceClosedLoop.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+simulink/runRLCResonanceOptimizationDemo.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+stream/OptunaRunner.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'matlab/+radia/+topopt/optimizeHexSheetTopology.m': 'Existing executable runtime scratch defaults; phase 2, unchanged in phase 1.',
    'packages/cubit-mesh-export/src/cubit_mesh_export/cubit_gui/toolbar_probe.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/pyproject.toml': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/bem/server.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/bibliography/plans/T14_canonical.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/document_meta/tools.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/grant_writing/tools.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/matlab/optuna_boundary.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/matlab/runtime.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/meta/bug_patterns.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/meta/catalog.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/paper_writing/_digest_lints.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/paper_writing/_tex_resolver.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/paper_writing/_undefined_variables.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/paper_writing/plans/T10.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/paper_writing/plans/T14.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/paper_writing/plans/T2.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/poster/_textutil.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/poster/plans/T14_qr_audit.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/poster/plans/T15_print_readiness_audit.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/poster/plans/T6_typography_lints.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/poster/tools.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/src/radia_mcp/presentation/tools.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'packages/radia-mcp/tests/test_axifem_signature_execution.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'packages/radia-mcp/tests/test_chroma_multilingual.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'packages/radia-mcp/tests/test_gmsh_knowledge.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'packages/radia-mcp/tests/test_install_deploy_current.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'packages/radia-mcp/tests/test_matlab_radia_mex_contract.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'packages/radia-mcp/tests/test_policy_lint.py': 'Privacy detector implementation and deliberate rejection fixtures.',
    'packages/radia-mcp/tests/test_radia_ngsolve_curve_vol_rule.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'packages/radia-mcp/tools/policy_lint.py': 'Privacy detector implementation and deliberate rejection fixtures.',
    'pyproject.toml': 'Existing operational build/runtime scratch routing; phase 2.',
    'src/core/rad_hacapk_hdiv.cpp': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'src/cubit_plugin/CubitMeshExportPlugin.cpp': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'src/cubit_plugin/MeshData.cpp': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'src/cubit_plugin/compact_netgen/dllmain_debug.cpp': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'src/cubit_plugin/cubit_build.ps1': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'src/cubit_plugin/debug_init.cpp': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'src/cubit_plugin/test_ho_sphere.jou': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'src/cubit_plugin/test_o2.jou': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'src/cubit_plugin/test_o2_fallback.jou': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'src/cubit_plugin/test_o2_only.jou': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'src/radia/panels/samples/em/em_elf_quarter_xz.jou': 'Existing executable runtime path/default interop; phase 2, no result writer exemption.',
    'tests/axifem/_vol_mesh.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/axifem/test_q2_curved.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/matlab/test_hcurl_topology_optimization.m': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/matlab/test_ih_geometry_update.m': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/matlab/test_ih_python_selection.m': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/matlab/test_ltspice_hysteresis_coupling.m': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/matlab/test_ltspice_rlc_resonance_optimization.m': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/matlab/test_optuna_table.m': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/matlab/test_peec_coil_geometry_resonance.m': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/matlab/test_rlc_resonance_optimization_demo.m': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/matlab/test_setup_python_selection.m': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/matlab/test_topology_optimization.m': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/mcp_integration/test_matlab_optuna_quality.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/mcp_server/fixtures/bad_radia_script.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/mcp_server/test_radia_rules.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_beam_transfer_mex_pybind_validation.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_ci_execution_policy.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_ci_preflight_mdx.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_ci_preflight_scope.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_cubit_installers.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_docs_notebook_contract.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_editable_intent.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_filaments_from_shape.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_hysteresis.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_mdx_runner_pool.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_public_runtime_privacy.py': 'Privacy detector implementation and deliberate rejection fixtures.',
    'tests/test_radia_promotion_gate.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_release_acceptance_hosts_shared.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_release_cubit_dual.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_release_quad_editable_intent.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_release_quad_editable_source.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_release_quad_optuna_candidate.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_release_quad_phase9_contract.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_release_quad_state.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_release_temp_shadows.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_stream_function.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tests/test_validation_compute_host.py': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tools/audit_application_block_contract.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/audit_pyside6_only.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/ci_preflight.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/ci_preflight_mdx.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/download_release_asset.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/editable_intent.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/README.md': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/build/accept_release.ps1': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/build/analyze_last_operation_log.ps1': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/build/enable_isolated_font_dumps.ps1': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/build/test_background.ps1': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/build/test_external_paste.ps1': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/build/test_font_session.ps1': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/docs/CANONICAL_OPERATION.md': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/docs/PALETTE_INTENT.md': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/src/eqnedt64_app.cpp': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/eqnedit64/tests/font_lifecycle_probe.cpp': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tools/eqnedit64/tests/read_office_wire.ps1': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tools/eqnedit64/tests/test_office_bit_parity.cjs': 'Existing operational scratch/route fixtures; phase 2, no result publication.',
    'tools/install_git_hooks.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/policy_lint.py': 'Privacy detector implementation and deliberate rejection fixtures.',
    'tools/release_acceptance.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/release_cubit_dual.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/release_quad.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'tools/release_temp_shadows.py': 'Existing operational build/runtime scratch routing; phase 2.',
    'validation_test/accelerator/benchmark_beam_transfer_mex_python.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/axifem/research/prototypes/nmr_validate.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/axifem/research/prototypes/test_p2_kelvin_sphere.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/axifem/research/prototypes/test_periodic_h1henrotte.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/bem/test_sibc_hacapk_match_ngsbem.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/common_observation_tube/radia_reference_py38.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/conftest.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/cubit/build_curved_hex_bdm2_cylinder.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/cubit/run_vol_multi_geometry_validation.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/cubit/test_export_jou.jou': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/cubit/test_export_no_phantom_block.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/cubit/test_p_convergence_regression.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/cubit_mesh_export/geometric_refit_benchmark.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/cubit_mesh_export/paper_capacitor_benchmark.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/cubit_mesh_export/paper_sphere_benchmark.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/equation/compare_boxes_with_tex.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/equation/score_against_tex.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/equivalence_source/phase3_e2e_cubit_to_sol.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/esrf_three_engine/results/surface_overlap_20260913/esrf6_surface_overlap_controls_v3.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/esrf_three_engine/results/surface_overlap_20260913/esrf6_surface_overlap_precision_controls.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/feec/bench_hdiv_rt2.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/feec/bench_hdiv_rt2_cross_topology.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/feec/test_hdiv_vim_cross_element.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/feec/tune_hdiv_hmatrix_optuna.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/ffag_topopt/build_ffag_cyclic_yoke_hex.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/ffag_topopt/validation_ffag_cubit_cyclic_nonlinear_yoke.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/force_validation/validation_axisymmetric_to_3d_vol_force.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/force_validation/validation_magnetic_material_pair_vol_force.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/hcurl_response_compression/evrs_pn_convergence.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/hdiv_vim/reduced_omega_retirement/reproduce_historical.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/hysteresis/bench_hysteresis_step.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/hysteresis/play_vs_energy/compare_play_vs_energy.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/hysteresis/regenerate_binput_play_fixture.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/hysteresis/verify_bqm_table.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/hysteresis/verify_coenergy_feasibility.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/hysteresis/verify_play_to_energy.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/induction_heating/beak_sibc_validity_map.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/induction_heating/prepare_ih_axisym_mex.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/induction_heating/team36_axisymmetric/run_radia_ih.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/ltspice/validate_interval_state_handoff.m': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/ngsolve_matlab_parity/run_sparsesolv_parity.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/optimization/runLtspiceCoilGeometryResonanceOptimization.m': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/optimization/validate_ltspice_coil_geometry_resonance.m': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/optimization/validate_ltspice_rlc_resonance_ac.m': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/optimization/validate_matlab_adjoint_quality.m': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/optimization/validate_rlc_resonance_simulink.m': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/build_test_motor_mesh.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/test_fem_coilmesh_esim_golden.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/test_heat_radiation.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/test_regcoil_iron_differentiator_golden.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/test_regcoil_parity_deliverable_golden.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/test_streamfunction_golden.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/test_streamfunction_greedy_golden.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/test_streamfunction_iron_golden.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/test_streamfunction_levels_golden.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/test_streamfunction_manufacture_e2e_golden.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/panels/test_streamfunction_pin_tiling_golden.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/peec_integration/verification/experiment_fem_convergence.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/peec_integration/verification/verify_ngsolve_inductance.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/radia_mcp/generate_annular_motor_dual_lane.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/radia_mcp/generate_field_study_production_artifacts.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/radia_mcp/mesh_quality_study/run_accuracy_per_dof.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/radia_mcp/test_vol2d_postprocess.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/radia_mcp/test_vol2d_scalar.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/sparsesolv/hiruma/bench_compact_ams.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/vim_coupled/validate_curved_hex_bdm2_cubit.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
    'validation_test/vim_coupled/validate_magnetic_conductor_disk.py': 'Executable validation admission/scratch routing; phase 2, result metadata is sanitized.',
}
PUBLIC_RUNTIME_OPERATIONAL_REASONS["tools/eqnedit64/CMakeLists.txt"] = "Existing native build scratch/default routing; phase 2."
PUBLIC_RUNTIME_OPERATIONAL_ALLOW = frozenset(PUBLIC_RUNTIME_OPERATIONAL_REASONS)
PUBLIC_RUNTIME_IDENTITY = re.compile(
    r"(?i)\bintel11\b|\bmdx[12]\b|\bhibino\b(?!20\d\d)|100号機|"
    r"\b(?:192\.168\.\d{1,3}\.\d{1,3}|10\.\d{1,3}\.\d{1,3}\.\d{1,3})\b|"
    r"(?<![A-Za-z])(?:(?-i:[WS]):[\\/]|C:[\\/]+(?:temp|Users)[\\/])|"
    r"(?-i:\bLAB\b)(?![ _-](?:style|colors|colour|space))")

PUBLIC_RUNTIME_FILENAME = re.compile(r"(?i)(?:^|/)(?:intel11|mdx[12]?|hibino|lab|100)(?:[_./-]|$)")

def _identifying_runtime_filename(path):
    return bool(PUBLIC_RUNTIME_IDENTITY.search(path) or
                (path.startswith(("validation_test/", "results/")) and
                 PUBLIC_RUNTIME_FILENAME.search(path)))

PUBLIC_RUNTIME_UNC = re.compile(r"(?:^|[\s=\x22\x27])\\\\[A-Za-z][A-Za-z0-9_.-]+\\[A-Za-z0-9_.-]")

def _public_strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield key
            yield from _public_strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from _public_strings(item)

def public_runtime_identity_hits(path, content):
    """Reject identifying metadata, including UNC roots, without TeX false hits."""
    if path in PUBLIC_RUNTIME_OPERATIONAL_ALLOW:
        return []
    if path == "packages/radia-mcp/src/radia_mcp/bibliography/data/references.bib":
        # Exact published citation keys; only the author-field name is masked.
        # All other fields and entries retain host/path/IP checks.
        citation_keys = ('article', 'cefc2026_hmatrix', 'compumag2025', 'HIBINO2018128', 'HIBINO2024', 'hibino2026aca', 'Ida2014')
        def mask_author(match):
            entry = match.group()
            if match.group(1) in citation_keys:
                return re.sub(r"(?im)^(\s*author\s*=.*)$", lambda line: re.sub(r"(?i)\bHibino\b", "PublishedAuthor", line.group()), entry)
            return entry
        content = re.sub(r"(?ms)^@\w+\{([^,]+),.*?(?=^@|\Z)", mask_author, content)
    texts = [content]
    structured = False
    if path.endswith((".json", ".ipynb")):
        try:
            document = json.loads(content)
            structured = True
            texts = list(_public_strings(document))
            def contains_host(value):
                if isinstance(value, dict):
                    return any((key.lower() in ("host", "hostname", "validation_host", "host_name", "computer_name", "machine_name") and isinstance(item, str) and item.lower() not in ("windows", "linux", "darwin", "runtime")) or (key.lower() == "machine" and isinstance(item, str) and item.lower() in ("lab", "intel11", "mdx1", "mdx2", "hibino", "100")) or contains_host(item) for key, item in value.items())
                if isinstance(value, list):
                    return any(contains_host(item) for item in value)
                return False
            if contains_host(document):
                return [f"{path}: identifying host metadata"]
        except ValueError:
            pass  # Standalone test strings and invalid JSON remain checked.
    hits = []
    for number, text in enumerate(texts, 1):
        if PUBLIC_RUNTIME_IDENTITY.search(text) or PUBLIC_RUNTIME_UNC.search(text) or (structured and re.search(r"(?i)C:[\\/]+(?:temp|Users)\b", text)):
            hits.append(f"{path}:{number}: identifying runtime provenance")
    if _identifying_runtime_filename(path):
        hits.append(f"{path}: identifying filename")
    return hits

def check_public_runtime_privacy():
    candidates = set()
    for pattern in (PUBLIC_RUNTIME_IDENTITY.pattern, PUBLIC_RUNTIME_UNC.pattern):
        for record in _git_grep(pattern, [], extra=("-P", "-I")):
            candidates.add(record.split(":", 1)[0])
    candidates.update(path for path in _git_ls_files()
                      if _identifying_runtime_filename(path) or path.endswith((".json", ".ipynb")))
    hits = []
    for path in sorted(candidates):
        if path not in PUBLIC_RUNTIME_OPERATIONAL_ALLOW:
            with open(os.path.join(REPO, path), encoding="utf-8") as source:
                hits.extend(public_runtime_identity_hits(path, source.read()))
    return hits

def check_all():
    """Return list of (policy_name, ok, detail) for all 8 policies."""
    results = []

    # 1: the removed unit-switching API must not return to executable Radia.
    hits = _git_grep("FldUnits", ["src/radia/*.py", "src/radia/**/*.py"])
    results.append(("Policy 1: no FldUnits() in Radia source", not hits,
                    hits[0] if hits else "ok"))

    # 2: no tracked binaries
    bins = _git_ls_files("*.pyd", "*.dll", "*.so", "*.lib", "*.exe", "*.obj")
    results.append(("Policy 2: no tracked binaries", not bins,
                    str(bins[:3]) if bins else "ok"))

    # 3: no Helmholtz wave kernel in src/core (hodge / decomposition excepted)
    hits = _git_grep("helmholtz", ["src/core/*.cpp", "src/core/*.h",
                                   "src/core/**/*.cpp", "src/core/**/*.h"],
                     extra=["-i"]) or []
    bad = [h for h in hits if not any(
        x in h.lower() for x in ("helmholtz-hodge", "helmholtz hodge",
                                 "helmholtz-decomposition", "helmholtz decomposition"))]
    results.append(("Policy 3: no Helmholtz wave kernel in core", not bad,
                    bad[0] if bad else "ok"))

    # 4: no CblasColMajor in src/core (LAPACK/HACApK allowlisted)
    hits = _git_grep("CblasColMajor", ["src/core/*.cpp", "src/core/*.h",
                                       "src/core/**/*.cpp", "src/core/**/*.h"]) or []
    bad = [h for h in hits if not any(a in h for a in CBLAS_COLMAJOR_ALLOW)]
    results.append(("Policy 4: no CblasColMajor in core", not bad,
                    bad[0] if bad else "ok"))

    # 5: no TRACKED generated files at repo root
    root_gen = [f for f in _git_ls_files("*.msh", "*.vtu", "*.vtk", "*.vol",
                                         "*.vts") if "/" not in f]
    results.append(("Policy 5: no tracked generated files at root", not root_gen,
                    str(root_gen) if root_gen else "ok"))

    # 6: forbid the legacy import path (use src/radia).  The needle is built
    # from fragments below so this module never self-matches -- Policy 6 greps
    # for that token on a line that also references a sys.path insertion.
    _legacy = "src/" + "python"
    hits = _git_grep(_legacy, ["*.py", "**/*.py"]) or []
    bad = [h for h in hits if "sys.path" in h]
    results.append((f"Policy 6: no legacy {_legacy} import path", not bad,
                    bad[0] if bad else "ok"))

    # 7: the examples tier is retired. Public demonstrations belong in docs,
    # numerical evidence in validation_test, and regressions in tests.
    examples = _git_ls_files("examples")
    results.append(("Policy 7: no tracked retired examples tier", not examples,
                    str(examples[:3]) if examples else "ok"))

    # 8: Piola-mapped HDiv field and geometry orders are dimension-dependent.
    # Keep one explicit table instead of reintroducing a tempting but incorrect
    # global relation such as ``geometry_order <= hdiv_order + 1``.
    cap = os.path.join(REPO, "src", "radia", "vim", "_capabilities.py")
    required = {
        "_vim.py": "validate_hdiv_configuration",
        "_vim2d.py": "validate_hdiv_configuration",
        "_solve.py": "validate_hdiv_configuration",
    }
    bad = []
    try:
        with open(cap, encoding="utf-8") as f:
            cap_text = f.read()
        for token in ('HDivCapability(2, "quad", 2, (1, 2, 3), 3)',
                      'HDivCapability(3, "hex", 2, (1, 2), 2)',
                      'HDivCapability(3, "wedge", 2, (1, 2), 2)'):
            if token not in cap_text:
                bad.append(f"_capabilities.py missing {token}")
    except OSError as exc:
        bad.append(f"cannot read _capabilities.py: {exc}")

    vim_dir = os.path.join(REPO, "src", "radia", "vim")
    relation = re.compile(
        r"(?:geometry_order|curve_order)\s*(?:<=|>=|<|>)\s*"
        r"(?:self\.)?(?:order|p)(?:\s*[+-]\s*\d+)?|"
        r"(?:self\.)?(?:order|p)(?:\s*[+-]\s*\d+)?\s*"
        r"(?:<=|>=|<|>)\s*(?:geometry_order|curve_order)")
    for name in sorted(os.listdir(vim_dir)):
        if not name.endswith(".py") or name == "_capabilities.py":
            continue
        path = os.path.join(vim_dir, name)
        try:
            with open(path, encoding="utf-8") as f:
                text = f.read()
        except OSError as exc:
            bad.append(f"cannot read {name}: {exc}")
            continue
        if name in required and required[name] not in text:
            bad.append(f"{name} does not use the central HDiv capability validator")
        match = relation.search(text)
        if match:
            line = text.count("\n", 0, match.start()) + 1
            bad.append(f"{name}:{line}: ad-hoc HDiv geometry/order relation")
    results.append(("Policy 8: central HDiv geometry/order capabilities", not bad,
                    bad[0] if bad else "ok"))

    # 9: a scrubbed commercial-solver brand must not return to Radia's own
    # identifiers, labels, source, tests or stored result JSONs.  Prose that
    # CITES published work keeps its brand words, which is why this checks
    # code and result data rather than Markdown.
    hits = _git_grep("|".join(SCRUBBED_BRANDS),
                     ["src/*", "tests/*", "matlab/*", "tools/*",
                      "validation_test/*.py", "validation_test/**/*.py",
                      "validation_test/**/*.json", "packages/**/*.py"],
                     extra=["-E", "-w", "-i"]) or []
    bad = [h for h in hits
           if not any(h.startswith(a + ":") for a in BRAND_CITATION_ALLOW)]
    results.append((f"Policy 9: no scrubbed commercial brand "
                    f"({'/'.join(SCRUBBED_BRANDS)}) in Radia code or results",
                    not bad, bad[0] if bad else "ok"))

    # 10: the published radia-mcp catalog must not advertise an MCP server for
    # a commercial solver.  Naming the reference tools is one thing; shipping
    # a public directory entry with their install path and repository URL
    # publishes the servers themselves, which the boundary forbids outright.
    catalog = os.path.join(REPO, "packages", "radia-mcp", "src", "radia_mcp",
                           "meta", "catalog.py")
    bad = []
    try:
        with open(catalog, encoding="utf-8") as f:
            entries = re.findall(r'^    "([A-Za-z0-9_-]+)": \{\n(.*?)^    \},$',
                                 f.read(), re.S | re.M)
        for name, body in entries:
            named = [b for b in COMMERCIAL_SOLVERS
                     if re.search(rf"\b{b}\b", name + body, re.I)]
            if named:
                bad.append(f"catalog entry {name!r} advertises {'/'.join(named)}")
    except OSError as exc:
        bad.append(f"cannot read the radia-mcp catalog: {exc}")
    results.append(("Policy 10: the public MCP catalog advertises no "
                    "commercial-solver server", not bad, bad[0] if bad else "ok"))

    identity_hits = check_public_runtime_privacy()
    results.append(("Policy 11: public runtime privacy", not identity_hits, "\n".join(identity_hits[:20]) if identity_hits else "ok"))
    return results


def main(argv=None):
    quiet = "--quiet" in (argv if argv is not None else sys.argv[1:])
    try:
        results = check_all()
    except (LintExecutionError, OSError) as exc:
        # A lint that could not run has not passed.  Say so and fail the
        # run; silently reporting PASS here is how a broken git once made
        # every policy green.
        print(f"FAIL  policy lint could not run: {exc}", file=sys.stderr)
        return 1
    nfail = 0
    for name, ok, detail in results:
        if ok:
            if not quiet:
                print(f"PASS  {name}")
        else:
            nfail += 1
            print(f"FAIL  {name}: {detail}")
    if nfail:
        print(f"\n{nfail}/{len(results)} policies FAILED", file=sys.stderr)
        return 1
    if not quiet:
        print(f"\nall {len(results)} policies pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
