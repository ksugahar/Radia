"""
bug_patterns.py -- learned catalog of bug patterns observed in real
Radia / radia-mcp / cubit-mesh-export incidents.

Purpose
=======
Make the recurrent bug classes that bite the lab session-after-session
explicit + queryable.  Claude (and any human contributor) should query
this catalog BEFORE writing new code that touches the affected area,
so the same bug never ships twice.

How to query
============
From an MCP client:
    bug_patterns_lookup()                          # all entries
    bug_patterns_lookup(topic="panel")             # filter
    bug_patterns_lookup(topic="release")
    bug_patterns_lookup(topic="cubit-license")
    bug_patterns_lookup(severity="high")
    bug_patterns_lookup(recent_days=14)            # only patterns
                                                   # observed in the
                                                   # last 14 days

How to add a new entry
======================
When debugging a NEW bug class (not already in the catalog), append
an entry to PATTERNS below.  Fields:

  id         : short stable slug, lowercase-with-dashes
  title      : one-sentence headline
  topics     : list of keyword tags for filter
  severity   : "high" | "medium" | "low"
  first_seen : YYYY-MM-DD (first incident date)
  last_seen  : YYYY-MM-DD (latest known recurrence; update on repeat)
  what       : 1-3 sentence plain description of what the user sees
  root_cause : technical explanation (file:line if possible)
  detection  : what test / audit / CI gate catches it now
  prevention : the rule to follow (links to skill / doc)
  related    : ["memory/<file>.md", "tools/<script>.py", ...]

KEEP entries SHORT.  This file ships to PyPI in the radia-mcp wheel;
verbose narrative belongs in memory/ or docs/.
"""

from __future__ import annotations

PATTERNS: list[dict] = [
    # =====================================================
    # PANEL / SUBPROCESS BUGS
    # =====================================================
    {
        "id": "late-taskmanager-import-unboundlocal",
        "title": "Late `from ngsolve import TaskManager` causes "
                 "UnboundLocalError on earlier `with TaskManager():`.",
        "topics": ["panel", "ngsolve", "taskmanager", "scoping"],
        "severity": "high",
        "first_seen": "2026-05-30",
        "last_seen": "2026-05-30",
        "what": "Cubit Verify-Mesh dialog popped up 'NGSolve verification "
                "failed' with a Python traceback ending in "
                "UnboundLocalError on TaskManager.",
        "root_cause": "Python compiles a function body as a whole.  An "
                      "`from ngsolve import TaskManager` LATER in the "
                      "function (e.g. inside a try: block) promotes "
                      "TaskManager to a LOCAL for the whole function, "
                      "and the EARLIER `with TaskManager():` then "
                      "raises UnboundLocalError on the local that "
                      "hasn't been assigned yet.  Real site: "
                      "calc_verify_vol.py:38 vs :83 (keiko 100号機).",
        "detection": "tests/panels/test_taskmanager_scoping.py (AST sweep) "
                     "+ tools/audit_new_panel_contract.py rule C6.",
        "prevention": "Put TaskManager in the TOP-of-function import "
                      "line alongside Mesh / Integrate / CF / BND.  "
                      "Never inside a try: block AFTER first use.",
        "related": ["tests/panels/test_taskmanager_scoping.py",
                    "tools/audit_new_panel_contract.py"],
    },
    {
        "id": "ih-log-truncated-by-super-then-append",
        "title": "Panel `.log` ends before the IH-specific summary "
                 "because the override called super first.",
        "topics": ["panel", "persistence-policy", "ih", "log"],
        "severity": "medium",
        "first_seen": "2026-05-30",
        "last_seen": "2026-05-30",
        "what": "User attaches the .log to a bug report; the .log "
                "ends at the generic _append_standard_summary block "
                "and is missing ~150 lines of IH-specific output (L_coil, "
                "P_workpiece, T_max, file paths, auto-open banner).",
        "root_cause": "IHWindow._on_finished called super()._on_finished "
                      "FIRST.  The base class's _persist_output_log() "
                      "ran at that point.  Then the override appended "
                      "150 more lines to _output that never made it "
                      "to the .log.",
        "detection": "Code review of any AnalysisWindow._on_finished "
                     "override.  No static gate today; consider adding.",
        "prevention": "Any subclass that overrides _on_finished and "
                      "appends additional output AFTER super() must "
                      "ALSO call self._persist_output_log() at the END "
                      "of the override.  See radia_ih.py:1784 for the "
                      "fix template.",
        "related": ["src/radia/radia_ih.py:1784",
                    "src/radia/radia_gui_base.py:1252"],
    },
    {
        "id": "calc-result-key-misnamed-time",
        "title": "calc_*.py emits `t_solve` instead of `t_solve_s`; "
                 "panel summary silently skips the entry.",
        "topics": ["panel", "result-output-policy", "naming"],
        "severity": "medium",
        "first_seen": "2026-05-30",
        "last_seen": "2026-05-30",
        "what": "Panel Output window shows ne / DoF but the 'Compute "
                "time' block is empty even though the calc ran for "
                "minutes.",
        "root_cause": "_append_standard_summary matches every result-"
                      "dict key against the pattern `t_*_s`.  A key "
                      "without `_s` suffix (e.g. `t_solve`, `t_mesh`) "
                      "or per-iteration only (not aggregated to top "
                      "level) is silently dropped.",
        "detection": "Result Output Policy (CLAUDE.md 2026-05-29); "
                     "no static gate yet -- relies on developer "
                     "discipline + Output-window eye check.",
        "prevention": "Always use `t_<step>_s` suffix.  For multi-"
                      "iteration solvers, aggregate per-iter `t_solve` "
                      "into a top-level `t_solve_s` via "
                      "`sum(h['t_solve'] for h in history)`.",
        "related": ["CLAUDE.md: Result Output Policy 2026-05-29"],
    },
    {
        "id": "calc-non-time-prefixed-t",
        "title": "`t_end_s` / `t_ext_C` get classified as compute "
                 "time by the panel summary helper.",
        "topics": ["panel", "result-output-policy", "naming"],
        "severity": "low",
        "first_seen": "2026-05-30",
        "last_seen": "2026-05-30",
        "what": "Panel shows '`t_ext_C 20.00`' in the Compute Time "
                "block.  20°C is the ambient temperature, not a solve "
                "time of 20 seconds.",
        "root_cause": "`_append_standard_summary` matches any key "
                      "shaped `t_*_s`.  `t_ext_C` and `t_end_s` "
                      "(sim-end time) accidentally fit the pattern.",
        "detection": "Visual review of panel summary table; no "
                     "automatic gate.",
        "prevention": "Don't prefix non-time fields with `t_`.  Use "
                      "`T_ext_C` (capital T = Temperature) and "
                      "`sim_end_s` instead.",
        "related": [],
    },
    {
        "id": "gmsh-file-key-mismatch",
        "title": "Open-GMSH button stays disabled because the calc's "
                 "result dict uses a key the panel doesn't recognize.",
        "topics": ["panel", "gmsh", "result-output-policy"],
        "severity": "medium",
        "first_seen": "2026-04-01",
        "last_seen": "2026-05-30",
        "what": "User finished a Run, .msh exists on disk, but the "
                "Open GMSH button greys out.",
        "root_cause": "AnalysisWindow._on_finished scans the result "
                      "dict for one of FOUR keys: gmsh_file, "
                      "field_gmsh_file, msh_output, msh_file.  Other "
                      "names (e.g. `output_msh`, `mesh_path`) are "
                      "silently ignored.",
        "detection": "tests/panels/test_open_gmsh_button.py (manual "
                     "result-dict tests).",
        "prevention": "Always use one of the 4 recognized keys.  "
                      "Document the convention in calc_*.py docstring.",
        "related": ["tests/panels/test_open_gmsh_button.py",
                    "src/radia/radia_gui_base.py: _on_finished"],
    },

    # =====================================================
    # CUBIT PLUGIN / .CCM BUGS
    # =====================================================
    {
        "id": "phantom-entity-from-meshexport-interface",
        "title": "iface->id_from_handle() materialises a phantom "
                 "block/sideset/nodeset in Cubit's database.",
        "topics": ["cubit", "ccm", "mesh-export", "side-effect"],
        "severity": "high",
        "first_seen": "2026-05-30",
        "last_seen": "2026-05-31",
        "what": "After `export netgen`, Cubit shows a new "
                "phantom block (id K+1) that the user never created.  "
                "Subsequent parse_cubit_list('volume', 'in block K+1') "
                "errors with 'No block with ID K+1 was found'.  "
                "Reported by keiko 100号機 2026-05-30 on a 6-turn "
                "loft coil journal.",
        "root_cause": "MeshExportInterface::get_block_list (and the "
                      "sideset / nodeset variants) returns a 'default' "
                      "handle covering elements / faces / nodes not "
                      "assigned to any user-defined group.  Calling "
                      "iface->id_from_handle() on that default handle "
                      "MATERIALISES the group as a real database entry "
                      "with id ~ K+1.  Source: MeshData.cpp::"
                      "extract_elements / _sidesets / _nodesets.",
        "detection": "tests/cubit/test_export_no_phantom_block.py "
                     "(state-snapshot diff before/after export). "
                     "cubit-smoke-test also exercises the round-trip.",
        "prevention": "BEFORE calling iface->get_*_list, snapshot the "
                      "USER-defined ids via "
                      "CubitInterface::parse_cubit_list(\"block\"/"
                      "\"sideset\"/\"nodeset\", \"all\").  Skip any "
                      "handle whose id_from_handle is NOT in the "
                      "snapshot.  See dce6ee3b for the sideset+nodeset "
                      "fix mirroring da1bfc72 for blocks.",
        "related": ["src/cubit_plugin/MeshData.cpp",
                    "tests/cubit/test_export_no_phantom_block.py"],
    },

    # =====================================================
    # CUBIT LICENSE BUGS
    # =====================================================
    {
        "id": "cubit-2025-8-logout-only-clears-local-cache",
        "title": "rlm_activate --logout no longer releases the server-"
                 "side seat in Cubit 2025.8+ (Web portal Deactivate "
                 "required).",
        "topics": ["cubit", "license", "rlm", "2025.8"],
        "severity": "high",
        "first_seen": "2026-05-30",
        "last_seen": "2026-05-30",
        "what": "After logging off LAB + 100号機 + mdx, all three "
                "still get NoAvailableSeats on the next --login.  Seat "
                "stays stuck for hours.",
        "root_cause": "Cubit 2025.8 introduced a new licensing portal.  "
                      "rlm_activate --logout exit-0 'Logged out' is "
                      "now a LOCAL cache clear ONLY -- the server "
                      "still considers the (user, machine) pair as "
                      "holding the activation.",
        "detection": "rlm_activate --login returns 'NoAvailableSeats' "
                     "from every machine in the pool, despite all "
                     "having run --logout cleanly.",
        "prevention": "Use the WEB PORTAL Deactivate at "
                      "https://coreform.com/account/ -> license -> "
                      "[...] -> Deactivate.  Alternative CLI trick: "
                      "log in from any machine with a DIFFERENT account "
                      "(displaces the previous activation row).  See "
                      "S:\\coreformcubit\\CoreformCubit.ps1 "
                      "-ReleaseSeat switch (2026-05-30 update).",
        "related": ["memory/feedback_coreform_http_500_login.md",
                    "memory/reference_rlm_activate_logoff_procedure.md"],
    },
    {
        "id": "cubit-rlm-per-user-per-machine-activation",
        "title": "Coreform RLM activation is keyed on (Windows user, "
                 "machine).  Administrator's --logout cannot release "
                 "keiko's seat even on the same machine.",
        "topics": ["cubit", "license", "rlm", "per-user"],
        "severity": "high",
        "first_seen": "2026-05-30",
        "last_seen": "2026-05-30",
        "what": "Admin ran rlm_activate --logout via SSH on 100号機.  "
                "Reported 'Logged out exit 0'.  Other users still "
                "see NoAvailableSeats.",
        "root_cause": "Server-side activations are keyed on (user, "
                      "machine, account) triple.  Administrator's "
                      "--logout only releases the (administrator, "
                      "INTEL11, 144576) row.  keiko's (keiko, "
                      "INTEL11, 144576) row is untouched.",
        "detection": "Inspect Cubit-Log.txt under EACH user profile, "
                     "not just one.",
        "prevention": "To release a specific user's seat: (a) log in "
                      "to coreform.com/account/ and Deactivate the "
                      "row for that (user, machine), OR (b) have that "
                      "user themselves run rlm_activate --logout in "
                      "their own Windows session, OR (c) use "
                      "runas/Start-Process -Credential to run as that "
                      "user (needs their Windows password).  Per-user "
                      "cache lives at %LOCALAPPDATA%\\Coreform\\.",
        "related": ["memory/reference_rlm_activate_logoff_procedure.md",
                    ".claude/skills/cubit-license/SKILL.md"],
    },

    # =====================================================
    # RELEASE / CI BUGS
    # =====================================================
    {
        "id": "tag-concurrency-race-3-tags-same-sha",
        "title": "Pushing 3 tags at the same SHA cancels the middle "
                 "queued CI run (concurrency-group=ci-${sha}).",
        "topics": ["release", "ci", "github-actions", "concurrency"],
        "severity": "medium",
        "first_seen": "2026-05-30",
        "last_seen": "2026-06-01",
        "what": "v4.85.0 tag CI was cancelled within 1 second of "
                "starting.  Release-publish workflow then never "
                "fired (workflow_run: success required).",
        "root_cause": "build-test.yml has "
                      "`concurrency: ci-${{ github.sha }}` with "
                      "cancel-in-progress=false (push events).  "
                      "Three simultaneous pushes on the same SHA: "
                      "1st runs, 2nd queues, 3rd CANCELS the queued "
                      "2nd and replaces it.  Middle tag loses.",
        "detection": "github.com/<org>/<repo>/actions shows the tag "
                     "CI as 'cancelled' (not 'failure').  ci-verify "
                     "would catch it via the missing junit XMLs.",
        "prevention": "Push tags ONE AT A TIME with >= 30 seconds "
                      "between pushes.  Or push max 2 at once.  Long-"
                      "term fix: change concurrency group to "
                      "`ci-${sha}-${ref}` so each tag has its own "
                      "group.",
        "related": ["memory/feedback_tag_push_concurrency_race.md",
                    ".github/workflows/build-test.yml"],
    },
    {
        "id": "init-py-version-mismatch-vs-pyproject",
        "title": "__init__.py version not bumped to match pyproject.toml; "
                 "release wheel verify rejects the package.",
        "topics": ["release", "version", "wheel"],
        "severity": "high",
        "first_seen": "2026-05-30",
        "last_seen": "2026-05-30",
        "what": "cme-v0.10.10 release workflow's 'Verify wheel' step "
                "failed with 'Version mismatch -- wheel=0.10.10 vs "
                "__init__.__version__=0.10.9'.  PyPI never got the "
                "new version.",
        "root_cause": "Bumped pyproject.toml but forgot to bump "
                      "src/<pkg>/__init__.py's __version__ in the "
                      "same commit.  Both are read by different tools.",
        "detection": "release.yml / release-cubit-mesh-export.yml / "
                     "release-radia-mcp.yml all check version "
                     "consistency in the 'Verify wheel' step.  "
                     "Catches at PUBLISH time, not at commit time.",
        "prevention": "Always bump BOTH files in lockstep.  release-"
                      "triple skill's Phase 2 lists them explicitly.  "
                      "Could add a pre-commit hook or a "
                      "test_version_consistency.py to catch earlier.",
        "related": ["tools/release_triple.py: Phase 2"],
    },
    {
        "id": "lab-editable-drift-after-pip-force-reinstall",
        "title": "`pip install --force-reinstall <package>` from PyPI "
                 "clobbers LAB's editable install pointer.",
        "topics": ["release", "lab", "editable"],
        "severity": "medium",
        "first_seen": "2026-04-28",
        "last_seen": "2026-06-01",
        "what": "After release-triple Phase 8, LAB's "
                "`pip show radia` no longer says "
                "'Editable project location: S:\\Radia\\01_GitHub'.  "
                "Source edits no longer flow to runtime; dev loop is "
                "broken.",
        "root_cause": "pip install --force-reinstall replaces the "
                      ".pth pointer with a regular install.  Easy "
                      "trap when Phase 8 deploy commands accidentally "
                      "run on LAB.",
        "detection": "tools/release_triple.py done's LAB-editable "
                     "gate (POLICY 2026-05-27).",
        "prevention": "Phase 8 deploy commands are for 100号機 + mdx "
                      "ONLY -- never run them locally on LAB.  After "
                      "any release, run `release_triple done` and "
                      "fix any DRIFT it reports.",
        "related": ["memory/project_ci_radia_mcp_editable_drift_fix.md",
                    "tools/release_triple.py: cmd_done"],
    },
    {
        "id": "ninja-stale-obj-after-netgen-upgrade",
        "title": "Ninja `#deps 0` masks include-order changes; .pyd "
                 "ships with mixed pybind11 ABI.",
        "topics": ["build", "ci", "ninja", "pybind11"],
        "severity": "high",
        "first_seen": "2026-05-30",
        "last_seen": "2026-06-01",
        "what": "3 intermittent main-push CI failures with 0 artifacts "
                "uploaded.  pytest crashed before junit XMLs were "
                "written because import radia segfaulted on ABI "
                "mismatch.",
        "root_cause": "After netgen pip upgrade, ngsolve started "
                      "shipping a bundled pybind11 with a different "
                      "ABI from the system pybind11.  CMakeLists.txt "
                      "was updated to `BEFORE PRIVATE NGSOLVE_BUNDLED_"
                      "INC` to use the bundled headers, but ninja's "
                      ".ninja_deps cache (#deps 0) didn't see the "
                      "include-order change, so old .obj files were "
                      "kept.",
        "detection": "L0+L1+L2 (pytest -r aR + Tee stdout + junit "
                     "annotation + LAB junit archive) makes the "
                     "crash visible without artifact dive.",
        "prevention": "Force `-Rebuild` always-on for main / PR / tag "
                      "(commits e6518268 + 7f46d995).  Adds ~15-25 min "
                      "to CI but eliminates the entire class.",
        "related": ["memory/feedback_ninja_stale_obj_after_netgen_upgrade.md",
                    ".github/workflows/build-test.yml"],
    },
    {
        "id": "ci-workflow-stale-binary-path-after-relocation",
        "title": "After relocating/renaming a bundled binary, the CI "
                 "workflow still sourced it from the OLD path/name -> "
                 "first push RED at the binaries-fetch step.",
        "topics": ["release", "ci", "tier2", "binary", "github-actions"],
        "severity": "high",
        "first_seen": "2026-06-01",
        "last_seen": "2026-06-01",
        "what": "First push after Tier-2 (Cubit plugin .pyd/.ccm moved "
                "OUT of src/radia into the cubit-mesh-export package) + "
                "the de-radia rename turned build-test CI RED at the "
                "'fetch cubit_mesh_curver.pyd from binaries' step (6 "
                "download attempts, exit 1).  Job died before tests -> "
                "no junit XML.  Build + tests themselves were fine.",
        "root_cause": "build-test.yml fetched the .pyd to src/radia and "
                      "the wheel-build step copied .pyd/.ccm FROM "
                      "src/radia -- but Build.ps1 now propagates them "
                      "into packages/cubit-mesh-export/src/cubit_mesh_"
                      "export.  The renamed asset (radia_cubit_mesh.pyd "
                      "-> cubit_mesh_curver.pyd) was also not yet in the "
                      "binaries release.",
        "detection": "ci-verify RED; runner _diag Worker log 'Failed to "
                     "download cubit_mesh_curver.pyd after 6 attempts'.",
        "prevention": "When relocating OR renaming a bundled binary, "
                      "grep ALL .github/workflows/*.yml (build-test + "
                      "release-*) for the old path/name and fix: fetch "
                      "dest, wheel-build source, pre-push hook upload "
                      "list, binaries-release asset name.  The runner has "
                      "the Cubit SDK so Build.ps1 supplies the .pyd "
                      "locally -- make the binaries fetch a FALLBACK "
                      "(skip-if-present), not mandatory.",
        "related": [".github/workflows/build-test.yml",
                    "memory/project_tier2_cme_sole_plugin_shipper_2026_06_01.md"],
    },
    {
        "id": "release-smoke-admin-license-fail-nonblocking",
        "title": "release_triple phase8 reports FAIL because cubit-smoke "
                 "over SSH runs as Administrator (no license) -- but the "
                 "deploy itself is verified OK.",
        "topics": ["release", "ci", "cubit", "license", "smoke"],
        "severity": "medium",
        "first_seen": "2026-06-01",
        "last_seen": "2026-06-01",
        "what": "phase8 on 100号機 + mdx fails: cubit-smoke 'export did "
                "not produce smoke.vol', Cubit exit=1 (0xC0000005); "
                "cubit.log says 'License Error: No license found'.  This "
                "stops `all` before phase8e/phase9 even though the "
                "binaries deployed fine.",
        "root_cause": "cubit-smoke-test over SSH runs as the Administrator "
                      "Windows profile, whose Coreform license is not "
                      "activated (renewals cache absent), so Cubit can't "
                      "start.  The DEPLOY is fine: cubit-plugin-install "
                      "--verify-only passes (sha256 match + compat OK + "
                      "old radia_cubit.* removed).",
        "detection": "cubit.log 'License Error: No license found' under "
                     "the smoke temp dir; --verify-only is green.",
        "prevention": "NON-BLOCKING: real lab users have their own "
                      "licenses.  `release_triple done` (preflight + "
                      "verify-editable + phase9) has NO smoke, so it "
                      "passes -- use it as the release gate, not phase8's "
                      "smoke.  Don't burn a Learn seat activating the "
                      "Administrator profile just for the smoke.",
        "related": [".claude/skills/cubit-license/SKILL.md",
                    "memory/reference_rlm_activate_logoff_procedure.md"],
    },
    {
        "id": "gitignored-skills-lag-code-renames",
        "title": ".claude/skills (gitignored, LAB-local) keep OLD "
                 "binary/command/cmake-target names after a code rename.",
        "topics": ["skills", "rename", "deploy", "release"],
        "severity": "medium",
        "first_seen": "2026-06-01",
        "last_seen": "2026-06-01",
        "what": "After renaming the Cubit plugin "
                "(radia_cubit.*->cubit_mesh_export.*/cubit_mesh_curver) + "
                "the command verb (radia_export->export, nastran-> "
                "jmag_nastran), the deploy/release-triple/build skills "
                "still listed the OLD cmake targets (radia_cubit_ccm), "
                "binary names, and radia_export commands -- they would "
                "fail if run verbatim.",
        "root_cause": ".claude/skills/ is gitignored (LAB-local dev "
                      "convenience), so a rename in TRACKED code does NOT "
                      "propagate to the skills; they must be swept by "
                      "hand.",
        "detection": "grep -r '<old-token>' .claude/skills after any "
                     "binary / command / target rename.",
        "prevention": "On ANY binary/command/cmake-target rename, grep "
                      ".claude/skills for the old token and sweep "
                      "(byte-level rename).  7 skills needed it this "
                      "time: deploy, release-triple, build, radia-plugin-"
                      "check, cubit-license, cubit-run, pyside6-health.",
        "related": ["memory/project_tier2_cme_sole_plugin_shipper_2026_06_01.md"],
    },

    # =====================================================
    # SILENT FALLBACK / NO-FALLBACKS POLICY VIOLATIONS
    # =====================================================
    {
        "id": "silent-fallback-on-unknown-enum",
        "title": "Setter accepts unknown enum value, silently falls "
                 "back to default (HACApK SetClusterStrategy).",
        "topics": ["api", "no-fallbacks", "validation"],
        "severity": "low",
        "first_seen": "2026-05-30",
        "last_seen": "2026-05-30",
        "what": "rad.SetClusterStrategy(2) accepted, getter returns "
                "2, but solver uses BBOX because the C-side dispatch "
                "is `== CHACAPK_CLUSTER_PCA` (==1) and unknown values "
                "fall through to BBOX.  User has no idea their setting "
                "had no effect.",
        "root_cause": "C setter `g_cluster_strategy = strategy` does "
                      "no validation.  Dispatch is an == comparison "
                      "that silently treats any other value as the "
                      "default branch.",
        "detection": "Code review.  No automated gate.",
        "prevention": "Validate at the setter (pybind layer or C "
                      "side): raise ValueError if outside the allowed "
                      "set.  Per CLAUDE.md 'No Fallbacks -- Fail Fast, "
                      "Fail Loud' policy.",
        "related": ["CLAUDE.md: No Fallbacks -- Fail Fast, Fail Loud"],
    },

    # =====================================================
    # TEST INFRASTRUCTURE BUGS
    # =====================================================
    {
        "id": "test-file-handle-leak-pytest-unraisable",
        "title": "`body = open(path).read()` leaks file handle; pytest "
                 "fails on PytestUnraisableExceptionWarning.",
        "topics": ["test", "python", "file-handle"],
        "severity": "low",
        "first_seen": "2026-05-30",
        "last_seen": "2026-05-30",
        "what": "CI panel-test job reports 3 failures: "
                "PytestUnraisableExceptionWarning 'Exception ignored "
                "in <_io.FileIO ...>'.",
        "root_cause": "open(...).read() never closes the file; GC "
                      "warns later.  pytest catches the unraisable "
                      "warning as a failure.",
        "detection": "CI panel-test runs pytest with default settings "
                     "which surface unraisable warnings.",
        "prevention": "Always use `with open(...) as f: body = f.read()` "
                      "in tests (and elsewhere).  Never the one-liner "
                      "form that leaks the handle.",
        "related": [],
    },
    {
        "id": "pardiso-mkl-thread-dll-fails-in-pytest-subprocess",
        "title": "A calc_*.py run with --solver pardiso FAILS in a "
                 "subprocess spawned under pytest on LAB "
                 "(mkl_intel_thread.dll cannot load).",
        "topics": ["test", "pytest", "pardiso", "mkl", "ngsolve", "lab"],
        "severity": "medium",
        "first_seen": "2026-06-02",
        "last_seen": "2026-06-02",
        "what": "A golden test that subprocess-runs a calc script with "
                "--solver pardiso fails (returncode 2) with 'Intel MKL "
                "FATAL ERROR: Cannot load mkl_intel_thread.dll'.  The SAME "
                "command run directly (no pytest) works fine.",
        "root_cause": "pytest's conftest MKL add_dll_directory shadow "
                      "poisons the DLL search path the subprocess "
                      "inherits; MKL's threading layer fails to load.  "
                      "sparsecholesky-based calcs are unaffected (no MKL "
                      "threads).  Same family as the pytest+PySide6 DLL "
                      "crash.",
        "detection": "Golden test red on LAB only; CI ignores tests/panels "
                     "so it does not run there either.",
        "prevention": "Guard pardiso golden tests: if returncode != 0 and "
                      "'MKL FATAL ERROR'/'Cannot load mkl' in output -> "
                      "pytest.skip; verify via a DIRECT (non-pytest) run.  "
                      "See tests/panels/test_fem_coilmesh_esim_golden.py.",
        "related": ["tests/panels/test_fem_coilmesh_esim_golden.py"],
    },

    # =====================================================
    # CI / GIT INFRASTRUCTURE BUGS
    # =====================================================
    {
        "id": "stale-index-lock-in-shared-clone",
        "title": "An interrupted background `git commit` in the shared "
                 "NAS clone leaves a stale .git/index.lock; the next "
                 "commit fails with 'a git process may have crashed'.",
        "topics": ["git", "release", "worktree", "lab"],
        "severity": "medium",
        "first_seen": "2026-06-02",
        "last_seen": "2026-06-02",
        "what": "git commit fails: 'a git process may have crashed in "
                "this repository earlier: remove the file manually to "
                "continue'.  HEAD does not advance; the staged index "
                "survives.",
        "root_cause": "The LAB main clone is shared with ~7 "
                      ".claude/worktrees/* sessions.  A backgrounded / "
                      "interrupted git commit can leave a 0-byte "
                      ".git/index.lock orphaned.  Worktree sessions use "
                      "their OWN index, so the main-clone index.lock is "
                      "safe to clear when no git process is running.",
        "detection": "Commit output shows the lock message; HEAD unchanged "
                     "and `git status` still shows the files staged (A).",
        "prevention": "Confirm `Get-Process git` is empty, `rm "
                      ".git/index.lock`, retry (staged index survives).  "
                      "Prefer FOREGROUND commits here; commit to a branch "
                      "you are not on via a dedicated worktree (its own "
                      "index avoids the collision entirely).",
        "related": [],
    },
    {
        "id": "policy-lint-helmholtz-hodge-false-positive",
        "title": "Policy Lint 'No Helmholtz in C++ core' false-positives "
                 "on the Helmholtz-Hodge DECOMPOSITION.",
        "topics": ["ci", "policy-lint", "greens-function", "hodge"],
        "severity": "medium",
        "first_seen": "2026-06-02",
        "last_seen": "2026-06-02",
        "what": "Policy Lint workflow red: 'Policy 3 failed: Laplace "
                "kernel only in C++ core', listing rad_application.h, "
                "rad_hacapk.{cpp,h}, rad_relaxation_methods.cpp.",
        "root_cause": "Policy 3 was a blunt `grep -i helmholtz` over "
                      "src/core.  The loop-projection-hodge merge "
                      "(v4.89.0) added 'Helmholtz-Hodge decomposition' "
                      "comments (SetLoopProjection ker(N) loop removal) "
                      "-- a vector-calculus DECOMPOSITION, NOT the "
                      "forbidden Helmholtz WAVE kernel (e^{-jkr}/r).  The "
                      "core is still Laplace-only.",
        "detection": ".github/workflows/policy-lint.yml Policy 3 step "
                     "(runs on push to main / PRs to main).",
        "prevention": "Policy 3 grep now excepts "
                      "'helmholtz[- ]?hodge|decomposition' while still "
                      "catching a real exp(-j*k*r) kernel (fixed on main "
                      "48d5205f).  When adding a forbidden-term grep, "
                      "scope it to the actual violation, not a substring "
                      "a legitimate concept shares.",
        "related": [".github/workflows/policy-lint.yml",
                    "CLAUDE.md: Green's Function: Laplace Kernel Only"],
    },
]


def lookup(topic: str | None = None,
           severity: str | None = None,
           recent_days: int | None = None) -> list[dict]:
    """Filter PATTERNS by topic / severity / age.

    Args:
        topic: substring match against entry topics tags + id + title.
            Case-insensitive.  None = no topic filter.
        severity: "high" | "medium" | "low".  None = any.
        recent_days: only entries with last_seen within N days of today.
            None = no age filter.

    Returns:
        List of matching entry dicts, ordered by last_seen desc then
        severity (high first).
    """
    from datetime import date, timedelta

    def _matches(p):
        if topic is not None:
            t = topic.lower()
            hay = (" ".join(p["topics"]) + " " + p["id"]
                   + " " + p["title"]).lower()
            if t not in hay:
                return False
        if severity is not None and p["severity"] != severity:
            return False
        if recent_days is not None:
            cutoff = date.today() - timedelta(days=recent_days)
            try:
                seen = date.fromisoformat(p["last_seen"])
            except ValueError:
                return False
            if seen < cutoff:
                return False
        return True

    sev_rank = {"high": 0, "medium": 1, "low": 2}
    out = [p for p in PATTERNS if _matches(p)]
    out.sort(key=lambda p: (p["last_seen"], -sev_rank.get(p["severity"], 99)),
             reverse=True)
    return out


def stats() -> dict:
    """Counts by topic + severity, for a quick health check."""
    from collections import Counter
    topic_count = Counter()
    sev_count = Counter()
    for p in PATTERNS:
        sev_count[p["severity"]] += 1
        for t in p["topics"]:
            topic_count[t] += 1
    return {
        "total": len(PATTERNS),
        "by_severity": dict(sev_count),
        "by_topic": dict(topic_count.most_common()),
    }
