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
        "detection": "Code review of any legacy AnalysisWindow._on_finished "
                     "override; current application blocks are covered by "
                     "tests/test_simulink_application.py and "
                     "tests/test_application_interface_manifest.py.",
        "prevention": "Prefer the Simulink application-runner artifact "
                      "contract; do not restore notebook-workbench adapters.",
        "related": ["tests/test_simulink_application.py",
                    "tests/test_application_interface_manifest.py"],
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
                      "public-safe curated corpus "
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
                      "qud Phase 2 lists them explicitly.  "
                      "Could add a pre-commit hook or a "
                      "test_version_consistency.py to catch earlier.",
        "related": ["tools/release_qud.py: Phase 2"],
    },
    {
        "id": "lab-editable-drift-after-pip-force-reinstall",
        "title": "`pip install --force-reinstall <package>` from PyPI "
                 "clobbers LAB's editable install pointer.",
        "topics": ["release", "lab", "editable"],
        "severity": "medium",
        "first_seen": "2026-04-28",
        "last_seen": "2026-06-01",
        "what": "After release-qud Phase 8, LAB's "
                "`pip show radia` no longer says "
                "'Editable project location: public-safe curated corpus'.  "
                "Source edits no longer flow to runtime; dev loop is "
                "broken.",
        "root_cause": "pip install --force-reinstall replaces the "
                      ".pth pointer with a regular install.  Easy "
                      "trap when Phase 8 deploy commands accidentally "
                      "run on LAB.",
        "detection": "tools/release_qud.py done's LAB-editable "
                     "gate (POLICY 2026-05-27).",
        "prevention": "Phase 8 deploy commands keep LAB + 100号機 editable, "
                      "deploy hibino from PyPI, and leave mdx to phase8e.  "
                      "After any release, run `python tools/release_qud.py done` and "
                      "fix any DRIFT it reports.",
        "related": ["memory/project_ci_radia_mcp_editable_drift_fix.md",
                    "tools/release_qud.py: cmd_done"],
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
                      "export.  The renamed curver asset was also not yet "
                      "in the binaries release.",
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
        "title": "release_qud phase8 reports FAIL because cubit-smoke "
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
                      "licenses.  `python tools/release_qud.py done` (preflight + "
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
                "jmag_nastran), the deploy/release-qud/build skills "
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
                      "time: deploy, release-qud, build, radia-plugin-"
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
    {
        "id": "mdx-pip-orphans-block-reinstall",
        "title": "pip interrupted-uninstall orphans (~adia, ~gsolve, ...) "
                 "in site-packages leave a package NOT importable after a "
                 "force-reinstall.",
        "topics": ["deploy", "pip", "mdx", "release", "lab"],
        "severity": "medium",
        "first_seen": "2026-06-02",
        "last_seen": "2026-06-02",
        "what": "release_qud phase8e reported radia not importable on "
                "mdx (ModuleNotFoundError) even though pip 'Successfully "
                "installed' earlier; cubit-smoke 'Cannot locate "
                "ih_bem_sample.jou -- install the radia package'.",
        "root_cause": "A previously interrupted pip uninstall left "
                      "tilde-prefixed orphan dirs (~adia, ~adia-X.dist-info, "
                      "~-mpy, ~pds, ~ryptography, ~ydantic_core) in "
                      "site-packages; pip then skips / half-writes the real "
                      "package and import fails.",
        "detection": "import <pkg> raises on the target; "
                     "Get-ChildItem site-packages -Filter '~*' lists them.",
        "prevention": "rm the site-packages '~*' orphan dirs, THEN "
                      "pip install --force-reinstall; the package then "
                      "imports cleanly.",
        "related": ["memory/project_radia_ih_june_completion_2026_06_02.md"],
    },
    {
        "id": "radia-ih-exe-launcher-lock-on-force-reinstall",
        "title": "pip --force-reinstall radia fails WinError 32 on "
                 "Scripts/radia-ih.exe when a radia-ih panel is running.",
        "topics": ["deploy", "pip", "windows", "release", "lab"],
        "severity": "low",
        "first_seen": "2026-06-02",
        "last_seen": "2026-06-02",
        "what": "pip ERROR: [WinError 32] ... 'radia-ih.exe' -> "
                "'radia-ih.exe.deleteme' (file in use).  The PACKAGE still "
                "installs + imports; only the entry-point launcher script "
                "rewrite is blocked.",
        "root_cause": "A running radia-ih panel process holds an open "
                      "handle on Scripts/radia-ih.exe; pip cannot replace "
                      "the launcher exe.",
        "detection": "pip exits non-zero with WinError 32 on a "
                     "Scripts/radia-*.exe at the end of the install.",
        "prevention": "Kill radia-* (+ coreform_cubit / mcp-server*) "
                      "processes, rm Scripts/radia-*.exe(.deleteme), then "
                      "reinstall so the launcher writes cleanly.",
        "related": [],
    },
    {
        "id": "tools-md-drift-wip-contamination",
        "title": "docs/TOOLS.md committed out-of-sync with radia_mcp code -- "
                 "regenerated WITH uncommitted WIP tools, or NOT regenerated "
                 "after adding tools -- radia-mcp matrix drift gate goes red.",
        "topics": ["ci", "release", "tools-md", "radia-mcp"],
        "severity": "high",
        "first_seen": "2026-06-04",
        "last_seen": "2026-06-05",
        "what": "The radia-mcp matrix 'TOOLS.md drift gate' (and the "
                "self-hosted CI 'Run basic tests' -> tests/mcp_server/"
                "test_tools_doc.py) fail: committed docs/TOOLS.md != "
                "gen_tools_doc.py regenerated from the COMMITTED code.  Two "
                "directions seen 2026-06-04/05: (a) a session regenerated "
                "TOOLS.md while ANOTHER session's force_validation tool was "
                "uncommitted WIP, baking a 339/force_validation row into the "
                "commit while committed code had only 338 (06003d40, 24f13628); "
                "(b) the SF/FEMM sessions ADDED femm_parity_documentation + "
                "force_validation + SF topics in committed code but never "
                "regenerated TOOLS.md, leaving it at the 338 snapshot "
                "(5c4b0216).",
        "root_cause": "gen_tools_doc.py imports the LIVE radia_mcp packages, so "
                      "regenerating on a dirty working tree captures UNCOMMITTED "
                      "tools; conversely, adding a tool without regenerating "
                      "leaves TOOLS.md behind.  CI regenerates from the clean "
                      "checked-out commit, so any mismatch with the committed "
                      "TOOLS.md fails the gate.",
        "detection": "python tools/ci_preflight.py (TOOLS.md drift gate is "
                     "WIP-aware: it WARNS when radia_mcp/src has uncommitted "
                     ".py changes); radia-mcp-matrix.yml 'TOOLS.md drift gate'; "
                     "tests/mcp_server/test_tools_doc.py.",
        "prevention": "Regenerate from CLEAN committed code: commit the "
                      "tool-adding code FIRST then regenerate, OR regenerate in "
                      "a detached worktree at HEAD (no WIP).  ALWAYS commit the "
                      "tool code + regenerated docs/TOOLS.md TOGETHER in one "
                      "commit.  Run `python tools/ci_preflight.py` before push.",
        "related": ["tools/ci_preflight.py",
                    "packages/radia-mcp/scripts/gen_tools_doc.py",
                    ".github/workflows/radia-mcp-matrix.yml",
                    "tests/mcp_server/test_tools_doc.py"],
    },
    {
        "id": "heavy-import-collection-break-minimal-dep-matrix",
        "title": "A new radia-mcp test imports ngsolve/netgen at MODULE level "
                 "-> pytest COLLECTION fails on the minimal-dep ubuntu matrix "
                 "-> the whole `pytest tests/` step goes red.",
        "topics": ["ci", "radia-mcp", "test-infrastructure", "ngsolve"],
        "severity": "high",
        "first_seen": "2026-06-05",
        "last_seen": "2026-06-05",
        "what": "radia-mcp matrix 'Pytest (meta health + ... + versioning)' "
                "failed on Python 3.10/3.11/3.12: the FEMM-parity + axisym FEM "
                "cross-validation tests (test_axi_*, test_planar_*, "
                "test_scalar_fem2d*, test_stranded, test_laminated_steel, "
                "test_nonlinear_magnet, test_femm_xcheck_*) do "
                "`from ngsolve import ...` / `from netgen.occ import ...` at "
                "module top.  The matrix installs only mcp + pytest + radia-mcp "
                "(--no-deps): no ngsolve/netgen, so COLLECTION raises "
                "ModuleNotFoundError and the entire suite errors out before any "
                "test runs (origin/main 5c4b0216).",
        "root_cause": "The radia-mcp matrix is intentionally lightweight "
                      "(knowledge servers must import without the heavy FEM "
                      "stack).  A heavy test imported at module level is "
                      "collected even when it cannot run, and a collection "
                      "ImportError is fatal to the whole pytest invocation.",
        "detection": "python tools/ci_preflight.py (the radia-mcp gate runs the "
                     "suite under RADIA_MCP_FORCE_MINIMAL=1, reproducing the "
                     "ubuntu minimal-dep collection on a full-env LAB box BEFORE "
                     "push); radia-mcp-matrix.yml 'Pytest' step.",
        "prevention": "packages/radia-mcp/tests/conftest.py skips collecting "
                      "any test module that imports ngsolve/netgen when ngsolve "
                      "is absent (or RADIA_MCP_FORCE_MINIMAL=1) -- a source-scan "
                      "covers future heavy tests automatically.  Keep heavy FEM "
                      "cross-validation behind that guard; run ci_preflight "
                      "before push.",
        "related": ["tools/ci_preflight.py",
                    "packages/radia-mcp/tests/conftest.py",
                    ".github/workflows/radia-mcp-matrix.yml"],
    },
    {
        "id": "flaky-test-rerun-masked-no-rootcause",
        "title": "Intermittent test absorbed by --reruns but never "
                 "root-caused -- a 1-attempt CI red looks identical to a real "
                 "regression dismissed as 'probably flaky'.",
        "topics": ["ci", "test-infrastructure", "flaky"],
        "severity": "medium",
        "first_seen": "2026-05-20",
        "last_seen": "2026-06-05",
        "what": "The self-hosted 'Run basic tests' step is the 2nd-most common "
                "CI failure (16 of the last 80 runs).  Some are genuine "
                "regressions; some are known flakes (test_omega_reduced_omega, "
                "test_B_accuracy_inside_iron) a rerun would have saved.  With "
                "no registry the two are indistinguishable, so a real "
                "regression gets waved off as 'flaky', OR a flake burns a "
                "fix-forward cycle.",
        "root_cause": "Iterative-solver tolerance near a golden band edge / "
                      "retired moment-path dipole-in-material sampling sensitivity makes a few "
                      "assertions non-deterministic.  --reruns 2 masks them "
                      "but records neither WHICH tests flake nor WHY.",
        "detection": "tests/known_flaky.md (the registry); the workflow log's "
                     "-r aR RERUN lines.  When CI reds in 'Run basic tests', "
                     "check known_flaky.md BEFORE assuming a regression.",
        "prevention": "Keep tests/known_flaky.md current; every listed test "
                      "MUST be under --reruns.  Listing is a STOPGAP -- "
                      "root-cause it (tighten tolerance, pin seed, stabilize "
                      "mesh) and delete the row.  NEVER list a "
                      "deterministically-failing test -- that is a bug.",
        "related": ["tests/known_flaky.md",
                    ".github/workflows/build-test.yml"],
    },
    # =====================================================
    # EXAMPLE BITROT / RADIA API DRIFT / NOTEBOOK PROMOTION
    # (surfaced during the examples/ -> docs/ consolidation, 2026-06)
    # =====================================================
    {
        "id": "radia-magnetization-tesla-not-apm",
        "title": "Magnetization passed as Tesla, but Radia uses A/m "
                 "(M = Br/mu_0) -> ~zero field, silent wrong result.",
        "topics": ["radia", "units", "magnetization", "example", "bitrot"],
        "severity": "high",
        "first_seen": "2026-06-26",
        "last_seen": "2026-06-26",
        "what": "ObjCylMag/ObjHexahedron/ObjWedge given mag like [0,0,1] "
                "('1 T') return |B| ~ 1e-4 mT (microtesla) instead of "
                "hundreds of mT; or a magpylib cross-check disagrees 100%.",
        "root_cause": "Radia magnetization is A/m, NOT Tesla. [0,0,1] = "
                      "1 A/m ~ 0. A permanent magnet needs M = Br/mu_0 "
                      "(Br=1.0 T -> 7.96e5; 1.05 T -> 8.36e5; 1.2 T -> "
                      "954930 A/m). magpylib uses polarization J in Tesla, "
                      "so M_radia = J/mu_0 for a cross-check to pass.",
        "detection": "Sanity-check the probe |B| is in mT-T, not micro-T. "
                     "Real sites: smco_array, simple_problems/compare_magpylib, "
                     "background_fields, chamfered_pole_piece.",
        "prevention": "Always M = Br/mu_0 in A/m. Never copy Br (Tesla) "
                      "straight into a mag vector. See CLAUDE.md "
                      "'Magnetization Units: A/m (NOT Tesla)'.",
        "related": ["packages/radia-mcp/src/radia_mcp/magnetic_materials/permanent_magnet_knowledge.py"],
    },
    {
        "id": "rad-solve-return-4-tuple",
        "title": "rad.Solve returns a 4-tuple [residual, _, _, iterations]; "
                 "indexing [3]/[4] for max|dM|/max|dH| is wrong.",
        "topics": ["radia", "solve", "api", "example", "bitrot"],
        "severity": "medium",
        "first_seen": "2026-06-26",
        "last_seen": "2026-06-26",
        "what": "Script crashes with IndexError on solve_result[4], or "
                "reads a wrong convergence value from [3].",
        "root_cause": "radia_pybind.cpp Solve returns make_tuple(D[0..3]) = "
                      "[residual, _, _, iterations] (4 elements, indices "
                      "0-3). Old scripts assumed [3]=max|dM|, [4]=max|dH|.",
        "detection": "grep example scripts for solve_result[3]/[4]. Real "
                     "site: background_fields/{quadrupole_analytical,"
                     "sphere_in_quadrupole,permeability_comparison}.py.",
        "prevention": "Use solve_result[0] (residual) for convergence and "
                      "int(solve_result[3]) for iteration count.",
        "related": ["src/lib/radia_pybind.cpp"],
    },
    {
        "id": "radia-constructor-arity-drift",
        "title": "Old example calls miss required Radia constructor args "
                 "(ObjArcCur 8-arg, ObjHexahedron/Wedge mag, MatSatIsoFrm list).",
        "topics": ["radia", "api", "arity", "example", "bitrot"],
        "severity": "high",
        "first_seen": "2026-06-26",
        "last_seen": "2026-06-26",
        "what": "TypeError: incompatible function arguments, or "
                "'Failed to generate ...', when running an older example.",
        "root_cause": "API signatures tightened: ObjArcCur needs 8 args "
                      "(center, radii, phi, h, nseg, man_auto, axis, j) -- "
                      "old 6-arg calls omit 'man','z'. ObjHexahedron/ObjWedge "
                      "now REQUIRE a magnetization arg (and the meshed-disk "
                      "core ring is a WEDGE, not a hex). MatSatIsoFrm takes a "
                      "SINGLE nested list [[ksi,ms],...], not 3 positional "
                      "lists.",
        "detection": "Run the example; check rad.<Ctor>.__doc__ for the "
                     "current signature. Sites: simple_problems/arc_current_*, "
                     "smco_array, background_fields/*.",
        "prevention": "Read the pybind __doc__ signature before porting an "
                      "old script; ObjArcCur('man','z'); ObjWedge(verts,[0,0,0]); "
                      "MatSatIsoFrm([[...],[...],[...]]).",
        "related": ["src/lib/radia_pybind.cpp",
                    "packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/radia.py"],
    },
    {
        "id": "objmltextrtg-nonplanar-faces",
        "title": "ObjMltExtRtg fails when stacked rectangles change in BOTH "
                 "in-plane dims (non-planar side faces).",
        "topics": ["radia", "geometry", "objmltextrtg", "example", "bitrot"],
        "severity": "medium",
        "first_seen": "2026-06-26",
        "last_seen": "2026-06-26",
        "what": "'Failed to generate convex polyhedron(s)' or 'vertex points "
                "... do not belong to one plane' from ObjMltExtRtg.",
        "root_cause": "ObjMltExtRtg connects consecutive axis-aligned "
                      "rectangles with quad side faces that must be PLANAR. "
                      "Chamfering both width AND thickness (or offset centers) "
                      "makes a doubly-ruled non-planar quad.",
        "detection": "Real site: simple_problems/chamfered_pole_piece.py "
                     "(dropped). See memory reference_objmltextrtg_planar_faces.",
        "prevention": "The public ObjMltExtRtg API is retired. Build tapered "
                      "geometry with Netgen or Cubit and pass the mesh through "
                      "the NGSolve-native path.",
        "related": ["packages/radia-mcp/src/radia_mcp/radia_ngsolve/knowledge/radia.py"],
    },
    {
        "id": "example-to-notebook-promotion-breakers",
        "title": "Embedding an example script in a docs notebook breaks on "
                 "__file__, cp932 stdout rewrap, or lab-private deps.",
        "topics": ["promotion", "notebook", "example", "consolidation", "cp932"],
        "severity": "medium",
        "first_seen": "2026-06-26",
        "last_seen": "2026-06-27",
        "what": "nbconvert --execute errors: NameError on __file__; "
                "AttributeError 'OutStream' has no attribute 'buffer'; "
                "ModuleNotFoundError (mcp_server_document / magpylib / a "
                "deleted sibling).",
        "root_cause": "Notebooks have no __file__ (sys.path / output_dir / "
                      "HERE built from it break). The cp932 console rewrap "
                      "`sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer)` "
                      "fails (Jupyter OutStream has no .buffer). LAB-private "
                      "styling (mcp_server_document) + optional libs aren't in "
                      "a clean env.",
        "prevention": "When folding a script into docs/<topic>/*.ipynb: "
                      "repoint __file__-derived paths to os.getcwd() (or the "
                      "kept corpus dir for sibling imports); strip the "
                      "codecs/stdout rewrap; strip Japanese (Repository "
                      "Language); for lab-private-styled plots, DISPLAY the "
                      "committed figures instead of re-running. Verify with "
                      "`jupyter nbconvert --execute` + an error-cell + "
                      "non-ASCII scan.",
        "related": [
            "AGENTS.md: No Development Cruft in SOURCE",
            "memory/docs_cruft_prune_2026_06_28.md",
        ],
    },

    # =====================================================
    # C++ SOLVER / CROSS-SOLVE CACHE BUGS
    # =====================================================
    {
        "id": "ngsolve-gettrafo-first-touch-garbage",
        "title": "NGSolve scalar tr(ip).point occasionally returns uninitialized memory on the "
                 "first process-wide extraction",
        "topics": ["ngsolve", "flaky", "hdiv-vim", "validation"],
        "severity": "high",
        "first_seen": "2026-07-03",
        "last_seen": "2026-07-03",
        "what": "The hex-BDM1 affine wiring golden (legacy filename says RT1) flaked "
                "(~2-8% of FRESH processes, bursty and "
                "machine-load-sensitive): the demag spectrum blew up to eig 1.9..9.9 while "
                "B/M_mass stayed bit-exact.  Single coordinates of a few lattice points came "
                "back as ~1e-310 denormals -- uninitialized memory from the scalar "
                "`mesh.GetTrafo(el)(ip).point` path (the SIMD mapped-rule Assemble path is "
                "unaffected).  Under heavy load the corruption window outlasted 8 back-to-back "
                "re-evaluations.",
        "root_cause": "NGSolve's scalar element-transformation evaluation can return "
                      "uninitialized coordinates on the first touch in a process; exact "
                      "upstream mechanism not yet isolated (likely the 0xc0000374 flake "
                      "sibling).  Any geometry extraction that trusts a single evaluation "
                      "silently poisons downstream quadrature.",
        "detection": "src/radia/vim/_vim.py::_trafo_lattice_nodes DETERMINISM CONTRACT: "
                     "re-evaluate until two CONSECUTIVE evaluations agree bit-for-bit and are "
                     "finite (max 16 tries with a 2 ms decorrelation sleep), else RAISE.  The "
                     "hex wiring golden also asserts the per-instance Gram state canary "
                     "(G.hex_state_check()).",
        "prevention": "Never trust a single scalar GetTrafo lattice evaluation for committed "
                      "geometry: route every extraction through _trafo_lattice_nodes (or an "
                      "equivalent two-consecutive-bit-identical contract).  Prefer assembled "
                      "(SIMD mapped-rule) quantities where possible.  See memory "
                      "ngsolve-gettrafo-first-touch-garbage for the full forensic chain.",
        "related": ["src/radia/vim/_vim.py",
                    "validation_test/feec/test_hdiv_vim_hex_rt1_wiring.py"],
    },
    {
        "id": "fem-reference-coil-polygon-current-deficit",
        "title": "FEM-reference coil disks meshed as inscribed polygons carry ~5% too little "
                 "current -- refinement ladders do NOT reveal it",
        "topics": ["validation", "ngsolve", "geometry", "mesh-export"],
        "severity": "high",
        "first_seen": "2026-07-03",
        "last_seen": "2026-07-03",
        "what": "A VIM-vs-FEM motor torque comparison sat at ~10% with BOTH sides perfectly "
                "mesh-converged.  The FEM reference represented coil wires as J-disks "
                "(J = I / (pi r^2)); netgen meshed each circle as an inscribed ~12-gon, so the "
                "integrated current was I * A_polygon / A_circle = -4.9%.  Torque (quadratic in "
                "drive) was off -9.7%.",
        "root_cause": "The disk boundary segment count comes from curvature heuristics, not "
                      "from the volume maxh being refined -- so h-refinement ladders leave the "
                      "polygon deficit CONSTANT and both sides look converged while disagreeing."
                      "  An uncurved mesh integrates J over the polygon, not the circle.",
        "detection": "The DRIVE-EQUIVALENCE probe: rerun the FEM box with the iron/scatterer "
                     "REMOVED (pure air, linear) and compare grad A at interior probe points "
                     "against the analytic source field.  A UNIFORM percentage deficit at every "
                     "probe = source-amplitude bug; position-dependent deviation = boundary "
                     "truncation.  Demonstrated executable: "
                     "docs/electric_machine/em_reference_audit.ipynb.",
        "prevention": "In ANY FEM reference with coil disks: mesh.Curve(order) ALWAYS, and "
                      "normalize J_k = I_k / Integrate(1, definedon=Materials(coil_k)) by the "
                      "MEASURED area so INT J dA == I_k exactly regardless of meshing (see "
                      "docs/electric_machine/planar_vim_motor_helpers.py::fem_reference_bar).  "
                      "Same failure class as the Mesh Export Consistency Check policy (mesh vs "
                      "CAD volume).",
        "related": ["docs/electric_machine/planar_vim_motor_helpers.py",
                    "docs/electric_machine/em_reference_audit.ipynb",
                    "memory/fem_reference_coil_polygon_deficit.md"],
    },
    {
        "id": "truncation-ladder-frozen-edge",
        "title": "A truncation-convergence ladder that never moves ONE boundary edge converges "
                 "to the wrong answer (2D log tails make it worse)",
        "topics": ["validation", "ngsolve", "open-boundary"],
        "severity": "high",
        "first_seen": "2026-07-03",
        "last_seen": "2026-07-03",
        "what": "A reduced eddy solve's collar Dirichlet ladder 'converged' (18.94 -> 18.75 -> "
                "18.74 W/m) while the true plate loss was ~12.2: every ladder step enlarged the "
                "SIDES and BOTTOM but the TOP edge stayed frozen at y=0.1 (squeezed under the "
                "source wires), forcing A_r = 0 exactly where the plate's reaction lobe reaches "
                "toward the coil (+55% loss).  Separately, an 8x8 all-in-one box biased the same "
                "loss -15%: the 2D coil-pair potential decays only ~1/r, so its value at the box "
                "boundary was the same order as at the plate.",
        "root_cause": "Anisotropic truncation: convergence in the moved edges says nothing "
                      "about the frozen edge.  In 2D the log/1-r tails put the far boundary "
                      "much further away than 3D intuition suggests.  (The frozen top came from "
                      "a FALSE constraint -- in the reduced split A_s enters only through the "
                      "sigma*(A_r+A_s) conductor term, so the collar MAY extend past the "
                      "sources/iron: the air region never evaluates A_s.)",
        "detection": "Move EVERY edge in the ladder (one edge per step if needed).  The "
                     "COIL-ONLY (no-iron) split test separates eddy-machinery error from "
                     "constitutive error; an independent closed form (thick-plate SIBC "
                     "loss/area = R_s |H_t|^2 / 2) arbitrates.  A magnetostatic twin lesson: a "
                     "finite outer Dirichlet reflects a scattered 2D dipole and suppresses beta "
                     "by exactly 1/(1 + beta/(B0 R^2)) -- closed-form-matched in "
                     "docs/electric_machine/em_reference_audit.ipynb.",
        "prevention": "Truncation ladders must scale ALL boundary edges; budget 2D boundaries "
                      "by the log/1-r tail, not 3D decay.  For pure-harmonic exteriors use the "
                      "exact open condition instead (n=1: dA/dr + A/R = 2 B0 cos(theta) Robin).",
        "related": ["docs/electric_machine/em_reference_audit.ipynb",
                    "memory/hdiv_vim_tri_quad_motor.md"],
    },
    {
        "id": "conjugate-potential-sign-branch-cut",
        "title": "2D conjugate potential A_z from magnetic charges: the -mu0 sign inverts the "
                 "coupling channel, and the atan2 formula's per-charge branch cuts poison "
                 "SURROUNDING evaluation sets",
        "topics": ["validation", "hdiv-vim", "hodge", "magnetization"],
        "severity": "high",
        "first_seen": "2026-07-03",
        "last_seen": "2026-07-03",
        "what": "Two distinct incidents in the iron->conductor A_s channel of the planar "
                "VIM<->FEM coupling: (1) the conjugate potential was coded as "
                "-mu0 q/(2 pi) atan2(dy,dx) -- the iron appeared to SHIELD the plate (loss "
                "12.4 -> 6.7 W/m) when it actually CONCENTRATES field onto it (12.4 -> 19.6, "
                "matching the FEM); the wrong sign produced 'same magnitude, opposite sign' "
                "contributions on the two sides of the comparison.  (2) with the correct sign, "
                "a bar ring SURROUNDING a rotor core crossed the -x branch ray of every charge "
                "(the atan2 cut), giving a flat ~43% torque error at every slip.",
        "root_cause": "The correct pair for a 2D charge q is (psi, A) = (-q/(2 pi) ln r, "
                      "+mu0 q/(2 pi) atan2(dy, dx)): dA/dy = mu0 H_x and -dA/dx = mu0 H_y.  As "
                      "a FORMULA the atan2 sum is discontinuous across each charge's -x ray "
                      "even though the total field A is single-valued (zero total charge); "
                      "evaluation sets that see the charges from ONE side (e.g. strictly below) "
                      "are safe, sets that SURROUND the body are not.",
        "detection": "Permanent startup gate: B(grad A_z) == mu0 * H_at by central differences "
                     "at a few exterior points (catches any sign/derivative error in the "
                     "A-channel).  The polar closure assert (going 2 pi around must return to "
                     "the anchor, Gauss/zero-total-charge) catches cut crossings.",
        "prevention": "Use radia.vim PlanarDemagBody.Az_at only for one-sided evaluation sets "
                      "(docstring caveat); for surrounding rings build the SINGLE-VALUED polar "
                      "construction dA/dphi = mu0 r H_r anchored on the cut-free +x axis "
                      "(docs/electric_machine/planar_vim_motor_helpers.py::"
                      "Az_iron_polar_voxel).  Never ship an A-channel without the gradient "
                      "consistency gate.",
        "related": ["src/radia/vim/_vim2d.py",
                    "docs/electric_machine/planar_vim_motor_helpers.py",
                    "docs/electric_machine/em_reference_audit.ipynb"],
    },
    {
        "id": "reference-secant-picard-oscillation",
        "title": "Per-element secant-nu Picard for a nonlinear FEM REFERENCE oscillates at deep "
                 "saturation -- use the closed-form nu(B) inversion + exact Newton instead",
        "topics": ["validation", "ngsolve", "solve", "nonlinear"],
        "severity": "high",
        "first_seen": "2026-07-03",
        "last_seen": "2026-07-23",
        "what": "The all-in-one nonlinear A_z FEM reference for the salient-bar motor case "
                "never converged at deep saturation: per-element secant nu updates oscillated "
                "(dA plateaued at ~0.1 for 80..300 iterations; relax 0.5 and 0.3, and an "
                "unsafeguarded Anderson, all failed) -- corner elements kept swinging across "
                "the BH knee.  Comparing against the UNCONVERGED reference produced a plausible "
                "but wrong ~10% 'disagreement'.  RECURRED 2026-07-23 on the saturable hodograph "
                "bend forward-check: undamped Picard nu(|grad A|) stalled at rel step 2.4e-1 "
                "after 60 iterations and reported a 142% wall-field error, which read as 'the "
                "hodograph design is wrong'.  Damping omega=0.35 converged in 42 iterations to "
                "9e-10 and the SAME design then matched to 0.3%.  Severity raised to high: the "
                "failure mode is not a slow solve, it is a WRONG CONCLUSION about another method.",
        "root_cause": "The successive-substitution nu(|B|) map loses contractivity when "
                      "element fields straddle the knee (nu spans 4 orders of magnitude inside "
                      "the body); per-element averaging amplifies the swing at corners.",
        "detection": "FAIL-LOUD the reference: raise when the Picard/nu residual has not "
                     "converged rather than returning the last iterate (never compare against "
                     "an unconverged reference).",
        "prevention": "When the constitutive law admits a closed-form B->H inversion (the "
                      "saturating law does: k H^2 + (1 + chi0 - k B/mu0) H - B/mu0 = 0), write "
                      "nu(|grad A|) as an exact CoefficientFunction and use "
                      "ngsolve.solvers.Newton on the semilinear form -- quadratic convergence, "
                      "6-9 iterations from a cold start at deep saturation "
                      "(docs/electric_machine/planar_vim_motor_helpers.py::fem_reference_bar).",
        "related": ["docs/electric_machine/planar_vim_motor_helpers.py",
                    "docs/electric_machine/em_reference_audit.ipynb",
                    "memory/hodograph_design_forward_verified.md"],
    },

    # =====================================================
    # NGSOLVE / NETGEN SETUP BUGS
    # =====================================================
    {
        "id": "numberspace-gauge-densifies-direct-solve",
        "title": "A NumberSpace mean-zero gauge on a Neumann solve makes the direct "
                 "factorization DENSE -- pin one DOF instead",
        "topics": ["ngsolve", "solve", "performance", "gauge"],
        "severity": "medium",
        "first_seen": "2026-07-23",
        "last_seen": "2026-07-23",
        "what": "A gradient-recovery (pure-Neumann) solve gauged with "
                "`H1(mesh, order=3) * NumberSpace(mesh)` ran for 10+ minutes at ~400 MB "
                "on a 2k-element 2-D mesh (~20k DOFs) and produced no output; it looked "
                "like an infinite loop.  The same solve takes 0.4 s once the gauge is "
                "changed.",
        "root_cause": "The `lam*v + u*mu` coupling terms put a FULL row and column into "
                      "the assembled matrix (the single number DOF touches every H1 DOF). "
                      "A sparse direct solver then fills in almost completely, so the "
                      "factorization degenerates to a dense O(N^3) / O(N^2)-memory solve.",
        "detection": "Wall-clock and memory blow up superlinearly with mesh refinement on "
                     "a problem that should be seconds; per-stage timing logs localize it "
                     "to the gauged solve.",
        "prevention": "Fix the gauge by pinning ONE degree of freedom, which keeps the "
                      "matrix sparse and SPD: "
                      "`free = BitArray(fes.FreeDofs()); free[0] = False` then "
                      "`gf.vec.data += a.mat.Inverse(free, inverse='sparsecholesky') * f.vec`. "
                      "Only the gauge constant changes, and for potentials recovered from a "
                      "gradient only differences are physical anyway.  Reserve NumberSpace "
                      "for iterative solvers or genuinely small systems.",
        "related": ["memory/hodograph_design_forward_verified.md"],
    },
    {
        "id": "duplicate-corner-points-stall-mesher",
        "title": "Concatenating sampled boundary curves duplicates each corner and stalls "
                 "the netgen mesher with no error",
        "topics": ["netgen", "geometry", "mesh", "hang"],
        "severity": "medium",
        "first_seen": "2026-07-23",
        "last_seen": "2026-07-23",
        "what": "`SplineGeometry(...).GenerateMesh()` on a closed outline built by "
                "concatenating four separately sampled boundary curves hung indefinitely "
                "with no exception and no progress output.",
        "root_cause": "Each corner is sampled TWICE (once by each adjacent curve).  The two "
                      "samples differ by a hair (different parameter offsets), so the loop "
                      "gains a tiny near-degenerate segment that properly crosses its "
                      "neighbour.  The outline is then not a simple polygon and the mesher "
                      "cannot resolve it.  A 1 um dedup tolerance was too tight to remove it.",
        "detection": "Run an explicit self-intersection check on the closed outline BEFORE "
                     "meshing; the corner pairs show up as crossings between the first/last "
                     "segments of adjacent curves.",
        "prevention": "Build the loop from `curve1 + curve2[1:] + curve3[::-1][1:] + "
                      "curve4[::-1][1:-1]` so shared endpoints appear once, dedup with a "
                      "geometric tolerance (~10 um, not 1 um), and RAISE on any remaining "
                      "self-intersection rather than handing a bad outline to the mesher.",
        "related": ["memory/hodograph_design_forward_verified.md"],
    },
    {
        "id": "terminal-corner-contaminates-designed-body",
        "title": "A designed duct verified with its end faces attached directly to the body "
                 "reads a large MESH-INDEPENDENT 'design error' that is a terminal artefact",
        "topics": ["validation", "hodograph", "design", "ngsolve", "geometry"],
        "severity": "medium",
        "first_seen": "2026-07-28",
        "last_seen": "2026-07-28",
        "what": "A hodograph-designed saturable IPM bridge was checked by nonlinear FEM on the "
                "designed outline.  The barrier-side wall was specified at a constant 1.900 T "
                "cap but the FEM read 1.728 T near the inlet -- a 9.1 % 'design failure' -- "
                "while the rest of the wall sat within 1 %.  Refining the mesh did NOT move it "
                "(9.07 % at h/8, 9.00 % at h/16), which is exactly what makes it easy to "
                "misread as a real formulation error rather than a boundary artefact.",
        "root_cause": "The design's terminals are MMF (Psi) equipotentials, so the flux "
                      "crosses them normally, while the walls are flux lines with the field "
                      "tangential.  Where the two meet, the field must reorganise over a "
                      "boundary layer whose reach is set by the channel width, not by h -- "
                      "here about 11 degrees of turn.  The inverse map is near-degenerate at "
                      "the same corner (min |J| 2.0e-8 against a 1.2e-6 median), so the "
                      "geometry there is the least trustworthy part of the design.",
        "detection": "Refine the mesh: an error that does NOT shrink is not discretisation.  "
                     "Plot the wall error against the along-path coordinate -- a terminal "
                     "artefact is monotone-decaying from one end, a formulation error is "
                     "spread over the body.  Report where min |J| of the inverse map sits; if "
                     "it is at a terminal corner the body is fine.",
        "prevention": "Design a LEAD-IN / LEAD-OUT (here 20 degrees of flat spec at each end, "
                      "which a real part has anyway where it merges into the core) and verify "
                      "only the body.  That alone took the worst body error from 9.1 % to "
                      "1.3 %.  Never trim the contaminated samples silently -- the trim length "
                      "must be justified by a measured decay, or the check is circular.",
        "related": ["validation_test/clebsch_legendre/verify_ipm_bridge_free_boundary.py",
                    "memory/hodograph_ipm_bridge_free_boundary.md"],
    },
    {
        "id": "hodograph-wall-cusp-speed-cancellation",
        "title": "A hodograph-designed wall cusps into self-intersection where the "
                 "wall-advance speed collapses (ramp kinks, excessive local contrast)",
        "topics": ["hodograph", "design", "geometry"],
        "severity": "medium",
        "first_seen": "2026-07-28",
        "last_seen": "2026-07-28",
        "what": "The recovered cap wall of a collecting-channel design curled into a "
                "~50 um hook and the outline self-intersected.  Two distinct causes, "
                "hit on consecutive iterations: (1) pinning the wall at the cap while "
                "the opposite wall still sat at the entry level drove the local "
                "contrast H(B_cap)/H(B_low) to ~17, which demands a ~50 um turning "
                "radius; (2) with linear-in-theta wall ramps, the wall speed dropped "
                "abruptly at the EXACT ramp-end angle (slowest sample at 70.0 deg for "
                "a 70-deg ramp) and the wall crawled into a cusp.",
        "root_cause": "Along an A = const wall the advance speed per unit field angle "
                      "is |Psi_theta + Psi_B * B'|/q: the local turning radius scales "
                      "like Phi/f(rho_local), so high wall-to-wall contrast means a "
                      "tiny radius; and a DISCONTINUOUS ramp slope B' makes the speed "
                      "jump at the kink, leaving a boundary layer where the wall "
                      "barely moves.",
        "detection": "Print the per-sample wall step |dr| along each recovered wall "
                     "and the location of min |J|: a step collapsing orders of "
                     "magnitude below its median at a specific angle IS the cusp, and "
                     "the angle names the offending profile feature.  Plot the "
                     "outline before meshing.",
        "prevention": "Keep rho_local = H(B_highwall)/H(B_lowwall) <= ~5 at every "
                      "theta (the exact-annulus chart's flat region, reused as a "
                      "local rule); make wall B-profiles C1 (sin(pi t/2T) ramps, "
                      "finite start slope, flat into the plateau); stagger profile "
                      "feature angles so no two coincide.",
        "related": ["validation_test/clebsch_legendre/verify_synrm_collector_design.py",
                    "memory/hodograph_synrm_channel_rungs.md"],
    },
    {
        "id": "outline-piece-reversed-phantom-self-intersections",
        "title": "A boundary piece concatenated in the wrong direction produces dozens "
                 "of phantom self-intersections that mimic a genuine design fold",
        "topics": ["geometry", "mesh", "validation", "hodograph"],
        "severity": "medium",
        "first_seen": "2026-07-28",
        "last_seen": "2026-07-28",
        "what": "Three debugging rounds attributed self-intersection failures of a "
                "designed outline to wall cusps / folding.  Real cusp issues existed "
                "and were fixed, but the crossings that kept firing came from the "
                "EXIT face being appended to the loop in its sampling order (low "
                "corner -> cap corner) when the circuit required the reverse: the "
                "loop jumped across the face and retraced it backward, crossing the "
                "adjacent wall repeatedly.",
        "root_cause": "Assembling a closed loop from separately sampled boundary "
                      "pieces without checking that consecutive pieces CONNECT.  A "
                      "reversed piece folds the circuit; the self-intersection "
                      "checker then reports crossings that look exactly like a "
                      "genuine geometry fold, sending the debugging at the physics.",
        "detection": "Look at WHERE the crossing pairs sit: phantom crossings "
                     "cluster near one piece junction or span two specific pieces; "
                     "a real fold sits inside one wall.  Check the seam gap between "
                     "each piece's last point and the next piece's first point "
                     "against the outline extent.",
        "prevention": "Assert seam continuity when concatenating (gap << outline "
                      "extent; the promoted collector driver raises with the piece "
                      "indices).  Sibling of duplicate-corner-points-stall-mesher: "
                      "that one is about shared corners appearing twice, this one "
                      "about a whole piece traversed backward.",
        "related": ["validation_test/clebsch_legendre/verify_synrm_collector_design.py",
                    "duplicate-corner-points-stall-mesher"],
    },
    {
        "id": "axis-aligned-manufactured-solution-fake-superconvergence",
        "title": "An axis-aligned manufactured solution on a structured "
                 "mesh fakes a huge accuracy win; the derivative "
                 "seminorm is the metric immune to it.",
        "topics": ["mesh", "verification", "fem", "benchmark", "ngsolve"],
        "severity": "high",
        "first_seen": "2026-08-07",
        "last_seen": "2026-08-07",
        "what": "A mesher/element comparison shows one route beating "
                "another by 10-70x in the L2 error -- far more than any "
                "real element difference -- while the curl/div/H1 "
                "seminorm on the same runs shows near parity. Measured: "
                "structured hex looked 30x better than tet in HCurl L2 "
                "and 10-18x in HDiv L2.",
        "root_cause": "The manufactured field was separable and aligned "
                      "with the coordinate axes (u = (sin ky, sin kz, "
                      "sin kx)) and the structured mesh was aligned with "
                      "the same axes, so each edge/face dof integrated "
                      "an exactly-constant integrand -- superconvergence "
                      "of the interpolant, not solver accuracy. Rotating "
                      "the FIELD by 30 deg z / 20 deg y (which preserves "
                      "curl curl u = k^2 u, because rotation commutes "
                      "with curl and div) collapsed hex/tet L2 from "
                      "0.017-0.031x to 0.39-0.44x (HCurl) and from "
                      "0.055-0.100x to 0.93-1.01x (HDiv -- a dead heat). "
                      "It also inflated hex's CG-iteration advantage "
                      "~7x for HDiv.",
        "detection": "Compare the L2 ratio against the derivative "
                     "seminorm ratio in the SAME run: alignment moves L2 "
                     "and leaves curl/div essentially untouched (0.658 "
                     "-> 0.657 and 0.973 -> 0.967 across the rotation). "
                     "An L2 advantage with no matching seminorm "
                     "advantage is the signature. "
                     "validation_test/radia_mcp/mesh_quality_study/"
                     "run_rotation_control.py is the standing control.",
        "prevention": "Never verify or benchmark a structured mesh with "
                      "an axis-aligned separable exact solution. Rotate "
                      "the field (cheaper and cleaner than rotating the "
                      "geometry: same meshes, same exact-solution "
                      "identity), or headline the derivative seminorm. "
                      "Treat any order-of-magnitude mesh-comparison gap "
                      "as an artifact until a control says otherwise.",
        "related": ["validation_test/radia_mcp/mesh_quality_study/run_rotation_control.py",
                    "validation_test/radia_mcp/mesh_quality_study/run_vector_elements.py",
                    "memory/mesh_quality_min_is_chaotic.md"],
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
