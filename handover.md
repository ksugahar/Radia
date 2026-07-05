# HANDOVER — purge pre-presentation manuscripts from git HISTORY

**To:** codex &nbsp; **From:** Claude &nbsp; **Date:** 2026-07-05

Transient coordination note (the previous MMMM-solver handover here was retired by
Sugahara as no longer needed). Delete this file once the purge is complete — it is
NOT a purge target, so it survives the rewrite and should be removed afterwards.

---

## 1. What to do

Sugahara decided: purge the pre-presentation paper **manuscripts** from git
**history** (not just HEAD). Execute a `git filter-repo` rewrite + force-push +
re-sync all clones. This is your domain (push / history / release / shared-.git);
Claude prepared the inventory + recipe below and did **not** run filter-repo or
force-push.

## 2. Current state (verified 2026-07-05)

- The HEAD-removal commits are already on public `origin/main`:
  `ca990656` (policy), `aee46758` (SF JA move), `a71ad521` (IGTE move),
  `4b664f5c` (URN move). These were a **normal push** — NOT a history rewrite.
- **The purge has NOT run.** All SHAs are unchanged and the manuscripts are still
  in `origin/main` history, retrievable from public:
  `git cat-file -s 4b664f5c^:docs/universal_relaxation_network/paper/urn_paper.tex`
  → `61284` bytes. So this is now a force-push **over already-public history**
  (it reduces future discoverability only; GitHub keeps orphan commits by SHA
  until GC / Support, and any fork/clone retains them).
- Include ALL in-flight commits in the pre-rewrite state (e.g. the hdiv work
  `4fc30f04` and anything newer) so the rewrite preserves them. This handover
  commit sits on top of that — clone the NAS local (below) to capture everything.

## 3. Purge target inventory (complete; historical paths across all renames)

Manuscripts moved through the examples→docs migration, the levitation→maglev
rename, and the radia-levitation dissolution, so every historical prefix must be
listed or remnants survive.

**Manuscript-ONLY directory subtrees — purge the whole subtree:**
```
docs/maglev/papers/
docs/paper/
docs/universal_relaxation_network/paper/
examples/Universal_Relaxation_Network/paper/
examples/universal_relaxation_network/paper/
examples/levitation/papers/
examples/maglev/papers/
packages/radia-levitation/papers/
```

**Loose manuscript FILES in MIXED dirs — per-file ONLY (do NOT purge the dir):**
```
examples/CLN/igte_symposium_2026.tex
examples/CLN/igte_symposium_2026.pdf
examples/CLN/igtesymp.cls
examples/CLN/A1_sibc_3stage_final.pdf
examples/CLN/A1_sibc_compact_with_AC.pdf
docs/stream_function/former_cad.md
docs/stream_function/paper_outline.md
docs/stream_function/paper_outline_sheet_metal.md
```
> ⚠️ `examples/CLN/` also holds the CLN research `.wls` scripts — purge only the
> IGTE files above, never the whole directory.
> ⚠️ CONFIRM: are `examples/CLN/A1_sibc_3stage_final.pdf` /
> `A1_sibc_compact_with_AC.pdf` paper figures (purge) or standalone CLN result
> figures referenced elsewhere (keep)? They also appear inside every
> `.../igte_symposium_2026/` manuscript dir, so they read as paper figures — but
> confirm before running.

**STAYS (NOT manuscripts — must NOT be purged):**
`docs/universal_relaxation_network/` scripts / notebooks / result-JSON /
`generate_paper_figures.py` / `results/*.tex`; `validation_test/maglev/research_cln/`
scripts; `docs/figures/lab_diagrams/*.tex`; `packages/radia-mcp/.../poster/templates/*.tex`.

## 4. Recipe (run in a FRESH clone, NEVER the shared NAS working tree)

```bash
# --- 0. make sure ALL intended commits are in the state to be rewritten ---
#     The clone below clones the NAS local working repo, which includes any
#     unpushed commits (4fc30f04 etc.) + this handover commit.

# --- 1. backup BEFORE anything (recover point) ---
git -C /s/Radia/01_GitHub bundle create /c/temp/radia_prepurge.bundle --all

# --- 2. fresh clone of the NAS local (includes unpushed commits) ---
rm -rf /c/temp/radia_purge
git clone --no-local /s/Radia/01_GitHub /c/temp/radia_purge
cd /c/temp/radia_purge

# --- 3. write the paths file ---
cat > /c/temp/purge_paths.txt <<'PATHS'
docs/maglev/papers/
docs/paper/
docs/universal_relaxation_network/paper/
examples/Universal_Relaxation_Network/paper/
examples/universal_relaxation_network/paper/
examples/levitation/papers/
examples/maglev/papers/
packages/radia-levitation/papers/
examples/CLN/igte_symposium_2026.tex
examples/CLN/igte_symposium_2026.pdf
examples/CLN/igtesymp.cls
examples/CLN/A1_sibc_3stage_final.pdf
examples/CLN/A1_sibc_compact_with_AC.pdf
docs/stream_function/former_cad.md
docs/stream_function/paper_outline.md
docs/stream_function/paper_outline_sheet_metal.md
PATHS

# --- 4. purge (git-filter-repo is installed on LAB) ---
git filter-repo --invert-paths --paths-from-file /c/temp/purge_paths.txt

# --- 5. VERIFY (all must be EMPTY / fail) ---
git log --all --oneline -- docs/universal_relaxation_network/paper/urn_paper.tex \
   docs/paper/urn_paper.tex 'examples/**/paper/urn_paper.tex'         # empty
git log --all --oneline -- docs/maglev/papers examples/CLN/igte_symposium_2026.tex # empty
git log --all --pretty=format: --name-only | \
   grep -iE 'igte_symposium_2026|urn_paper|/papers?/|former_cad|paper_outline'      # empty
git cat-file -p HEAD:docs/universal_relaxation_network/README.md | grep -n 'W:'      # README still points to W:
#   also: repo still builds + `python tools/ci_preflight.py` clean.
```

## 5. Publish (codex only)

```bash
# filter-repo drops 'origin' by design; re-add and force-push
cd /c/temp/radia_purge
git remote add origin git@github.com:ksugahar/Radia.git
git push --force-with-lease origin main
```
Then coordinate the re-sync (the disruptive part):
- **LAB + 100号機 share the NAS `.git`** — reset it to the rewritten `origin/main`
  (`git -C /s/Radia/01_GitHub fetch origin && git -C /s/Radia/01_GitHub reset --hard origin/main`),
  after confirming no un-committed / un-pushed work is stranded (esp. your own,
  and Claude may be mid-edit — coordinate the mutex).
- Re-`pip install -e` the editable packages on LAB + 100号機 if paths shifted.
- CI runner clone: reset or re-clone.
- mdx / hibino are PyPI consumers (not clones) — unaffected.

## 6. Rollback

If verification fails, do NOT force-push. Restore from
`/c/temp/radia_prepurge.bundle` (`git clone /c/temp/radia_prepurge.bundle`).

---
Inventory rationale + why-codex-not-Claude: `memory/pre_presentation_history_purge.md`
(Claude's private memory). Same recipe also at `C:\temp\HISTORY_PURGE_HANDOFF.md`
+ `C:\temp\purge_paths.txt`.
