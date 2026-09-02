"""TaskManager knowledge -- NGSolve work-stealing parallelism for the
Radia + NGSolve stack.

Read this when:
* Writing a new `calc_*.py` panel script that does any NGSolve solve.
* Reviewing an existing `calc_*.py` for parallelisation correctness.
* Confused by why MKL `mkl_set_num_threads` and ngsolve's TaskManager
  both seem to control thread count and yet disagree.
* Implementing a custom C++ kernel that should obey CLAUDE.md
  "Parallelisation: NGSolve TaskManager (CLAUDE.md TaskManager policy)".
* Adding a `--nthreads` CLI flag to a calc_*.py.

The MCP server exposes this via ``taskmanager(topic=...)``.  Topics:
overview, usage, set_num_threads, mkl_interaction, cpp_kernels,
audit_radia_ih, common_mistakes, all.
"""


TASKMANAGER_OVERVIEW = r"""
# NGSolve TaskManager — work-stealing parallelism for Radia

## Policy (CLAUDE.md)

> **Use NGSolve TaskManager for thread-level parallelisation, NOT raw
> OpenMP parallel regions.**

NGSolve's `TaskManager` provides work-stealing task-based parallelism
that integrates with MKL and avoids nested-OpenMP issues.  All new
parallel code in Radia should use it.

## What it does

`TaskManager` runs Python / C++ NGSolve operations multi-threaded:

- `BilinearForm.Assemble()` — assembly is parallelised over elements.
- `Mesh.Curve(...)` — curving is parallelised over curved elements.
- Krylov solvers (`CGSolver`, `GMRESSolver`) — MatVecs and inner
  products parallelise.
- `Integrate(...)` — domain quadrature parallelises over elements.
- `Set(...)` projection onto a GridFunction (L2 projection).
- ``add_ngsolve_python_module``-built C++ TUs that use
  `ngcore::ParallelFor` (e.g. `radia.sparsesolv_ngsolve`,
  `rad_equivalence_source.cpp`).

## Without TaskManager

Without the context manager, NGSolve falls back to **single thread**:
assembly, solve, integrate all run on one core.  On a 28-core LAB
machine you lose ~20× wall-clock for free.  This is the most common
"why is my solve so slow" cause in panel scripts.

## See also

- `taskmanager('usage')` — the `with TaskManager():` pattern + when
  to wrap.
- `taskmanager('set_num_threads')` — `ngsolve.SetNumThreads(n)` +
  `--nthreads` CLI convention.
- `taskmanager('mkl_interaction')` — when MKL and TaskManager both
  want to set thread count.
- `taskmanager('cpp_kernels')` — `ngcore::ParallelFor` from custom
  C++ TUs (compact_netgen, equivalence_source, etc.).
- `taskmanager('audit_radia_ih')` — concrete audit of the IH headless solver
  paths (status: all PASS as of 2026-05-27).
- `taskmanager('common_mistakes')` — assemble outside TaskManager,
  missing `--nthreads` CLI, double-wrapping, etc.
"""


TASKMANAGER_USAGE = r"""
# `with TaskManager():` pattern

## Canonical usage

```python
from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                    TaskManager, ...)
import ngsolve

# Optional: cap thread count for the whole script (default = all cores)
ngsolve.SetNumThreads(8)

# ---- wrap every solve in `with TaskManager():` ----
with TaskManager():
    bf = BilinearForm(fes)
    bf += nu_cf * grad(u) * grad(v) * dx
    bf.Assemble()                          # parallelised

    lf = LinearForm(fes)
    lf += J0 * v * dx("coil")
    lf.Assemble()                          # parallelised

    gfu = GridFunction(fes)
    gfu.vec.data = bf.mat.Inverse(
        fes.FreeDofs(), inverse="pardiso") * lf.vec   # parallelised
```

## What counts as "the solve" — wrap from BilinearForm() to .Set()

Empirically (verified 2026-05-27 audit):

- `BilinearForm` definition (the `+=` lines): cheap, can be inside or
  outside.  Convention: inside.
- `.Assemble()`: **MUST** be inside.  Parallelises over elements.
- `Mesh.Curve(...)`: **MUST** be inside.  Geometry projection per element.
- `mat.Inverse()` / `CGSolver` / `GMRESSolver`: **MUST** be inside.
  Direct solver internal parallelism + iterative solver MatVecs.
- `Integrate(...)` post-processing: **SHOULD** be inside (parallel
  domain integration); negligible if obs is point.
- `GridFunction.Set(...)` (L2 projection): **MUST** be inside.
- `gf.Save(path)`: outside is fine (pure I/O).
- VTK / GMSH export: outside is fine (pure I/O).

## Nesting is OK but redundant

`with TaskManager(): with TaskManager(): ...` is a no-op for the
inner context — NGSolve detects and reuses the outer thread pool.
No need to defensively re-wrap inside a helper.

## Single-thread debug mode

For reproducible debugging:

```python
ngsolve.SetNumThreads(1)
with TaskManager():
    bf.Assemble()
# Now Assemble runs serially -- numerical results are deterministic
# (assembly order matters for FP accumulation at the ULP level).
```
"""


TASKMANAGER_SET_NUM_THREADS = r"""
# `ngsolve.SetNumThreads(n)` and the `--nthreads` CLI convention

## The setter

```python
import ngsolve
ngsolve.SetNumThreads(8)        # use 8 threads
ngsolve.SetNumThreads(0)        # default (= all cores; n_cpu - 1 on Linux)
ngsolve.SetNumThreads(1)        # serial / debug
```

Call BEFORE entering `with TaskManager():`.  The thread pool size is
fixed when the context manager is entered.

## Panel script convention (`calc_*.py`)

Add a CLI flag so the panel UI / a benchmarking script can override:

```python
parser.add_argument(
    "--nthreads", type=int, default=0,
    help="TaskManager threads (0 = auto = all cores; 1 = serial)")

# ... later ...

if args.nthreads:
    ngsolve.SetNumThreads(args.nthreads)
with TaskManager():
    ...
```

**Convention** (verified in `calc_fem_kelvin.py` line 359-362):
- Default `0` = NGSolve picks (= all cores)
- Explicit `1` = serial (debugging / golden tests)
- Explicit `> 1` = exact thread cap

## Default thread count by machine

| Machine | Cores | NGSolve default |
|---------|-------|-----------------|
| LAB (Threadripper Pro) | 28 | 28 |
| 100号機 (multi-user)   | 16 | 16 (shared) |
| mdx (compute node)     | 28-56 | 28-56 |
| CI runner (GH Actions) | 4 | 4 |

On 100号機 (shared with other users), prefer `--nthreads 8` or `12` to
avoid starving other users.  LAB and mdx: default is fine.

## Reading the current thread cap

```python
import ngsolve
print(ngsolve.ngsglobals.numthreads)      # 0 means "default = all cores"
```

If you see `numthreads == 1` in production and didn't ask for it,
something earlier in the import chain called `SetNumThreads(1)` --
look for `os.environ["OMP_NUM_THREADS"] = "1"` in dependency setup
scripts.
"""


TASKMANAGER_MKL_INTERACTION = r"""
# TaskManager × Intel MKL — who controls what

## Two thread pools

NGSolve and MKL each have their own thread pool:

1. **NGSolve TaskManager** — wraps the solver's outer parallel loop
   (element assembly, MatVec dispatch, etc.).
2. **MKL** — runs the matrix kernels (DGEMM, sparse LU
   factorisation, Pardiso, BLAS-3) — controlled by
   `mkl_set_num_threads()`.

## The nesting problem

If TaskManager spawns 28 tasks AND each calls a 28-thread MKL DGEMM,
you get 28 × 28 = 784 threads competing for 28 cores → context
switching kills throughput.

## The convention NGSolve uses

NGSolve resolves this by setting MKL to **single-thread inside the
TaskManager region** and using its own outer parallel loop:

```
TaskManager: 28 workers
  worker i:  Assemble element i (calls MKL DGEMM in single-thread mode)
```

So inside `with TaskManager(): bf.Assemble()`, MKL sees
`mkl_set_num_threads(1)` — automatically.

## When MKL parallelism IS useful

For the **Pardiso solve** of a large factored system:

```python
with TaskManager():
    bf.Assemble()                          # 28 NGSolve workers, MKL=1
    gfu.vec.data = bf.mat.Inverse(
        fes.FreeDofs(), inverse="pardiso"  # ← MKL Pardiso, controlled
    ) * lf.vec                              # by mkl_set_num_threads
```

NGSolve's TaskManager hands off the factorisation to Pardiso, which
internally uses MKL parallelism.  Pardiso reads
`mkl_set_num_threads()` AT FACTORISATION TIME, so set it BEFORE the
`with TaskManager():` block:

```python
import mkl
mkl.set_num_threads(8)             # for the Pardiso factorise step
ngsolve.SetNumThreads(8)           # for the assemble step
with TaskManager():
    ...                            # both honour 8 threads
```

(In practice the two values can differ; the lab uses
`ngsolve.SetNumThreads(args.nthreads)` and lets MKL default to all
cores for Pardiso.)

## Diagnostics

If wall-clock is suspiciously high:

```python
import time, mkl, ngsolve
print(f"MKL: {mkl.get_max_threads()}")
print(f"NGSolve: {ngsolve.ngsglobals.numthreads}")
t0 = time.perf_counter()
with TaskManager():
    bf.Assemble()
print(f"Assemble: {time.perf_counter()-t0:.2f}s")
```

If Assemble takes 30 s on a 28-core machine for a 100k-DOF FEM,
something is wrong — likely TaskManager not active (= the wrap was
forgotten one level up).
"""


TASKMANAGER_CPP_KERNELS = r"""
# TaskManager in custom C++ — `ngcore::ParallelFor`

## Policy (CLAUDE.md)

> **All new parallel code in Radia C++ should use TaskManager.**
> Raw OpenMP `#pragma omp parallel for` is allowed only in legacy
> code that hasn't been migrated.

## Pattern

```cpp
#include <core/taskmanager.hpp>

void EvaluateStaticH(
    const double* centroids, ...,
    int n_obs,
    double* H_out,
    int n_threads)
{
    // Optionally cap threads (default = NGSolve global)
    int requested = (n_threads > 0)
        ? n_threads
        : ngcore::TaskManager::GetMaxThreads();
    if (requested <= 0) requested = 1;
    ngcore::RegionTaskManager rtm(requested);

    ngcore::ParallelFor(ngcore::IntRange(n_obs), [&](size_t j) {
        // ... compute H at obs[j] ...
        H_out[j*3 + 0] = ...;
        H_out[j*3 + 1] = ...;
        H_out[j*3 + 2] = ...;
    });
}
```

Reference implementation: `src/core/rad_equivalence_source.cpp`
(Phase A 2026-05-26, 50-120× over Python).

## Calling from Python

The C++ side accepts `int n_threads`.  Python wrapper passes 0
(NGSolve default) unless the user overrides:

```python
H = nfs.evaluate_static_H(obs, use_cpp=True, n_threads=0)
```

## Why NOT raw OpenMP

`#pragma omp parallel for` opens its own thread pool that:
- Doesn't honour `ngsolve.SetNumThreads(...)` from the calling
  Python.
- Nests badly with NGSolve's own TaskManager (28 × 28 = 784 threads).
- Doesn't integrate with MKL's automatic single-thread fallback.

`ngcore::ParallelFor` honours all of the above for free.

## When to release the GIL

If the C++ kernel is bound via pybind11 and runs for more than a few
microseconds, release the GIL so other Python threads (e.g. the
panel GUI) stay responsive:

```cpp
m.def("_EquivalenceSourceStaticH",
    [](py::array_t<double> centroids, ...) {
        // ... unpack args ...
        py::gil_scoped_release release;        // ← release GIL
        radia::eqsrc::EvaluateStaticH(...);    // C++ runs without GIL
        // GIL re-acquired automatically here
    });
```

Verified in `src/lib/radia_pybind.cpp` `_EquivalenceSourceStaticH`
binding (Phase A).
"""


TASKMANAGER_AUDIT_RADIA_IH = r"""
# Audit: TaskManager coverage in IH headless solvers

Audit date: **2026-05-27** (commit ~HEAD).
Scope: every `calc_*.py` script owned by `IHDesignSpec` and the IH Simulink
application.

## Result table

| Script | IH role | TaskManager status |
|---|---|---|
| `calc_inductance.py` | Four inductance modes | ✅ OK — BIE assembly is in helper `radia.bem_sibc_solver` (line 235 has `with TaskManager():`) |
| `calc_peec.py` | Vacuum-coil mode | ✅ OK — no NGSolve solve (pure scipy/numpy + PEEC topology); the selected BLAS owns its threads |
| `calc_fem_kelvin.py` | PEEC+FEM+Kelvin | ✅ OK — 5 × `with TaskManager():` wraps; `--nthreads` CLI; `SetNumThreads()` honoured |
| `calc_fem_coilmesh.py` | Full FEM | ✅ OK — 3 × `with TaskManager():` wraps |
| `calc_heat.py` | Thermal 3D | ✅ OK — `with TaskManager():` wrap |
| `calc_heat_axisym.py` | Thermal axisymmetric | ✅ OK — `with TaskManager():` wrap |
| `calc_heat_with_em_table.py` | Thermal/EM table validation | ✅ OK — `with TaskManager():` wrap |

## Conclusion

**All IH headless NGSolve paths are correctly TaskManager-parallelised.**

## Helper modules (also covered)

| Helper | TaskManager status |
|---|---|
| `src/radia/bem_sibc_solver.py` | ✅ line 235 `with TaskManager():` for BIE assembly |
| `src/radia/esim_coupled_solver.py` | n/a — pure scipy / numpy (cell problem) |
| `src/radia/esim_cell_problem.py` | n/a — pure scipy / numpy |

## How to re-audit

```bash
# From the repository root:
rg -n "TaskManager|\.Assemble\(\)|\.Inverse\(.*inverse" \
     src/radia/panels/calc_inductance.py \
     src/radia/panels/calc_fem_kelvin.py \
     src/radia/panels/calc_fem_coilmesh.py \
     src/radia/panels/calc_heat*.py \
     src/radia/panels/calc_peec.py
```

Any `Assemble()` line without a `with TaskManager():` ancestor in the
same function = parallelism bug. Re-run after changes to an IH solver path.
"""


TASKMANAGER_COMMON_MISTAKES = r"""
# Common TaskManager mistakes (and how to spot them)

## 1. Assemble outside TaskManager

```python
# WRONG — Assemble runs serially
bf.Assemble()
with TaskManager():
    gfu.vec.data = bf.mat.Inverse(...) * lf.vec

# CORRECT
with TaskManager():
    bf.Assemble()
    gfu.vec.data = bf.mat.Inverse(...) * lf.vec
```

How to spot: `grep -n "\.Assemble\(\)" calc_*.py` and confirm each
hit is inside a `with TaskManager():` block in the same function.

## 2. Forgetting `--nthreads` CLI

Without the flag, the script always uses `ngsolve.ngsglobals.numthreads`
which is set by whatever ran first in the process.  On 100号機 this
can be 1 (left over from another user's debug session).  Always:

```python
parser.add_argument("--nthreads", type=int, default=0,
                    help="TaskManager threads (0 = auto, 1 = serial)")
...
if args.nthreads:
    ngsolve.SetNumThreads(args.nthreads)
```

## 3. Wrapping Mesh() / mesh.Curve()

`Mesh(path)` itself is single-threaded I/O and shouldn't be wrapped
(no benefit, just noise).  `mesh.Curve(p)` IS parallel and MUST be
wrapped:

```python
mesh = Mesh(vol_path)                     # outside
with TaskManager():
    mesh.Curve(args.order)                # inside
    bf.Assemble()
    ...
```

## 4. Setting threads INSIDE the context

```python
# WRONG — SetNumThreads is ignored once you're inside the TM context
with TaskManager():
    ngsolve.SetNumThreads(8)              # no effect
    bf.Assemble()                         # uses whatever was set before

# CORRECT
ngsolve.SetNumThreads(8)                  # before
with TaskManager():
    bf.Assemble()
```

## 5. Multiple Pardiso solves not benefiting

Pardiso's REUSE-symbolic-factorisation path is faster on a second
solve with the same matrix structure.  TaskManager doesn't fix this
automatically — you need `inverse="pardiso"` + Sparskit reuse on
your own.  TaskManager only parallelises within ONE Pardiso call.

## 6. Mixing raw OpenMP and TaskManager in same C++ TU

If a custom C++ file has `#pragma omp parallel for` AND
`ngcore::ParallelFor` in the same call chain, the inner OMP region
fights with NGSolve's pool → context-switch storm.  Pick one
(NGSolve TaskManager per CLAUDE.md policy).

## 7. cubit-spawned subprocess inheriting wrong thread count

When Cubit's embedded Python launches `calc_*.py` via subprocess.run,
the child inherits Cubit's environment.  If Cubit sets
`OMP_NUM_THREADS=1` (it sometimes does on Windows), all child
NGSolve threads run serial.  Workaround in the spawning code:

```python
env = os.environ.copy()
env.pop("OMP_NUM_THREADS", None)        # let NGSolve decide
env.pop("MKL_NUM_THREADS", None)
subprocess.run([python, "calc_x.py", ...], env=env)
```
"""


TASKMANAGER_HELPER_VS_CALLER = r"""
# Policy: Caller Wraps, Helper Does NOT (2026-05-27)

## The rule

| Where the code lives | `with TaskManager():` ? |
|---|---|
| `src/radia/panels/calc_*.py` (CALLER) | YES — wrap the FEM block |
| `examples/**/*.py`           (CALLER) | YES — wrap the FEM block |
| `tests/**/*.py`              (CALLER) | YES — wrap the FEM block |
| `src/radia/**/*.py` (HELPER, NOT calc_*) | NO — caller wraps |
| C++ kernels using `ngcore::ParallelFor` | n/a — honour caller's context |

## Why "caller wraps, helper does not"

1. **Composability**: a caller running `helper_a()` → `helper_b()` →
   `helper_c()` opens ONE TM region for all three; no start/stop
   overhead × 3.
2. **Predictability**: the parallelism intent is visible at the call
   site in `calc_*.py`, not buried inside a helper module.
3. **Single audit point**: `tools/audit_taskmanager.py` only needs
   to grep callers for missing wraps + helpers for stray wraps; no
   call-graph traversal.
4. **No visual noise**: `with TaskManager(): with TaskManager(): ...`
   is a no-op for the inner context (NGSolve detects and reuses) but
   it muddies the source.

## Migration history (2026-05-27)

The 7 helper modules under `src/radia/` that historically had
internal `with TaskManager():` were converted to no-internal-wrap:

| Helper | Sites removed | Notes |
|---|---|---|
| `bem_sibc_solver.py` | 1 | BIE assembly — major helper |
| `bem/coil_inductance_ngsolve.py` | 1 | inductance from src/sink |
| `kelvin_solver.py` | 2 | `_assemble_and_solve` |
| `kelvin_validate.py` | 1 | reference-energy compare |
| `scalar_bie_sibc.py` | 1 | scalar BIE-SIBC class |
| `scalar_potential_solver.py` | 1 | Ω-reduced Ω scalar |
| `vector_potential_solver.py` | 4 | A-formulation HCurl |

Callers (`calc_*.py` + `examples/**.py` + `tests/**.py`) were swept
to add the caller-side wrap.  154 example files (mostly
adaptive-Kelvin convergence sweeps) gained a `with TaskManager():`
around `main()` or the top-level execution block.

## Helper diagnostic (optional, defensive)

Some helpers may add an assertion that yells if the caller forgot
to wrap.  This is OPTIONAL — `audit_taskmanager.py` catches the
missing-wrap case at lint time:

```python
def _assemble_and_solve(a_bf, f_lf, fes, inverse="pardiso"):
    \"\"\"Assemble and solve.  Caller MUST be inside `with TaskManager():`
    per CLAUDE.md "Caller Wraps, Helper Does NOT" (2026-05-27).
    \"\"\"
    a_bf.Assemble()
    f_lf.Assemble()
    gfu = GridFunction(fes)
    gfu.vec.data = a_bf.mat.Inverse(
        fes.FreeDofs(), inverse=inverse) * f_lf.vec
    return gfu
```

Note the docstring: explicit caller-contract.  No assertion code
(NGSolve's API doesn't expose "am I in TM" — only timing-based
heuristics are possible, which are noisy in CI).

## Caller pattern (canonical)

```python
import argparse
import ngsolve
from ngsolve import (..., TaskManager, ...)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--nthreads", type=int, default=0,
                    help="TaskManager threads (0=auto, 1=serial)")
    args = p.parse_args()

    if args.nthreads:
        ngsolve.SetNumThreads(args.nthreads)

    mesh = Mesh(args.vol)
    with TaskManager():
        mesh.Curve(args.order)
        # ... ALL helper calls, ALL FEM solves go here ...
        result = solve_workpiece(mesh, ...)            # helper, no wrap
        gfu_E.Set(result.E_cf, dual=True)              # parallel Set()
        # ...

    # Outside: I/O, JSON dump, etc.  Cheap, single-threaded OK.
    json.dump(result, sys.stdout)

if __name__ == "__main__":
    main()
```

## Audit tool

```bash
python tools/audit_taskmanager.py            # repo-wide
python tools/audit_taskmanager.py --quiet-ok # violations only
python tools/audit_taskmanager.py --json     # machine output
```

Exit 0 = clean.  Non-zero = list of violations with file:line and
the suggested fix.  Run after editing any solver helper or example.

## See also

- CLAUDE.md "TaskManager Wrap Policy: Caller Wraps, Helper Does NOT"
- `.claude/skills/taskmanager-wrap/SKILL.md` (skill alias for editors)
- `taskmanager('audit_radia_ih')` — IH panel audit result
- `taskmanager('common_mistakes')` — old pitfalls
- `taskmanager('only')` — 2026-05-27 expansion: TaskManager is the
  SOLE parallelization mechanism, even for not-yet-parallel NGSolve ops
"""


TASKMANAGER_ONLY = r"""
# Policy: TaskManager-Only (2026-05-27)

## The rule

**NGSolve's `TaskManager` is the SOLE parallelization mechanism
for any Radia code that touches NGSolve.**  Do NOT introduce
`threading`, `multiprocessing`, `concurrent.futures`, raw `#pragma
omp parallel`, or any other parallel runtime for NGSolve-driven
computation.  Radia C++ kernels reach the same threadpool via
`ngcore::ParallelFor`, which honours the active TaskManager.

**Even NGSolve operations that are NOT (yet) internally parallel
go through TaskManager.**  Wrap them anyway.

## What "go through TaskManager" covers (expanded list)

| Op                            | Currently parallel inside NGSolve? | Wrap anyway? |
|---|---|---|
| `BilinearForm.Assemble()`     | YES                                | YES |
| `LinearForm.Assemble()`       | YES                                | YES |
| `mat.Inverse(inverse="...")`  | YES (PARDISO etc.)                 | YES |
| `mesh.Curve(p)`               | YES                                | YES |
| `gf.Set(cf, ...)`             | YES (CF evaluation)                | YES |
| `Integrate(cf, mesh, ...)`    | YES                                | YES |
| `CGSolver/GMRESSolver/...`    | YES (per-iter mat-vec)             | YES |
| `BilinearForm(fes)` ctor      | NO (mostly serial setup)           | **YES** |
| `LinearForm(fes)` ctor        | NO                                 | **YES** |
| `H1 / HCurl / HDiv / L2 / VectorH1 / Periodic / Compress / ...` | partial | **YES** |
| `HDivSurface / HCurlSurface`  | NO                                 | **YES** |
| `GridFunction(fes)` ctor      | NO                                 | **YES** |
| `Mesh(geo.GenerateMesh(...))` | NO (Netgen mesher is single-thread) | **YES** |
| `Norm(...)`, `InnerProduct(...)` inside `Integrate` | YES | YES |

The "wrap anyway" column is the policy.  Why:

1. **Future-proof.**  When NGSolve adds parallelism to FES
   construction or quadrature in a future release, callers benefit
   without code changes.
2. **One coherent block.**  Reader sees "the NGSolve computation"
   bracketed by one `with TaskManager():`, not a patchwork of
   per-op micro-wraps.
3. **No "is this op parallel?" gotcha.**  Contributors don't have
   to memorize which 2026-05 NGSolve calls happen to be parallel.
4. **Aligns with NGSolve's own tutorials.**  Canonical NGSolve
   examples wrap the whole solve block, not individual lines.

## What is NOT NGSolve and does NOT need wrapping

- Pure `numpy` / `scipy` operations
- File I/O (`open`, `json.dump`, `pickle.dumps`, ...)
- `argparse`, `logging`, `print`
- `radia.*` Python calls that ultimately invoke Radia C++ kernels
  (those kernels use `ngcore::ParallelFor` which honours the active
  TaskManager — you still wrap them for the same future-proof reason,
  but they are not what triggers the audit)

## Forbidden alternatives

- `import threading` for any FEM-numeric loop
- `import multiprocessing` for assembly / solve / mesh
- `concurrent.futures.ThreadPoolExecutor` for NGSolve work
- `#pragma omp parallel` in new C++ that calls into NGSolve

If you need PROCESS-level parallelism (sweep over parameter sets,
each running its own NGSolve session), use `subprocess` / a job
scheduler / `joblib` only AT THE OUTERMOST level — each child
process opens its own `TaskManager` for its own NGSolve work.

## Audit (expanded 2026-05-27)

`tools/audit_taskmanager.py` flags the broader set of NGSolve ops
listed above.  76 caller violations were uncovered when the
expanded check went live; 73 were auto-fixed by
`tools/fix_taskmanager_callers.py`, 3 required manual care due to
nested `def` in the wrap region.

```bash
python tools/audit_taskmanager.py            # expect "All clean"
```

## Canonical script template

```python
import argparse, sys, json
import ngsolve
from ngsolve import *           # brings TaskManager + everything

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--nthreads", type=int, default=0,
                    help="0 = auto, 1 = serial")
    p.add_argument("--vol", required=True)
    args = p.parse_args()
    if args.nthreads:
        ngsolve.SetNumThreads(args.nthreads)

    with TaskManager():
        # everything NGSolve-touching goes here
        mesh = Mesh(args.vol)
        mesh.Curve(3)
        fes  = H1(mesh, order=3, dirichlet="outer")
        u, v = fes.TnT()
        a = BilinearForm(fes); a += grad(u)*grad(v)*dx; a.Assemble()
        f = LinearForm(fes);   f += 1*v*dx;             f.Assemble()
        gfu = GridFunction(fes)
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(),
                                      inverse="sparsecholesky") * f.vec
        err = Integrate((grad(gfu)-grad_u_exact)**2, mesh)
        result = {"err": err}

    # OUTSIDE TM: cheap I/O
    json.dump(result, sys.stdout)

if __name__ == "__main__":
    main()
```

## See also

- CLAUDE.md "TaskManager-Only Policy: Align with NGSolve, No
  Alternatives"
- `taskmanager('helper_vs_caller')` — original "caller wraps" rule
- `taskmanager('mkl_interaction')` — how TM and MKL share threads
- `taskmanager('cpp_kernels')` — `ngcore::ParallelFor` in custom C++
"""


_TOPIC_MAP = {
    "overview":         TASKMANAGER_OVERVIEW,
    "usage":            TASKMANAGER_USAGE,
    "with":             TASKMANAGER_USAGE,                   # alias
    "set_num_threads":  TASKMANAGER_SET_NUM_THREADS,
    "threads":          TASKMANAGER_SET_NUM_THREADS,         # alias
    "nthreads":         TASKMANAGER_SET_NUM_THREADS,         # alias
    "mkl_interaction":  TASKMANAGER_MKL_INTERACTION,
    "mkl":              TASKMANAGER_MKL_INTERACTION,         # alias
    "pardiso":          TASKMANAGER_MKL_INTERACTION,         # alias
    "cpp_kernels":      TASKMANAGER_CPP_KERNELS,
    "cpp":              TASKMANAGER_CPP_KERNELS,             # alias
    "parallel_for":     TASKMANAGER_CPP_KERNELS,             # alias
    "audit_radia_ih":   TASKMANAGER_AUDIT_RADIA_IH,
    "audit":            TASKMANAGER_AUDIT_RADIA_IH,          # alias
    "radia_ih":         TASKMANAGER_AUDIT_RADIA_IH,          # alias
    "common_mistakes":  TASKMANAGER_COMMON_MISTAKES,
    "mistakes":         TASKMANAGER_COMMON_MISTAKES,         # alias
    "pitfalls":         TASKMANAGER_COMMON_MISTAKES,         # alias
    "helper_vs_caller": TASKMANAGER_HELPER_VS_CALLER,        # NEW 2026-05-27
    "policy":           TASKMANAGER_HELPER_VS_CALLER,        # alias
    "caller_wraps":     TASKMANAGER_HELPER_VS_CALLER,        # alias
    "helper":           TASKMANAGER_HELPER_VS_CALLER,        # alias
    "only":             TASKMANAGER_ONLY,                    # NEW 2026-05-27 (expanded)
    "taskmanager_only": TASKMANAGER_ONLY,                    # alias
    "tm_only":          TASKMANAGER_ONLY,                    # alias
    "no_alternatives":  TASKMANAGER_ONLY,                    # alias
    "broad":            TASKMANAGER_ONLY,                    # alias
    "future_proof":     TASKMANAGER_ONLY,                    # alias
}


def get_taskmanager_knowledge(topic: str = "overview") -> str:
    """Dispatch TaskManager knowledge topics.

    Topics:
        overview         - Policy + what TaskManager does (DEFAULT)
        usage            - `with TaskManager():` pattern, when to wrap
        set_num_threads  - `SetNumThreads()` + `--nthreads` CLI convention
        mkl_interaction  - TaskManager vs MKL thread pool nesting
        cpp_kernels      - `ngcore::ParallelFor` in custom C++ TUs
        helper_vs_caller - "Caller wraps, helper does NOT" (2026-05-27)
        only             - TaskManager-only policy: SOLE parallel path,
                           wrap even not-yet-parallel NGSolve ops (2026-05-27)
        audit_radia_ih   - 2026-05-27 audit of IH panel solver scripts
        common_mistakes  - Wrap-outside-Assemble, missing --nthreads, etc.
        all              - Everything
    """
    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join([
            TASKMANAGER_OVERVIEW,
            TASKMANAGER_USAGE,
            TASKMANAGER_SET_NUM_THREADS,
            TASKMANAGER_MKL_INTERACTION,
            TASKMANAGER_CPP_KERNELS,
            TASKMANAGER_HELPER_VS_CALLER,
            TASKMANAGER_ONLY,
            TASKMANAGER_AUDIT_RADIA_IH,
            TASKMANAGER_COMMON_MISTAKES,
        ])
    if topic in _TOPIC_MAP:
        return _TOPIC_MAP[topic]
    return (f"Unknown topic '{topic}'. Available: "
            "overview, usage, set_num_threads, mkl_interaction, "
            "cpp_kernels, helper_vs_caller, only, "
            "audit_radia_ih, common_mistakes, all.")
