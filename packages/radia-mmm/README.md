# radia-mmm — legacy collocation MMM/MSC demag

The **legacy** home for Radia's collocation demag solvers and their examples:

- **MMM** — Magnetic Moment Method (tetrahedra, volume magnetization `M` as DOF).
- **MSC** — Magnetic Surface Charge (hexahedra/wedges, surface charge `σ` as DOF; includes the
  Yano-Sugahara distortion elements).

## Status: legacy — HDiv-VIM is the default path

Radia's demag core is **migrating to the FEEC HDiv-VIM** (`radia.hdiv_vim`), which is the **default
path for new work** (symmetric `N = BᵀGB`, loop modes field-null by construction, μ_r-independent
convergence, no hand-crafted loop-star, robust stiff/saturating Newton on the analytic field). The
collocation MMM/MSC here is the **legacy** path, homed in this package so it can be maintained,
deprecated, and eventually removed **independently of the new core**.

- The MMM/MSC **C++ kernels stay in the `radia` wheel** (this package `depends on radia`); only the
  Python API namespace, examples, and tests live here.
- The Yano-Sugahara MSC backend is **deprecated** — see `radia.set_demag_backend("yano" | "hdiv")`.
- For new work use `radia.hdiv_vim` (`build_demag` / `DemagOperator` / `solve_demag_newton`).

## Install

```bash
pip install -e packages/radia-mmm   # editable (LAB); depends on `radia`
```

```python
import radia_mmm as rm
pm = rm.ObjHexahedron(vertices, [0, 0, 954930])   # same C++ element as radia.ObjHexahedron
B  = rm.Fld(pm, "b", [0, 0, 0.1])
```

## Examples

Migrated here (pure-Radia MMM/MSC):

- `examples/smco_magnet_array/` — SmCo permanent-magnet array (PM, no solve).
- `examples/hantila_solver/` — Hantila polarization MMM (constant-LHS, factor-once) solver.

**Migration in progress** (staged to verify each, avoid breakage). Still under the top-level
`examples/` pending review/move:

- pure MMM/MSC: `background_fields`, `simple_problems`.
- MMM-centric with an NGSolve comparison/mesh dependency (move with care): `mmm_eigenvalue_study`,
  `cube_uniform_field`, `nodal_force`, `c_type_electromagnet`.

(NGSolve / HDiv / PEEC / Kelvin / induction-heating examples stay under the top-level `examples/` —
they are not MMM/MSC.)
