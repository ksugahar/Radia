"""radia_mcp.motor knowledge: the 2D planar machine-modelling stack in Radia.

Radia's planar soft-iron path is HDiv-VIM plus shared material and
post-processing helpers.  This file is intentionally free of retired solver
names so MCP answers describe the current public API only.
"""
from __future__ import annotations

SECTIONS: dict[str, str] = {
    "overview": """\
# 2D Planar Machine Modelling In Radia

Radia uses the HDiv-VIM planar route for soft-iron demag and keeps the
surrounding planar helpers method-neutral:

- `radia.vim._vim2d` / `radia.vim.Solve`: planar HDiv-VIM soft-iron demag.
- `radia.planar_charges`: exterior field, vector potential, Maxwell torque,
  inter-body force, and fixed permanent-magnet source fields.
- `radia.planar_materials`: B-H tables, per-region material laws, and
  anisotropic susceptibility tensors.
- `radia.planar_aniso`: dense direct anisotropic solve
  `(I - X N) M = X H0` on the shared planar-charge kernel.
- `radia.planar_hysteresis`: play-hysteresis demag on the shared direct-N
  operator.

Use NGSolve `TaskManager` around planar solves.  Keep exploratory scripts in
`C:\\temp`; promote reusable checks to `validation_test/feec` or result-bearing
docs notebooks only after they are stable.
""",
    "eddy_coupling": """\
# Staggered Eddy-Current Coupling

`radia.planar_eddy.couple(...)` couples an arbitrary soft-iron solve callback
with an NGSolve reduced-potential complex `A_z` eddy FEM in a separate
conductor.  The iron solve callback has signature
`H_ext_complex (nEl, 2) -> M (nEl, 2)`.

The conductor FEM matrix is independent of the staggered iteration and should
be assembled/factored once.  The iron field is injected as the shared
`M.n` log-charge cloud rendered as an `atan2` coefficient function.
""",
    "pm_motor": """\
# Permanent Magnets

Permanent magnets are rigid fixed-M sources.  In planar HDiv workflows use the
`magnets=[(pm_mesh, M_fixed)]` route for a separate PM body.  The PM field is
added to the applied field seen by the soft iron; the PM itself is not solved.
""",
    "nonlinear": """\
# Nonlinear Soft Iron

The planar material layer parses B-H tables once and exposes monotone
anhysteretic laws for HDiv and dense planar helpers.  Nonlinear solves must
fail loudly on non-convergence and record the tolerance used; engineering
notebooks should avoid over-tight tolerances when mesh/discretization error is
already larger.
""",
    "anisotropic": """\
# Anisotropic GO-Steel Demag

`radia.planar_aniso.solve_anisotropic_demag(...)` supports scalar or
per-region `chi_par`, `chi_perp`, and `easy_deg`.  It assembles the dense
planar demag operator `N` and solves `(I - X N) M = X H0` directly.  This is the
right route for tensor susceptibility because a scalar Picard iteration is not
well-conditioned for realistic `chi`.
""",
    "hysteresis": """\
# 2D Play-Hysteresis Demag

`radia.planar_materials.PlayHysteresis` and
`radia.planar_hysteresis.solve_hysteresis_demag(...)` implement a
Prandtl-Ishlinskii style play operator on the shared direct-N demag operator.
The incremental susceptibility is non-negative, so Newton on the shared
operator stays well-conditioned; failures raise instead of falling back.
""",
    "api": """\
# API Quick Reference

```python
import radia.planar_charges as pc
import radia.planar_eddy as pe
from radia.vim import Solve

with ngsolve.TaskManager():
    r = Solve(iron_mesh, mu_r=200.0, H_ext=ngsolve.CF((H0, 0.0)))
    H = pc.exterior_field(iron_mesh, r["M"], points)
    torque = pc.maxwell_torque(iron_mesh, r["M"], Rc, H_ext=(H0, 0.0))
```
""",
    "validation": """\
# Validation Ladder

- `validation_test/feec/test_hdiv_vim_2d_magnets.py`: planar fixed-M PM source.
- `validation_test/feec/test_planar_materials.py`: canonical constitutive laws.
- `validation_test/feec/test_planar_aniso.py`: tensor susceptibility and
  embedded fixed-M PM contract.
- `validation_test/feec/test_planar_hysteresis.py`: play-hysteresis loop and
  fail-loud Newton behavior.
""",
}


def get_planar_coupling(topic: str = "overview") -> str:
    """2D planar HDiv-VIM demag + shared eddy/PM coupling knowledge."""
    t = (topic or "overview").strip().lower()
    if t in ("all", "*"):
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return "Unknown topic %r. Valid topics: %s (or 'all')." % (topic, valid)
    return SECTIONS[t]
