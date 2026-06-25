# axifem — Henrotte / Meeker axisymmetric FEM examples

Examples that exercise the `radia.axifem` submodule (FEMM-style
axisymmetric finite elements with `{1, r², z}`-linear shape functions —
Henrotte 1993 / Meeker FEMM convention) on representative IH-relevant
problems.

The C++ implementation lives at `src/ext/axifem/` and ships in the
radia wheel as `src/radia/axifem.pyd`.  Element-matrix unit
tests live at `tests/axifem/`.

## Examples

| Directory | What it shows |
|---|---|
| [`disk_convergence/`](disk_convergence/) | Convergence of the first eddy-current relaxation time `τ_1` for a Cu disk vs an independent Mathematica BEM-Foster reference (208.32 µs). Documents the P1 plateau at ~5–6 % above ref and the Q2 path that closes the gap. |
| [`nmr_validation/`](nmr_validation/) | Permanent-magnet NMR axisymmetric reproduction; bundles FEMM and NGSolve mixed-formulation references so the example runs anywhere. `B_z(r)` agrees within ~0.3 % with both references except at the magnet edge discontinuity. |

Each subdirectory is self-contained: it bundles the small pure-Python
prototype it uses (axifem_core.py, sigma_mass.py) plus any reference
data files, so `python script.py` works after `pip install radia`.

## Why axifem matters for IH

Induction heating coils generate **circumferential (φ-direction)**
eddy currents in axisymmetric workpieces (cylinders, disks).  An
axisymmetric formulation reduces the problem from 3-D to 2-D `(r, z)`
without losing physics, which is the most efficient route for IH
analysis.  The Henrotte / Meeker FEMM trick of using `{1, r², z}` shape
functions removes the standard P1 axis singularity in `B_r ∝ 1/r`
and gives clean per-element `B_z = const`, `B_r ∝ 1/r` — perfect for
time-constant (CLN) extraction and Foster-form network identification.

The lab's IH research stack therefore depends on `radia.axifem`
being installable on every machine (LAB / 100号機 / mdx) so the 21 lab
users can run axisymmetric IH analyses through the standard
`pip install radia` path.
