# Single-stroke chain construction

A real wound coil must be **one continuous conductor** driven by a
single current source.  The iso-contours of the stream function are
*N closed loops* — to manufacture, they must be connected into one
chain.  The connection segments carry the full current and produce a
**parasitic field** that was not part of the original SF design.

This page documents the three single-stroke methods shipped + the
empirical complexity tier framework that bounds their performance.

## The three chain methods

All three are in
[`demo_sf_to_peec_gx.py`](../../examples/stream_function/demo_sf_to_peec_gx.py),
selectable via `--chain-method {greedy, lobe, kuijpers}`.

### `greedy`: global nearest-neighbour

For each polyline, find the closest unvisited polyline (in periodic
`(φ, z)` distance) and connect via a cylinder-surface geodesic.  Each
polyline is opened at the point closest to the previous chain end.

**Pros**: shortest total wire length.
**Cons**: chain criss-crosses the cylinder, visually obvious "wasted"
arcs (the user's flag 2026-05-30: 「等高線同士を無駄に接続している」).

Measured (planar Gx fingerprint baseline): RMS 21.78 %, x-axis
nonlin 11.40 %, length 22.65 m, 1350 segments.

### `lobe`: 4-quadrant Maxwell-pair topology

Classify each closed contour into one of four saddle quadrants
`(sign(x_centroid), sign(z_centroid))`.  Within each quadrant, sort by
polyline arc length (outer first) and chain via greedy NN.  Inter-
quadrant transitions follow a fixed order `(+x,+z) → (+x,-z) → (-x,-z)
→ (-x,+z)` with alternating direction so the exit of one quadrant is
near the entry of the next.

**Pros**: only 3 inter-quadrant arcs (vs N−1 for greedy).
**Cons**: slightly longer than greedy because forced quadrant ordering.

Measured: RMS 23.98 %, x-axis nonlin 9.67 %, length 24.67 m, 1350
segments.

### `kuijpers`: per-lobe cut line with straight rungs *(default, recommended)*

After Kuijpers, Jansen, Lomonova (Compumag 2023 [525]) Method-1.  Each
lobe is given a **fixed cut φ** (`+x` lobes at φ=0, `-x` at φ=π).
Every contour in a lobe is opened at the point closest to its lobe's
cut φ.  Contours within a lobe are sorted by the z coordinate of their
cut point; adjacent contours are connected by **one straight (φ, z)
"rung" at the cut**.  Inter-lobe transitions use the same straight
blend.

This is the Fig.4 rung pattern of Kuijpers et al.

**Pros**: cleanest visual structure (no criss-cross), best DSV RMS, 23%
fewer segments (rungs are 1 segment vs 8-segment geodesic helix).
**Cons**: slightly longer wire than greedy.

Measured: RMS **16.24 %** (-25 % vs greedy, -32 % vs lobe),
nonlin 9.73 %, length 23.99 m, **1040 segments**.

Confirms the paper's observation "deviation region = field error
region" — minimising the deviation minimises the field error.

## Complexity tier framework

A coil's REACHABLE design quality is bounded by its TOPOLOGY:

| Tier   | Topology                  | Baseline    | + Path-A     | Path-A effectiveness         |
| ------ | ------------------------- | ----------- | ------------ | ---------------------------- |
| EASY   | axisym. / planar uniform  | RMS 2–3 %   | < 1 %        | useful (basis-loop)          |
|        |   + FE-direct H¹ ψ        | 2 %         | **0.47 %**   | **MONOTONE convergence**     |
| MEDIUM | cylindrical Gz            | already OK  | redundant    | smooth helix natural         |
| HARD   | cylindrical Gx fingerprint| 16 %        | 15 %         | tier-bounded                 |
| HARDER | shielded / biplanar / 3D  | —           | —            | needs FE-direct or D-path    |

Past EASY needs FE-direct continuous ψ.  Past HARD needs B-spline SFD
(Kuijpers Methods 2/3) or multivalued-potential reformulation (Path D
in `aca_tsvd(topic=single_stroke)`).

## Dead-end variants (do NOT re-try)

Six attempts on 2026-05-30 that all *reduced* field accuracy:

**Chain-orientation traps** (memory entry
`feedback_single_stroke_chain_orientation_traps`):

  - Sort by sign(centroid_x): RMS 16 → **43 %** (lobe reversal flips
    current direction for half the loops).
  - Force top/bottom cut on closed contours: RMS 16 → **30 %** (same
    flip).
  - Pair-aware reverse (`pos_y` outer-to-inner then `neg_y` reversed
    inner-to-outer): RMS 16 → **40 %**.
  - 8-sub-lobe (sx, sy, sz): RMS 16 → 35 %.

Common cause: any per-contour traversal-direction reversal flips
matplotlib's natural CCW orientation → current direction inverts →
field broken.

  **SAFE knob** = global within-lobe SORT KEY.
  **UNSAFE** = per-contour TRAVERSAL DIRECTION change.

**Path-A on basis-loop ψ** (memory entry
`feedback_path_a_naive_picard_negative`):

  - α=1.0 full step: 16 → 55 % (divergence).
  - α=0.5 / 0.1: oscillates 40–55 % (no contraction).
  - α=0.05 + N=60 iters: oscillates, best 15 % by accident.

Common cause: matplotlib contour topology jumps under small ψ
perturbations → `B_c(ψ)` not smooth → Picard does not contract.

## Why FE-direct ψ unsticks Path-A

`H¹` GridFunction is a CONTINUOUS function of its DOF vector → contour
family deforms SMOOTHLY → matplotlib emits a consistent contour count
and ordering → `B_c(ψ)` smooth → naive Picard CONTRACTS.

Empirical proof (planar uniform Bz):

  iter 40: residual 0.62 %
  iter 41: 0.60 %
  iter 42: 0.58 %
  ...
  iter 47: **0.49 %**       ← monotone, no backtracking

Same iteration on basis-loop ψ oscillates 1.7 % → 15 % around the best.

## Visualising the chain

```bash
python examples/stream_function/view_sf_coil_gx_gmsh.py --mode contours
python examples/stream_function/view_sf_coil_gx_gmsh.py --mode chain    # = kuijpers
python examples/stream_function/view_sf_coil_gx_gmsh.py --mode step
```

`--mode contours`: shows the SF design's raw closed contours (no
connection arcs).  Each contour is a separate Physical Group in GMSH
so you can toggle them individually.

`--mode chain`: shows the single-stroke chain (= what PEEC sees).  Use
to inspect the rung pattern and identify any criss-cross artefacts.

`--mode step`: opens the multi-piece loft chain STEP from
`demo_sf_to_peec_gx.py --with-peec`.

(Pitfall: GMSH may restore the window to an off-screen monitor; the
viewer scripts force `General.GraphicsPositionX/Y` to (100, 100).  If
you still see no window, see `gmsh_usage(topic=pitfalls)` #9.)

## Cross-reference

  - Math: [theory.md](theory.md) Path-A section
  - Open extension paths: [paper_outline.md](paper_outline.md) Path A/B/D
  - MCP topic: `aca_tsvd(topic=single_stroke)`,
    `aca_tsvd(topic=session_2026_05_30)` sections 2 (chain methods) +
    3 (Path-A representation) + 4 (dead-end variants).
