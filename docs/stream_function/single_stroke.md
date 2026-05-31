# Single-stroke chain construction

A real wound coil must be **one continuous conductor** driven by a
single current source.  The iso-contours of the stream function are
*N closed loops* — to manufacture, they must be connected into one
chain.  The connection segments carry the full current and produce a
**parasitic field** that was not part of the original SF design.

This page documents the three single-stroke methods shipped + the
empirical complexity tier framework that bounds their performance.

## The chain methods

All are in
[`demo_sf_to_peec_gx.py`](../../examples/stream_function/demo_sf_to_peec_gx.py),
selectable via
`--chain-method {field_aware, kuijpers, lobe, greedy, nn_blend}`.

> **Use the [`single-stroke-chain`](../../.claude/skills/single-stroke-chain/SKILL.md)
> skill** when joining a new coil's contours.  The single-stroke
> connection has no clean closed-form optimum — it is a reason-and-verify
> task (pick a method that respects the current-sign order → build the
> chain → measure DSV RMS → keep only if it beats the baseline → escalate
> to Path-A).  Five deterministic "improvements" have been tried and most
> made the field WORSE; the skill encodes the traps so you do not repeat
> them.

### `field_aware`: sign-order + azimuthal-min cuts *(default, recommended)*

Uses the `kuijpers` lobe/current-sign visiting ORDER (the dominant factor),
but instead of snapping every cut to a fixed φ, chooses each contour's cut
by coordinate descent to **minimise the azimuthal arc to its chain
neighbours** (axial `dz` is free — it carries no stray field).  This lets
the rungs' stray fields cancel more symmetrically over the DSV.

**Pros**: best DSV RMS found; beats `kuijpers` in every tested config.
**Cons**: slightly more, slightly scattered rungs than `kuijpers`.

Measured (default nphi=24 nz=40 nlevels=12): RMS **9.29 %**, x-axis
nonlin **7.20 %** — vs `kuijpers` 16.24 % / 9.73 %.  Robust: 30–54 % lower
RMS than `kuijpers` across nlevels 10/12/16 and nphi32/nz48.

**Key subtlety**: the azimuthal rung *total* is NOT the predictor —
`field_aware` has MORE azimuthal arc than `kuijpers` (12522 vs 7444 mm) yet
LOWER RMS, because field impact is the symmetric *cancellation* of the rung
stray fields, not their summed length.  And the geometry-only
`nn_blend` (12847 mm azimuthal, near-equal) gives RMS 0.65 because its
order interleaves the lobe signs.  The lever is the sign-respecting ORDER,
then symmetric cut placement — neither reducible to a single scalar metric
(hence the skill).

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

### `kuijpers`: per-lobe cut line with straight rungs *(prior best, robust fallback)*

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

## Avoiding the long bridges entirely: multi-wire

For a multi-lobe coil (Gx fingerprint) the long inter-lobe "rungs" the
single stroke needs are **structurally irreducible** — one continuous
wire MUST traverse all four lobes, so it must bridge between them.  Two
attempts to remove them both FAILED (2026-05-31):

  - **Reorder (boustrophedon / snake)** so consecutive lobes meet at the
    saddle: did NOT shorten the worst rung (still 291 mm) and made the
    field WORSE (9.3 % → 19.2 %) — the rung arrangement's symmetric
    stray-cancellation breaks.
  - **Cut the chain at the long rungs** into sub-wires: produces OPEN
    current paths (current can't start/stop mid-air; `div J ≠ 0`), so the
    Biot-Savart field is unphysical (9.3 % → 21 %).

The only physical way to have NO bridges is for every independent wire to
be a **closed loop** — i.e. drive each closed contour as its own
conductor (multiple current feeds).  That is exactly
[`demo_coil_design_gx.py`](../../examples/stream_function/demo_coil_design_gx.py):
independent saddle-shaped closed loops, **no bridges, DSV RMS 0.81 %** —
an order of magnitude better than any single stroke.

| Design | wires / feeds | bridges | DSV RMS |
|--------|---------------|---------|---------|
| single-stroke (`demo_sf_to_peec_gx.py`, field_aware) | 1 | yes (irreducible) | 9.3 % |
| multi-wire (`demo_coil_design_gx.py`, independent loops) | N (one per loop) | **none** | **0.81 %** |

**RULE**: if the wasteful long connections are unacceptable, switch to
the independent-closed-loop (multi-wire) design — do NOT try to cut or
reorder a single stroke to remove them.  The single stroke buys ONE feed
at the cost of the bridges + ~10× worse field; the multi-wire design buys
the best field + no bridges at the cost of N feeds.

## Complexity tier framework

A coil's REACHABLE design quality is bounded by its TOPOLOGY:

| Tier   | Topology                  | Baseline    | + Path-A     | Path-A effectiveness         |
| ------ | ------------------------- | ----------- | ------------ | ---------------------------- |
| EASY   | axisym. / planar uniform  | RMS 2–3 %   | < 1 %        | useful (basis-loop)          |
|        |   + FE-direct H¹ ψ        | 2 %         | **0.47 %**   | **MONOTONE convergence**     |
| MEDIUM | cylindrical Gz            | already OK  | redundant    | smooth helix natural         |
| HARD   | cylindrical Gx fingerprint| 16 % (kuijpers) → **9.3 % (field_aware)** | **8.1 % (field_aware + Path-A)** | chain method matters; Path-A composes |
| HARDER | shielded / biplanar / 3D  | —           | —            | needs FE-direct or D-path    |

Past EASY needs FE-direct continuous ψ.  **The HARD tier "16 % ceiling"
was a `kuijpers`-method artifact, not a fundamental bound**: the
`field_aware` chain (sign-order + azimuthal-min cuts, 2026-05-31) reaches
9.3 % on the same SF design with no Path-A — so the connection METHOD, not
just the SF design, gates HARD-tier quality.  The two best ideas COMPOSE:
`field_aware` + Path-A (`--compensated-iter 40 --compensated-step 0.3`)
reaches **8.1 %**, roughly half the old `kuijpers` baseline.  Going lower
still likely needs B-spline SFD (Kuijpers Methods 2/3) or multivalued-
potential reformulation (Path D in `aca_tsvd(topic=single_stroke)`).

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
