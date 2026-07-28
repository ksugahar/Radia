"""SynRM topology optimization knowledge for radia_mcp.motor.

Distilled from:

- Kishi-Wakao-Murata-Makino-Takeuchi-Matsushita, "Multi-Objective
  Topology Optimization of Synchronous Reluctance Motors Using
  Autoencoder-Estimated Flux Barrier Shapes", IEEJ Trans. PE, 2025.
- Liu Xinyao 2025 thesis (public-safe curated corpus)
  applying the Wakao method end-to-end on a 4-pole SynRM with
  ONELAB cross-validation.

Cross-references the general framework in radia_mcp.topology_optimization
(shape derivative, topology derivative) and adds the **SynRM-specific
auto-encoder bootstrap** technique that bypasses the slow LS-method
warm-up.
"""

from __future__ import annotations


WAKAO_AUTOENCODER_LS = """\
## Wakao 2025 — Autoencoder + Level-Set hybrid for SynRM

Reference: Kishi-Wakao-Murata-Makino-Takeuchi-Matsushita, IEEJ Trans.
on Power and Energy 2025.

### The problem

Synchronous reluctance motors (SynRM) achieve torque by reluctance
variation between d-axis (low-reluctance) and q-axis (high-reluctance)
paths.  The **flux barriers** (air slots) inside the rotor shape this
reluctance distribution.  Multi-objective design targets:

- maximize average torque  T_ave(shape)
- minimize torque ripple   T_rip(shape) = (T_max - T_min) / T_ave

with binary material distribution (iron vs air per element).

Standard approach: **level-set (LS) method** (Allaire-Jouve-Toader
2004) — represent the iron/air boundary by the zero level-set of a
scalar field φ(x), evolve φ by `∂φ/∂t = -V_n |∇φ|` with shape gradient
V_n.  Robust but **slow to converge** from a uniform initial guess
because the level-set has to nucleate many disconnected flux barriers
from scratch.

### The Wakao idea

Use a **convolutional autoencoder (AE)** trained on a library of
already-good SynRM rotor shapes to suggest a near-feasible initial
shape, then refine with LS.  The autoencoder is bottlenecked through
a 2D latent space (`z_1`, `z_2`) which the authors discovered to
naturally align with the two Pareto objectives:

| Latent axis | Encoded objective |
|-------------|-------------------|
| `z_1`       | Average torque T_ave (negative correlation) |
| `z_2`       | Torque ripple T_rip (positive correlation) |

So the entire Pareto frontier in (T_ave, T_rip) space is parameterized
by a continuous 2D walk in latent space.

### Training data

- 543 SynRM rotor shapes (4-pole, 36 stator slots)
- Each shape: 100×100 binary bitmap (1 = iron, 0 = air)
- Pre-computed via JMAG / ONELAB → (T_ave, T_rip) labels per shape

### Architecture

```
Encoder:  100×100 binary → Conv(32) → Conv(64) → Conv(128) → Dense(2)
Latent:   z = (z_1, z_2) ∈ ℝ²
Decoder:  Dense(2) → ConvT(128) → ConvT(64) → ConvT(32) → 100×100 sigmoid
Loss:     binary cross-entropy + λ |z|² (regularization)
```

The **clean** correlation between latent axes and (T_ave, T_rip) is
not enforced explicitly during training — it emerges from the dataset
distribution and the bottleneck.  Tested by Kishi-Wakao by walking a
grid in z-space and labeling decoded shapes with JMAG-computed T_ave,
T_rip.

### Optimization workflow

1. **Pick** target Pareto point `(T_ave*, T_rip*)`.
2. **Decode** initial shape from `z = (z_1*, z_2*)` matched by the
   latent-to-objective regression.
3. **Refine** with level-set method (Allaire-Jouve-Toader) — only
   needs ~10-30 iterations vs ~200 for cold start.
4. **Evaluate** with high-fidelity FEA (ONELAB/JMAG/NGSolve).
5. **Loop** over Pareto front.

### Reported results

For the test case (4-pole 36-slot SynRM, 8.5 kW), the AE-LS hybrid
achieves the same Pareto front as cold-start LS in **~15% of the
compute time**.  Notably, the **shape diversity** is greater — the
LS warm-start from AE produces multiple basins (3-, 4-, 5-barrier
designs) that pure LS misses because of the fixed nucleation pattern.

### Implementation in radia_mcp.motor

The MCP-server tool `motor_topology_optimization(topic="wakao_ae_ls")`
exposes:

1. AE architecture template (PyTorch)
2. Training-set generation script using `onelab_knowledge`
3. Latent-to-objective regression recipe
4. LS-method refinement using `radia_mcp.topology_optimization.
   shape_optimization_knowledge` (gradient via adjoint)

The training-set generation is itself an MCP-orchestrated workflow:
write 543 `pmsm.geo` variants → run ONELAB driver → collect torque
waveforms → label (T_ave, T_rip) → train AE.
"""

LIU_THESIS_APPLICATION = """\
## Liu Xinyao 2025 — Wakao method applied end-to-end

Location: `public-safe curated corpus`.

### Contribution

Liu Xinyao's senior thesis re-applied the Wakao 2025 autoencoder + LS
method to a different SynRM geometry (Kindai University spec) and
verified the **cross-tool reproducibility**:

1. Train the AE on JMAG-labeled shapes.
2. Decode candidate Pareto points from latent space.
3. Cross-check decoded shapes' (T_ave, T_rip) with **ONELAB** instead
   of JMAG.
4. Verify ONELAB and JMAG agree on the Pareto front position.

### Key finding

JMAG and ONELAB produce **the same Pareto front** (within ~1%
T_ave, ~5% T_rip).  This is non-trivial because:

- JMAG uses Galerkin-type edge-element FEA (proprietary mesh).
- ONELAB uses the Whitney 1-form A-formulation with sliding
  air-gap moving-band.

The agreement validates that the discretization choice does not
contaminate the optimization landscape.  The **shape sensitivity
curve** (∂T_ave/∂shape, ∂T_rip/∂shape) is therefore topology-method-
agnostic and can be computed with NGSolve (open source, BSD license)
without losing fidelity vs proprietary tools.

### Recipe for replication in NGSolve

1. Generate 200-500 random SynRM rotor shapes via the AE-decoder
   (using a fixed seed grid in z-space).
2. For each shape: write Netgen OCC geometry script that subtracts
   the binary-bitmap air pockets from a solid rotor disk.
3. Mesh with Netgen, label per ONELAB convention
   (`rotor_iron`, `air_gap`, `stator_iron`, `stator_ind_*`).
4. Run NGSolve time-stepped magnetodynamic A-formulation
   (see `onelab_knowledge.ngsolve_xlate`).
5. Post-process torque waveform → (T_ave, T_rip).
6. Compare with the JMAG / ONELAB labels — agreement < 5%?  Good.
7. Train a per-NGSolve regressor mapping z → (T_ave, T_rip) and
   re-do step 2 of the AE workflow with the NGSolve labels.

### Open question

What is the **adjoint-mode shape derivative** for T_rip?  The standard
result for T_ave uses the magnetostatic adjoint (one extra solve per
shape).  T_rip is a **functional of the full time history** of T(t)
and requires an unrolled time-adjoint (memory cost = #timesteps ·
#DOF).  Liu Xinyao avoids this by using the AE for the dominant
contribution and only LS for the local refinement (where T_rip
sensitivity is well-approximated by finite differences).

This open question is what makes the AE warm-start so valuable for
T_rip: the AE encodes T_rip dependence directly without needing an
adjoint.
"""

PARETO_NAVIGATION = """\
## Pareto navigation in latent space

Once the autoencoder is trained, the entire Pareto front is just a
1D curve in 2D latent space.  Practical navigation:

### Step 1 — Calibrate

Walk a coarse grid (5×5) in z-space, decode each, evaluate with FEA,
fit:

  T_ave  =  a + b·z_1 + c·z_2  + ...
  T_rip  =  d + e·z_1 + f·z_2  + ...

Quadratic terms catch the local Pareto curvature.

### Step 2 — Pareto-front sampling

The Pareto front is the set of z* satisfying:

  ∂T_ave/∂z₁ · ∂T_rip/∂z₂  -  ∂T_ave/∂z₂ · ∂T_rip/∂z₁  = 0

(the Jacobian-determinant zero condition for trade-off optimality).
Solve numerically — typically a 1D root-finding along a parametric
curve.

### Step 3 — Refine with LS

Each Pareto-optimal z* gives a decoded shape; refine with LS for
local optimality (corrects for AE reconstruction error).

### Step 4 — Validate

Run a high-fidelity NGSolve time-stepped solve with **fine mesh**
(>200k elements) to confirm T_ave and T_rip on the refined shape.

### Code stub (planned)

```python
import torch
from radia_mcp.motor.topology_opt_knowledge import load_synrm_ae

ae = load_synrm_ae("models/synrm_4pole_36slot.pt")
z_grid = torch.linspace(-3, 3, 21).reshape(-1, 1)
z = torch.cat([z_grid.repeat_interleave(21, dim=0),
               z_grid.repeat(21, 1)], dim=1)
shapes = ae.decode(z)                # 441 × 100 × 100
# Then drive NGSolve / ONELAB on each shape...
```
"""

SATURABLE_BRIDGE_HODOGRAPH = """\
## Saturable rotor bridge / rib: shape the iron with the field-plane hodograph

VERIFIED 2026-07-28 --
`validation_test/clebsch_legendre/verify_ipm_bridge_free_boundary.py`
(golden-banded, exits nonzero on violation; sidecar
`results_ipm_bridge_free_boundary.json`).

### When this applies

The IPM bridge / rib is the one place in a machine where saturation is the
DESIGN INTENT, not a nuisance: the iron neck is meant to saturate so leakage
flux is capped and predictable.  The engineering constraint is therefore a CAP

    |B| must nowhere exceed the knee B_knee

and the cheapest iron is the iron that sits EXACTLY at the cap on its
most-loaded surface -- no hot spot, no wasted margin.

In physical space that surface is a FREE BOUNDARY (an unknown curve on which a
field condition holds), so locating it means an outer shape loop whose every
iteration is a nonlinear solve.  In the field-plane (Chaplygin) hodograph
`(B, theta)` the SAME condition is the COORDINATE LINE `B = B_knee`, i.e. a
fixed Dirichlet edge, and the design collapses to ONE LINEAR SOLVE.  The
saturating law enters as a known coefficient of the independent variable:

    d/dB( a A_B ) + b A_thth = 0,   a = B mu_d / mu_s^2,   b = 1/(mu_s B)

with `mu_s(B) = B/H` the measured SECANT curve (no table inversion) and
`mu_d = dB/dH` the differential.  In log-polar `u = ln B` the operator is
Laplace stretched by exactly `mu_d / mu_s` -- linear material gives 1
(conformal); saturation IS the departure from conformality.  For silicon steel
near 1.9 T that ratio is ~0.21, a domain stretch of only ~2.2, so the solve is
well conditioned.  Ellipticity never fails: `(mu q)' = dB/dH > 0` for any
monotone B-H curve, so there is no magnetic analogue of the sonic line.

### What the verification measured

90-degree turn around a barrier tip, inner (barrier-side) wall pinned at the
1.900 T cap, outer wall ramping 0.90 -> 1.70 T (the funnel closing), leakage
flux 1.1e-3 Wb/m, representative steel `mu_r(B) = 1 + 6999/(1+(B/1 T)^4)`.
Designed by one linear hodograph solve, then checked by an INDEPENDENT
nonlinear FEM on the designed outline:

| quantity | hodograph design | naive baseline |
|---|---|---|
| body inner wall abs(B) vs the 1.900 T cap | 1.875--1.902 T | 1.329--2.156 T |
| cap overshoot | +0.10 % | +13.50 % |
| inner-wall spread | 1.40 % | 43.56 % |
| body iron area, same flux | 0.7590 mm^2 | 0.8049 mm^2 (+6.0 %) |
| MMF (independent global check) | 2.4447 A design vs 2.4476 A FEM (0.12 %) | -- |

The naive baseline is a circular centreline through the same two body end
mid-points, same turn, same two widths, same lead-in/out, same flux -- what a
competent engineer draws, not a straw man.

### Two facts you can build on

1. **Flux scale-freedom is EXACT.**  The equation is linear in `A` and `Psi` is
   linear in `A`, so a fixed FIELD spec determines the geometry up to one scale
   that is exactly proportional to the flux (measured deviation 2e-12 at half
   flux).  Practical consequence: solve ONCE, then scale to the
   mechanically-set throat width and read the leakage flux straight off.  A
   flux sweep costs nothing.
2. **Terminals contaminate the body -- always design a lead-in / lead-out.**
   With the end faces attached directly to the body, the inlet corner corrupts
   the first ~11 degrees of the inner wall by up to 9.1 %, MESH-INDEPENDENTLY
   (9.07 % at h/8, 9.00 % at h/16).  A 20-degree flat lead-in/out at each end
   (which a real bridge has where it merges into the core) drops the worst body
   error to 1.3 %.  See bug pattern `terminal-corner-contaminates-designed-body`.

### Scope -- read before quoting this

- This is a **LOCAL** design kernel.  The bridge's boundary data (flux, turn
  angle, terminal field levels) comes from the surrounding rotor solve, and the
  throat width is a MECHANICAL input (centrifugal stress), not an output.  The
  intended architecture is: global nonlinear FEM -> extract the local boundary
  data -> hodograph designs the local shape -> substitute -> re-solve.
- The claim is NOT that the hodograph beats FEM.  A nonlinear shape loop
  converges to the same shape.  The claim is that the CAP becomes a boundary
  condition instead of an outcome to check afterwards, and the loop of
  nonlinear solves collapses to one linear solve.
- **Where it does NOT apply**: slots (current-carrying, so not source-free);
  magnet interiors (fixed `M`, so `B` and `H` are not collinear and there is no
  `mu(|H|)`); the air gap (`mu = mu0`, where the hodograph degenerates to a
  conformal map -- that is classical Schwarz-Christoffel / Carter territory,
  nothing new); and the cross-section AS A WHOLE, because field nulls on the
  d/q symmetry lines are singular points of `(B, theta)` and flux barriers make
  the iron multiply connected so `theta` is not single-valued.  One tooth or
  one bridge, away from nulls, is the right unit.

See also `radia_mcp.electromagnet` topic `clebsch_hodograph` for the
formulation and the 90-degree bend case it generalises.
"""

FLUX_CHANNEL_HODOGRAPH = """\
## SynRM flux channel: saturated sizing chart + free-form collector design

VERIFIED 2026-07-28 -- two golden-banded drivers in
`validation_test/clebsch_legendre/`:
`verify_synrm_channel_annulus_lock.py` (part 1) and
`verify_synrm_collector_design.py` (part 2), with committed JSON sidecars
and the outline figure `synrm_collector_outlines.png`.

### Part 1 -- the pure turn is exactly solvable; use the chart, not FEM loops

A 90-degree turning channel with flux-line walls and equipotential terminals
is exactly solvable for ANY material law: curl H = 0 with azimuthal H forces
H = C/r regardless of mu(B).  The cap-binding family reduces to quadrature:

    rho = H_cap/H(B_out),  f(rho) = INT_1^rho B_of_H(H_cap/s) ds,
    Phi = r_in f(rho),     body area = (Theta/2) r_in^2 (rho^2 - 1)

Numbers for the representative steel (cap 1.90 T; percentages are
material-model-dependent, the (B,H) samples are in the JSON):

- area/Phi^2 minimum at rho* = 5.83, FLAT within 5% over rho in [3.8, 9.5]
  -- pick rho by layout, the iron penalty is small.
- radial width w/Phi is nearly aspect-independent (0.55..0.67 per T): the
  q-axis budget cost of a channel barely depends on the turn aspect; what
  varies (5x) is the iron AREA, i.e. loss volume.
- The LINEAR-designed annulus (B ~ 1/r sizing, linear-optimal rho = 2.22)
  wastes 26.8% of the cap (peak 1.391 T) and uses 2.67x the optimal iron
  (sizing-with-the-real-curve x0.518, then aspect x0.722).  The linear 1/r
  width rule overestimates channel width by +18% (rho 1.5) to +59% (rho 3).
  The practitioner's B_avg ~ 1.5 T experience rule is far better (~+10% at
  rho = 3); the chart replaces the experience constant with a guaranteed
  cap-binding value.

The hodograph adds NOTHING here (the constant-wall solution IS the annulus);
its role in part 1 is to be LOCKED against this exact nonlinear reference
(r_in/rho/MMF to 2.6e-6..4.6e-5, circularity ~1e-8) -- the strongest
validation of the design machinery in the lane.

### Part 2 -- the collecting channel: where free-form shaping pays

The real SynRM channel collects flux DISTRIBUTED along the gap-side face.
That kills H = C/r (no quadrature) and the hodograph becomes the only
linear-cost tool.  New BC class: entry face = the hodograph segment
B = B_e with the Dirichlet ramp A(theta) = Phi theta/theta_c; barrier-side
wall pinned at the cap through the carrying turn.  Verified by independent
nonlinear FEM: flat-cap region mean 0.032% / max 0.068%, peak 1.900 T;
MMF agrees at the |B| level x the saturation slope (see rule 4).

Payoff: the compass-drawn baseline through the SAME entry face and SAME
exit (circular-arc walls -- what an engineer sketches) uses +78% iron
(3.29 vs 1.84 mm^2) and peaks at only 1.371 T (28% of the cap unused).
With fixed terminals, circular walls cannot follow the accumulating flux.
SynRM translation: same d-axis flux, same terminals, -44% channel iron ->
directly thicker barriers.

### Four reusable design rules (each cost one failed iteration)

1. **Wall assignment is forced by non-crossing in (B, theta)**: the cap wall
   is the LONG wall attached at the theta = 0 corner of the entry face.
   The reverse assignment makes the two wall images cross (domain pinch).
2. **Local contrast rule**: keep rho_local = H(B_capwall)/H(B_lowwall) <= ~5
   at EVERY theta (the part-1 chart reused as a local rule).  At ~17 the
   wall demands a ~50 um turning radius and cusps into self-intersection.
3. **Ramps must be C1**: the wall-advance speed along an A = const wall is
   |Psi_theta + Psi_B B'|/q, so a discontinuous ramp slope B' puts a cusp
   exactly at the kink.  Use sin(pi t/2T) ramps; stagger feature angles so
   no two profile features coincide.
4. **MMF comparisons amplify B errors by mu_s/mu_d**: on the saturating
   curve dH/H = (mu_s/mu_d) dB/B (~4.3 at 1.45 T).  Band MMF checks at the
   |B| band times this slope, or a healthy design reads as a failure.

### Scope -- read before quoting

- LOCAL design kernel: entry level B_e, flux, and turn come from the
  surrounding rotor solve; matching a specific stator MMF harmonic content
  at the face is future work (the design's required gap loading dA/ds is
  reported: near-uniform ~ B_e here).
- The baseline is ONE compass construction (arcs through the same
  anchors); the claim is the mechanism (circular walls cannot track
  accumulating flux), not superiority over every hand method.
- Embedding in a full rotor and measuring L_d/L_q/saliency is the open
  next rung; until then the -44% iron is a channel-level, not a
  machine-level, number.

See also topic `saturable_bridge_hodograph` (the free-boundary cap trick on
the leakage bridge) and `radia_mcp.electromagnet` topic `clebsch_hodograph`
(formulation and the verified 90-degree bend).
"""

SECTIONS = {
    "wakao_ae_ls": WAKAO_AUTOENCODER_LS,
    "liu_thesis_application": LIU_THESIS_APPLICATION,
    "pareto_navigation": PARETO_NAVIGATION,
    "saturable_bridge_hodograph": SATURABLE_BRIDGE_HODOGRAPH,
    "flux_channel_hodograph": FLUX_CHANNEL_HODOGRAPH,
}


def get_topology_opt_knowledge(topic: str = "wakao_ae_ls") -> str:
    """Return SynRM topology-optimization knowledge.

    Args:
        topic: section name; "all" returns everything.
    """
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
