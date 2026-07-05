"""Stream-function (SF) coil-design knowledge for the radia-streamfunction MCP
server.

The detailed, self-contained knowledge lives in THIS server's own
``radia_mcp.streamfunction.knowledge.aca_tsvd`` (overview / method / api /
kernel_agnostic / performance / cmaes / validation / literature / workflow /
single_stroke / regularized -- the last carries the folded-Tikhonov + Pareto
front + sheet-metal lever section).  It was MOVED here from
``radia_mcp.radia_ngsolve`` in 2026-06 because SF coil design is not a general
NGSolve usage and was making radia_ngsolve too large.  This module adds a
dedicated SF OVERVIEW + topic map and dispatches every other topic to
``get_aca_tsvd_knowledge`` (which has ~80 aliases).
"""

from .knowledge.aca_tsvd import get_aca_tsvd_knowledge


SF_OVERVIEW = r"""
================================================================
Stream-Function (SF) Coil Design -- radia-streamfunction
================================================================

The SF method designs a surface coil that produces a TARGET field over a
design volume (DSV).  A surface current ``K`` is written from a scalar
stream function ``psi`` on the conductor surface,

    K = n_hat x grad psi          (planar: K = z_hat x grad psi)

so |K| = |grad psi| is the surface current density and the iso-contours of
psi at equal increments ARE the wire paths.  Designing the field reduces to
solving the (massively underdetermined) least-norm system

    A psi = B_target        (M target constraints, N >> M surface DOFs)

with the (ACA+)+TSVD kernel-agnostic solver, then choosing WHICH least-norm
psi via a regularisation.

WHAT THIS FRAMEWORK PROVIDES
----------------------------
  * Kernel-agnostic (ACA+)+TSVD least-norm solver -- the matrix entry
    A(i,j) is a callback, so the SAME solver drives coils (Biot-Savart) OR
    magnets (Radia MMM material kernel).  [topic: kernel_agnostic, method]
  * FE-direct psi as a continuous H1 GridFunction on ANY surface
    (plane / cylinder / sphere / conformal / 3D-printed former) -- the case
    the structured basis-loop grid cannot represent.  [topic: single_stroke]
    MULTI-SURFACE works with no special-casing: a BIPLANAR coil (two plates in
    ONE mesh) designs through the same path; abe groups EACH disconnected
    component's edges into its own free constant (ndof_free < ndof) and the
    contours close on every plate.  Locked by test_streamfunction_biplanar.
  * Regularisation / DESIGN OBJECTIVE menu (--regularize): L2 (min |psi|) /
    H1 (min |grad psi|, a smoothness proxy) / INDUCTANCE (min 1/2 psi^T L psi
    -- the PHYSICAL min-stored-energy gradient-coil objective, Turner/Forbes;
    L = mu0 C^T SL C from ngsolve.bem LaplaceSL, K = n x grad psi; validated
    torus -0.6 %, Nagaoka solenoid 0.78; dense, moderate N) + L-inf peak cap
    (pareto lever).  All fold onto ONE ACA factorisation via
    RegularizedTSVD.from_stiffness(base, S).  [topic: regularized]
  * INDUCTANCE: MINIMISE (gradient coil, fast slew) OR TARGET A VALUE (IH
    resonance).  --target-inductance L_H (or --resonance-cap C with --peec-freq
    f -> L_target = 1/((2 pi f)^2 C)) SEARCHES nlevels (turns, L_coil ~ N^2; the
    field design is nlevels-independent) by bisection of single-stroke -> PEEC
    L_coil for the turn count that resonates; reports `resonance` {nlevels,
    achieved L, resonance_freq_Hz, L_range_H, in_range}.  Same BEM L machinery,
    opposite goal.  LAB: C=22nF f=200kHz -> 28.8uH -> nlevels 13 -> 30.3uH ->
    coil resonates 195kHz; connects SF designer to radia-ih.  WITH --greedy-turns
    (few-turn coil) the turns are FIXED by greedy (no nlevels search; L_coil ~ N^2
    pins N), so it instead reports required_cap_F = 1/((2 pi f)^2 L_coil) -- the
    tank cap to resonate the DELIVERED low-turn coil at --peec-freq (+ for
    --resonance-cap, the resonance_freq_Hz that coil+cap reach).  [topic: panel,
    low_turn]
  * Folded TIKHONOV (the "+ alpha I" core) -> the (field homogeneity, PEAK
    current density) PARETO FRONT swept at ~50 us/point, with FOUR stackable
    levers: Tikhonov alpha (L-curve), L-inf minimax, geometry (former size /
    cylinder length), and sheet-metal surface forming.  [topic: pareto]
  * Single-stroke chain -- connect the contour family into ONE
    continuous wire with least stray-field impact (field_aware / Kuijpers).
    [topic: single_stroke]
  * Sheet-metal, TWO distinct kinds:
      - WIRE distortion: bend the manufactured single wire (psi fixed) to
        cancel the single-stroke residual -- one feed, no extra shims.
      - SURFACE forming: reshape the conductor SURFACE + re-solve psi to lower
        the PEAK current density on the Pareto front (genuine bending,
        standoff-vs-bending decomposition).  Lever direction is
        geometry-dependent: planar out-of-surface, cylinder in-surface (axial
        for low-m / Z2 targets, +azimuthal for m>=2 ellipse C2 / S2 shims).
    [topic: pareto, single_stroke]
  * Surface-deformation outer loop (CMA-ES / NSGA-II) for accuracy.

PIPELINE
--------
  target field  ->  A psi = B   ((ACA+)+TSVD, regularised)
                ->  psi on the surface  ->  equal-increment iso-contours
                ->  single-stroke chain (one wire)  ->  CAD / PEEC / manufacture
  Optional outer loops: regularisation-shape (sigma), geometry / sheet-metal
  forming, surface deformation -- to push the (homogeneity, peak-J) front.

TOPIC MAP  (query: streamfunction("<topic>"))
---------------------------------------------
  overview          this page
  theory  / method  SFM + (ACA+)+TSVD math
  api               radia.stream_function API (aca_tsvd, RegularizedTSVD, ...)
  kernel_agnostic   the callback matrix-entry contract (coils OR magnets)
  regularized       regularisation menu + folded Tikhonov closed form
  pareto            (homogeneity, peak-J) Pareto front + 4 levers + sheet-metal
  material_aware    SF design WITH IRON (--iron-vol): Kelvin-FEM reaction folded
                    into the whole pipeline; obs-adjoint; exact-quad source
  low_turn          low/many-turn discrete refinement: greedy constructive
                    (--greedy-turns, monotone) + connector-aware; --pin-tiling
                    driven array; --optimize-levels; --distort; single-pass
  single_stroke     one-wire chain, FE-direct on arbitrary formers, wire forming
  cmaes             SA-25-020 CMA-ES outer loop
  performance       ACA+ amortisation numbers
  validation        analytic-benchmark checks
  literature        SFM lineage (Turner / Peeren / current potential)
  workflow          end-to-end demo recipes
  panel             the design/pareto/manufacture PANEL + calc_streamfunction.py
  boundary_conditions  --confine off/on/abe (Abe edge-equipotential BC)
  contour           contour=flux-line; --contour-sub order-p + --flux-plot bubble
  harmonics         spherical-harmonic target (--target-harmonic Z2/X/..) +
                    achieved-field (l,m) decomposition: purity / contamination
  shielding         active shielding (--shield-vol): primary+shield joint
                    design (Turner 1986), ~86x stray reduction; coverage lesson
  fusion            stellarator Stage-2 (REGCOIL/NESCOIL/FOCUS): winding-surface
                    current potential, net current, force/stress, VMEC, winding-shape;
                    + radia-SF-vs-REGCOIL SUPERSET (parity + deliverable/iron, goldens a/b/c)

DOCS + DEMOS
------------
  docs/stream_function/  (README, theory, regularization, single_stroke,
    deformation, examples, api, benchmarks)
  SHOWCASE NOTEBOOKS (executed, outputs embedded -- the method run live):
    docs/stream_function/demo_gallery.ipynb    -- public demo gallery for the
       remaining human-facing examples; validation/benchmark runners are split
       to validation_test/stream_function, with reusable computation promoted
       to src/radia APIs
    docs/stream_function/theory.ipynb          -- build_fem_matrix on plane /
       cylinder / sphere + achievable RMS; Path-A residual-vs-iteration (honest
       non-monotone with auto-levels); complexity-tier table from the 4 demos
    docs/stream_function/regularization.ipynb  -- the {l2,h1,h1_sigma,
       inductance_diag,linf} menu (refreshes the comparison table); H1 order
       sweep p=1..5 (p=3 sweet spot, non-monotone); folded RegularizedTSVD
       checks (S=I == pseudo_inverse ~1e-15; alpha>0 == dense Tikhonov)
    docs/stream_function/deformation.ipynb     -- run_deformation_search live:
       penalty form (--minimize-reg) + NSGA-II Pareto front (--pareto);
       cost-vs-trial + (RMS, psi^T S psi) front; flat/accuracy/reg-min table
    docs/stream_function/examples_catalog.ipynb -- source/result migration
       ledger for transitional demo scripts, synchronized with JSON sidecars.
    (the .md keeps the prose/derivations; the .ipynb embeds the produced
    numbers + figures; validation and benchmarks execute from validation_test)
  Demo records include the Pareto + sheet-metal cases:
    demo_pareto_tikhonov_aca / demo_pareto_geometry_nsga /
    demo_pareto_cylinder / demo_pareto_deform / demo_pareto_cylinder_deform

The detailed knowledge lives in this server's own knowledge module
(``radia_mcp.streamfunction.knowledge.aca_tsvd``, moved from radia-ngsolve
in 2026-06) -- this server is the SF-focused front door over it.
"""


SF_PANEL = r"""
================================================================
SF coil-design PANEL + FE-direct calc  (calc_streamfunction.py)
================================================================

The GUI panel (radia_streamfunction.py, Layer 3) wraps the headless calc
``src/radia/panels/calc_streamfunction.py`` (Layer 4).  ONE argparser drives
both, with THREE modes (--method):

Standalone launch:

  radia-streamfunction path/to/coil.vol

The Cubit ``Solve -> Radia-NGSolve`` launcher and the standalone entry point
use the same Layer-3 window and Layer-4 calc scripts.  The command-line
argument fills the panel's coil/conductor ``.vol`` browse row; each mode keeps
its own additional inputs such as ``--eval-vol``.

  design       target -> A psi = B (folded-Tikhonov RegularizedTSVD) -> psi,
               field homogeneity over the eval region, peak surface current.
  pareto       (homogeneity, peak current density) front via --pareto-lever:
                 alpha     Tikhonov L-curve (one factorisation, swept alpha)
                 linf      minimax IRLS trajectory (peak down at ~const homo)
                 geometry  eval-region (DSV) scale sweep (dual of former size)
  manufacture  psi iso-contours -> orientation-consistent equal-current turns
               -> field-aware single-stroke wire -> (sheet-metal --distort)
               -> CAD STEP (--step-output) -> PEEC L,R (--peec).
               TURN-COUNT / few-turn levers (--greedy-turns greedy constructive,
               --pin-tiling driven array, --optimize-levels): see 'low_turn'.

MATERIAL-AWARE design (iron yoke / shield / core): --iron-vol folds a Kelvin-FEM
iron reaction into design + pareto + manufacture so the whole pipeline accounts
for the iron.  See the 'material_aware' topic.

I/O: --coil-vol (a STANDALONE 2D surface .vol; psi = H1 on it, Setup-B
``definedon=Boundaries('.*')`` + ``grad(v).Trace()`` + ``ds``), --eval-vol
(surface OR volume), --target-cf (a CoefficientFunction expr of x,y,z;
scalar -> Bz, 3-vector -> full B).

CURRENT-CONFINEMENT BOUNDARY CONDITION  (--confine {off, on, abe})
-----------------------------------------------------------------
On a FINITE former the contours run off the edges; closing them with a rim
chord injects a spurious edge current (LAB short cylinder Gx: single-current
rms 0.54, 42/42 contours open).  Confine the current to the patch:

  off   no constraint (default; fine when contours close on their own, e.g.
        full-ring solenoid / long former).
  on    psi = 0 on every boundary edge (H1 dirichlet_bbnd).  Simple gradient/
        shim BC; BREAKS solenoid-type targets (the two ends are forced equal).
  abe   the CANONICAL Abe edge-equipotential BC (M. Abe, IEEE Trans. Magn.,
        DUCAS; Appendix eq.6 T = R.T_IN, A-1/A-3): each PHYSICAL boundary
        component gets ONE FREE constant + one ground.  Closes the contours on
        ANY former AND works for gradient AND solenoid (the two ends take
        different free constants).  Implemented as a DOF-reduction matrix R
        (off/on are the column-select special case); physical edges vs a CAD
        seam are told apart by element adjacency (a boundary mesh-edge borders
        ONE surface element, a seam two).

abe is the best DESIGN + SEPARATE-TURN + GENERAL choice: LAB short cylinder
Gx -> n_open 0, separate-loops single-current 0.022 (vs on's 0.149), and does
NOT break uniform (vs on which degrades it).  CAVEAT: abe is NOT automatically
best for the SINGLE-STROKE WIRE -- its edge equipotential makes a contour hug
the boundary, so the field-aware chain can connect it worse than `on`; abe's
value is the canonical BC + generality + design/separate-turn accuracy.

CONTOUR DRAWING = FLUX-LINE DRAWING (same principle)
----------------------------------------------------
The psi iso-contours are drawn by the magnetic-flux-line rule: between two
adjacent lines flows a FIXED amount (current for psi, flux for A_z) -- Abe's
"between nodes i,j flows T_i - T_j".  Equal-psi-interval contouring therefore
already gives wire density ~ |grad psi| = |K|, the same density rule the
flux-line bubble system (Hirahatake/Noguchi/Igarashi/Yamashita, bubble
r ~ 1/sqrt|B|) enforces.  Two manufacture refinements:

  --contour-sub N   order-p contour: the default marches on the VERTEX
                    (order-1) psi, dropping an order-2/3 design's edge DOFs.
                    sub=3 subdivides each triangle 3x3 and evaluates the
                    FULL-order psi via mesh.GetTrafo(el) + gfu(trafo(ip))
                    (the element-trafo MeshPoint dodges the boundary-point-
                    eval-returns-0 quirk) -- the FE analogue of the analytical
                    flux-line trace.  LAB Gx o2: loops_homo 1.32e-4 -> 1.15e-4.
  --flux-plot p.png  bubble-system flux-line view of the DESIGNED coil's B
  --flux-plane {x,y,z}   field on a cut-plane, bubble-seeded (density ~ |B|) +
                    matplotlib streamplot.  Physical check (the four-lobe Gx
                    gradient saddle renders correctly).
  --steps-plot p.png  per-step manufacturing 2x2 3D view: (1) equal-current
                    contours (N = --nlevels turns -- this sets the line count),
                    (2) single-stroke wire, (3) sheet-metal --distort wire,
                    (4) wire WITH thickness (--wire-diam, twist-free
                    parallel-transport tube) + distortion.
  --scan-map        MEASUREMENT TWIN: the 3-component (Bx,By,Bz) field of the
                    DELIVERED single-stroke wire (distorted if --distort) on a
                    planar raster -- the simulation counterpart of a 3-axis
                    magnetometer (TMR / fluxgate) scan (fixed sensor ARRAY or
                    one sensor on a moving stage sample the SAME grid).  Reports
                    Bnormal (the along-normal MAIN uniform field) AND Bt (the
                    in-plane TRANSVERSE leakage a 3-axis sensor catches -- the
                    single-stroke connector stray a 1-axis probe misses), plus
                    homogeneity; the full grid is saved to committed JSON
                    (artifact "sf_scan_map_3component"; Data Persistence).
                    --scan-plane {x,y,z} normal, --scan-standoff (default = the
                    eval/DSV height), --scan-half (default = DSV footprint),
                    --scan-n (n x n; No-Fallback: >=2), --scan-current (default =
                    the best-fit design current; set the bench current -> real
                    Tesla; field is LINEAR in current).  Compare against a
                    measured +-I current-reversal scan (coil-only, Earth + sensor
                    offset cancelled) -> the sheet-metal-PoC verification base.

PRINTABLE FORMER (3D-print the out-of-surface groove the wire threads through)
-----------------------------------------------------------------------------
The DELIVERED wire (distorted if --distort) is realised as a base plate with the
wire CHANNEL cut out (plate - swept tube), water-soluble-support printed.  Two
outputs, shared --former-clearance/-margin/-wall/-decimate knobs:
  --former-stl f.stl   RECOMMENDED.  A ROBUST MESH boolean (trimesh + manifold3d)
                    -> STL (the native 3D-print / slicer format).  Handles the 3D
                    SELF-CROSSING single-stroke wire WATERTIGHT in ~1 s (channel
                    ratio ~0.93).  Needs `pip install trimesh manifold3d` (lazy
                    import, fails loud if absent).
  --former-output f.step   OCCT (netgen.occ cylinder-Fuse) -> STEP.  OK for the
                    PLANAR demonstrator (channel ~0.96) but SILENTLY under-cuts /
                    segfaults on a 3D self-crossing wire.
CAD-KERNEL LESSON: an arbitrary-path SWEEP of a thin-profile / large-path /
self-crossing wire is fragile in the OCCT kernel ITSELF -- netgen.occ `Pipe`
segfaults, `Glue` silently no-ops, `Fuse` degenerates, and build123d's dedicated
`swept()` ALSO degenerates / raises `StdFail_NotDone` (same kernel).  So for a
robust boolean of that geometry use the MESH boolean (--former-stl); and author
CAD generally with build123d or Cubit (ACIS -- a different, robust kernel), NOT
ad-hoc netgen.occ (which is a STEP->mesh I/O library).  The full former-CAD
research note remains a private lab artifact; public radia-mcp keeps only this
general CAD-kernel lesson.

CHAIN (--chain {field_aware, nn})
---------------------------------
The "extra lines" in the single-stroke view are the inter-loop CONNECTORS
(rungs): chaining N contour loops into one wire needs N-1 bridges, and they
carry the series current, so their stray field is real (not cosmetic).
field_aware (default) keeps that field small two ways:
  (1) CUT placement -- each loop's entry/exit cut by coordinate descent to
      minimise the FULL one-current wire error min_I ||I (loops+connectors) - B||
      (NOT ||connectors|| alone -- worse on open contours; the connectors are
      routed to CANCEL the rim-chord residual, not merely to be short).
  (2) VISIT ORDER -- the same wire-error objective is minimised over a small
      candidate set {nearest-neighbour, 2-opt-shortened}, keeping whichever the
      cut-opt drives lowest.  The 2-opt shortens the long "jump to a far lobe
      and back" rungs (LAB Gx: max 372->289 mm, delivered +19..+70 %) but a
      length-optimal reorder can break the rungs' symmetric stray-cancellation
      and HURT some cases (abe nl=16: -78 %) -- the documented "shorter rungs
      != better field" trap.  Selecting the lower-wire-error order makes the
      2-opt GUARANTEED never worse than nearest-neighbour while capturing the
      real gains.  (Pure --chain nn keeps NN order + naive closing, no cut-opt,
      so it never sees the 2-opt -- which would only hurt a cut-opt-free chain.)
  (3) RESOLUTION -- the cut-opt's candidate count (--chain-ncut, cuts/loop) and
      descent rounds (--chain-passes) AUTO-scale with the loop count when not
      given (~1.1 cuts/loop capped [24,64]; passes 4..8).  More loops -> more
      connectors -> the JOINT cut-opt needs more candidates/rounds to drive
      their net field to zero.  More PASSES at a fixed n_cut monotonically lower
      the objective (coordinate descent re-scans the same grid); a higher n_cut
      samples a FINER but not strictly-nested cut grid, so auto resolution is a
      STRONG NET improvement (NOT a per-geometry guarantee -- a marginal case can
      shift a few %).  It matters most for CLUSTERED levels (--optimize-levels
      packs contours to sharpen the loop-set field): the OLD fixed default (24,4)
      under-resolved them and left a 45-loop Z2 wire ~0.20 vs a 0.0098 loop-set
      floor; auto (~50,8) recovers it to ~0.067 (3x, then --distort ~0.052).
      Override with --chain-ncut / --chain-passes.  (Bench across 6 formers: 4
      better 2-28%, 2 marginally worse 3-5% -- the cut-grid non-nesting.)
field_aware is never worse than nn.  For a SINGLE-region (l=1) psi with closed
contours (--confine abe) it reaches the separate-turns floor with no --distort.
For a MULTI-region (saddle, l>=2) psi the inter-region connectors leave a
residual gap above the loop-set floor (clustered Z2: ~0.067 delivered vs ~0.0098
loop-set) that auto-resolution SHRINKS but does not close -- splitting into one
wire per current group is WORSE (verified): the cut-opt nulls connectors GLOBALLY
by cross-cancelling the +/- groups' rungs, a freedom per-group chaining removes.

END-TO-END VALIDATION vs an INDEPENDENT codebase
------------------------------------------------
validation_test/stream_function/verify_coil_field_independent.py designs a coil
(MRI-gradient scale: cylinder r=0.15 m, L=0.5 m, DSV r=0.05 m) and checks the
field TWO ways: the numpy straight-segment Biot-Savart used in the designer AND
Radia's C++ rad.ObjFlmCur + rad.Fld (a separate codebase).  They agree to
8-11 digits (uniform 3.5e-11, Gx 1.1e-8); the abe-confined Gx coil reaches
1.0 % nonlinearity on the short former, cross-validated.  Locked by
validation_test/panels/test_notebook_workbench.py plus the streamfunction
golden validation lane.
"""


SF_FUSION = r"""
================================================================
SF for FUSION -- stellarator Stage-2 (NESCOIL / REGCOIL / FOCUS)
================================================================

The SAME surface stream function (current potential psi on a conductor
surface, K = n x grad psi) IS the object the stellarator-coil community calls
the WINDING-SURFACE CURRENT POTENTIAL -- the unknown of NESCOIL (Merkel 1987),
REGCOIL (Landreman 2017), and the surface part of FOCUS (Zhu 2018).  The
"Stage-2" coil problem is mathematically identical to MRI-gradient / IH coil
design:

    given a target NORMAL field B.n on the PLASMA boundary,
    find psi on a WINDING SURFACE around it whose Biot-Savart n.B reproduces it
    (iso-contours of psi = the coils)

The design-matrix rows are just the plasma-normal component n.B of the
winding-surface Biot-Savart kernel we already assemble (A3 = 3-component A at
the plasma points; A_n = einsum('mc,mcj->mj', plasma_normal, A3)).  No
fusion-specific solver code.

DOCS / VALIDATION SURFACES
--------------------------
  docs/stream_function/fusion.md
  docs/stream_function/demo_gallery.ipynb
  validation_test/stream_function/regcoil_fusion_helpers.py

demo_regcoil_fusion.py
  1. FORWARD MAP IS EXACT.  Producible targets -- uniform vertical (PF /
     equilibrium coil) and non-axisym sin(th)cos(2ph) (stellarator-like) --
     reproduce to B.n residual ~2e-9 (machine precision).
  2. REGCOIL L-CURVE.  On a hard, NOT-cheaply-producible target sin(3th)cos(5ph)
     (decays across the plasma-coil gap), sweeping alpha traces the classic
     (B.n residual, peak|grad psi|) trade-off; knee at alpha_rel ~ 2e-2, peak
     current density saturates at the surface representation limit.
  3. NET CURRENT = the multivalued / SECULAR term.  A single-valued psi carries
     ZERO net current through each torus hole; the full current potential is
     Psi = psi + (G/2pi) zeta + (I/2pi) theta.  The TWO extra DOFs are the first
     cohomology generators of the winding surface; their COUNT is the surface
     Betti number b1 = 2 (genus 1), CONFIRMED gmsh-free via the Euler
     characteristic of the triangulated torus (b1 = 2 - chi = 2g).  The
     volume T-Omega cohomology CUT engine src/radia/cohomology_cut.py is
     itself gmsh-free (radia.cohomology); on the torus the surface
     generators are analytic: grad(zeta), grad(theta) are single-valued.
     K_zeta = n x grad(zeta) = the net-POLOIDAL-current (TF)
     sheet; VERIFIED to give the textbook toroidal field B_tor*R = const inside
     the tube (Ampere 1/R, to 0.2 %) and ~0 outside.  KEY PHYSICS: the TF field
     is TANGENT to the plasma, so its B.n footprint is >1000x smaller than the
     net-toroidal generator's -> the net poloidal current is NOT fitted from
     B.n, it is a PRESCRIBED engineering parameter (1 T on axis @ R=0.3 -> 1.5
     MA).  This is exactly why REGCOIL takes net_poloidal_current as an INPUT.
  4. VMEC-SHAPED BOUNDARY.  A non-axisym rotating-ellipse boundary in the VMEC
     Fourier form R = sum RBC cos(m th - n NFP ph), Z = sum ZBS sin(...), with
     analytic normals; --wout wout_*.nc reads a real equilibrium boundary
     (netCDF4 rmnc/zmns/xm/xn/nfp, stellarator-symmetric only -- raises on
     lasym=T).

demo_regcoil_fusion_advanced.py
  A. COIL FORCE / STRESS.  Lorentz force per area f = K x B_avg (B_avg = the
     +/-eps average of the coil self-field across the sheet).  For the TF coil
     |f| = magnetic pressure B_tor^2/(2 mu0) (ratio ~0.99) and CONCENTRATES on
     the INBOARD leg (~5x outboard) -- why a TF coil is inboard-stress-limited.
     Honest: the magnetic force per area (stress DRIVER, N/m^2), not a
     structural hoop-stress model.
  B. REAL EQUILIBRIUM.  Designs against the li383 (NCSX-like, NFP=3, QA)
     reference wout from simsopt (MIT); --wout for any VMEC output, else fetches
     li383 (121 kB).  B.n on the genuine 25-mode boundary to ~4e-8.
  C. FOCUS STANDOFF.  Coil complexity peak|grad psi| is MONOTONIC in the winding
     gap (closer = simpler), so the DISTANCE optimum is CONSTRAINT-BOUND (push to
     the minimum engineering standoff d_min).
  D. FOCUS WINDING-SHAPE (the core FOCUS contribution).  _surface_mesh_from_grid
     builds an NGSolve surface mesh from an ARBITRARY (theta,phi) point grid
     (manual netgen Element2D + FaceDescriptor; SF design machine-precision,
     matches the OCC torus), so the winding can be CONFORMAL to a shaped plasma.
     For an elongated (kappa=2) plasma, blending the winding circular (varying
     gap) -> conformal (uniform gap) at the same min standoff cuts coil
     complexity by ~34 %.  The winding SHAPE, not just distance, is a real lever.

HONEST SCOPE
------------
  * Parts 1-2 single-valued psi (PF/RMP/shaping); part 3 the multivalued secular
    term (generator COUNT computed, generators analytic ON THE TORUS -- a general
    winding surface uses cohomology_cut.py).
  * Force is magnetic force/area, not structural stress.
  * The winding-shape study sweeps a 1-parameter circular->conformal blend; a
    full Fourier-mode winding-surface optimiser is the named next step.
  * Real VMEC: the reader is round-trip-verified against the wout schema; we do
    NOT run VMEC (no simsopt/desc here) -- the default rotating-ellipse is an
    analytic MODEL, --wout drops in a real equilibrium.

Locked by tests/panels/test_streamfunction_golden.py
(test_regcoil_fusion_* : forward machine precision, b1==2, TF Ampere 1/R +
tangency, VMEC non-axisym, wout round-trip + lasym-raise, force=pressure,
FOCUS monotonic + conformal<circular).

radia-SF vs NESCOIL / REGCOIL / FOCUS -- A SUPERSET (2026-06-21)
---------------------------------------------------------------
radia ships the SAME object (winding-surface current potential psi,
K = n x grad psi) but is a strict SUPERSET on deliverable + physics.
Precise claim: "REGCOIL is a special case (vacuum, toroidal, design-only)
of radia's SF; design AT PARITY, deliverable + iron a SUPERSET; a complete
win for single-conductor coils."

  DESIGN -- AT PARITY.  The demos above reproduce the REGCOIL workflow
    (current potential, Tikhonov L-curve, net-current secular term, VMEC
    boundary, force/stress, FOCUS standoff + winding-shape); B.n ~2e-9 on
    producible targets.  No fusion-specific solver code.
  SCALE / METHOD -- ACA + ridge-TSVD.  The design matrix is the dense
    Biot-Savart coupling COMPRESSED by ACA (an H-matrix) and inverted by a
    ridge (Tikhonov) TRUNCATED-SVD pseudo-inverse (radia.stream_function
    aca_tsvd; topics 'regularized' / 'performance').  The ridge IS REGCOIL's
    lambda, but ACA + the FE-direct surface make it a SCALABLE regularised
    least-squares on an ARBITRARY meshed winding surface -- not REGCOIL's
    dense Fourier least-squares on a parameterised torus.
  OPTIMISE -- fast TSVD -> Pareto + surface deformation.  The TSVD core is
    FOLDED ONCE, so each regularised re-solve is a ~50 us k x k core solve
    (RegularizedTSVD; demo_pareto_tikhonov_aca.py).  That cheap re-solve makes
    a whole front affordable, so radia sweeps a MULTI-OBJECTIVE PARETO over
    BOTH the regularisation alpha (misfit-vs-energy + an L-inf IRLS
    misfit-vs-PEAK-current front) AND the winding-surface SHAPE
    (geometry/geom_scale lever + sheet-metal deform) -- richer than REGCOIL's
    single-alpha L-curve.  (topic 'pareto'.)
  DELIVERABLE -- design-to-manufacture, not design-only.  NESCOIL/REGCOIL/
    FOCUS stop at psi.  radia continues: contours -> SINGLE-STROKE wire
    (grad-psi winding orientation, so l>=2 saddle shims chain without the
    common series current cancelling) -> sheet-metal distort -> STEP CAD
    (OCC WriteStep) -> PEEC.  The PEEC step is a full circuit-extraction
    SOLVER (L,R,C,M + SPICE, MMM coupling), not just an inductance number.
    (topics 'single_stroke' / 'low_turn'.)
  PHYSICS -- iron.  NESCOIL/REGCOIL/FOCUS are free-space (vacuum Biot-Savart).
    radia's material-aware kernel (Kelvin-FEM DtN, M = M_free + M_react,
    --iron-vol) designs coils WITH iron yoke/shield/core -- domains REGCOIL
    structurally cannot enter.  (topic 'material_aware'.)
  HONEST NUANCE.  The single-stroke win is decisive for SINGLE-CONDUCTOR coils
    (MRI gradient/shim, IH).  Stellarator MODULAR coils are intentionally
    SEPARATE coils, so "one wire" is not their manufacturing step; there radia
    still adds STEP CAD + PEEC L beyond REGCOIL's psi, and its distort is
    analogous to FOCUS filament-shape opt.  Open gaps: no SIMSOPT/STELLOPT
    Stage-1+2 integration; head-to-head at community mode counts unmeasured.

EARNED BY MEASUREMENT (the three goldens, Repository-First)
----------------------------------------------------------
  (a) VACUUM PARITY + DELIVERABLE -- test_regcoil_parity_deliverable_golden.py:
      a uniform-vertical B.n design on the torus winding surface hits
      B.n_rel ~4.9e-9
      (parity) AND the SAME run emits a ~954 kB STEP + PEEC L ~3.09 uH
      ("design equal, deliverable beyond").
  (b) IRON DIFFERENTIATOR -- test_regcoil_iron_differentiator_golden.py
      (act8_03_general_iron_design SCENARIO B): on a target behind iron the
      FREE-SPACE design MISSES by >20% while the MATERIAL-AWARE design HITS to
      <1e-2 (ratio >10x) -- the physics REGCOIL cannot do.
  (c) MANUFACTURE END-TO-END -- test_streamfunction_manufacture_e2e_golden.py:
      ONE --method manufacture run takes target x -> psi (design ~0.042%) ->
      single-stroke wire (~0.197%) -> sheet-metal distort (~0.175%, no regress)
      -> STEP CAD -> PEEC L ~81.8 uH, on the cylinder + DSV fixture.

REFERENCES
----------
  P. Merkel, Nucl. Fusion 27, 867 (1987)               -- NESCOIL
  M. Landreman, Nucl. Fusion 57, 046003 (2017)         -- REGCOIL
  C. Zhu et al., Nucl. Fusion 58, 016008 (2018)        -- FOCUS
  docs/stream_function/fusion.md
"""


SF_HARMONICS = r"""
SPHERICAL-HARMONIC TARGET + FIELD DECOMPOSITION (MRI gradient / shim basis)
==========================================================================
In a current-free region Bz is harmonic, so it expands in REAL REGULAR SOLID
HARMONICS R_l^m(x,y,z) -- homogeneous harmonic polynomials (Laplace nabla^2
R = 0).  These ARE the named MRI gradients/shims; the SF designer speaks them
natively for BOTH the target and the analysis of the achieved field.

NAMED BASIS (l <= 4 -> 1+3+5+7+9 = 25 harmonics; one poly-string table =
the single source of truth):
  l=0  Z0 = 1
  l=1  Z = z (Gz),  X = x (Gx),  Y = y (Gy)
  l=2  Z2 = z^2 - (x^2+y^2)/2,  ZX = xz,  ZY = yz,  C2 = x^2-y^2,  S2 = xy
  l=3  Z3 = z^3 - 1.5 z(x^2+y^2),  Z2X, Z2Y, ZC2, ZS2, C3, S3
  l=4  Z4, Z3X, Z3Y, Z2C2, Z2S2, ZC3, ZS3, C4 (= x^4-6x^2y^2+y^4), S4
  m>0 = cos(m.phi) (C);  m<0 = sin(|m|.phi) (S).

TARGET -- ``--target-harmonic`` (alternative to ``--target-cf``):
  a name, ``l=L,m=M``, or ``(L,M)``, optionally weighted and summed:
    --target-harmonic X            Gx gradient
    --target-harmonic Z2           pure 2nd-order shim
    --target-harmonic Z4           4th-order zonal shim
    --target-harmonic "Z2:1.0,Z:0.1"   Z2 with a Z offset
    --target-harmonic "l=4,m=-4"   == S4  (the comma inside l=L,m=M /
                                    (L,M) is rejoined, not split into terms)
  It generates the solid-harmonic ``--target-cf`` polynomial (so the whole
  pipeline -- design / pareto / manufacture / single-stroke -- is unchanged).
  Give ``--target-cf`` OR ``--target-harmonic`` (loud error on both).

ANALYSIS -- the achieved Bz over the DSV is decomposed in DESIGN mode
(``result["harmonics"]``, depth ``--harmonic-lmax`` default 3, max 4):
  * spectrum   per-(l,m) field RMS over the DSV (the harmonic content, in T),
               sorted; each with its LSQ coefficient and field fraction
  * residual_fraction   how completely harmonics up to lmax capture the field
  * purity     (with a --target-harmonic) the target-harmonic field fraction
               -- the standard "gradient purity" gradient-coil quality metric
  * max_contaminant   the largest non-target harmonic (e.g. a Gx coil's Z2X)

VERIFIED (cylinder fixture, order 2, confine abe; gallery record
demo_shim_coil_purity.py): --target-harmonic X ->
dominant X, purity 1.000, residual 1.5e-4, Z2X contaminant 7e-5; Z2 ->
purity 1.000; the 4th-order Z4 shim -> purity 0.99983, residual 1.5e-2,
named Z3 contaminant 9.5e-3 (high-l shims are harder: an l-th harmonic's
field scales as r^l over a fixed DSV).  The named-basis harmonicity
(Laplacian 0, all 25 entries) and the target<->decompose round-trip are
golden-locked (tests/panels/test_streamfunction_golden.py
test_harmonic_basis_is_harmonic / test_harmonic_l4_forms_and_decompose).
The panel auto-generates the two flags (cli-diff clean).
"""


SF_SHIELDING = r"""
ACTIVE SHIELDING -- primary + shield coil (Turner 1986 / Mansfield-Chapman)
==========================================================================
An actively-shielded gradient coil adds a second, OUTER cylindrical surface
(the SHIELD) whose current cancels the STRAY field outside the assembly --
essential in MRI so the switching gradients do not induce eddy currents in
the cryostat.

CLI (calc_streamfunction.py, design mode only):
  --shield-vol SHIELD.vol        outer shield coil surface mesh
  --shield-eval-vol EXT.vol      external region where stray field -> 0
                                 (a LARGE annular shell OUTSIDE the shield)
  --shield-weight w  (def 1.0)   weight on stray-null rows vs DSV target rows

METHOD: primary + shield designed JOINTLY in one stacked LS system
  A_stack = [ A_dsv_p  | A_dsv_s  ]   (DSV target rows)
            [ w*A_ext_p| w*A_ext_s ]  (stray-null rows)
  B_stack = [B_target; 0]
  block-diagonal seminorm diag(S_primary, S_shield); RegularizedTSVD
  unchanged.  psi_f = [primary_free | shield_free] is split after solving.
  (_build_shielded_problem in calc_streamfunction.py.)

REPORTED (shielded run JSON):
  homogeneity_rms   DSV field error (target preserved)
  stray_rms         HONEST combined field at INDEPENDENT external measure
                    points (e_mpts = external-mesh centroids), relative to
                    the DSV target -- the generalising shielding metric
  stray_fit_rms     circular fit residual at the constraint VERTICES (info
                    only -- NOT the shielding quality)
  peak_J_primary / peak_J_shield

COVERAGE IS EVERYTHING (measured, the key lesson): the external null region
must COVER the exterior, not a thin mid-plane slice.
  * LARGE full-length shell (r=0.235-0.27, spans the coil) -> GENUINE
    ~86x (39 dB) stray reduction, 100-1700x in the far field, DSV
    homogeneity preserved (gallery record demo_active_shield.py).
  * TOO-SMALL thin slice -> point-sampled nulling OVERFITS the slice (~4x
    locally) and makes the field WORSE at larger z.
  * the un-sampled GAP between shield (r=0.20) and null region (r=0.235)
    stays unshielded (~2-3x).
An analytic external-multipole-moment constraint (vs point sampling) would
remove the sampling dependence -- documented next step, not done.

VERIFIED: test_streamfunction_active_shielding (golden, large-shell fixture
make_shielded_coil_vol.py).  Panel auto-emits the 3 flags (cli-diff clean).

HONESTY NOTE: the first implementation reported stray at the CONSTRAINT
points (circular -> ~machine-zero "-184 dB"); verify-first caught it and the
headline is now the independent-point measure.  Repository-first: a real
~86x beats a perfect-looking circular artifact.
"""


SF_CLEBSCH_3D = r"""
================================================================
3D STREAM FUNCTION + COHOMOLOGY (Clebsch volume extension)  [IN PROGRESS]
================================================================
STATUS 2026-06-11: verified-conceptual + math; implementation Phase 1-3 gated on
golden tests (two derivation/verification workflows in flight).  Honestly marked:
the GENERAL non-axisymmetric 3D inverse designer is a non-convex + helicity-
obstructed FRONTIER, NOT a shipped capability.  Repository-first: a clearly
marked frontier is a real result, not a failure.

WHY (CEFC 2026 WA-O1, Thessaloniki -- corrected attribution)
------------------------------------------------------------
  * Talk #1 (presented by Z. Ren as proxy): "Dual Formulations of Multi-Port EM
    Field Problem Using Cohomology Modeling", Zhou/Yan/Xu/Ren (IEE-CAS + GeePs
    Sorbonne).  de Rham cohomology for non-local port excitations; thick cuts /
    thick links; H1<->MMF, H2<->flux; dual h/b formulations + complementary
    energy bounds.  = "Ren's cohomology".
  * Talk #2 (Tampere -- the interesting one, NOT Ren, NOT cohomology):
    "Bidirectional Coordinate Transformation ... 2-D Magnetic Field Problems",
    Dervisha/Marjamaki/Rasilo/Tarhasaari.  Exterior-calculus "potential
    coordinates" (A, phi): A = vector-potential z-component = the 2D stream/flux
    function; bidirectional forward (x,y)->(A,phi) and inverse (A,phi)->(x,y)
    (redraw geometry in potential coordinates).
  * Cohomology-cut citation stack: Kotiuga (theory) -> Pellikka 2013 (the Gmsh
    algorithm cohomology_cut.py calls) -> Ren 2002 / talk#1 (T-Omega dual, 2ndary).

THE UNIFICATION (the payoff)
----------------------------
Everything here lives on ONE NGSolve de Rham / FEEC complex:
    H1 --grad--> H(curl) --curl--> H(div) --div--> L2
  * H1 (0-forms): stream/Clebsch scalars psi, chi; scalar potential phi; 2D flux
    function A.  (The single-stroke deformation field is VectorH1.)
  * H(curl) (1-forms): cohomology generators h_k; grad of the multivalued chi;
    H field.  H1 cohomology = curl-free-but-not-gradient (harmonic) subspace.
  * H(div) (2-forms): B, and J = grad(psi) x grad(chi) = grad(lambda) x grad(mu).
    H2 cohomology = div-free-but-not-curl subspace.
Consequence: Clebsch, cohomology, the stream function and the single-stroke
deformation are ALL native NGSolve/Netgen FEEC objects -> the stack can drop
Gmsh.  Cohomology generators are computable natively as the Hodge-Laplacian
zero-eigenspace on H(curl) (replaces Gmsh computeHomology; reuses the Compact HX
H(curl) regular decomposition).

FIELD-CLEBSCH vs CURRENT-CLEBSCH (do not confuse)
-------------------------------------------------
  * FIELD-Clebsch  B = grad(psi) x grad(chi) -- the field rep
    (clebsch_potential.py ClebschSolver, forward analysis).  BILINEAR in
    (psi,chi).
  * CURRENT-Clebsch J = grad(lambda) x grad(mu) -- the design bridge.
    div J = 0 identically.  wires = {lambda=const} INT {mu=const}; per-tube
    current dI = Dlambda * Dmu (signed flux = Jacobian; valid where the
    (lambda,mu) chart is non-degenerate, i.e. grad lambda not || grad mu).

ACA+TSVD ON CLEBSCH (the key result)
------------------------------------
ACA+TSVD solves LINEAR least-norm A x = b.
  * FIELD-Clebsch is bilinear -> ACA+TSVD does NOT directly apply to inverse
    design.
  * CURRENT-Clebsch with the foliation mu FIXED is LINEAR in lambda:
    lambda = sum_j lambda_j phi_j -> J = sum_j lambda_j (grad(phi_j) x grad(mu)),
    so B = A lambda is linear and the EXISTING kernel-agnostic ACA+TSVD engine
    (radia.stream_function aca_tsvd / RegularizedTSVD) applies UNCHANGED -- only
    the per-DOF basis current changes from the surface n x grad(phi_j) to the
    volume grad(phi_j) x grad(mu).
  * The surface SF is the SINGLE-SURFACE special case (mu = one surface); the
    volume version = a STACK / foliation = a true 3D bulk winding.  The ONLY
    nonconvex part is choosing the foliation mu (outer loop); the surface method
    dodges it by mu = the winding surface.  ACA+'s relative win is LARGER in 3D
    (volume Biot-Savart kernel is more expensive; ACA+ speed is kernel-cost-
    driven, measured 1.6x->10.8x as N grows for the surface case).

CLEBSCH <-> COHOMOLOGY (the connection)
---------------------------------------
  * axisym Clebsch (psi = r*A_phi, chi = phi) IS a stream-function rep; chi =
    atan2(y,x) has grad = (-y/r^2, x/r^2, 0) = a de Rham 1-cohomology GENERATOR
    (= 2*pi * the Gmsh/Pellikka h_1: SAME class, DIFFERENT representative + 2*pi
    normalization, NOT literally equal).  The atan2 branch cut IS the cohomology
    cut -- a topological feature, not a numerical bug.
  * Net current => one Clebsch scalar MUST be multivalued = a cohomology
    generator (REGCOIL-style: prescribe net current as the secular term, fit the
    single-valued remainder).  See the 'fusion' topic for the surface b1=2 case
    already shipped.
  * TWO DISTINCT global obstructions, do NOT conflate: (a) domain topology b1
    (net current); (b) field helicity / Hopf (linked field lines, no global
    single-valued Clebsch under B.n=0).  Axisym is helicity-zero, so the repo has
    only hit (a).

BIDIRECTIONAL POTENTIAL-COORDINATE MAP (Tampere -> 3D)  [VERIFIED sympy 2026-06-11]
----------------------------------------------------------------------------------
  * 2D inverse-map PDE (Tampere): d/dA(mu dx/dA) + d/dphi((1/mu) dx/dphi) = 0
    (and for y).  Forward = geometry->potentials; inverse = potentials->geometry
    -- the design lever ("design in flux coordinates, manufacture in physical").
  * 2D Jacobian: d(x,y)/d(A,phi) = +-1/(|B||H|) (Lagrange identity, VERIFIED).
  * 3D flux coordinates (psi, chi, phi): SEMI-orthogonal -- grad(phi) perp
    grad(psi) AND grad(phi) perp grad(chi) (from grad psi x grad chi = -mu grad
    phi); (psi,chi) skew within the flux surface.  (VERIFIED)
  * 3D Jacobian (VERIFIED sympy): det[grad psi; grad chi; grad phi]
    = -(1/mu)|B|^2 = -|B||H| -> inverse = -1/(|B||H|), UNIFYING 2D and 3D under
    |B||H| != 0; degeneracy locus = field nulls.

VECTOR-T CONVEX INVERSE + HELICITY OBSTRUCTION (VERIFIED sympy 2026-06-11)
-------------------------------------------------------------------------
  * "Problem with T = lambda*grad(mu)": t=lambda d mu has t^dt=0 (Frobenius),
    i.e. T.curl T = 0 pointwise (VERIFIED) -> CURRENT helicity H = int T.curl T = 0
    for any Clebsch current; with J.n=0 on dOmega it is gauge-invariant.
    CONTRAPOSITIVE: a NONZERO-helicity current (linked/knotted lines, force-free-
    like) has NO global Clebsch J=grad(lambda)x grad(mu).  Ordinary coils (nested
    loops / solenoid) are helicity ~ 0 so Clebsch is fine; the obstruction bites
    only for topologically nontrivial (linked) current designs.
  * REFRAME: non-convexity is SPECIFIC to the Clebsch (bilinear) form, NOT to
    volume inverse design.  With the FULL VECTOR T (HCurl), J=curl T is LINEAR in
    T -> A T = B_target is a CONVEX least-norm the existing ACA+TSVD solves
    UNCHANGED; the gauge null space (T=grad chi -> curl=0, VERIFIED) is in ker(A)
    and TSVD truncates it automatically (NO tree-cotree).  Helicity does NOT
    obstruct the vector-T form.
  * TRADE-OFF: vector-T = convex / no-direct-wires / no-helicity-limit;
    Clebsch mu-free = nonconvex / wires / zero-helicity; Clebsch mu-fixed
    (foliated) = convex / wires / zero-helicity.  The REAL frontier is (i) WIRE
    EXTRACTION from a general T/J, (ii) the Clebsch FOLIATION choice -- NOT
    "non-convex helicity".  Stage A (buildable) = vector-T convex inverse golden;
    Stage B (frontier) = wire extraction / foliation.

SINGLE-STROKE ERROR CORRECTION via DEFORM (already shipped; VectorH1 upgrade)
----------------------------------------------------------------------------
  * _distort_wire (calc_streamfunction.py; inner closure literally deform())
    bends the ONE single-stroke wire (psi/levels fixed, one current) to cancel
    the connector residual.  Framed as "control-grid bilinear = the discrete
    realisation of an NGSolve VectorH1 mesh deformation".  Measured: plane
    12290->340 ppm; cylinder Gx 8.5-9.3%->1.4%; sphere Z2 4.3%->0.36%.  See
    'single_stroke'.
  * Clean upgrade: replace the control-grid bilinear with a proper VectorH1
    deformation field on the FE-direct surface mesh (H1-regularised bends,
    manufacturability via BC, arbitrary surfaces).  VectorH1 is the deformation
    FIELD REPRESENTATION + shape derivative; the field eval stays Biot-Savart
    (NOT a SetDeformation-then-PDE-solve loop).

IMPLEMENTATION STATUS (tests/ -> examples/ -> panels/)  [Stage A DONE 2026-06-11]
--------------------------------------------------------------------------------
  Phase 0  record the unification (memory/clebsch_cohomology_streamfunction_
           unification.md) -- DONE.
  Phase 1  2D bidirectional map (Tampere) -- DONE: validation_test/feec/vim_legacy/
           bidirectional_map_2d.py + tests/feec/test_bidirectional_map_2d.py (5).
           w=zeta^2: forward A,phi 1e-15; orthogonality 2.6e-14; Jacobian
           det=|B||H| 1.6e-16; inverse x,y harmonic in (A,phi) 1.6e-8; round-trip
           6e-8.  Pure NGSolve, no Gmsh.
  #2       vector-T convex inverse (general form) -- DONE: validation_test/feec/vim_legacy/
           vector_t_inverse.py + tests/feec/test_vector_t_inverse.py (5).
           J=curl T (HCurl), ACA+TSVD min-norm: convex fit 6e-16; div J=5.7e-16;
           GAUGE T=grad(chi)->field 2.2e-20, TSVD truncates (k_aca=9<=M<<ndof=3360).
  Phase 3  foliated-Clebsch ACA+TSVD solenoid -- DONE:
           validation_test/feec/vim_legacy/foliated_clebsch_solenoid.py +
           validation_test/feec/test_foliated_clebsch_solenoid.py
           (4).  mu=r fixed -> J linear in lambda -> radia.stream_function
           aca_tsvd UNCHANGED; uniform Bz to 3.3e-5, weak div J 8.5e-6.
  Phase 2  axisym 3D probe (psi,phi) flux coords in (r,z) -- the 2D map covers
           the meridian case; not separately built.
  STAGE B  WIRE EXTRACTION (continuous volume current -> windable wires) -- the
           BUILDABLE part DONE:
    B1     helicity diagnostic -- validation_test/feec/vim_legacy/helicity_diagnostic.py +
           validation_test/feec/test_helicity_diagnostic.py (4).  H_rel=|int T.curlT|/
           (||T|| ||curlT||) in [0,1]: axial 9e-5, ABC/Beltrami 1.000.  Gates
           whether a clean Clebsch/level-set extraction is even possible.
    B2     wire extractor -- validation_test/feec/vim_legacy/foliated_solenoid_wires.py +
           validation_test/feec/test_foliated_solenoid_wires.py (5).  Foliated solenoid
           lambda -> per-cylinder equal-Delta-lambda contours -> 59 equal-current
           wires (I=dlam*dmu) -> Biot-Savart reproduces uniform Bz to 4.7%, AGREES
           WITH RADIA rad.ObjFlmCur+rad.Fld to 3.4e-10 (two-codebase invariant).
  FRONTIER F1/F2/F3 -- DEMONSTRATED as fundamental walls (golden-locked, honest
           partly-negative results; NOT fake-completed, NOT fake-failed):
    F1     foliation_choice_wall.py + test (4): foliation choice = bilevel
           CONVEX-inner/NON-convex-outer.  Field fit ~1e-15 at every tilt, so the
           non-convexity is in COIL COMPLEXITY (||J||, REGCOIL measure), invisible
           to the field residual.  g(t): single-axis 1 min (convex), two-axis 2
           mins + 18% barrier.  = stellarator winding-surface choice (Merkel/
           Landreman/Kaptanoglu arXiv:2408.08267).  Remedy = global opt or vector-T.
    F2     clebsch_recovery_wall.py + test (5): recovery of (lambda,mu) from J is
           ILL-POSED.  WALL1 non-unique (transport J.grad(mu)=0: mu0=x vs x^2 both
           valid, mu_b=mu_a^2 -> infinite-dim area-preserving-diffeo gauge); WALL2
           cond ~1/|J|^2 slope -4.00 at field nulls (Yoshida 2009, Qin 2018).
    F3     helical_current_no_clebsch.py + test (6): nonzero helicity = topological
           obstruction.  ABC target H_rel~1; convex Clebsch fit FORCED to H_rel~0
           -> helicity GAP ~1.000 (foliation+mesh-independent; the L2 residual is
           the WRONG metric).  Clebsch streamline closes (2e-3), ABC does not (143
           ring-lengths, ergodic).  Moffatt 1969, Enciso-Peralta-Salas 2020.  Only
           honest outputs: vector-T bulk distribution OR multi-patch atlas + cuts.
  Session golden total: Stage A 14 + Stage B 9 + frontier 15 = 38, all green.
  PANELIZED (examples -> panels, 2026-06-11): Stage A/B shipped as the
           radia_streamfunction "Volume 3D" mode.  shipped pipeline
           src/radia/streamfunction_volume.py (design_volume_coil); headless
           calc src/radia/panels/calc_streamfunction_volume.py (conductor .vol +
           --target-bz -> equal-current wires + GMSH wire overlay); PySide6
           _Volume3DPanel in radia_streamfunction.py; golden
           tests/panels/test_streamfunction_volume_golden.py reproduces n_wires=59,
           field 4.7%, two-codebase 3.4e-10 on a frozen tube fixture.  Panel
           covers the CLEAN regime only; F1/F2/F3 stay research demos by design.

Files: src/radia/clebsch_potential.py (ClebschSolver, AxisymStreamFunctionSolver),
src/radia/cohomology_cut.py (gmsh-free cohomology cut via radia.cohomology), src/radia/
stream_function.py + panels/calc_streamfunction.py (ACA+TSVD engine),
memory/clebsch_cohomology_streamfunction_unification.md (full record).
"""


SF_MATERIAL_AWARE = r"""
================================================================
SF coil design WITH MAGNETIC MATERIAL (iron yoke / shield / core)
================================================================
The free-space Biot-Savart kernel BREAKS the moment iron is in the picture: a
coil designed in vacuum produces the wrong field once a yoke / shield / core
re-routes the flux.  ``--iron-vol`` makes the ENTIRE existing pipeline
(design / pareto / manufacture / single-stroke / distort) MATERIAL-AWARE, so the
designed psi already accounts for the iron reaction.  (This was the Track-B
"SF with iron" handoff goal; now shipped in the production solver.)

CLI (calc_streamfunction.py):
  --iron-vol IRON.vol     a Kelvin open-boundary FE mesh of the iron region
                          (materials: the iron block + air + a `kelvin` shell;
                          Periodic kelvin_int<->kelvin_ext + a GND vertex --
                          the same .vol convention as the Kelvin DtN demos)
  --mu-r MU               iron relative permeability (linear)
  --iron-mat NAME         iron material label in the .vol (default "iron")
  --iron-exact-source     [OPT-IN] evaluate the basis-current source at the iron
                          QUADRATURE points (FE-direct, no P1-vertex interp)
  --iron-quad-order N     quadrature order for --iron-exact-source

METHOD (reduced scalar potential, mu0 cancels by linearity)
-----------------------------------------------------------
  H = H_s - grad(Omega)
    H_s   = Biot-Savart of the SF basis current K_j = n x grad(phi_j)
    Omega = the IRON REACTION, a Kelvin open-boundary FE solve driven by H_s
  The design transfer splits  M = M_free + M_react ; design = invert M.  The
  iron reaction is assembled as A_react and ADDED to the free-space design matrix
  (Ac/Am += reaction) behind --iron-vol, so nothing downstream changes.

VERIFIED (test_streamfunction_iron_golden.py + the DtN act8_10_current_sheet_design_matrix it was ported
from): designing WITH M HITS the target (isolation 1.3e-15) while the free-space
M MISSES by ~12-13%; iron vs free-space peak_J differs ~3.8% on the full CLI.
iron_changes_the_design / iron_aware_design / iron_exact_source /
iron_missing_material_raises are golden-locked.  (commits: iron 3f16aff5,
obs-adjoint 7baa5ae9, exact-source 7603e053.)

SCALABILITY -- obs-adjoint (commit 7baa5ae9)
--------------------------------------------
The reaction is built with ONE back-substitution per OBS ROW (not one per
psi-DOF): A_react = einsum(E_iron, B_iron), E_iron = D^T Ainv G^T read on the
iron source DOFs, G[row] = the -grad(.)(obs)[comp] covector (element-local basis
grads via mesh(*p).nr -> ElementId(VOL,nr)).  Verified BIT-IDENTICAL to the
forward per-DOF loop (1.1e-14); the cap rose 3000 -> 20000 DOF (now B_iron-memory
bound, not solve-bound).  --iron-exact-source uses the adjoint identity
A_react[i,j] = (mu_r-1) integral_iron B_s^j . grad(X_i) with B_s^j at the
quadrature points; exact vs P1 differ ~14% on a coarse shell (the removed
source-interpolation plateau).  Default stays the fast P1 (opt-in exact).

WHY NOT FMM FOR THE DESIGN MATRIX (recurring question)
------------------------------------------------------
FMM (ngsolve.bem.BiotSavartCF) takes FILAMENT segments.  The SF DESIGN matrix
needs the field PER psi-DOF, and each basis current n x grad(phi_j) has NO clean
filament representation, so FMM CANNOT assemble the design matrix.  FMM only
fits a KNOWN / finished psi (filament-izable into contour loops) for apply-once
field eval (viz / verify / a fixed source).  The SF design matrix is FE-direct
Biot-Savart (O(N^2)) + ACA+TSVD (HACApK) compression -- this is exactly why the
accuracy upgrade became exact-quad-adjoint, NOT FMM.

No-Fallback: --iron-vol and --shield-vol together RAISE (an active shield coil
is not iron); free-space is the default when neither is given.

See also: 'shielding' (active shield COIL, a different lever), the maintained
Kelvin/open-boundary routing in docs/kelvin/ARCHIVE_RETIREMENT.md, and
[[sf-coil-model-current-sheet-vs-dirichlet]].
"""


SF_LOW_TURN = r"""
================================================================
LOW-TURN / MANY-TURN MANUFACTURE -- discrete refinement of psi
================================================================
CONCEPTUAL FRAME (user-approved): the stream function is the CONVEX CONTINUOUS
relaxation of the coil; a finite-turn coil is that relaxation PLUS discrete
ROUNDING / REFINEMENT (the topology-opt SIMP analogy); single-pass (one wire)
is the N=1 extreme (path optimisation).  So "few turns" is not a separate
paradigm -- it is the discrete back-end that completes SF for ALL turn counts.

CONTOUR ORIENTATION = grad-psi winding (the single-wire l>=2 fix, production)
---------------------------------------------------------------------------
The single-wire / equal-current manufacture path orients each contour loop by
the surface-current winding K = n_hat x grad_s(psi) (a KDTree majority vote of
segment-direction . K per loop; global sense aligned to the target), NOT by the
old per-loop field alignment sign(f.B).  WHY it matters: the marching-triangles
polyline order is arbitrary, so its winding does not carry the current sign.
sign(f.B) ('flip every loop to help the target') recovers it ONLY for a MONOTONE
psi (l=1 gradient); for a SADDLE psi (l>=2 shim: Z2, Z3, C2/S2 ...) the natural
winding has MIXED senses (the single wire reverses through the saddle) and
sign(f.B) FLATTENS them, so one common series current cancels itself.  Measured
(L=1.0 former): Z2 single-wire 88% (sign(f.B)) -> 5.6% grad-psi loop-set (uniform
N=20) / 0.98% with --optimize-levels -> ~6.7% DELIVERED single-stroke wire
(grad-psi + auto-resolution field_aware chain).  Gx (l=1) is IDENTICAL either way
(grad-psi is a strict generalisation).  This was an orientation BUG, not a physics
limit -- single-wire l>=2 shims DO work.  `_gradpsi_signs` / `_build_gradpsi_kdt`
in calc_streamfunction.py; verified by
validation_test/stream_function/verify_gradpsi_orientation.py.
  DELIVERED gap (loop-set ~0.01 vs single-stroke wire ~0.067) = the inter-region
  CONNECTOR (bridge) field, NOT a loop-set error.  It was mostly chain UNDER-
  RESOLUTION: the old fixed cut-opt (ncut=24,passes=4) left clustered-level Z2 at
  ~0.20; auto-scaling resolution with loop count (see 'panel' CHAIN (3)) recovers
  ~0.067 (3x).  CHARACTER of the remaining ~7x gap (decomposed, 2026-06-20): the
  unit-current connector field is TINY (~7e-4 of ||B||) but the coil runs at high
  current, which AMPLIFIES it to the ~6.6% residual.  It is ROBUST to visit-ORDER:
  a lobe-aware reorder (group by current sense, cross at the saddle to shorten the
  long inter-lobe rung) does NOT help -- the saddle sits at the DSV centre, so a
  SHORT crossover there contaminates the target MORE than a long one routed away
  (the "shorter rungs != better field" trap); the cut-opt already finds the field-
  best of {nn, 2-opt}.  Splitting into one wire PER current group is also WORSE
  (~0.85): the cut-opt nulls connectors GLOBALLY by cross-cancelling the +/- rungs.
  LEVERS that DO work: --distort bends the one wire against the discrete Biot-Savart
  (Z2 0.067 -> 0.052 grid3/iter5 ~5mm, -> 0.045 grid4/iter10 ~7mm; sphere Z2
  4.3%->0.36%) -- a field/manufacturability trade (more bend), so it is a user knob,
  NOT auto-scaled.  For sub-1% l>=2 use INDEPENDENT feeds (a driven ARRAY:
  loops_homogeneity_rms ~ 0, the --pin-tiling regime), NOT one series wire.

(1) EQUAL-deltaI CONTOURS  (default, --nlevels N)
    N iso-contours at equal psi increments = N equal-current turns; the
    asymptotic large-N continuum.  KEY PITFALL: the contour-quantise turn count
    is NON-MONOTONIC (N=2 equal-current residual 0.24 can BEAT N=20's 0.43), so
    a "bisect for the minimum N" search is UNSOUND -- use greedy (3) for a
    monotone min-turns answer.

(2) --optimize-levels  (commit f9a9121a)
    Optimise the N ordered contour LEVELS (vs the uniform equal-deltaI default)
    to minimise the EQUAL-CURRENT discrete-loop residual.  Coordinate descent on
    the ordered levels (each 1-D bounded between its neighbours, per-level field
    cached); keeps clean, manufacturable contours.  Gain (nlevels=5, Gx):
    0.68->0.48 .. 0.77->0.24.  Complements --distort (position) + --chain.

(3) --greedy-turns N  (greedy constructive, the MONOTONE low-turn strategy;
    commit 39780894)
    "Put up a pin, lay the next turn where it helps most."  Matching pursuit with
    a COMMON series current: each step adds the candidate loop (with sign) that
    best ALIGNS the unit combined field with the target.  MONOTONE by
    construction -> the rms-vs-turns trace directly gives the fewest turns for
    any spec.  --greedy-target R stops once the residual <= R.
      --greedy-dict {contour, pin, bubble}
        contour : dense psi iso-contour loops (the manufacturable default)
        pin     : ON-SURFACE k-ring loops at surface pins (mesh verts at graph
                  distance k, angle-ordered in the tangent plane; k=1 1-ring,
                  k=2,3 adaptive size)
        bubble  : bubble-packed pins, density ~ |grad psi| (= the (4) dictionary)
      --greedy-connector-weight w   (commit 9fd546e1)
        CONNECTOR-AWARE selection: score = field_residual + w * gap/coil_diag,
        gap = nearest distance from the candidate's centroid to a loop already
        laid (0 for the first/stacked turn).  Keeps only field-improving turns
        so the trace stays MONOTONE; biases toward SHORT single-stroke connectors
        ("rungs") -> a manufacturable wire.  Verified: w=0.6 cuts connector
        length 0.502->0.447 m at an unchanged field residual.  w=0 (default) =
        pure-field.
      --greedy-plot p.png   turn-budget UI: the monotone rms-vs-turns curve +
        the --greedy-target spec line + the min-turns mark (no in-figure title).
    HONEST FINDING: greedy with ONE common current SATURATES fast (single-
    direction alignment).  The pin dictionary's poor chained wire_homo (~0.98 vs
    contour ~0.28) is from chaining FEW ISOLATED loops (long connectors), NOT
    loop size -- connector-aware (w>0) is the fix; the contour dictionary stays
    the manufacturable low-turn choice.

(4) --pin-tiling  (dense bubble-tiling pin coil; commit 9fd546e1)
    Equal-current k-ring loops packed over the surface with DENSITY proportional
    to |K| = |grad psi| (Hirahatake/Noguchi/Igarashi/Yamashita bubble system,
    clearance r ~ 1/sqrt|K|, k bisected to --pin-tiling-pins) and each loop sized
    to its bubble (--pin-tiling-frac).  HONEST REGIME (reported in the result +
    pin_tiling_note): isolated equal-current bubble loops reproduce the target
    under INDEPENDENT currents (loops_homogeneity_rms ~ 4.7e-15 -- the driven pin
    / shim ARRAY regime) but NOT as one series wire (wire_homogeneity_rms ~ 0.85),
    because isolated loops do NOT tile-cancel into the surface current K the way
    psi iso-contours do.  USE --pin-tiling for a DRIVEN-ARRAY layout (shim /
    matrix coil, independent currents); use iso-contours (default) or
    --greedy-turns for a single wound wire.

(5) --distort  (sheet-metal WIRE distortion; already shipped)
    Bend the ONE chained single-stroke wire (psi / levels fixed, one current)
    against the ACTUAL discrete-wire Biot-Savart to cancel the single-stroke
    residual without extra feeds.  Measured: plane 12290->340 ppm; cylinder Gx
    8.5-9.3%->1.4%; sphere Z2 4.3%->0.36%.  (See 'single_stroke'.)

(6) SINGLE-PASS (N=1 extreme)  -- demos sp1/sp2/sp3
    Inverse design through the iron-aware map (sp1); WITH feed leads at a
    near-coil workpiece (sp2, IH geometry -- lead field is first-order); Fourier
    lead compensation is FUNDAMENTALLY LIMITED (sp3: 46%->39%) -> geometric
    cancellation (bifilar / coaxial leads) or MATERIAL (flux concentrator over
    the leads, via --iron-vol) is the lever, not coil reshaping.

Metrics in the manufacture result: homogeneity_rms (continuous psi design ref),
loops_homogeneity_rms (N separate turns, independent currents -- best
discretisation / driven-array regime), wire_homogeneity_rms (single-stroke one
series current -- the delivered wire; the primary `rms`).  See 'panel' for
--chain / --confine and the connector ("rung") physics.

Goldens: test_streamfunction_greedy_golden.py (contour/pin/target/plot +
connector-aware), test_streamfunction_pin_tiling_golden.py (array regime /
density / bubble dict), test_streamfunction_levels_golden.py.
"""


TOPICS = {
    "overview": "SF coil-design framework + pipeline + topic map (this server's front door)",
    "material_aware": "SF coil design WITH iron (--iron-vol): Kelvin-FEM reaction folded into the whole pipeline, obs-adjoint scalability, exact-quad source, why-not-FMM",
    "low_turn": "low/many-turn manufacture: greedy constructive (--greedy-turns, monotone) + dicts contour/pin/bubble + connector-aware + turn-budget plot; --pin-tiling driven array; --optimize-levels; --distort; single-pass",
    "harmonics": "spherical-harmonic target (--target-harmonic Z2/X/..) + achieved-field (l,m) decomposition: purity / contamination (MRI gradient/shim basis)",
    "shielding": "active shielding (--shield-vol): primary+shield joint design (Turner 1986), ~86x stray reduction; coverage-is-everything lesson",
    "panel": "the design/pareto/manufacture PANEL + calc_streamfunction.py CLI",
    "boundary_conditions": "current confinement BC --confine off/on/abe (Abe edge-equipotential)",
    "contour": "contour=flux-line principle; --contour-sub order-p + --flux-plot bubble view",
    "theory": "SFM + (ACA+)+TSVD math (least-norm A psi = B, K = n x grad psi)",
    "method": "alias of theory",
    "api": "radia.stream_function API: aca_tsvd, RegularizedTSVD, pseudo_inverse_solve",
    "kernel_agnostic": "callback matrix-entry contract -- coils (Biot-Savart) OR magnets (MMM)",
    "regularized": "regularisation menu (L2/H1/sigma/inductance/L-inf) + folded Tikhonov closed form",
    "pareto": "(homogeneity, peak current density) Pareto front + 4 levers + sheet-metal",
    "sheet_metal": "sheet-metal: wire distortion + surface forming (-> pareto / single_stroke)",
    "single_stroke": "one-wire chain, FE-direct on arbitrary formers, wire distortion",
    "deformation": "surface-deformation outer loop (accuracy) -- see also pareto for the peak front",
    "fe_direct": "FE-direct H1 psi on plane / cylinder / sphere / conformal formers",
    "cmaes": "SA-25-020 CMA-ES outer loop",
    "performance": "(ACA+)+TSVD amortisation numbers",
    "validation": "analytic-benchmark validation",
    "literature": "SFM lineage (Turner / Peeren / current potential method)",
    "workflow": "end-to-end demo recipes",
    "fusion": "stellarator Stage-2 (NESCOIL/REGCOIL/FOCUS): winding-surface current potential -> plasma B.n, net-current secular term (cohomology), coil force/stress, real VMEC boundary, FOCUS winding-distance + winding-SHAPE; PLUS the radia-SF-vs-REGCOIL SUPERSET positioning (design parity; ACA+ridge-TSVD scale; fast-TSVD Pareto over alpha+shape; STEP+PEEC deliverable; iron) earned by goldens (a)(b)(c)",
    "clebsch_3d": "[IN PROGRESS] 3D stream function + cohomology: NGSolve de Rham unification (Clebsch/cohomology/SF/deform on H1->HCurl->HDiv); ACA+TSVD on CURRENT-Clebsch J=grad(lambda)xgrad(mu) with mu fixed; Tampere bidirectional (A,phi) map -> 3D flux coords, Jacobian=1/(|B||H|); honest non-convex+helicity frontier",
}


# SF-side remap for topic names that have no direct aca_tsvd topic/alias.
_REMAP = {
    "theory": "method",
    "deformation": "single_stroke",   # surface-deformation / sheet-metal content
    "deform": "single_stroke",
    "benchmarks": "literature",
}


def get_streamfunction_documentation(topic: str = "overview") -> str:
    """Return SF coil-design knowledge for ``topic``.

    ``overview`` returns this server's dedicated front-door page; every other
    topic dispatches to the local ``knowledge.aca_tsvd`` module (~80 aliases).
    """
    t = (topic or "overview").strip().lower()
    if t in ("overview", "index", "", "sf", "streamfunction", "stream_function",
             "stream-function", "front_door", "home"):
        return SF_OVERVIEW
    if t in ("panel", "calc", "calc_streamfunction", "gui", "modes",
             "design", "manufacture", "boundary_conditions", "boundary",
             "bc", "confine", "confinement", "abe", "edge_equipotential",
             "contour", "contours", "contour_sub", "order_p", "flux",
             "flux_line", "flux_lines", "bubble", "bubble_system",
             "cross_codebase", "validation_panel", "chain"):
        return SF_PANEL
    if t in ("material_aware", "material", "iron", "iron_vol", "iron_aware",
             "with_iron", "yoke", "core", "shield_iron", "reaction",
             "iron_reaction", "kelvin_reaction", "obs_adjoint", "exact_source",
             "iron_exact_source", "mu_r"):
        return SF_MATERIAL_AWARE
    if t in ("low_turn", "low_turns", "few_turn", "few_turns", "turn_count",
             "turns", "greedy", "greedy_turns", "greedy_dict", "greedy_target",
             "greedy_plot", "greedy_connector", "connector_aware", "connector",
             "pin", "pin_tiling", "pin_coil", "bubble", "bubble_tiling",
             "shim_array", "driven_array", "optimize_levels", "level_opt",
             "discrete_refinement", "single_pass", "one_pass", "1pass",
             "distort", "sheet_metal_wire"):
        return SF_LOW_TURN
    if t in ("harmonics", "harmonic", "spherical_harmonic", "spherical_harmonics",
             "spherical", "solid_harmonic", "solid_harmonics", "target_harmonic",
             "gradient", "gradients", "shim", "shims", "purity", "contamination",
             "mri", "z2", "decompose", "decomposition", "lm", "ylm"):
        return SF_HARMONICS
    if t in ("shielding", "shield", "active_shield", "active_shielding",
             "shielded", "primary_shield", "primary_and_shield", "stray",
             "stray_field", "turner", "mansfield", "shield_vol", "cryostat",
             "eddy_shield"):
        return SF_SHIELDING
    if t in ("fusion", "regcoil", "nescoil", "focus", "stellarator",
             "winding_surface", "winding_shape", "secular", "net_current",
             "cohomology", "vmec", "wout", "coil_force", "coil_stress",
             "force", "stress", "plasma", "stage2", "stage_2",
             "current_potential"):
        return SF_FUSION
    if t in ("clebsch_3d", "clebsch", "clebsch_potential", "3d", "3d_sf",
             "3d_stream_function", "3d_streamfunction", "volume_clebsch",
             "current_clebsch", "field_clebsch", "bidirectional",
             "bidirectional_map", "potential_coordinates", "flux_coordinates",
             "flux_coords", "de_rham", "derham", "feec", "feec_unification",
             "unification", "aca_clebsch", "aca_tsvd_clebsch", "tampere",
             "clebsch_cohomology", "cohomology_3d", "foliation",
             "foliated_clebsch", "helicity"):
        return SF_CLEBSCH_3D
    return get_aca_tsvd_knowledge(_REMAP.get(t, t))
