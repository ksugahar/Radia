"""Stream-function (SF) coil-design knowledge for the radia-streamfunction MCP
server.

The detailed, self-contained knowledge lives in THIS server's own
``radia_mcp.streamfunction.knowledge.aca_tsvd`` (overview / method / api /
kernel_agnostic / performance / cmaes / validation / literature / workflow /
single_stroke / regularized -- the last carries the folded-Tikhonov + Pareto
front + sheet-metal (板金) lever section).  It was MOVED here from
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
    coil resonates 195kHz; connects SF designer to radia-ih.  [topic: panel]
  * Folded TIKHONOV (the "+ alpha I" core) -> the (field homogeneity, PEAK
    current density) PARETO FRONT swept at ~50 us/point, with FOUR stackable
    levers: Tikhonov alpha (L-curve), L-inf minimax, geometry (former size /
    cylinder length), and SHEET-METAL (板金) surface forming.  [topic: pareto]
  * Single-stroke ("一筆書き") chain -- connect the contour family into ONE
    continuous wire with least stray-field impact (field_aware / Kuijpers).
    [topic: single_stroke]
  * Sheet-metal (板金), TWO distinct kinds:
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
  single_stroke     one-wire chain, FE-direct on arbitrary formers, wire 板金
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
  fusion            stellarator Stage-2 (REGCOIL/NESCOIL/FOCUS): winding-surface
                    current potential, net current, force/stress, VMEC, winding-shape

DOCS + DEMOS
------------
  docs/stream_function/  (README, theory, regularization, single_stroke,
    deformation, examples, api, benchmarks)
  examples/stream_function/demo_*.py  -- incl. the Pareto + sheet-metal
    demos: demo_pareto_tikhonov_aca / demo_pareto_geometry_nsga /
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

  design       target -> A psi = B (folded-Tikhonov RegularizedTSVD) -> psi,
               field homogeneity over the eval region, peak surface current.
  pareto       (homogeneity, peak current density) front via --pareto-lever:
                 alpha     Tikhonov L-curve (one factorisation, swept alpha)
                 linf      minimax IRLS trajectory (peak down at ~const homo)
                 geometry  eval-region (DSV) scale sweep (dual of former size)
  manufacture  psi iso-contours -> orientation-consistent equal-current turns
               -> field-aware single-stroke wire -> (sheet-metal --distort)
               -> CAD STEP (--step-output) -> PEEC L,R (--peec).

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
                    (2) single-stroke (一筆書き) wire, (3) sheet-metal (板金)
                    --distort wire, (4) wire WITH thickness (太さ, --wire-diam,
                    twist-free parallel-transport tube) + distortion.

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
field_aware is never worse than nn; with closed contours (--confine abe) it
reaches the separate-turns floor with no --distort.

END-TO-END VALIDATION vs an INDEPENDENT codebase
------------------------------------------------
examples/stream_function/verify_coil_field_independent.py designs a coil
(MRI-gradient scale: cylinder r=0.15 m, L=0.5 m, DSV r=0.05 m) and checks the
field TWO ways: the numpy straight-segment Biot-Savart used in the designer AND
Radia's C++ rad.ObjFlmCur + rad.Fld (a separate codebase).  They agree to
8-11 digits (uniform 3.5e-11, Gx 1.1e-8); the abe-confined Gx coil reaches
1.0 % nonlinearity on the short former, cross-validated.  Locked by
tests/panels/test_streamfunction_golden.py and ..._panel_qt.py.
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

DEMOS
-----
  examples/stream_function/demo_regcoil_fusion.py            (4 parts)
  examples/stream_function/demo_regcoil_fusion_advanced.py   (4 more parts)

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
     Betti number b1 = 2 (genus 1), CONFIRMED via Gmsh homology
     (addHomologyRequest('Cohomology',[pg],[],[1]); computeHomology() -> 2) --
     the engine wrapped by src/radia/cohomology_cut.py (which is a VOLUME
     scalar-potential solver, so the honest integration is "same engine,
     analytic generators on the torus": grad(zeta), grad(theta) are
     single-valued).  K_zeta = n x grad(zeta) = the net-POLOIDAL-current (TF)
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

VERIFIED (cylinder fixture, order 2, confine abe;
examples/stream_function/demo_shim_coil_purity.py): --target-harmonic X ->
dominant X, purity 1.000, residual 1.5e-4, Z2X contaminant 7e-5; Z2 ->
purity 1.000; the 4th-order Z4 shim -> purity 0.99983, residual 1.5e-2,
named Z3 contaminant 9.5e-3 (high-l shims are harder: an l-th harmonic's
field scales as r^l over a fixed DSV).  The named-basis harmonicity
(Laplacian 0, all 25 entries) and the target<->decompose round-trip are
golden-locked (tests/panels/test_streamfunction_golden.py
test_harmonic_basis_is_harmonic / test_harmonic_l4_forms_and_decompose).
The panel auto-generates the two flags (cli-diff clean).
"""


TOPICS = {
    "overview": "SF coil-design framework + pipeline + topic map (this server's front door)",
    "harmonics": "spherical-harmonic target (--target-harmonic Z2/X/..) + achieved-field (l,m) decomposition: purity / contamination (MRI gradient/shim basis)",
    "panel": "the design/pareto/manufacture PANEL + calc_streamfunction.py CLI",
    "boundary_conditions": "current confinement BC --confine off/on/abe (Abe edge-equipotential)",
    "contour": "contour=flux-line principle; --contour-sub order-p + --flux-plot bubble view",
    "theory": "SFM + (ACA+)+TSVD math (least-norm A psi = B, K = n x grad psi)",
    "method": "alias of theory",
    "api": "radia.stream_function API: aca_tsvd, RegularizedTSVD, pseudo_inverse_solve",
    "kernel_agnostic": "callback matrix-entry contract -- coils (Biot-Savart) OR magnets (MMM)",
    "regularized": "regularisation menu (L2/H1/sigma/inductance/L-inf) + folded Tikhonov closed form",
    "pareto": "(homogeneity, peak current density) Pareto front + 4 levers + sheet-metal (板金)",
    "sheet_metal": "sheet-metal (板金): wire distortion + surface forming (-> pareto / single_stroke)",
    "single_stroke": "one-wire (一筆書き) chain, FE-direct on arbitrary formers, wire distortion",
    "deformation": "surface-deformation outer loop (accuracy) -- see also pareto for the peak front",
    "fe_direct": "FE-direct H1 psi on plane / cylinder / sphere / conformal formers",
    "cmaes": "SA-25-020 CMA-ES outer loop",
    "performance": "(ACA+)+TSVD amortisation numbers",
    "validation": "analytic-benchmark validation",
    "literature": "SFM lineage (Turner / Peeren / current potential method)",
    "workflow": "end-to-end demo recipes",
    "fusion": "stellarator Stage-2 (NESCOIL/REGCOIL/FOCUS): winding-surface current potential -> plasma B.n, net-current secular term (cohomology), coil force/stress, real VMEC boundary, FOCUS winding-distance + winding-SHAPE",
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
    if t in ("harmonics", "harmonic", "spherical_harmonic", "spherical_harmonics",
             "spherical", "solid_harmonic", "solid_harmonics", "target_harmonic",
             "gradient", "gradients", "shim", "shims", "purity", "contamination",
             "mri", "z2", "decompose", "decomposition", "lm", "ylm"):
        return SF_HARMONICS
    if t in ("fusion", "regcoil", "nescoil", "focus", "stellarator",
             "winding_surface", "winding_shape", "secular", "net_current",
             "cohomology", "vmec", "wout", "coil_force", "coil_stress",
             "force", "stress", "plasma", "stage2", "stage_2",
             "current_potential"):
        return SF_FUSION
    return get_aca_tsvd_knowledge(_REMAP.get(t, t))
