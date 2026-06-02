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
  * Regularisation menu: L2 / H1 (min surface-current energy) / sigma-weighted
    (ohmic) / inductance-diagonal / L-inf (peak cap).  All fold onto ONE ACA
    factorisation via RegularizedTSVD.  [topic: regularized]
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
field_aware (default) chooses each loop's entry/exit CUT by coordinate descent
to minimise the FULL one-current wire error min_I ||I (loops+connectors) - B||
(NOT ||connectors|| alone -- that is worse on open contours).  Never worse than
nn; with closed contours it reaches the separate-turns floor with no --distort.

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


TOPICS = {
    "overview": "SF coil-design framework + pipeline + topic map (this server's front door)",
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
    return get_aca_tsvd_knowledge(_REMAP.get(t, t))
