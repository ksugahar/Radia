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


TOPICS = {
    "overview": "SF coil-design framework + pipeline + topic map (this server's front door)",
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
    return get_aca_tsvd_knowledge(_REMAP.get(t, t))
