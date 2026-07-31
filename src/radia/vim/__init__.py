"""radia.vim -- HDiv-type VIM demag operator.

The FEEC H(div) BDM production demag path: a SYMMETRIC demag
operator N = B^T G B whose loop modes are field-null by construction (de Rham), giving mu_r-independent
convergence with no hand-crafted loop-star.

This is the PRODUCTION home for the validated HDiv-VIM core.  Canonical docs:
docs/hdiv_vim/README.md; roadmap:
docs/hdiv_vim/PRODUCTIONIZATION.md.

Public API (NGSolve-aligned validated solve primitives):
  Solve(mesh, mu_r=/bh_table=, H_ext=.., B_r=..)
      -> the production one-call HDiv-VIM demag solve.  This is the preferred method-level entry,
         matching the NGSolve convention that solve/operator objects use short CamelCase names inside
         their owning namespace.  Scalar recoil mu_r plus a constant/spatial B_r selects the level-2
         linear permanent-magnet law.
  HDivSolver(mesh, order=1|2, ...)
      -> persistent 3D solve object for load sweeps, history, and coupled
         bodies.  It owns the HDiv space and geometry-only ChargeGram; result
         dictionaries contain no reusable-operator handle.
  SolveHysteresis(mesh, h_steps, play=(K, eta, f_k_tables), order=1|2)
      -> quasi-static B-input hysteresis stepping on the SAME BDM charge Gram: the chi-free
         H-matrix is built ONCE and reused by every step / nonlinear iteration.  BDM1 keeps
         element-average states; BDM2 uses NGSolve IntegrationRuleSpace quadrature states.
         Committed material states advance only after each converged step.
         EnergyStopMaterial supplies the production C++ vector B-input Stop law for hard magnets;
         explicit state supports irreversible demagnetization, recoil, restart, and persistent
         final-field evaluation.
  DemagOperator(HDiv(mesh, order=1), intorder=, eps=)
      -> an ngsolve.bem-style operator.  `.mat` is the H-matrix-backed NGSolve BaseMatrix
         N = B^T G B, which composes with NGSolve solvers / BlockMatrix exactly like ngsolve.bem's
         SingleLayerPotentialOperator.  `.DemagFactor(M_cf)` -> the demag factor (~1/3).
  ChargeGram(HDiv(mesh, order=1), ...)
      -> (B, G, M_mass), the charge map, charge-Gram H-matrix, and HDiv mass used by DemagOperator.
  MagnetizationSource(mesh, M_given, order=1)
      -> a fixed/given magnetization L2-projected into its own HDiv space, with a persistent native C++
         H-field CoefficientFunction for coupling to a separate soft-iron Solve.
  SolveCoupled([CoupledBody(...), ...])
      -> block Gauss-Seidel coupling for independent linear-recoil PM and linear/nonlinear
         iron spaces.  Each geometry-only ChargeGram is built once and reused.
  SolveCoupledHysteresis(history_bodies, bodies, h_steps)
      -> one or more stateful EnergyStop/Play PMs plus independent
         linear/nonlinear bodies.  Each outer trial restarts every body from
         its committed history; only the converged all-body state commits.

Permanent-magnet levels are canonical: (1) MagnetizationSource fixed/given M,
(2) Solve with recoil mu_r+B_r, (3) PlayHysteresisMaterial, and
(4) EnergyStopMaterial.  They are physical model levels, not compatibility aliases.
  MeshSoftIron(mesh, mu_r=/bh_table=, order=1|2) / VolSoftIron(path, mu_r=/bh_table=, order=1|2)
      -> method-layer constructors for mesh-backed Radia soft iron.  For ordinary user code prefer the
         user-intent API `rad.SoftIron(geometry, mu_r=...).solve(...)`.  When `rad.Solve(..., image=...)`
          is used on a MeshSoftIron, `rad.Fld(iron, ...)` evaluates the solved BDM field and its reflected
          IMA contributions directly; `M_avg_reduced` is the
         reduced-domain diagnostic and `M_avg` is the physical full-domain average.  Unconstrained explicit
         full-solve `rad.Fld` parity is a separate 10-eps validation target, not a percent-level tolerance.
  PlanarSolve(...) and PlanarDemagBody(...)
      -> the 2D planar tri/quad layer.

NOTE: the production HDiv-VIM solve requires Radia's C++/pybind core for the
charge Gram, persistent matrix/solvers, prescribed-source field CF, and batch
field evaluator; Python declares the NGSolve spaces and forms.
"""
import inspect as _inspect

from . import _nonlinear  # noqa: F401
from ._vim import (  # noqa: F401  (ngsolve.bem-style operator + .mat)
    DemagOperator,
    build_charge_gram as _charge_gram_impl,
)
from ._solve import hdiv_demag_solve as _solve_impl  # noqa: F401  (production demag solve)
from ._solver import HDivSolver  # noqa: F401
from ._vim2d import (  # noqa: F401  (2D planar motor-cross-section layer; vim.Solve dispatches here)
    PlanarDemagBody,
    maxwell_torque_circle,
    solve_planar_demag as _solve_planar_demag,
)
from ._radsolve import (  # noqa: F401  (.vol/mesh -> both-backend iron)
    soft_iron_from_mesh as _mesh_soft_iron_impl,
    soft_iron_from_vol as _vol_soft_iron_impl,
)
from ._hysteresis import (  # noqa: F401  (B-input hysteresis stepping: ONE Gram build, per-step W-CG)
    EnergyStopMaterial,
    PlayHysteresisMaterial,
    SolveHysteresis as _solve_hysteresis_impl,
)
from ._field_batch import (  # noqa: F401  (batch exterior field of the BDM1/BDM2 solution)
    field_coefficient_from_solution as _field_coefficient_from_solution_impl,
    field_from_solution as _field_from_solution_impl,
)
from ._magnetization_source import MagnetizationSource  # noqa: F401
from ._coupled import (  # noqa: F401
    CoupledBody,
    CoupledHistoryBody,
    field_from_coupled_hysteresis as _field_from_coupled_hysteresis_impl,
    field_from_coupled_solution as _field_from_coupled_solution_impl,
    solve_coupled_hysteresis as _solve_coupled_hysteresis_impl,
    solve_coupled as _solve_coupled_impl,
)
from ._shapes import soft_iron_box, soft_iron_hex, magnet_box, magnet_hex  # noqa: F401  (mesh-less-SHAPE intent constructors: soft iron -> HDiv-VIM; PM -> analytic)
from ._capabilities import HDivCapability, hdiv_capabilities  # noqa: F401
from ._eddy_hybrid import (  # noqa: F401  (reduced HCurl/T + surface-Omega/SIBC eddy-current VIM primitives)
    MU0,
    SkinDepth,
    EddySIBCApplicability,
    EddyTracePolynomialDim,
    EddyParentOrderLedger,
    SampledCurrentBasis,
    SampledMagnetizationBasis,
    VolumeCurrentBasis,
    MagnetizationBasis,
    SurfaceOmegaBasis,
    EddyFaceTopology,
    EddyConductorGraphEdge,
    EddyConductorCycle,
    EddyConductorGraph,
    EddyMeshTopology,
    EddyDofPolicy,
    EddyReductionPlan,
    EddyBubbleDecomposition,
    EddyBubbleHCurlBasis,
    HCurlCellFamilyInventory,
    EddyBubbleReduction,
    NgsolveHCurlCellFamilies,
    ClassifyNgsolveEddyTopology,
    NgsolveEddyDofPolicy,
    NgsolveEddyBubbleReduction,
    NgsolveEddyBubbleHCurlBasis,
    NgsolveBridgeCycleCurrentBasis,
    SampleNgsolveVectorCFs,
    NgsolveVolumeCurrentBasis,
    NgsolveHCurlVectorPotentialPort,
    NgsolveMagnetizationBasis,
    NgsolveHDivMagnetizationBasis,
    HDivMultipolePortSet,
    PlanarHarmonicPortSet,
    NgsolveHDivRegularSolidHarmonicPorts,
    NgsolvePlanarHarmonicPorts,
    NgsolveHDivExternalFieldRHS,
    HDivMMMReducedModel,
    NgsolveHDivMMMReduction,
    NgsolveHDivMMMResponseReduction,
    NgsolveBDMHDivMMMResponseReduction,
    PlanarHDivMMMReducedSolution,
    PlanarHDivMMMReducedModel,
    NgsolvePlanarHDivMMMResponseReduction,
    NgsolveHCurlCurlBasis,
    NgsolveSurfaceOmegaBasis,
    NgsolveMatrixToDense,
    NgsolveVectorToArray,
    NgsolveCouplingDofMasks,
    ResponseBasis,
    EVRSBasis,
    CompressHCurlResponseInCurrentGram,
    BlockKrylovBasis,
    NgsolveBlockKrylovBasis,
    NgsolveOperatorBlockKrylovBasis,
    NgsolveStaticCondensedBlockKrylovBasis,
    SampledLaplaceInteraction,
    HACApKSampledLaplaceInteraction,
    HACApKSampledPlanarLogInteraction,
    SampledHACApKOperator,
    ReducedInteractionMatrix,
    NGSolveProjectedOperator,
    NGSolveProjectedInteraction,
    CoupledReducedOperator,
    CurrentMagneticFluxDensitySamples,
    MagnetizationCurrentCoupling,
    EVRSTMethodAlgebra,
    ReducedPortAdmittance,
    ReducedPortImpedance,
    SharedMeshMaterialModel,
    HCurlVIMHDivMMMSolution,
    CoupledHDivEVRSSystem,
    CoupledHDivHybridVIMSystem,
    MixedGalerkinHDivHybridVIMSystem,
    HCurlVIMHDivMMMSystem,
    CoupleHDivMagnetizationToEVRS,
    CoupleHCurlVIMWithHDivMMM,
    CoupleHybridVIMWithHDivMMM,
    CoupleEddyBubbleHCurlBasisWithHDivMMM,
    MixedGalerkinOrthogonalization,
    SurfaceImpedanceGram,
    LocalESIMSurfaceLUT,
    LocalESIMSurfaceModel,
    LocalESIMSurfaceSolution,
    HybridVIMSystem,
    HCurlEddyCLNModel,
    HCurlEddyCLNFromVIM,
    AssembleHybridVIM,
    AssembleSurfaceImpedanceGram,
    BuildLocalESIMSurfaceLUT,
    ValidateLocalESIMSurfaceLUT,
    SolveLocalESIMSurfaceVIM,
    TopologyAwareHybridVIM,
    NgsolveTopologyAwareHybridVIM,
    NgsolveEddyBubbleHybridVIM,
    NgsolveHCurlVIMHDivMMM,
    NgsolveBDMEddyBubbleVIM,
    SkinImpedance,
    SIBCAdmittanceTail,
    SIBCSchurTerminationImpedance,
    SIBCSchurTerminationAdmittance,
    ExternalVectorPotentialRHS,
)
from ._hcurl_tet_interaction import (  # noqa: F401
    HCurlHMatrixOperator,
    HCurlTetVolumeInteraction,
    HCurlCellVolumeInteraction,
    NgsolveHCurlTetVolumeInteraction,
    NgsolveHCurlCellVolumeInteraction,
    SampleNgsolveHCurlCellSubtetVelocities,
)
from ._hcurl_planar_interaction import (  # noqa: F401
    HCurlPlanarVolumeInteraction,
    NgsolveHCurlPlanarVolumeInteraction,
)
from ._matlab_bridge import (  # noqa: F401
    ExportHCurlEddyCLNJSON,
    ExportHCurlEddyCLNFamilyJSON,
)
def Solve(*args, **kwargs):
    """NGSolve-style production HDiv-VIM one-call solve.
    """
    if any(name.startswith("_") for name in kwargs):
        raise TypeError(
            "vim.Solve does not accept private operator-cache arguments; "
            "use vim.HDivSolver for repeated solves")
    return _solve_impl(*args, **kwargs)


def SolveHysteresis(*args, **kwargs):
    """Quasi-static B-input history solve on one HDiv geometry."""
    if any(name.startswith("_") for name in kwargs):
        raise TypeError(
            "vim.SolveHysteresis does not accept private operator-cache arguments; "
            "use vim.HDivSolver.SolveHysteresis for continuation")
    return _solve_hysteresis_impl(*args, **kwargs)


def ChargeGram(*args, **kwargs):
    """NGSolve-style charge-Gram builder for an HDiv finite element space.

    Returns ``(B, G, M_mass)``.
    """
    return _charge_gram_impl(*args, **kwargs)


def MeshSoftIron(*args, **kwargs):
    """Build a mesh-backed Radia soft iron container from an NGSolve mesh.

    User-facing workflows should usually use ``rad.SoftIron(mesh, ...)``.
    """
    return _mesh_soft_iron_impl(*args, **kwargs)


def VolSoftIron(*args, **kwargs):
    """Build a mesh-backed Radia soft iron container from a Netgen ``.vol`` file.

    User-facing workflows should usually use ``rad.SoftIron(path, ...)``.
    """
    return _vol_soft_iron_impl(*args, **kwargs)


def PlanarSolve(*args, **kwargs):
    """NGSolve-style alias for the 2D planar HDiv-VIM solve."""
    return _solve_planar_demag(*args, **kwargs)


def FieldFromSolution(*args, **kwargs):
    """Batch demagnetizing H (A/m) from the full BDM1/BDM2 HDiv solution.

    Pass ``vim.Solve``'s result dict.  The solve owns a persistent C++ source
    evaluator; ordinary batches use exact flat/curved element leaves, while very
    large target-source work may use the accuracy- and timing-probed treecode.
    """
    return _field_from_solution_impl(*args, **kwargs)


def FieldCoefficientFromSolution(*args, **kwargs):
    """Persistent NGSolve demagnetizing-field CoefficientFunction."""
    return _field_coefficient_from_solution_impl(*args, **kwargs)


def SolveCoupled(*args, **kwargs):
    """Mutually coupled independent HDiv bodies with cached ChargeGrams."""
    return _solve_coupled_impl(*args, **kwargs)


def FieldFromCoupledSolution(*args, **kwargs):
    """Batch magnetization field summed over a SolveCoupled result."""
    return _field_from_coupled_solution_impl(*args, **kwargs)


def SolveCoupledHysteresis(*args, **kwargs):
    """Stateful PM plus independently meshed linear/nonlinear HDiv bodies."""
    return _solve_coupled_hysteresis_impl(*args, **kwargs)


def FieldFromCoupledHysteresis(*args, **kwargs):
    """Batch final field summed over a SolveCoupledHysteresis result."""
    return _field_from_coupled_hysteresis_impl(*args, **kwargs)


for _new, _old in [
    (Solve, _solve_impl),
    (SolveHysteresis, _solve_hysteresis_impl),
    (ChargeGram, _charge_gram_impl),
    (MeshSoftIron, _mesh_soft_iron_impl),
    (VolSoftIron, _vol_soft_iron_impl),
    (PlanarSolve, _solve_planar_demag),
    (FieldFromSolution, _field_from_solution_impl),
    (FieldCoefficientFromSolution, _field_coefficient_from_solution_impl),
    (SolveCoupled, _solve_coupled_impl),
    (FieldFromCoupledSolution, _field_from_coupled_solution_impl),
    (SolveCoupledHysteresis, _solve_coupled_hysteresis_impl),
    (FieldFromCoupledHysteresis, _field_from_coupled_hysteresis_impl),
]:
    _signature = _inspect.signature(_old)
    _new.__signature__ = _signature.replace(parameters=[
        parameter for parameter in _signature.parameters.values()
        if not parameter.name.startswith("_")
    ])

__all__ = [
    "Solve", "HDivSolver", "DemagOperator", "ChargeGram",
    "MeshSoftIron", "VolSoftIron", "PlanarSolve",
    "PlanarDemagBody", "maxwell_torque_circle",
    "soft_iron_box", "soft_iron_hex",
    "magnet_box", "magnet_hex",
    "SolveHysteresis", "EnergyStopMaterial", "PlayHysteresisMaterial",
    "FieldFromSolution", "FieldCoefficientFromSolution", "MagnetizationSource",
    "CoupledBody", "SolveCoupled", "FieldFromCoupledSolution",
    "CoupledHistoryBody", "SolveCoupledHysteresis", "FieldFromCoupledHysteresis",
    "HDivCapability", "hdiv_capabilities",
    "MU0", "SkinDepth", "EddySIBCApplicability",
    "EddyTracePolynomialDim", "EddyParentOrderLedger",
    "SampledCurrentBasis", "SampledMagnetizationBasis",
    "VolumeCurrentBasis", "MagnetizationBasis", "SurfaceOmegaBasis",
    "EddyFaceTopology", "EddyConductorGraphEdge", "EddyConductorCycle",
    "EddyConductorGraph", "EddyMeshTopology", "EddyDofPolicy",
    "EddyReductionPlan", "EddyBubbleDecomposition", "EddyBubbleHCurlBasis",
    "HCurlCellFamilyInventory",
    "EddyBubbleReduction",
    "ClassifyNgsolveEddyTopology", "NgsolveEddyDofPolicy", "NgsolveEddyBubbleReduction",
    "NgsolveHCurlCellFamilies",
    "NgsolveEddyBubbleHCurlBasis",
    "NgsolveBridgeCycleCurrentBasis", "SampleNgsolveVectorCFs", "NgsolveVolumeCurrentBasis",
    "NgsolveHCurlVectorPotentialPort",
    "NgsolveMagnetizationBasis", "NgsolveHDivMagnetizationBasis",
    "HDivMultipolePortSet", "PlanarHarmonicPortSet",
    "NgsolveHDivRegularSolidHarmonicPorts", "NgsolvePlanarHarmonicPorts",
    "NgsolveHDivExternalFieldRHS", "HDivMMMReducedModel", "NgsolveHDivMMMReduction",
    "NgsolveHDivMMMResponseReduction", "NgsolveBDMHDivMMMResponseReduction",
    "PlanarHDivMMMReducedSolution", "PlanarHDivMMMReducedModel",
    "NgsolvePlanarHDivMMMResponseReduction",
    "NgsolveHCurlCurlBasis", "NgsolveSurfaceOmegaBasis",
    "NgsolveMatrixToDense", "NgsolveVectorToArray", "NgsolveCouplingDofMasks",
    "ResponseBasis", "EVRSBasis", "CompressHCurlResponseInCurrentGram",
    "BlockKrylovBasis", "NgsolveBlockKrylovBasis",
    "NgsolveOperatorBlockKrylovBasis", "NgsolveStaticCondensedBlockKrylovBasis",
    "SampledLaplaceInteraction", "HACApKSampledLaplaceInteraction",
    "SampledHACApKOperator", "ReducedInteractionMatrix", "HybridVIMSystem",
    "SurfaceImpedanceGram", "AssembleSurfaceImpedanceGram",
    "LocalESIMSurfaceLUT", "LocalESIMSurfaceModel", "LocalESIMSurfaceSolution",
    "BuildLocalESIMSurfaceLUT", "ValidateLocalESIMSurfaceLUT",
    "SolveLocalESIMSurfaceVIM",
    "HACApKSampledPlanarLogInteraction", "NGSolveProjectedOperator",
    "NGSolveProjectedInteraction", "CoupledReducedOperator",
    "HCurlHMatrixOperator", "HCurlTetVolumeInteraction", "HCurlCellVolumeInteraction",
    "NgsolveHCurlTetVolumeInteraction", "NgsolveHCurlCellVolumeInteraction",
    "SampleNgsolveHCurlCellSubtetVelocities",
    "HCurlPlanarVolumeInteraction", "NgsolveHCurlPlanarVolumeInteraction",
    "HCurlEddyCLNModel", "HCurlEddyCLNFromVIM",
    "ExportHCurlEddyCLNJSON", "ExportHCurlEddyCLNFamilyJSON",
    "CurrentMagneticFluxDensitySamples", "MagnetizationCurrentCoupling",
    "EVRSTMethodAlgebra", "ReducedPortAdmittance", "ReducedPortImpedance",
    "SharedMeshMaterialModel", "HCurlVIMHDivMMMSolution",
    "CoupledHDivEVRSSystem", "CoupledHDivHybridVIMSystem",
    "MixedGalerkinHDivHybridVIMSystem",
    "HCurlVIMHDivMMMSystem",
    "CoupleHDivMagnetizationToEVRS", "CoupleHCurlVIMWithHDivMMM",
    "CoupleHybridVIMWithHDivMMM",
    "CoupleEddyBubbleHCurlBasisWithHDivMMM",
    "AssembleHybridVIM", "TopologyAwareHybridVIM", "NgsolveTopologyAwareHybridVIM",
    "NgsolveEddyBubbleHybridVIM", "NgsolveHCurlVIMHDivMMM",
    "NgsolveBDMEddyBubbleVIM",
    "SkinImpedance", "SIBCAdmittanceTail",
    "SIBCSchurTerminationImpedance", "SIBCSchurTerminationAdmittance",
    "ExternalVectorPotentialRHS",
    "_nonlinear", "_vim", "_solve", "_radsolve", "_hysteresis",
]
