"""Two-stage HDiv-MMM topology optimization policy.

``MMM-topology`` has exactly two optimization stages:

1. whole-element iron generation/removal proposed by ACA--QR--TSVD and
   accepted only after the exact Schur/full HDiv-MMM checks owned by
   :func:`radia.topology_optimization.grow_hdiv_mmm_by_superposition`;
2. topology-preserving NGSolve ``GetTrafo`` deformation with a complete
   physical re-solve for every accepted trial.

Coreform Cubit is the primary CAD/mesh rebuild route when a GetTrafo quality
gate requires a new mesh.  This facade does not select Sculpt.  Sculpt remains
available to separate voxel, volume-fraction, CT/microstructure, difficult-STL,
and independent remesh-validation workflows.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from .sheet_metal_optimization import (
    TopologyPreservingShapeResult,
    TopologyPreservingShapeState,
    optimize_topology_preserving_shape,
)
from .topology_optimization import (
    HDivMMMGenerationResult,
    grow_hdiv_mmm_by_superposition,
)


@dataclass(frozen=True)
class MMMTopologyPolicy:
    """Auditable policy label carried by every two-stage result."""

    name: str = "MMM-topology"
    lego_stage: str = "aca-qr-tsvd-exact-schur-binary-lego"
    shape_stage: str = "ngsolve-gettrafo-full-resolve"
    primary_cad_mesher: str = "coreform-cubit"
    fallback_mesher: str | None = None
    sculpt_in_core_loop: bool = False


DEFAULT_MMM_TOPOLOGY_POLICY = MMMTopologyPolicy()


@dataclass(frozen=True)
class MMMTopologyOptimizationResult:
    """Results from the binary Lego stage followed by the GetTrafo stage.

    ``final_target_accepted`` is deliberately optional.  Stage termination is
    not evidence that the requested engineering bands were met.  A caller that
    owns those bands supplies ``final_acceptance`` and receives an explicit
    Boolean; otherwise this field remains ``None``.
    """

    policy: MMMTopologyPolicy
    generation: HDivMMMGenerationResult
    shape: TopologyPreservingShapeResult
    final_target_accepted: bool | None

    @property
    def active_elements(self) -> np.ndarray:
        return self.generation.active_elements

    @property
    def final_evaluation(self):
        return self.shape.state.evaluation


def optimize_mmm_topology(
    *,
    lego_options,
    build_shape_state,
    linearize_shape_step,
    deformation_factory,
    rebuild_shape_model,
    evaluate_shape_model,
    move_limit,
    shape_options=None,
    final_acceptance=None,
    policy=DEFAULT_MMM_TOPOLOGY_POLICY,
) -> MMMTopologyOptimizationResult:
    """Run the canonical two-stage ``MMM-topology`` workflow.

    ``lego_options`` is passed to the production HDiv-MMM whole-element
    generation/removal driver.  ``build_shape_state`` converts its exact final
    binary state into the topology-fixed shape model and parameterization.
    The remaining callbacks are the existing GetTrafo shape-driver contract.

    A Cubit backend may be supplied through ``shape_options['cubit_backend']``
    when an application owns a CAD-quality rebuild.  No fallback mesher is
    installed by this facade.
    """
    if not isinstance(policy, MMMTopologyPolicy):
        raise TypeError("policy must be an MMMTopologyPolicy")
    if policy.sculpt_in_core_loop:
        raise ValueError("Sculpt is not part of the MMM-topology core loop")
    if not isinstance(lego_options, Mapping):
        raise TypeError("lego_options must be a mapping")
    if shape_options is None:
        shape_options = {}
    if not isinstance(shape_options, Mapping):
        raise TypeError("shape_options must be a mapping")

    generation = grow_hdiv_mmm_by_superposition(**dict(lego_options))
    if not isinstance(generation, HDivMMMGenerationResult):
        raise RuntimeError(
            "the Lego stage did not return an HDivMMMGenerationResult")

    initial_shape_state = build_shape_state(generation)
    if not isinstance(initial_shape_state, TopologyPreservingShapeState):
        raise TypeError(
            "build_shape_state must return TopologyPreservingShapeState")

    reserved = {
        "linearize_step",
        "deformation_factory",
        "rebuild_model",
        "evaluate_model",
        "move_limit",
    }
    overlap = reserved.intersection(shape_options)
    if overlap:
        raise TypeError(
            "shape_options cannot override the MMM-topology stage contract: "
            + ", ".join(sorted(overlap)))

    shape = optimize_topology_preserving_shape(
        initial_shape_state,
        linearize_step=linearize_shape_step,
        deformation_factory=deformation_factory,
        rebuild_model=rebuild_shape_model,
        evaluate_model=evaluate_shape_model,
        move_limit=move_limit,
        **dict(shape_options),
    )
    if not isinstance(shape, TopologyPreservingShapeResult):
        raise RuntimeError(
            "the GetTrafo stage did not return TopologyPreservingShapeResult")

    accepted = None
    if final_acceptance is not None:
        value = final_acceptance(generation, shape)
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError("final_acceptance must return a Boolean")
        accepted = bool(value)

    return MMMTopologyOptimizationResult(policy, generation, shape, accepted)


__all__ = [
    "MMMTopologyPolicy",
    "DEFAULT_MMM_TOPOLOGY_POLICY",
    "MMMTopologyOptimizationResult",
    "optimize_mmm_topology",
]
