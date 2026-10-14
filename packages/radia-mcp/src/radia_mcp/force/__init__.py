"""Shared electromagnetic-force knowledge and MCP entry points.

``radia_mcp.force`` is the application-neutral front door used by motor,
maglev, differential-forms, and radia-ngsolve.  Solver-specific field
assembly remains in the owning solver packages; this package owns method
selection, common post-processing entry points, and validation guidance.
"""

from .gates import (
    electromagnetic_force_method_selection_gate,
    force_action_reaction_gate,
    force_torque_method_agreement_gate,
    force_weight_equilibrium_gate,
)
from .knowledge import (
    TOPICS,
    get_force_extras,
    get_force_knowledge,
    get_force_methods,
    get_force_recipe,
    get_force_validation,
)

__all__ = [
    "TOPICS",
    "electromagnetic_force_method_selection_gate",
    "force_action_reaction_gate",
    "force_torque_method_agreement_gate",
    "force_weight_equilibrium_gate",
    "get_force_extras",
    "get_force_knowledge",
    "get_force_methods",
    "get_force_recipe",
    "get_force_validation",
]
