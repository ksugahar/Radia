"""MCP Server: radia_mcp.maglev  (magnetic levitation)

Magnetic levitation knowledge -- UNIFIED.  Covers both the maglev
SYSTEMS (EMS/EDS trains, SCMaglev, Halbach/Inductrack, magnetic wheels,
PM/SC bearings) AND the levitation FORCE physics (induction/eddy-current
lift, EML melting, active magnetic bearings, superconducting
Meissner/pinning, diamagnetic, Earnshaw + loopholes, force computation).

Headline content is the lab's own Radia-based maglev research line
(CAE-AI Lab, Yano + Sugahara): Radia IEM <-> reduced-potential FEM weak
coupling for moving-magnet eddy-current force (topic radia_iem_fem) and
Cauer Ladder Network model-order reduction for real-time control-coupled
maglev (topic cln_mor_control).  The former linear-drive material
(LIM/LSM, end effects) was removed; the former separate
radia_mcp.levitation server was consolidated into this one.

Cross-references:
- `radia_mcp.ih` — induction heating (EML = levitation + IH together)
- `radia_mcp.mor` (mor_cln) — Cauer Ladder Network MOR theory
- `radia_mcp.fem` (potential_formulations) — A-phi / T-Omega / A-T gauges
- `radia_mcp.team_benchmark.force_motion.problem_28` — TEAM 28 levitation benchmark
- `radia_mcp.motor` — analogous rotary motor knowledge

Usage:
    mcp-server-maglev              # stdio
    mcp-server-maglev --selftest   # self-test
"""
import json
import sys
from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool
from .knowledge import get_knowledge, TOPICS
from .periodic_settling_gate import (
    rotating_conductor_periodic_settling_gate as _rotating_conductor_periodic_settling_gate,
)
from .team28_dynamic_gate import (
    team28_cycle_averaged_motion_gate as _team28_cycle_averaged_motion_gate,
)

mcp = FastMCP("mcp-server-maglev")


@mcp.tool()
def maglev(topic: str = "overview") -> str:
    """
    Magnetic levitation knowledge -- maglev systems + levitation force physics.

    Args:
        topic: One of:
          maglev systems + the lab's Radia/CLN research --
            "overview"            - Unified maglev landscape (DEFAULT)
            "radia_iem_fem"       - Radia IEM <-> reduced-potential FEM weak coupling (Yano)
            "cln_mor_control"     - Cauer Ladder Network MOR for control-coupled maglev (Yano)
            "team28_dynamic_scope"- 50 Hz cycle-average mechanical motion vs full EM transient
            "physical_tensor_rom" - Physical polarizability tensor alpha(s) as a passive LTI (AAA+NNLS; Kameari+Kelvin breakdown)
            "pm_maglev_zero_power"- Passive PM levitation, Maxwell-Earnshaw
            "eddy_current_maglev" - Eddy-current EDS, Kansai 2D model, Arago
            "sumitomo_heavy_industrial" - JP 7-327337 PM bearing + JP 2007-215264 mover
            "kansai_research"     - Saiki/Fujii magnetic-wheel lineage
            "scmaglev_eds"        - SCMaglev (Chuo Shinkansen) -- SC-EDS levitation
            "halbach_arrays"      - Halbach + Inductrack
          levitation FORCE physics --
            "induction_levitation"- Eddy-current (AC) lift; jumping/Thomson ring
            "eml_melting"         - Electromagnetic levitation melting (ties to IH)
            "magnetic_bearings"   - Active magnetic bearings (AMB), flywheels
            "superconducting"     - Meissner vs flux pinning, HTS bulk
            "diamagnetic"         - grad(B^2) levitation (graphite / water-frog)
            "earnshaw_stability"  - Earnshaw's theorem + its 5 loopholes
            "force_computation"   - Maxwell stress / virtual work / time-avg J x B
            "benchmarks"          - TEAM 28, jumping-ring analytic, EML lift coeff
            "all"                 - Everything
    """
    return get_knowledge(topic)


@mcp.tool()
def rotating_conductor_periodic_settling_gate(
    response: list[float],
    steps_per_period: int,
    angle_step_deg: float,
    reference_response: list[float] | None = None,
    maximum_final_period_relative_l2: float = 2.0e-3,
    maximum_contraction_factor: float = 0.35,
    reference_relative_l2_tolerance: float = 1.0e-8,
) -> str:
    """Gate full-turn convergence of a rotating-conductor eddy response."""

    try:
        result = _rotating_conductor_periodic_settling_gate(
            response,
            steps_per_period=steps_per_period,
            angle_step_deg=angle_step_deg,
            reference_response=reference_response,
            maximum_final_period_relative_l2=maximum_final_period_relative_l2,
            maximum_contraction_factor=maximum_contraction_factor,
            reference_relative_l2_tolerance=reference_relative_l2_tolerance,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "rotating_conductor_periodic_settling_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)


@mcp.tool()
def team28_cycle_averaged_motion_gate(
    summary: dict,
    claim_scope: str = "cycle_averaged_mechanical_motion",
    expected_frequency_hz: float = 50.0,
) -> str:
    """Gate TEAM 28 cycle-averaged motion evidence and reject transient overclaims."""

    try:
        result = _team28_cycle_averaged_motion_gate(
            summary,
            claim_scope=claim_scope,
            expected_frequency_hz=expected_frequency_hz,
        )
    except (TypeError, ValueError) as exc:
        result = {
            "policy": "team28_cycle_averaged_motion_gate_v1",
            "status": "invalid_input",
            "error": str(exc),
        }
    return json.dumps(result, indent=2, sort_keys=True)




register_status_tool(
    mcp,
    server_name='mcp-server-maglev',
    description='Magnetic levitation, UNIFIED: maglev systems (EMS/EDS/PM/SC/Halbach) + levitation FORCE physics (induction/EML/AMB/superconducting/diamagnetic/Earnshaw/force-computation). Lab research line: Radia IEM<->FEM weak coupling + Cauer Ladder Network MOR for control-coupled maglev (Yano, CAE-AI).',
    subpackage='radia_mcp.maglev',
    related_servers=["mor", "motor", "ih"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-maglev',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("maglev MCP server self-test:")
        print(f"  knowledge: {len(get_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
