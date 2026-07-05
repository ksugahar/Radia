"""MCP Server: radia_mcp.team_benchmark

TEAM Workshop benchmark problems reference layer for EM solver validation.

Covers 30 problems organized by physics class:
- Magnetostatic (problems 6, 13)
- Eddy current (1a/1b, 2, 3, 4, 5, 7, 9, 21, 24)
- Force / motion / levitation (17, 20, 23, 28, 33b) ★ lab core
- NDT (8, 15, 27)
- Inverse / optimization (22, 25, 31)
- Hysteresis (32) ★ lab core
- HF cavity / motor (18, 19, 29, 30b, 34)

★ Problems 13, 20, 23, 32, 33b are heavily used in lab validation.

Distilled from 45 files / 449 MB in
public-safe curated corpus

Use cases:
- "Which TEAM problem validates Maxwell stress force?" → problem 20
- "Where's the geometry for problem 32 (vector hysteresis)?" → catalog
- "How does Hollaus MSFEM compare on problem 24?" → eddy current

Usage:
    mcp-server-team-benchmark              # stdio
    mcp-server-team-benchmark --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .catalog_knowledge import get_catalog_knowledge
from .magnetostatic_knowledge import get_magnetostatic_knowledge
from .eddy_current_knowledge import get_eddy_current_knowledge
from .force_motion_knowledge import get_force_motion_knowledge
from .ndt_inverse_knowledge import get_ndt_inverse_knowledge
from .special_knowledge import get_special_knowledge


mcp = FastMCP("mcp-server-team-benchmark")


@mcp.tool()
def team_catalog(topic: str = "overview") -> str:
    """
    TEAM Workshop benchmark catalog.

    Args:
        topic: One of:
            "overview"         - All 30 problems by physics class (DEFAULT)
            "validation_map"   - Which TEAM problem validates technique X
            "all"              - Everything
    """
    return get_catalog_knowledge(topic)


@mcp.tool()
def team_magnetostatic(topic: str = "problem_06") -> str:
    """
    TEAM magnetostatic problems: 6 (sphere), 13 (nonlinear yoke).

    Args:
        topic: One of:
            "problem_06"  - Sphere in Uniform Field (analytic ref) (DEFAULT)
            "problem_13"  - ★ 3-D Non-Linear Magnetostatic (lab core)
            "all"         - Both
    """
    return get_magnetostatic_knowledge(topic)


@mcp.tool()
def team_eddy_current(topic: str = "overview") -> str:
    """
    TEAM eddy current problems: 1a/1b, 2, 3, 4, 5, 7, 9, 21, 24.

    Args:
        topic: One of:
            "overview"     - All eddy current problems summary (DEFAULT)
            "problem_07"   - ★ Asymmetrical Conductor with Hole (canonical 3D)
            "problem_09"   - Velocity Effects (moving conductor)
            "problem_24"   - Nonlinear Time-Transient Rig
            "all"          - Everything
    """
    return get_eddy_current_knowledge(topic)


@mcp.tool()
def team_force_motion(topic: str = "overview") -> str:
    """
    TEAM force / motion / levitation: 17, 20, 23, 28, 33b ★ lab core.

    Args:
        topic: One of:
            "overview"     - All force problems summary (DEFAULT)
            "problem_20"   - ★ 3-D Static Force (lab core)
            "problem_23"   - ★ Forces in Permanent Magnets (lab core)
            "problem_33b"  - ★ Electric Local Force (per-node validation)
            "problem_28"   - Electrodynamic Levitation
            "all"          - Everything
    """
    return get_force_motion_knowledge(topic)


@mcp.tool()
def team_ndt_inverse(topic: str = "ndt_problems") -> str:
    """
    TEAM NDT and inverse / optimization problems.

    Args:
        topic: One of:
            "ndt_problems"  - Problems 8, 15, 27 (NDT) (DEFAULT)
            "inverse_opt"   - Problems 22, 25, 31 (inverse / optimization)
            "all"           - Both
    """
    return get_ndt_inverse_knowledge(topic)


@mcp.tool()
def team_special(topic: str = "hysteresis_32") -> str:
    """
    TEAM special problems: hysteresis (32), motors (30b/34), HF (18/19/29).

    Args:
        topic: One of:
            "hysteresis_32"   - ★ Problem 32 vector hysteresis (lab core)
            "motors_30b_34"   - Motor analytic references
            "high_freq"       - Problems 18, 19, 29 (HF, out of lab scope)
            "references"      - Key reference documents in _references/
            "all"             - Everything
    """
    return get_special_knowledge(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def find_validation_problem(feature: str) -> str:
    """Find the right TEAM problem to validate a solver feature."""
    guidance = {
        "maxwell_stress_force": (
            "Maxwell stress force computation:\n"
            "→ team_force_motion('problem_20')\n"
            "Problem 20 is THE canonical static force benchmark.\n"
            "All 7 lab-implemented force methods validate here.\n"
        ),
        "pm_force": (
            "Permanent magnet force:\n"
            "→ team_force_motion('problem_23')\n"
            "Problem 23 has multiple PM arrangements (attractive, repulsive,\n"
            "parallel, perpendicular) with measured forces.\n"
        ),
        "local_force_per_node": (
            "Per-node local force formulations:\n"
            "→ team_force_motion('problem_33b')\n"
            "Problem 33b validates: Kameari per-node, Iino-Okamoto two-stage,\n"
            "Eggshell, Sensitivity-VWP (Pile 2018).\n"
        ),
        "nonlinear_magnetostatic": (
            "Nonlinear magnetostatic FEM:\n"
            "→ team_magnetostatic('problem_13')\n"
            "Problem 13 is the canonical 3D nonlinear C-yoke test.\n"
            "Multiple published reference solutions (Picard, Newton).\n"
        ),
        "3d_eddy_current": (
            "3D harmonic eddy current solver:\n"
            "→ team_eddy_current('problem_07')\n"
            "Problem 7 (Asymmetrical Conductor with Hole) is the most-cited\n"
            "3D eddy current benchmark. 50 Hz, plate with off-center hole.\n"
        ),
        "vector_hysteresis": (
            "Vector hysteresis model (Energy-Based, Play, Stop):\n"
            "→ team_special('hysteresis_32')\n"
            "Problem 32 — rotating B excitation on T-core, measure H locus.\n"
            "Cross-reference: radia_mcp.magnetic_materials.hysteresis_models\n"
        ),
        "moving_conductor_eddy": (
            "Moving conductor eddy current (v × B):\n"
            "→ team_eddy_current('problem_09')\n"
            "Problem 9 — axial-motion bar in coil field.\n"
        ),
        "nonlinear_transient": (
            "Nonlinear time-transient FEM:\n"
            "→ team_eddy_current('problem_24')\n"
            "Problem 24 — actual rotational test rig at IGTE.\n"
            "Cross-reference: radia.calc_motor_transient\n"
        ),
        "ndt_eddy_current": (
            "Eddy current NDT solver:\n"
            "→ team_ndt_inverse('ndt_problems')\n"
            "Problem 8 (coil above crack), Problem 15 (slot), Problem 27 (deep flaws).\n"
            "Cross-reference: radia_mcp.ndt.eddy_current_NDT\n"
        ),
        "topology_optimization": (
            "Topology / shape optimization:\n"
            "→ team_ndt_inverse('inverse_opt')\n"
            "Problem 25 (Die Press) for topology opt.\n"
            "Problem 22 (SMES) for parametric shape opt.\n"
        ),
        "magnetic_levitation": (
            "Magnetic levitation device:\n"
            "→ team_force_motion('problem_28')\n"
            "Problem 28 — eddy current Lorentz levitation.\n"
        ),
    }
    return guidance.get(feature, (
        f"Unknown feature '{feature}'. Available:\n"
        "  maxwell_stress_force, pm_force, local_force_per_node,\n"
        "  nonlinear_magnetostatic, 3d_eddy_current, vector_hysteresis,\n"
        "  moving_conductor_eddy, nonlinear_transient, ndt_eddy_current,\n"
        "  topology_optimization, magnetic_levitation.\n"
        "For full catalog, see team_catalog('overview').\n"
    ))




register_status_tool(
    mcp,
    server_name='mcp-server-team-benchmark',
    description='TEAM Workshop benchmark problems reference layer (30 problems × physics class). ★ Lab core: 13, 20, 23, 32, 33b',
    subpackage='radia_mcp.team_benchmark',
    related_servers=["fem", "bem", "motor"],
)


def main():
    """Entry point for mcp-server-team-benchmark."""
    if "--selftest" in sys.argv:
        print("team_benchmark MCP server self-test:")
        print(f"  catalog: {len(get_catalog_knowledge('all'))} chars")
        print(f"  magnetostatic: {len(get_magnetostatic_knowledge('all'))} chars")
        print(f"  eddy_current: {len(get_eddy_current_knowledge('all'))} chars")
        print(f"  force_motion: {len(get_force_motion_knowledge('all'))} chars")
        print(f"  ndt_inverse: {len(get_ndt_inverse_knowledge('all'))} chars")
        print(f"  special: {len(get_special_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
