"""MCP Server: radia_mcp.magnetic_materials

Comprehensive magnetic material knowledge: hysteresis models,
iron loss, silicon steel database, permanent magnets,
demagnetization factor, Radia implementation status.

Usage:
    mcp-server-magnetic-materials              # stdio
    mcp-server-magnetic-materials --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .hysteresis_models_knowledge import get_hysteresis_models_knowledge
from .iron_loss_knowledge import get_iron_loss_knowledge
from .silicon_steel_knowledge import get_silicon_steel_knowledge
from .permanent_magnet_knowledge import get_permanent_magnet_knowledge
from .demagnetization_knowledge import get_demagnetization_knowledge
from .radia_status_knowledge import get_radia_status_knowledge


mcp = FastMCP("mcp-server-magnetic-materials")


@mcp.tool()
def magnetic_materials_hysteresis(topic: str = "lab_core") -> str:
    """
    Hysteresis model catalog & decision tree.

    ★ DEFAULT topic = 'lab_core' (B-input Stop based energy model),
    which is the Sugahara Lab's PRIMARY production hysteresis method.

    Covers 13+ models documented in W:/磁気特性/ヒステリシス/:
    Jiles-Atherton, Play, Stop, Energy-Based (Henrotte/Bergqvist),
    Preisach, Chua, Chan, Bouc-Wen, E&S, Lee, Potter-Schmulian, Zirka,
    LLG, plus adjacent (Cellular Automaton, Inverse Distribution,
    Thermal extension, Friction).

    Args:
        topic: One of:
            "lab_core"       - ★ Sugahara Lab CORE method: B-input Stop
                                + Energy formulation (DEFAULT, start here)
            "catalog"        - 13-model overview table + paper references
            "decision_tree"  - When to use which (decision tree)
            "jiles_atherton" - JA model formulation + parameters
            "play"           - Play model (B-input alternative)
            "energy_based"   - Energy-Based formulation theory
            "preisach"       - Preisach formalism (Everett function, FORC)
            "vector"         - Vector hysteresis (E&S, vector Preisach, RSST)
            "all"            - Everything
    """
    return get_hysteresis_models_knowledge(topic)


@mcp.tool()
def magnetic_materials_iron_loss(topic: str = "decision") -> str:
    """
    Iron loss models: Steinmetz family, Bertotti 3-term, Carstensen,
    non-sinusoidal corrections (MSE, iGSE).

    Args:
        topic: One of:
            "steinmetz_family" - Steinmetz + MSE + iGSE for non-sinusoid
            "bertotti"         - 3-term decomp + JIS steel constants
            "carstensen"       - AC copper + per-layer iron loss
            "waveform"         - PWM / DC bias corrections + decision tree
            "decision"         - Default = waveform topic
            "all"              - Everything
    """
    return get_iron_loss_knowledge(topic)


@mcp.tool()
def magnetic_materials_silicon_steel(topic: str = "grades") -> str:
    """
    JIS silicon steel grade database + processing/handling notes.

    Source: W:/磁気特性/00_教科書/電磁鋼板の磁気特性と取扱方法.pdf
    (94 MB, JFE Steel handbook).

    Args:
        topic: One of:
            "grades"   - JIS naming + material constants table
            "handling" - Processing effects on iron loss + Building Factor
            "all"      - Both
    """
    return get_silicon_steel_knowledge(topic)


@mcp.tool()
def magnetic_materials_permanent_magnet(topic: str = "families") -> str:
    """
    Permanent magnet datasheets: NdFeB, SmCo, Ferrite, AlNiCo
    + temperature derating + demag curves + BHmax + knee.

    Args:
        topic: One of:
            "families"    - Family comparison + selection guide
            "demag_curve" - B-H vs J-H, BHmax, knee, operating point
            "all"         - Both
    """
    return get_permanent_magnet_knowledge(topic)


@mcp.tool()
def magnetic_materials_demagnetization(topic: str = "overview") -> str:
    """
    Demagnetization factor N (反磁場係数): Osborn 1945 closed-form
    for ellipsoids + non-ellipsoid approximations.

    Args:
        topic: One of:
            "overview"      - Why N matters, sum rule, FE caveats
            "osborn"        - Ellipsoid closed-form + reference table
            "non_ellipsoid" - Cylinder, prism, cube approximations
            "all"           - Everything
    """
    return get_demagnetization_knowledge(topic)


@mcp.tool()
def magnetic_materials_radia_status(topic: str = "overview") -> str:
    """
    Radia magnetic material implementation status (Mat classes).

    Maps the hysteresis_models / iron_loss / silicon_steel /
    permanent_magnet theoretical knowledge into running Radia code.

    Args:
        topic: One of:
            "overview"      - Mat class table + .hys format + state mgmt
            "how_to_choose" - Decision tree for picking a Mat class
            "todo_list"     - Research-grade additions still missing
            "all"           - Everything
    """
    return get_radia_status_knowledge(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def new_magnetic_material_simulation(material_type: str) -> str:
    """Set up a new simulation requiring magnetic material data."""
    guidance = {
        "soft_iron": (
            "Soft iron (silicon steel) simulation:\n"
            "1. Identify JIS grade (e.g., 35JN230)\n"
            "   → magnetic_materials_silicon_steel('grades')\n"
            "2. DC or AC? Linear or nonlinear?\n"
            "   → magnetic_materials_radia_status('how_to_choose')\n"
            "3. Iron loss method (Steinmetz / Bertotti / iGSE)?\n"
            "   → magnetic_materials_iron_loss('waveform')\n"
            "4. Laminated stack? Use Hollaus MSFEM\n"
            "   → motor_hollaus_eddy('effective_material')\n"
        ),
        "permanent_magnet": (
            "Permanent magnet simulation:\n"
            "1. Family selection (NdFeB / SmCo / Ferrite / AlNiCo)?\n"
            "   → magnetic_materials_permanent_magnet('families')\n"
            "2. Demag curve operating point + safety margin?\n"
            "   → magnetic_materials_permanent_magnet('demag_curve')\n"
            "3. Temperature derating?\n"
            "   → Use NdFeB-SH / EH grades for >100°C\n"
            "4. Force on PM?\n"
            "   → calc_em_force.py --pm-magnetization\n"
            "   → motor_em_force_extras('permanent_magnet_force')\n"
        ),
        "hysteresis": (
            "Hysteresis-aware simulation:\n"
            "1. Which model? (Play / Energy-Based / JA / Preisach)\n"
            "   → magnetic_materials_hysteresis('decision_tree')\n"
            "2. Calibration data available? (major loop / FORC)\n"
            "   → magnetic_materials_hysteresis('catalog')\n"
            "3. .hys file format for Radia Play/Energy:\n"
            "   → magnetic_materials_radia_status('overview')\n"
            "4. Transient solver:\n"
            "   → calc_motor_transient.py + MatPlayHysteresis state mgmt\n"
        ),
        "demag_sizing": (
            "Demagnetization-factor-based sizing:\n"
            "1. Approximate body as ellipsoid? Use Osborn 1945\n"
            "   → magnetic_materials_demagnetization('osborn')\n"
            "2. Non-ellipsoid (cylinder, prism)? Use tabulated approximation\n"
            "   → magnetic_materials_demagnetization('non_ellipsoid')\n"
            "3. Final design: let FE handle demag, use Osborn as cross-check\n"
        ),
    }
    return guidance.get(material_type, (
        f"Unknown material_type '{material_type}'. Available: soft_iron, "
        "permanent_magnet, hysteresis, demag_sizing.\n"
        "For overview, see magnetic_materials_radia_status('overview').\n"
    ))




register_status_tool(
    mcp,
    server_name='mcp-server-magnetic-materials',
    description='Magnetic materials: hysteresis (Play/Energy lab core), iron loss (Bertotti/Steinmetz/iGSE), JIS silicon steel, PM datasheets, Osborn...',
    subpackage='radia_mcp.magnetic_materials',
    related_servers=["ih", "motor", "electromagnet"],
)


def main():
    """Entry point for mcp-server-magnetic-materials."""
    if "--selftest" in sys.argv:
        print("magnetic_materials MCP server self-test:")
        print(f"  hysteresis: {len(get_hysteresis_models_knowledge('all'))} chars")
        print(f"  iron_loss: {len(get_iron_loss_knowledge('all'))} chars")
        print(f"  silicon_steel: {len(get_silicon_steel_knowledge('all'))} chars")
        print(f"  permanent_magnet: {len(get_permanent_magnet_knowledge('all'))} chars")
        print(f"  demagnetization: {len(get_demagnetization_knowledge('all'))} chars")
        print(f"  radia_status: {len(get_radia_status_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
