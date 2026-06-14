"""MCP Server: radia_mcp.pcb

Wireless Power Transfer knowledge layer for EV / robot / bearingless
motor / dynamic / FOD applications.

Covers: WPT landscape + decision tree, coil design + compensation
topology (SS/LCC/LCL), efficiency (Q+kQ) + safety + IEC/SAE standards,
Foreign Object Detection (★ lab core), applications (dynamic EV, robot,
bearingless motor), alternatives (capacitive, microwave/rectenna,
metamaterial).

Distilled from 269 files / 2.1 GB in
W:/03_文献・論文/00_電磁界解析/99_アプリケーション/05_ワイヤレス給電/

Cross-references:
- `radia_mcp.peec` — PEEC for coil L, M, k computation
- `radia_mcp.ih` — analogous SIBC / eddy current physics
- `radia_mcp.bem` — open-boundary BEM for coil field
- `radia_mcp.motor` — bearingless motor design

Usage:
    mcp-server-pcb              # stdio
    mcp-server-pcb --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .overview_knowledge import get_overview_knowledge
from .coil_compensation_knowledge import get_coil_compensation_knowledge
from .efficiency_safety_knowledge import get_efficiency_safety_knowledge
from .fod_knowledge import get_fod_knowledge
from .applications_knowledge import get_applications_knowledge
from .alternatives_knowledge import get_alternatives_knowledge


mcp = FastMCP("mcp-server-pcb")


@mcp.tool()
def pcb_overview(topic: str = "decision_tree") -> str:
    """
    WPT landscape: regimes, decision tree, lab focus.

    Args:
        topic: One of:
            "landscape"        - WPT regimes (inductive/resonant/microwave)
            "decision_tree"    - Which topic for which design Q (DEFAULT)
            "history"          - WPT history + lab lineage
            "all"              - Everything
    """
    return get_overview_knowledge(topic)


@mcp.tool()
def pcb_coil_compensation(topic: str = "coil_design") -> str:
    """
    Coil design + compensation topology + resonance matching.

    Args:
        topic: One of:
            "coil_design"            - Coil geometry, Litz, ferrite (DEFAULT)
            "compensation_topology"  - SS/LCC/LCL/CLLLC choice
            "resonance_matching"     - Resonant condition, Z reflection
            "all"                    - Everything
    """
    return get_coil_compensation_knowledge(topic)


@mcp.tool()
def pcb_efficiency_safety(topic: str = "q_factor_kQ") -> str:
    """
    Efficiency (Q, k, kQ) + safety + IEC/SAE standards.

    Args:
        topic: One of:
            "q_factor_kQ"             - Q, k, kQ product, η_max (DEFAULT)
            "coupling_k"              - k computation + misalignment
            "safety_human_exposure"   - ICNIRP, pacemaker, implant safety
            "standards_IEC_SAE"       - IEC 61980, SAE J2954, Qi
            "all"                     - Everything
    """
    return get_efficiency_safety_knowledge(topic)


@mcp.tool()
def pcb_fod(topic: str = "overview") -> str:
    """
    ★ Foreign Object Detection (FOD) — lab core research.

    Args:
        topic: One of:
            "overview"               - 6 FOD method families (DEFAULT)
            "sensorless"             - V/I/freq monitoring methods
            "sensor_based"           - Search coil array + camera + radar
            "permeability_patents"   - Permeability + patents + IEC 61980-3
            "all"                    - Everything
    """
    return get_fod_knowledge(topic)


@mcp.tool()
def pcb_applications(topic: str = "dynamic_ev") -> str:
    """
    WPT applications: dynamic EV, robot, bearingless motor.

    Args:
        topic: One of:
            "dynamic_ev"           - EV running, segmented primary (DEFAULT)
            "robot_bearingless"    - ★ Lab specialty
            "pcb_lab_lineage"      - Sugahara lab paper lineage
            "all"                  - Everything
    """
    return get_applications_knowledge(topic)


@mcp.tool()
def pcb_alternatives(topic: str = "capacitive") -> str:
    """
    Alternative WPT: capacitive, microwave/rectenna, metamaterial.

    Args:
        topic: One of:
            "capacitive"            - E-field WPT (DEFAULT)
            "microwave_rectenna"    - Microwave / rectenna / drone
            "metamaterial"          - Metamaterial-enhanced WPT
            "all"                   - Everything
    """
    return get_alternatives_knowledge(topic)


# ============================================================
# MCP Prompts
# ============================================================

@mcp.prompt()
def design_pcb_system(scenario: str) -> str:
    """Walk through WPT system design for a scenario."""
    guidance = {
        "ev_stationary_85khz": (
            "EV stationary WPT @ 85 kHz (SAE J2954 compliant):\n"
            "1. Coil: DD or DDQ for misalignment tolerance\n"
            "   → pcb_coil_compensation('coil_design')\n"
            "2. Compensation: LCC-LCC topology (constant-I source)\n"
            "   → pcb_coil_compensation('compensation_topology')\n"
            "3. Target kQ > 10 for >90% efficiency\n"
            "   → pcb_efficiency_safety('q_factor_kQ')\n"
            "4. FOD: required per IEC 61980-3\n"
            "   → pcb_fod('overview')\n"
            "5. Standards: SAE J2954 + IEC 61980-1\n"
            "   → pcb_efficiency_safety('standards_IEC_SAE')\n"
            "6. Code: PEEC for L/M/k, Radia MMM for ferrite plate\n"
            "   → radia_mcp.peec\n"
        ),
        "robot_multi_joint": (
            "Multi-joint arm robot WPT:\n"
            "1. Chain of LCC-LCC stages, one per rotation axis\n"
            "   → pcb_applications('robot_bearingless')\n"
            "2. Use ディスクリピーター (disc repeater) for compact size\n"
            "3. Litz wire for high Q\n"
            "   → pcb_coil_compensation('coil_design')\n"
            "4. Cross-reference lab papers (WPT2016-23, 26, 39)\n"
            "5. Code: PEEC + RNA for joint-by-joint analysis\n"
            "   → radia_mcp.peec, radia_mcp.bem.mmm_msc.rna_mmm\n"
        ),
        "bearingless_motor": (
            "Bearingless motor + WPT integration (★ lab specialty):\n"
            "1. Stator: stationary coil for both levitation + power\n"
            "2. Rotor: receiver coil + rectifier + load\n"
            "   → pcb_applications('robot_bearingless')\n"
            "3. LCC-S compensation (constant-I to rotor)\n"
            "   → pcb_coil_compensation('compensation_topology')\n"
            "4. Cross-reference: 整流コイルを用いたベアリングレスモータ papers\n"
            "5. Coupling with motor design\n"
            "   → radia_mcp.motor\n"
        ),
        "fod_design": (
            "FOD subsystem design:\n"
            "1. Pick detection family by cost / sensitivity\n"
            "   → pcb_fod('overview')\n"
            "2. For EV: sensorless freq sweep + search coil array (dual-use)\n"
            "   → pcb_fod('sensorless')\n"
            "   → pcb_fod('sensor_based')\n"
            "3. For living object: capacitive sensor\n"
            "4. Comply with IEC 61980-3 test items\n"
            "   → pcb_fod('permeability_patents')\n"
            "5. Lab actively researching multi-method fusion\n"
        ),
        "dynamic_charging": (
            "Dynamic EV charging (running):\n"
            "1. Segmented primary coil array\n"
            "   → pcb_applications('dynamic_ev')\n"
            "2. LCC-LCC (k-insensitive)\n"
            "3. Position detection via dual-use search coil\n"
            "   → pcb_fod('sensor_based')\n"
            "4. Feed-forward Q control for transitions\n"
            "5. Code: PEEC for segmented primary, transient sim\n"
        ),
        "drone_microwave": (
            "Drone microwave WPT (NOT lab core but reference):\n"
            "1. 28 GHz mmWave with rectenna\n"
            "   → pcb_alternatives('microwave_rectenna')\n"
            "2. Pencil beam with steering\n"
            "3. Rectenna efficiency: 70-80% small power, 30-50% high power\n"
            "4. Cross-link: ngsolve.bem for high-freq full-wave\n"
        ),
    }
    return guidance.get(scenario, (
        f"Unknown scenario '{scenario}'. Available:\n"
        "  ev_stationary_85khz, robot_multi_joint, bearingless_motor,\n"
        "  fod_design, dynamic_charging, drone_microwave.\n"
        "For overview, see pcb_overview('decision_tree').\n"
    ))




register_status_tool(
    mcp,
    server_name='mcp-server-pcb',
    description='Wireless Power Transfer: coil + compensation (SS/LCC/LCL), efficiency, IEC 61980 / SAE J2954, FOD, dynamic EV / robot / bearingless...',
    subpackage='radia_mcp.pcb',
    related_servers=["peec", "litz-transmission", "maglev"],
)


def main():
    """Entry point for mcp-server-pcb."""
    if "--selftest" in sys.argv:
        print("wpt MCP server self-test:")
        print(f"  overview: {len(get_overview_knowledge('all'))} chars")
        print(f"  coil_compensation: {len(get_coil_compensation_knowledge('all'))} chars")
        print(f"  efficiency_safety: {len(get_efficiency_safety_knowledge('all'))} chars")
        print(f"  fod: {len(get_fod_knowledge('all'))} chars")
        print(f"  applications: {len(get_applications_knowledge('all'))} chars")
        print(f"  alternatives: {len(get_alternatives_knowledge('all'))} chars")
        print("OK")
        return
    mcp.run()


if __name__ == "__main__":
    main()
