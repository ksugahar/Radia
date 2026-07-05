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

    Covers 13+ models documented in public-safe curated corpus
    30_磁気特性/02_ヒステリシスモデル/ (1184 files):
    Jiles-Atherton, Play, Stop, Energy-Based (Bergqvist/Henrotte/
    Francois-Lavet/Jacques/Egger lineage), Preisach, Chua, Chan,
    Bouc-Wen, E&S, Lee, Potter-Schmulian, Zirka, LLG, plus adjacent
    (Cellular Automaton, Inverse Distribution, Thermal extension,
    Friction).

    Args:
        topic: One of:
            "lab_core"       - ★ Sugahara Lab CORE method: B-input Stop
                                + Energy formulation (DEFAULT, start here)
            "catalog"        - 13-model overview table + paper references
            "decision_tree"  - When to use which (decision tree)
            "jiles_atherton" - JA model formulation + parameters
            "play"           - Play model (B-input alternative)
            "energy_based"   - Energy-Based formulation theory
                                (Bergqvist/Henrotte/Francois-Lavet/Jacques)
            "efficient_eval" - Egger 2025 regularization + Schur complement
                                → O(K) inverse operator (the computational
                                recipe behind MatEnergyHysteresis)
            "preisach"       - Preisach formalism (Everett function, FORC)
            "vector"         - Vector hysteresis (E&S, vector Preisach, RSST)
            "input_dependent_shape" - Matsuo trilogy: w(B)·g_0 product form
                                       to remove "equal vertical chords"
                                       restriction; required for silicon
                                       steel; 5 identification methods
            "lamination_homogenization" - Henrotte-Steentjes-Hameyer-Geuzaine
                                          2015 IET 2-step approach:
                                          1D mesoscale → algebraic H(B,Ḃ,p_k)
                                          for macroscale motor FE
            "cellular_automaton" - Saito 兆古 (Hosei) Preisach ≡ 2D CA;
                                    Bitter-method domain visualization;
                                    1/f / Barkhausen connection
            "mathematical_foundations" - Krasnoselskii-Pokrovskii 1989 +
                                          Visintin 1994 (Applied Math
                                          Sci 111): formal play/stop
                                          variational inequalities,
                                          Hilpert's inequality,
                                          Prandtl-Ishlinskii equivalence,
                                          Preisach homogenization theorem
            "rheological_models" - Krasnoselskii §39 mechanical spring
                                    + friction analogy that motivates
                                    Bergqvist-Henrotte energy hysteresis;
                                    graph representation; transducers M/W
            "lab_core_production" - ★★ Production deep dive for the lab
                                    本命: Henrotte 2006 (h = h_r + h_i,
                                    elementary 4-param + N+1 fraction
                                    combined) + Francois-Lavet 2013
                                    (variational Ω = u - h·J + χ|J-J_p|,
                                    atanh saturation, M250-50A concrete
                                    params, Picard FE algorithm) +
                                    Egger 2025 Schur inverse.  Use this
                                    when implementing / calibrating.
            "stop_vs_play_silicon_steel" - Matsuo 2003 IEEE TMAG 39:1361
                                    silicon steel violates "equal vertical
                                    chords" so pure stop fails on
                                    asymmetric dc-biased loops.  Play+stop
                                    combination model C_i alternating
                                    corrections recovers measured loops
                                    precisely.  Lab production path:
                                    Energy-based > combination >
                                    input-dependent shape.
            "forward_inverse_variational" - Egger 2025 IEEE TMAG 61:7300207
                                    Rigorous theoretical foundation for
                                    Energy-based forward operator B=∂_H w*
                                    (scalar potential FE) AND inverse
                                    H=∂_B w (vector potential FE) via
                                    convex duality.  4 Assertions
                                    proven; TEAM 32 benchmark validated;
                                    O(N_χ) forward vs O(N_χ²) inverse
                                    (without Schur trick).
            "tp_eec_steady_state" - Kitao 2012 IEEE TMAG 48:3375 Novel
                                    Time-Periodic Explicit Error
                                    Correction for play model.  Correct
                                    UNKNOWNS AND HYSTERON STATES to
                                    recover symmetric steady state after
                                    inrush current.  Drastic convergence
                                    speedup for long transients with
                                    dc bias.  Ring core circuit equations.
            "h_input_b_input_combination" - Ito 2013 IEEE TMAG 49:1985
                                    H-input AND B-input play models BOTH
                                    needed for ferrite equivalent circuit
                                    (NiZn ring core 100kHz-2MHz).
                                    Decompose: H_AC=h_fast+R_1·dB/dt,
                                    B_DC=b_slow+b_fast.  Newton-Raphson
                                    inversion between H and B play.
            "bobbio_1997_unification" - Bobbio 1997 IEEE TMAG 33:4417
                                    Foundational PTM/STM unification.
                                    Series-PTM ≡ Parallel-PTM ≡
                                    Bobbio-Marrucci 1993.  Preisach ⊃
                                    PIM ⊃ PTM/STM hierarchy.  Numerical
                                    PTM ≈ STM^(-1) inversion.  Vitrovac
                                    7505Z + ferroresonance validation.
            "asp_model"      - Lee 2014 IEEE TMAG 50:7300104 Asymmetric
                                    transition Probability model.
                                    Function-based (NO hysterons), BH
                                    curve = FORC compressed along B.
                                    Beats 60-hysteron play on minor
                                    loops.  For MRI iron cores, memory
                                    motors.  ("ASPモデル" in Matsuo 2014.)
            "zirka_hdhm"     - Zirka 2004 IEEE TMAG 40:390 congruency-
                                    based History-Dependent Hysteresis
                                    Model.  Transplantation method:
                                    higher-order curves from internal
                                    FORC segments.  Negative-slope-free
                                    by construction.  Madelung's 5 rules.
                                    For PWM transients; HIHM (history-
                                    independent) 18% loss error on PWM.
            "chua_model"     - Chua-Bass 1972 dynamic circuit-oriented
                                    ODE dy/dt = w(dx/dt)·h(y)·g(x-f(y)).
                                    w-weight controls freq behaviour
                                    (loop widening type-c / narrowing
                                    type-d for fluorescent lamps).
                                    SPICE/ferroresonance; Saito Hosei
                                    lineage (= CLN MOR group).
            "other_models"   - Niche: Chan (SPICE transformer, 3-param
                                    hyperbolic branches), Bouc-Wen
                                    (structural/piezo ODE), Potter-
                                    Schmulian (recording media), LLG
                                    (nm micromagnetics ground-truth).
            "jacques_monograph" - ★★ Jacques 2018 ULiege PhD (254p),
                                    THE canonical Energy-Based reference.
                                    Full Clausius-Duhem thermodynamic
                                    derivation; single/multi-cell +
                                    homogenization; THREE discrete FE
                                    implementations (Vector Play /
                                    Variational / Angle-searching) with
                                    numerical pitfalls; Newton-Raphson
                                    inversion w/ analytical permeability
                                    tensor.  The WHY behind lab_core.
            "peeling_identification" - Analytic 剥ぎ取り shape-function
                                    identification of a B-input play/stop
                                    model from the descending FORC ladder
                                    by a difference scheme (Everett-density
                                    mu); NO least-squares (per-hysteron LSQ
                                    on the play variable is ill-conditioned
                                    because all p_k share B).  Sugahara/
                                    Ahagon playmodel.py shapeFunction +
                                    prototype_step0-4 calibration pipeline.
            "vector_play_model" - Canonical Sugahara/Ahagon B-input vector
                                    play model H=sum Hfunc_i(p_i) with play
                                    operator p_i=clip(p_i,B-zeta_i,B+zeta_i),
                                    equidistant zeta, Hfunc via peeling.
                                    Congruency (matches H-axis-congruent
                                    silicon steel); rotational BQM variants;
                                    complementary to the energy STOP (p=B-s).
            "all"            - Everything (27 topics)
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
            "minor_loops_dc"   - Taitoda 2015 DC-biased k_h(B_m, B_max)
                                  + apparent frequency for general
                                  hysteresis with minor loops
                                  (17% → 6% error vs conventional)
            "decision"         - Default = waveform topic
            "all"              - Everything
    """
    return get_iron_loss_knowledge(topic)


@mcp.tool()
def magnetic_materials_silicon_steel(topic: str = "grades") -> str:
    """
    JIS silicon steel grade database + processing/handling notes.

    Source: public-safe curated corpus
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
