"""
MOR (Model Order Reduction) MCP Server (radia_mcp.mor)

Knowledge layer for eddy-current FEM model order reduction, centered
on the **Cauer Ladder Network (CLN)** method — a LAB SPECIALTY of
菅原研 (Sugahara Lab) where Sugahara is co-author on the canonical
foundational papers (2018-2020).

Usage:
    mcp-server-mor              # Start MCP server (stdio)
    mcp-server-mor --selftest   # Self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .cln_knowledge import get_cln_documentation
from .systematic_knowledge import get_systematic_mor_knowledge

# CLN deep-dive modules from the 2026-05-25 absorption of
# public-safe curated corpus (~5238 lines / 57 topics across 5 modules).
# See CHANGELOG.md 0.74.0 entry.
from .cln_practice_knowledge import get_cln_practice_documentation
from .cln_multiport_knowledge import get_cln_multiport_documentation
from .cln_advanced_knowledge import get_cln_advanced_documentation
from .cln_specialty_knowledge import get_cln_specialty_documentation
from .cln_collab_knowledge import get_cln_collab_documentation

try:
    from .bibliography_index_knowledge import get_bibliography_index
except ImportError:
    def get_bibliography_index(query: str = "") -> str:
        return "MOR bibliography index not yet generated."


mcp = FastMCP("mcp-server-mor")


@mcp.tool()
def mor_cln(topic: str = "all") -> str:
    """Cauer Ladder Network (CLN) -- the lab-specialty MOR for eddy-current FEM.

    Topics:
      "all"
      "overview"         - CLN big picture, Sugahara Lab's role
      "recursion"        - Basic A-formulation recursion + Cauer circuit
      "multiple"         - Multiple expansion points (Kuriyama 2019)
      "nonlinear"        - Nonlinear ferromagnetic extension
      "applications"     - Industrial inductors, WPT, hybrid twin
      "pgd_positioning"  - CLN as a PGD instance (Köster-König-Bíró 2021,
                           IEEE TMag 57(6)); implications for the
                           base-agnostic Theorem 1 Schur augmentation
                           and for Multi-K / POD positioning
    """
    return get_cln_documentation(topic)


@mcp.tool()
def mor_systematic(topic: str = "mor_taxonomy") -> str:
    """
    Systematic MOR knowledge -- distilled from the deGruyter 3-volume
    MOR Handbook (1221 pages) + Kiss-Orosz 2024 Energies rotating-
    machines review + Ioan Vol3 Ch5 EM-specific MOR.

    Args:
        topic: One of:
            "mor_taxonomy"          - Three families: projection / system / data
            "projection_pod_rb"     - POD + Reduced Basis (Vol 2 Ch 2, 4)
            "projection_krylov"     - Krylov / Arnoldi / Lanczos (Vol 1 Ch 3)
            "pgd"                   - Proper Generalized Decomposition (Vol 2 Ch 3)
            "hyperreduction"        - DEIM / EIM for nonlinear (Vol 2 Ch 5)
            "system_theoretic_bt"   - Balanced Truncation, H2/H_inf (Vol 1 Ch 2)
            "data_driven_dmd_oi"    - DMD, OI, Loewner (Vol 1 Ch 6 + Vol 2 Ch 7)
            "parametric_pmor"       - Parametric MOR (Vol 2 Ch 1)
            "cln_sugahara"          - CLN positioned in MOR taxonomy
            "em_specific_ioan"      - Ioan Vol 3 Ch 5 (56p EM-specific)
            "rotating_machines_kiss_orosz_2024" - 2024 review for motors
            "software_lab"          - pyMOR / MORLab / lab tools
            "lab_recommendation"    - Decision guide for radia + NGSolve
            "all"                   - Everything
    """
    return get_systematic_mor_knowledge(topic)


@mcp.tool()
def mor_bibliography(query: str = "") -> str:
    """Search the MOR bibliography catalog (87 papers in lab library)."""
    return get_bibliography_index(query)


# ============================================================
# CLN deep-dive tools -- 5238 lines / 57 topics, absorbed
# 2026-05-25 from public-safe curated corpus (~500 .m / .docx / .pdf
# files across 16 topic folders + 6 root references). Group out by
# theme so a topic query maps to a tractable knowledge chunk.
# ============================================================

@mcp.tool()
def mor_cln_practice(topic: str = "all") -> str:
    """CLN MATLAB+COMSOL practice corpus -- foundations + 2020_11_04 lab
    practice + A-phi / 2D-rethink / Bessel root references.

    Source: public-safe curated corpus, 02_境界条件,
    09_級数展開, 2020_11_04_線形のCLNの練習} + A-phi.pdf + 2D-rethink
    docx + Bessel-inequality docx.

    Topics (12):
      "overview", "matlab_class" (CLN.m 71-line full class),
      "e0_mode" (HelmholtzEquation+withsol recursion idiom +
      Legendre analytical formulas to n=9), "h1_mode",
      "boundary_conditions" (Robin vs Infinite vs Kelvin -- the
      _NG.m failure mode!), "robin_bc", "series_expansion",
      "a_phi", "quasi_static_2d", "lab_lessons" (11 pitfalls),
      "key_files", "citations", "all".
    """
    return get_cln_practice_documentation(topic)


@mcp.tool()
def mor_cln_multiport(topic: str = "all") -> str:
    """CLN multi-port + multi-expansion-point + 3D extensions.

    Source: public-safe curated corpus, 04_展開点,
    10_タイル, 11_3次元}.

    Topics (10):
      "overview", "multiport_theory" (Z_ij matrix from scalar Z(f),
      Sugahara 2017 voltage-source / Matsuo 2018c matrix-form),
      "multiport_code" (FreeFEM++ Multi-turnLadderSeries.edp +
      MATLAB transformer/CLN_scalar.m quoted),
      "multi_expansion_theory" (Kuriyama 2019 K = C^T nu C + s_0 sigma
      operator shift, 4 variants A/T/3D/AK),
      "multi_expansion_code" (One Conductor_CLN_s0.edp quoted),
      "periodic_tiled" (tile graphical ladder pedagogy vs
      Bloch-Floquet periodic), "three_d" (A_phi_Gridap.jl HCurl/H1
      saddle-point, 4-stage VTK output), "lab_lessons" (10),
      "key_files", "citations", "all".
    """
    return get_cln_multiport_documentation(topic)


@mcp.tool()
def mor_cln_advanced(topic: str = "all") -> str:
    """CLN advanced physics extensions -- HF, nonlinear, circuit,
    Bloch, fixed-point (FP) CLN.

    Source: public-safe curated corpus, 06_非線形,
    07_電気力学, 14_ブロッホ方程式, 16_FP_CLN, 2020_12_07_高周波CLN
    の練習}.

    Topics (12):
      "overview", "high_frequency" (Helmholtz regime symbolic
      prototype, J_{n+1} uses epsilon not sigma),
      "high_frequency_code", "nonlinear" (4 generations 2017-2023,
      Tobita's jw method = FP-CLN predecessor), "nonlinear_code",
      "electrodynamic_circuit" (Shindo electromagnet shindo.edp,
      437 lines -- canonical CLN-as-SPICE-block pattern),
      "bloch_equation" (honestly thin: 1 design doc 2019),
      "fp_cln" (★ Fixed-Point CLN with constant virtual mu_FP
      linearising LHS, CEFC 2024 Sugahara-Tobita-Matsuo-Takahashi,
      MATLAB OO + NGSolve coupling), "fp_cln_code",
      "lab_lessons", "key_files", "citations", "all".
    """
    return get_cln_advanced_documentation(topic)


@mcp.tool()
def mor_cln_specialty(topic: str = "all") -> str:
    """CLN lab-signature techniques -- termination, Hiruma method,
    BEM+FEM coupling, Nagamine error theory.

    Source: public-safe curated corpus, 12_比留間法,
    13_BEM+FEM, 15_誤差論@長嶺氏}.

    Topics (11):
      "overview", "termination_theory" (frequency-dependent
      Z_term(s) vs open/short truncation; "more stages > better
      terminator" empirical finding 2016-12-23),
      "termination_code" (cf.m + termination.m quoted),
      "hiruma_method" (★ Shingo Hiruma, Hokkaido Igarashi Lab ->
      Kyoto Matsuo Lab: non-symmetric Lanczos producing Cauer ladder
      directly from algebraic (G+sC)x=b, unifies CLN with
      PVL/SyPVL/PRIMA), "hiruma_code" (hiruma_scalar.m +
      hiruma_matrix.m + hiruma3.m expansion variants),
      "bem_fem_coupled" (Kameari-CLN interior + BEM exterior;
      Dirichlet->Robin BC change makes lambda_n s-dependent;
      TSVD truncation O(M*N_m) -> K=5-15 ports),
      "bem_fem_code" (gen_svd_mode.m + MoM_A.m Bessel reference),
      "nagamine_error_theory" (★ Hideaki Nagamine, Kyoto Matsuo
      Lab: 4 PDFs culminating CEFC 2022; verified-arithmetic exact
      R_1..R_18 / L_1..L_17 for square section; mesh-adequacy rule
      delta_n >= 10*Delta_x from Foster cut-off frequency),
      "lab_lessons", "key_files", "citations", "all".
    """
    return get_cln_specialty_documentation(topic)


@mcp.tool()
def mor_cln_collab(topic: str = "all") -> str:
    """CLN external-collaborator work + Cauer-form conversion +
    application examples.

    Source: public-safe curated corpus,
    2022_10_09_遠藤@法政, 2023_08_08_松本@法政,
    2026_04_01_長方形CLN} + 2017_05_17_インバータ回路.docx.

    Topics (10):
      "overview", "cauer_i_vs_ii" (continued-fraction expansion of
      Z(s) around s=0 (CauerI = series R, parallel L) vs around
      s=infty (CauerII = series L, parallel R); two-matrix Lanczos
      in K-inner-product, N<=7 stability before loss of
      orthogonality), "cauer_form_conversion_code",
      "endo_hosei" (Endo @ Hosei Univ 2022: 2D vector-potential CLN
      on 4-square + 1-cylinder Al, COMSOL LiveLink-MATLAB sweep
      5 freq x 9 offset x 4 nStage), "matsumoto_hosei" (Matsumoto
      handoff 2023: 1 PDF + 2 mesh dumps, year-2 contract),
      "rectangular_cln_tbd" (2026_04_01: empty placeholder for
      rectangular-cross-section CLN, litz/motor slot),
      "inverter_application" (2017 design memo: CLN as SPICE
      subcircuit + H-bridge inverter, prototype for 07_電気力学),
      "lab_lessons", "key_files", "citations", "all".
    """
    return get_cln_collab_documentation(topic)


register_status_tool(
    mcp,
    server_name='mcp-server-mor',
    description='Model Order Reduction: PRIMA, Cauer Ladder Network, hyperreduction (DEIM)',
    subpackage='radia_mcp.mor',
    related_servers=["radia-ngsolve", "rna-mec"],
)


def main():
    if "--selftest" in sys.argv:
        print("MOR MCP server self-test:")
        print(f"  cln('all'): {len(mor_cln('all'))} chars")
        from .systematic_knowledge import SECTIONS as S
        for k in S:
            r = mor_systematic(k)
            assert len(r) > 500, f"{k} short ({len(r)})"
            print(f"  systematic({k!r}): {len(r)} chars")
        # CLN deep-dive modules (2026-05-25 absorption)
        cln_deep = [
            ("cln_practice", mor_cln_practice),
            ("cln_multiport", mor_cln_multiport),
            ("cln_advanced", mor_cln_advanced),
            ("cln_specialty", mor_cln_specialty),
            ("cln_collab", mor_cln_collab),
        ]
        for name, fn in cln_deep:
            r = fn("all")
            assert len(r) > 2000, f"{name}('all') short ({len(r)})"
            print(f"  {name}('all'): {len(r)} chars")
            r_ov = fn("overview")
            assert len(r_ov) > 200, f"{name}('overview') short"
        print("  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
