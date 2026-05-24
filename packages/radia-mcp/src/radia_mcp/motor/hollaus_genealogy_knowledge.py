"""Hollaus genealogy: visualization of the TU Wien MSFEM research lineage.

The Karl Hollaus / TU Wien MSFEM stack has 30+ papers (1992-2026) that
build on each other in a complex dependency graph.  This module renders
the genealogy as a text-format DAG with edges labeled by what each
paper contributes downstream.

The graph is meant for navigation: "given I'm reading paper X, which
papers does it build on, and which papers extend it?"

Usage:
    motor_hollaus_genealogy("dot")        # Graphviz DOT format
    motor_hollaus_genealogy("ascii")      # ASCII art tree
    motor_hollaus_genealogy("chrono")     # Chronological list
    motor_hollaus_genealogy("by_topic")   # Grouped by research thread
"""
from __future__ import annotations


# ============================================================
# Paper Database
# ============================================================

# Each entry: short_key -> (year, authors_short, title_short,
#                            parents=[keys], topics=[set], doi)
PAPERS = {
    # --- Foundational FE theory (pre-MSFEM, used as baseline) ---
    "hiptmair_xu_2007": (
        2007,
        "Hiptmair-Xu",
        "Nodal auxiliary space preconditioning H(curl)/H(div)",
        [],
        {"AMS", "preconditioning"},
        "10.1137/060660588",
    ),

    # --- Hollaus genealogy starts ~2012 ---
    "hollaus_2012_msfem_intro": (
        2012,
        "Hollaus-Schoberl",
        "MSFEM for laminated cores (first paper)",
        [],
        {"MSFEM", "2D-1D", "TEAM21"},
        None,
    ),
    "hannukainen_hollaus_schoberl_2014": (
        2014,
        "Hannukainen-Hollaus-Schoberl",
        "Two-scale homogenization (foundational nonlinear)",
        ["hollaus_2012_msfem_intro"],
        {"MSFEM", "two_scale", "nonlinear"},
        "10.1007/978-3-642-41834-8_13",
    ),
    "hollaus_msfem_2d1d_2014": (
        2014,
        "Hollaus",
        "MSFEM 2D-1D T-formulation (rotating machines)",
        ["hollaus_2012_msfem_intro"],
        {"MSFEM", "2D-1D", "T_formulation"},
        None,
    ),
    "hollaus_t_phi_open_2015": (
        2015,
        "Hollaus",
        "T-Phi-Phi MSFEM for open circuits (motors)",
        ["hollaus_msfem_2d1d_2014"],
        {"MSFEM", "T_Phi", "motors"},
        None,
    ),
    "hollaus_effective_material_2016": (
        2016,
        "Hollaus",
        "Effective material (EM) homogenization workhorse",
        ["hannukainen_hollaus_schoberl_2014"],
        {"EM", "effective_material", "homogenization"},
        None,
    ),
    "hollaus_error_estimator_2017": (
        2017,
        "Hollaus-Schoberl",
        "Equilibrated error estimator (Prager-Synge)",
        ["hollaus_effective_material_2016"],
        {"error_estimator", "MSFEM"},
        None,
    ),
    "schobinger_hollaus_tsukerman_2020": (
        2020,
        "Schobinger-Hollaus-Tsukerman",
        "Nonasymptotic homogenization (two mu_eff defs, resonance)",
        ["hannukainen_hollaus_schoberl_2014"],
        {"homogenization", "mu_eff", "resonance"},
        "10.1109/TMAG.2019.2949000",
    ),
    "schobinger_hollaus_2021_hierarchical": (
        2021,
        "Schobinger-Hollaus",
        "Hierarchical error estimator for 3D MSFEM",
        ["hollaus_error_estimator_2017"],
        {"error_estimator", "hierarchical", "3D"},
        None,
    ),
    "hollaus_schobinger_2024_msfem_mor_deim": (
        2024,
        "Hollaus-Schobinger",
        "MSFEM + MOR + structure-preserving DEIM",
        ["hollaus_effective_material_2016",
         "schobinger_hollaus_tsukerman_2020"],
        {"MSFEM", "MOR", "DEIM"},
        "10.1109/TMAG.2023.3334562",
    ),
    "hanser_2021_diploma": (
        2021,
        "Hanser (V.)",
        "Diploma thesis: MSFEM eddy currents + vector Preisach",
        ["hollaus_effective_material_2016",
         "hannukainen_hollaus_schoberl_2014"],
        {"MSFEM", "Preisach", "diploma_thesis", "TU_Wien"},
        "10.34726/hss.2021.87241",
    ),
    "hanser_msfem_preisach_2022": (
        2022,
        "Hanser-Schobinger-Hollaus",
        "Mixed MSFEM with T,Phi-Phi + vector hysteresis (COMPEL)",
        ["hanser_2021_diploma",
         "hollaus_t_phi_open_2015"],
        {"MSFEM", "T_Phi", "Preisach", "hysteresis"},
        "10.1108/COMPEL-07-2021-0240",
    ),
    "hanser_schobinger_hollaus_2025": (
        2025,
        "Hanser-Schobinger-Hollaus",
        "Schur-complement circuit coupling (300x-5000x speedup)",
        ["hollaus_schobinger_2024_msfem_mor_deim",
         "hanser_msfem_preisach_2022"],
        {"circuit_coupling", "Schur", "A_formulation"},
        "10.1108/COMPEL-04-2024-0163",
    ),
    "frljic_2026_perpendicular_flux": (
        2026,
        "Frljic-Trkulja-Zilic",
        "Anisotropic mu_eff for perpendicular flux (open-type cores)",
        ["hollaus_effective_material_2016",
         "schobinger_hollaus_tsukerman_2020"],
        {"EM", "anisotropic", "perpendicular_flux"},
        "10.1109/TMAG.2025.3543210",
    ),

    # --- Adjacent / referenced lineages ---
    "lange_henrotte_hameyer_2009": (
        2009,
        "Lange-Henrotte-Hameyer",
        "Temporary linearization for field-circuit coupling",
        [],
        {"FEMM", "linearization", "circuit"},
        "10.1109/TMAG.2009.2012585",
    ),
    "henrotte_hameyer_2004": (
        2004,
        "Henrotte-Hameyer",
        "EM force density (energy method)",
        [],
        {"force", "Henrotte"},
        "10.1109/TMAG.2004.825293",
    ),
    "carstensen_2007_phd": (
        2007,
        "Carstensen",
        "187p PhD: SRM winding eddy current losses",
        [],
        {"Carstensen", "AC_copper_loss", "SRM"},
        None,
    ),
}


# ============================================================
# Renderers
# ============================================================

def _render_dot() -> str:
    """Render as Graphviz DOT (for graphviz `dot -Tpng`)."""
    lines = ["digraph hollaus_genealogy {",
             "  rankdir=LR;",
             "  node [shape=box, fontname=\"Arial\", fontsize=10];",
             ""]
    # Group by year for left-to-right time axis
    for k, (year, auths, title, parents, _, doi) in sorted(
            PAPERS.items(), key=lambda kv: kv[1][0]):
        # Color by year
        col = "lightgray"
        if year >= 2024:
            col = "lightgreen"
        elif year >= 2020:
            col = "lightyellow"
        elif year >= 2014:
            col = "lightblue"
        label = f"{year}\\n{auths}\\n{title[:50]}"
        lines.append(f'  "{k}" [label="{label}", fillcolor="{col}", '
                     'style=filled];')
    lines.append("")
    for k, (_, _, _, parents, _, _) in PAPERS.items():
        for p in parents:
            if p in PAPERS:
                lines.append(f'  "{p}" -> "{k}";')
    lines.append("}")
    return "\n".join(lines)


def _render_ascii_tree() -> str:
    """ASCII art tree by lineage depth."""
    lines = ["Hollaus MSFEM Genealogy (text-tree view)",
             "=" * 60, ""]
    # Find roots (papers with no Hollaus parents)
    children = {k: [] for k in PAPERS}
    for k, (_, _, _, parents, _, _) in PAPERS.items():
        for p in parents:
            if p in children:
                children[p].append(k)
    roots = [k for k, (_, _, _, parents, _, _) in PAPERS.items()
             if not parents]

    def _show(key, depth=0):
        year, auths, title, _, topics, doi = PAPERS[key]
        prefix = "  " * depth + ("|-- " if depth > 0 else "")
        tag = ",".join(sorted(topics))[:35]
        out = [f"{prefix}[{year}] {auths}: {title[:55]}"]
        out.append(" " * (len(prefix) + 2) + f"({tag})")
        for c in sorted(children.get(key, []), key=lambda x: PAPERS[x][0]):
            out.extend(_show(c, depth + 1))
        return out

    for r in sorted(roots, key=lambda x: PAPERS[x][0]):
        lines.extend(_show(r))
        lines.append("")
    return "\n".join(lines)


def _render_chrono() -> str:
    """Chronological flat list."""
    lines = ["Hollaus papers (chronological)",
             "=" * 60, ""]
    by_year = sorted(PAPERS.items(), key=lambda kv: kv[1][0])
    for k, (year, auths, title, parents, topics, doi) in by_year:
        lines.append(f"{year}  {auths}")
        lines.append(f"      {title}")
        if doi:
            lines.append(f"      DOI: {doi}")
        if parents:
            lines.append(f"      builds on: {', '.join(parents)}")
        lines.append("")
    return "\n".join(lines)


def _render_by_topic() -> str:
    """Grouped by main topic/thread."""
    threads = {
        "Foundational MSFEM (2012-2014)": [
            "hollaus_2012_msfem_intro",
            "hannukainen_hollaus_schoberl_2014",
            "hollaus_msfem_2d1d_2014",
        ],
        "Open-circuit motors (T-Phi, 2015-)": [
            "hollaus_t_phi_open_2015",
        ],
        "Effective Material workhorse (2016-)": [
            "hollaus_effective_material_2016",
            "frljic_2026_perpendicular_flux",
        ],
        "Error estimators (2017, 2021)": [
            "hollaus_error_estimator_2017",
            "schobinger_hollaus_2021_hierarchical",
        ],
        "Nonasymptotic homogenization (2020)": [
            "schobinger_hollaus_tsukerman_2020",
        ],
        "MOR + DEIM + circuit coupling (2024-2025)": [
            "hollaus_schobinger_2024_msfem_mor_deim",
            "hanser_schobinger_hollaus_2025",
        ],
        "Hanser Diploma + Preisach hysteresis chain": [
            "hanser_2021_diploma",
            "hanser_msfem_preisach_2022",
        ],
        "Adjacent: FEMM/Henrotte/Carstensen": [
            "lange_henrotte_hameyer_2009",
            "henrotte_hameyer_2004",
            "carstensen_2007_phd",
        ],
    }
    lines = ["Hollaus genealogy grouped by research thread",
             "=" * 60, ""]
    for thread, keys in threads.items():
        lines.append(f"## {thread}")
        for k in keys:
            if k in PAPERS:
                y, auths, title, _, _, _ = PAPERS[k]
                lines.append(f"  [{y}] {auths}: {title}")
        lines.append("")
    return "\n".join(lines)


def get_hollaus_genealogy(view: str = "by_topic") -> str:
    """Render the genealogy in the requested view.

    Args:
        view: "dot" (Graphviz), "ascii" (tree), "chrono" (timeline),
              "by_topic" (research-thread groups), "all" (all four).
    """
    if view == "dot":
        return _render_dot()
    if view == "ascii":
        return _render_ascii_tree()
    if view == "chrono":
        return _render_chrono()
    if view == "by_topic":
        return _render_by_topic()
    if view == "all":
        return ("\n\n" + "=" * 60 + "\n").join([
            _render_by_topic(),
            _render_chrono(),
            _render_ascii_tree(),
            "## Graphviz DOT (paste into https://dreampuf.github.io/GraphvizOnline/)\n\n"
            + _render_dot(),
        ])
    return f"Unknown view: {view}.  Try 'by_topic', 'chrono', 'ascii', 'dot', or 'all'."
