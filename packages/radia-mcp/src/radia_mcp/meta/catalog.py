"""Authoritative catalog of all radia_mcp servers.

Source of truth: this dict. .mcp.json (in the Radia monorepo root)
should mirror it.

The catalog includes the `meta` server itself (self-referential entry)
so that `radia_mcp_overview` reports a complete picture and
`radia_mcp_by_tag("meta")` finds the entry-point.

Each entry:
    "<short_name>": {
        "subpackage": "radia_mcp.<name>",
        "entry_point": "mcp-server-<name>",
        "description": str,
        "primary_tools": [str, ...],
        "related": [other short names that pair well],
        "audit_command": optional heavy validation / repo-audit command,
        "selftest_command": auto-added from entry_point by list/get helpers,
        "tags": [str, ...]   for filtered queries
    }
"""

from __future__ import annotations
from typing import Any

# Order roughly follows logical groupings (general → application)
CATALOG: dict[str, dict[str, Any]] = {
    # ============================================================
    # CAD / mesh authoring
    # ============================================================
    "cubit": {
        "subpackage": "radia_mcp.cubit",
        "entry_point": "mcp-server-cubit",
        "description": "Cubit mesh scripting, hex/tet workflow, export formats",
        "primary_tools": ["cubit_exec", "cubit_mesh_auto", "cubit_docs",
                          "cubit_audit_summary"],
        "related": ["build123d", "gmsh", "radia-ngsolve"],
        "audit_command": "mcp-server-cubit --selftest --audit-examples",
        "tags": ["cad", "mesh"],
    },
    "build123d": {
        "subpackage": "radia_mcp.build123d",
        "entry_point": "mcp-server-build123d",
        "description": "build123d STEP authoring (CAD-as-code) + Cubit interop",
        "primary_tools": ["build123d_api", "execute_build123d",
                            "build123d_to_cubit_hex"],
        "related": ["cubit", "gmsh", "radia-ngsolve"],
        "tags": ["cad"],
    },
    "gmsh": {
        "subpackage": "radia_mcp.gmsh",
        "entry_point": "mcp-server-gmsh",
        "description": "GMSH MSH v4.1 inspect/validate/convert/write_node_data",
        "primary_tools": ["gmsh_usage", "gmsh_reference", "gmsh_audit_summary",
                          "gmsh_numsubedges_remediation_plan",
                          "gmsh_mesh_generation_remediation_plan"],
        "related": ["build123d", "cubit", "radia-ngsolve"],
        "audit_command": "mcp-server-gmsh --selftest --audit-examples",
        "tags": ["mesh"],
    },
    # ============================================================
    # FEM / BEM / specialty solvers
    # ============================================================
    "radia-ngsolve": {
        "subpackage": "radia_mcp.radia_ngsolve",
        "entry_point": "mcp-server-radia-ngsolve",
        "description": "Radia + NGSolve: Kelvin / sparsesolv / CLN / PEEC / "
                       "analytical formulas / lint",
        "primary_tools": ["kelvin_transformation", "ngsolve_usage",
                            "analytical_formulas", "peec_inductance"],
        "related": [
            "bem",
            "build123d",
            "cubit",
            "fem",
            "gmsh",
            "mathematica",
            "matrix-solvers",
            "mor",
            "peec",
            "radia-streamfunction",
        ],
        "tags": ["fem", "solver"],
    },
    "radia-streamfunction": {
        "subpackage": "radia_mcp.streamfunction",
        "entry_point": "mcp-server-radia-streamfunction",
        "description": "Stream-function coil design: (ACA+)+TSVD least-norm, "
                       "FE-direct psi, regularisation / folded-Tikhonov Pareto "
                       "front, single-stroke chain, sheet-metal levers",
        "primary_tools": ["streamfunction"],
        "related": [
            "radia-ngsolve",
            "bayesian-opt",
            "nmr-mri",
            "electromagnet",
            "peec",
        ],
        "tags": ["optimization", "solver"],
    },
    "fem": {
        "subpackage": "radia_mcp.fem",
        "entry_point": "mcp-server-fem",
        "description": "FEM formulations theory layer (A-Omega / T-Omega / H / "
                       "Reduced / Darwin, edge / HO / XFEM / IGA / DG, gauging "
                       "+ Kelvin, MSFEM, Schur circuit coupling, NGSolve hierarchical)",
        "primary_tools": ["fem_usage"],
        "related": [
            "bem",
            "differential-forms",
            "gnn",
            "matrix-solvers",
            "pinn",
            "radia-ngsolve",
            "team-benchmark",
        ],
        "tags": ["fem", "theory"],
    },
    "bem": {
        "subpackage": "radia_mcp.bem",
        "entry_point": "mcp-server-bem",
        "description": "MoM/BEM theory: RWG, EFIE/MFIE/CFIE/PMCHWT, "
                       "Loop-Star, Calderon, Radia MMM/MSC, HACApK, FEM-BEM",
        "primary_tools": ["bem_usage"],
        "related": [
            "fem",
            "peec",
            "radia-ngsolve",
            "team-benchmark",
        ],
        "tags": ["fem", "theory"],
    },
    "matrix-solvers": {
        "subpackage": "radia_mcp.matrix_solvers",
        "entry_point": "mcp-server-matrix-solvers",
        "description": "Sparse solver theory + decision tree: Krylov "
                       "(CG/BiCGSTAB/GMRES/COCG/COCR/IDR), preconditioners "
                       "(AMG, Hiptmair-Xu AMS), Biro-Preis A-V, tree-cotree",
        "primary_tools": ["matrix_solvers_usage"],
        "related": ["radia-ngsolve", "fem"],
        "tags": ["solver", "theory"],
    },
    "mor": {
        "subpackage": "radia_mcp.mor",
        "entry_point": "mcp-server-mor",
        "description": "Model Order Reduction: PRIMA, Cauer Ladder Network, "
                       "hyperreduction (DEIM)",
        "primary_tools": ["mor_usage"],
        "related": ["data-assimilation", "radia-ngsolve", "rna-mec"],
        "tags": ["solver"],
    },
    # ============================================================
    # Physics applications
    # ============================================================
    "ih": {
        "subpackage": "radia_mcp.ih",
        "entry_point": "mcp-server-ih",
        "description": "Induction heating: SIBC, ESIM, Karl iteration, "
                       "workpiece coupling",
        "primary_tools": ["induction_heating", "ih_sibc", "ih_esim"],
        "related": [
            "litz-transmission",
            "magnetic-materials",
            "ndt",
            "panel-review",
            "peec",
            "rna-mec",
        ],
        "tags": ["application"],
    },
    "peec": {
        "subpackage": "radia_mcp.peec",
        "entry_point": "mcp-server-peec",
        "description": "PEEC filament/panel, FastHenry parser, HOIBC, "
                       "Carstensen AC copper loss",
        "primary_tools": ["peec_usage", "peec_hoibc", "peec_carstensen_ac_loss"],
        "related": [
            "bem",
            "ih",
            "litz-transmission",
            "radia-ngsolve",
            "radia-streamfunction",
            "pcb",
        ],
        "tags": ["application"],
    },
    "electromagnet": {
        "subpackage": "radia_mcp.electromagnet",
        "entry_point": "mcp-server-electromagnet",
        "description": "Accelerator electromagnet: CoilBuilder, Hantila, "
                       "Play/Energy hysteresis",
        "primary_tools": ["electromagnet_usage"],
        "related": [
            "accelerator",
            "fusion-reactor",
            "magnetic-materials",
            "motor",
            "nmr-mri",
            "panel-review",
            "radia-streamfunction",
        ],
        "tags": ["application"],
    },
    "motor": {
        "subpackage": "radia_mcp.motor",
        "entry_point": "mcp-server-motor",
        "description": "Motor analysis: ONELAB transient, Hollaus effective "
                       "material (lamination), Wakao autoencoder topology, "
                       "Kaimori-Mifune Darwin TD",
        "primary_tools": ["motor_usage"],
        "related": [
            "electromagnet",
            "maglev",
            "magnetic-materials",
            "panel-review",
            "team-benchmark",
            "topology-optimization",
        ],
        "tags": ["application"],
    },
    "accelerator": {
        "subpackage": "radia_mcp.accelerator",
        "entry_point": "mcp-server-accelerator",
        "description": "Accelerator physics: beam optics, dipole/quad/sext "
                       "magnets, undulator/wiggler",
        "primary_tools": ["accelerator_usage"],
        "related": ["electromagnet", "fusion-reactor", "nmr-mri"],
        "tags": ["application"],
    },
    "fusion-reactor": {
        "subpackage": "radia_mcp.fusion_reactor",
        "entry_point": "mcp-server-fusion-reactor",
        "description": "Fusion reactor magnets: tokamak ITER + stellarator "
                       "LHD/W7-X/heliotron lineage",
        "primary_tools": ["fusion_reactor"],
        "related": ["accelerator", "data-assimilation", "electromagnet"],
        "tags": ["application"],
    },
    # ============================================================
    # Materials / hysteresis / losses
    # ============================================================
    "magnetic-materials": {
        "subpackage": "radia_mcp.magnetic_materials",
        "entry_point": "mcp-server-magnetic-materials",
        "description": "Magnetic materials: hysteresis (Play/Energy lab core), "
                       "iron loss (Bertotti/Steinmetz/iGSE), JIS silicon steel, "
                       "PM datasheets, Osborn demag factor",
        "primary_tools": ["magnetic_materials_usage"],
        "related": [
            "electromagnet",
            "ih",
            "motor",
            "ndt",
            "rna-mec",
        ],
        "tags": ["application"],
    },
    "litz-transmission": {
        "subpackage": "radia_mcp.litz_transmission",
        "entry_point": "mcp-server-litz-transmission",
        "description": "Litz wire AC loss (Dowell, homogenization, magnetic-"
                       "plated wire) + multiconductor transmission line theory",
        "primary_tools": ["litz_transmission"],
        "related": [
            "ih",
            "metamaterial",
            "peec",
            "pcb",
        ],
        "tags": ["application"],
    },
    "rna-mec": {
        "subpackage": "radia_mcp.rna_mec",
        "entry_point": "mcp-server-rna-mec",
        "description": "RNA / Magnetic Equivalent Circuit. ★ Lab specialty: "
                       "dynamic hysteresis MEC (Play + Cauer)",
        "primary_tools": ["rna_mec"],
        "related": ["magnetic-materials", "mor", "ih"],
        "tags": ["application"],
    },
    # ============================================================
    # ML / optimization
    # ============================================================
    "topology-optimization": {
        "subpackage": "radia_mcp.topology_optimization",
        "entry_point": "mcp-server-topology-optimization",
        "description": "Topology optimization: SIMP, level set, ON/OFF, MMA, "
                       "Wakao autoencoder+LS SynRM",
        "primary_tools": ["topology_optimization_usage"],
        "related": [
            "evolutionary",
            "motor",
            "bayesian-opt",
        ],
        "tags": ["ml", "optimization"],
    },
    "bayesian-opt": {
        "subpackage": "radia_mcp.bayesian_opt",
        "entry_point": "mcp-server-bayesian-opt",
        "description": "BO + GP regression + FMQA + surrogate models (57 lab "
                       "files; ARD kernel, PI-GP, multi-fidelity)",
        "primary_tools": ["bayesian_opt"],
        "related": ["evolutionary", "pinn", "radia-streamfunction",
                    "topology-optimization"],
        "tags": ["optimization"],
    },
    "evolutionary": {
        "subpackage": "radia_mcp.evolutionary",
        "entry_point": "mcp-server-evolutionary",
        "description": "GA / DE / PSO / CMA-ES / Immune / NSGA-II for EM",
        "primary_tools": ["evolutionary"],
        "related": ["bayesian-opt", "topology-optimization"],
        "tags": ["optimization"],
    },
    "data-assimilation": {
        "subpackage": "radia_mcp.data_assimilation",
        "entry_point": "mcp-server-data-assimilation",
        "description": "Kalman / EnKF / 4D-Var for EM state estimation + "
                       "sensor fusion",
        "primary_tools": ["data_assimilation"],
        "related": ["mor", "fusion-reactor"],
        "tags": ["optimization"],
    },
    "gnn": {
        "subpackage": "radia_mcp.gnn",
        "entry_point": "mcp-server-gnn",
        "description": "Graph Neural Networks for PDE/EM. Physics-Embedded "
                       "GNN, E(n)-GNN / NequIP / MACE",
        "primary_tools": ["gnn"],
        "related": ["pinn", "fem"],
        "tags": ["ml"],
    },
    "pinn": {
        "subpackage": "radia_mcp.pinn",
        "entry_point": "mcp-server-pinn",
        "description": "Physics-Informed Neural Networks + Gaussian Processes "
                       "for EM",
        "primary_tools": ["pinn_usage"],
        "related": ["gnn", "bayesian-opt", "fem"],
        "tags": ["ml"],
    },
    # ============================================================
    # Domain-specific
    # ============================================================
    "pcb": {
        "subpackage": "radia_mcp.pcb",
        "entry_point": "mcp-server-pcb",
        "description": "Wireless Power Transfer: coil + compensation (SS/LCC/"
                       "LCL), efficiency, IEC 61980 / SAE J2954, FOD, dynamic EV / "
                       "robot / bearingless motor, capacitive / microwave / metamaterial",
        "primary_tools": ["pcb_usage"],
        "related": [
            "litz-transmission",
            "maglev",
            "metamaterial",
            "peec",
        ],
        "tags": ["application"],
    },
    "ndt": {
        "subpackage": "radia_mcp.ndt",
        "entry_point": "mcp-server-ndt",
        "description": "Non-destructive testing: eddy current testing, "
                       "magnetic flux leakage, MFL signal analysis",
        "primary_tools": ["ndt_usage"],
        "related": ["ih", "magnetic-materials"],
        "tags": ["application"],
    },
    "metamaterial": {
        "subpackage": "radia_mcp.metamaterial",
        "entry_point": "mcp-server-metamaterial",
        "description": "Metamaterials: homogenization, effective medium, "
                       "periodic structures",
        "primary_tools": ["metamaterial_usage"],
        "related": ["pcb", "litz-transmission"],
        "tags": ["application"],
    },
    "nmr-mri": {
        "subpackage": "radia_mcp.nmr_mri",
        "entry_point": "mcp-server-nmr-mri",
        "description": "NMR/MRI: gradient coils, B0 shimming, RF coils, "
                       "field uniformity",
        "primary_tools": ["nmr_mri_usage"],
        "related": ["electromagnet", "accelerator", "radia-streamfunction"],
        "tags": ["application"],
    },
    "maglev": {
        "subpackage": "radia_mcp.maglev",
        "entry_point": "mcp-server-maglev",
        "description": "Magnetic levitation, UNIFIED: maglev systems "
                       "(EMS/EDS/PM/SC/Halbach) + levitation FORCE physics "
                       "(induction/EML/AMB/superconducting/diamagnetic/"
                       "Earnshaw/force-computation). Lab research: Radia "
                       "IEM<->FEM weak coupling + Cauer Ladder Network MOR "
                       "for control-coupled maglev (Yano, CAE-AI).",
        "primary_tools": ["maglev"],
        "related": ["motor", "pcb"],
        "tags": ["application"],
    },
    # ============================================================
    # Knowledge / meta / theory
    # ============================================================
    "team-benchmark": {
        "subpackage": "radia_mcp.team_benchmark",
        "entry_point": "mcp-server-team-benchmark",
        "description": "TEAM Workshop benchmark problems reference layer "
                       "(30 problems × physics class). ★ Lab core: 13, 20, 23, 32, 33b",
        "primary_tools": ["team_benchmark"],
        "related": ["fem", "bem", "motor"],
        "tags": ["application"],
    },
    "differential-forms": {
        "subpackage": "radia_mcp.differential_forms",
        "entry_point": "mcp-server-differential-forms",
        "description": "Differential forms / exterior calculus for EM: "
                       "de Rham complex, cohomology, EM forces theory",
        "primary_tools": ["differential_forms_usage", "forces_knowledge"],
        "related": ["fem", "mathematica"],
        "tags": ["theory"],
    },
    "mathematica": {
        "subpackage": "radia_mcp.mathematica",
        "entry_point": "mcp-server-mathematica",
        "description": "Mathematica recipes: vector calc, Kelvin transform, "
                       "symbolic Maxwell, evaluation pipeline",
        "primary_tools": ["mathematica_recipes", "mathematica_status"],
        "related": ["differential-forms", "radia-ngsolve", "figure",
                      "md2html"],
        "tags": ["theory"],
    },
    "md2html": {
        "subpackage": "radia_mcp.md2html",
        "entry_point": "mcp-server-md2html",
        "description": "Markdown -> self-contained HTML with MathJax v3 + "
                       "tables + fenced code + codehilite + base64 image "
                       "embed + [N] reference linking + UTF-8 / cp932 "
                       "fallback for legacy Japanese files.  Promoted "
                       "from mcp-server-document.",
        "primary_tools": ["md2html_convert"],
        "related": ["mathematica", "figure", "literature-index"],
        "tags": ["meta"],
    },
    "figure": {
        "subpackage": "radia_mcp.figure",
        "entry_point": "mcp-server-figure",
        "description": "Sugahara Lab publication-figure toolkit "
                       "(Times New Roman default): beamer/slide + "
                       "IEEE/IEEJ paper profiles, MATLAB + Matplotlib "
                       "snippets, lab style rules (units in parentheses, "
                       "no in-figure title).  Promoted from "
                       "mcp-server-document.",
        "primary_tools": ["figure_style_guide", "figure_size_for_target",
                            "paper_figure_profiles", "paper_figure_recipe",
                            "paper_figure_quality_rules"],
        "related": ["mathematica", "literature-index", "chart2d",
                      "md2html"],
        # 'meta' tag = cross-cutting utility usable by any paper /
        # digest, not bound to a single solver domain.
        "tags": ["meta"],
    },
    "chart2d": {
        "subpackage": "radia_mcp.chart2d",
        "entry_point": "mcp-server-chart2d",
        "description": "22 paper-quality 2D charts as MCP tools: line / "
                       "loglog / semilog / step / errorbar / fill_between"
                       " / bode / histogram / bar / box / violin / ecdf /"
                       " contour / contourf / pcolormesh / quiver / "
                       "streamplot / imshow / polar / scatter / phase "
                       "(Nyquist).  Each accepts return_mode='recipe' "
                       "(Python text) | 'image' (MCP Image inline) | "
                       "'both'.  Inherits radia_mcp.figure profile + "
                       "gate stack.",
        "primary_tools": ["chart2d_catalog", "chart2d_line", "chart2d_bode",
                            "chart2d_scatter", "chart2d_contourf"],
        "related": ["figure"],
        "tags": ["meta"],
    },
    "literature-index": {
        "subpackage": "radia_mcp.literature_index",
        "entry_point": "mcp-server-literature-index",
        "description": "★ Meta-MCP: full-text search across 3,889 lab "
                       "literature files in W:/00_電磁界解析",
        "primary_tools": ["literature_search", "literature_by_folder",
                            "literature_folder_tree", "literature_stats",
                            "literature_semantic_search"],
        "related": ["meta", "figure", "md2html"],
        "tags": ["meta"],
    },
    "document-meta": {
        "subpackage": "radia_mcp.document_meta",
        "entry_point": "mcp-server-document-meta",
        "description": "Cross-cutting document/repo helpers: deadline, "
                       "version diff, templates, lint-all, result-saving "
                       "notebook audits, examples->docs/validation_test "
                       "promotion audits, and root-level panels migration "
                       "impact checks.",
        "primary_tools": [
            "document_meta_notebook_result_audit",
            "document_meta_examples_notebook_audit",
            "document_meta_panel_layout_audit",
            "document_meta_lint_all",
        ],
        "related": ["meta"],
        "tags": ["meta"],
    },

    # ============================================================
    # Meta (the catalog itself — recommended first call)
    # ============================================================
    "meta": {
        "subpackage": "radia_mcp.meta",
        "entry_point": "mcp-server-radia-meta",
        "description": "★ RECOMMENDED FIRST CALL. Cross-server catalog "
                       "of all radia_mcp.* servers — answers \"which "
                       "tool covers concept X?\" without trial-and-error.",
        "primary_tools": ["radia_mcp_overview", "radia_mcp_get",
                            "radia_mcp_by_tag", "radia_mcp_related",
                            "radia_mcp_health"],
        "related": ["literature-index", "document-meta"],
        "tags": ["meta"],
    },

    # ============================================================
    # Panel review (Radia GUI panel skill chain)
    # ============================================================
    "panel-review": {
        "subpackage": "radia_mcp.panel_review",
        "entry_point": "mcp-server-panel-review",
        "description": "Radia GUI panel review skill chain "
                       "(panel-cli-diff / panel-review / panel-qt-test / "
                       "panel-preview / panel-smoke) + bug catalogue. "
                       "Surfaces the 13-check list and known panel bug "
                       "patterns through one queryable tool.",
        "primary_tools": ["panel_review"],
        "related": ["electromagnet", "ih", "motor"],
        "tags": ["meta"],
    },
}


# ============================================================
# External lab MCP packages — not in CATALOG (different PyPI dists,
# different repos), but reported alongside for cross-team discovery.
# Listed here so radia_mcp_overview() surfaces them as "you may also
# want to install ...". Not subject to the per-server status-tool
# policy (different repos, different conventions).
# ============================================================
EXTERNAL_PACKAGES: dict[str, dict[str, Any]] = {
    "optuna-mcp": {
        "pypi": "optuna-mcp",
        "github": "https://github.com/optuna/optuna-mcp",
        "install": "pip install --upgrade optuna optuna-mcp",
        "entry_point": "optuna-mcp",
        "description": "Official Optuna MCP server for Study, Trial, "
                       "Visualization, and Dashboard operations. Use "
                       "`--storage sqlite:///C:/temp/optuna_mcp.db` for "
                       "persistent lab sessions.",
        "related": [
            "bayesian-opt",
            "evolutionary",
            "topology-optimization",
            "radia-streamfunction",
        ],
    },
    "elf": {
        "pypi": "mcp-server-elf",
        "github": "https://github.com/ksugahar/mcp-server-elf",
        "install": "pip install mcp-server-elf",
        "entry_point": "mcp-server-elf",
        "description": "ELF600 BEM electromagnetic analysis docs "
                       "(MAGIC magnetostatic / ELFIN electrostatic / "
                       "BEAM particle tracking) — 16 tools + 1 prompt "
                       "over help HTM, examples, vendor wiki, ctypes API.",
        "related": ["radia-ngsolve", "motor", "mor"],
    },
    "comsol": {
        "pypi": "(fork-only, not on PyPI yet)",
        "github": "https://github.com/ksugahar/COMSOL_Multiphysics_MCP",
        "install": "git clone + pip install -e .",
        "entry_point": "comsol-mcp",
        "description": "COMSOL Multiphysics MCP — Sugahara lab fork of "
                       "wjc9011/COMSOL_Multiphysics_MCP with multilingual "
                       "embedding, Kelvin-transform tutorial, Cubit-to-"
                       "COMSOL bridge, Japanese chapter detection.",
        "related": ["fem", "matrix-solvers", "differential-forms"],
    },
    "mcp-server-document": {
        "pypi": "(LAB-private, not on PyPI)",
        "github": "(internal: S:/mcp-server)",
        "install": "pip install -e S:/mcp-server",
        "entry_point": "mcp-server-document",
        "description": "Paper/grant/poster/presentation writing helpers "
                       "(grant_writing_*, paper_writing_*, poster_*, "
                       "presentation_*, bibliography_*, pdf_*, ocr_*) — "
                       "thousands of skill-chain tools for academic "
                       "document QA. Lab-only.",
        "related": [],
    },
}


# Aliases: CLI script names (mcp-server-<X>) and underscore / hyphen
# variants resolve to the catalog short name.  Without these, a user
# who remembered the CLI binary name ('mcp-server-radia-meta') would
# call `radia_mcp_get('radia-meta')` and get "Unknown server", because
# the catalog uses short names like 'meta' (the 'radia-' prefix is only
# in the CLI script name, not the catalog key).
#
# Rule of thumb: any string a user might reasonably type to refer to a
# server -- the CLI suffix, the underscore subpackage name, the
# hyphen short name -- should resolve to the same CATALOG entry.
_ALIASES = {
    # CLI-name -> catalog key
    "radia-meta": "meta",
    "radia_meta": "meta",
}
# Auto-generate underscore variants for every hyphenated catalog key
# (e.g. 'magnetic-materials' resolves from both 'magnetic-materials' and
# 'magnetic_materials').  Done once at module import to keep `get()` O(1).
for _k in list(CATALOG.keys()):
    if "-" in _k:
        _ALIASES.setdefault(_k.replace("-", "_"), _k)
del _k  # don't leak the loop variable


def _resolve(name: str) -> str | None:
    """Resolve a user-typed name to a canonical CATALOG key.

    Accepts: catalog short names ('meta'), CLI-script suffixes
    ('radia-meta'), and underscore variants ('radia_meta',
    'magnetic_materials').  Returns None if no match.
    """
    if name in CATALOG:
        return name
    if name in _ALIASES:
        return _ALIASES[name]
    # Last-ditch: try the mcp-server-<name> prefix-stripped form
    if name.startswith("mcp-server-"):
        stripped = name[len("mcp-server-"):]
        if stripped in CATALOG:
            return stripped
        if stripped in _ALIASES:
            return _ALIASES[stripped]
    return None


def _resolve_external(name: str) -> str | None:
    """Resolve a user-typed name to an EXTERNAL_PACKAGES key."""
    if name in EXTERNAL_PACKAGES:
        return name
    if name.startswith("mcp-server-"):
        stripped = name[len("mcp-server-"):]
        if stripped in EXTERNAL_PACKAGES:
            return stripped
    return None


def _entry(name: str, info: dict[str, Any]) -> dict[str, Any]:
    """Return a public catalog entry with uniform runtime commands."""
    out = {"name": name, **info}
    entry_point = info.get("entry_point")
    if entry_point:
        out.setdefault("selftest_command", f"{entry_point} --selftest")
    return out


def list_all() -> list[dict]:
    """Return catalog entries as a flat list."""
    return [_entry(n, info) for n, info in CATALOG.items()]


def list_external() -> list[dict]:
    """External lab MCP packages (not in CATALOG; different repos)."""
    return [{"name": n, **info} for n, info in EXTERNAL_PACKAGES.items()]


def get(name: str) -> dict | None:
    """Look up a single server by short name, CLI-name, or alias.

    Resolves aliases first so callers can pass any of:
        'meta'           - catalog short name
        'radia-meta'     - CLI script suffix (mcp-server-radia-meta)
        'radia_meta'     - underscore-variant
        'mcp-server-meta' - full CLI name
    """
    key = _resolve(name)
    return _entry(key, CATALOG[key]) if key else None


def find_by_tag(tag: str) -> list[dict]:
    """All servers tagged with `tag` (e.g. 'optimization')."""
    return [
        _entry(n, info)
        for n, info in CATALOG.items()
        if tag in info.get("tags", [])
    ]


def find_related(name: str) -> list[dict]:
    """Servers/packages listed as `related` of `name` (alias-aware).

    External packages are not part of ``CATALOG`` and are therefore not
    imported or self-tested here.  They are still returned from related lookup
    so the meta server can steer agents from a radia-mcp knowledge lane to the
    official external MCP that should own a workflow, e.g. Optuna study/trial
    operation.
    """
    key = _resolve(name)
    ext_key = _resolve_external(name)
    if key is None and ext_key is None:
        return []
    info = CATALOG.get(key) if key is not None else EXTERNAL_PACKAGES.get(ext_key)
    if not info:
        return []

    related = []
    for r in info.get("related", []):
        if r in CATALOG:
            related.append(_entry(r, CATALOG[r]))
        elif r in EXTERNAL_PACKAGES:
            related.append({"name": r, "external": True, **EXTERNAL_PACKAGES[r]})

    if key is not None:
        for ext_name, ext_info in EXTERNAL_PACKAGES.items():
            if key in ext_info.get("related", []):
                related.append({"name": ext_name, "external": True, **ext_info})

    return related
