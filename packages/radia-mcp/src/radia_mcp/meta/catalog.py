"""Authoritative catalog of all 37 radia_mcp servers.

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
        "primary_tools": ["cubit_exec", "cubit_mesh_auto", "cubit_docs"],
        "related": ["build123d", "gmsh", "interop"],
        "tags": ["cad", "mesh"],
    },
    "build123d": {
        "subpackage": "radia_mcp.build123d",
        "entry_point": "mcp-server-build123d",
        "description": "build123d STEP authoring (CAD-as-code) + Cubit interop",
        "primary_tools": ["build123d_api", "execute_build123d",
                            "build123d_to_cubit_hex"],
        "related": ["cubit", "interop"],
        "tags": ["cad"],
    },
    "interop": {
        "subpackage": "radia_mcp.interop",
        "entry_point": "mcp-server-radia-interop",
        "description": "Cross-CAD interop (STEP/IGES/CadQuery <-> Cubit/Netgen)",
        "primary_tools": ["any_step_to_cubit_hex", "freecad_to_cubit_hex"],
        "related": ["cubit", "build123d"],
        "tags": ["cad"],
    },
    "gmsh": {
        "subpackage": "radia_mcp.gmsh",
        "entry_point": "mcp-server-gmsh",
        "description": "GMSH MSH v4.1 inspect/validate/convert/write_node_data",
        "primary_tools": ["gmsh_usage", "gmsh_reference"],
        "related": ["cubit"],
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
            "fem",
            "mathematica",
            "matrix-solvers",
            "mor",
            "peec",
        ],
        "tags": ["fem", "solver"],
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
            "wpt",
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
            "fusion",
            "magnetic-materials",
            "motor",
            "nmr-mri",
            "panel-review",
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
            "maglev-linear",
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
        "related": ["electromagnet", "fusion", "nmr-mri"],
        "tags": ["application"],
    },
    "fusion": {
        "subpackage": "radia_mcp.fusion",
        "entry_point": "mcp-server-fusion",
        "description": "Fusion reactor magnets: tokamak ITER + stellarator "
                       "LHD/W7-X/heliotron lineage",
        "primary_tools": ["fusion"],
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
            "wpt",
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
            "mcmc",
            "motor",
            "optuna",
        ],
        "tags": ["ml", "optimization"],
    },
    "optuna": {
        "subpackage": "radia_mcp.optuna",
        "entry_point": "mcp-server-optuna",
        "description": "Optuna black-box optimization (Sano-Akiba-Imamura "
                       "2023 textbook)",
        "primary_tools": ["optuna_usage", "optuna_algorithm",
                            "optuna_lab_applications"],
        "related": [
            "bayesian-opt",
            "evolutionary",
            "mcmc",
            "topology-optimization",
        ],
        "tags": ["optimization"],
    },
    "bayesian-opt": {
        "subpackage": "radia_mcp.bayesian_opt",
        "entry_point": "mcp-server-bayesian-opt",
        "description": "BO + GP regression + FMQA + surrogate models (57 lab "
                       "files; ARD kernel, PI-GP, multi-fidelity)",
        "primary_tools": ["bayesian_opt"],
        "related": ["optuna", "mcmc", "pinn"],
        "tags": ["optimization"],
    },
    "evolutionary": {
        "subpackage": "radia_mcp.evolutionary",
        "entry_point": "mcp-server-evolutionary",
        "description": "GA / DE / PSO / CMA-ES / Immune / NSGA-II for EM",
        "primary_tools": ["evolutionary"],
        "related": ["optuna", "mcmc", "topology-optimization"],
        "tags": ["optimization"],
    },
    "mcmc": {
        "subpackage": "radia_mcp.mcmc",
        "entry_point": "mcp-server-mcmc",
        "description": "MCMC + MCTS + SPM for EM. ★ Hokkaido Sato 2023 MCTS "
                       "PM motor + Yin 2024 inductor + Saotome 1995 SPM",
        "primary_tools": ["mcmc_algorithms", "mcmc_libraries",
                            "mcmc_em_applications", "mcmc_mcts_lab", "mcmc_spm"],
        "related": [
            "bayesian-opt",
            "data-assimilation",
            "evolutionary",
            "optuna",
            "topology-optimization",
        ],
        "tags": ["optimization"],
    },
    "data-assimilation": {
        "subpackage": "radia_mcp.data_assimilation",
        "entry_point": "mcp-server-data-assimilation",
        "description": "Kalman / EnKF / 4D-Var for EM state estimation + "
                       "sensor fusion",
        "primary_tools": ["data_assimilation"],
        "related": ["mcmc", "mor", "fusion"],
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
    "wpt": {
        "subpackage": "radia_mcp.wpt",
        "entry_point": "mcp-server-wpt",
        "description": "Wireless Power Transfer: coil + compensation (SS/LCC/"
                       "LCL), efficiency, IEC 61980 / SAE J2954, FOD, dynamic EV / "
                       "robot / bearingless motor, capacitive / microwave / metamaterial",
        "primary_tools": ["wpt_usage"],
        "related": [
            "litz-transmission",
            "maglev-linear",
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
        "related": ["wpt", "litz-transmission"],
        "tags": ["application"],
    },
    "nmr-mri": {
        "subpackage": "radia_mcp.nmr_mri",
        "entry_point": "mcp-server-nmr-mri",
        "description": "NMR/MRI: gradient coils, B0 shimming, RF coils, "
                       "field uniformity",
        "primary_tools": ["nmr_mri_usage"],
        "related": ["electromagnet", "accelerator"],
        "tags": ["application"],
    },
    "maglev-linear": {
        "subpackage": "radia_mcp.maglev_linear",
        "entry_point": "mcp-server-maglev-linear",
        "description": "Maglev (EMS/EDS/SCMaglev/Halbach/bearingless ★) + "
                       "linear drives (LIM/LSM). Lab specialty: bearingless + WPT",
        "primary_tools": ["maglev_linear"],
        "related": ["motor", "wpt"],
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
        "related": ["differential-forms", "radia-ngsolve"],
        "tags": ["theory"],
    },
    "literature-index": {
        "subpackage": "radia_mcp.literature_index",
        "entry_point": "mcp-server-literature-index",
        "description": "★ Meta-MCP: full-text search across 3,889 lab "
                       "literature files in W:/00_電磁界解析",
        "primary_tools": ["literature_search", "literature_by_folder",
                            "literature_folder_tree", "literature_stats",
                            "literature_semantic_search"],
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
        "related": ["literature-index"],
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


def list_all() -> list[dict]:
    """Return catalog entries as a flat list."""
    return [{"name": n, **info} for n, info in CATALOG.items()]


def list_external() -> list[dict]:
    """External lab MCP packages (not in CATALOG; different repos)."""
    return [{"name": n, **info} for n, info in EXTERNAL_PACKAGES.items()]


def get(name: str) -> dict | None:
    """Look up a single server by short name."""
    return CATALOG.get(name)


def find_by_tag(tag: str) -> list[dict]:
    """All servers tagged with `tag` (e.g. 'optimization')."""
    return [
        {"name": n, **info}
        for n, info in CATALOG.items()
        if tag in info.get("tags", [])
    ]


def find_related(name: str) -> list[dict]:
    """Servers listed as `related` of `name`."""
    info = CATALOG.get(name)
    if not info:
        return []
    return [
        {"name": r, **CATALOG[r]}
        for r in info.get("related", [])
        if r in CATALOG
    ]
