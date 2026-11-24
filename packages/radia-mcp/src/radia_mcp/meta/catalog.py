"""Authoritative catalog of all radia_mcp servers.

Source of truth: this dict together with the package entry points. Generated
documentation may project this catalog, but no workstation-local MCP client
configuration is part of the repository contract.

The catalog includes the `meta` server itself (self-referential entry)
so that `radia_mcp_overview` reports a complete picture and
`radia_mcp_by_tag("meta")` finds the entry-point.

Each entry:
    "<short_name>": {
        "subpackage": "radia_mcp.<name>",
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
    "build123d-mcp": {
        "subpackage": "radia_mcp.build123d",
        "description": "build123d STEP authoring (CAD-as-code) + Cubit interop",
        "primary_tools": ["build123d_api", "build123d_run", "execute_build123d",
                            "build123d_render_views", "build123d_diff_shapes",
                            "build123d_to_cubit_hex"],
        "related": ["cubit-mcp", "gmsh-mcp", "onshape-mcp", "radia-ngsolve"],
        "tags": ["cad"],
    },
    "ngsolve-mcp.age": {
        "subpackage": "radia_mcp.age",
        "description": "Air-Gap Element (AGE) analytical/harmonic rotor-stator coupling for NGSolve",
        "primary_tools": ["age_usage"],
        "related": [],
        "tags": ["fem"],
    },
    "ngsolve-mcp.dtn": {
        "subpackage": "radia_mcp.dtn",
        "description": "Dirichlet-to-Neumann (DtN) open-boundary operator and coarse-mesh coupling for NGSolve",
        "primary_tools": ["dtn_usage"],
        "related": [],
        "tags": ["fem"],
    },
    "ngsolve-mcp.gauging": {
        "subpackage": "radia_mcp.gauging",
        "description": "Electromagnetic gauging (Tree-Cotree, Coulomb/Lorenz, open-boundary)",
        "primary_tools": ["gauging_usage"],
        "related": [],
        "tags": ["fem"],
    },
    "ngsolve-mcp.kelvin": {
        "subpackage": "radia_mcp.kelvin",
        "description": "Kelvin transformation open-boundary exterior inversion for NGSolve/Radia",
        "primary_tools": ["kelvin_usage"],
        "related": [],
        "tags": ["fem"],
    },
    "ngsolve-mcp.weak-form": {
        "subpackage": "radia_mcp.weak_form",
        "description": "NGSolve weak-form variational formulations (HCurl/HDiv/H1, BilinearForm, SymbolicBFI)",
        "primary_tools": ["weak_form_usage"],
        "related": [],
        "tags": ["fem"],
    },
    "onshape-mcp": {
        "subpackage": "radia_mcp.onshape",
        "description": "Research-only, version-pinned Onshape STEP, Parasolid and ACIS export",
        "primary_tools": ["onshape_export_cad"],
        "related": ["build123d-mcp"],
        "tags": ["cad"],
    },
    "gmsh-mcp": {
        "subpackage": "radia_mcp.gmsh",
        "description": "GMSH MSH v4.2 inspect/validate + post-display "
                       "launch artifacts + visualization-only policy lint",
        "primary_tools": ["gmsh_inspect_msh", "gmsh_validate_msh",
                          "gmsh_validate_geo", "gmsh_exec",
                          "gmsh_render",
                          "gmsh_write_post_launch_artifact",
                          "gmsh_usage", "gmsh_reference"],
        "related": ["blender-mcp", "build123d-mcp", "cubit-mcp", "radia-ngsolve"],
        "audit_command": "radia-mcp serve gmsh-mcp --selftest --audit-repo",
        "tags": ["mesh"],
    },
    "blender-mcp": {
        "subpackage": "radia_mcp.blender",
        "description": "Solver-faithful 3-D model animation through an optional "
                       "headless Blender worker (write + separate reopen check)",
        "primary_tools": ["blender_field_animation", "blender_usage"],
        "related": ["gmsh-mcp"],
        "tags": ["mesh"],
    },
    # ============================================================
    # FEM / BEM / specialty solvers
    # ============================================================
    "radia-ngsolve": {
        "subpackage": "radia_mcp.radia_ngsolve",
        "description": "Radia + NGSolve: Kelvin / sparsesolv / PEEC / "
                       "analytical formulas / lint",
        "primary_tools": ["kelvin_transformation", "ngsolve_usage",
                            "analytical_formulas", "peec_inductance"],
        "related": [
            "bem",
            "build123d-mcp",
            "cubit-mcp",
            "fem",
            "force",
            "gmsh-mcp",
            "radia-acoustic",
            "matlab-mcp-lab",
            "acoustic-fembem",
            "radia-mcp.mathematica",
            "matrix-solvers",
            "mor",
            "peec",
            "streamfunction",
        ],
        "tags": ["fem", "solver"],
    },
    "force": {
        "subpackage": "radia_mcp.force",
        "description": "Shared electromagnetic-force and torque method "
                       "selection, static/phasor sample integration, virtual "
                       "work, air-gap torque, and validation gates",
        "primary_tools": [
            "force",
            "force_result",
            "force_recipe",
            "force_lorentz",
            "force_maxwell_surface",
            "force_time_average_lorentz",
            "force_time_average_maxwell_surface",
            "force_coenergy_torque",
            "force_air_gap_torque_samples",
            "force_method_agreement_gate",
            "force_validation_guide",
        ],
        "related": [
            "differential-forms",
            "maglev",
            "motor",
            "radia-ngsolve",
            "team-benchmark",
        ],
        "tags": ["solver", "theory"],
    },
    "matlab-mcp-lab": {
        "subpackage": "radia_mcp.matlab",
        "description": "Official MATLAB MCP composition, Radia/NGSolve MEX "
                       "capability contract, radia-optuna MATLAB/Simulink "
                       "difference support, and 43 generic ML/RL gates",
        "primary_tools": ["matlab_agent_guide",
                            "matlab_extension_contract",
                            "matlab_official_server_config",
                            "matlab_official_server_status",
                            "matlab_radia_mex_contract",
                            "matlab_optuna_mcp_route",
                            "matlab_optuna_simulink_contract",
                            "matlab_radia_acoustic_interface_contract"],
        "related": ["acoustic-fembem", "radia-ngsolve"],
        "tags": ["ml", "optimization", "solver"],
    },
    "acoustic-fembem": {
        "subpackage": "radia_mcp.acoustic_fembem",
        "description": "Readable MATLAB P1 acoustic FEM/BEM education solver, "
                       "convolution quadrature, and validation gates",
        "primary_tools": ["acoustic_fembem_agent_guide",
                            "acoustic_fembem_extension_contract",
                            "acoustic_fembem_server_config"],
        "related": ["matlab-mcp-lab", "radia-ngsolve", "bem", "fem"],
        "tags": ["fem", "theory"],
    },
    "radia-acoustic": {
        "subpackage": "radia_mcp.radia_acoustic",
        "description": "Production NGSolve/ngsolve.bem acoustics, FSI, and CQ",
        "primary_tools": ["radia_acoustic_usage", "radia_acoustic_fsi_preflight", "radia_acoustic_cq_grid"],
        "related": ["radia-ngsolve", "bem"],
        "tags": ["fem", "solver"],
    },
    "streamfunction": {
        "subpackage": "radia_mcp.streamfunction",
        "description": "Stream-function coil design: (ACA+)+TSVD least-norm, "
                       "FE-direct psi, regularisation / folded-Tikhonov Pareto "
                       "front, single-stroke chain, sheet-metal levers; also "
                       "NMR/MRI gradient and shim coil design",
        "primary_tools": ["streamfunction"],
        "related": [
            "radia-ngsolve",
            "bayesian-opt",
            "electromagnet",
            "peec",
        ],
        "tags": ["optimization", "solver"],
    },
    "fem": {
        "subpackage": "radia_mcp.fem",
        "description": "FEM formulations theory layer (A-Omega / T-Omega / H / "
                       "Reduced / Darwin, edge / HO / XFEM / IGA / DG, gauging "
                       "+ Kelvin, MSFEM, Schur circuit coupling, NGSolve hierarchical)",
        "primary_tools": ["fem_overview"],
        "related": [
            "bem",
            "differential-forms",
            "gnn",
            "matrix-solvers",
            "pinn",
            "radia-ngsolve",
            "team-benchmark",
            "acoustic-fembem",
        ],
        "tags": ["fem", "theory"],
    },
    "bem": {
        "subpackage": "radia_mcp.bem",
        "description": "MoM/BEM theory: RWG, EFIE/MFIE/CFIE/PMCHWT, "
                       "Loop-Star, Calderon, Radia HDiv-VIM, HACApK, FEM-BEM",
        "primary_tools": ["bem_overview"],
        "related": [
            "fem",
            "peec",
            "radia-ngsolve",
            "team-benchmark",
            "acoustic-fembem",
            "radia-acoustic",
        ],
        "tags": ["fem", "theory"],
    },
    "matrix-solvers": {
        "subpackage": "radia_mcp.matrix_solvers",
        "description": "Sparse solver theory + decision tree: Krylov "
                       "(CG/BiCGSTAB/GMRES/COCG/COCR/IDR), preconditioners "
                       "(AMG, Hiptmair-Xu AMS), Biro-Preis A-V, tree-cotree",
        "primary_tools": ["matrix_solvers_overview"],
        "related": ["radia-ngsolve", "fem"],
        "tags": ["solver", "theory"],
    },
    "mor": {
        "subpackage": "radia_mcp.mor",
        "description": "Model Order Reduction: POD, PRIMA/Krylov, "
                       "hyperreduction (DEIM)",
        "primary_tools": ["mor_systematic"],
        "related": ["data-assimilation", "radia-ngsolve", "rna-mec"],
        "tags": ["solver"],
    },
    # ============================================================
    # Physics applications
    # ============================================================
    "ih": {
        "subpackage": "radia_mcp.ih",
        "description": "Induction heating: SIBC, ESIM, Karl iteration, "
                       "workpiece coupling",
        "primary_tools": ["induction_heating", "ih_sibc", "ih_esim"],
        "related": [
            "litz-transmission",
            "magnetic-materials",
            "peec",
            "rna-mec",
        ],
        "tags": ["application"],
    },
    "peec": {
        "subpackage": "radia_mcp.peec",
        "description": "PEEC filament/panel, FastHenry parser, HOIBC, "
                       "Carstensen AC copper loss",
        "primary_tools": ["peec_usage", "peec_hoibc", "peec_carstensen_ac_loss"],
        "related": [
            "bem",
            "ih",
            "litz-transmission",
            "radia-ngsolve",
            "streamfunction",
            "pcb",
        ],
        "tags": ["application"],
    },
    "electromagnet": {
        "subpackage": "radia_mcp.electromagnet",
        "description": "Accelerator electromagnet: CoilBuilder, Hantila, "
                       "Play/Energy hysteresis",
        "primary_tools": ["electromagnet_usage"],
        "related": [
            "accelerator",
            "fusion-reactor",
            "magnetic-materials",
            "motor",
            "streamfunction",
        ],
        "tags": ["application"],
    },
    "motor": {
        "subpackage": "radia_mcp.motor",
        "description": "Motor analysis: ONELAB transient, Hollaus effective "
                       "material (lamination), Wakao autoencoder topology, "
                       "Kaimori-Mifune Darwin TD",
        "primary_tools": ["motor_validation_router", "motor_field_quick_check"],
        "related": [
            "electromagnet",
            "force",
            "maglev",
            "magnetic-materials",
            "team-benchmark",
            "topology-optimization",
        ],
        "tags": ["application"],
    },
    "accelerator": {
        "subpackage": "radia_mcp.accelerator",
        "description": "Accelerator physics: beam optics, dipole/quad/sext "
                       "magnets, undulator/wiggler",
        "primary_tools": ["accelerator"],
        "related": ["electromagnet", "fusion-reactor"],
        "tags": ["application"],
    },
    "fusion-reactor": {
        "subpackage": "radia_mcp.fusion_reactor",
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
        "description": "Magnetic materials: hysteresis (Play/Energy lab core), "
                       "iron loss (Bertotti/Steinmetz/iGSE), JIS silicon steel, "
                       "PM datasheets, Osborn demag factor",
        "primary_tools": ["magnetic_materials_hysteresis", "magnetic_materials_iron_loss"],
        "related": [
            "electromagnet",
            "ih",
            "motor",
            "rna-mec",
        ],
        "tags": ["application"],
    },
    "litz-transmission": {
        "subpackage": "radia_mcp.litz_transmission",
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
        "description": "RNA / Magnetic Equivalent Circuit. ★ Lab specialty: "
                       "Play-model hysteresis MEC",
        "primary_tools": ["rna_mec"],
        "related": ["magnetic-materials", "mor", "ih"],
        "tags": ["application"],
    },
    # ============================================================
    # ML / optimization
    # ============================================================
    "topology-optimization": {
        "subpackage": "radia_mcp.topology_optimization",
        "description": "Topology optimization: SIMP, level set, ON/OFF, MMA, "
                       "Wakao autoencoder+LS SynRM",
        "primary_tools": ["topology_opt_topology_derivative"],
        "related": [
            "evolutionary",
            "motor",
            "bayesian-opt",
        ],
        "tags": ["ml", "optimization"],
    },
    "bayesian-opt": {
        "subpackage": "radia_mcp.bayesian_opt",
        "description": "BO + GP regression + FMQA + surrogate models (57 lab "
                       "files; ARD kernel, PI-GP, multi-fidelity)",
        "primary_tools": ["bayesian_opt"],
        "related": ["evolutionary", "pinn", "streamfunction",
                    "topology-optimization"],
        "tags": ["optimization"],
    },
    "evolutionary": {
        "subpackage": "radia_mcp.evolutionary",
        "description": "GA / DE / PSO / CMA-ES / Immune / NSGA-II for EM",
        "primary_tools": ["evolutionary"],
        "related": ["bayesian-opt", "topology-optimization"],
        "tags": ["optimization"],
    },
    "data-assimilation": {
        "subpackage": "radia_mcp.data_assimilation",
        "description": "Kalman / EnKF / 4D-Var for EM state estimation + "
                       "sensor fusion",
        "primary_tools": ["data_assimilation"],
        "related": ["mor", "fusion-reactor"],
        "tags": ["optimization"],
    },
    "gnn": {
        "subpackage": "radia_mcp.gnn",
        "description": "Graph Neural Networks for PDE/EM. Physics-Embedded "
                       "GNN, E(n)-GNN / NequIP / MACE",
        "primary_tools": ["gnn"],
        "related": ["pinn", "fem"],
        "tags": ["ml"],
    },
    "pinn": {
        "subpackage": "radia_mcp.pinn",
        "description": "Physics-Informed Neural Networks + Gaussian Processes "
                       "for EM",
        "primary_tools": ["pinn"],
        "related": ["gnn", "bayesian-opt", "fem"],
        "tags": ["ml"],
    },
    # ============================================================
    # Domain-specific
    # ============================================================
    "pcb": {
        "subpackage": "radia_mcp.pcb",
        "description": "Wireless Power Transfer: coil + compensation (SS/LCC/"
                       "LCL), efficiency, IEC 61980 / SAE J2954, FOD, dynamic EV / "
                       "robot / bearingless motor, capacitive / microwave / metamaterial",
        "primary_tools": ["pcb_overview"],
        "related": [
            "litz-transmission",
            "maglev",
            "metamaterial",
            "peec",
        ],
        "tags": ["application"],
    },
    "metamaterial": {
        "subpackage": "radia_mcp.metamaterial",
        "description": "Metamaterials: homogenization, effective medium, "
                       "periodic structures",
        "primary_tools": ["metamaterial"],
        "related": ["pcb", "litz-transmission"],
        "tags": ["application"],
    },
    "maglev": {
        "subpackage": "radia_mcp.maglev",
        "description": "Magnetic levitation, UNIFIED: maglev systems "
                       "(EMS/EDS/PM/SC/Halbach) + levitation FORCE physics "
                       "(induction/EML/AMB/superconducting/diamagnetic/"
                       "Earnshaw/force-computation). Lab research: Radia "
                       "IEM<->FEM weak coupling + Foster modal reduced "
                       "models for control-coupled maglev.",
        "primary_tools": ["maglev"],
        "related": ["force", "motor", "pcb"],
        "tags": ["application"],
    },
    # ============================================================
    # Knowledge / meta / theory
    # ============================================================
    "team-benchmark": {
        "subpackage": "radia_mcp.team_benchmark",
        "description": "TEAM Workshop benchmark problems reference layer "
                       "(30 problems × physics class). ★ Lab core: 13, 20, 23, 32, 33b",
        "primary_tools": ["team_catalog"],
        "related": ["fem", "bem", "force", "motor"],
        "tags": ["application"],
    },
    "differential-forms": {
        "subpackage": "radia_mcp.differential_forms",
        "description": "Visual differential geometry / exterior calculus "
                       "for EM: Cartan, de Rham, cohomology, EM forces theory",
        "primary_tools": ["differential_forms_visual_geometry",
                          "differential_forms_geometry_gate",
                          "differential_forms_maxwell"],
        "related": ["fem", "force", "radia-mcp.mathematica"],
        "tags": ["theory"],
    },
    "radia-mcp.mathematica": {
        "subpackage": "radia_mcp.mathematica",
        "description": "Mathematica recipes: vector calc, Kelvin transform, "
                       "symbolic Maxwell, evaluation pipeline",
        "primary_tools": [
            "mathematica_verification_guide",
            "mathematica_check_identities",
            "mathematica_run_script",
            "mathematica_status",
        ],
        "related": ["differential-forms", "radia-ngsolve"],
        "tags": ["theory"],
    },

    # ============================================================
    # Meta (the catalog itself — recommended first call)
    # ============================================================
    "radia-mcp.meta": {
        "subpackage": "radia_mcp.meta",
        "description": "★ RECOMMENDED FIRST CALL. Cross-server catalog "
                       "of all radia_mcp.* servers — answers \"which "
                       "tool covers concept X?\" without trial-and-error.",
        "primary_tools": ["radia_mcp_overview", "radia_mcp_get",
                            "radia_mcp_by_tag", "radia_mcp_related",
                            "radia_mcp_health",
                              "radia_mcp_validate_solver_artifact",
                              "radia_mcp_operations", "tests_design_guide"],
        "related": [ ],
        "tags": ["meta"],
    },

    # ============================================================
    # Panel review (Radia Simulink application-block contract)
    # ============================================================
}


# ============================================================
# External lab MCP packages — not in CATALOG (different PyPI dists,
# different repos), but reported alongside for cross-team discovery.
# Listed here so radia_mcp_overview() surfaces them as "you may also
# want to install ...". Not subject to the per-server status-tool
# policy (different repos, different conventions).
# ============================================================
EXTERNAL_PACKAGES: dict[str, dict[str, Any]] = {
    "cubit-mcp": {
        "pypi": "cubit-mesh-export",
        "github": "https://github.com/ksugahar/cubit-mesh-export",
        "install": "pip install cubit-mesh-export",
        "subpackage": "cubit_mesh_export.mcp",
        "selftest_command": "mcp-server-cubit --selftest",
        "audit_command": "mcp-server-cubit --selftest --audit-repo",
        "description": (
            "Independent Cubit APREPRO/headless meshing, export and MCP; no Radia required. "
            "Direct mesh export uses the Cubit mesh. "
            "Netgen/NGSolve are used for Python curving and solver-side mesh checks."
        ),
        "related": ["build123d-mcp", "gmsh-mcp", "radia-ngsolve", "topology-optimization"],
    },
    "optuna-mcp": {
        "pypi": "optuna-mcp",
        "github": "https://github.com/optuna/optuna-mcp",
        "install": "pip install --upgrade optuna optuna-mcp",
        "entry_point": "optuna-mcp",
        "description": "Official Optuna MCP server for Study, Trial, "
                       "Visualization, and Dashboard operations. Use "
                       "`--storage sqlite:///C:/temp/optuna_mcp.db` for "
                       "persistent lab sessions.",
        "owns": "Every shared Optuna operation present in the live MCP tools/list response.",
        "radia_difference_tool": "matlab-mcp-lab:matlab_optuna_mcp_route",
        "license": "MIT",
        "validation_storage": "fresh per-run temporary SQLite database",
        "automated_dashboard": False,
        "related": [
            "bayesian-opt",
            "evolutionary",
            "topology-optimization",
            "streamfunction",
            "matlab-mcp-lab",
        ],
    },
    "elf-magic-mcp": {
        "pypi": "mcp-server-elf",
        "github": "https://github.com/ksugahar/mcp-server-elf",
        "install": "pip install mcp-server-elf",
        "description": "ELF600 BEM electromagnetic analysis docs "
                       "(MAGIC magnetostatic / ELFIN electrostatic / "
                       "BEAM particle tracking) — 16 tools + 1 prompt "
                       "over help HTM, examples, vendor wiki, ctypes API.",
        "related": ["radia-ngsolve", "motor", "mor"],
    },
}


# Lookup is exact on catalog keys; retired names and underscore
# variants have no aliases (user decision 2026-11-16).
from ..capability_packs import PACKS, modules_for

for _pack, _spec in PACKS.items():
    CATALOG[_pack] = {
        "subpackage": _spec["subpackage"],
        "description": _spec["description"],
        "primary_tools": ["capability_pack_status"],
        "related": [],
        "tags": ["meta"],
        "kind": "capability-pack",
        "profiles": ["all", *[p for p in _spec["profiles"] if p != "all"]],
        "members": list(modules_for(_pack)),
    }

# Keep old entry points during migration; expose recommended replacements.
for _name, _info in CATALOG.items():
    if _name in PACKS:
        continue
    _module = _info["subpackage"].removeprefix("radia_mcp.")
    _packs = [p for p in PACKS if _module in modules_for(p)]
    if _packs:
        _info["capability_packs"] = _packs

_PACK_MEMBERS = {m for _p in PACKS for m in modules_for(_p)}
for _name, _info in CATALOG.items():
    _module = _info["subpackage"].removeprefix("radia_mcp.")
    if _info.get("kind") == "capability-pack" or _module not in _PACK_MEMBERS:
        _info["entry_point"] = f"radia-mcp serve {_name}"
    else:
        _info["kind"] = "domain"
del _name, _info, _module

def _resolve(name: str) -> str | None:
    """Return ``name`` when it is a catalog key, else None.

    Lookup is exact: no underscore variants and no retired names
    (user decision 2026-11-16: no aliases).
    """
    return name if name in CATALOG else None


def _resolve_external(name: str) -> str | None:
    """Resolve a user-typed name to an EXTERNAL_PACKAGES key."""
    if name in EXTERNAL_PACKAGES:
        return name
    return None


def _entry(name: str, info: dict[str, Any]) -> dict[str, Any]:
    """Return a public catalog entry with uniform runtime commands."""
    out = {"name": name, **info}
    out.setdefault("python_name", info.get("subpackage"))
    entry_point = info.get("entry_point")
    if entry_point:
        out.setdefault("server_id", name)
        out.setdefault("command", entry_point)
        out.setdefault("selftest_command", f"{entry_point} --selftest")
    return out


def list_all() -> list[dict]:
    """Return catalog entries as a flat list."""
    return [_entry(n, info) for n, info in CATALOG.items()]


def list_external() -> list[dict]:
    """Separately distributed MCP packages, including monorepo siblings."""
    return [{"name": n, **info} for n, info in EXTERNAL_PACKAGES.items()]


def get(name: str) -> dict | None:
    """Look up a single server by its exact catalog key."""
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
