"""Declarative, startup-selected MCP capability packs."""

PACKS = {
    "radia-mcp-analysis": {
        "subpackage": "radia_mcp.radia_analysis",
        "description": "Field formulations, numerical methods and independent verification",
        "profiles": {
            "field": ("radia_ngsolve", "fem", "bem"),
            "verification": ("force", "team_benchmark"),
            "theory": ("matrix_solvers", "differential_forms", "mor"),
        },
    },
    "radia-mcp-learning": {
        "subpackage": "radia_mcp.radia_learning",
        "description": "Surrogate models, physics-informed learning and data assimilation",
        "profiles": {
            "all": ("gnn", "pinn", "data_assimilation"),
        },
    },
    "radia-mcp-optimize": {
        "subpackage": "radia_mcp.radia_optimize",
        "description": "Optimization: solver-neutral diagnostics, Bayesian/evolutionary search, topology and shape optimization, metamaterial design and stream-function coil design",
        "profiles": {
            "search": ("optimization", "bayesian_opt", "evolutionary"),
            "topology": ("topology_optimization", "metamaterial"),
            "coil": ("streamfunction", "fusion_reactor"),
        },
    },
    "radia-mcp-ih": {
        "subpackage": "radia_mcp.radia_ih",
        "description": "Induction heating: SIBC, ESIM, Karl iteration and screening physics",
        "profiles": {
            "all": ("ih", "force", "mor"),
        },
    },
    "radia-mcp-maglev": {
        "subpackage": "radia_mcp.radia_maglev",
        "description": "Magnetic levitation and motion: guidance, propulsion and force workflows",
        "profiles": {
            "all": ("maglev", "electromagnet", "force"),
        },
    },
    "radia-mcp-electromagnet": {
        "subpackage": "radia_mcp.radia_electromagnet",
        "description": "Electromagnets and accelerator magnets",
        "profiles": {
            "all": ("electromagnet", "accelerator", "magnetic_materials"),
        },
    },
    "radia-mcp-motor": {
        "subpackage": "radia_mcp.radia_motor",
        "description": "Electric machines: motor design, circuit coupling and materials",
        "profiles": {
            "all": ("motor", "electromagnet", "magnetic_materials"),
        },
    },
    "radia-mcp-pcb": {
        "subpackage": "radia_mcp.radia_pcb",
        "description": "PCB coils, wireless power, PEEC and conductor losses",
        "profiles": {
            "all": ("pcb", "peec", "litz_transmission"),
        },
    },
    "radia-mcp-acoustic": {
        "subpackage": "radia_mcp.radia_acoustic_workflows",
        "description": "Production and educational acoustic routes with explicit ownership",
        "profiles": {
            "production": ("radia_acoustic",),
            "education": ("acoustic_fembem",),
        },
    },
    "radia-mcp-rna": {
        "subpackage": "radia_mcp.radia_rna",
        "description": "Reluctance network analysis / magnetic equivalent circuits and hysteresis modeling",
        "profiles": {
            "all": ("rna_mec", "magnetic_materials"),
        },
    },
}


def modules_for(pack: str, profile: str = "all") -> tuple[str, ...]:
    """Resolve a profile before importing any optional domain dependencies."""
    profiles = PACKS[pack]["profiles"]
    if profile == "all":
        return tuple(dict.fromkeys(m for members in profiles.values() for m in members))
    if profile not in profiles:
        raise ValueError(f"Unknown profile {profile!r}; choose all or {', '.join(profiles)}")
    return profiles[profile]
