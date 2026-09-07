"""Declarative, startup-selected MCP capability packs."""

PACKS = {
    "radia-design": {
        "description": "Optimization, surrogate models and data assimilation",
        "profiles": {
            "optimization": ("optimization", "topology_optimization", "bayesian_opt", "evolutionary"),
            "learning": ("gnn", "pinn", "data_assimilation"),
        },
    },
    "radia-electrical": {
        "description": "PEEC, wireless power and conductor losses",
        "profiles": {"all": ("peec", "pcb", "litz_transmission")},
    },
    "radia-motion": {
        "description": "Motor, MagLev, electromagnet, accelerator and rotating IH workflows",
        "profiles": {
            "motor": ("motor", "electromagnet", "magnetic_materials"),
            "maglev": ("maglev", "electromagnet", "force"),
            "ih": ("ih", "force", "mor"),
            "accelerator": ("accelerator", "electromagnet"),
        },
    },
    "radia-analysis": {
        "description": "Field formulations, numerical methods and independent verification",
        "profiles": {
            "field": ("radia_ngsolve", "fem", "bem"),
            "verification": ("force", "team_benchmark"),
            "theory": ("matrix_solvers", "differential_forms", "mor"),
        },
    },
    "radia-acoustic-workflows": {
        "description": "Production and educational acoustic routes with explicit ownership",
        "profiles": {
            "production": ("radia_acoustic",),
            "education": ("acoustic_fembem",),
        },
    },
    "document-ops": {
        "description": "Document inspection, conversion and project checks",
        "profiles": {
            "inspect": ("pdf",),
            "convert": ("doc_convert", "md2html"),
            "project": ("document_meta", "research_project"),
        },
    },
    "radia-publication": {
        "description": "Paper, slide, poster, figure and bibliography workflows",
        "profiles": {
            "paper": ("paper_writing", "bibliography"),
            "poster": ("paper_writing", "chart2d", "bibliography"),
            "figures": ("chart2d",),
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
