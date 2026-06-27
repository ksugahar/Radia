"""Universal Relaxation Network public API."""

try:
    from .universal_relaxation_network import (
        URNConfig,
        UniversalRelaxationNetwork,
        generate_spice_netlist,
        train_urn,
    )
except ModuleNotFoundError as exc:
    if exc.name == "torch":
        raise ModuleNotFoundError(
            "radia.urn requires PyTorch. Install it with `pip install radia[urn]`."
        ) from exc
    raise

__all__ = [
    "URNConfig",
    "UniversalRelaxationNetwork",
    "generate_spice_netlist",
    "train_urn",
]
