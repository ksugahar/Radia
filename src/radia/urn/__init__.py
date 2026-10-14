"""Universal Relaxation Network public API."""

try:
    from .universal_relaxation_network import (
        URNConfig,
        UniversalRelaxationNetwork,
        generate_spice_netlist,
        train_urn,
    )
    from .reduction import reduce_y_admittance_urn, redundant_y_bases
    from .y_admittance_urn import (
        StackedYURN,
        YAdmittanceURN,
        YAdmittanceURNActiveBasis,
        YAdmittanceURNConfig,
        active_basis_refit_config,
        refit_y_admittance_active_bases,
        s_domain_rmse,
        train_stacked_y_urn,
        train_y_admittance_urn,
    )
except ModuleNotFoundError as exc:
    if exc.name == "torch":
        raise ModuleNotFoundError(
            "radia.urn requires PyTorch. Install it with `pip install radia[urn]`."
        ) from exc
    raise

__all__ = [
    "URNConfig",
    "StackedYURN",
    "UniversalRelaxationNetwork",
    "YAdmittanceURN",
    "YAdmittanceURNActiveBasis",
    "YAdmittanceURNConfig",
    "active_basis_refit_config",
    "generate_spice_netlist",
    "refit_y_admittance_active_bases",
    "reduce_y_admittance_urn",
    "redundant_y_bases",
    "s_domain_rmse",
    "train_stacked_y_urn",
    "train_urn",
    "train_y_admittance_urn",
]
