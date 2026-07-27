"""Universal Relaxation Network public API."""

try:
    from .universal_relaxation_network import (
        URNConfig,
        UniversalRelaxationNetwork,
        generate_spice_netlist,
        train_urn,
    )
    from .cauer_ladder_urn import (
        CauerLadderURN,
        CauerLadderURNConfig,
        CauerRationalFit,
        fit_rational_pole_zero,
        train_cauer_ladder_alternating,
        train_cauer_ladder_progressive,
        train_cauer_ladder_tail_then_polish,
        train_cauer_ladder_urn,
    )
    from .cln_peeling_urn import (
        CLNBranchFit,
        CLNPeelingConfig,
        CLNPeelingNetwork,
        CLNPeelingStage,
        train_cln_peeling_urn,
    )
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
    "CauerLadderURN",
    "CauerLadderURNConfig",
    "CauerRationalFit",
    "CLNBranchFit",
    "CLNPeelingConfig",
    "CLNPeelingNetwork",
    "CLNPeelingStage",
    "UniversalRelaxationNetwork",
    "YAdmittanceURN",
    "YAdmittanceURNActiveBasis",
    "YAdmittanceURNConfig",
    "active_basis_refit_config",
    "fit_rational_pole_zero",
    "generate_spice_netlist",
    "refit_y_admittance_active_bases",
    "s_domain_rmse",
    "train_cauer_ladder_alternating",
    "train_cauer_ladder_progressive",
    "train_cauer_ladder_tail_then_polish",
    "train_cauer_ladder_urn",
    "train_cln_peeling_urn",
    "train_stacked_y_urn",
    "train_urn",
    "train_y_admittance_urn",
]
