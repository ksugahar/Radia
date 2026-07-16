"""Public HDiv-VIM API naming contract.

The public ``radia.vim`` surface follows NGSolve-style names.  The old
implementation-shaped snake_case names are intentionally not compatibility
aliases.
"""

import radia.vim as vim


def test_ngsolve_style_hdiv_vim_names_are_public():
    for name in (
        "Solve",
        "DemagOperator",
        "ChargeGram",
        "MeshSoftIron",
        "VolSoftIron",
        "PlanarSolve",
        "SolveHysteresis",
        "EnergyStopMaterial",
        "PlayHysteresisMaterial",
        "MagnetizationSource",
        "FieldCoefficientFromSolution",
        "CoupledBody",
        "SolveCoupled",
        "FieldFromCoupledSolution",
        "CoupledHistoryBody",
        "SolveCoupledHysteresis",
        "FieldFromCoupledHysteresis",
    ):
        assert callable(getattr(vim, name))
        assert name in vim.__all__


def test_permanent_magnet_api_exposes_all_four_production_levels():
    import inspect

    assert "B_r" in inspect.signature(vim.Solve).parameters
    assert callable(vim.MagnetizationSource)       # fixed/given M
    assert callable(vim.Solve)                     # linear-recoil B_r + mu_r
    assert callable(vim.PlayHysteresisMaterial)    # simplified Play
    assert callable(vim.EnergyStopMaterial)        # full B-input EnergyStop
    assert vim.MagnetizationSource.permanent_magnet_level == 1
    assert vim.PlayHysteresisMaterial.permanent_magnet_level == 3
    assert vim.EnergyStopMaterial.permanent_magnet_level == 4


def test_legacy_hdiv_vim_names_are_not_public_aliases():
    for name in (
        "hdiv_demag_solve",
        "build_charge_gram",
        "soft_iron_from_mesh",
        "soft_iron_from_vol",
        "solve_planar_demag",
    ):
        assert not hasattr(vim, name)
        assert name not in vim.__all__
