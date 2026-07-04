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
        "ChargeGramGauss",
        "MeshSoftIron",
        "VolSoftIron",
        "PlanarSolve",
    ):
        assert callable(getattr(vim, name))
        assert name in vim.__all__


def test_legacy_hdiv_vim_names_are_not_public_aliases():
    for name in (
        "hdiv_demag_solve",
        "build_charge_gram",
        "build_charge_gauss",
        "soft_iron_from_mesh",
        "soft_iron_from_vol",
        "solve_planar_demag",
    ):
        assert not hasattr(vim, name)
        assert name not in vim.__all__
