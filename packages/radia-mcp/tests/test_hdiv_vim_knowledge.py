from radia_mcp.radia_ngsolve.knowledge.hdiv_vim import (
    get_hdiv_vim_documentation,
)


def test_hdiv_vim_knowledge_distinguishes_fixed_and_evolving_pm_paths():
    implementation = " ".join(get_hdiv_vim_documentation("implementation").split())
    assert "MagnetizationSource" in implementation
    assert "fixed-M source" in implementation
    assert "EnergyStopMaterial" in implementation
    assert "initial_b_path" in implementation
    assert "initial_state" in implementation
    assert "mutually evolving PM plus nonlinear iron" in implementation


def test_hdiv_vim_knowledge_exposes_the_four_level_permanent_magnet_ladder():
    implementation = " ".join(get_hdiv_vim_documentation("implementation").split())
    assert "four-level model ladder" in implementation
    assert "MagnetizationSource(mesh, M_given)" in implementation
    assert "mu_r=mu_rec, B_r=B_r" in implementation
    assert "simplified Play" in implementation
    assert "full B-input EnergyStop" in implementation
    assert "rigid `mu_rec=1` limit" in implementation


def test_hdiv_vim_knowledge_lists_irreversible_pm_validation_gates():
    verification = " ".join(get_hdiv_vim_documentation("verification").split())
    assert "positive-gamma proximal" in verification
    assert "stationarity" in verification
    assert "reverse-field remanence loss" in verification
    assert "split-run restart parity" in verification
