from radia_mcp.radia_ngsolve.knowledge.hdiv_vim import (
    get_hdiv_vim_documentation,
)


def test_hdiv_vim_knowledge_distinguishes_fixed_and_evolving_pm_paths():
    implementation = " ".join(get_hdiv_vim_documentation("implementation").split())
    assert "MagnetizationSource" in implementation
    assert "fixed-M source" in implementation
    assert "EnergyStopMaterial" in implementation
    assert "HDivSolver" in implementation
    assert "initial_b_path" in implementation
    assert "initial_state" in implementation
    assert "SolveCoupledHysteresis([history_pm, ...]" in implementation
    assert "all states commit together" in implementation


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


def test_hdiv_vim_knowledge_identifies_the_ngsolve_hdiv_family():
    family = " ".join(get_hdiv_vim_documentation("family").split())
    overview = " ".join(get_hdiv_vim_documentation("overview").split())

    assert "HDiv(mesh, order=p)" in family
    assert "BDM, the NGSolve default" in family
    assert "HDiv(mesh, order=p, RT=True)" in family
    assert "Raviart--Thomas" in family
    assert "BDM1/BDM2" in overview
    assert "requires `RT=True`" in overview


def test_hdiv_vim_knowledge_exposes_bdm_eddy_bubble_production_path():
    coupled = " ".join(get_hdiv_vim_documentation("eddy_bubble").split())

    assert "NgsolveBDMEddyBubbleVIM" in coupled
    assert "NgsolveBDMHDivMMMResponseReduction" in coupled
    assert 'parent_family="BDM"' in coupled
    assert "Only conductor faces touching air are SIBC faces" in coupled
    assert "eddy_flux_density(points)" in coupled


def test_hdiv_vim_knowledge_exposes_h1_hodge_mixed_validation_boundary():
    coupled = " ".join(get_hdiv_vim_documentation("eddy_bubble").split())
    implementation = " ".join(
        get_hdiv_vim_documentation("implementation").split()
    )
    verification = " ".join(
        get_hdiv_vim_documentation("verification").split()
    )

    assert "H1HodgeDemagOperator" in coupled
    assert "hdiv.FreeDofs()" in coupled
    assert "snapshot residual below `3e-12`" in coupled
    assert "does not repair the open-boundary charge kernel" in coupled
    assert "hdiv_definedon" in coupled
    assert "demag_operator_factory" in coupled
    assert "exact restricted BDM space" in coupled
    assert "standard unit H1 stiffness metric" in coupled
    assert "weighted or Kelvin-transformed metric" in coupled
    assert "1.0916977441" in coupled
    assert "4.65%" in coupled
    assert "0.98%" in coupled
    assert "strict h/p ladder" in coupled
    assert "H1HodgeDemagOperator" in implementation
    assert "C.T K^-1 C" in implementation
    assert "material-restricted HDiv spaces" in verification


def test_hdiv_vim_knowledge_rejects_single_port_as_a_completeness_proof():
    coupled = " ".join(get_hdiv_vim_documentation("eddy_bubble").split())

    assert "h-refinement as proof" in coupled
    assert "single uniform-field vector-potential port" in coupled
    assert "8 and 16 Krylov steps" in coupled
    assert "A, r^2 A, r^4 A, z^2 A" in coupled
    assert "select the mode" in coupled
    assert "port residue" in coupled
    assert "largest generalized eigenvalue" in coupled
    assert "both h- and p-refinement" in coupled


def test_hdiv_vim_knowledge_distinguishes_hdiv_star_deficiency_from_hcurl_current():
    coupled = " ".join(get_hdiv_vim_documentation("eddy_bubble").split())

    assert "HDiv local-response completeness gate" in coupled
    assert "not extra current modes" in coupled
    assert "J=curl(T)" in coupled
    assert "NgsolveHDivLocalPolynomialTrainingPorts" in coupled
    assert "30 normalized vector-polynomial probes" in coupled
    assert "training-only" in coupled
    assert "rejects them as physical `external_fields`" in coupled
    assert "degree two against degree three" in coupled
    assert "within 1.02% in Joule loss" in coupled
    assert "universal accuracy claim" in coupled


def test_hdiv_vim_knowledge_requires_frequency_appropriate_sibc_routing():
    coupled = " ".join(get_hdiv_vim_documentation("eddy_bubble").split())

    assert "EddySIBCApplicability" in coupled
    assert "assembled model identity" in coupled
    assert "within the same volumetric or SIBC regime" in coupled
    assert "solve_frequency` rejects" in coupled
    assert "Rebuild the topology-aware system" in coupled
    assert "invalid extrapolation" in coupled
    assert "machine precision" in coupled


def test_hdiv_vim_knowledge_exposes_coupled_local_esim_without_overclaiming():
    coupled = " ".join(get_hdiv_vim_documentation("eddy_bubble").split())

    assert "SIBC basis does not by itself imply ESIM" in coupled
    assert "solve_frequency_local_esim" in coupled
    assert "CoupledHDivHCurlLocalESIMSolution" in coupled
    assert "exactly one physical excitation" in coupled
    assert "50/1000/5000 A/m" in coupled
    assert "fixed-Gram replay to machine precision" in coupled
    assert "nonlinear **skin** impedance" in coupled
    assert "bulk nonlinear B-H operator" in coupled
    assert "is not yet implemented" in coupled


def test_hdiv_vim_knowledge_bounds_magnetic_conductor_claims():
    coupled = " ".join(get_hdiv_vim_documentation("eddy_bubble").split())

    assert "Thin magnetic-conductor adjudication" in coupled
    assert "axisymmetric Q2 volume solve" in coupled
    assert "full 3-D HCurl A-form" in coupled
    assert "mapped-HEX BDM1" in coupled
    assert "BDM2 HEX ladder is non-monotone" in coupled
    assert "not an accuracy oracle" in coupled
    assert "5.38% to 3.34% to 1.90%" in coupled
    assert "sampled-interaction" in coupled
    assert "native Q2 isoparametric HEX Gram" in coupled
    assert "hex_geometry_backend" in coupled
    assert "splitting a trilinear HEX" in coupled
    assert "strict 32/96/384-HEX h ladder" in coupled
    assert "8.03e-7" in coupled
    assert "matrix-free restricted operator" in coupled
    assert "1819.7 s to 339.8 s" in coupled
    assert "18432-element fine axisymmetric Q2 reference" in coupled
    assert "2.19%" in coupled
    assert "frozen reject" in coupled
    assert "universal solver-superiority claim" in coupled
