from radia_mcp.radia_ngsolve.knowledge.esim import get_esim_documentation


def test_evrs_coupling_topic_records_production_contract() -> None:
    text = get_esim_documentation("evrs_coupling")

    assert "SurfaceImpedanceGram" in text
    assert "AssembleSurfaceImpedanceGram" in text
    assert "LocalESIMSurfaceModel" in text
    assert "BuildLocalESIMSurfaceLUT" in text
    assert "ValidateLocalESIMSurfaceLUT" in text
    assert "LocalESIMSurfaceLUT" in text
    assert "cell_solve_count = 0" in text
    assert "extrapolation is forbidden" in text
    assert "SolveLocalESIMSurfaceVIM" in text
    assert "solve_frequency_local_esim" in text
    assert "CoupledHDivHCurlLocalESIMSolution" in text
    assert "still an integration task" in text
    assert "A_bb^-1 f_b" in text
    assert "3557 active parent-HCurl DoFs" in text
    assert "25 retained eddy coordinates" in text
    assert "conductor-air/exterior" in text
    assert "Retain volume EVRS/VIM" in text


def test_full_esim_documentation_includes_evrs_coupling() -> None:
    text = get_esim_documentation("all")

    assert "Production HCurl EVRS + local ESIM-SIBC coupling" in text
    assert "local-ESIM port difference" in text


def test_evrs_coupling_topic_does_not_use_cln_naming() -> None:
    text = get_esim_documentation("evrs_coupling")

    assert "CLN" not in text
    # Generic Cauer/Dowell ladder descriptions are not CLN reduction.
