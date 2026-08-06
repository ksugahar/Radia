from radia_mcp.radia_ngsolve.knowledge.cln_sibc_orthogonal import (
    get_cln_sibc_orthogonal_documentation,
    get_cln_sibc_orthogonal_section,
)


def test_evrs_esim_section_records_production_contract() -> None:
    text = get_cln_sibc_orthogonal_section("evrs_esim")

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


def test_full_cln_sibc_documentation_includes_evrs_esim() -> None:
    text = get_cln_sibc_orthogonal_documentation()

    assert "Production HCurl EVRS + local ESIM-SIBC mixed Galerkin" in text
    assert "local-ESIM port difference" in text


def test_evrs_esim_section_is_listed() -> None:
    assert "evrs_esim" in get_cln_sibc_orthogonal_section("list")
