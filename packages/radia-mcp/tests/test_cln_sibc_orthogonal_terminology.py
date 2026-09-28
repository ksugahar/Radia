from radia_mcp.radia_ngsolve.knowledge.cln_sibc_orthogonal import (
    get_cln_sibc_orthogonal_section,
)


def test_current_basis_terminology_precedes_historical_sections() -> None:
    for name in ("overview", "all"):
        text = get_cln_sibc_orthogonal_section(name)
        assert text.startswith("## Current terminology")
        assert text.count("## Current terminology") == 1
        assert "**Foster + SIBC**" in text
        assert "**CLN + SIBC**" in text
        assert "retained for compatibility" in text


def test_evrs_esim_contract_moved_to_the_esim_topic() -> None:
    assert "evrs_esim" not in get_cln_sibc_orthogonal_section("list")
    assert "Production HCurl EVRS" not in get_cln_sibc_orthogonal_section("all")
