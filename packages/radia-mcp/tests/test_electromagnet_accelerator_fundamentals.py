"""Regression tests for the accelerator textbook knowledge in electromagnet."""

from __future__ import annotations

from radia_mcp.electromagnet.accelerator_fundamentals_knowledge import (
    DOCUMENTS,
    SOURCE_CATALOG,
    get_accelerator_source_guide,
)
from radia_mcp.electromagnet.accelerator_fundamentals_knowledge import (
    TOPICS as ACCELERATOR_TOPICS,
)
from radia_mcp.electromagnet.em_knowledge import (
    TOPICS,
    get_electromagnet_documentation,
)
from radia_mcp.electromagnet.server import (
    electromagnet_accelerator_sources,
    electromagnet_usage,
    new_electromagnet_simulation,
)

EXPECTED_TOPICS = {
    "accelerator_fundamentals",
    "beam_optics_contract",
    "accelerator_magnet_types",
    "accelerator_magnet_design",
    "rapid_cycling_magnets",
    "superconducting_accelerator_magnets",
    "accelerator_magnet_measurement",
    "accelerator_model_boundaries",
    "accelerator_sources",
}


def test_accelerator_topics_are_registered_and_dispatchable() -> None:
    assert EXPECTED_TOPICS == set(ACCELERATOR_TOPICS)
    assert EXPECTED_TOPICS <= set(TOPICS)
    assert EXPECTED_TOPICS == set(DOCUMENTS)

    for topic in sorted(EXPECTED_TOPICS):
        direct = get_electromagnet_documentation(topic)
        through_tool = electromagnet_usage(topic)
        assert direct == through_tool
        assert len(direct) > 500, topic
        assert not direct.startswith("Unknown topic"), topic


def test_foundations_preserve_the_lattice_to_field_contract() -> None:
    fundamentals = electromagnet_usage("accelerator_fundamentals")
    assert "3.33564095198152" in fundamentals
    assert "theta = integral(B_perp ds) / (B rho)" in fundamentals
    assert "static magnetic" in fundamentals
    assert "does no work" in fundamentals
    assert "Classical cyclotron" in fundamentals
    assert "Synchrotron" in fundamentals
    assert "Storage ring" in fundamentals
    assert "FFAG" in fundamentals
    assert "particle species" in fundamentals
    assert "integrated normal and skew multipoles" in fundamentals

    optics = electromagnet_usage("beam_optics_contract")
    assert "beta * gamma - alpha^2 = 1" in optics
    assert "sigma_x^2 = beta_x*epsilon_x + D_x^2*sigma_delta^2" in optics
    assert "Delta Q_x,y = xi_x,y * delta" in optics
    assert "integrated multipoles" in optics
    assert "Maxwell-consistent 3-D map" in optics


def test_engineering_topics_cover_static_ramped_and_superconducting_magnets() -> None:
    design = electromagnet_usage("accelerator_magnet_design")
    assert "N I approximately B_gap * g / mu0" in design
    assert "G*r_p^2/(2*mu0)" in design
    assert "B''*r_p^3/(6*mu0)" in design
    assert "2-D body" in design
    assert "full 3-D magnet" in design
    assert "m_dot = P/(c_p Delta T)" in design

    ramped = electromagnet_usage("rapid_cycling_magnets")
    assert "skin_depth = sqrt(2/(omega*mu*sigma))" in ramped
    assert "Laminate the yoke" in ramped
    assert "stranded or transposed conductors" in ramped
    assert "stray capacitance to ground" in ramped
    assert "B-I" in ramped

    superconducting = electromagnet_usage(
        "superconducting_accelerator_magnets"
    )
    assert "J_c(B,T,strain)" in superconducting
    assert "Fine filaments" in superconducting
    assert "coupling-current time constants" in superconducting
    assert "minimum propagation zone" in superconducting
    assert "hot-spot temperature" in superconducting


def test_measurement_and_scope_fail_loud_instead_of_overclaiming() -> None:
    measurement = electromagnet_usage("accelerator_magnet_measurement")
    for instrument in ("Rotating/harmonic coil", "Hall probe", "NMR probe"):
        assert instrument in measurement
    assert "feed-down" in measurement
    assert "current history" in measurement

    boundaries = electromagnet_usage("accelerator_model_boundaries")
    for coupled_topic in (
        "RF phase stability",
        "space-charge tune shift",
        "vacuum-chamber impedance",
        "dynamic aperture",
        "injection/extraction efficiency",
    ):
        assert coupled_topic in boundaries
    assert "field solve alone" in boundaries


def test_source_guide_tracks_all_twelve_pdf_sources_without_local_paths() -> None:
    assert len(SOURCE_CATALOG) == 12
    assert sum(int(entry["pages"]) for entry in SOURCE_CATALOG) == 1844
    assert len({entry["id"] for entry in SOURCE_CATALOG}) == 12
    assert len({entry["filename"] for entry in SOURCE_CATALOG}) == 12

    guide = get_accelerator_source_guide()
    assert "12 PDFs, 1844 pages" in guide
    for filename in (
        "EddyCurrents.pdf",
        "OHO_txt-1984-Ⅲ.pdf",
        "OHO_ビーム輸送の基礎.pdf",
        "ZGOUBI_Understanding the Physics of Particle Accelerators.pdf",
        "パリティ物理学_加速器科学.pdf",
        "多芯線と導体.pdf",
    ):
        assert filename in guide

    backslash = chr(92)
    private_prefixes = (
        "W:" + backslash,
        "S:" + backslash,
        "C:" + backslash + "Users" + backslash,
    )
    for private_prefix in private_prefixes:
        assert private_prefix not in guide


def test_source_search_returns_relevant_books_and_handles_misses() -> None:
    rapid = electromagnet_accelerator_sources("rapid cycling")
    assert "Magnets for High-Intensity Proton Synchrotrons" in rapid
    assert "OHO_txt-2001-9.pdf" in rapid

    space_charge = electromagnet_accelerator_sources("space charge")
    assert "Beam Instability and Space-Charge Limits" in space_charge
    assert "Space-Charge Effects in High-Intensity Accelerators" in space_charge

    superconducting = electromagnet_accelerator_sources("superconducting")
    assert "Multifilamentary Wires and Conductors" in superconducting

    missing = electromagnet_accelerator_sources("no-such-topic-7f9c")
    assert missing.startswith("No accelerator textbook sources match query")


def test_new_simulation_prompt_starts_from_beam_requirements() -> None:
    prompt = new_electromagnet_simulation("quadrupole")
    assert "0. Translate lattice requirements through magnetic rigidity" in prompt
    assert "accelerator_fundamentals" in prompt
    assert "beam_optics_contract" in prompt
    assert "Gradient: G = dBy/dx" in prompt
