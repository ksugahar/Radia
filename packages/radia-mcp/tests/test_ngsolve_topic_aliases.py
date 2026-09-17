"""Keep ambiguous NGSolve topic aliases deterministic and domain-specific."""

from radia_mcp.radia_ngsolve.knowledge.ngsolve import get_ngsolve_documentation
from radia_mcp.radia_ngsolve.server import ngsolve_usage


def test_usage_defaults_to_bounded_unique_topic_index():
    index = ngsolve_usage()
    assert index == get_ngsolve_documentation("index")
    assert len(index) < 12_000
    assert index.count("- coax") == 1
    assert "transmission_line" in index
    assert "NGSolve documentation topics" in index


def test_explicit_topic_and_legacy_all_remain_available():
    assert ngsolve_usage("coax") == ngsolve_usage("transmission_line")
    assert len(ngsolve_usage("all")) > len(ngsolve_usage())


def test_machine_topic_aliases_select_their_specific_guidance():
    assert "rotor-sweep" in get_ngsolve_documentation("reluctance_torque")
    assert "back-EMF" in get_ngsolve_documentation("torque_constant")
    assert "MAGNETIZING (AIR-GAP) INDUCTANCE" in get_ngsolve_documentation(
        "magnetising_inductance"
    )


def test_generic_leakage_and_drum_aliases_keep_their_existing_public_meaning():
    assert "Slot-leakage inductance" in get_ngsolve_documentation(
        "leakage_inductance"
    )
    assert "vibro-acoustic teaching lane" in get_ngsolve_documentation("drum")
    assert "magnetic circuit" in get_ngsolve_documentation(
        "magnetic_circuit_leakage"
    ).lower()
