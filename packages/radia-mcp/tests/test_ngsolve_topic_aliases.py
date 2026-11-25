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


def test_direct_solver_guidance_is_explicit_and_keeps_mkl_scope_separate():
    doc = get_ngsolve_documentation("solvers")
    assert 'inverse="sparsecholesky"' in doc
    assert 'inverse="pardiso"' not in doc
    assert 'inverse="umfpack"' not in doc
    assert "no automatic backend fallback" in doc
    assert "MATLAB MEX" in doc
    assert "BLAS/LAPACK and FFT" in doc
    assert "true-residual validation" in doc

def test_consolidated_fem_usage_and_inline_ih_recipe():
    from ngsolve_mcp.gauging.server import gauging_usage
    from ngsolve_mcp.kelvin.server import kelvin_usage
    from ngsolve_mcp.weak_form.server import weak_form_usage
    from radia_mcp.fem.server import fem_gauge_open_boundary
    from radia_mcp.ih.server import ih_sibc
    from radia_mcp.radia_ngsolve.server import dtn_coarse_mesh, kelvin_transformation

    recipe = ih_sibc("ngsolve_recipes")
    for topic in ("eddy_current", "sibc_robin", "ih_reference",
                  "magnetostatics_a", "symbolic_bfi", "hdiv_flux"):
        assert ngsolve_usage(topic) == weak_form_usage(topic)
    for topic in ("eddy_current", "sibc_robin"):
        assert weak_form_usage(topic) in recipe
    verify_first = kelvin_transformation("verify_first")
    assert verify_first == kelvin_usage("verify_first")
    assert verify_first in recipe
    coil = fem_gauge_open_boundary("coil_av_tree")
    assert gauging_usage("coil_av_tree") in coil
    assert gauging_usage("coil_av_tree") in recipe
    assert "within 1%" in coil and "30% lower" in coil
    assert all(f"# Check {number}:" in recipe for number in (1, 4, 5))
    assert "Unknown topic" not in dtn_coarse_mesh("spherical_harmonics")
    assert "Unknown topic" not in ngsolve_usage("taskmanager")
