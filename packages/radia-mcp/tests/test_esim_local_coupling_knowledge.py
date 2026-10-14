from radia_mcp.radia_ngsolve.knowledge.esim import get_esim_documentation


def test_esim_knowledge_prefers_local_gram_and_bounds_bulk_nonlinearity():
    overview = " ".join(get_esim_documentation("overview").split())
    iteration = " ".join(get_esim_documentation("karl_iteration").split())
    api = " ".join(get_esim_documentation("module_api").split())

    assert "Skin depth" in overview or "Skin depth" in overview.capitalize()
    assert "LocalESIMSurfaceModel" in iteration
    assert "AssembleSurfaceImpedanceGram" in iteration
    assert "solve_frequency_local_esim" in iteration
    assert "Scalar legacy" in iteration
    assert "simultaneous ordinary bulk nonlinear B-H" in iteration
    assert "solve_frequency_local_esim" in api
    assert "exactly one physical excitation" in api
    assert "not yet a simultaneous bulk nonlinear B-H solve" in api
