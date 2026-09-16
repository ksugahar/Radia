from radia_mcp.radia_ngsolve.server import kelvin_transformation


def test_constitutive_sampling_is_separate_from_iteration_convergence():
    text = kelvin_transformation('source_in_omega_form')
    for phrase in ('B_h-B(H_h)', 'additional interior quadrature points',
                   'not an independent A-formulation', 'rigorous field-error bound'):
        assert phrase in text


def test_mcp_boundary_guidance_distinguishes_surface_from_gauge():
    text = kelvin_transformation('source_in_omega_form')
    for phrase in ('surface_dirichlet', 'g - phi_s', 'Reduced zero is not total zero',
                   'only a BBBND point gauge', 'wrong-region',
                   'not proof of arbitrary-mesh', 'full-reassembly parity'):
        assert phrase in text
    assert 'same explicit' in text
    assert 'spatially varying source' in text
    assert 'not a\n  free-interface validation' in text
    assert 'remove internal' in text
    assert 'nonzero on the' in text


def test_cache_comparison_freezes_source_and_initial_state():
    text = kelvin_transformation('source_in_omega_form')
    for phrase in ('Reuse one mesh and one', 'source projection',
                   'disabled before coil construction', 'identical initial material states',
                   'a speedup does not establish external-solver accuracy'):
        assert phrase in text


def test_nonlinear_equivalence_requires_material_range_and_source_error_audit():
    text = kelvin_transformation('source_in_omega_form')
    for phrase in ('actual material H range', 'match the high-field tail',
                   'reject the equivalence claim', 'solved internal field',
                   'diagnostic change, not proof of a fix'):
        assert phrase in text
    assert 'not a ground-truth error bound' in text
    assert 'independent mesh/order convergence' in text
