from radia_mcp.radia_ngsolve.server import kelvin_transformation


def test_mcp_boundary_guidance_distinguishes_surface_from_gauge():
    text = kelvin_transformation('source_in_omega_form')
    for phrase in ('surface_dirichlet', 'g - phi_s', 'Reduced zero is not total zero',
                   'only a BBBND point gauge', 'wrong-region',
                   'not proof of arbitrary-mesh', 'full-reassembly parity'):
        assert phrase in text
    assert 'same explicit' in text
    assert 'spatially varying source' in text
