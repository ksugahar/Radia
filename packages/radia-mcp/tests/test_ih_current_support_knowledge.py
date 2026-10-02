from radia_mcp.ih.ih_knowledge import get_induction_heating_documentation
from radia_mcp.ih.esim_knowledge import ESIM_USAGE_OVERVIEW


def test_current_sibc_and_heat_contracts_are_explicit():
    text = get_induction_heating_documentation('all')
    assert 'lumped_surface_p1' in text
    assert 'Unmapped points fail loudly' in text
    assert 'nonuniform arrays fail loudly' in text.lower()
    assert 'K_gamma' in text
    assert 'peak current squared' in text
    assert 'native thermal preview is linear' in text
    assert 'must not rotate the accepted temperature a second time' in text


def test_per_panel_parser_option_is_not_advertised_as_production_support():
    assert '--esim-per-panel' in ESIM_USAGE_OVERVIEW
    assert 'not a working production coupling' in ESIM_USAGE_OVERVIEW
    assert 'strong path requires linear SIBC' in ESIM_USAGE_OVERVIEW


def test_sibc_default_is_current_architecture_not_redirect_stub():
    from radia_mcp.ih.sibc_knowledge import get_ih_sibc_documentation
    text = get_ih_sibc_documentation('peec_fem')
    assert 'BDDC_AMS_COARSE_CYCLES = 3' in text
    assert 'edge-only wirebasket' in text
    assert 'variable nodal Z_s' in text
    assert 'calc_peec.py' not in text
    assert 'v4.6.0' not in text
    assert text in get_ih_sibc_documentation()

def test_esim_is_not_described_as_unimplemented():
    from radia_mcp.ih.esim_knowledge import get_ih_esim_documentation
    text = get_induction_heating_documentation('all') + get_ih_esim_documentation('all')
    assert 'ESIM WIP' not in text
    assert 'ESIM-WIP' not in text
    assert '"esim_converged": false' in text

def test_3d_kelvin_reluctivity_factor_is_squared():
    text = get_induction_heating_documentation('all')
    assert "(r'/a)^4" not in text
    assert "nu' = nu0 * (r'/R)^2 in exterior sphere (3D HCurl A" in text
