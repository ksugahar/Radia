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
