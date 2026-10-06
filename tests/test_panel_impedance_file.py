import json

import numpy as np
import pytest

from radia.surface_impedance import read_panel_impedance, write_panel_impedance


@pytest.fixture
def surface():
    import ngsolve as ng
    from netgen.occ import Box, OCCGeometry, Pnt, Glue
    body = Box(Pnt(0, 0, 0), Pnt(.01, .01, .01))
    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(Glue(list(body.faces))).GenerateMesh(maxh=.006))


def test_file_roundtrip_and_rejects_wrong_frequency(surface, tmp_path):
    import ngsolve as ng
    values = np.linspace(1, 2, surface.GetNE(ng.BND))*(.001+.002j)
    path = tmp_path/'zs.json'
    write_panel_impedance(path, surface, values, frequency_hz=1000)
    np.testing.assert_array_equal(read_panel_impedance(path, surface, frequency_hz=1000).values, values)
    with pytest.raises(ValueError, match='frequency mismatch'):
        read_panel_impedance(path, surface, frequency_hz=2000)


@pytest.mark.parametrize('change', ['mesh', 'shape', 'active', 'nan', 'unit'])
def test_bad_or_stale_panel_file_fails(surface, tmp_path, change):
    import ngsolve as ng
    path = tmp_path/'zs.json'
    write_panel_impedance(path, surface, np.ones(surface.GetNE(ng.BND))*(1+1j), frequency_hz=1000)
    data = json.loads(path.read_text())
    if change == 'mesh':
        data['surface_sha256'] = 'wrong'
    elif change == 'shape':
        data['imag_ohm'] = data['imag_ohm'][:-1]
    elif change == 'active':
        data['real_ohm'][0] = -1
    elif change == 'nan':
        data['real_ohm'][0] = float('nan')
    else:
        data['unit'] = 'milliohm'
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        read_panel_impedance(path, surface, frequency_hz=1000)
