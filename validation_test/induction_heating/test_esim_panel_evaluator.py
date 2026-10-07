"""Interpolation contracts; final certification is deliberately never cached."""
import numpy as np
import pytest

from radia.esim_panel_evaluator import PanelESIMEvaluator


class Cell:
    def __init__(self, law=None):
        self.calls = []
        self.law = law or (lambda h: 2 + .1*np.log(h) + 1j*(3 + .2*np.log(h)))

    def solve(self, h):
        self.calls.append(h)
        return dict(Z=self.law(h), converged=True)


def evaluator(cell, **kwargs):
    return PanelESIMEvaluator(cell, configuration=lambda: {'law': 'test'}, **kwargs)


@pytest.mark.parametrize('mode', ['direct', 'table'])
def test_all_panels_directly_certified_even_duplicates(mode):
    cell = Cell()
    ev = evaluator(cell, mode=mode)
    fields = np.array([.001, .02, 1., 1., 10.])
    expected = np.array([cell.law(h) for h in fields])
    np.testing.assert_allclose(ev.evaluate(fields), expected, rtol=1e-14)
    before = len(cell.calls)
    np.testing.assert_allclose(ev.certify(fields), expected, rtol=1e-14)
    assert len(cell.calls) - before == len(fields)
    assert ev.diagnostics()['certification_cell_calls'] == len(fields)


def test_table_cache_reuse_and_floor():
    cell = Cell()
    ev = evaluator(cell, mode='table')
    ev.evaluate([0, .0001])
    assert cell.calls == [.001]
    ev.evaluate([.001])
    assert len(cell.calls) == 1


def test_table_range_extension_no_extrapolation_and_component_errors():
    cell = Cell(lambda h: (1+.1*np.log1p(h)**2)*(1+2j))
    ev = evaluator(cell, mode='table', interpolation_tol=1e-4)
    for fields in (np.geomspace(.01, 10, 23), np.geomspace(.001, 100, 43)):
        actual = ev.evaluate(fields)
        exact = np.array([cell.law(h) for h in fields])
        assert np.max(np.abs(actual-exact)/np.abs(exact)) < 2e-4
        assert ev.knots[0] <= min(fields) and ev.knots[-1] >= max(fields)
        assert np.all(actual.real >= 0)
    assert ev.diagnostics()['max_accepted_probe_error'] <= ev.tolerance


def test_complete_identity_mutation_rejected():
    state = {'sigma': 1}
    ev = PanelESIMEvaluator(Cell(), mode='table', configuration=lambda: state)
    ev.evaluate([1.])
    state['sigma'] = 2
    with pytest.raises(RuntimeError, match='configuration changed'):
        ev.certify([1.])


@pytest.mark.parametrize('fields', [[-1], [np.nan], [np.inf], [[1]]])
def test_invalid_fields(fields):
    with pytest.raises(ValueError, match='finite nonnegative 1D'):
        evaluator(Cell()).evaluate(fields)


@pytest.mark.parametrize('value', [-1+1j, complex(np.nan, 1)])
def test_bad_direct_cells_rejected(value):
    with pytest.raises(RuntimeError, match='nonfinite or nonpassive'):
        evaluator(Cell(lambda h: value)).certify([1.])


def test_unconverged_direct_cell_rejected():
    class Failed:
        def solve(self, h):
            return dict(Z=1j, converged=False)
    with pytest.raises(RuntimeError, match='converg'):
        evaluator(Failed()).evaluate([1.])


def test_table_budget_fails_loudly():
    ev = evaluator(Cell(lambda h: 1+np.log(h)**2+1j), mode='table',
                   interpolation_tol=1e-10, max_table_cells=3)
    with pytest.raises(RuntimeError, match='budget exhausted'):
        ev.evaluate([.1, 1.])


def test_empirical_probe_is_not_used_as_final_certification():
    # Endpoint/midpoint sampling misses this passive interior oscillation.
    cell = Cell(lambda h: 1 + .1*np.sin(2*np.pi*np.log(h)) + 1j)
    ev = evaluator(cell, mode='table', interpolation_tol=1e-4)
    fields = np.exp([0., .25, 1.])
    interpolated = ev.evaluate(fields)
    certified = ev.certify(fields)
    assert np.max(abs(interpolated-certified)/abs(certified)) > .01
    assert ev.diagnostics()['certification_cell_calls'] == len(fields)


def test_real_finite_cell_configuration_and_constant_reduction():
    from radia.esim_cell_problem import ESIMFiniteSlabSolver
    cell = ESIMFiniteSlabSolver(.003, sigma=5.8e7, frequency=1000,
                               mu_r=1., geometry='cylinder')
    ev = PanelESIMEvaluator(cell, mode='table', interpolation_tol=1e-4)
    fields = np.geomspace(.001, 1e4, 31)
    np.testing.assert_allclose(ev.evaluate(fields), ev.certify(fields), rtol=1e-10)
    cell.sigma *= 2
    with pytest.raises(RuntimeError, match='configuration changed'):
        ev.evaluate(fields)


@pytest.mark.parametrize('kwargs', [dict(mode='unknown'), dict(interpolation_tol=0),
    dict(interpolation_tol=np.nan), dict(max_table_cells=2), dict(max_table_cells=True)])
def test_invalid_controls(kwargs):
    with pytest.raises(ValueError):
        evaluator(Cell(), **kwargs)


@pytest.mark.parametrize('mode', ['direct', 'table'])
def test_missing_production_cell_identity_names_attributes(mode):
    with pytest.raises(ValueError, match='missing attributes: half_thickness'):
        PanelESIMEvaluator(Cell(), mode=mode)
