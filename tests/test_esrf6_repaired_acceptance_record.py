"""Recompute the committed nominal ESRF6 acceptance, without a solver runtime."""
import hashlib
import json
from pathlib import Path

import numpy as np


def test_ci_wheel_nominal_three_engine_record():
    directory = (Path(__file__).resolve().parents[1] / 'validation_test' /
                 'esrf_three_engine/results/candidate_59b094d8')
    path = directory / 'three_engine_case6_bdm1_bonus12.json'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == (
        'b00ba5d7a6532721e23eefe6e8b1c0b8889b9f52e1f4240bbd2f2a99570d722f')
    result = json.loads(path.read_bytes())
    assert result['passed'] is True and result['nonlinear_converged'] is True
    runtime = result['runtime_identity']
    assert runtime['editable'] is False and runtime['installed_import'] is True
    assert runtime['direct_url']['archive_info']['hashes']['sha256'] == (
        '46e3aa1ddd89010e21419e1d28f6e44e403cf95014d7e5486109eca4232145b9')
    assert result['implementation_sha256']['radia._radia_pybind'] == (
        'cbc9bbbc61b5cd0ee7035c642197522b2615a871e3d4c49c2887184f4f9040ec')
    hdiv = result['engines']['hdiv_mmm']['nonlinear_stats']
    assert hdiv['nonlinear_line_search_exhausted'] is False
    assert hdiv['nonlinear_converged_final_stage'] is True
    assert 0 <= hdiv['nonlinear_final_relative_residual'] <= 2e-5
    for name, key in [('reduced_a', 'final_relative_change'),
                      ('mixed_total_reduced_omega', 'relative_B_change')]:
        stats = result['engines'][name]['nonlinear_stats']
        assert stats['converged'] is True
        assert 0 <= stats[key] <= 2e-5
    assert result['provenance']['engine_settings']['mixed_total_reduced_omega']['bonus_intorder'] == 12
    points = np.asarray(result['observation_points_m'])
    assert points.shape == (45, 3) and np.isfinite(points).all()
    names = ('hdiv_mmm', 'reduced_a', 'mixed_total_reduced_omega')
    assert set(result['fields_T']) == set(names)
    for field in result['fields_T'].values():
        values = np.asarray(field)
        assert values.shape == (45, 3) and np.isfinite(values).all()
    assert set(result['pairwise_core']) == {
        names[0] + '__vs__' + names[1],
        names[0] + '__vs__' + names[2],
        names[1] + '__vs__' + names[2],
    }
    core = abs(points[:, 0]) <= .02
    assert core.sum() == 27
    errors = []
    for pair, recorded in result['pairwise_core'].items():
        left, right = pair.split('__vs__')
        a, b = (np.asarray(result['fields_T'][name])[core] for name in (left, right))
        error = np.linalg.norm(a-b) / np.linalg.norm(a)
        assert np.isfinite(error) and error <= .03
        np.testing.assert_allclose(error, recorded['relative_rms'], rtol=1e-13)
        errors.append(error)
    assert len(errors) == 3
    np.testing.assert_allclose(max(errors), result['maximum_core_pairwise_relative_rms'], rtol=1e-13)


def test_repaired_bdm1_nominal_three_engine_record():
    directory = (Path(__file__).resolve().parents[1] / 'validation_test' /
                 'esrf_three_engine/results/candidate_4d85e72cc')
    hdiv_path = directory / 'case6_bdm1_mass_riesz.json'
    hdiv = json.loads(hdiv_path.read_bytes())
    result = json.loads((directory / 'case6_fem_repaired_acceptance.json').read_bytes())
    assert result['completed'] is True and result['passed'] is True
    assert result['hdiv_result_sha256'] == hashlib.sha256(hdiv_path.read_bytes()).hexdigest()
    assert result['implementation'] == hdiv['implementation']
    assert result['source'] == hdiv['source']
    np.testing.assert_array_equal(result['bh_table'], hdiv['bh_table'])
    assert result['fem_mesh_sha256'] == (
        'dfc12b84f80db1fa17fb5012b6f87c072e8aa1fb3c085b61269f047a879a54ac')
    assert result['observation_half_width_m'] == hdiv['observation_half_width_m']
    assert result['threads'] == 8 and result['fem_order'] == 2
    assert result['nonlinear_tolerance'] == 2e-5
    assert hdiv['options']['order'] == 1
    assert hdiv['options']['preconditioner'] == 'mass-riesz'
    assert hdiv['options']['image'] is None
    stats = hdiv['nonlinear_stats']
    residual = stats['nonlinear_final_relative_residual']
    assert np.isfinite(residual) and 0 <= residual <= hdiv['options']['nl_tol']
    assert stats['nonlinear_converged_final_stage'] is True
    assert hdiv['completed'] is True and hdiv['accepted'] is True
    points = np.asarray(result['observation_points_m'])
    np.testing.assert_array_equal(points, hdiv['observation_points_m'])
    selector = abs(points[:, 0]) <= .02
    assert np.count_nonzero(selector) == 27
    names = ('hdiv_mmm', 'reduced_a', 'mixed_total_reduced_omega')
    assert set(result['fields_T']) == set(names)
    fields = {name: np.asarray(result['fields_T'][name]) for name in names}
    np.testing.assert_array_equal(fields['hdiv_mmm'], hdiv['B_T'])
    for field in fields.values():
        assert field.shape == (45, 3) and np.all(np.isfinite(field))
    for name in names[1:]:
        fem_stats = result['engines'][name]['nonlinear_stats']
        assert fem_stats['converged'] is True
        # FEM Picard changes are not energy-Newton equation residuals.
        key = 'final_relative_change' if name == 'reduced_a' else 'relative_B_change'
        change = fem_stats[key]
        assert np.isfinite(change) and 0 <= change <= fem_stats['tolerance']
    errors = []
    for index, left in enumerate(names):
        for right in names[index + 1:]:
            a, b = fields[left][selector], fields[right][selector]
            error = np.linalg.norm(a - b) / np.linalg.norm(a)
            assert np.isfinite(error) and error <= .03
            recorded = result['pairwise_core'][left + '__vs__' + right]['relative_rms']
            np.testing.assert_allclose(recorded, error, rtol=1e-13, atol=0)
            errors.append(error)
    np.testing.assert_allclose(result['maximum_core_pairwise_relative_rms'], max(errors),
                               rtol=1e-13, atol=0)
