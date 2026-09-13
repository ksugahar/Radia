"""Recompute the committed nominal ESRF6 acceptance, without a solver runtime."""
import hashlib
import json
from pathlib import Path

import numpy as np


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
