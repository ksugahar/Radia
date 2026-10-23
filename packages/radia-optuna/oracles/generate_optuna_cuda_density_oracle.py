"""Regenerate density values from the pinned upstream Optuna, never MATLAB."""
from pathlib import Path
import json
import numpy as np
import optuna
from optuna.samplers._tpe import probability_distributions as p

assert optuna.__version__ == "5.0.0"
weights = np.array([0.0, 0.3, 0.7])
mu, sigma = np.array([-0.7, 0.2, 1.1]), np.array([0.4, 0.8, 1.3])
cases = []
for name, cls, low, high, step, samples in [
    ("continuous", p._BatchedTruncNormDistributions, -2, 3, None, [-2, -.5, 0, 1, 3]),
    ("log", p._BatchedTruncLogNormDistributions, .1, 10, None, [.1, .5, 1, 4, 10]),
    ("discrete", p._BatchedDiscreteTruncNormDistributions, -2, 3, 1, [-2, -1, 0, 1, 3]),
    ("log_discrete", p._BatchedDiscreteTruncLogNormDistributions, 1, 9, 1, [1, 3, 5, 7, 9]),
]:
    distribution = cls(mu, sigma, low, high) if step is None else cls(mu, sigma, low, high, step)
    expected = p._MixtureOfProductDistribution(weights, [distribution]).log_pdf(np.array(samples)[:, None])
    cases.append(dict(name=name, samples=samples, weights=weights.tolist(), mu=mu.tolist(),
        sigma=sigma.tolist(), low=distribution.adapted_low, high=distribution.adapted_high,
        log_scale=distribution.is_log, step=step, expected=expected.tolist()))
probabilities = np.array([[.2, .3, .5], [.1, .6, .3], [.7, .2, .1]])
samples = np.array([0, 1, 2, 0, 2])
expected = p._MixtureOfProductDistribution(weights, [p._BatchedCategoricalDistributions(probabilities)]).log_pdf(samples[:, None])
cases.append(dict(name="categorical", samples=(samples+1).tolist(), weights=weights.tolist(),
    probabilities=probabilities.tolist(), expected=expected.tolist()))
for entry in cases:
    for key, value in dict(mu=[], sigma=[], low=0, high=1, log_scale=False, step=None, probabilities=[]).items():
        entry.setdefault(key, value)
Path(__file__).with_name("optuna_cuda_density_oracle.json").write_text(json.dumps(dict(
    schema="radia.optuna.cuda-density-oracle.v1", optuna_version=optuna.__version__, cases=cases),
    indent=2)+"\n", encoding="utf-8")
