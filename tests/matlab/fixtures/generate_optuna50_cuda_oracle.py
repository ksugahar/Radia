"""Generate CUDA density/selection truth from the pinned upstream implementation."""
from pathlib import Path
import hashlib
import json

import numpy as np
import optuna
from optuna.samplers._tpe import probability_distributions as pd


def generate():
    if optuna.__version__ != "5.0.0":
        raise RuntimeError("This fixture requires Optuna 5.0.0")
    below_weights = np.array([0.2, 0.3, 0.5])
    above_weights = np.array([0.4, 0.6])

    def numerical(kind, upper):
        mu = np.array([-0.5, 0.7] if upper else [-1.0, 0.0, 1.0])
        sigma = np.array([0.7, 1.1] if upper else [0.6, 0.9, 1.3])
        if kind == "continuous":
            return pd._BatchedTruncNormDistributions(mu, sigma, -3.0, 3.0)
        if kind == "log":
            return pd._BatchedTruncLogNormDistributions(mu, sigma, np.exp(-3.0), np.exp(3.0))
        if kind == "discrete":
            return pd._BatchedDiscreteTruncNormDistributions(mu, sigma, -3.0, 3.0, 1.0)
        if kind == "log_discrete":
            return pd._BatchedDiscreteTruncLogNormDistributions(mu + 1.5, sigma, 1.0, 16.0, 1.0)
        return pd._BatchedCategoricalDistributions(np.array(
            [[0.1, 0.7, 0.2], [0.6, 0.2, 0.2]] if upper else
            [[0.2, 0.3, 0.5], [0.1, 0.7, 0.2], [0.8, 0.1, 0.1]]))

    def descriptor(d, weights):
        result = {"weights": weights.tolist(), "categorical": isinstance(d, pd._BatchedCategoricalDistributions)}
        if result["categorical"]:
            result["probabilities"] = d.weights.tolist()
        else:
            result.update(mu=d.mu.tolist(), sigma=d.sigma.tolist(),
                          internal_low=float(d.adapted_low), internal_high=float(d.adapted_high),
                          low=float(d.low), high=float(d.high), log=bool(d.is_log),
                          step=float(d.step) if d.step else None)
        return result

    cases = []
    for case_index, kinds in enumerate([[k] for k in ["continuous", "log", "discrete", "log_discrete", "categorical"]]
                                       + [["continuous", "categorical", "log", "discrete", "log_discrete"]]):
        lower = pd._MixtureOfProductDistribution(below_weights, [numerical(k, False) for k in kinds])
        upper = pd._MixtureOfProductDistribution(above_weights, [numerical(k, True) for k in kinds])
        seed, count = 201 + case_index, 64
        rng = np.random.RandomState(seed)
        samples = lower.sample(rng, count)
        below_pdf, above_pdf = lower.log_pdf(samples), upper.log_pdf(samples)
        best = int(np.argmax(below_pdf - above_pdf))
        matlab_samples = samples.copy()
        for dim, kind in enumerate(kinds):
            if kind == "categorical":
                matlab_samples[:, dim] += 1
        cases.append(dict(name="_".join(kinds), seed=seed, count=count,
                          below=[descriptor(d, below_weights) for d in lower.distributions],
                          above=[descriptor(d, above_weights) for d in upper.distributions],
                          samples=matlab_samples.tolist(), below_pdf=below_pdf.tolist(), above_pdf=above_pdf.tolist(),
                          best_index=best + 1, best=matlab_samples[best].tolist(), next_random=rng.rand(8).tolist()))
    return dict(schema="radia.optuna.cuda-oracle.v1", optuna_version=optuna.__version__,
                numpy_version=np.__version__, generator_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
                cases=cases)


if __name__ == "__main__":
    Path(__file__).with_name("optuna50_cuda_oracle.json").write_text(
        json.dumps(generate(), indent=2, allow_nan=False) + "\n", encoding="utf-8")
