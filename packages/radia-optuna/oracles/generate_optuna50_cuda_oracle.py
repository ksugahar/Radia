"""Generate CUDA density/selection truth from the pinned upstream implementation."""
from pathlib import Path
import hashlib
import json

import numpy as np
import optuna
from optuna.samplers._tpe import probability_distributions as pd
from optuna.samplers._tpe.parzen_estimator import _ParzenEstimator, _ParzenEstimatorParameters
from optuna.samplers._tpe.sampler import default_weights


def model_cases():
    space = {"x": optuna.distributions.FloatDistribution(-2, 2),
             "log": optuna.distributions.FloatDistribution(float(np.exp(-2)), float(np.exp(2)), log=True),
             "integer": optuna.distributions.IntDistribution(0, 4096),
             "category": optuna.distributions.CategoricalDistribution(["a", "b", "c"])}
    cases = []
    for n in [0, 1, 2, 7, 33, 257, 513]:
        for endpoints in [False, True]:
            seed = 700 + len(cases)
            def observations(offset):
                order = np.random.RandomState(seed + offset).permutation(n)
                return {"x": np.linspace(-1.5, 1.5, n)[order],
                        "log": np.exp(np.linspace(-1.3, 1.1, n))[order],
                        "integer": (order % 4 if n <= 33 else order).astype(float),
                        "category": (order % 3).astype(float)}
            good, bad = observations(0), observations(1)
            bad["x"] += 0.2
            bad["log"] *= np.exp(0.1)
            bad["integer"] += 5 if n <= 33 else 1024
            bad["category"] = (bad["category"] + 1) % 3
            parameters = _ParzenEstimatorParameters(1.0, True, endpoints, default_weights, True)
            lower, upper = (_ParzenEstimator(o, space, parameters) for o in [good, bad])
            rng = np.random.RandomState(seed)
            samples = lower.sample(rng, 24)
            selected = int(np.argmax(lower.log_pdf(samples) - upper.log_pdf(samples)))
            def matrix(o):
                result = np.column_stack(list(o.values()))
                result[:, -1] += 1
                return result.tolist()
            cases.append(dict(name=f"rows_{n}_endpoints_{endpoints}", seed=seed,
                              endpoints=endpoints, good=matrix(good), bad=matrix(bad),
                              lows=[d.low if hasattr(d, "low") else 0 for d in space.values()],
                              highs=[d.high if hasattr(d, "high") else 0 for d in space.values()],
                              best=[float(samples[k][selected]) + (k == "category") for k in space],
                              next_random=rng.rand(8).tolist()))
    return cases


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
                cases=cases, model_cases=model_cases())


if __name__ == "__main__":
    Path(__file__).with_name("optuna50_cuda_oracle.json").write_text(
        json.dumps(generate(), indent=2, allow_nan=False) + "\n", encoding="utf-8")
