# MATLAB Optuna MOTPE Parzen migration

## Retired implementation

The first MATLAB `MOTPESampler` used one scalar Gaussian KDE bandwidth per
parameter and formed its candidate set from every good observation plus uniform
random points.  That path was deterministic and useful while the Study/Trial
compatibility layer was being built, but it was not Optuna's MOTPE estimator.

On the five-seed ZDT1 quality lane (120 trials), the retired path produced a
median Pareto-front error of 0.316553.  The result exposed a problem in the
estimator and candidate-generation semantics, not in TPE as an algorithm.

## Canonical implementation

MOTPE now uses the shared `ParzenEstimator` used by single-objective TPE:

- hypervolume-contribution weights are attached to below observations;
- the prior remains an explicit mixture component;
- candidates are sampled from the below estimator;
- selection maximizes `log l(x) - log g(x)`;
- categorical parameters use the same weighted mixture contract.

The scalar KDE helpers are deleted.  Do not restore them as a fallback.  New
MOTPE work should improve the shared Parzen path and retain the deterministic
random/TPE/NSGA-II validation comparison.

## Multivariate note

Joint TPE uses a shared mixture component and Optuna's multivariate bandwidth
rule.  The table schema now persists versioned numeric and categorical
distribution metadata, so ordinary scalar suggestions participate in the
automatically inferred intersection search space when `Multivariate=true`.
Conditional or changed distributions stay outside that intersection and use
independent sampling.  `Trial.suggestVector` remains only as an explicit
numeric convenience; it is no longer required for joint sampling.
