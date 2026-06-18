"""Global, derivative-free population optimization -- Differential Evolution (Storn-Price
DE/rand/1/bin), the GLOBAL complement to the package's LOCAL optimizers (Nelder-Mead in the
`outer_loop` knowledge, Levenberg-Marquardt in :mod:`nonlinear_lsq`) and to the LINEAR inverse
solvers (:mod:`linear_inverse`). Use it when the objective is MULTIMODAL / non-convex (many local
minima) so a descent method started anywhere would stall -- e.g. magnet/coil layout objectives
with competing harmonics, or any boxed design space.

DE/rand/1/bin: for each member x_i form a mutant v = a + F (b - c) from three distinct random
members, binomially cross v with x_i (rate CR, at least one gene forced), and keep the trial only
if it lowers the objective. The population drifts to the global basin without gradients.

SCALING CAVEAT (learned): the budget must grow with dimension -- popsize~15/dim & a few hundred
generations solve 2-D multimodal problems, but a 5-D Rastrigin needs popsize~25/dim and ~2000
generations to reliably reach the global optimum (smaller budgets stall one Rastrigin "ring" away).
"""
import numpy as np


def differential_evolution(func, bounds, popsize=15, F=0.7, CR=0.9, maxiter=400,
                           tol=1e-12, seed=0):
    """Minimize ``func(x)`` over the box ``bounds`` (list of (lo, hi) per dimension) by Differential
    Evolution. ``popsize`` is members-per-dimension (population NP = popsize*dim, min 5); ``F`` the
    differential weight, ``CR`` the crossover rate; ``seed`` makes the run deterministic (numpy
    default_rng). Stops at ``maxiter`` generations or when the population objective spread < ``tol``.
    Returns ``{x, fun, nfev, nit}``."""
    rng = np.random.default_rng(seed)
    bounds = np.asarray(bounds, dtype=float)
    lo, hi = bounds[:, 0], bounds[:, 1]
    d = len(bounds)
    NP = max(5, popsize * d)
    pop = lo + rng.random((NP, d)) * (hi - lo)
    fit = np.array([float(func(x)) for x in pop])
    nfev = NP
    nit = 0
    for nit in range(1, maxiter + 1):
        for i in range(NP):
            choices = [j for j in range(NP) if j != i]
            a, b, c = pop[rng.choice(choices, 3, replace=False)]
            mutant = np.clip(a + F * (b - c), lo, hi)
            cross = rng.random(d) < CR
            cross[rng.integers(d)] = True                 # force >=1 gene from the mutant
            trial = np.where(cross, mutant, pop[i])
            ft = float(func(trial)); nfev += 1
            if ft < fit[i]:
                pop[i], fit[i] = trial, ft
        if fit.max() - fit.min() < tol:
            break
    k = int(np.argmin(fit))
    return {"x": pop[k], "fun": float(fit[k]), "nfev": nfev, "nit": nit}
