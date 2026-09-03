"""ih_optimize.py

LLM-driven optimization loop harness for IH coil design.

IHOptimizer takes a pluggable "proposer" (a callable that maps
history -> next CoilSpec) and drives it through a fixed number of
iterations against a cached IHWorkpieceContext. Each trial is
recorded with its spec, workpiece metrics, objective value and
timing. The best-so-far trial is returned.

Radia ships only ``RandomProposer`` as a deterministic baseline. AI-guided
proposal belongs to the MCP/agent layer: pass any callable with the contract
below instead of launching a vendor-specific CLI from the solver package.

Custom proposers are trivial: any callable with signature
    proposer(history: list[Trial]) -> CoilSpec
works.
"""

import dataclasses
import time
from typing import Callable, Optional

import numpy as np

from radia.coil_spec import CoilSpec, build_coil_from_spec


@dataclasses.dataclass
class Trial:
    iteration: int
    spec: CoilSpec
    metrics: dict           # empty when evaluation failed
    objective: float        # -inf / +inf when eval failed
    t_eval_ms: float
    failure: Optional[str] = None


class IHOptimizer:
    """Sequential optimizer: propose -> evaluate -> record -> repeat.

    Args:
        ctx: IHWorkpieceContext (BEM solver cached once up-front).
        proposer: callable(history) -> CoilSpec. First call gets [].
        objective_fn: callable(metrics: dict) -> float. Default:
            lambda m: m['P_total'] -- maximize heating power.
        maximize: True to maximize objective_fn, False to minimize.
        nw, nh: filament grid density passed to ctx.evaluate (default 5, 8).
        **kwargs: forwarded to ctx.evaluate.
    """

    def __init__(self, ctx, proposer: Callable,
                 objective_fn: Optional[Callable[[dict], float]] = None,
                 maximize: bool = True, nw: int = 5, nh: int = 8,
                 **kwargs):
        self.ctx = ctx
        self.proposer = proposer
        if objective_fn is None:
            objective_fn = lambda m: m['P_total']
        self.objective_fn = objective_fn
        self.maximize = maximize
        self.nw = nw
        self.nh = nh
        self.kwargs = kwargs
        self.history: list[Trial] = []

    def _bad_objective(self) -> float:
        return float('-inf') if self.maximize else float('inf')

    def _is_better(self, a: float, b: float) -> bool:
        return (a > b) if self.maximize else (a < b)

    def run(self, n_iter: int, verbose: bool = True) -> Optional[Trial]:
        for i in range(n_iter):
            spec = self.proposer(self.history)
            t0 = time.perf_counter()
            failure = None
            metrics: dict = {}
            obj = self._bad_objective()
            try:
                coil = build_coil_from_spec(spec)
                metrics = self.ctx.evaluate(
                    coil, nw=self.nw, nh=self.nh, **self.kwargs)
                obj = float(self.objective_fn(metrics))
            except Exception as e:
                failure = f"{type(e).__name__}: {e}"
            t_eval = (time.perf_counter() - t0) * 1e3

            trial = Trial(iteration=i, spec=spec, metrics=metrics,
                          objective=obj, t_eval_ms=t_eval, failure=failure)
            self.history.append(trial)

            if verbose:
                self._print_trial(trial)

        return self.best()

    def _print_trial(self, trial: Trial):
        best = self.best()
        is_best = best is not None and best.iteration == trial.iteration
        flag = "  *BEST" if is_best else ""
        s = trial.spec
        params = (f"r={s.r_wire*1e3:.2f}mm" if s.profile_type == 'circle'
                  else f"{s.wire_w*1e3:.2f}x{s.wire_h*1e3:.2f}mm")
        head = (f"  [{trial.iteration:2d}] {params}  R_maj={s.R_major*1e3:.1f}mm  "
                f"arc={s.arc_deg:.0f}deg  ")
        if trial.failure:
            print(head + f"FAIL: {trial.failure[:60]}")
            return
        m = trial.metrics
        print(head + f"H={m['H_t_rms']:.3f}  P={m['P_total']:.3e}  "
              f"obj={trial.objective:.4e}  t={trial.t_eval_ms:.0f}ms{flag}")

    def best(self) -> Optional[Trial]:
        valid = [t for t in self.history if t.failure is None]
        if not valid:
            return None
        key = (lambda t: t.objective) if self.maximize else \
              (lambda t: -t.objective)
        return max(valid, key=key)

    def summary(self) -> str:
        if not self.history:
            return "no trials"
        ok = [t for t in self.history if t.failure is None]
        fail = [t for t in self.history if t.failure is not None]
        best = self.best()
        lines = [
            f"trials: {len(self.history)} total, "
            f"{len(ok)} ok, {len(fail)} failed",
            f"best  : iter {best.iteration}  obj={best.objective:.4e}",
            f"spec  : {best.spec.to_json()}",
            f"H_t_rms={best.metrics['H_t_rms']:.4f}  "
            f"P_total={best.metrics['P_total']:.4e}  "
            f"Q_total={best.metrics['Q_total']:.4e}",
        ]
        return "\n".join(lines)


class RandomProposer:
    """Uniform random sampling baseline for round-wire torus designs.

    ranges is a dict of named parameters; expected keys for circle
    torus: 'r_wire_mm', 'R_major_mm', 'arc_deg'. Values are (lo, hi).
    """

    def __init__(self, ranges: dict, profile_type: str = 'circle',
                 current: float = 1.0, seed: Optional[int] = 0):
        self.ranges = ranges
        self.profile_type = profile_type
        self.current = current
        self.rng = np.random.default_rng(seed)

    def __call__(self, history: list) -> CoilSpec:
        r = self.ranges
        R_major = self.rng.uniform(*r['R_major_mm']) * 1e-3
        arc_deg = float(self.rng.uniform(*r['arc_deg']))
        if self.profile_type == 'circle':
            r_wire = self.rng.uniform(*r['r_wire_mm']) * 1e-3
            # Clamp to keep inner radius positive (avoid immediate fail)
            r_wire = min(r_wire, 0.9 * R_major)
            return CoilSpec(profile_type='circle', r_wire=r_wire,
                            R_major=R_major, arc_deg=arc_deg,
                            current=self.current)
        wire_w = self.rng.uniform(*r['wire_w_mm']) * 1e-3
        wire_h = self.rng.uniform(*r['wire_h_mm']) * 1e-3
        return CoilSpec(profile_type='rect', wire_w=wire_w, wire_h=wire_h,
                        R_major=R_major, arc_deg=arc_deg,
                        current=self.current)
