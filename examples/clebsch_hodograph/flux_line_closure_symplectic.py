r"""Flux-line closure: the dynamical face of the Clebsch / de Rham structure.

A magnetic flux line is an integral curve of B:  dx/ds = B(x).  In 2-D,
`B = grad(A_z) x z_hat = (dA_z/dy, -dA_z/dx)`, so the flux-line ODE IS Hamilton's
equations with **A_z as the Hamiltonian** -- i.e. A_z (= the Clebsch potential,
`hdiv_vim_clebsch_2d_az.py`) is conserved along a flux line, so flux lines are
A_z level sets and **close**.

That closure has TWO requirements, and this script isolates both:

  (1) the FIELD must be a CLOSED 2-form  (`div B = 0`, the Clebsch/loop part) --
      otherwise no global A_z exists and the flux line cannot close.  Adding a
      charge/star part `eps*grad(phi)` (the de Rham complement) destroys the
      conserved A_z and the line SPIRALS.  [the de Rham / edge-FE requirement --
      flux lines computed from an edge (H(curl)) potential `B = curl A` are
      exactly divergence-free; cf. Noguchi, "Flux-line computation from
      hexahedral edge-finite-element results / bubble placement", IEEJ (JP).]

  (2) the INTEGRATOR must be SYMPLECTIC (conserve the Hamiltonian A_z) --
      a non-symplectic step (RK4) drifts A_z and the line spirals even for a
      perfectly closed field.  [the accelerator-tracking requirement -- Sugahara
      2020, "Implicit symplectic flux-line tracking", which noted the analogy to
      circular-accelerator beam-orbit tracking.]

A_z = the Clebsch potential = the flux-line-flow Hamiltonian is the single object
behind both faces -- the field (FEEC / Clebsch) face and the dynamical
(symplectic) face.

run:  python flux_line_closure_symplectic.py
"""
import os
import numpy as np


# Analytic 2-D flux function A_z (the Clebsch potential = the Hamiltonian).  Off-centre + anisotropic
# so the closed flux line is a non-trivial (elliptic, shifted) loop, not a bare circle.
def A_z(x, y):
    return -0.5 * (x * x + 1.6 * y * y) - 0.2 * x


def B_loop(x, y):
    """B = grad(A_z) x z_hat = (dA_z/dy, -dA_z/dx) -- a CLOSED 2-form (div B = 0)."""
    return np.array([-1.6 * y, x + 0.2])


def grad_phi(x, y):
    """the charge / star complement grad(phi), phi = (x^2+y^2)/2 -> div = 2 (NOT a closed 2-form)."""
    return np.array([x, y])


def _trace_forward_euler(vel, x0, ds, n):
    """non-symplectic FORWARD (explicit) Euler -- SAME 1st order as the symplectic Euler, so the
    contrast is the symplecticity itself (bounded vs secularly growing), not the order."""
    xs = np.zeros((n + 1, 2)); xs[0] = x0
    for i in range(n):
        xs[i + 1] = xs[i] + ds * vel(*xs[i])
    return xs


def _trace_symplectic(x0, ds, n):
    """symplectic (semi-implicit Euler) trace of the Hamiltonian flow x'=dA/dy, y'=-dA/dx.
    Uses the updated x in the y-step -> preserves a nearby invariant (A_z bounded, no secular drift)."""
    xs = np.zeros((n + 1, 2)); xs[0] = x0
    for i in range(n):
        x, y = xs[i]
        xn = x + ds * (-1.6 * y)            # x' = dA_z/dy = -1.6 y
        yn = y + ds * (xn + 0.2)            # y' = -dA_z/dx = x + 0.2  (uses xn -> symplectic)
        xs[i + 1] = [xn, yn]
    return xs


def _trace_symplectic_charged(x0, ds, n, eps):
    """The SAME semi-implicit (symplectic-style) step, but on the CHARGED field B_loop + eps*grad(phi)
    (div != 0, NOT a closed 2-form).  Isolates the FIELD requirement: same integrator, but no global
    A_z exists, so the line spirals regardless of the integrator."""
    xs = np.zeros((n + 1, 2)); xs[0] = x0
    for i in range(n):
        x, y = xs[i]
        xn = x + ds * (-1.6 * y + eps * x)
        yn = y + ds * ((xn + 0.2) + eps * y)
        xs[i + 1] = [xn, yn]
    return xs


def _metrics(xs, x0):
    Az = np.array([A_z(p[0], p[1]) for p in xs])
    drift = float(np.max(np.abs(Az - Az[0])) / (abs(Az[0]) + 1e-30))
    # closure: smallest return distance to the start after leaving a neighbourhood
    d = np.linalg.norm(xs - np.asarray(x0), axis=1)
    left = np.where(d > 0.5 * d.max())[0]
    ret = float(d[left[0]:].min()) if len(left) else float(d.min())
    return drift, ret, Az


def analyze(turns=25, steps_per_turn=300):
    x0 = np.array([2.0, 0.0])
    period = 2 * np.pi / np.sqrt(1.6)          # linear oscillator period (x'=-1.6y, y'=x)
    ds = period / steps_per_turn
    n = int(turns * steps_per_turn)
    sym = _trace_symplectic(x0, ds, n)
    fe = _trace_forward_euler(lambda x, y: B_loop(x, y), x0, ds, n)    # same order, non-symplectic
    eps = 0.03
    mix = _trace_symplectic_charged(x0, ds, n, eps)                    # symplectic step, but charged field
    drift_sym, ret_sym, Az_sym = _metrics(sym, x0)
    drift_fe, ret_fe, Az_fe = _metrics(fe, x0)
    drift_mix, ret_mix, Az_mix = _metrics(mix, x0)
    return {
        "turns": turns, "n": n,
        "drift_sym": drift_sym, "ret_sym": ret_sym,
        "drift_fe": drift_fe, "ret_fe": ret_fe,
        "drift_mix": drift_mix, "ret_mix": ret_mix, "eps": eps,
        "sym": sym, "fe": fe, "mix": mix, "Az_sym": Az_sym, "Az_fe": Az_fe, "Az_mix": Az_mix,
    }


def _plot(r):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.4), dpi=150)
    # left: the INTEGRATOR contrast on the SAME closed field, matched scale
    ax.plot(r["fe"][:, 0], r["fe"][:, 1], "C2-", lw=0.7, label="forward Euler (spirals out)")
    ax.plot(r["sym"][:, 0], r["sym"][:, 1], "C0-", lw=1.3, label="symplectic (closes)")
    ax.plot(2.0, 0.0, "k.", ms=10)
    lim = 1.15 * np.max(np.abs(r["fe"]))
    ax.set_aspect("equal"); ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_title("same 1st order, closed field:\nsymplectic closes vs forward Euler spirals")
    ax.legend(fontsize=8, loc="lower left")
    # right: A_z drift (log) -- both requirements at once
    s = np.arange(len(r["Az_sym"]))
    ax2.semilogy(s, np.abs(r["Az_sym"] / r["Az_sym"][0] - 1) + 1e-16, "C0-", lw=1, label="symplectic, closed (CLOSES)")
    ax2.semilogy(s, np.abs(r["Az_fe"] / r["Az_fe"][0] - 1) + 1e-16, "C2-", lw=1, label="forward Euler (non-sympl)")
    ax2.semilogy(s, np.abs(r["Az_mix"] / r["Az_mix"][0] - 1) + 1e-16, "C3-", lw=1, label="+ charge (not closed)")
    ax2.set_xlabel("step"); ax2.set_ylabel("$|A_z/A_z(0) - 1|$  (Hamiltonian drift)")
    ax2.set_title("A_z drift = loss of closure")
    ax2.legend(fontsize=8, loc="lower right")
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.tight_layout(); fig.savefig(png, bbox_inches="tight"); plt.close(fig)
    print(f"  figure saved: {png}")


def main():
    print("Flux-line closure: the dynamical face of the Clebsch / de Rham structure\n")
    r = analyze()
    print(f"  {r['turns']} turns, {r['n']} steps.  A_z = Clebsch potential = flux-line Hamiltonian.")
    print(f"  (1) CLOSED 2-form + SYMPLECTIC   : A_z drift={r['drift_sym']:.2e}  return dist={r['ret_sym']:.2e}  -> CLOSES")
    print(f"  (2) closed field, FORWARD Euler  : A_z drift={r['drift_fe']:.2e}  return dist={r['ret_fe']:.2e}  -> spirals out")
    print(f"  (3) + charge (not a closed 2-form): A_z drift={r['drift_mix']:.2e}  return dist={r['ret_mix']:.2e}  -> SPIRALS")
    print(f"\n  ratio drift(forward Euler)/drift(symplectic) = {r['drift_fe']/r['drift_sym']:.1f}x  (same order -> the symplecticity)")
    print(f"  ratio drift(charge)/drift(symplectic)        = {r['drift_mix']/r['drift_sym']:.1f}x  (the field is not a closed 2-form)")
    print("\n  => flux-line closure needs BOTH a closed 2-form field (de Rham / edge-FE -- Noguchi) AND")
    print("     a symplectic A_z-conserving integrator (accelerator tracking -- Sugahara 2020).  A_z =")
    print("     the Clebsch potential = the flux-line-flow Hamiltonian is the single object behind both.")
    _plot(r)


if __name__ == "__main__":
    main()
