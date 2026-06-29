r"""When does the per-station 2-D twist design break?  The fast-twist leaf coupling.

THE QUESTION (the next rung after the Frenet sweep)
---------------------------------------------------
twisting_quadrupole_pole.py and combined_function_frenet_sweep.py design a
twisting magnet as a STACK of independently-rotated 2-D cross-sections -- the
SLOW-TWIST (adiabatic) limit d phi / ds -> 0.  When the twist is FAST (the pole
orientation phi(s) = k s turns rapidly along the beam, k = d phi / ds), adjacent
"leaves" couple: the true 3-D field is NOT the 2-D stack.  This file quantifies
the coupling and finds the threshold -- the twist analogue of rung-1's
foliate-and-perturb (leaf_coupling_perturbation_3d.py), where the straight
magnet's coupling was ~ gap / L.

THE EXACT MODEL (no FEM -- the analytic helical multipole)
---------------------------------------------------------
A magnet whose order-n multipole twists with constant rate k has HELICAL symmetry
Phi(r, theta, s) = Phi(r, theta - k s).  The current-free (Laplace) harmonic with
that symmetry is the HELICAL MULTIPOLE (standard for helical undulators / Siberian
snakes / twisted quadrupoles):

    Phi_n = C * I_n(n k r) * sin( n (theta - k s) ),

I_n the modified Bessel function (Laplace + helical symmetry => the modified
Bessel equation).  The field (B = -grad Phi, cylindrical):

    B_r =  -C n k I_n'(n k r) sin(n psi)
    B_th =  -C (n/r) I_n(n k r) cos(n psi)        psi = theta - k s
    B_s =  +C n k I_n(n k r) cos(n psi)           <-- LONGITUDINAL, absent in 2-D

As k -> 0, I_n(n k r) -> (n k r / 2)^n / n!, B_s -> 0 and B_perp -> the pure 2-D
multipole rotated by phi(s) = k s -- EXACTLY the per-station 2-D stack (the
Frenet-sweep design).  The leaf coupling is the deviation of the full helical
field from that 2-D stack, controlled by the dimensionless TWIST-PER-APERTURE

    ka = k * a = 2 pi a / P          (a = aperture radius, P = 2 pi / k = pitch).

WHAT IS MEASURED (analytic, exact)
----------------------------------
On the aperture circle r = a:
  - eps(ka) = || B_perp(3-D) - B_perp(2-D stack) || / || B_perp(2-D stack) ||
              the transverse focusing-quality error -> scales as (ka)^2;
  - B_s / B_perp = the longitudinal field fraction -> scales as ka (first order);
  - the THRESHOLD ka* where eps = 1 percent -> the validity ratio P/a = 2 pi / ka*.

THE RESULT (the bridge to rung-1)
---------------------------------
eps ~ (ka)^2 (slope 2), B_s/B_perp ~ ka (slope 1), and eps = 1 percent at
ka* ~ 0.14, i.e. PITCH / APERTURE ~ 46.  The per-station 2-D twist design holds
when the pitch exceeds the aperture by ~ several tens -- the SAME "longitudinal
scale >> transverse scale by ~40x" rule as rung-1's straight magnet (L/gap ~ 40);
the twist version replaces gap/L by a/P.

run:  python twist_rate_leaf_coupling.py            # the coupling sweep + threshold
      python twist_rate_leaf_coupling.py --fig        # + figure
"""
import argparse
import json
import math
import os

import numpy as np


def _coupling(ka, n=2, n_theta=512):
    """Helical order-n multipole at the aperture r=a (set a=1, k=ka), s=0.
    Returns (eps_transverse, bs_over_bperp): the deviation of the full helical
    transverse field from the 2-D stack (per-station rotation), and the
    longitudinal-field fraction."""
    from scipy.special import iv, ivp
    th = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    arg = n * ka                                               # n k r at r=1
    lead = (n * ka / 2.0) ** n / math.factorial(n)            # I_n leading term
    # full helical field (C = 1)
    Br3 = -n * ka * ivp(n, arg) * np.sin(n * th)
    Bt3 = -n * iv(n, arg) * np.cos(n * th)
    Bs3 = +n * ka * iv(n, arg) * np.cos(n * th)
    # 2-D stack (the per-station rotation = the k->0 transverse shape)
    Br2 = -n * lead * np.sin(n * th)
    Bt2 = -n * lead * np.cos(n * th)
    dperp = np.sqrt((Br3 - Br2) ** 2 + (Bt3 - Bt2) ** 2)
    perp2 = np.sqrt(Br2 ** 2 + Bt2 ** 2)
    perp3 = np.sqrt(Br3 ** 2 + Bt3 ** 2)
    eps = float(np.sqrt(np.mean(dperp ** 2)) / np.sqrt(np.mean(perp2 ** 2)))
    bs = float(np.sqrt(np.mean(Bs3 ** 2)) / np.sqrt(np.mean(perp3 ** 2)))
    return eps, bs


def twist_sweep(kas=(0.05, 0.1, 0.2, 0.4, 0.8), n=2, tol=0.01):
    """Sweep the twist-per-aperture ka; fit the power laws and find the threshold
    ka* where the transverse error eps = tol."""
    kas = np.asarray(kas, dtype=float)
    eps = np.array([_coupling(k, n)[0] for k in kas])
    bs = np.array([_coupling(k, n)[1] for k in kas])
    eps_slope = float(np.polyfit(np.log(kas), np.log(eps), 1)[0])
    bs_slope = float(np.polyfit(np.log(kas), np.log(bs), 1)[0])
    kfine = np.linspace(0.02, float(kas[-1]), 600)
    efine = np.array([_coupling(k, n)[0] for k in kfine])
    ka_star = float(np.interp(tol, efine, kfine))
    return {
        "n": int(n), "tol": float(tol),
        "ka": kas.tolist(), "eps_transverse": eps.tolist(),
        "bs_over_bperp": bs.tolist(),
        "eps_slope": eps_slope, "bs_slope": bs_slope,
        "ka_star": ka_star, "pitch_over_aperture_star": float(2.0 * np.pi / ka_star),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2, help="multipole order (2 = quad)")
    ap.add_argument("--tol", type=float, default=0.01, help="eps threshold")
    ap.add_argument("--fig", action="store_true")
    args = ap.parse_args()

    print("=" * 76)
    print("Fast-twist leaf coupling -- when does the per-station 2-D twist break?")
    print("=" * 76)
    sw = twist_sweep(n=args.n, tol=args.tol)
    print(f"helical order n = {sw['n']}  (the exact analytic helical multipole "
          f"Phi = I_n(n k r) sin n(theta - k s))")
    print(f"  {'ka=2pi a/P':<12}{'eps (transv)':<15}{'B_s/B_perp':<12}")
    for ka, e, b in zip(sw["ka"], sw["eps_transverse"], sw["bs_over_bperp"]):
        print(f"  {ka:<12.3f}{e:<15.3e}{b:<12.3e}")
    print(f"  transverse error eps ~ ka^{sw['eps_slope']:.2f}  (the focusing-quality, "
          f"2nd order)")
    print(f"  longitudinal B_s/B_perp ~ ka^{sw['bs_slope']:.2f}  (the 3-D field, "
          f"1st order)")
    print(f"  => THRESHOLD: eps = {sw['tol']*100:.0f}% at ka* = {sw['ka_star']:.3f}, "
          f"i.e. PITCH/APERTURE = {sw['pitch_over_aperture_star']:.1f}")
    print("  => per-station 2-D twist design holds for pitch >> aperture (~tens) --")
    print("     the SAME longitudinal>>transverse rule as rung-1 (L/gap ~ 40),")
    print("     the twist version replacing gap/L by a/P (k -> 0 = the 2-D stack).")

    jpath = os.path.splitext(os.path.abspath(__file__))[0] + ".json"
    with open(jpath, "w") as f:
        json.dump(sw, f, indent=2)
    print(f"\nsaved {jpath}")

    if args.fig:
        _figure(sw)


def _figure(sw):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ka = np.array(sw["ka"])
    eps = np.array(sw["eps_transverse"])
    bs = np.array(sw["bs_over_bperp"])
    fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.2))
    # LEFT: eps and B_s/B_perp vs ka (log-log, slopes 2 and 1)
    ax[0].loglog(ka, eps, "o-", color="C3",
                 label=f"transverse error eps ~ ka^{sw['eps_slope']:.2f}")
    ax[0].loglog(ka, bs, "s-", color="C0",
                 label=f"B_s / B_perp ~ ka^{sw['bs_slope']:.2f}")
    ax[0].loglog(ka, eps[0] * (ka / ka[0]) ** 2, "k--", lw=0.8, label="slope 2")
    ax[0].loglog(ka, bs[0] * (ka / ka[0]), "k:", lw=0.8, label="slope 1")
    ax[0].axhline(sw["tol"], color="0.5", lw=0.8)
    ax[0].set_xlabel("twist per aperture  ka = 2*pi*a/P")
    ax[0].set_ylabel("relative coupling")
    ax[0].set_title("Fast-twist leaf coupling: transverse (2nd order)\n"
                    "+ longitudinal B_s (1st order)")
    ax[0].legend(fontsize=8)
    # RIGHT: the validity map -- eps vs pitch/aperture, with the threshold
    pa = 2.0 * np.pi / ka
    ax[1].loglog(pa, eps * 100, "o-", color="C3")
    ax[1].axhline(sw["tol"] * 100, color="0.5", lw=0.9, ls="--",
                  label=f"{sw['tol']*100:.0f}% tolerance")
    ax[1].axvline(sw["pitch_over_aperture_star"], color="C2", lw=1.2, ls=":",
                  label=f"P/a = {sw['pitch_over_aperture_star']:.0f} (threshold)")
    ax[1].set_xlabel("pitch / aperture  P/a")
    ax[1].set_ylabel("transverse error eps [%]")
    ax[1].set_title("Validity of the per-station 2-D twist design:\n"
                    "holds for P/a >> 1 (like rung-1's L/gap ~ 40)")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    png = os.path.splitext(os.path.abspath(__file__))[0] + ".png"
    fig.savefig(png, dpi=130)
    print(f"saved {png}")


if __name__ == "__main__":
    main()
