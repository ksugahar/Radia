"""maglev_sphere_force.py -- induced-dipole AC levitation force on a
conducting sphere, with the coefficient PINNED by the analytic perfect-
conductor limit and the frequency response REDUCED by CLN/Cauer.

This is the isotropic levitation-FORCE primitive of the 3D CLN-SIBC route
(docs/maglev/demos/cuboid/cln_sibc_cuboid_3d.py).  A sphere
is isotropic, so the scalar polarizability already ported applies directly
-- no anisotropic alpha tensor is needed to demonstrate (and verify) a
levitation force.  The cuboid a!=b!=c anisotropy is a separable refinement
-- and note it is SHAPE anisotropy (the dimensions a,b,c), NOT material
anisotropy: copper stays a scalar sigma/mu.  A field along z drives eddy
currents in the x-y cross-section (a x b), along x in the y-z cross-section
(b x c), etc., so a!=b!=c gives three different eddy time constants and
hence a direction-split alpha = diag(alpha_x, alpha_y, alpha_z).  It is the
AC/eddy-current generalization of the magnetostatic demagnetizing-factor
tensor (sphere: isotropic 1/3; ellipsoid/brick: direction-dependent, from
shape alone).

Physics (Landau-Lifshitz, Electrodynamics of Continuous Media, sec. 59):
  A conducting sphere radius a, conductivity sigma, in a uniform AC field
  develops an induced magnetic moment  m = (4 pi a^3 / mu0) G(x) B0  with
      x    = (1+i) a / delta,   delta = sqrt(2/(omega mu0 sigma)),
      G(x) = -1/2 [ 1 - 3/x^2 + (3/x) cot x ].
  Limits:  G(0)=0 (DC, field penetrates, no eddy response);
           G(inf)=-1/2 (perfect-conductor flux exclusion, diamagnet).

Time-averaged force on the induced dipole in a field gradient (B0 = peak
amplitude, real in space):
      <F> = (pi a^3 / mu0) Re[G(x)] grad(B0^2).
  Re[G] < 0  =>  force toward weaker field  =>  LIFT away from the coil.
  Coefficient cross-checked against the perfect-conductor sphere energy
  U = -1/2 m.B  ->  <F> -> -(pi a^3 / 2 mu0) grad(B0^2)  at high freq.

CLN/Cauer reduction: the sphere has the Foster form
      G(s) = -1/2 + sum_n (3/(n pi)^2) / (1 - s tau_n),  tau_n = mu0 sigma a^2/(n pi)^2,
  i.e. a ladder of relaxation modes -- the same modal structure the cuboid
  CLN-SIBC port reduces.  A Lanczos/Cauer reduction of this modal system
  converges to G with a handful of stages.

Run:  python maglev_sphere_force.py     (a few seconds; pure numpy)
"""
import json
import os
import sys

import numpy as np

mu0 = 4 * np.pi * 1e-7
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from _validation_output import validation_output  # noqa: E402

# ----- sphere -----
a = 5e-3          # radius 5 mm
sigma = 5.8e7     # copper
NMODES = 1000     # modal truncation for the Foster reference


def delta(omega):
    return np.sqrt(2.0 / (omega * mu0 * sigma))


def _cot(x):
    """cot(x) = i (e^{2ix}+1)/(e^{2ix}-1); stable for Im(x)>0 (no tan overflow)."""
    e = np.exp(2j * x)
    return 1j * (e + 1.0) / (e - 1.0)


def G_exact(omega):
    """Closed-form sphere magnetic response G(x), x=(1+i) a/delta."""
    x = (1 + 1j) * a / delta(omega)
    return -0.5 * (1.0 - 3.0 / x**2 + 3.0 / x * _cot(x))


# ----- Foster (modal) form -----
_n = np.arange(1, NMODES + 1)
TAU = mu0 * sigma * a**2 / (_n * np.pi) ** 2          # relaxation times
W = 3.0 / (_n * np.pi) ** 2                            # modal residues (sum W = 1/2)


def G_modal(omega, N=NMODES):
    s = 1j * omega
    n = slice(0, N)
    return -0.5 + np.sum(W[n] / (1.0 - s * TAU[n]))


def G_cln(omega, nstage):
    """Lanczos/Cauer reduction of the modal system b^T (I - sT)^{-1} b.

    T = diag(tau_n), b_n = sqrt(W_n); G+1/2 = b^T (I - sT)^{-1} b exactly.
    Lanczos on (T, b) builds a tridiagonal T_m; the reduced transfer
    function is ||b||^2 e1^T (I - s T_m)^{-1} e1 -- the Cauer continued
    fraction truncated at nstage relaxation stages.
    """
    b = np.sqrt(W)
    beta = np.linalg.norm(b)
    q = b / beta
    Q, alphas, betas = [q], [], []
    qprev = np.zeros_like(q)
    bprev = 0.0
    for _ in range(nstage):
        z = TAU * q                      # T is diagonal
        al = q @ z
        alphas.append(al)
        z = z - al * q - bprev * qprev
        for u in Q:                      # full reorthogonalization (small m)
            z = z - (u @ z) * u
        bn = np.linalg.norm(z)
        if bn < 1e-14 or len(alphas) == nstage:
            break
        betas.append(bn)
        qprev, bprev = q, bn
        q = z / bn
        Q.append(q)
    Tm = np.diag(alphas) + np.diag(betas, 1) + np.diag(betas, -1)
    s = 1j * omega
    m = Tm.shape[0]
    e1 = np.zeros(m); e1[0] = 1.0
    y = np.linalg.solve(np.eye(m) - s * Tm, e1)
    return -0.5 + beta**2 * (e1 @ y)


def levitation_force(omega, B0, gradB2):
    """Time-averaged vertical lift force [N].  gradB2 = d(B0^2)/dz < 0 going
    away from the coil; returns the magnitude of the upward lift."""
    ReG = G_exact(omega).real
    Fz = (np.pi * a**3 / mu0) * ReG * gradB2      # (pi a^3/mu0) ReG grad(B0^2)
    return Fz


def main():
    print(f"Cu sphere a={a*1e3:.1f} mm, sigma={sigma:.2e};  modal modes={NMODES}")

    # ----- check 1: closed form vs Foster modal sum -----
    print(f"\n[check 1] closed-form G(x)  vs  Foster modal sum ({NMODES} modes)")
    print("   f[Hz]   a/delta   Re G_exact  Re G_modal   |diff|")
    for f in (50, 500, 5e3, 5e4, 5e5):
        w = 2 * np.pi * f
        ge, gm = G_exact(w), G_modal(w)
        print(f"  {f:7.0f}  {a/delta(w):7.3f}   {ge.real:+9.5f}  {gm.real:+9.5f}"
              f"   {abs(ge - gm):.2e}")

    # ----- check 2: physical limits of G -----
    print("\n[check 2] limits of G(x)")
    g_lo = G_exact(2 * np.pi * 1e-2)
    g_hi = G_exact(2 * np.pi * 1e9)
    print(f"   low  freq (f=0.01 Hz): G = {g_lo.real:+.3e} {g_lo.imag:+.3e}j"
          f"   (Re->0, no eddy response)")
    print(f"   high freq (f=1 GHz)  : G = {g_hi.real:+.6f} {g_hi.imag:+.3e}j"
          f"   (Re->-0.5 perfect-conductor exclusion)")
    assert abs(g_hi.real + 0.5) < 1e-3, "high-freq G must -> -1/2"
    assert abs(g_lo.real) < 1e-3 and g_lo.real <= 0, "low-freq Re G -> 0^-"

    # Re G < 0 across the whole band -> force is a LIFT at every frequency
    fs = np.logspace(0, 8, 200)
    ReGs = np.array([G_exact(2 * np.pi * f).real for f in fs])
    assert np.all(ReGs <= 1e-9), "Re G must be <= 0 (always a lift)"
    print(f"   Re G(f) < 0 for all f in [1 Hz, 100 MHz]  -> lift at every freq  OK")

    # ----- check 3: CLN/Cauer reduction converges -----
    # the CLN reduces the N-mode modal operator; its target is G_modal(N),
    # not the closed form (the Foster<->closed-form gap is check 1).
    print("\n[check 3] CLN/Cauer reduction of the sphere modal system (f=5 kHz)")
    w = 2 * np.pi * 5e3
    gtgt = G_modal(w)
    print("   stages   Re G_CLN     rel.err vs full modal")
    for m in (1, 2, 3, 4, 6, 8):
        gc = G_cln(w, m)
        print(f"     {m:2d}     {gc.real:+9.6f}    {abs(gc - gtgt)/abs(gtgt)*100:8.4f} %")
    assert abs(G_cln(w, 8) - gtgt) / abs(gtgt) < 1e-3, "CLN must reduce the modal system"

    # ----- check 4: perfect-conductor force coefficient -----
    # gradient field: B0=0.1 T, scale length L -> grad(B0^2) ~ -B0^2/L (away from coil)
    B0, L = 0.1, 0.05
    gradB2 = -B0**2 / L
    print(f"\n[check 4] levitation force, B0={B0} T, grad(B0^2)={gradB2:+.3f} T^2/m")
    F_pc = (np.pi * a**3 / (2 * mu0)) * abs(gradB2)   # analytic high-freq limit
    print(f"   perfect-conductor analytic lift = (pi a^3/2 mu0)|gradB2| = {F_pc*1e3:.4f} mN")
    print("   f[Hz]    a/delta    Re G     lift[mN]   /F_pc")
    sweep = {}
    for f in (50, 500, 5e3, 5e4, 5e5, 5e6, 5e7, 5e8):
        w = 2 * np.pi * f
        Fz = levitation_force(w, B0, gradB2)
        lift = abs(Fz)
        print(f"  {f:9.0f}  {a/delta(w):8.3f}  {G_exact(w).real:+7.4f}"
              f"  {lift*1e3:8.4f}  {lift/F_pc:6.4f}")
        sweep[f] = lift
    f_hi = 1e9
    assert abs(abs(levitation_force(2*np.pi*f_hi, B0, gradB2)) / F_pc - 1.0) < 1e-3, \
        "high-freq lift must approach the perfect-conductor coefficient"
    print("   high-freq lift -> perfect-conductor coefficient (ratio -> 1)  OK")

    # ----- persist -----
    fs_plot = np.logspace(1, 7.5, 120)
    lift_plot = [abs(levitation_force(2 * np.pi * f, B0, gradB2)) for f in fs_plot]
    data = {
        "sphere": {"a_m": a, "sigma": sigma},
        "field": {"B0_T": B0, "scale_L_m": L, "gradB2_T2_per_m": gradB2},
        "perfect_conductor_lift_N": F_pc,
        "force_formula": "F = (pi a^3/mu0) Re[G(x)] grad(B0^2)",
        "sweep_lift_N": {str(int(f)): v for f, v in sweep.items()},
        "freq_Hz": fs_plot.tolist(),
        "lift_N": lift_plot,
    }
    output = validation_output("maglev_sphere_force_results.json", HERE)
    with open(output, "w") as fp:
        json.dump(data, fp, indent=2, allow_nan=False)
    print(f"\n wrote {output}")
    plot(fs_plot, np.array(lift_plot), F_pc)


def plot(fs, lift, F_pc):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "font.family": "serif", "axes.linewidth": 0.8})
    fig, ax = plt.subplots(figsize=(8 / 2.54, 6 / 2.54))
    ax.semilogx(fs, lift * 1e3, "C0-", lw=1.2, label="induced-dipole + CLN")
    ax.axhline(F_pc * 1e3, color="0.6", ls=":", lw=0.8)
    ax.text(fs[0], F_pc * 1e3, " perfect conductor", va="bottom", ha="left",
            fontsize=8, color="0.4")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("levitation force (mN)")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout(pad=0.3)
    out = os.path.join(HERE, "maglev_sphere_force_vs_freq.png")
    fig.savefig(out, dpi=300)
    print(f" wrote {out}")


if __name__ == "__main__":
    main()
