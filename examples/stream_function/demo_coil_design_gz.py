"""Cylindrical z-gradient (Gz) coil design via the stream function method.

End-to-end demonstration of the stream function method built on the
(ACA+)+TSVD least-norm solver:

    target field  ->  surface stream function psi(z)  ->  wire rings  ->  verify

A Gz gradient (target Bz = Gz*z on the axis) is axisymmetric, so the surface
current on the cylinder (radius a) is purely azimuthal and the stream function
reduces to psi(z).  The natural basis is therefore a set of full azimuthal
current RINGS at z_1..z_N (radius a).  The coupling A(i,j) = Bz at field point i
from unit ring j is Radia's Biot-Savart (ObjFlmCur ring + radia.Fld via
radia_field_kernel).  We solve the underdetermined system  A I = Gz*z  with the
TSVD-regularised pseudo-inverse (ACA+ accelerated) for the continuous ring
current density I(z); the stream function psi(z) = cumulative integral of I.

Contouring psi(z) at equal current spacing gives discrete equal-current wire
rings (a generalised Maxwell pair); we build them as ObjFlmCur rings and verify
the on-axis gradient linearity.

(Using full rings -- the correct reduced model for the axisymmetric Gz problem
-- also keeps every field evaluation in Radia's reliable full-ring path.)

Run:  python demo_coil_design_gz.py [--nrings N] [--k MODES] [--nwires N]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import radia as rad
from radia.stream_function import aca_tsvd, pseudo_inverse_solve, radia_field_kernel

NSEG = 120  # polygon segments per ring (>=72 keeps polygon error < 1e-3)


def ring(a, z0, current):
    t = np.linspace(0.0, 2.0 * np.pi, NSEG + 1)
    pts = [[a * np.cos(tt), a * np.sin(tt), z0] for tt in t]
    return rad.ObjFlmCur(pts, current)


def contour(psi, z, n_wires):
    """Equal-current contour of psi(z) -> [(z_wire, current)], current=+/-dI."""
    lo, hi = float(psi.min()), float(psi.max())
    dI = (hi - lo) / n_wires
    if dI <= 0:
        return []
    out = []
    for lev in lo + (np.arange(n_wires) + 0.5) * dI:
        for j in range(len(z) - 1):
            a0, a1 = psi[j] - lev, psi[j + 1] - lev
            if a0 == 0.0:
                a0 = 1e-300
            if a0 * a1 < 0.0:
                frac = a0 / (a0 - a1)
                zc = z[j] + frac * (z[j + 1] - z[j])
                slope = psi[j + 1] - psi[j]
                out.append((zc, dI * (1.0 if slope > 0 else -1.0)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nrings", type=int, default=80)
    ap.add_argument("--k", type=int, default=0, help="TSVD modes (0 = auto)")
    ap.add_argument("--nwires", type=int, default=16)
    args = ap.parse_args()

    a = 0.15           # cylinder radius [m]
    L = 0.50           # cylinder length [m]
    Gz = 1.0           # target gradient dBz/dz [T/m]
    dsv = 0.10         # DSV half-length [m] (target region |z| <= dsv)

    rad.UtiDelAll()
    try:
        # --- ring basis + on-axis target -------------------------------------
        z_ring = np.linspace(-L / 2, L / 2, args.nrings)
        basis = [ring(a, zr, 1.0) for zr in z_ring]
        N = len(basis)

        zt = np.linspace(-dsv, dsv, 25)
        obs = np.column_stack([np.zeros(zt.size), np.zeros(zt.size), zt])
        M = obs.shape[0]
        B_target = Gz * zt

        entry = radia_field_kernel(obs, basis, component=2, field="b")
        modes = min(M, N) if args.k <= 0 else min(args.k, M, N)
        res = aca_tsvd(M, N, entry, modes=modes, kmax=min(M, N),
                       aca_eps=1.0e-8, method=3)
        k_use = res.modes if args.k <= 0 else min(args.k, res.modes)
        I = pseudo_inverse_solve(res, B_target, k_mode=k_use)   # ring currents

        # continuous-coil residual (sanity: A I should match Gz*z)
        rad.UtiDelAll()
        cont = rad.ObjCnt([ring(a, zr, ii) for zr, ii in zip(z_ring, I)])
        Bc = np.asarray(rad.Fld(cont, "b", obs.tolist()),
                        dtype=float).reshape(-1, 3)[:, 2]
        cont_res = float(np.linalg.norm(Bc - B_target) / np.linalg.norm(B_target))
        rad.UtiDelAll()

        print("=== Cylindrical Gz gradient coil (stream function + (ACA+)+TSVD) ===")
        print(f"cylinder a={a*1e3:.0f} mm, L={L*1e3:.0f} mm | {N} candidate rings"
              f" | on-axis target M={M}, Gz={Gz} T/m")
        print(f"ACA+ k_aca = {res.k_aca} (of min(M,N)={min(M, N)}),  "
              f"TSVD modes used = {k_use}")
        print(f"continuous-current-density coil: on-axis ||A I - Gz z||/||Gz z|| "
              f"= {cont_res:.3e}")

        # --- contour psi(z) = cumulative current -> equal-current wire rings --
        psi = np.concatenate([[0.0], np.cumsum(0.5 * (I[1:] + I[:-1]))])  # cumtrapz
        wires = contour(psi, z_ring, args.nwires)
        print(f"discretised into {len(wires)} equal-current wire rings "
              f"(I/ring = {abs(wires[0][1]) if wires else 0:.4g} A-turns)")

        # --- verify the DISCRETE wire coil on the axis -----------------------
        rad.UtiDelAll()
        coil = rad.ObjCnt([ring(a, zc, cur) for zc, cur in wires])
        zc = np.linspace(-dsv, dsv, 41)
        axis = np.column_stack([np.zeros_like(zc), np.zeros_like(zc), zc])
        Bz = np.asarray(rad.Fld(coil, "b", axis.tolist()),
                        dtype=float).reshape(-1, 3)[:, 2]
        G_fit = float(np.dot(zc, Bz) / np.dot(zc, zc))
        nonlin = float(np.max(np.abs(Bz - G_fit * zc)) / (np.max(np.abs(Bz)) + 1e-30))
        print(f"discrete wire coil: fitted dBz/dz = {G_fit:.4g} (per unit psi)")
        print(f"on-axis gradient nonlinearity over DSV = {nonlin:.3e}")
        rad.UtiDelAll()

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
            ax1.plot(z_ring * 1e3, I / np.abs(I).max(), "o-", ms=3,
                     label="ring current I(z) (norm.)")
            ax1.plot(z_ring * 1e3, psi / np.abs(psi).max(), "s-", ms=3,
                     label="stream fn psi(z) (norm.)")
            for zc0, cur in wires:
                ax1.axvline(zc0 * 1e3, color="r" if cur > 0 else "b", alpha=0.3, lw=1)
            ax1.set_xlabel("z [mm]"); ax1.set_title("Current density + wire rings")
            ax1.legend(); ax1.grid(alpha=0.3)
            ax2.plot(zc * 1e3, Bz, "o-", ms=3, label="discrete coil Bz")
            ax2.plot(zc * 1e3, G_fit * zc, "k--", label="linear fit")
            ax2.set_xlabel("z [mm]"); ax2.set_ylabel("Bz [T]")
            ax2.set_title(f"On-axis gradient (nonlin {nonlin:.1e})")
            ax2.legend(); ax2.grid(alpha=0.3)
            out = HERE / "demo_coil_design_gz.png"
            fig.tight_layout(); fig.savefig(out, dpi=110)
            print(f"Saved plot: {out.name}")
        except ImportError:
            print("(matplotlib not installed; skipped plot)")
    finally:
        rad.UtiDelAll()


if __name__ == "__main__":
    main()
