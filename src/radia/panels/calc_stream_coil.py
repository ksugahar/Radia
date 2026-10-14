"""Stream-function cylindrical Gz gradient-coil designer (Stage-2 CLI).

Layer 4 computation script (no GUI, no Cubit).  Designs a cylindrical
z-gradient (Gz) gradient coil via the stream function method built on the
(ACA+)+TSVD least-norm solver:

    target field  ->  surface stream function psi(z)  ->  wire rings  ->  verify

A Gz gradient (target Bz = Gz*z) is axisymmetric, so the surface current on
the cylinder (radius a) is purely azimuthal and the stream function reduces to
psi(z).  The natural basis is a set of full azimuthal current RINGS at
z_1..z_N (radius a).  The coupling A(i,j) = Bz at field point i from unit ring
j is Radia's Biot-Savart (rad.ObjFlmCur ring + radia.Fld via
radia_field_kernel).  We solve the underdetermined system  A I = Gz*z  with the
TSVD-regularised pseudo-inverse (ACA+ accelerated) for the continuous ring
current density I(z); the stream function psi(z) = cumulative integral of I.

Contouring psi(z) at equal current spacing gives discrete equal-current wire
rings (a generalised Maxwell pair); we build them as rad.ObjFlmCur rings and
verify the on-axis gradient linearity.

Off-axis volume target: an on-axis-only target leaves the off-axis (volume)
gradient uniformity unconstrained.  Because the basis is FULL CIRCULAR RINGS,
the produced field is ALWAYS axisymmetric, so adding target points at a SINGLE
azimuth (rho, 0, z) for several rho > 0 is sufficient AND safe -- it enforces
Bz = Gz*z over the whole DSV cylinder without breaking symmetry.

Using full rings (the correct reduced model for the axisymmetric Gz problem)
also keeps every field evaluation in Radia's reliable full-ring path; the known
rad.ObjFlmCur orientation bug only affects small TILTED loops, not full rings.

Called as a subprocess (Layer 3 panel / test) -- prints a JSON object as the
last stdout line.  Example::

    python calc_stream_coil.py --radius 0.15 --length 0.5 --gradient 1.0

JSON keys: k_aca, modes, n_wires, fitted_dBdz, gradient_nonlinearity,
continuous_residual, wire_z, wire_I.
"""

import argparse
import os
import sys

import numpy as np

NSEG = 120  # polygon segments per ring (>=72 keeps polygon error < 1e-3)


def _ring(rad, a, z0, current):
    """Build a full azimuthal current ring (rad.ObjFlmCur) at z=z0, radius a."""
    t = np.linspace(0.0, 2.0 * np.pi, NSEG + 1)
    pts = [[a * np.cos(tt), a * np.sin(tt), z0] for tt in t]
    return rad.ObjFlmCur(pts, current)


def _contour(psi, z, n_wires):
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


def design_stream_coil(radius, length, gradient, dsv, n_rings, n_wires,
                       modes, aca_eps, method):
    """Design a cylindrical Gz gradient coil via stream function + (ACA+)+TSVD.

    Args:
        radius: Cylinder radius a [m].
        length: Cylinder length L [m] (rings span z in [-L/2, +L/2]).
        gradient: Target gradient dBz/dz Gz [T/m].
        dsv: DSV half-length [m] (target region |z| <= dsv).
        n_rings: Number of candidate rings (basis size N).
        n_wires: Number of equal-current wire rings to contour into.
        modes: TSVD modes (0 = auto = min(M, N), clamped to ACA+ rank).
        aca_eps: ACA+ stopping tolerance.
        method: TSVD recompression method (2 or 3).

    Returns:
        dict with keys k_aca, modes, n_wires, fitted_dBdz,
        gradient_nonlinearity, continuous_residual, wire_z, wire_I.
    """
    import radia as rad
    from radia.stream_function import (
        aca_tsvd, pseudo_inverse_solve, radia_field_kernel,
    )

    a = float(radius)
    L = float(length)
    Gz = float(gradient)
    dsv_h = float(dsv)

    rad.UtiDelAll()
    try:
        # --- ring basis -----------------------------------------------------
        z_ring = np.linspace(-L / 2, L / 2, n_rings)
        basis = [_ring(rad, a, zr, 1.0) for zr in z_ring]
        N = len(basis)

        # --- on-axis + OFF-AXIS volume target -------------------------------
        # The ring basis produces an axisymmetric field, so off-axis target
        # points at a SINGLE azimuth (here phi=0 -> y=0) fully constrain the
        # field over the DSV cylinder without breaking symmetry.  Sample
        # rho in {0, 0.4a, 0.7a} x z in the DSV and require Bz = Gz*z.
        zt = np.linspace(-dsv_h, dsv_h, 25)
        rho_t = np.array([0.0, 0.4 * a, 0.7 * a])    # target radii (phi = 0)
        RHO, ZT = np.meshgrid(rho_t, zt, indexing="ij")
        rho_flat = RHO.ravel()
        z_flat = ZT.ravel()
        obs = np.column_stack([rho_flat, np.zeros_like(rho_flat), z_flat])
        M = obs.shape[0]
        B_target = Gz * z_flat                        # Bz = Gz*z (rho-indep.)

        entry = radia_field_kernel(obs, basis, component=2, field="b")
        n_modes = min(M, N) if modes <= 0 else min(modes, M, N)
        res = aca_tsvd(M, N, entry, modes=n_modes, kmax=min(M, N),
                       aca_eps=aca_eps, method=method)
        k_use = res.modes if modes <= 0 else min(modes, res.modes)
        I = pseudo_inverse_solve(res, B_target, k_mode=k_use)  # ring currents

        # continuous-coil residual (sanity: A I should match Gz*z over the
        # full on-axis + off-axis target set)
        rad.UtiDelAll()
        cont = rad.ObjCnt([_ring(rad, a, zr, ii)
                           for zr, ii in zip(z_ring, I)])
        Bc = np.asarray(rad.Fld(cont, "b", obs.tolist()),
                        dtype=float).reshape(-1, 3)[:, 2]
        cont_res = float(np.linalg.norm(Bc - B_target)
                         / np.linalg.norm(B_target))
        rad.UtiDelAll()

        # --- contour psi(z) = cumulative current -> equal-current wires -----
        psi = np.concatenate([[0.0],
                              np.cumsum(0.5 * (I[1:] + I[:-1]))])  # cumtrapz
        wires = _contour(psi, z_ring, n_wires)

        # --- verify the DISCRETE wire coil on the axis ----------------------
        rad.UtiDelAll()
        coil = rad.ObjCnt([_ring(rad, a, zc, cur) for zc, cur in wires])
        zc = np.linspace(-dsv_h, dsv_h, 41)
        axis = np.column_stack([np.zeros_like(zc), np.zeros_like(zc), zc])
        Bz = np.asarray(rad.Fld(coil, "b", axis.tolist()),
                        dtype=float).reshape(-1, 3)[:, 2]
        G_fit = float(np.dot(zc, Bz) / np.dot(zc, zc))
        nonlin = float(np.max(np.abs(Bz - G_fit * zc))
                       / (np.max(np.abs(Bz)) + 1e-30))
        rad.UtiDelAll()

        return {
            "k_aca": int(res.k_aca),
            "modes": int(k_use),
            "n_wires": int(len(wires)),
            "fitted_dBdz": G_fit,
            "gradient_nonlinearity": nonlin,
            "continuous_residual": cont_res,
            "wire_z": [float(zw) for zw, _ in wires],
            "wire_I": [float(cur) for _, cur in wires],
        }
    finally:
        rad.UtiDelAll()


def main():
    _panels_dir = os.path.dirname(os.path.abspath(__file__))
    if _panels_dir not in sys.path:
        sys.path.insert(0, _panels_dir)
    from calc_common import calc_main

    parser = argparse.ArgumentParser(
        description="Stream-function cylindrical Gz gradient-coil designer")
    parser.add_argument("--radius", type=float, default=0.15,
                        help="Cylinder radius a [m] (default 0.15)")
    parser.add_argument("--length", type=float, default=0.5,
                        help="Cylinder length L [m] (default 0.5)")
    parser.add_argument("--gradient", type=float, default=1.0,
                        help="Target gradient dBz/dz Gz [T/m] (default 1.0)")
    parser.add_argument("--dsv", type=float, default=0.10,
                        help="DSV half-length [m] (target |z| <= dsv, "
                             "default 0.10)")
    parser.add_argument("--n-rings", type=int, default=80,
                        help="Number of candidate rings / basis size "
                             "(default 80)")
    parser.add_argument("--n-wires", type=int, default=16,
                        help="Number of equal-current wire rings (default 16)")
    parser.add_argument("--modes", type=int, default=0,
                        help="TSVD modes (0 = auto = min(M, N), "
                             "clamped to ACA+ rank)")
    parser.add_argument("--aca-eps", type=float, default=1e-8,
                        help="ACA+ stopping tolerance (default 1e-8)")
    parser.add_argument("--method", type=int, default=3, choices=(2, 3),
                        help="TSVD recompression method 2 or 3 (default 3)")

    def run(args):
        return design_stream_coil(
            radius=args.radius, length=args.length, gradient=args.gradient,
            dsv=args.dsv, n_rings=args.n_rings, n_wires=args.n_wires,
            modes=args.modes, aca_eps=args.aca_eps, method=args.method)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
