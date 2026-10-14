"""End-to-end workflow: stream function -> single-stroke CAD coil -> PEEC -> field.

Closes the design loop for a cylindrical Gz gradient coil:

  [1] stream function method   target Bz=Gz*z -> ring currents I(z) (ACA+TSVD)
  [2] equal-current contour     psi(z)=int I -> wire rings [(z_k, +/-I_w)]
  [3] single-stroke (1-pen)     coaxial rings -> ONE continuous helix polyline
                                (reverse winding handedness at current sign flips)
  [4] CAD                       sweep a round wire section along the helix -> STEP
  [5] PEEC                      STEP -> coil_from_cad filaments -> L, R, current
  [6] field                     exact Biot-Savart of the filaments -> B(r)
  [7] verify (loop closed)      manufactured-coil Bz(0,0,z) vs the design Gz*z

The single-stroke connection is what turns an abstract current density into a
real coil wound from one continuous conductor.  Field evaluation uses the exact
biot_savart.h_segments_batch (NOT ObjFlmCur, which has an orientation bug).

Run:  python demo_sf_to_peec_gz.py [--nrings N] [--nwires N] [--with-peec]
"""
import argparse
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

import radia as rad
from radia.stream_function import aca_tsvd, pseudo_inverse_solve, radia_field_kernel
from radia.biot_savart import h_segments_batch, MU0

A_CYL = 0.15       # cylinder radius [m]
L_CYL = 0.50       # cylinder length [m]
GZ = 1.0           # target gradient dBz/dz [T/m]
DSV = 0.10         # DSV half-length [m]
NSEG_RING = 120    # field-basis ring segments


def _ring(a, z0, cur):
    t = np.linspace(0.0, 2.0 * np.pi, NSEG_RING + 1)
    pts = [[a * np.cos(tt), a * np.sin(tt), z0] for tt in t]
    return rad.ObjFlmCur(pts, cur)


def design_rings(n_rings, n_wires, k_modes=0):
    """[1][2] SF solve + equal-current contour -> wire rings [(z, +/-I_w)]."""
    z_ring = np.linspace(-L_CYL / 2, L_CYL / 2, n_rings)
    basis = [_ring(A_CYL, zr, 1.0) for zr in z_ring]
    N = len(basis)
    zt = np.linspace(-DSV, DSV, 25)
    obs = np.column_stack([np.zeros(zt.size), np.zeros(zt.size), zt])
    M = obs.shape[0]
    entry = radia_field_kernel(obs, basis, component=2, field="b")
    modes = min(M, N) if k_modes <= 0 else min(k_modes, M, N)
    res = aca_tsvd(M, N, entry, modes=modes, kmax=min(M, N), aca_eps=1e-8)
    I = pseudo_inverse_solve(res, GZ * zt, k_mode=res.modes)
    rad.UtiDelAll()
    psi = np.concatenate([[0.0], np.cumsum(0.5 * (I[1:] + I[:-1]))])
    lo, hi = float(psi.min()), float(psi.max())
    dI = (hi - lo) / n_wires
    wires = []
    for lev in lo + (np.arange(n_wires) + 0.5) * dI:
        for j in range(len(z_ring) - 1):
            a0, a1 = psi[j] - lev, psi[j + 1] - lev
            if a0 == 0.0:
                a0 = 1e-300
            if a0 * a1 < 0.0:
                frac = a0 / (a0 - a1)
                zc = z_ring[j] + frac * (z_ring[j + 1] - z_ring[j])
                slope = psi[j + 1] - psi[j]
                wires.append((zc, dI * (1.0 if slope > 0 else -1.0)))
    return wires, dI, res.k_aca, min(M, N)


def single_stroke_helix(wires, a, nseg_per_turn=160, ramp_frac=0.18):
    """[3] Connect coaxial equal-current rings into ONE continuous polyline.

    Each ring is a near-flat helical turn at its z (so the field matches the SF
    design); the last `ramp_frac` of every turn is a smooth axial ramp
    (a *blended* crossover) to the next turn's z -- no sharp kink, so the CAD
    sweep stays clean.  Winding handedness encodes the current sign (CCW=+phi,
    CW=-phi); it reverses where the current sign flips.  A single current I_w
    flows along the whole path -> one continuous conductor.
    """
    wires = sorted(wires, key=lambda w: w[0])
    z_arr = [w[0] for w in wires]
    pts = []
    for k, (zk, ck) in enumerate(wires):
        d = 1.0 if ck > 0 else -1.0
        z_next = z_arr[k + 1] if k + 1 < len(wires) else zk
        for s in np.linspace(0.0, 1.0, nseg_per_turn, endpoint=False):
            phi = d * 2.0 * np.pi * s
            if s < 1.0 - ramp_frac:
                z = zk
            else:                                  # smooth axial ramp to z_next
                u = (s - (1.0 - ramp_frac)) / ramp_frac
                z = zk + (z_next - zk) * (0.5 - 0.5 * np.cos(np.pi * u))  # cosine blend
            pts.append([a * np.cos(phi), a * np.sin(phi), z])
    path = np.array(pts)
    I_w = abs(wires[0][1])
    return path, I_w


def field_on_axis(path, I_w, z):
    """[6] Exact Biot-Savart Bz(0,0,z) of the continuous wire path."""
    segs = np.stack([path[:-1], path[1:]], axis=1)        # (Nseg, 2, 3)
    obs = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])
    H = h_segments_batch(segs, obs, current=I_w)
    return MU0 * H[:, 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nrings", type=int, default=80)
    ap.add_argument("--nwires", type=int, default=16)
    ap.add_argument("--with-peec", action="store_true",
                    help="also run STEP export + coil_from_cad + PEEC")
    args = ap.parse_args()

    print("=== SF -> single-stroke coil -> (CAD/PEEC) -> field workflow (Gz) ===")

    # [1][2] stream function design -> equal-current wire rings
    wires, dI, k_aca, mn = design_rings(args.nrings, args.nwires)
    print(f"[1][2] SF design: {args.nrings} rings, ACA+ k_aca={k_aca} (of {mn}); "
          f"{len(wires)} equal-current turns, I_w={dI:.4g} A-turns")

    # [3] single-stroke helix
    path, I_w = single_stroke_helix(wires, A_CYL)
    length = float(np.sum(np.linalg.norm(np.diff(path, axis=0), axis=1)))
    print(f"[3] single-stroke path: {len(path)} points, {len(path)-1} segments, "
          f"wire length = {length:.3f} m (one continuous conductor)")

    # [6][7] field of the single-stroke coil vs the design target
    zc = np.linspace(-DSV, DSV, 41)
    Bz = field_on_axis(path, I_w, zc)
    G_fit = float(np.dot(zc, Bz) / np.dot(zc, zc))
    nonlin = float(np.max(np.abs(Bz - G_fit * zc)) / (np.max(np.abs(Bz)) + 1e-30))
    print(f"[6][7] single-stroke coil on axis: fitted dBz/dz = {G_fit:.4g}, "
          f"gradient nonlinearity over DSV = {nonlin:.3e}")
    print(f"       (compare: the abstract ring set is the SF design target Gz*z)")

    if args.with_peec:
        run_peec_chain(path, I_w)

    # plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams["pdf.fonttype"] = 42
        fig = plt.figure(figsize=(11, 4))
        ax = fig.add_subplot(1, 2, 1, projection="3d")
        ax.plot(path[:, 0] * 1e3, path[:, 1] * 1e3, path[:, 2] * 1e3, lw=0.6)
        ax.set_xlabel("x"); ax.set_ylabel("y")
        ax.text2D(0.02, 0.96, "(a)", transform=ax.transAxes, fontsize=10)
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.plot(zc * 1e3, Bz, "o-", ms=3, label="single-stroke coil Bz")
        ax2.plot(zc * 1e3, G_fit * zc, "k--", label="linear fit")
        ax2.set_xlabel("$z$ (mm)"); ax2.set_ylabel("$B_z$ (T)")
        ax2.text(0.02, 0.96, "(b)", transform=ax2.transAxes, fontsize=10)
        ax2.legend(); ax2.grid(alpha=0.3)
        out = HERE / "demo_sf_to_peec_gz.png"
        fig.tight_layout(); fig.savefig(out, dpi=110)
        print(f"Saved plot: {out.name}")
    except ImportError:
        print("(matplotlib not installed; skipped plot)")


def run_peec_chain(path, I_w):
    """[4] CAD STEP (sweep round wire along the helix) + [5] PEEC (L, R)."""
    import os
    import time

    # subsample for CAD/PEEC (keep the OCC spline + PEEC matrix tractable)
    stride = max(1, len(path) // 700)
    cpath = path[::stride]
    if not np.allclose(cpath[-1], path[-1]):
        cpath = np.vstack([cpath, path[-1]])

    # [4] CAD: sweep a circular wire cross-section along the smooth helix
    step_path = str(HERE / "sf_coil_gz.step")
    try:
        import build123d as b3d
        from radia.coil_from_cad import export_step_with_labels
        mm = 1000.0
        pmm = cpath * mm
        wire = b3d.Spline(*[b3d.Vector(float(x), float(y), float(z)) for x, y, z in pmm])
        # wire radius < half the min turn spacing, else adjacent tubes self-
        # intersect and OCC MakeSolid fails on the swept helix.
        dz_min = float(np.min(np.abs(np.diff(sorted({round(z, 6) for z in cpath[:, 2]})))) * mm) \
            if len({round(z, 6) for z in cpath[:, 2]}) > 1 else 4.0
        r_wire = max(0.2, min(1.0, 0.35 * dz_min))
        tan = pmm[1] - pmm[0]
        tan = tan / np.linalg.norm(tan)
        pl = b3d.Plane(origin=b3d.Vector(*pmm[0]), z_dir=b3d.Vector(*tan))
        sec = pl.location * b3d.Circle(radius=r_wire).face()
        t0 = time.time()
        try:
            solid = b3d.sweep(sec, path=wire, is_frenet=True)
            try:
                solid.label = "coil"
            except Exception:
                pass
            export_step_with_labels([solid], step_path)
            print(f"[4] CAD: round wire r={r_wire:.2f}mm swept along helix -> "
                  f"{os.path.basename(step_path)} ({os.path.getsize(step_path)//1024} KB, "
                  f"vol={solid.volume:.0f} mm^3, {time.time()-t0:.1f}s)")
        except Exception as se:
            # tight multi-turn helix can defeat OCC pipe MakeSolid; export the
            # centerline as a STEP WIRE (the windable conductor path) instead.
            b3d.export_step(wire, step_path)
            print(f"[4] CAD: solid sweep failed ({type(se).__name__}); exported "
                  f"centerline WIRE -> {os.path.basename(step_path)} "
                  f"({os.path.getsize(step_path)//1024} KB)")
    except Exception as e:
        print(f"[4] CAD STEP export FAILED: {type(e).__name__}: {str(e)[:120]}")

    # [5] PEEC: filaments built directly from the design centerline (we know it,
    # so no STEP centerline re-extraction needed) -> port impedance -> L, R.
    try:
        from radia.peec_matrices import PEECBuilder
        from radia.peec_topology import PEECCircuitSolver
        b = PEECBuilder()
        nodes = [b.add_node_at(float(p[0]), float(p[1]), float(p[2])) for p in cpath]
        w = hh = 4e-3
        for i in range(len(nodes) - 1):
            b.add_connected_segment(nodes[i], nodes[i + 1], w, hh, sigma=5.8e7)
        b.add_port(nodes[0], nodes[-1])
        topo = b.build_topology()
        solver = PEECCircuitSolver(topo)
        freq = 1.0e3
        t0 = time.time()
        Z = solver.compute_port_impedance(freq)
        Zc = complex(np.asarray(Z).ravel()[0]) if np.ndim(Z) else complex(Z)
        L = Zc.imag / (2 * np.pi * freq)
        R = Zc.real
        length = float(np.sum(np.linalg.norm(np.diff(cpath, axis=0), axis=1)))
        print(f"[5] PEEC ({len(nodes) - 1} filament segs, {length:.2f} m @ {freq:.0f} Hz, "
              f"{time.time() - t0:.1f}s): L = {L * 1e6:.3f} uH, R = {R * 1e3:.3f} mOhm")
        print(f"    (low freq -> ~uniform current I_w; the on-axis field above is "
              f"the PEEC-coil field at the design current)")
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"[5] PEEC FAILED: {type(e).__name__}: {str(e)[:120]}")


if __name__ == "__main__":
    main()
