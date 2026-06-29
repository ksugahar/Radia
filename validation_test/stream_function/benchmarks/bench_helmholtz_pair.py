"""Benchmark vs analytical Helmholtz pair (uniform Bz).

Target: uniform Bz over a spherical DSV of radius r centred on the axis.
Reference: two coaxial circular loops of radius a separated by distance
a (= Helmholtz pair).  Analytical Bz on the axis:

    Bz(z) = (mu_0 * I * a^2 / 2) * [
        1 / (a^2 + (z + a/2)^2)^{3/2}
      + 1 / (a^2 + (z - a/2)^2)^{3/2}
    ]

The Maxwell-pair geometry was designed to make d^2 Bz/dz^2 = 0 at the
midplane, giving uniform Bz to leading order over a small DSV.

We use this as the LITERATURE-EQUIVALENT BASELINE for our planar
uniform-Bz SF coil design.  The DSV uniformity of a Helmholtz pair is
well known (~1% over a sphere of radius a/4).  Our planar SF design
should reach the SAME OR BETTER uniformity at the same target spec.

This is the FIRST shipped benchmark.  Other industry benchmarks are
documented in docs/stream_function/benchmarks.md until they have runnable
validation_test implementations.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "src"))  # for radia
from benchmark_framework import (
    Benchmark, BenchmarkSpec, add_common_args,
)
from radia.biot_savart import h_segments_batch, MU0


def helmholtz_pair_bz_axis(I, a, z_axis):
    """Analytical Bz on the axis of a Helmholtz pair.

    Two loops at z = +/- a/2, radius a, current I.
    """
    z_plus = z_axis + a / 2.0
    z_minus = z_axis - a / 2.0
    r2_plus = a * a + z_plus * z_plus
    r2_minus = a * a + z_minus * z_minus
    Bz = (MU0 * I * a * a / 2.0) * (
        1.0 / r2_plus ** 1.5 + 1.0 / r2_minus ** 1.5
    )
    return Bz


def helmholtz_pair_bz_3d(I, a, obs_points):
    """Helmholtz pair Bz at arbitrary 3D points by Biot-Savart over the
    two ring sources.  Used to compute the 3D DSV uniformity."""
    n_seg = 256
    t = np.linspace(0, 2 * np.pi, n_seg + 1)
    Bz_total = np.zeros(obs_points.shape[0])
    for sign in (-1.0, +1.0):
        ring_pts = np.column_stack([
            a * np.cos(t), a * np.sin(t),
            (sign * a / 2.0) * np.ones_like(t),
        ])
        segs = np.stack([ring_pts[:-1], ring_pts[1:]], axis=1)
        H = h_segments_batch(segs, obs_points, current=I)
        Bz_total = Bz_total + MU0 * H[:, 2]
    return Bz_total


class HelmholtzPairBenchmark(Benchmark):
    """Compare our planar SF uniform-Bz design vs a Helmholtz pair."""

    def __init__(self, a=0.10, dsv_r=0.025, B0=1.0e-3,
                 n_target_per_axis=5):
        """
        Parameters
        ----------
        a : Helmholtz pair radius [m] (also half the loop separation)
        dsv_r : DSV radius [m] -- uniformity is evaluated here
        B0 : target Bz [T] (we scale Helmholtz current to match)
        n_target_per_axis : sphere is approximated by a cube grid
        """
        self._a = a
        self._dsv_r = dsv_r
        self._B0 = B0
        self._n_target = n_target_per_axis
        # Helmholtz pair current to give Bz = B0 at origin
        Bz_at_origin = helmholtz_pair_bz_axis(1.0, a, np.array([0.0]))[0]
        self._I_helm = B0 / Bz_at_origin
        # DSV grid: 3D cube cropped to sphere
        g = np.linspace(-dsv_r, dsv_r, n_target_per_axis)
        Xt, Yt, Zt = np.meshgrid(g, g, g, indexing="ij")
        pts = np.column_stack([Xt.ravel(), Yt.ravel(), Zt.ravel()])
        rad = np.linalg.norm(pts, axis=1)
        self._dsv_obs = pts[rad <= dsv_r + 1e-12]

    @property
    def spec(self) -> BenchmarkSpec:
        return BenchmarkSpec(
            name="helmholtz_pair",
            reference="Maxwell (1873) / Helmholtz (classic); analytical baseline",
            target_kind="uniform_Bz",
            target_params={
                "B0_T": self._B0,
                "dsv_radius_m": self._dsv_r,
                "n_target": int(self._dsv_obs.shape[0]),
            },
            geometry={
                "type": "helmholtz_pair",
                "loop_radius_a_m": self._a,
                "loop_separation_m": self._a,
                "current_A": self._I_helm,
            },
        )

    @property
    def baseline(self) -> dict:
        """Helmholtz pair: compute uniformity over our same DSV."""
        Bz_helm = helmholtz_pair_bz_3d(
            self._I_helm, self._a, self._dsv_obs)
        target = self._B0 * np.ones_like(Bz_helm)
        Bz_mean = float(np.mean(Bz_helm))
        rms = float(np.linalg.norm(Bz_helm - target) /
                     (np.linalg.norm(target) + 1e-30))
        p2p = (float(np.max(Bz_helm)) - float(np.min(Bz_helm))) / \
              (abs(Bz_mean) + 1e-30)
        # Helmholtz pair wire length = 2 * 2*pi*a
        wire_length = float(2.0 * 2.0 * np.pi * self._a)
        return {
            "rms": rms,
            "p2p_mean": p2p,
            "Bz_mean_T": Bz_mean,
            "wire_length_m": wire_length,
            "n_contours": 2,
            "current_A": self._I_helm,
        }

    def solve(self) -> dict:
        """Run our planar SF design with the same DSV target."""
        # Import the SF design pipeline
        from demo_planar_uniform_fem_psi_advanced import (
            solve_and_evaluate,
        )

        target_args = {
            "half": self._dsv_r,
            "z": 1.5 * self._a,        # target plane above the source
            "n": self._n_target,
        }
        # Source plane large enough to embrace the projection of the DSV
        plane_half = 2.5 * self._a
        result = solve_and_evaluate(
            plane_half=plane_half,
            maxh=0.025,
            order=3,                    # session-discovered sweet spot
            target_args=target_args,
            B0=self._B0,
            regularize="h1",
            sigma_cf_expr=None,
            jmax=None,
            surface_z_fn=None,
            n_sample=81,
            nlevels=12,
        )
        return {
            "rms": result["rms"],
            "p2p_mean": result["p2p_mean"],
            "Bz_mean_T": result["Bz_mean"],
            "wire_length_m": result["length"],
            "n_contours": result["n_wire"],
            "current_A": result["I_w"],
        }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=float, default=0.10,
                    help="Helmholtz pair radius [m]")
    ap.add_argument("--dsv-r", type=float, default=0.025,
                    help="DSV radius [m]")
    ap.add_argument("--B0", type=float, default=1.0e-3,
                    help="target Bz [T]")
    ap.add_argument("--n-target", type=int, default=5)
    add_common_args(ap)
    args = ap.parse_args()

    bench = HelmholtzPairBenchmark(
        a=args.a, dsv_r=args.dsv_r, B0=args.B0,
        n_target_per_axis=args.n_target)

    print(f"=== Benchmark: {bench.spec.name} ===")
    print(f"  reference: {bench.spec.reference}")
    print(f"  target: uniform Bz = {args.B0*1e3:.2f} mT over"
          f" DSV r = {args.dsv_r*1e3:.0f} mm")
    print(f"  Helmholtz pair: a = {args.a*1e3:.0f} mm,"
          f" current = {bench._I_helm:.3g} A")
    print()

    result = bench.run(json_out=args.json)
    print(f"=== Baseline (Helmholtz pair) ===")
    for k, v in result.published_baseline.items():
        print(f"  {k}: {v}")
    print()
    print(f"=== Our SF design ===")
    for k, v in result.our_result.items():
        print(f"  {k}: {v}")
    print()
    print(f"=== Ratio (us / them) ===")
    if result.ratio:
        for k, v in result.ratio.items():
            verdict = "BETTER" if v < 1.0 else "WORSE"
            print(f"  {k}: {v:.3f}  ({verdict})")
    print(f"  elapsed: {result.elapsed_s:.1f} s")


if __name__ == "__main__":
    main()
