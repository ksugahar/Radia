r"""Electric-machine demo: torque ripple tau(theta), physical scaling, and SKEW
cancellation -- stitching radia-ngsolve's open motor-FEA capabilities into one
runnable, analytic-gated showcase.

A diametric PM rotor in a 2-fold SALIENT iron stator gives a reluctance torque
tau(theta) ~ T0 sin(2 theta).  The rotor is rotated by rotating the magnetisation
(magnet phi_deg, no remesh).  This script then demonstrates the three 2D->3D motor
post-processors added from the open ONELAB/GetDP study:

  1. tau(theta) via the weighted-Maxwell-stress (eggshell) torque on an air-gap band.
  2. MachineScaling : 2D per-unit-depth torque -> physical N*m for a real stack
     length (ONELAB AxialLength x SymmetryFactor).
  3. skew_average   : the torque a continuously-SKEWED machine produces -- on the
     REAL swept curve.  A harmonic of order n is reduced by sinc(n*skew/2) =
     skew_factor(n, skew) (validated here, FE <-> analytic), and a half-period skew
     NULLS the order-2 ripple.

Pure radia-ngsolve (open-source); analytic / self-consistency gated -- no commercial
tool needed to run or to validate.

    python cogging_skew_demo.py
"""
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ngsolve import Mesh, CoefficientFunction, grad, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, MoveTo, Glue, X, Y

from radia_mcp.radia_ngsolve.solve import (reluctivity, solve_planar_magnetostatic,
                                           skew_average, skew_factor)
from radia_mcp.radia_ngsolve.force import eggshell_torque_2d
from radia_mcp.radia_ngsolve.machine_scaling import MachineScaling

MU0 = 4e-7 * math.pi
HC = 1.0 / MU0                       # Br ~ 1 T rigid PM (mu_r = 1 recoil)
Rr, Rs_in, Rs_out, BOX = 10.0, 13.0, 22.0, 70.0       # mm
RB_IN, RB_OUT = 10.6, 12.4           # eggshell band inside the air gap [mm]
HERE = Path(__file__).resolve().parent
VALIDATION_RESULT = (
    HERE.parents[1]
    / "validation_test"
    / "electric_machine"
    / "cogging_skew_demo_results.json"
)


def mesh_salient(Wx=44.0, Wy=26.0):
    """Diametric-PM rotor in a 2-fold salient iron stator (iron on the +-x sides)."""
    rotor = WorkPlane().Circle(0, 0, Rr).Face(); rotor.faces.name = "rotor"
    disk_in = WorkPlane().Circle(0, 0, Rs_in).Face()
    airgap = disk_in - rotor; airgap.faces.name = "air"
    rect = MoveTo(-Wx / 2, -Wy / 2).Rectangle(Wx, Wy).Face()
    iron = rect - disk_in; iron.faces.name = "stator"
    box = MoveTo(-BOX, -BOX).Rectangle(2 * BOX, 2 * BOX).Face()
    outer = box - rect; outer.faces.name = "air"
    for sel in (outer.edges.Max(X), outer.edges.Min(X),
                outer.edges.Max(Y), outer.edges.Min(Y)):
        sel.name = "outer"
    rotor.faces.maxh = 1.5; airgap.faces.maxh = 0.7; iron.faces.maxh = 2.0
    return Mesh(OCCGeometry(Glue([rotor, airgap, iron, outer]), dim=2).GenerateMesh(maxh=6.0))


def torque_curve(mesh, thetas_deg):
    """2D eggshell torque tau(theta) [per-unit-depth, mm units] over the rotor sweep."""
    nu = reluctivity(mesh, {"stator": 1000.0})        # linear iron; air/rotor -> nu0
    out = []
    with TaskManager():
        for th in thetas_deg:
            A = solve_planar_magnetostatic(mesh, nu, magnets={"rotor": (HC, th)}, order=2)
            B = CoefficientFunction((grad(A)[1], -grad(A)[0]))
            out.append(eggshell_torque_2d(B, mesh, (0.0, 0.0), RB_IN, RB_OUT))
    return np.array(out)


def _runtime_versions():
    versions = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
    }
    for name in ("ngsolve", "netgen", "radia_mcp"):
        try:
            module = __import__(name)
            versions[name] = getattr(module, "__version__", "unknown")
        except Exception as exc:  # pragma: no cover - diagnostic only
            versions[name] = f"unavailable: {exc}"
    return versions


def run_demo(verbose=False):
    """Run the analytic-gated cogging/skew demo and return serializable results."""
    thetas = list(range(0, 360, 30))                  # 12 rotor angles, full revolution
    tau2d = torque_curve(mesh_salient(), thetas)

    # 1) shape: dominant order 2, zero mean (reluctance ripple of a 2-fold stator)
    C = np.fft.rfft(tau2d - tau2d.mean())
    order = int(np.argmax(np.abs(C[1:]))) + 1
    mean_frac = abs(tau2d.mean()) / np.abs(tau2d).max()
    if verbose:
        print(f"1) tau(theta): dominant order = {order} (expect 2), "
              f"mean/amp = {100*mean_frac:.1f}% (expect ~0)")
    assert order == 2 and mean_frac < 0.15

    # 2) physical scaling: 2D mm-mesh per-depth torque -> N*m for a 50 mm stack, full machine
    sc = MachineScaling(stack_length=0.050, symmetry_factor=1.0, length_unit_m=1e-3)
    amp_Nm = sc.torque(np.abs(tau2d).max())
    if verbose:
        print(f"2) MachineScaling(50 mm stack): peak |tau| = {amp_Nm*1e3:.4f} mN*m "
              f"(= 2D x length_unit^2 x stack)")

    # 3) skew: the FE skew_average must match the analytic skew_factor on the REAL curve
    if verbose:
        print("3) skew_average(tau) vs analytic skew_factor (order-2 ripple):")
    amp0 = math.sqrt(2.0 * np.mean((tau2d - tau2d.mean()) ** 2))
    skew_rows = []
    for skew_deg in (30.0, 90.0, 180.0):
        sk = math.radians(skew_deg)
        ts = skew_average(np.radians(thetas), tau2d - tau2d.mean(), sk, n_slices=0)
        amp = math.sqrt(2.0 * np.mean(ts ** 2))
        predicted = abs(skew_factor(2, sk)) * amp0
        tag = "  <- NULLS the order-2 ripple" if skew_deg == 180.0 else ""
        ratio = amp / amp0
        predicted_ratio = abs(skew_factor(2, sk))
        if verbose:
            print(f"   skew {skew_deg:5.0f} deg : amp {ratio:6.3f} x  "
                  f"(skew_factor predicts {predicted_ratio:6.3f}){tag}")
        assert math.isclose(amp, predicted, rel_tol=2e-2, abs_tol=1e-9 * amp0)
        skew_rows.append({
            "skew_deg": skew_deg,
            "amplitude_ratio": float(ratio),
            "predicted_ratio": float(predicted_ratio),
            "amplitude_per_depth": float(amp),
        })

    result = {
        "schema": "radia.validation.electric_machine.cogging_skew_demo.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "versions": _runtime_versions(),
        "checks": {
            "dominant_order_is_2": bool(order == 2),
            "mean_fraction_lt_0p15": bool(mean_frac < 0.15),
            "skew_matches_analytic": bool(all(
                math.isclose(
                    r["amplitude_ratio"],
                    r["predicted_ratio"],
                    rel_tol=2e-2,
                    abs_tol=1e-12,
                )
                for r in skew_rows
            )),
        },
        "theta_deg": thetas,
        "torque_2d_per_depth": [float(v) for v in tau2d],
        "dominant_order": int(order),
        "mean_fraction": float(mean_frac),
        "peak_torque_Nm": float(amp_Nm),
        "peak_torque_mNm": float(amp_Nm * 1e3),
        "skew_rows": skew_rows,
    }
    if verbose:
        print("\n[OK] electric-machine demo: tau(theta) reluctance ripple (order 2, zero mean), "
              "MachineScaling -> physical N*m, and skew_average == analytic skew_factor on the "
              "real FE curve (a half-period skew nulls the ripple). Pure radia-ngsolve, analytic-gated.")
    return result


def write_results_json(result, path=VALIDATION_RESULT):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def main():
    result = run_demo(verbose=True)
    out = write_results_json(result)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
