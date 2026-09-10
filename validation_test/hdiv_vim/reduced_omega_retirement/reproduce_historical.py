#!/usr/bin/env python
"""Is the TOSCA total/reduced split a pure change of variables?

HISTORICAL RESEARCH SCRIPT -- requires the original assets and pre-retirement
solver (for example source at f2f2fd799). The following hypothesis was NOT
confirmed: the exterior conventions differ. This is not a production recipe
or an upper-bound certification. See the adjacent README and result JSON.

Reading of the implementation: inside the total materials the mixed engine
solves for ``phi_total`` with

    H = H_s + grad(Phi_s) - grad(phi_total)      (total region)
    H = H_s - grad(phi_reduced)                  (reduced region)
    phi_total - phi_reduced = Phi_s              (on the interface)

If ``Phi_s`` lives in the SAME polynomial space as ``phi_total`` -- which is
what ``source_projection_order == order`` buys -- then substituting
``psi := phi_total - Phi_s`` gives ``H = H_s - grad(psi)`` in the total region
with ``psi = phi_reduced`` on the interface.  The pair therefore glues into ONE
globally continuous scalar of the same order, so the admissible set is
IDENTICAL to the plain reduced Omega--Omega formulation and the discrete
minimum must be the same number.

Consequences if true:
  - the 2.5-3 % relative harmonic norm of the iron Hodge split costs NOTHING;
    it is a gauge choice, absorbed by the unknown.
  - the split buys conditioning only (no H_s / grad(Omega) cancellation in
    high-permeability iron), not accuracy.
  - W_Omega is then a genuine variational upper bound, so on a FIXED mesh it
    must fall monotonically with p.

This script tests the claim directly: same mesh, same order, same shifted
energy, plain reduced Omega versus mixed total/reduced Omega.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import ngsolve as ng
import radia as rad

from radia.coil_builder import CoilBuilder
from radia.kelvin_identify_ngsolve import detect_kelvin_offset
from radia.kelvin_material import make_kelvin_aware_radia_H_s_cf
from radia.kelvin_solver import solve_magnetostatic_reduced_omega_kelvin
from radia.static_electromagnet import (
    StaticElectromagnetMixedDomain,
    solve_static_electromagnet_mixed_total_reduced_omega,
)

MU0 = 4.0e-7 * math.pi
NU0 = 1.0 / MU0

MIXED_DOMAIN = StaticElectromagnetMixedDomain(
    reduced_materials=("air",),
    total_materials=("iron", "kelvin"),
    nonlinear_materials=("iron",),
    reduced_total_interface="iron_air_interface",
)


def build_coil() -> int:
    radius = 0.0225
    straight_x = 0.050
    straight_y = 0.0625
    centre = np.array([0.0, 0.13125, 0.0])
    start = centre + np.array([0.5 * straight_x + radius, -0.5 * straight_y, 0.0])
    builder = (
        CoilBuilder(-2000.0)
        .set_start(start)
        .set_cross_section(0.035, 0.105)
        .add_straight(straight_y).add_arc(radius, 90.0)
        .add_straight(straight_x).add_arc(radius, 90.0)
        .add_straight(straight_y).add_arc(radius, 90.0)
        .add_straight(straight_x).add_arc(radius, 90.0)
    )
    if not builder.is_closed or builder.gap > 1e-12:
        raise RuntimeError("C-type CoilBuilder path is open by %.3e m" % builder.gap)
    return rad.ObjCnt(builder.to_radia(arc_max_segment_length=0.004))


def shifted_omega_energy(mesh, h_cf, mu_cf, h_source, kelvin_h_source,
                         mu_r, quadrature_orders):
    """W - 1/2 int nu_0 |B_s|^2, subtracted pointwise (same as the bracket)."""
    air = mesh.Materials("air")
    iron = mesh.Materials("iron")
    exterior = mesh.Materials("kelvin")
    rows = {}
    with ng.TaskManager():
        for q in quadrature_orders:
            omega_air = 0.5 * float(ng.Integrate(
                MU0 * ng.InnerProduct(h_cf - h_source, h_cf + h_source), mesh,
                definedon=air, order=q))
            omega_iron = 0.5 * float(ng.Integrate(
                mu_cf * ng.InnerProduct(h_cf, h_cf)
                - MU0 * ng.InnerProduct(h_source, h_source), mesh,
                definedon=iron, order=q))
            omega_kelvin = 0.5 * float(ng.Integrate(
                mu_cf * ng.InnerProduct(h_cf - kelvin_h_source,
                                        h_cf + kelvin_h_source), mesh,
                definedon=exterior, order=q))
            rows[str(q)] = {
                "omega_air": omega_air,
                "omega_iron": omega_iron,
                "omega_kelvin": omega_kelvin,
                "W_omega_hat_J": omega_air + omega_iron + omega_kelvin,
            }
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="plain reduced Omega versus mixed total/reduced Omega")
    parser.add_argument("--level", type=Path,
                        default=Path("C:/temp/radia-ctype-family/coarse"))
    parser.add_argument("--orders", nargs="+", type=int, default=[2])
    parser.add_argument("--mu-r", type=float, default=1000.0)
    parser.add_argument("--quadrature-orders", nargs="+", type=int,
                        default=[16, 22])
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    options = parser.parse_args()
    if options.threads > 0:
        ng.SetNumThreads(options.threads)

    report = json.loads((options.level / "mesh_result.json").read_text(encoding="utf-8"))
    mesh = ng.Mesh(str(options.level / "kelvin_domain.vol"))
    offset = np.asarray(detect_kelvin_offset(mesh), dtype=float)
    physical_center = np.asarray(report["kelvin_physical_center_m"], dtype=float)
    kelvin_center = tuple((physical_center + offset).tolist())
    radius = float(report["kelvin_radius_m"])

    rad.UtiDelAll()
    coil = build_coil()
    h_source = rad.RadiaField(coil, "h")
    kelvin_h_source = rad.KelvinRadiaFieldStrength(
        coil, kelvin_center, radius, tuple(physical_center.tolist()))
    h_source_global = make_kelvin_aware_radia_H_s_cf(
        mesh, coil, radius, kelvin_center,
        phys_center=tuple(physical_center.tolist()))

    # The Kelvin material's mu' is supplied by the transform itself and must
    # not be overridden here.
    mu_r_by_material = {"iron": float(options.mu_r), "air": 1.0}
    rows = []
    for order in options.orders:
        print("=== order %d %s" % (order, datetime.now().strftime("%H:%M:%S")),
              flush=True)

        started = time.perf_counter()
        with ng.TaskManager():
            plain = solve_magnetostatic_reduced_omega_kelvin(
                mesh, h_source_global, radius, kelvin_center,
                mu_r_by_material=mu_r_by_material, order=int(order),
                dirichlet_bbbnd="GND")
        plain_seconds = time.perf_counter() - started
        plain_h = h_source_global - ng.grad(plain["phi"])
        plain_rows = shifted_omega_energy(
            mesh, plain_h, plain["mu_cf"], h_source, kelvin_h_source,
            options.mu_r, options.quadrature_orders)

        started = time.perf_counter()
        with ng.TaskManager():
            mixed = solve_static_electromagnet_mixed_total_reduced_omega(
                mesh, h_source, MIXED_DOMAIN, radius, kelvin_center,
                order=int(order),
                linear_mu_r_by_material={"iron": float(options.mu_r)},
                bh_table=None,
                source_projection_order=int(order),
                source_potential_contract="total_hodge",
                kelvin_source_h=kelvin_h_source)
        mixed_seconds = time.perf_counter() - started
        mixed_rows = shifted_omega_energy(
            mesh, mixed["H_cf"], mixed["mu_cf"], h_source, kelvin_h_source,
            options.mu_r, options.quadrature_orders)

        contract = mixed["static_electromagnet_contract"]["source_trace"]
        q = str(options.quadrature_orders[-1])
        w_plain = plain_rows[q]["W_omega_hat_J"]
        w_mixed = mixed_rows[q]["W_omega_hat_J"]
        print("  plain  %.9f   mixed  %.9f   relative %.3e"
              % (w_plain, w_mixed, abs(w_mixed - w_plain) / abs(w_plain)),
              flush=True)
        print("  iron relative harmonic norm %.4f %%  (plain %.0f s, mixed %.0f s)"
              % (100.0 * float(contract.get("iron_relative_harmonic_norm", float("nan"))),
                 plain_seconds, mixed_seconds), flush=True)
        rows.append({
            "order": int(order),
            "plain_reduced_omega": plain_rows,
            "mixed_total_reduced_omega": mixed_rows,
            "relative_difference": abs(w_mixed - w_plain) / abs(w_plain),
            "iron_relative_harmonic_norm": contract.get(
                "iron_relative_harmonic_norm"),
            "source_projection_order": int(contract["projection_order"]),
            "plain_runtime_s": plain_seconds,
            "mixed_runtime_s": mixed_seconds,
        })

    options.output.write_text(json.dumps({
        "schema": "radia.research.reduced-vs-mixed-omega.v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "radia_version": rad.__version__,
        "level": str(options.level),
        "mesh_elements": int(mesh.ne),
        "mu_r": float(options.mu_r),
        "claim": "with source_projection_order == order the two admissible sets "
                 "coincide, so the shifted Omega energies must agree",
        "orders": rows,
    }, indent=2), encoding="utf-8")
    print("wrote %s" % options.output, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
