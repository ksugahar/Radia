"""
Coupled axisymmetric induction heating (volumetric eddy currents + heat).

For a workpiece of revolution heated by coaxial coils this solves the
time-harmonic A_phi eddy-current problem on an (r, z) mesh of workpiece,
coils and air with sigma(T) and mu_r(T) per element, deposits the Joule
heat in the workpiece volume and advances the temperature with the enthalpy
Newton integrator; the EM problem is re-solved from the current temperature
every --em-every steps.  See radia.ih_axisym_coupled for the formulation.

This resolves the skin depth, so it remains valid across the Curie band,
where a surface source with a frozen |H_t| is not (see
validation_test/induction_heating/results/coupled_curie_cylinder_frozen_ht
.json).  It covers axisymmetric coil/workpiece arrangements; a coil with
leads or a non-axisymmetric part needs the 3D route.

Inputs
------
``--mesh``            2D (r, z) Netgen .vol, r = x >= 0, with the workpiece,
                      coil and air regions as materials.
``--workpiece``       Workpiece material name.
``--coil NAME=AMPS``  Coil region and its total peak current (ampere-turns);
                      repeat for several coils.
``--em-dirichlet``    Boundary expression with A_phi = 0; must include the
                      axis.  Other boundaries get H_tangential = 0.
``--em-material-table`` CSV ``T_C,sigma_S_m,mu_r`` (or --sigma/--mu-r).
Thermal options are those of calc_heat_axisym (constant or
--material-table, latent heat, convection and radiation boundaries of the
workpiece).

Output
------
JSON with the power and temperature histories, every EM solve's power
balance and skin depths, the thermal energy audit, thermal exposure and
optional case depth; the final temperature field with its sidecar.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time

import numpy as np

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import calc_main, progress, setup_paths  # noqa: E402
from calc_heat import (  # noqa: E402
    THERMAL_PRESETS,
    _resolve_boundary_role,
    _resolve_material,
)


def _log(msg):
    progress("IH_AXI", msg)


def _parse_coils(items):
    coils = {}
    for item in items or ():
        if "=" not in item:
            raise ValueError(f"--coil expects NAME=AMPS, got {item!r}")
        name, amps = item.split("=", 1)
        coils[name.strip()] = float(amps)
    if not coils:
        raise ValueError("at least one --coil NAME=AMPS is required")
    return coils


def solve_ih_axisym_coupled(mesh_path, *, workpiece, coils, frequency,
                            em_dirichlet, em_material_table="", sigma=None,
                            mu_r=None, allow_em_extrapolation=False,
                            material="steel", rho=None, cp=None, k=None,
                            material_table="", latent_heat=0.0,
                            latent_range=None, allow_table_extrapolation=False,
                            h_conv=0.0, t_ext=20.0, convection_boundaries="",
                            emissivity=0.0, radiation_boundaries="",
                            t_initial=20.0, dt=0.1, t_end=1.0, em_every=1,
                            em_order=2, thermal_order=2, newton_tol=1e-3,
                            exposure_thresholds=(), temperature_limit=None,
                            depth_threshold=None, depth_boundaries="",
                            depth_span=None, temperature_output="",
                            _mesh=None):
    setup_paths()
    t0 = time.perf_counter()
    from ngsolve import Mesh
    import ih_heat_transient
    import ih_thermal
    import ih_thermal_material
    import ih_thermal_post
    from radia import ih_axisym_coupled as C

    if _mesh is None and not os.path.isfile(mesh_path):
        return {"error": f"--mesh not found: {mesh_path}"}
    mesh = Mesh(mesh_path) if _mesh is None else _mesh
    if mesh.dim != 2:
        return {"error": "the coupled axisymmetric solver needs a 2D (r, z) "
                         "mesh"}
    try:
        if em_material_table:
            em_mat = C.EMMaterialTable.from_csv(
                em_material_table, allow_extrapolation=allow_em_extrapolation)
        else:
            if sigma is None or mu_r is None:
                raise ValueError("give --em-material-table or both --sigma "
                                 "and --mu-r")
            em_mat = C.EMMaterialTable.constant(float(sigma), float(mu_r))
        rho_v, cp_v, k_v = _resolve_material(material, rho, cp, k)
        if material_table:
            th_mat = ih_thermal_material.ThermalMaterial.from_csv(
                material_table, rho=rho_v, latent_heat=latent_heat,
                latent_range=latent_range,
                allow_extrapolation=allow_table_extrapolation)
        else:
            if latent_heat:
                raise ValueError("--latent-heat needs --material-table")
            th_mat = ih_thermal_material.ThermalMaterial.constant(
                rho_v, cp_v, k_v)
        conv_sel, conv_names = _resolve_boundary_role(
            mesh, convection_boundaries, "--convection-boundaries",
            required=float(h_conv) != 0.0)
        rad_sel, rad_names = _resolve_boundary_role(
            mesh, radiation_boundaries, "--radiation-boundaries",
            required=float(emissivity) != 0.0)
    except (ValueError, OSError) as exc:
        return {"error": str(exc)}
    _log(f"MESH:{mesh.ne} elements materials={sorted(set(mesh.GetMaterials()))}"
         f" coils={coils} f={frequency:g} Hz")
    boundaries = ih_heat_transient.HeatBoundaryTerms(
        convection=conv_sel if float(h_conv) else "", h_conv=float(h_conv),
        t_ext=float(t_ext), radiation=rad_sel if float(emissivity) else "",
        emissivity=float(emissivity))

    def on_step(step, entry):
        _log(f"STEP:{step} t={entry['t_s']:.4g}s P={entry['P_W']:.4e}W "
             f"T_max={entry['T_max_nodal_C']:.1f}C")

    try:
        result, gfT, em = C.run_coupled(
            mesh, frequency=frequency, workpiece=workpiece, coils=coils,
            dirichlet=em_dirichlet, em_material=em_mat,
            thermal_material=th_mat, boundaries=boundaries, dt=dt,
            t_end=t_end, t_initial=t_initial, em_every=em_every,
            em_order=em_order, thermal_order=thermal_order,
            newton_tol_K=newton_tol, on_step=on_step)
    except (ValueError, RuntimeError) as exc:
        return {"error": str(exc)}

    thresholds = sorted(set(float(t) for t in (exposure_thresholds or ())))
    if temperature_limit is not None:
        thresholds = sorted(set(thresholds) | {float(temperature_limit)})
    result["thermal_exposure"] = _exposure_on_workpiece(
        mesh, gfT, thresholds, workpiece, ih_thermal_post)
    if temperature_limit is not None:
        result["temperature_limit"] = ih_thermal_post.limit_check(
            result["thermal_exposure"], temperature_limit)
    if depth_threshold is not None:
        names = [b for b in str(depth_boundaries).split("|") if b.strip()]
        if not names:
            return {"error": "--depth-threshold needs --depth-boundaries"}
        result["case_depth"] = ih_thermal_post.case_depth(
            mesh, gfT, float(depth_threshold), boundary_names=names,
            span=depth_span, region=workpiece)
    if temperature_output:
        if _mesh is not None and not os.path.isfile(mesh_path):
            return {"error": "--temperature-output needs the mesh on disk"}
        gfT.Save(temperature_output)
        ih_thermal.write_field_sidecar(
            temperature_output, mesh_path=mesh_path, mesh=mesh,
            fes_order=int(thermal_order),
            quantity=ih_thermal.TEMPERATURE_QUANTITY,
            unit=ih_thermal.TEMPERATURE_UNIT,
            extra={"producer": "calc_ih_axisym_coupled",
                   "t_end_s": float(t_end)}, definedon=workpiece)
        result["temperature_output"] = temperature_output
    result["t_total_s"] = round(time.perf_counter() - t0, 2)
    return result


def _exposure_on_workpiece(mesh, gfT, thresholds, workpiece, post):
    """thermal_exposure restricted to the workpiece region."""
    return post.thermal_exposure(mesh, gfT, thresholds, axisymmetric=True,
                                 region=workpiece)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mesh", required=True)
    p.add_argument("--workpiece", required=True)
    p.add_argument("--coil", action="append", default=[],
                   help="NAME=AMPS (total peak ampere-turns); repeatable")
    p.add_argument("--frequency", type=float, required=True)
    p.add_argument("--em-dirichlet", required=True,
                   help="boundaries with A_phi = 0 (must include the axis)")
    p.add_argument("--em-material-table", default="")
    p.add_argument("--sigma", type=float, default=None)
    p.add_argument("--mu-r", type=float, default=None)
    p.add_argument("--allow-em-extrapolation", action="store_true")
    p.add_argument("--material", default="steel",
                   choices=list(THERMAL_PRESETS) + ["custom"])
    p.add_argument("--rho", type=float, default=None)
    p.add_argument("--cp", type=float, default=None)
    p.add_argument("--k", type=float, default=None)
    p.add_argument("--material-table", default="")
    p.add_argument("--latent-heat", type=float, default=0.0)
    p.add_argument("--latent-range", default="")
    p.add_argument("--allow-table-extrapolation", action="store_true")
    p.add_argument("--h-conv", type=float, default=0.0)
    p.add_argument("--t-ext", type=float, default=20.0)
    p.add_argument("--convection-boundaries", default="")
    p.add_argument("--emissivity", type=float, default=0.0)
    p.add_argument("--radiation-boundaries", default="")
    p.add_argument("--t-initial", type=float, default=20.0)
    p.add_argument("--dt", type=float, default=0.1)
    p.add_argument("--t-end", type=float, default=1.0)
    p.add_argument("--em-every", type=int, default=1)
    p.add_argument("--em-order", type=int, default=2)
    p.add_argument("--thermal-order", type=int, default=2)
    p.add_argument("--newton-tol", type=float, default=1e-3)
    p.add_argument("--exposure-thresholds", default="")
    p.add_argument("--temperature-limit", type=float, default=None)
    p.add_argument("--depth-threshold", type=float, default=None)
    p.add_argument("--depth-boundaries", default="")
    p.add_argument("--depth-span", type=float, default=None)
    p.add_argument("--temperature-output", default="")

    def run(a):
        try:
            coils = _parse_coils(a.coil)
            thr = [float(t) for t in a.exposure_thresholds.split(",")
                   if t.strip()]
            lr = (tuple(float(t) for t in a.latent_range.split(","))
                  if a.latent_range else None)
        except ValueError as exc:
            return {"error": str(exc)}
        return solve_ih_axisym_coupled(
            a.mesh, workpiece=a.workpiece, coils=coils,
            frequency=a.frequency, em_dirichlet=a.em_dirichlet,
            em_material_table=a.em_material_table, sigma=a.sigma,
            mu_r=a.mu_r, allow_em_extrapolation=a.allow_em_extrapolation,
            material=a.material, rho=a.rho, cp=a.cp, k=a.k,
            material_table=a.material_table, latent_heat=a.latent_heat,
            latent_range=lr,
            allow_table_extrapolation=a.allow_table_extrapolation,
            h_conv=a.h_conv, t_ext=a.t_ext,
            convection_boundaries=a.convection_boundaries,
            emissivity=a.emissivity,
            radiation_boundaries=a.radiation_boundaries,
            t_initial=a.t_initial, dt=a.dt, t_end=a.t_end,
            em_every=a.em_every, em_order=a.em_order,
            thermal_order=a.thermal_order, newton_tol=a.newton_tol,
            exposure_thresholds=thr, temperature_limit=a.temperature_limit,
            depth_threshold=a.depth_threshold,
            depth_boundaries=a.depth_boundaries, depth_span=a.depth_span,
            temperature_output=a.temperature_output)

    calc_main(run, p)


if __name__ == "__main__":
    main()
