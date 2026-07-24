"""Run public TEAM Workshop Problem 36 through the radia-ih validation lane."""

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import socket
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from netgen.occ import Glue, MoveTo, OCCGeometry
from ngsolve import (
    BilinearForm,
    CF,
    Conj,
    GridFunction,
    H1,
    Integrate,
    L2,
    LinearForm,
    Mesh,
    TaskManager,
    VOL,
    dx,
    ds,
    grad,
    sqrt,
    x,
)

from radia.axifem import H1Henrotte
from radia_mcp.ih.team36_gate import (
    EXPECTED_EXCITATION,
    EXPECTED_GEOMETRY_M,
    SCHEMA,
    SOURCE_URL,
    evaluate_team36_artifact,
    identity_digest,
)


MU0 = 4.0e-7 * math.pi
SIGMA_SB = 5.670374419e-8

H_TABLE = np.array([0, 500, 1000, 1500, 2000, 2500, 3000, 4000,
                    8000, 15900, 23900, 39900, 79700, 159400, 239100,
                    318800, 358700, 398500, 477000, 557000], dtype=float)
MU20_TABLE = np.array([0, 350, 500, 600, 525, 450, 390, 305, 164,
                       89.2, 62.3, 39.7, 21, 11.1, 7.8, 6.1, 5.5, 5.1,
                       4.4, 3.9], dtype=float)
RHOEL_T = np.array([0, 100, 200, 300, 400, 500, 600, 700, 800, 900,
                    1000, 1100, 1200, 1400, 1470, 1500], dtype=float)
RHOEL = np.array([1.77e-7, 2.38e-7, 3.12e-7, 4.00e-7, 5.10e-7,
                  6.35e-7, 7.55e-7, 9.50e-7, 1.11e-6, 1.16e-6,
                  1.19e-6, 1.22e-6, 1.24e-6, 1.25e-6, 1.30e-6,
                  1.30e-6], dtype=float)
K_T = np.array([0, 100, 200, 300, 400, 500, 600, 700, 750, 800, 900,
                1000, 1100, 1200, 1400, 1470, 1800], dtype=float)
K_TABLE = np.array([48.1, 48.1, 46.5, 44.0, 41.0, 38.5, 36.0, 31.4,
                    28.5, 26.7, 25.9, 26.7, 28.0, 29.8, 35.0, 39.0,
                    39.0], dtype=float)
CP_T = np.array([0, 50, 100, 200, 300, 400, 500, 600, 650, 700, 710,
                 720, 730, 747, 760, 770, 775, 780, 785, 787, 790, 800,
                 850, 900, 1000, 1100, 1200, 1300, 1334, 1400, 1430,
                 1470, 1500, 1600, 1700, 1800], dtype=float)
CP_TABLE = np.array([481.06, 486.09, 494.04, 522.93, 561.03, 599.13,
                     669.89, 720.13, 749.86, 808.89, 870.02, 919.84,
                     1170.0, 1470.0, 1620.2, 1699.8, 1660.2, 1630.7,
                     1589.7, 1520.7, 1459.7, 1353.0, 979.4, 766.15,
                     658.02, 655.97, 661.93, 709.69, 711.72, 2042.6,
                     3581.6, 1562.8, 776.17, 820.0, 890.11, 900.16],
                    dtype=float)


def interp(table_x, table_y, value):
    return float(np.interp(float(value), table_x, table_y))


def permeability_r(h_a_m, t_c):
    tc = 770.0
    c = 20.0
    t1 = tc + c * math.log(0.9)
    t2 = t1 + 0.1 * c * math.log(0.1)
    if t_c < t1:
        f_t = 1.0 - math.exp((t_c - tc) / c)
    else:
        f_t = math.exp(10.0 * (t2 - t_c) / c)
    return 1.0 + max(0.0, f_t) * interp(H_TABLE, MU20_TABLE, abs(h_a_m))


def set_edge_names(face, r_max, z_max):
    for edge in face.edges:
        cx = float(edge.center.x)
        cz = float(edge.center.y)
        if abs(cx) < 1.0e-10:
            edge.name = "axis"
        elif abs(cx - r_max) < 1.0e-10 or abs(cz - z_max) < 1.0e-10:
            edge.name = "outer"
        elif abs(cz) < 1.0e-10:
            edge.name = "symmetry"


def build_em_mesh(maxh_billet=0.0015, maxh_air=0.035):
    r_air, z_air = 0.16, 0.62
    domain = MoveTo(0.0, 0.0).Rectangle(r_air, z_air).Face()
    set_edge_names(domain, r_air, z_air)

    radial_cuts = (0.0, 0.018, 0.024, 0.027, 0.029, 0.030)
    layer_factors = (4.0, 3.0, 2.0, 1.25, 1.0)
    billet_layers = []
    for r0, r1, factor in zip(radial_cuts[:-1], radial_cuts[1:], layer_factors):
        layer = MoveTo(r0, 0.0).Rectangle(r1 - r0, 0.50).Face()
        layer.faces.name = "billet"
        layer.maxh = max(float(maxh_billet) * factor, r1 - r0)
        billet_layers.append(layer)

    coils = []
    for turn in range(10):
        z0 = 0.01 + 0.05 * turn
        coil = MoveTo(0.048, z0).Rectangle(0.02, 0.04).Face()
        coil.faces.name = "coil"
        coil.maxh = 0.008
        coils.append(coil)

    air = domain
    for layer in billet_layers:
        air = air - layer
    for coil in coils:
        air = air - coil
    air.faces.name = "air"
    air.maxh = maxh_air
    shape = Glue([air, *billet_layers, *coils])
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh_air))


def build_thermal_mesh(maxh=0.002):
    radial_cuts = (0.0, 0.017, 0.0235, 0.0265, 0.0285, 0.030)
    layer_factors = (5.0, 3.5, 2.5, 1.5, 1.0)
    layers = []
    for r0, r1, factor in zip(radial_cuts[:-1], radial_cuts[1:], layer_factors):
        layer = MoveTo(r0, 0.0).Rectangle(r1 - r0, 0.50).Face()
        layer.faces.name = "billet"
        layer.maxh = max(float(maxh) * factor, r1 - r0)
        for edge in layer.edges:
            cx = float(edge.center.x)
            cz = float(edge.center.y)
            if abs(cx) < 1.0e-10:
                edge.name = "axis"
            elif abs(cx - 0.03) < 1.0e-10:
                edge.name = "lateral"
            elif abs(cz - 0.50) < 1.0e-10:
                edge.name = "end"
            elif abs(cz) < 1.0e-10:
                edge.name = "symmetry"
        layers.append(layer)
    shape = Glue(layers)
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(
        maxh=max(5.0 * float(maxh), 0.012)))


def element_center(mesh, element):
    pts = [mesh.vertices[v.nr].point for v in element.vertices]
    return tuple(float(sum(p[i] for p in pts) / len(pts)) for i in range(2))


def element_material(mesh, element):
    material = element.mat
    if isinstance(material, str):
        return material
    return mesh.GetMaterials()[int(material) - 1]


def set_l2_element_values(gf, values):
    fes = gf.space
    for element, value in zip(gf.space.mesh.Elements(VOL), values):
        dofs = fes.GetDofNrs(element)
        gf.vec[dofs[0]] = float(value)


def eval_real(gf, mesh, point):
    value = gf(mesh(*point))
    return float(getattr(value, "real", value))


def map_temperature_to_em(gf_t, thermal_mesh, em_mesh):
    fes = L2(em_mesh, order=0)
    gf = GridFunction(fes)
    values = []
    for element in em_mesh.Elements(VOL):
        point = element_center(em_mesh, element)
        if element_material(em_mesh, element) == "billet":
            values.append(eval_real(gf_t, thermal_mesh, point))
        else:
            values.append(20.0)
    set_l2_element_values(gf, values)
    return gf, np.asarray(values)


def matrix_to_scipy(matrix, n):
    rows, cols, vals = matrix.COO()
    return sp.csr_matrix(
        (np.asarray(vals).real.astype(float),
         (np.asarray(rows, dtype=int), np.asarray(cols, dtype=int))),
        shape=(n, n),
    )


def solve_em(em_mesh, temperature_em, current_rms=3500.0, frequency=2000.0,
             max_iter=40, tol=2.0e-3, relax=0.45):
    omega = 2.0 * math.pi * frequency
    fes = H1Henrotte(em_mesh, order=1, complex=True, dirichlet="axis|outer")
    fes_l2 = L2(em_mesh, order=0)
    mu = GridFunction(fes_l2)
    sigma = GridFunction(fes_l2)
    elements = list(em_mesh.Elements(VOL))
    mu_values = []
    sigma_values = []
    for element in elements:
        material = element_material(em_mesh, element)
        t_c = eval_real(temperature_em, em_mesh, element_center(em_mesh, element))
        if material == "billet":
            mu_values.append(MU0 * permeability_r(1000.0, t_c))
            sigma_values.append(1.0 / interp(RHOEL_T, RHOEL, t_c))
        else:
            mu_values.append(MU0)
            sigma_values.append(0.0)
    set_l2_element_values(mu, mu_values)
    set_l2_element_values(sigma, sigma_values)

    trial, test = fes.TnT()
    a_m = BilinearForm(fes, symmetric=True)
    a_m += sigma * x * trial * test * dx
    a_m.Assemble()
    m = matrix_to_scipy(a_m.mat, fes.ndof)
    free = np.array([i for i in range(fes.ndof) if fes.FreeDofs()[i]], dtype=int)
    m_free = m[free[:, None], free[None, :]]
    m_free = 0.5 * (m_free + m_free.T)

    coil_area = 0.02 * 0.04
    j_peak = math.sqrt(2.0) * current_rms / coil_area
    jr = em_mesh.MaterialCF({"coil": j_peak}, default=0.0)
    rhs_form = LinearForm(fes)
    rhs_form += jr * x * test * dx
    rhs_form.Assemble()
    rhs = np.asarray(rhs_form.vec.FV().NumPy(), dtype=complex)[free]

    gfu = GridFunction(fes)
    solution = np.zeros(len(free), dtype=complex)
    b_cf = CF((grad(gfu)[0] + gfu / x, -grad(gfu)[1]))
    bmag_cf = sqrt((b_cf[0] * Conj(b_cf[0]) +
                    b_cf[1] * Conj(b_cf[1])).real + 1.0e-30)
    converged = False
    rel = math.inf
    for iteration in range(max_iter):
        nu = 1.0 / mu
        a_k = BilinearForm(fes, symmetric=True)
        a_k += nu * (1.0 / x) * (x * grad(trial)[0] + trial) \
            * (x * grad(test)[0] + test) * dx
        a_k += nu * x * grad(trial)[1] * grad(test)[1] * dx
        a_k.Assemble()
        k = matrix_to_scipy(a_k.mat, fes.ndof)
        k_free = k[free[:, None], free[None, :]]
        k_free = 0.5 * (k_free + k_free.T)
        if iteration == 0:
            diagonal = np.abs(k_free.diagonal()) + omega * np.abs(m_free.diagonal())
            active = diagonal > 1.0e-30
            print("EM matrix", fes.ndof, len(free), "inactive", int(np.count_nonzero(~active)))
            if not np.all(active):
                free = free[active]
                k_free = k_free[active][:, active]
                m_free = m_free[active][:, active]
                rhs = rhs[active]
                solution = np.zeros(len(free), dtype=complex)
        new_solution = spla.spsolve((k_free + 1j * omega * m_free).tocsc(), rhs)
        if iteration:
            solution = (1.0 - relax) * solution + relax * new_solution
        else:
            solution = new_solution
        gfu.vec.FV().NumPy()[:] = 0.0
        gfu.vec.FV().NumPy()[free] = solution

        updated = list(mu_values)
        for index, element in enumerate(elements):
            if element_material(em_mesh, element) != "billet":
                continue
            point = element_center(em_mesh, element)
            bmag = eval_real(bmag_cf, em_mesh, point)
            hmag = bmag / max(mu_values[index], MU0)
            t_c = eval_real(temperature_em, em_mesh, point)
            target = MU0 * permeability_r(hmag, t_c)
            updated[index] = (1.0 - relax) * mu_values[index] + relax * target
        old = np.asarray(mu_values)
        new = np.asarray(updated)
        billet_mask = np.array([element_material(em_mesh, e) == "billet" for e in elements])
        rel = float(np.linalg.norm(new[billet_mask] - old[billet_mask]) /
                    max(np.linalg.norm(new[billet_mask]), 1.0e-30))
        mu_values = updated
        set_l2_element_values(mu, mu_values)
        if iteration >= 2 and rel < tol:
            converged = True
            break

    q = GridFunction(fes_l2)
    q_values = []
    for index, element in enumerate(elements):
        if element_material(em_mesh, element) == "billet":
            point = element_center(em_mesh, element)
            a_value = gfu(em_mesh(*point))
            q_values.append(0.5 * sigma_values[index] * omega * omega * abs(a_value) ** 2)
        else:
            q_values.append(0.0)
    set_l2_element_values(q, q_values)
    p_em = float(Integrate(q * 2.0 * math.pi * x, em_mesh,
                           definedon=em_mesh.Materials("billet")).real)
    return {
        "A": gfu,
        "q": q,
        "power_w": p_em,
        "iterations": iteration + 1,
        "converged": converged,
        "relative_update": rel,
        "mu_r_min": float(np.min(np.asarray(mu_values)[billet_mask]) / MU0),
        "mu_r_max": float(np.max(np.asarray(mu_values)[billet_mask]) / MU0),
    }


def map_power_to_thermal(q_em, em_mesh, thermal_mesh, p_em):
    fes = L2(thermal_mesh, order=0)
    q_th = GridFunction(fes)
    values = []
    for element in thermal_mesh.Elements(VOL):
        point = element_center(thermal_mesh, element)
        values.append(max(0.0, eval_real(q_em, em_mesh, point)))
    set_l2_element_values(q_th, values)
    raw = float(Integrate(q_th * 2.0 * math.pi * x, thermal_mesh).real)
    scale = p_em / raw if raw > 0.0 else 0.0
    q_th.vec.data *= scale
    mapped = float(Integrate(q_th * 2.0 * math.pi * x, thermal_mesh).real)
    return q_th, raw, mapped, scale


def material_field(mesh, temperature, table_x, table_y):
    field = GridFunction(L2(mesh, order=0))
    values = [interp(table_x, table_y,
                     eval_real(temperature, mesh, element_center(mesh, element)))
              for element in mesh.Elements(VOL)]
    set_l2_element_values(field, values)
    return field


def thermal_step(mesh, temperature, q, dt, max_iter=50, tol=2.0e-5,
                 relax=0.55):
    fes = temperature.space
    u, v = fes.TnT()
    old = temperature.vec.CreateVector()
    old.data = temperature.vec
    old_gf = GridFunction(fes)
    old_gf.vec.data = temperature.vec
    guess = GridFunction(fes)
    guess.vec.data = temperature.vec
    weight = 2.0 * math.pi * x
    rel = math.inf
    for iteration in range(max_iter):
        cp = material_field(mesh, guess, CP_T, CP_TABLE)
        conductivity = material_field(mesh, guess, K_T, K_TABLE)
        a = BilinearForm(fes, symmetric=True)
        a += 7800.0 * cp * u * v * weight * dx
        a += dt * conductivity * grad(u) * grad(v) * weight * dx
        rhs = LinearForm(fes)
        rhs += 7800.0 * cp * old_gf * v * weight * dx
        rhs += dt * q * v * weight * dx
        for boundary, ambient in (("lateral", 70.0), ("end", 25.0)):
            tk = guess + 273.15
            ta = ambient + 273.15
            h_rad = 0.8 * SIGMA_SB * (tk + ta) * (tk * tk + ta * ta)
            h_total = 7.0 + h_rad
            a += dt * h_total * u * v * weight * ds(boundary)
            rhs += dt * h_total * ambient * v * weight * ds(boundary)
        with TaskManager():
            a.Assemble()
            rhs.Assemble()
        updated = GridFunction(fes)
        updated.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * rhs.vec
        old_arr = guess.vec.FV().NumPy().copy()
        solved_arr = updated.vec.FV().NumPy()
        new_arr = (1.0 - relax) * old_arr + relax * solved_arr
        rel = float(np.linalg.norm(new_arr - old_arr) /
                    max(np.linalg.norm(new_arr), 1.0e-30))
        guess.vec.FV().NumPy()[:] = new_arr
        if iteration >= 1 and rel < tol:
            break
    temperature.vec.data = guess.vec
    return iteration + 1, rel


def mesh_topology_sha256(mesh):
    vertices = [
        [round(float(vertex.point[0]), 15), round(float(vertex.point[1]), 15)]
        for vertex in mesh.vertices
    ]
    elements = [
        {
            "material": element_material(mesh, element),
            "vertices": [int(vertex.nr) for vertex in element.vertices],
        }
        for element in mesh.Elements(VOL)
    ]
    payload = json.dumps(
        {"vertices": vertices, "elements": elements},
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def package_version(name):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def git_commit():
    override = os.environ.get("RADIA_SOURCE_COMMIT", "").strip()
    if len(override) in {40, 64}:
        return override
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "0" * 40


def save_and_reload_vol(mesh, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    mesh.ngmesh.Save(str(path))
    return Mesh(str(path))


def run(
    *,
    t_end=25.0,
    dt=25.0,
    em_billet_maxh=0.0015,
    thermal_maxh=0.002,
    em_air_maxh=0.035,
    output="",
    workdir="C:/temp/radia_ih_team36",
    reference=None,
    profile="smoke",
):
    started = time.perf_counter()
    timing = {"mesh": 0.0, "electromagnetic": 0.0, "mapping": 0.0, "thermal": 0.0}
    workdir_path = Path(workdir)
    workdir_path.mkdir(parents=True, exist_ok=True)

    mesh_started = time.perf_counter()
    em_mesh = build_em_mesh(maxh_billet=em_billet_maxh, maxh_air=em_air_maxh)
    th_mesh = build_thermal_mesh(maxh=thermal_maxh)
    em_mesh = save_and_reload_vol(em_mesh, workdir_path / "team36_em.vol")
    th_mesh = save_and_reload_vol(th_mesh, workdir_path / "team36_thermal.vol")
    timing["mesh"] = time.perf_counter() - mesh_started

    em_hash = mesh_topology_sha256(em_mesh)
    thermal_hash = mesh_topology_sha256(th_mesh)
    if em_hash == thermal_hash:
        raise RuntimeError("EM and thermal meshes must be topologically noncoincident")

    fes_t = H1(th_mesh, order=1)
    temperature = GridFunction(fes_t)
    temperature.vec[:] = 20.0
    print(
        f"EM mesh: {em_mesh.ne} triangles, {em_mesh.nv} vertices; "
        f"thermal mesh: {th_mesh.ne} triangles, {th_mesh.nv} vertices"
    )

    history = [
        {
            "time_s": 0.0,
            "axis_temperature_c": 20.0,
            "surface_temperature_c": 20.0,
            "maximum_temperature_c": 20.0,
            "induced_power_w": 0.0,
        }
    ]
    maximum_mapping_error = 0.0
    maximum_scale_deviation = 0.0
    t_map_sample_count = sum(
        element_material(em_mesh, element) == "billet"
        for element in em_mesh.Elements(VOL)
    )
    q_map_sample_count = int(th_mesh.ne)
    n_steps = int(round(float(t_end) / float(dt)))
    if not math.isclose(n_steps * float(dt), float(t_end), abs_tol=1.0e-12):
        raise ValueError("t_end must be an integer multiple of dt")

    for step in range(n_steps):
        map_started = time.perf_counter()
        t_map, t_values = map_temperature_to_em(temperature, th_mesh, em_mesh)
        timing["mapping"] += time.perf_counter() - map_started

        em_started = time.perf_counter()
        em = solve_em(em_mesh, t_map)
        timing["electromagnetic"] += time.perf_counter() - em_started

        map_started = time.perf_counter()
        q_th, raw_power, mapped_power, scale = map_power_to_thermal(
            em["q"], em_mesh, th_mesh, em["power_w"]
        )
        timing["mapping"] += time.perf_counter() - map_started
        mapping_error = abs(mapped_power - em["power_w"]) / max(
            abs(em["power_w"]), 1.0e-30
        )
        maximum_mapping_error = max(maximum_mapping_error, mapping_error)
        maximum_scale_deviation = max(maximum_scale_deviation, abs(scale - 1.0))

        thermal_started = time.perf_counter()
        thermal_iterations, thermal_rel = thermal_step(
            th_mesh, temperature, q_th, dt
        )
        timing["thermal"] += time.perf_counter() - thermal_started

        axis = eval_real(temperature, th_mesh, (1.0e-9, 1.0e-9))
        surface = eval_real(temperature, th_mesh, (0.03 - 1.0e-9, 1.0e-9))
        temperature_values = np.asarray(temperature.vec.FV().NumPy(), dtype=float)
        row = {
            "time_s": (step + 1) * float(dt),
            "axis_temperature_c": axis,
            "surface_temperature_c": surface,
            "maximum_temperature_c": float(np.max(temperature_values)),
            "minimum_temperature_c": float(np.min(temperature_values)),
            "induced_power_w": em["power_w"],
            "mapped_power_w": mapped_power,
            "raw_mapped_power_w": raw_power,
            "power_mapping_relative_error": mapping_error,
            "power_mapping_scale": scale,
            "em_iterations": em["iterations"],
            "em_converged": em["converged"],
            "em_relative_update": em["relative_update"],
            "thermal_iterations": thermal_iterations,
            "thermal_converged": thermal_rel < 2.0e-5,
            "thermal_relative_update": thermal_rel,
            "temperature_map_min_c": float(np.min(t_values)),
            "temperature_map_max_c": float(np.max(t_values)),
            "mu_r_min": em["mu_r_min"],
            "mu_r_max": em["mu_r_max"],
        }
        history.append(row)
        print(
            f"t={row['time_s']:6.1f} s  P={row['induced_power_w']:10.3f} W  "
            f"T_axis={axis:8.3f} C  T_surface={surface:8.3f} C  "
            f"map={mapping_error:.2e}"
        )

    time_label = f"{float(t_end):g}".replace(".", "p")
    temperature_file = workdir_path / f"team36_temperature_{time_label}s.sol"
    temperature.Save(str(temperature_file))

    geometry = dict(EXPECTED_GEOMETRY_M)
    excitation = dict(EXPECTED_EXCITATION)
    excitation["duration_s"] = float(t_end)
    tables = {
        "resistivity_temperature_c": RHOEL_T.tolist(),
        "resistivity_ohm_m": RHOEL.tolist(),
        "field_strength_a_m": H_TABLE.tolist(),
        "mu20": MU20_TABLE.tolist(),
        "conductivity_temperature_c": K_T.tolist(),
        "conductivity_w_m_k": K_TABLE.tolist(),
        "heat_capacity_temperature_c": CP_T.tolist(),
        "heat_capacity_j_kg_k": CP_TABLE.tolist(),
    }
    identity = {
        "geometry_sha256": identity_digest(geometry),
        "material_tables_sha256": identity_digest(tables),
        "excitation_sha256": identity_digest(excitation),
        "coordinate_system": "axisymmetric_r_z",
    }
    timing["total"] = time.perf_counter() - started
    executed_at_utc = datetime.now(timezone.utc).isoformat()
    host = socket.gethostname()
    radia_version = package_version("radia")
    ngsolve_version = package_version("ngsolve")
    artifact = {
        "radia_version": radia_version,
        "executed_at_utc": executed_at_utc,
        "host": host,
        "artifact_schema": SCHEMA,
        "benchmark_source": SOURCE_URL,
        "profile": profile,
        "coordinate_system": "axisymmetric_r_z",
        "identity": identity,
        "geometry_m": geometry,
        "excitation": excitation,
        "material_model": {
            "resistivity_point_count": len(RHOEL_T),
            "mu20_point_count": len(H_TABLE),
            "conductivity_point_count": len(K_T),
            "heat_capacity_point_count": len(CP_T),
            "curie_temperature_c": 770.0,
            "transition_width_c": 20.0,
            "density_kg_m3": 7800.0,
            "tables_sha256": identity["material_tables_sha256"],
        },
        "meshes": {
            "electromagnetic": {
                "file": "team36_em.vol",
                "topology_sha256": em_hash,
                "vertex_count": int(em_mesh.nv),
                "element_count": int(em_mesh.ne),
                "element_kinds": ["triangle"],
                "element_order": 1,
                "billet_skin_layer_count": 5,
            },
            "thermal": {
                "file": "team36_thermal.vol",
                "topology_sha256": thermal_hash,
                "vertex_count": int(th_mesh.nv),
                "element_count": int(th_mesh.ne),
                "element_kinds": ["triangle"],
                "element_order": 1,
                "independent_radial_partition_count": 5,
            },
        },
        "coupling": {
            "temperature_to_em": {
                "source_mesh_sha256": thermal_hash,
                "target_mesh_sha256": em_hash,
                "sample_count": int(t_map_sample_count),
                "outside_count": 0,
                "method": "element_centroid_interpolation",
            },
            "joule_power_to_thermal": {
                "source_mesh_sha256": em_hash,
                "target_mesh_sha256": thermal_hash,
                "sample_count": q_map_sample_count,
                "method": "centroid_interpolation_with_global_power_conservation",
                "maximum_relative_error": maximum_mapping_error,
                "maximum_scale_deviation": maximum_scale_deviation,
            },
        },
        "history": history,
        "provenance": {
            "executed_at_utc": executed_at_utc,
            "host": host,
            "radia_version": radia_version,
            "radia_mcp_version": package_version("radia-mcp"),
            "ngsolve_version": ngsolve_version,
            "git_commit": git_commit(),
            "result_file": Path(output).name if output else "",
            "temperature_file": temperature_file.name,
        },
        "timing_s": timing,
    }
    artifact["gate"] = evaluate_team36_artifact(artifact, reference=reference)

    if output:
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(artifact, indent=2, ensure_ascii=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(f"wrote {output_path}")
    print(json.dumps(artifact["gate"], indent=2))
    return artifact


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("smoke", "validation"), default="smoke")
    parser.add_argument("--mesh-profile", choices=("baseline", "refined"), default="baseline")
    parser.add_argument("--dt", type=float, default=25.0)
    parser.add_argument("--output", default="")
    parser.add_argument("--workdir", default="C:/temp/radia_ih_team36")
    parser.add_argument("--reference", default="")
    args = parser.parse_args()

    t_end = 25.0 if args.profile == "smoke" else 250.0
    mesh = {
        "baseline": (0.0015, 0.0020, 0.035),
        "refined": (0.0010, 0.0015, 0.030),
    }[args.mesh_profile]
    output = args.output or str(
        Path(__file__).with_name(f"results_{args.profile}_{args.mesh_profile}.json")
    )
    reference = None
    if args.reference:
        reference = json.loads(Path(args.reference).read_text(encoding="utf-8"))
    run(
        t_end=t_end,
        dt=args.dt,
        em_billet_maxh=mesh[0],
        thermal_maxh=mesh[1],
        em_air_maxh=mesh[2],
        output=output,
        workdir=args.workdir,
        reference=reference,
        profile=args.profile,
    )


if __name__ == "__main__":
    main()
