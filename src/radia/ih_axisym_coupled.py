"""Coupled axisymmetric induction heating: volumetric eddy currents + heat.

For a body of revolution heated by coaxial coils the magnetic vector
potential has a single component ``A_phi(r, z)``.  This module solves the
time-harmonic eddy-current problem on an (r, z) mesh of workpiece, coils and
air with element-wise ``sigma(T)`` and ``mu_r(T)``, deposits the
time-averaged Joule density ``Q = sigma |omega A|^2 / 2`` in the workpiece
volume, and advances the temperature with the enthalpy Newton integrator
(:class:`radia.ih_heat_transient.NonlinearHeatStepper`).  The EM problem is
re-solved from the current temperature every ``em_every`` time steps
(staggered coupling; ``em_every = 1`` is first order in ``dt``).

Unlike a surface-impedance model this resolves the skin depth, so it stays
valid when the skin depth grows to millimetres above the Curie band -- the
regime where a frozen surface source was shown to fail
(validation_test/induction_heating/results/coupled_curie_cylinder_frozen_ht
.json).

Weak form (A = A_phi, test v, both in H1 with A = 0 on the axis)::

    int nu [ (r dA/dr + A)(r dv/dr + v) / r + r dA/dz dv/dz ] dr dz
      + j omega int_wp sigma r A v dr dz  =  int_coil J r v dr dz

with ``nu = 1 / (mu0 mu_r)``.  The natural boundary condition is
``H_tangential = 0``; Dirichlet boundaries are chosen by name.

Every EM solve is checked twice.  The discrete identity
``P_joule = P_coil = Re(j omega int_coil J* A 2 pi r) / 2`` holds for any
Galerkin solution, so it checks the linear solve and the quadrature, not
accuracy.  Accuracy is gated by the skin-depth resolution: every workpiece
element within three local skin depths of the workpiece boundary must be
no larger than ``skin_resolution`` times the local skin depth (default 1).
Against the Bessel solution of an infinite cylinder, order 2 at a 0.25 mm
skin depth gives +5.8 % at h/delta = 2, +0.29 % at 1 and +1e-4 at 0.5.

The permeability is linear per element (no saturation).  Material data come
from the user; nothing here is alloy data.
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field

import numpy as np

MU0 = 4.0e-7 * math.pi


# ---------------------------------------------------------------------------
# EM material
# ---------------------------------------------------------------------------

@dataclass
class EMMaterialTable:
    """Electrical conductivity and linear relative permeability vs T."""
    T: np.ndarray
    sigma: np.ndarray
    mu_r: np.ndarray
    allow_extrapolation: bool = False
    source: str = "table"
    max_excursion_C: float = 0.0

    @classmethod
    def constant(cls, sigma, mu_r):
        return cls(T=np.array([-273.15, 1.0e4]), sigma=np.array([sigma] * 2),
                   mu_r=np.array([mu_r] * 2), source="constant")

    @classmethod
    def from_csv(cls, path, *, allow_extrapolation=False):
        with open(path, encoding="utf-8-sig", newline="") as handle:
            rows = [r for r in csv.reader(handle)
                    if r and not r[0].lstrip().startswith("#")]
        header = [h.strip().lower() for h in rows[0]]
        want = ("t_c", "sigma_s_m", "mu_r")
        if tuple(header[:3]) != want:
            raise ValueError(f"{path}: header must start with "
                             f"{','.join(want)}, got {rows[0]}")
        d = np.asarray([[float(v) for v in r[:3]] for r in rows[1:]])
        return cls(T=d[:, 0], sigma=d[:, 1], mu_r=d[:, 2],
                   allow_extrapolation=allow_extrapolation, source=str(path))

    def __post_init__(self):
        self.T = np.asarray(self.T, float)
        self.sigma = np.asarray(self.sigma, float)
        self.mu_r = np.asarray(self.mu_r, float)
        if not (len(self.T) == len(self.sigma) == len(self.mu_r) >= 2):
            raise ValueError("EM material table needs at least two rows")
        if not np.all(np.diff(self.T) > 0):
            raise ValueError("EM material temperatures must increase")
        if np.any(self.sigma <= 0) or np.any(self.mu_r < 1.0):
            raise ValueError("need sigma > 0 and mu_r >= 1")

    def _T(self, T):
        T = np.asarray(T, float)
        lo, hi = self.T[0], self.T[-1]
        over = max(0.0, float(T.max()) - hi, lo - float(T.min()))
        if over > 0:
            if not self.allow_extrapolation:
                raise ValueError(
                    f"workpiece temperature {T.min():.1f}..{T.max():.1f} C "
                    f"leaves the EM material table {lo:.1f}..{hi:.1f} C; "
                    "extend it or allow extrapolation explicitly")
            self.max_excursion_C = max(self.max_excursion_C, over)
        return np.clip(T, lo, hi)

    def evaluate(self, T):
        Tc = self._T(T)
        return (np.interp(Tc, self.T, self.sigma),
                np.interp(Tc, self.T, self.mu_r))

    def audit(self):
        return {"source": self.source,
                "T_range_C": [float(self.T[0]), float(self.T[-1])],
                "sigma_range_S_m": [float(self.sigma.min()),
                                    float(self.sigma.max())],
                "mu_r_range": [float(self.mu_r.min()), float(self.mu_r.max())],
                "allow_extrapolation": self.allow_extrapolation,
                "max_excursion_C": self.max_excursion_C}


# ---------------------------------------------------------------------------
# Eddy-current solver
# ---------------------------------------------------------------------------

def _check_axis_dirichlet(mesh, pts, dirichlet):
    """Every boundary segment on r = 0 must be in the Dirichlet set."""
    import re
    from ngsolve import BND
    if not dirichlet:
        raise ValueError("a Dirichlet set containing the r = 0 axis is "
                         "required (A_phi = 0 there)")
    pat = re.compile(str(dirichlet))
    scale = max(float(np.max(np.abs(pts[:, 0]))), 1e-300)
    on_axis, missing = 0, set()
    for el in mesh.Elements(BND):
        r = [abs(float(pts[v.nr, 0])) for v in el.vertices]
        if max(r) <= 1e-12 * scale:
            on_axis += 1
            if not pat.fullmatch(el.mat):
                missing.add(el.mat)
    if on_axis == 0:
        raise ValueError("the mesh has no boundary on the r = 0 axis")
    if missing:
        raise ValueError(
            f"boundary segments on the r = 0 axis named {sorted(missing)} are "
            f"not in the Dirichlet set {dirichlet!r}; A_phi must vanish on the "
            "whole axis")


class AxisymEddyCurrent:
    """Time-harmonic A_phi on an (r, z) mesh with element-wise materials.

    ``coils`` maps a mesh material name to the total peak current [A]
    (ampere-turns) through that region; the current density is uniform
    over the region.
    """

    def __init__(self, mesh, *, frequency, workpiece, coils, dirichlet,
                 order=2, linear_solver="pardiso"):
        from ngsolve import CF, GridFunction, H1, Integrate, L2
        mats = set(mesh.GetMaterials())
        if workpiece not in mats:
            raise ValueError(f"workpiece material {workpiece!r} is not in "
                             f"the mesh ({sorted(mats)})")
        for name in coils:
            if name not in mats:
                raise ValueError(f"coil material {name!r} is not in the mesh")
        pts = np.asarray([v.point for v in mesh.vertices])
        if np.min(pts[:, 0]) < -1e-12:
            raise ValueError("an axisymmetric mesh needs r = x >= 0")
        _check_axis_dirichlet(mesh, pts, dirichlet)
        self.mesh = mesh
        self.workpiece = workpiece
        self.omega = 2.0 * math.pi * float(frequency)
        self.frequency = float(frequency)
        self.fes = H1(mesh, order=int(order), complex=True,
                      dirichlet=str(dirichlet))
        self.p0 = L2(mesh, order=0)
        self.gf_sigma = GridFunction(self.p0)
        self.gf_mu = GridFunction(self.p0)
        self.gf_mu.vec[:] = 1.0
        self.gfA = GridFunction(self.fes)
        self.linear_solver = linear_solver
        area = np.asarray(Integrate(CF(1.0), mesh, element_wise=True), float)
        self._wp = np.asarray([el.nr for el in mesh.Elements()
                               if el.mat == workpiece])
        self._area = area
        dens = {}
        for name, current in coils.items():
            a = float(Integrate(CF(1.0), mesh,
                                definedon=mesh.Materials(name)).real)
            dens[name] = float(current) / a
        self.J = mesh.MaterialCF(dens, default=0.0)
        self.coils = dict(coils)

    @property
    def workpiece_elements(self):
        return self._wp

    def solve(self, sigma_wp, mu_wp):
        """Solve with per-workpiece-element sigma, mu_r; return a record."""
        from ngsolve import (BilinearForm, Integrate, InnerProduct,
                             LinearForm, dx, grad, x as r)
        s = np.zeros(self.mesh.ne)
        m = np.ones(self.mesh.ne)
        s[self._wp] = sigma_wp
        m[self._wp] = mu_wp
        self.gf_sigma.vec.FV().NumPy()[:] = s
        self.gf_mu.vec.FV().NumPy()[:] = m
        from ngsolve import ET, IntegrationRule
        nu = 1.0 / (MU0 * self.gf_mu)
        u, v = self.fes.TnT()
        # The Joule and source terms are polynomial (P0 material, r times
        # two order-p functions): integrate them exactly, so the discrete
        # power identity P_joule = P_coil holds to solver precision.
        p = self.fes.globalorder
        exact = {et: IntegrationRule(et, 2 * p + 1)
                 for et in {el.type for el in self.mesh.Elements()}}
        # The stiffness carries 1/r: its Gauss rule has interior points
        # only, so no quadrature point sits on the r = 0 axis.
        interior = {et: IntegrationRule(et, 2 * p + 2)
                    for et in {el.type for el in self.mesh.Elements()}}
        a = BilinearForm(self.fes, symmetric=True)
        a += nu / r * (r * grad(u)[0] + u) * (r * grad(v)[0] + v) \
            * dx(intrules=interior)
        a += nu * r * grad(u)[1] * grad(v)[1] * dx(intrules=interior)
        a += 1j * self.omega * self.gf_sigma * r * u * v \
            * dx(definedon=self.mesh.Materials(self.workpiece),
                 intrules=exact)
        f = LinearForm(self.fes)
        f += self.J * r * v * dx(intrules=exact)
        # Parallelism follows the caller's ngsolve.TaskManager().
        a.Assemble()
        f.Assemble()
        self.gfA.vec.data = a.mat.Inverse(
            self.fes.FreeDofs(), inverse=self.linear_solver) * f.vec
        if not np.all(np.isfinite(np.asarray(self.gfA.vec.FV().NumPy()))):
            raise RuntimeError("the eddy-current solution is not finite")
        wp = self.mesh.Materials(self.workpiece)
        P_joule = float(Integrate(
            self.heat_density() * 2 * math.pi * r, self.mesh,
            definedon=wp, order=2 * self.fes.globalorder + 3).real)
        S = Integrate(1j * self.omega * self.gfA * self.J * r, self.mesh,
                      order=2 * self.fes.globalorder + 3)
        P_coil = float(0.5 * (2 * math.pi) * complex(S).real)
        return {"P_joule_W": P_joule, "P_coil_W": P_coil,
                "power_balance_relative_error":
                    (P_joule - P_coil) / max(abs(P_coil), 1e-300)}

    def heat_density(self):
        """Time-averaged Joule density [W/m^3] as a coefficient function."""
        from ngsolve import InnerProduct
        return (0.5 * self.gf_sigma * self.omega ** 2
                * InnerProduct(self.gfA, self.gfA).real)

    def skin_depth(self, sigma, mu_r):
        return np.sqrt(2.0 / (self.omega * MU0 * np.asarray(mu_r)
                              * np.asarray(sigma)))

    def skin_resolution(self, delta_wp):
        """Largest h/delta over workpiece elements within three local skin
        depths of the workpiece boundary; returns (ratio, element nr)."""
        from scipy.spatial import cKDTree
        if not hasattr(self, "_geom"):
            pts = np.asarray([v.point for v in self.mesh.vertices])
            wp_v, cent, size = set(), [], []
            for el in self.mesh.Elements():
                if el.mat != self.workpiece:
                    continue
                vv = [v.nr for v in el.vertices]
                wp_v.update(vv)
                P = pts[vv]
                cent.append(P.mean(axis=0))
                size.append(max(np.linalg.norm(P[i] - P[j])
                                for i in range(len(vv))
                                for j in range(i + 1, len(vv))))
            # The eddy currents flow under the workpiece surface that faces
            # the field: its interface with non-conducting regions.  Edges
            # on the mesh boundary (the axis, symmetry or truncation cuts)
            # are not such a surface.
            wp_edges, other_edges = set(), set()
            for el in self.mesh.Elements():
                vv = [v.nr for v in el.vertices]
                es = {tuple(sorted((vv[i], vv[(i + 1) % len(vv)])))
                      for i in range(len(vv))}
                (wp_edges if el.mat == self.workpiece else other_edges) \
                    .update(es)
            interface = wp_edges & other_edges
            if not interface:
                raise ValueError("the workpiece has no interface with a "
                                 "non-conducting region")
            bpts = [0.5 * (pts[a] + pts[b]) for a, b in interface]
            self._geom = (np.asarray(cent), np.asarray(size),
                          cKDTree(np.asarray(bpts)))
        cent, size, tree = self._geom
        depth, _ = tree.query(cent)
        band = depth <= 3.0 * delta_wp
        if not np.any(band):
            return 0.0, -1
        ratio = np.where(band, size / delta_wp, 0.0)
        i = int(np.argmax(ratio))
        return float(ratio[i]), int(self._wp[i])


# ---------------------------------------------------------------------------
# Coupled run
# ---------------------------------------------------------------------------

def element_mean_temperature(mesh, gfT, elements):
    """Revolved-volume mean of gfT over each listed element."""
    from ngsolve import Integrate, x as r
    num = np.asarray(Integrate(gfT * r, mesh, element_wise=True), float)
    den = np.asarray(Integrate(r, mesh, element_wise=True), float)
    return num[elements] / den[elements]


def run_coupled(mesh, *, frequency, workpiece, coils, dirichlet,
                em_material, thermal_material, boundaries, dt, t_end,
                t_initial=20.0, em_every=1, em_order=2, thermal_order=2,
                power_tolerance=1e-7, newton_tol_K=1e-3, max_newton=25,
                max_halvings=0, skin_resolution=1.0, energy_tolerance=1e-4,
                on_step=None, linear_solver_em="pardiso",
                linear_solver_heat="sparsecholesky"):
    """Staggered EM-thermal transient.  Returns ``(result, gfT, em)``.

    ``em_every`` (>= 1) re-solves the EM problem every that many steps; the
    validation record shows +3.8 % (2) and +12 % (5) in the peak temperature
    across the Curie band against every step, so values above 1 trade
    accuracy for speed and are recorded.
    """
    from ngsolve import CF, GridFunction, H1, x as r
    from . import ih_heat_transient as iht

    if isinstance(em_every, bool) or int(em_every) != em_every or em_every < 1:
        raise ValueError(f"em_every must be an integer >= 1, got {em_every!r}")
    if not skin_resolution > 0:
        raise ValueError("skin_resolution must be positive")

    em = AxisymEddyCurrent(mesh, frequency=frequency, workpiece=workpiece,
                           coils=coils, dirichlet=dirichlet, order=em_order,
                           linear_solver=linear_solver_em)
    fesT = H1(mesh, order=int(thermal_order), definedon=workpiece)
    gfT = GridFunction(fesT)
    gfT.Set(CF(float(t_initial)), definedon=mesh.Materials(workpiece))
    wp_el = em.workpiece_elements
    Q = em.heat_density()
    stepper = iht.NonlinearHeatStepper(
        gfT, thermal_material, boundaries, weight=2 * math.pi * r,
        volume_source=Q, region=workpiece,
        linear_solver=linear_solver_heat, newton_tol_K=newton_tol_K,
        max_newton=max_newton, max_halvings=max_halvings)

    def em_update():
        # the table must cover every temperature in the elements, not only
        # their means
        em_material.evaluate(np.asarray(stepper._sampled_range()))
        Tel = element_mean_temperature(mesh, gfT, wp_el)
        sig, mu = em_material.evaluate(Tel)
        delta = em.skin_depth(sig, mu)
        ratio, worst = em.skin_resolution(delta)
        if ratio > skin_resolution:
            raise ValueError(
                f"the workpiece mesh does not resolve the skin depth: element "
                f"{worst} is {ratio:.2f} skin depths across (limit "
                f"{skin_resolution:g}; min skin depth "
                f"{float(delta.min()) * 1e3:.3f} mm). Refine the surface "
                "layer of the workpiece mesh.")
        rec = em.solve(sig, mu)
        if not abs(rec["power_balance_relative_error"]) <= power_tolerance:
            raise RuntimeError(
                "eddy-current solve inconsistent: Joule "
                f"{rec['P_joule_W']:.6e} W vs coil {rec['P_coil_W']:.6e} W")
        rec.update({"T_elem_max_C": float(Tel.max()),
                    "mu_r_min": float(mu.min()),
                    "skin_depth_min_m": float(delta.min()),
                    "skin_depth_max_m": float(delta.max()),
                    "skin_resolution_h_over_delta": ratio})
        return rec

    n = int(round(float(t_end) / float(dt)))
    if abs(n * float(dt) - float(t_end)) > 1e-9 * max(1.0, float(t_end)):
        raise ValueError("t_end must be a multiple of dt")
    history = []
    em_history = []
    rec = em_update()
    em_history.append({"t_s": 0.0, **rec})
    for step in range(1, n + 1):
        if step > 1 and (step - 1) % int(em_every) == 0:
            rec = em_update()
            em_history.append({"t_s": (step - 1) * float(dt), **rec})
        stepper.advance(float(dt))
        entry = {"t_s": step * float(dt), "P_W": rec["P_joule_W"],
                 "T_max_nodal_C": float(np.max(stepper._nodal()))}
        history.append(entry)
        if on_step is not None:
            on_step(step, entry)
    stepper.check_energy(energy_tolerance)
    result = {
        "model": "axisymmetric volumetric eddy current + enthalpy heat, "
                 "staggered",
        "skin_resolution_limit": float(skin_resolution),
        "frequency_Hz": float(frequency), "coils_A": dict(coils),
        "em_every": int(em_every), "dt_s": float(dt), "t_end_s": float(t_end),
        "em_order": int(em_order), "thermal_order": int(thermal_order),
        "em_material": em_material.audit(),
        "thermal_material": thermal_material.audit(),
        "history": history, "em_history": em_history,
        "nonlinear_transient": stepper.audit.as_dict(),
    }
    return result, gfT, em
