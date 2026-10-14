"""Common physical-volume energy quadrature for current-excited comparisons.

This is a diagnostic, not a solver objective. An explicit Kelvin map includes
the exterior; without it the result is a finite-domain partial integral.
Permanent magnets, hysteresis and anisotropy need different energy contracts.
"""
from __future__ import annotations

import hashlib
import json
import numpy as np

MU0 = 4e-7 * np.pi


def integrate_energy(B, H, weights, regions, laws):
    """Integrate independent W(B) and W'(H) with the same material law.

    ``laws`` maps material labels to relative permeability or [H,B] tables.
    A nonzero Fenchel gap measures constitutive inconsistency; we never force
    W+W'=H.B by constructing one observable from the other returned field.
    """
    from radia.scalar_potential_solver import _build_bh_coenergy_interpolator
    from radia.vector_potential_solver import _build_nu_of_b_interpolator

    B, H = np.asarray(B), np.asarray(H)
    w = np.asarray(weights, dtype=float)
    labels = np.asarray(regions)
    if (np.iscomplexobj(B) or np.iscomplexobj(H) or B.shape != H.shape
            or B.shape != (len(w), 3) or labels.shape != w.shape
            or not len(w) or not np.isfinite(B).all() or not np.isfinite(H).all()
            or not np.isfinite(w).all() or np.any(w <= 0)):
        raise ValueError("finite real B/H vectors and positive physical volume weights required")
    if set(labels) - set(laws):
        raise ValueError("missing material energy law")
    b, h = np.linalg.norm(B, axis=1), np.linalg.norm(H, axis=1)
    rows = {}
    for label in sorted(set(labels)):
        take = labels == label
        law = np.asarray(laws[label], dtype=float)
        if law.ndim == 0:
            if not np.isfinite(law) or law <= 0:
                raise ValueError("relative permeability must be positive")
            mu = MU0 * float(law)
            energy, coenergy = b[take]**2 / (2 * mu), mu * h[take]**2 / 2
        else:
            if (law.ndim != 2 or law.shape[1] != 2 or len(law) < 2
                    or not np.isfinite(law).all() or np.any(law[0] != 0)
                    or np.any(np.diff(law, axis=0) <= 0)):
                raise ValueError("single-valued strictly increasing [H,B] table from (0,0) required")
            primitive = _build_bh_coenergy_interpolator(law)
            reluctivity, _ = _build_nu_of_b_interpolator(law)
            hb = np.array([reluctivity(x) * x for x in b[take]])
            energy = hb * b[take] - np.array([primitive(x) for x in hb])
            coenergy = np.array([primitive(x) for x in h[take]])
        hb_density = np.einsum("ij,ij->i", H[take], B[take])
        rows[str(label)] = {
            "volume_m3": float(w[take].sum()),
            "energy_J": float(w[take] @ energy),
            "coenergy_J": float(w[take] @ coenergy),
            "H_dot_B_J": float(w[take] @ hb_density),
            "fenchel_gap_J": float(w[take] @ (energy + coenergy - hb_density)),
        }
    return {"regions": rows, **{
        key: sum(row[key] for row in rows.values())
        for key in ("volume_m3", "energy_J", "coenergy_J", "H_dot_B_J", "fenchel_gap_J")
    }}


def compare_energy(records):
    """Report scalar ordering without asserting a variational bound."""
    names = ("hdiv_mmm", "reduced_a", "mixed_total_reduced_omega")
    if set(records) != set(names):
        raise ValueError("three converged energy records required")
    contracts = [records[n]["contract"] for n in names]
    if any(c != contracts[0] for c in contracts[1:]):
        raise ValueError("energy quadrature/material/domain contracts differ")
    out = {"contract": contracts[0], "ordering_is_bound": False}
    for key in ("energy_J", "coenergy_J"):
        values = np.array([records[n][key] for n in names], dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("non-finite energy")
        scale = max(float(np.max(np.abs(values))), np.finfo(float).tiny)
        out[key] = {"values": dict(zip(names, values.tolist())),
                    "relative_spread": float(np.ptp(values) / scale),
                    "hdiv_between": bool(min(values[1:]) <= values[0] <= max(values[1:]))}
    return out


class PhysicalVolumeEnergy:
    """Common FEM quadrature with an optional physical exterior pullback."""
    def __init__(self, mesh, bh_table, order, *, iron="iron", exterior="kelvin",
                 kelvin_center=None, kelvin_radius=None, physical_center=(0., 0., 0.)):
        import ngsolve as ng
        from ngsolve.comp import IntegrationRuleSpace
        if type(order) is not int or order < 1:
            raise ValueError("positive integration order required")
        names = tuple(mesh.GetMaterials())
        if set(names) - {iron, exterior, "air"}:
            raise ValueError("energy observer supports only iron/air/Kelvin current-excited meshes")
        with ng.TaskManager():
            space = IntegrationRuleSpace(mesh, order=order)
            rule = ng.dx(intrules=space.GetIntegrationRules())
            weight = ng.LinearForm(space)
            weight += space.TestFunction() * rule
            weight.Assemble()
            coords = ng.GridFunction(ng.VectorValued(space, dim=3))
            coords.Interpolate(ng.CF((ng.x, ng.y, ng.z)))
            ids = ng.GridFunction(space)
            ids.Interpolate(mesh.MaterialCF({name: i for i, name in enumerate(names)}))
        labels = np.array(names)[np.rint(ids.vec.FV().NumPy()).astype(int)]
        include_exterior = kelvin_center is not None and kelvin_radius is not None
        if (kelvin_center is None) != (kelvin_radius is None):
            raise ValueError("Kelvin centre and radius must be supplied together")
        if include_exterior and exterior not in names:
            raise ValueError("explicit exterior map requires a Kelvin mesh region")
        keep = np.ones(len(labels), dtype=bool) if include_exterior else labels != exterior
        self.computational_points = coords.vec.FV().NumPy().reshape(3, space.ndof).T[keep].copy()
        self.points = self.computational_points.copy()
        self.weights = weight.vec.FV().NumPy()[keep].copy()
        self.regions = labels[keep]
        self.exterior = self.regions == exterior
        self.kelvin_center = kelvin_center
        self.kelvin_radius = kelvin_radius
        self.physical_center = physical_center
        if include_exterior:
            from radia.kelvin_source import kelvin_computational_to_physical
            if not np.isfinite(kelvin_radius) or kelvin_radius <= 0:
                raise ValueError("positive Kelvin radius required")
            delta = self.points[self.exterior] - np.asarray(kelvin_center)
            rho = np.linalg.norm(delta, axis=1)
            if np.any(rho <= 0) or np.any(rho > kelvin_radius*(1+1e-10)):
                raise ValueError(f"Kelvin quadrature radii [{rho.min()}, {rho.max()}] outside (0, {kelvin_radius}]")
            self.points[self.exterior] = kelvin_computational_to_physical(
                self.points[self.exterior], kelvin_center, physical_center, kelvin_radius)
            self.weights[self.exterior] *= (kelvin_radius / rho)**6
        self.iron = iron
        self.laws = {iron: np.asarray(bh_table).tolist(), "air": 1.0, exterior: 1.0}
        digest = hashlib.sha256()
        for array in (self.points, self.weights):
            digest.update(np.asarray(array, dtype="<f8").tobytes())
        digest.update(json.dumps(self.regions.tolist()).encode())
        self.contract = {
            "schema": "radia.validation.physical-volume-energy.v1",
            "scope": ("physical_and_kelvin_exterior" if np.any(labels != exterior) else "kelvin_exterior_only")
                     if include_exterior else "physical_mesh_volume_only",
            "all_space_energy": bool(include_exterior and np.any(labels != exterior)),
            "quadrature_convergence_verified": False,
            "excluded": ([] if np.any(labels != exterior) else ["physical_interior"])
                        if include_exterior else [exterior, "exterior_tail"],
            "kelvin_center_m": list(kelvin_center) if include_exterior else None,
            "physical_center_m": list(physical_center), "kelvin_radius_m": kelvin_radius,
            "quadrature_order": order,
            "quadrature_sha256": digest.hexdigest(), "points": len(self.points),
            "laws": self.laws, "law_interpolation": "pchip_with_vacuum_slope_tail",
        }

    def _record(self, B, H):
        return {"contract": self.contract, **integrate_energy(
            B, H, self.weights, self.regions, self.laws)}

    def fem(self, mesh, B, H):
        if np.any(self.exterior):
            raise ValueError("choose the formulation-specific exterior field adapter")
        mapped = mesh(*self.points.T)
        return self._record(np.asarray(B(mapped)), np.asarray(H(mapped)))

    def _exterior_form(self, values, form):
        from radia.kelvin_source import kelvin_solution_to_physical
        return kelvin_solution_to_physical(
            values, self.points[self.exterior], kelvin_center=self.kelvin_center,
            physical_center=self.physical_center, radius=self.kelvin_radius, form=form)

    def omega(self, mesh, B, H):
        mapped = mesh(*self.computational_points.T)
        b, h = np.asarray(B(mapped)).copy(), np.asarray(H(mapped)).copy()
        if np.any(self.exterior):
            b[self.exterior] = self._exterior_form(b[self.exterior], "flux_density")
            h[self.exterior] = self._exterior_form(h[self.exterior], "field_strength")
        return self._record(b, h)

    def reduced_a(self, solver, coil):
        import ngsolve as ng
        import radia as rad
        mapped = solver.mesh(*self.computational_points.T)
        b = np.asarray(solver.get_B()(mapped)).copy()
        h = np.asarray(solver.get_H()(mapped)).copy()
        if np.any(self.exterior):
            # Only the reaction curl is in the Kelvin frame. Add the actual
            # source at physical exterior points after transforming it.
            reaction = np.asarray(ng.curl(solver.get_A())(mapped))[self.exterior]
            b[self.exterior] = self._exterior_form(reaction, "flux_density")
            b[self.exterior] += np.asarray(rad.Fld(coil, "b", self.points[self.exterior]))
            h[self.exterior] = b[self.exterior] / MU0
        return self._record(b, h)

    def hdiv(self, result, coil):
        import ngsolve as ng
        import radia as rad
        from radia import vim
        if result.get("image") is not None:
            raise ValueError("energy observation of image-reduced iron needs explicit material mapping")
        with ng.TaskManager():
            H = np.asarray(vim.FieldFromSolution(result, self.points, algorithm="direct"))
        H += np.asarray(rad.Fld(coil, "h", self.points))
        B = MU0 * H
        selected = self.regions == self.iron
        gf = result["gfM"]
        B[selected] += MU0 * np.asarray(gf(gf.space.mesh(*self.points[selected].T)))
        return self._record(B, H)
