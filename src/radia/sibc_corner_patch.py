"""Matched-patch correction for a surface impedance near tips and corners.

A Leontovich impedance is the leading term of an expansion in the ratio of the
skin depth to the radius of curvature of the surface.  Where that ratio stops
being small -- a rounded tip only a skin depth or two across, or a genuine
corner, where the radius is zero and no term of the expansion applies --
refining the surface discretisation does not help: the error is in the model.
Measured on the beak fin at 150 kHz, a 2-D SIBC solve converged to 1024
perimeter panels loses ``17%`` of the section's dissipation outright and is
``+16%`` on the beak loss share and ``-15%`` on the tip share.  That section's
tip radius is 0.25 mm against a skin depth of 0.171 mm, and its beak root and
far end are square corners.

This implements the overlapping-patch idea of Proekt, Yuferev, Tsukerman and
Ida (2002): keep the global surface-impedance solution where it is valid,
resolve the conductor's interior in a neighbourhood of the offending features,
and match the two where both hold.

It is a matched decomposition rather than a literal overlap.  The patch and
the outer band partition the conductor at a chord, and the matching is one
way: the patch is told what the outer solution says on that chord, and the
outer solution is never told what the patch found.  That is a real limitation,
and it is what the residual is made of, so the chord is placed deep in the
valid region and the driver sweeps it instead of asserting a position.

A section usually has more than one offending feature.  Patching the
interesting one and leaving the rest to the outer model recovers the patched
region and leaves the rest wrong, which then contaminates every ratio taken
against the whole: on this fixture one patch fixes the fin's own dissipation
from ``-21%`` to ``+2%`` while the unpatched far end stays ``-13%`` out.
``side`` exists so both ends can be patched and the outer model keeps only the
flat band it is actually valid on.

Where the patch is cut matters more than how the matching is written.  The cut
is placed across the *thick body*, where the skin depth is a small fraction of
the local thickness and of every radius of curvature, so the outer solution is
accurate there in the sense its own expansion claims.  On that chord the outer
solution is not merely a surface current: to leading order its interior field
is the flat half-space profile ``exp(-(1+j) d / delta)`` measured from each
face, and that is what the patch is given.  A patch cut through the beak
instead would be matched where the outer model is already failing, which is the
one place the method must not be used.

    laplace(A) - j omega mu0 sigma A = -mu0 sigma E0   in the conductor
    laplace(A)                       =  0              in the air

``E0`` is the uniform axial electric field the global solve already produces as
its Lagrange multiplier, so the patch is driven by the same source as the
global problem rather than by an independent excitation.  On the surface the
patch's own current ``J = sigma (E0 - j omega A)`` reduces to the global
model's ``(E0 - Zs K) / Zs`` wherever the thin-skin limit holds, which is what
makes the two comparable on the overlap band rather than merely adjacent.

Dauge, Dular, Krahenbuhl, Peron, Perrussel and Poignard (2014) put the corner
layer of the eddy-current model at the scale of ``delta``; the patch has to
contain that layer, which it does by construction here, since it contains the
whole fin and its root corners.  What the cut position then tests is the
opposite question -- how far back the matching must sit before the answer stops
depending on it -- and that is measured rather than asserted.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MU0 = 4e-7 * np.pi


def surface_potential(panel_current, panel_ds, *, surface_impedance,
                      axial_E_field, omega):
    """``A_z`` at each panel midpoint, from the global solve's own equation.

    The surface-impedance panel equation reads ``Zs K + j omega A = E0``, so
    the global solution already carries its own surface potential exactly; it
    does not have to be recovered by re-integrating the singular kernel over
    the sheet the evaluation point lies on.
    """
    current = np.asarray(panel_current)
    ds = np.asarray(panel_ds, dtype=float)
    k_surface = current / ds
    return (complex(axial_E_field)
            - complex(surface_impedance) * k_surface) / (1j * omega)


def _balanced_sum(terms):
    """Sum CoefficientFunctions pairwise, so the tree stays shallow.

    A left fold over a thousand panels builds a tree a thousand deep, which is
    slow to construct and needlessly deep to evaluate.
    """
    items = list(terms)
    if not items:
        raise ValueError("nothing to sum")
    while len(items) > 1:
        items = [items[i] + items[i + 1] if i + 1 < len(items) else items[i]
                 for i in range(0, len(items), 2)]
    return items[0]


def exterior_potential_cf(panel_xy, panel_ds, panel_current):
    """``A_z`` in the air, as a CoefficientFunction over the panel currents.

    This is the same logarithmic kernel the global solve is built on, so the
    boundary data handed to the patch is the global solution itself rather than
    an independent approximation of it.  The distance is floored at half a
    panel width, which is the usual regularisation and keeps a boundary point
    that grazes the sheet bounded; :func:`check_exterior_potential` measures
    what that costs on the surface, where it is worst.
    """
    import ngsolve as ng

    xy = np.asarray(panel_xy, dtype=float)
    ds = np.asarray(panel_ds, dtype=float)
    current = np.asarray(panel_current)
    terms = []
    for (px, py), dsk, ik in zip(xy, ds, current):
        r2 = (ng.x - float(px)) ** 2 + (ng.y - float(py)) ** 2
        floor = (0.5 * float(dsk)) ** 2
        terms.append(ng.CF(complex(ik))
                     * ng.log(1.0 / ng.sqrt(ng.IfPos(r2 - floor, r2, floor))))
    return (MU0 / (2.0 * np.pi)) * _balanced_sum(terms)


def check_exterior_potential(panel_xy, panel_ds, panel_current, *,
                             surface_impedance, axial_E_field, omega,
                             sample=64):
    """Worst relative disagreement of the kernel CF with the exact surface A.

    The kernel is evaluated where it is least accurate -- on the sheet itself,
    at panel midpoints -- and compared with the potential the global equation
    gives there directly.  Anywhere else on the patch boundary the agreement is
    better, so this bounds the error in the boundary data.
    """
    xy = np.asarray(panel_xy, dtype=float)
    ds = np.asarray(panel_ds, dtype=float)
    current = np.asarray(panel_current)
    exact = surface_potential(current, ds, surface_impedance=surface_impedance,
                              axial_E_field=axial_E_field, omega=omega)
    take = np.linspace(0, len(xy) - 1, int(sample)).astype(int)
    offset = xy[take][:, None, :] - xy[None, :, :]
    distance = np.linalg.norm(offset, axis=2)
    floor = 0.5 * ds[None, :]
    distance = np.maximum(distance, floor)
    approx = (MU0 / (2.0 * np.pi)) * (np.log(1.0 / distance) @ current)
    scale = np.max(np.abs(exact))
    return float(np.max(np.abs(approx - exact[take])) / scale)


@dataclass(frozen=True)
class PatchSolution:
    """One resolved neighbourhood of the conductor, and how it was matched."""

    x_cut_m: float | None
    side: str
    skin_depth_m: float
    ndof: int
    mesh: object
    gf_A: object
    conductor: str
    sigma_S_per_m: float
    omega_rad_per_s: float
    axial_E_field_V_per_m: complex
    boundary_data_error: float | None

    def current_density(self):
        """``J = sigma (E0 - j omega A)`` inside the conductor."""
        import ngsolve as ng

        return self.sigma_S_per_m * (
            ng.CF(self.axial_E_field_V_per_m)
            - 1j * self.omega_rad_per_s * self.gf_A)

    def loss_density(self):
        """``|J|^2 / (2 sigma)``: time-average watts per cubic metre."""
        import ngsolve as ng

        j = self.current_density()
        return (j * ng.Conj(j)).real / (2.0 * self.sigma_S_per_m)

    def loss_beyond(self, x_m):
        """Dissipation per unit length in the conductor beyond ``x_m``."""
        import ngsolve as ng

        return float(ng.Integrate(
            self.loss_density() * ng.IfPos(ng.x - float(x_m), 1.0, 0.0),
            self.mesh, definedon=self.mesh.Materials(self.conductor),
            order=10).real)

    def total_loss(self):
        """Dissipation per unit length in the whole patched conductor."""
        import ngsolve as ng

        return float(ng.Integrate(
            self.loss_density(), self.mesh,
            definedon=self.mesh.Materials(self.conductor), order=10).real)

    def total_current(self):
        """Axial current carried by the patched conductor, in amperes."""
        import ngsolve as ng

        return complex(ng.Integrate(
            self.current_density(), self.mesh,
            definedon=self.mesh.Materials(self.conductor), order=10))


def solve_cut_patch(section_face, *, x_cut, panel_xy, panel_ds, panel_current,
                    axial_E_field, surface_impedance, omega, sigma,
                    maxh_conductor, maxh_air=None, air_pad=None, order=3,
                    curve_order=3, conductor="conductor",
                    exterior_potential=None, side="high"):
    """Resolve the conductor's interior beyond ``x_cut``, matched to the outer.

    Args:
        section_face: ``netgen.occ`` face of the conductor cross-section, in
            metres, with the fin pointing along +x.
        x_cut: matching plane, in metres.  It must lie in the thick body, where
            the outer solution is valid; the caller is expected to sweep it and
            show the answer does not depend on where it sits.  ``None``
            resolves the whole section instead, which is the control the sweep
            is measured against rather than a patch.
        panel_xy, panel_ds, panel_current: the global surface-impedance solve.
        axial_E_field, surface_impedance: that solve's multiplier and ``Zs``.
        omega, sigma: angular frequency and conductivity.
        maxh_conductor: element size in the conductor; it has to resolve the
            skin, not the geometry.
        maxh_air: element size in the air (default: twenty times the skin).
        air_pad: air margin beyond the section (default: ten skin depths).
        order, curve_order: polynomial and geometry order.
        exterior_potential: optional CoefficientFunction replacing the panel
            kernel on the air boundary, for a geometry whose exterior field is
            known in closed form. Its accuracy cannot be inferred from the
            panel data: a warning is emitted and boundary_data_error is None.
        side: ``"high"`` keeps the conductor beyond ``x_cut``, ``"low"`` keeps
            the part before it.  A section with a failing feature at each end
            needs one of each, with the outer model left the band between.

    Returns:
        PatchSolution.
    """
    import ngsolve as ng
    from netgen.occ import Glue, MoveTo, OCCGeometry, X, Y

    delta = float(np.sqrt(2.0 / (omega * MU0 * sigma)))
    if air_pad is None:
        air_pad = 10.0 * delta
    if maxh_air is None:
        maxh_air = 20.0 * delta

    # Geometry comes from the panel polygon, not from the CAD bounding box:
    # OCC pads a bounding box by a gap that is a large fraction of a skin
    # depth at this scale, and the matching data has to sit on the surface the
    # global solve actually discretised.
    xy = np.asarray(panel_xy, dtype=float)
    y_lo, y_hi = float(np.min(xy[:, 1])), float(np.max(xy[:, 1]))
    x_lo, x_hi = float(np.min(xy[:, 0])), float(np.max(xy[:, 0]))
    thickness = y_hi - y_lo

    if side not in ("high", "low"):
        raise ValueError("side must be 'high' (keep x > x_cut) or 'low'")
    keep_high = side == "high"

    if x_cut is None:
        # No cut: the whole section is resolved, and the only approximation
        # left is the truncation of the open air region.  This is the control
        # the cut sweep is measured against, not a patch.
        solid = section_face
        x_left, x_right = x_lo - air_pad, x_hi + air_pad
    else:
        if not x_lo < x_cut < x_hi:
            raise ValueError("x_cut must fall inside the section")
        # The cut has to be a full-height chord of the body: that is what
        # makes the half-space profile below the right matching data.  A cut
        # through the taper would be matched where the outer model is already
        # failing.
        probe = _chord_extent(xy, x_cut)
        slack = 1e-3 * thickness
        if abs(probe[0] - y_lo) > slack or abs(probe[1] - y_hi) > slack:
            raise ValueError(
                f"the cut at x={x_cut:.4e} m spans "
                f"y={probe[0]:.4e}..{probe[1]:.4e} but the section spans "
                f"{y_lo:.4e}..{y_hi:.4e}: place the cut in the "
                "full-thickness body, not in the fin")
        if thickness < 6.0 * delta:
            raise ValueError(
                f"the body is {thickness / delta:.1f} skin depths thick at "
                "the cut; the outer solution is not valid there")
        # Trim in x only: the window clears the section in y, so the boolean
        # has no edge coincident with the conductor's own faces.
        x_left = x_cut if keep_high else x_lo - air_pad
        x_right = x_hi + air_pad if keep_high else x_cut
        keep = MoveTo(x_left, y_lo - thickness).Rectangle(
            x_right - x_left, 3.0 * thickness).Face()
        solid = section_face * keep

    solid.faces.name = conductor
    solid.maxh = float(maxh_conductor)
    if x_cut is not None:
        (solid.edges.Min(X) if keep_high else solid.edges.Max(X)).name = "cut"

    air = MoveTo(x_left, y_lo - air_pad).Rectangle(
        x_right - x_left, thickness + 2.0 * air_pad).Face()
    air.faces.name = "air"
    air.maxh = float(maxh_air)
    outside = air - solid
    outside.faces.name = "air"
    for selector in (outside.edges.Min(X), outside.edges.Max(X),
                     outside.edges.Min(Y), outside.edges.Max(Y)):
        selector.name = "outer"

    shape = Glue([outside, solid])
    mesh = ng.Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=float(maxh_air)))
    if curve_order > 1:
        mesh.Curve(int(curve_order))

    if exterior_potential is None:
        outer_cf = exterior_potential_cf(panel_xy, panel_ds, panel_current)
        data_error = check_exterior_potential(
            panel_xy, panel_ds, panel_current,
            surface_impedance=surface_impedance,
            axial_E_field=axial_E_field, omega=omega)
    else:
        # A caller-supplied exterior, for cases whose outer field is known in
        # closed form -- which is what makes an analytic anchor possible.
        outer_cf = exterior_potential
        import warnings
        warnings.warn(
            "caller-supplied exterior_potential has unverified boundary data; "
            "validate it independently before accepting corrected losses",
            RuntimeWarning, stacklevel=2)
        data_error = None
    boundary = {"outer": outer_cf}

    if x_cut is not None:
        a_surface = surface_potential(
            panel_current, panel_ds, surface_impedance=surface_impedance,
            axial_E_field=axial_E_field, omega=omega)
        a_top = _nearest_panel_value(panel_xy, a_surface, (x_cut, y_hi))
        a_bot = _nearest_panel_value(panel_xy, a_surface, (x_cut, y_lo))
        # Deep inside the conductor the axial field vanishes, so A tends to
        # E0 / (j omega) rather than to zero; it is the DEVIATION from that
        # value that decays through the skin.  Matching the raw exponential
        # instead leaves J = sigma E0 flowing through the cross-section.
        a_deep = complex(axial_E_field) / (1j * omega)
        decay = (1.0 + 1j) / delta
        boundary["cut"] = (
            ng.CF(a_deep)
            + ng.CF(complex(a_top) - a_deep) * ng.exp(-decay * (y_hi - ng.y))
            + ng.CF(complex(a_bot) - a_deep) * ng.exp(-decay * (ng.y - y_lo)))

    dirichlet = "|".join(boundary)
    fes = ng.H1(mesh, order=order, complex=True, dirichlet=dirichlet)
    u, v = fes.TnT()
    sigma_cf = mesh.MaterialCF({conductor: float(sigma)}, default=0.0)
    form = ng.BilinearForm(fes, symmetric=True)
    form += ng.grad(u) * ng.grad(v) * ng.dx
    form += 1j * omega * MU0 * sigma_cf * u * v * ng.dx
    rhs = ng.LinearForm(fes)
    rhs += MU0 * sigma_cf * complex(axial_E_field) * v * ng.dx

    gf = ng.GridFunction(fes)
    gf.Set(mesh.BoundaryCF(boundary, default=0.0),
           definedon=mesh.Boundaries(dirichlet))
    form.Assemble()
    rhs.Assemble()
    residual = rhs.vec.CreateVector()
    residual.data = rhs.vec - form.mat * gf.vec
    gf.vec.data += form.mat.Inverse(freedofs=fes.FreeDofs(),
                                    inverse="umfpack") * residual

    return PatchSolution(
        x_cut_m=None if x_cut is None else float(x_cut), side=side,
        skin_depth_m=delta, ndof=int(fes.ndof),
        mesh=mesh, gf_A=gf, conductor=conductor, sigma_S_per_m=float(sigma),
        omega_rad_per_s=float(omega),
        axial_E_field_V_per_m=complex(axial_E_field),
        boundary_data_error=data_error)


def _chord_extent(panel_xy, x_cut):
    """Lowest and highest y where the outline crosses ``x = x_cut``.

    Read off the closed panel polygon rather than the CAD, so the chord is the
    one the global solve's own discretisation sees.  A cut that is not a simple
    chord -- fewer or more than two crossings -- is rejected by the caller
    through the full-thickness check.
    """
    xy = np.asarray(panel_xy, dtype=float)
    a = xy
    b = np.roll(xy, -1, axis=0)
    straddles = (a[:, 0] - x_cut) * (b[:, 0] - x_cut) < 0.0
    hits = list(xy[np.isclose(xy[:, 0], x_cut, rtol=0.0,
                              atol=1e-12 * max(abs(x_cut), 1e-6)), 1])
    if np.any(straddles):
        a, b = a[straddles], b[straddles]
        t = (x_cut - a[:, 0]) / (b[:, 0] - a[:, 0])
        hits.extend(a[:, 1] + t * (b[:, 1] - a[:, 1]))
    if not hits:
        raise ValueError(f"the outline does not cross x={x_cut:.4e} m")
    return float(np.min(hits)), float(np.max(hits))


def _nearest_panel_value(panel_xy, values, point):
    xy = np.asarray(panel_xy, dtype=float)
    target = np.asarray(point, dtype=float)
    return np.asarray(values)[int(np.argmin(
        np.linalg.norm(xy - target[None, :], axis=1)))]
