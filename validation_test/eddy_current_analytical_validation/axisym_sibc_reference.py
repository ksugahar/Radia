"""Independent axisymmetric A-form FEM-SIBC reference, exp(+i omega t).

Validation-only in-memory OCC geometry. Scalar impedance and impressed
uniform-J circular coil are explicit inputs; no Radia BEM operator is reused.
"""
import math
import numpy as np
import ngsolve as ng
from netgen.occ import WorkPlane, Glue, OCCGeometry


def solve(*, bore=0.0, outer=.15, maxh=.006, boundary_h=.001,
          coil_a=.0005, order=3, frequency=1000, sigma=5.8e7):
    mu0 = 4e-7 * math.pi
    omega = 2 * math.pi * frequency
    zs = (1+1j) * math.sqrt(omega*mu0/(2*sigma))
    box = WorkPlane().MoveTo(0, -outer).Rectangle(outer, 2*outer).Face()
    for edge in box.edges:
        edge.name = "axis" if abs(edge.center.x) < 1e-12 else "outer"
    body = WorkPlane().MoveTo(bore, -.0125).Rectangle(.025-bore, .025).Face()
    for edge in body.edges:
        edge.name = "sibc"
        edge.maxh = boundary_h
    coil = WorkPlane().MoveTo(.03+coil_a, 0).Direction(0, 1).Arc(coil_a, 180).Arc(coil_a, 180).Face()
    coil.name = "coil"
    coil.maxh = coil_a/2
    air = box-body-coil
    air.name = "air"
    mesh = ng.Mesh(OCCGeometry(Glue([air, coil]), dim=2).GenerateMesh(maxh=maxh))
    mesh.Curve(order)
    fes = ng.H1(mesh, order=order, complex=True, dirichlet="axis|outer")
    u, v = fes.TnT()
    r = ng.x
    form = ng.BilinearForm(fes)
    form += ng.grad(u)*ng.grad(v)/(mu0*r)*ng.dx
    # Outward normal of the AIR points into the conductor. Positive j is
    # essential: the opposite sign gives an active, not passive, boundary.
    form += 1j*omega/zs*u*v/r*ng.ds("sibc")
    load = ng.LinearForm(fes)
    load += v/(math.pi*coil_a**2)*ng.dx("coil")
    form.Assemble()
    load.Assemble()
    field = ng.GridFunction(fes)
    field.vec.data = form.mat.Inverse(fes.FreeDofs(), inverse="pardiso")*load.vec
    pairing = np.dot(load.vec.FV().NumPy(), field.vec.FV().NumPy())
    reaction = -math.pi*omega*pairing.imag
    loss = math.pi*omega**2*zs.real/abs(zs)**2 * ng.Integrate(
        (field*ng.Conj(field)).real/r, mesh, ng.BND,
        definedon=mesh.Boundaries("sibc"), order=2*order+4)
    # On the outer cylindrical metal surface, H_z=+j omega A_phi/Zs.
    z = np.linspace(-.010, .010, 21)
    hz = np.array([1j*omega/zs*field(mesh(.025, float(zi)))/.025 for zi in z])
    return dict(P_surface=float(loss), P_reaction=float(reaction),
                power_balance_relative_error=float(abs(loss-reaction)/loss),
                ndof=fes.ndof, bore=bore, outer=outer, boundary_h=boundary_h,
                order=order, coil_a=coil_a,
                z=z.tolist(), H_z_real=hz.real.tolist(), H_z_imag=hz.imag.tolist())
