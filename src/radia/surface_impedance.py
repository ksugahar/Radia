"""Source-side P1 surface forms for piecewise-constant passive impedance.

Panel order is exactly ``mesh.Elements(BND)`` order. A tagged value prevents
accidentally interpreting a vertex array as a panel array. Caller owns TaskManager.
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PanelSurfaceImpedance:
    """One peak-phasor impedance [ohm] per BND triangle, in mesh order."""

    values: np.ndarray

    def __post_init__(self):
        values = np.array(self.values, dtype=complex, copy=True)
        if (values.ndim != 1 or not values.size
                or not np.all(np.isfinite(values)) or np.any(values.real < 0)):
            raise ValueError("Panel Z_s must be a nonempty finite passive vector")
        values.flags.writeable = False
        object.__setattr__(self, "values", values)

    def coefficient(self, fes, *, imaginary=False):
        from ngsolve import BND, GridFunction, SurfaceL2
        mesh = fes.mesh
        elements = list(mesh.Elements(BND))
        if (fes.ndof != mesh.nv or len(elements) != len(self.values)
                or any(len(el.vertices) != 3 for el in elements)):
            raise ValueError("Panel Z_s requires P1 surface triangles and one value per BND element")
        space = SurfaceL2(mesh, order=0)
        coefficient = GridFunction(space)
        values = self.values.imag if imaginary else self.values.real
        for el, value in zip(elements, values):
            dofs = space.GetDofNrs(el)
            if len(dofs) != 1 or dofs[0] < 0:
                raise ValueError("Expected one constant surface DOF per panel")
            coefficient.vec[dofs[0]] = float(value)
        return coefficient


def assemble_panel_impedance_stiffness(fes, impedance):
    """K_Z[i,j] = integral Z_s grad_s(phi_i).grad_s(phi_j) dS.

    NGSolve owns mapped gradients and quadrature. This uses the same weighted
    Gram principle as VIM, in the BIE's scalar potential basis (not VIM modes).
    """
    from ngsolve import BilinearForm, InnerProduct, ds, grad
    from scipy.sparse import coo_matrix
    u, v = fes.TnT()
    matrix = np.zeros((fes.ndof, fes.ndof), dtype=complex)
    for imaginary, factor in ((False, 1), (True, 1j)):
        z = impedance.coefficient(fes, imaginary=imaginary)
        form = BilinearForm(fes)
        form += z * InnerProduct(grad(u).Trace(), grad(v).Trace()) * ds
        form.Assemble()
        rows, cols, values = form.mat.COO()
        matrix += factor * coo_matrix((values, (rows, cols)), shape=matrix.shape).toarray()
    return matrix


def panel_tangential_field_rms(fes, phi_vec):
    """Per-panel sqrt(mean |grad_s phi|^2), with peak phasor coefficients."""
    from ngsolve import BND, CF, GridFunction, InnerProduct, Integrate, grad
    phi = np.asarray(phi_vec, dtype=complex)
    if phi.shape != (fes.ndof,) or not np.all(np.isfinite(phi)):
        raise ValueError("Invalid surface potential")
    real, imag = GridFunction(fes), GridFunction(fes)
    real.vec.FV().NumPy()[:] = phi.real
    imag.vec.FV().NumPy()[:] = phi.imag
    h2 = InnerProduct(grad(real), grad(real)) + InnerProduct(grad(imag), grad(imag))
    area = np.asarray(Integrate(CF(1), fes.mesh, BND, element_wise=True))
    integral = np.asarray(Integrate(h2, fes.mesh, BND, element_wise=True))
    if np.any(area <= 0):
        raise ValueError("Surface panels require positive area")
    return np.sqrt(np.maximum(integral / area, 0))
