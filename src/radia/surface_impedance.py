"""Source-side P1 surface forms for piecewise-constant passive impedance.

Panel order is exactly ``mesh.Elements(BND)`` order. A tagged value prevents
accidentally interpreting a vertex array as a panel array. Caller owns TaskManager.
"""
from dataclasses import dataclass

import numpy as np


def require_legacy_vertex_layout(record):
    """Legacy visualization/export must never reinterpret triangle fields as vertices."""
    if record.get('esim_impedance_layout') not in (None, 'legacy-vertex-order'):
        raise ValueError('This legacy vertex consumer cannot read BND-element-order fields; '
                         'use the current solver qsurf/Ht artifacts instead')


def panel_mesh_identity(mesh):
    """Hash ordered triangle coordinates; vertex renumbering alone is harmless."""
    import hashlib
    from ngsolve import BND
    triangles = np.asarray([[tuple(mesh.vertices[v.nr].point) for v in el.vertices]
                            for el in mesh.Elements(BND)], dtype='<f8')
    if triangles.ndim != 3 or triangles.shape[1:] != (3, 3):
        raise ValueError('Panel impedance requires triangular 3D surfaces')
    return hashlib.sha256(triangles.tobytes()).hexdigest(), triangles.mean(axis=1)


def write_panel_impedance(path, mesh, values, *, frequency_hz):
    """Write mesh-bound JSON in ohms; values use peak-phasor convention."""
    import json
    from pathlib import Path
    impedance = PanelSurfaceImpedance(values)
    identity, centroids = panel_mesh_identity(mesh)
    if len(centroids) != len(impedance.values):
        raise ValueError('One impedance is required per surface triangle')
    if not np.isfinite(frequency_hz) or frequency_hz <= 0:
        raise ValueError('frequency_hz must be finite and positive')
    payload = dict(schema='radia.panel_surface_impedance.v1', layout='BND-element-order',
                   unit='ohm', frequency_hz=float(frequency_hz), surface_sha256=identity,
                   centroids_m=centroids.tolist(), real_ohm=impedance.values.real.tolist(),
                   imag_ohm=impedance.values.imag.tolist())
    Path(path).write_text(json.dumps(payload, indent=2, allow_nan=False), encoding='utf-8')


def read_panel_impedance(path, mesh, *, frequency_hz):
    """Reject stale mesh order, frequency, invalid shape and active/nonfinite Zs."""
    import json
    from pathlib import Path
    data = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    identity, centroids = panel_mesh_identity(mesh)
    if (data.get('schema') != 'radia.panel_surface_impedance.v1'
            or data.get('layout') != 'BND-element-order' or data.get('unit') != 'ohm'
            or data.get('surface_sha256') != identity):
        raise ValueError('Panel Zs schema, units or surface mesh/order mismatch')
    frequency = float(data.get('frequency_hz', float('nan')))
    if not np.isfinite(frequency) or not np.isclose(frequency, frequency_hz, rtol=1e-12, atol=0):
        raise ValueError('Panel Zs frequency mismatch')
    real, imag = np.asarray(data['real_ohm'], float), np.asarray(data['imag_ohm'], float)
    if real.shape != (len(centroids),) or imag.shape != real.shape:
        raise ValueError('One real/imag impedance is required per surface triangle')
    return PanelSurfaceImpedance(real + 1j*imag)


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
    from scipy.sparse import coo_matrix, csr_matrix
    u, v = fes.TnT()
    matrix = csr_matrix((fes.ndof, fes.ndof), dtype=complex)
    for imaginary, factor in ((False, 1), (True, 1j)):
        z = impedance.coefficient(fes, imaginary=imaginary)
        form = BilinearForm(fes)
        form += z * InnerProduct(grad(u).Trace(), grad(v).Trace()) * ds
        form.Assemble()
        rows, cols, values = form.mat.COO()
        matrix += factor * coo_matrix((values, (rows, cols)), shape=matrix.shape).tocsr()
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


def main(argv=None):
    """Export an editable, mesh-bound panel table for the Simulink assembler."""
    import argparse
    from ngsolve import Mesh, BND
    from radia.panels.calc_inductance import _extract_bnd_only_inline
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('workpiece')
    parser.add_argument('--label', default='sibc')
    parser.add_argument('--frequency', required=True, type=float)
    parser.add_argument('--real-ohm', required=True, type=float)
    parser.add_argument('--imag-ohm', required=True, type=float)
    parser.add_argument('--output', required=True)
    args = parser.parse_args(argv)
    mesh = Mesh(args.workpiece)
    if args.label not in mesh.GetBoundaries():
        raise ValueError('Template label must name a workpiece boundary')
    surface = _extract_bnd_only_inline(mesh, args.label)
    values = np.full(surface.GetNE(BND), complex(args.real_ohm, args.imag_ohm))
    write_panel_impedance(args.output, surface, values, frequency_hz=args.frequency)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
