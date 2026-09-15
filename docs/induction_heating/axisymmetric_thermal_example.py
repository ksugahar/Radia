"""Local subroutine for the executed axisymmetric thermal notebook.

Diagnostic in-memory geometry only; not the production VOL/CAD entry point.
"""
import hashlib
import importlib.util
import sys

import ngsolve as ng
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


def reproduce(root):
    from validation_test.induction_heating._axisym_test_mesh import make_axisymmetric_mesh

    source = root / 'src/radia/simulink/ih_operator_assembly.py'
    spec = importlib.util.spec_from_file_location('ih_notebook_assembly', source)
    assembly = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = assembly
    spec.loader.exec_module(assembly)
    mesh = make_axisymmetric_mesh()
    radius, k, rho, cp, t0, a, dt = .025, 46.6, 7800., 467., 293.15, 1e4, .1
    q = 2*k*a*radius
    weight = 2*np.pi*ng.x
    options = assembly.IHOperatorAssemblyOptions(
        thermal_order=2, workpiece_label='outer', convection_W_per_m2K=0)
    with ng.TaskManager():
        # Use the checkout's production assembly kernel, without replacing
        # an installed Radia package or pretending this is a VOL-label gate.
        op = assembly._assemble_p2_thermal(
            mesh, np.full(mesh.nv, q), options, weight, 'outer')
        constant = op.constant_coefficients
        load = np.asarray(op.heat_to_temperature_projection).reshape(
            op.n_temperature, -1) @ op.unit_heat_density_W_per_m3
        power = q*2*np.pi*radius*.025
        power_error = abs(constant @ load/power-1)
        capacity_error = abs(constant @ op.temperature_cell_weights_J_per_K /
                             (rho*cp*np.pi*radius**2*.025)-1)
        assert power_error < 1e-12 and capacity_error < 1e-12

        fes = ng.H1(mesh, order=2)
        u, v = fes.TnT()
        mass = ng.BilinearForm(fes)
        mass += weight*rho*cp*u*v*ng.dx
        mass.Assemble()
        step = ng.BilinearForm(fes)
        step += weight*(rho*cp*u*v+dt*k*ng.grad(u)*ng.grad(v))*ng.dx
        step.Assemble()
        rhs = ng.LinearForm(fes)
        rhs += weight*q*v*ng.ds('outer')
        rhs.Assemble()
        temperature = ng.GridFunction(fes)
        temperature.Set(t0+a*ng.x**2)
        coefficients = temperature.vec.FV().NumPy().copy()
        inverse = step.mat.Inverse(fes.FreeDofs(), inverse='sparsecholesky')
        errors = []
        for n in range(1, 101):
            b = temperature.vec.CreateVector()
            b.data = mass.mat*temperature.vec+dt*rhs.vec
            temperature.vec.data = inverse*b
            exact = t0+a*ng.x**2+4*k*a*n*dt/(rho*cp)
            errors.append(float(np.sqrt(
                ng.Integrate(weight*(temperature-exact)**2, mesh) /
                ng.Integrate(weight*(exact-t0)**2, mesh))))
        assert max(errors) < 1e-9

        # Independently consume the production CSR/mixed-load ABI in SciPy.
        # This is not presented as a fresh MATLAB/MEX execution.
        m = csr_matrix((op.mass_value, op.mass_col, op.mass_row_ptr))
        stiffness = csr_matrix((op.stiffness_value, op.stiffness_col, op.stiffness_row_ptr))
        for _ in range(100):
            coefficients = spsolve(m+dt*stiffness, m@coefficients+dt*load)
        difference = ng.GridFunction(fes)
        difference.vec.FV().NumPy()[:] = coefficients-temperature.vec.FV().NumPy()
        rms = float(np.sqrt(ng.Integrate(weight*difference**2, mesh) /
                            ng.Integrate(weight, mesh)))
        assert rms < 1e-7
    return mesh, temperature, {
        'NGSolve': ng.__version__, 'assembly_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'example_sha256': hashlib.sha256(__file_bytes()).hexdigest(),
        'temperature_DOFs': fes.ndof, 'steps': 100, 'heat_power_W': power,
        'power_relative_error': power_error, 'capacity_relative_error': capacity_error,
        'max_exact_relative_rise_L2': max(errors), 'matrix_vs_independent_RMS_K': rms,
        'scope': 'LAB in-memory algebra reproduction, not production CAD/VOL acceptance',
    }


def __file_bytes():
    from pathlib import Path
    return Path(__file__).read_bytes()
