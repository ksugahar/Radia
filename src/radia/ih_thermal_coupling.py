"""IH nonlinear heat law on caller-assembled quadrature maps, SI/kelvin.

NGSolve owns the space, basis, gradients, measures and geometry. This module
only contracts supplied evaluation maps with Radia material laws. It neither
constructs an FE space nor substitutes nodal temperatures for coefficients.
The dense CPU Newton route is an explicit bounded reference, not a GPU fallback.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

SIGMA = 5.670374419e-8
ENERGY_RELATIVE_LIMIT = 1e-6
POWER_FLOOR_W = 1e-12


def port_heat_balance(current_peak_A, voltage_peak_V, coil_resistance_ohm, heat_power_W):
    """Independent port observation versus integrated heat, peak phasors.

    Voltage must be obtained from the physical reaction/port operator, not
    inferred from the heat under test. Resistive coil power is accounted once.
    """
    i, v = complex(current_peak_A), complex(voltage_peak_V)
    r, heat = float(coil_resistance_ohm), float(heat_power_W)
    if not all(np.isfinite(x) for x in (i.real, i.imag, v.real, v.imag, r, heat)) or r < 0:
        raise ValueError("Finite port observation and nonnegative coil resistance required")
    power, coil = .5*(v*i.conjugate()).real, .5*r*abs(i)**2
    defect = power-coil-heat
    relative = abs(defect)/max(abs(power), abs(coil)+abs(heat), POWER_FLOOR_W)
    result = dict(input_power_W=power, coil_loss_W=coil, heat_power_W=heat,
        defect_W=defect, relative_defect=relative, limit=ENERGY_RELATIVE_LIMIT)
    if relative > ENERGY_RELATIVE_LIMIT:
        raise ValueError("Independent port/heat energy balance exceeds 1e-6")
    return result


def _array(value, shape, name):
    a = np.array(value, dtype=float, copy=True)
    if a.shape != shape or not np.all(np.isfinite(a)):
        raise ValueError(f"Invalid {name} shape or finite values")
    a.setflags(write=False)
    return a


@dataclass(frozen=True)
class ThermalLaw:
    """Per-volume-quadrature cp/k and per-boundary-quadrature emissivity.

    Table rows are quadrature samples, columns increasing absolute kelvin.
    Enthalpy integrates the piecewise-linear cp exactly, including partial
    intervals. Density is constant mass density, not electrical resistivity.
    """
    kelvin: np.ndarray
    cp: np.ndarray
    kappa: np.ndarray
    emissivity: np.ndarray

    def __post_init__(self):
        axis = np.asarray(self.kelvin, float)
        if axis.ndim != 1 or len(axis) < 2 or not np.all(np.isfinite(axis)) or axis[0] <= 0 or np.any(np.diff(axis) <= 0):
            raise ValueError("Strictly increasing positive Kelvin axis required")
        axis = _array(axis, axis.shape, 'kelvin')
        cp = np.asarray(self.cp, float)
        if cp.ndim != 2 or cp.shape[1] != len(axis) or not cp.shape[0]:
            raise ValueError("cp table must be volume-quadrature by temperature")
        cp = _array(cp, cp.shape, 'cp')
        k = _array(self.kappa, cp.shape, 'kappa')
        e = np.asarray(self.emissivity, float)
        if e.ndim != 2 or e.shape[1] != len(axis):
            raise ValueError("emissivity table must be boundary-quadrature by temperature")
        e = _array(e, e.shape, 'emissivity')
        if np.any(cp <= 0) or np.any(k <= 0) or np.any((e < 0) | (e > 1)):
            raise ValueError("Positive cp/kappa and emissivity in [0,1] required")
        for key, value in [('kelvin', axis), ('cp', cp), ('kappa', k), ('emissivity', e)]:
            object.__setattr__(self, key, value)

    def evaluate(self, table, temperature):
        t = np.asarray(temperature, float)
        if t.shape != (len(table),) or not np.all(np.isfinite(t)) or np.any(t < self.kelvin[0]) or np.any(t > self.kelvin[-1]):
            raise ValueError("Temperature outside declared Kelvin material table")
        j = np.clip(np.searchsorted(self.kelvin, t, side='right') - 1, 0, len(self.kelvin)-2)
        row = np.arange(len(t))
        slope = (table[row, j+1]-table[row, j]) / (self.kelvin[j+1]-self.kelvin[j])
        return table[row, j] + slope*(t-self.kelvin[j]), slope

    def enthalpy(self, temperature):
        t = np.asarray(temperature, float)
        cp, _ = self.evaluate(self.cp, t)
        j = np.clip(np.searchsorted(self.kelvin, t, side='right')-1, 0, len(self.kelvin)-2)
        cumulative = np.c_[np.zeros(len(t)), np.cumsum(.5*(self.cp[:, :-1]+self.cp[:, 1:])*np.diff(self.kelvin), axis=1)]
        row = np.arange(len(t))
        return cumulative[row, j] + .5*(self.cp[row, j]+cp)*(t-self.kelvin[j])


@dataclass(frozen=True)
class ThermalQuadrature:
    basis: np.ndarray
    gradient: np.ndarray
    volume_weights: np.ndarray
    density: np.ndarray
    boundary_basis: np.ndarray
    boundary_weights: np.ndarray
    convection: np.ndarray
    constant_coefficients: np.ndarray
    law: ThermalLaw

    def __post_init__(self):
        e = np.asarray(self.basis, float)
        if e.ndim != 2 or not 0 < e.shape[1] <= 256 or not 0 < e.shape[0] <= 4096:
            raise ValueError("Bounded quadrature reference requires 1..256 coefficients and 1..4096 volume samples")
        q, n = e.shape
        g = np.asarray(self.gradient, float)
        if g.ndim != 3 or g.shape[0] != q or g.shape[2] != n or not 1 <= g.shape[1] <= 3:
            raise ValueError("Gradient must be volume-quadrature by dimension by coefficient")
        b = np.asarray(self.boundary_basis, float)
        if b.ndim != 2 or b.shape[1] != n or b.shape[0] > 4096:
            raise ValueError("Boundary evaluation shape mismatch")
        arrays = dict(basis=(e, e.shape), gradient=(g, g.shape), boundary_basis=(b, b.shape),
            volume_weights=(self.volume_weights, (q,)), density=(self.density, (q,)),
            boundary_weights=(self.boundary_weights, (len(b),)), convection=(self.convection, (len(b),)),
            constant_coefficients=(self.constant_coefficients, (n,)))
        for key, (value, shape) in arrays.items():
            object.__setattr__(self, key, _array(value, shape, key))
        if np.any(self.volume_weights <= 0) or np.any(self.density <= 0) or np.any(self.boundary_weights <= 0) or np.any(self.convection < 0):
            raise ValueError("Positive quadrature measures/density and nonnegative convection required")
        c = self.constant_coefficients
        if not np.allclose(e@c, 1, rtol=0, atol=1e-12) or not np.allclose(b@c, 1, rtol=0, atol=1e-12) or not np.allclose(g@c, 0, rtol=0, atol=1e-12):
            raise ValueError("Constant field is not preserved by supplied FE maps")
        if self.law.cp.shape[0] != q or self.law.emissivity.shape[0] != len(b):
            raise ValueError("Material rows must follow quadrature numbering")

    def residual(self, candidate, accepted, load, ambient, dt):
        n = self.basis.shape[1]
        x = _array(candidate, (n,), 'candidate'); old = _array(accepted, (n,), 'accepted')
        f = _array(load, (n,), 'thermal load W')
        if not np.isfinite(dt) or dt <= 0 or not np.isfinite(ambient) or ambient <= 0:
            raise ValueError("Positive time step and ambient Kelvin required")
        e, g, b = self.basis, self.gradient, self.boundary_basis
        tq, oldq, tb = e@x, e@old, b@x
        cp, _ = self.law.evaluate(self.law.cp, tq)
        k, dk = self.law.evaluate(self.law.kappa, tq)
        eps, deps = self.law.evaluate(self.law.emissivity, tb)
        de = self.volume_weights*self.density*(self.law.enthalpy(tq)-self.law.enthalpy(oldq))/dt
        grad = g@x
        boundary = self.convection*(tb-ambient) + eps*SIGMA*(tb**4-ambient**4)
        dboundary = self.convection + SIGMA*(deps*(tb**4-ambient**4)+4*eps*tb**3)
        conduction = np.einsum('qdi,qd,q->i', g, grad, self.volume_weights*k)
        storage_vector = e.T@de
        boundary_vector = b.T@(self.boundary_weights*boundary)
        r = storage_vector + conduction + boundary_vector - f
        residual_scale = max(float(np.linalg.norm(v)) for v in
            (f, storage_vector, conduction, boundary_vector))
        residual_scale = max(residual_scale, POWER_FLOOR_W)
        j = e.T@((self.volume_weights*self.density*cp/dt)[:, None]*e)
        j += np.einsum('qdi,qdj,q->ij', g, g, self.volume_weights*k)
        j += np.einsum('qdi,qd,qj,q->ij', g, grad, e, self.volume_weights*dk)
        j += b.T@((self.boundary_weights*dboundary)[:, None]*b)
        power = float(self.constant_coefficients@f)
        storage = float(np.sum(de))
        loss = float(self.boundary_weights@boundary)
        defect = power-storage-loss
        scale = max(abs(power), abs(storage)+abs(loss), POWER_FLOOR_W)
        balance = dict(input_power_W=power, storage_power_W=storage, boundary_loss_W=loss,
            defect_W=defect, relative_defect=abs(defect)/scale,
            limit=ENERGY_RELATIVE_LIMIT, power_floor_W=POWER_FLOOR_W,
            residual_scale_W=residual_scale)
        return r, j, balance

    def step(self, accepted, load, ambient, dt, *, tolerance=1e-10, max_iterations=30):
        """Explicit CPU reference; transactional return, no mutation on failure."""
        if not 0 < tolerance <= 1e-6 or not 1 <= max_iterations <= 100:
            raise ValueError("Invalid nonlinear tolerance/iteration bound")
        x = np.array(accepted, float, copy=True)
        for iteration in range(max_iterations):
            r, j, balance = self.residual(x, accepted, load, ambient, dt)
            relative = float(np.linalg.norm(r))/balance['residual_scale_W']
            if relative <= tolerance and balance['relative_defect'] <= ENERGY_RELATIVE_LIMIT:
                return x, dict(iterations=iteration, true_residual=relative, balance=balance,
                    scheme='backward-Euler enthalpy; EM load evaluated at accepted previous temperature')
            increment = np.linalg.solve(j, -r)
            linear_residual = np.linalg.norm(j@increment+r)/max(np.linalg.norm(r), POWER_FLOOR_W)
            if not np.isfinite(linear_residual) or linear_residual > 1e-6:
                raise ValueError("Thermal Newton linear true residual exceeds 1e-6")
            x += increment
        raise ValueError("Thermal Newton failed; accepted temperature unchanged")

    def native_config(self):
        """Explicit quadrature contract for IH native kernels; no FE inference."""
        return dict(thermal_material_mode='quadrature-enthalpy-v1',
            thermal_temperature_axis_K=self.law.kelvin.tolist(),
            thermal_basis=self.basis.tolist(), thermal_gradient=self.gradient.reshape(-1, self.basis.shape[1]).tolist(),
            thermal_dimension=self.gradient.shape[1], thermal_volume_weights_m3=self.volume_weights.tolist(),
            thermal_density_kg_per_m3=self.density.tolist(), thermal_cp_J_per_kgK=self.law.cp.tolist(),
            thermal_kappa_W_per_mK=self.law.kappa.tolist(), thermal_boundary_basis=self.boundary_basis.tolist(),
            thermal_boundary_weights_m2=self.boundary_weights.tolist(), thermal_convection_W_per_m2K=self.convection.tolist(),
            thermal_emissivity=self.law.emissivity.tolist(), thermal_nonlinear_tolerance=1e-10,
            thermal_energy_relative_limit=ENERGY_RELATIVE_LIMIT, thermal_energy_power_floor_W=POWER_FLOOR_W,
            thermal_max_iterations=30)
