# -*- coding: utf-8 -*-
"""2D->3D machine-global scaling for radia-ngsolve motor workflows.

A 2D (x, y) machine cross-section is solved PER UNIT AXIAL LENGTH and often on a
1/N rotational SECTOR.  To report physical, whole-machine SI globals (torque,
flux linkage, force, loss) the raw FE-primitive outputs must be scaled by the
axial stack length and the symmetry factor -- exactly ONELAB/GetDP's
``AxialLength`` and ``SymmetryFactor`` constants (machine_magstadyn_a.pro)::

    Flux               = SymmetryFactor * AxialLength * ... * CompZ[a]
    GlobalTerm Dof{Ub} / SymmetryFactor ...          (a sector carries 1/N)
    Rb[]               = Factor_R_3DEffects * AxialLength * FillFactor * N^2 / SurfCoil / sigma

This module gives those scalings a single, tested home instead of the ad-hoc
``MM_TORQUE_TO_NM = 1e-6 * depth`` constant that was hand-rolled per test.  It
complements the analytic winding / leakage layer in
:mod:`radia_mcp.radia_ngsolve.solve` (``skew_factor``,
``slot_leakage_inductance``, ``magnetizing_inductance_per_phase``).

Verified against the ONELAB ElectricMachines source (2026-06-09); see
``radia_mcp.motor.onelab_knowledge`` section ``twod_corrections`` (#1, #2, #4).
ONELAB is open-source (ULiege-ULB) -- the Rb[] form below is ported verbatim and
attributed, not re-derived.
"""
from __future__ import annotations

from dataclasses import dataclass

MU0 = 4e-7 * 3.141592653589793


@dataclass(frozen=True)
class MachineScaling:
    """Convert per-unit-depth, per-sector 2D FE outputs to whole-machine SI globals.

    Parameters
    ----------
    stack_length : float
        Axial active length l_stk [m] (ONELAB ``AxialLength``).
    symmetry_factor : float, default 1.0
        Whole machine / modelled sector, N (ONELAB ``SymmetryFactor``).  A
        rotational (anti)periodic sector solve yields 1/N of each extensive
        global, so the physical value is N x the sector value.  1 = full machine.
    length_unit_m : float, default 1.0
        Mesh length unit in metres.  Use ``1e-3`` when the geometry is meshed in
        millimetres (the convention in the radia-motor cogging / cogging-dq tests).

    Notes
    -----
    :meth:`torque` carries ``length_unit_m ** 2`` (the 2D Maxwell-stress torque
    density integrates force/length x lever-arm = two powers of length); generic
    extensive globals use :pyattr:`global_scale` = ``stack_length*symmetry_factor``
    and the caller supplies any remaining unit power specific to their quantity.
    """

    stack_length: float
    symmetry_factor: float = 1.0
    length_unit_m: float = 1.0

    @property
    def global_scale(self) -> float:
        """``AxialLength * SymmetryFactor`` -- the bare extensive-global multiplier."""
        return self.stack_length * self.symmetry_factor

    def torque(self, tau_2d: float) -> float:
        """Physical torque [N*m] from a 2D ``eggshell_torque_2d`` output.

        ``radia_ngsolve.force.eggshell_torque_2d`` returns the weighted
        Maxwell-stress torque PER UNIT AXIAL LENGTH in (mesh-length)^2 field
        units; the physical whole-machine torque is
        ``tau_2d * length_unit_m**2 * stack_length * symmetry_factor``.  With a
        millimetre mesh (``length_unit_m=1e-3``) and ``symmetry_factor=1`` this
        reproduces the ``1e-6 * l_stk`` constant used in the motor regression tests.
        """
        return tau_2d * self.length_unit_m ** 2 * self.global_scale

    def force(self, f_2d: float) -> float:
        """Physical force [N] from a 2D per-unit-depth ``eggshell_force_2d`` output.

        Force density integrates one lever-free power less than torque, so it
        carries a single ``length_unit_m`` (mesh-unit) factor in addition to the
        ``AxialLength * SymmetryFactor`` machine scaling.
        """
        return f_2d * self.length_unit_m * self.global_scale

    def flux_linkage(self, lambda_per_depth: float) -> float:
        """Whole-machine flux linkage [Wb] from a per-unit-depth, per-sector value.

        Applies the ONELAB ``SymmetryFactor * AxialLength`` factor only.  The flux
        value must already be in SI per unit depth (the A_z mesh-length-unit
        handling depends on the solve's source / reluctivity units, so it stays
        with the caller; only the machine scaling is applied here).
        """
        return lambda_per_depth * self.global_scale


def winding_resistance_onelab(turns, surf_coil, sigma, stack_length,
                              factor_3d_effects=1.0, fill_factor=1.0):
    """Coil/phase DC resistance replicating ONELAB's ``Rb[]`` (machine_magstadyn_a.pro)::

        Rb = Factor_R_3DEffects * AxialLength * FillFactor * N^2 / SurfCoil / sigma

    ``factor_3d_effects`` (>= 1) inflates the in-slot 2D resistance to account for
    the END-TURN copper length the 2D model cannot see; ``fill_factor`` follows
    ONELAB's convention (it multiplies, as written above).  Ported verbatim from
    the open ONELAB reference, not re-derived -- see ``twod_corrections`` #4.
    """
    return factor_3d_effects * stack_length * fill_factor * turns * turns / (surf_coil * sigma)
