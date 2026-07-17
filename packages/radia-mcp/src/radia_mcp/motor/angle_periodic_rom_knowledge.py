"""HCurl Eddy Bubble + HDiv-MMM angle-periodic motor ROM knowledge."""
from __future__ import annotations


SECTIONS: dict[str, str] = {
    "architecture": """\
# HCurl Eddy Bubble + HDiv-MMM Motor Architecture

Use BDM1/BDM2 HDiv-MMM for magnetization and HCurl high-order parent spaces
for divergence-free eddy current.  Reduction is role-aware: an air/exterior
conductor face may carry a surface/SIBC trace, a conductor-conductor face must
preserve loop transport, and ordinary low-response interior coordinates may
be removed by mixed-Galerkin response compression.  Do not select p=6 as a
universal answer; increase the parent order until the retained response,
corner/edge observables, and passivity gates stop changing.
""",
    "face_policy": """\
# Face-Adaptive Eddy Reduction

Classify every conductive face before eliminating a DoF.  The production
roles are conductor-air, conductor-exterior, conductor-insulator,
conductor-conductor, and conductive-interface.  Only conductor-air/exterior
is an SIBC half-space face.  Conductor-insulator may need a non-SIBC trace;
conductor-conductor and conductive-interface faces form the conductive graph.
Reduce their many bridge coordinates to a cycle basis, never delete them as
local bubbles, because eddy current must remain able to circulate globally.
""",
    "angle_rom": """\
# Angle-Periodic Energy Law

The generalized current order is physical phase currents followed by internal
eddy coordinates.  Build odd, uniform, endpoint-exclusive angle tables for
L(theta), R(theta), permanent-magnet flux, v-cross-B motion flux, and scalar
cogging coenergy.  Fourier interpolation gives analytic angle derivatives;
continuous skew multiplies each harmonic by its sinc factor.  One coenergy law
drives reluctance, PM, cogging, and virtual-work torque, while the same motion
flux vector appears in speed voltage and Lorentz torque to preserve power.
""",
    "time_domain": """\
# Time Domain, SIBC, Hysteresis, and Temperature

Advance flux implicitly with mechanical angle and speed.  Temperature scales
the resistance by a positive diagonal congruence, preserving positive
semidefiniteness; optional end-winding L/R and a lumped thermal state remain
explicit model parameters.  A frequency-domain sqrt(s) DtN/SIBC term must be
realized as a positive-real CLN state system before time stepping and must
never be replaced silently by a constant resistance.  Hysteresis uses a pure
trial callback; commit its restart state only after the coupled step succeeds.
""",
    "ports": """\
# Simulink, C ABI, and FMI Boundary

`radia.motor_rom_export.SaveMotorROMBundle` writes canonical NPZ arrays, a MAT
mirror and MATLAB loader, a JSON port/unit manifest, and an FMI 3.0.2 variable
fragment.  `src/core/rad_motor_rom_c.h` is C ABI version 1 with row-major
tables and an opaque handle.  The bundle is an FMI source boundary, not a
packaged FMU.  A wrapper must retain phase/eddy ordering, units, the pure-trial
hysteresis rule, and the discrete energy-balance output.
""",
    "mesh_gate": """\
# Curved Mesh Gate

Before a solve, sample the actual high-order element transformation, not only
the linear corner tetrahedron.  Reject non-positive Jacobians, near-collapsed
scaled Jacobians, missing material/boundary labels, and SIBC labels attached
to a conductor-insulator or conductor-conductor face.  The production
`check-vol` path exposes curve order, mapping samples, face-role counts, SIBC
candidates, and loop-bridge counts.  Curved air-gap geometry is part of the
torque discretization, not cosmetic postprocessing.
""",
    "validation": """\
# Current Numerical Locks

The curved 8-pole/24-physical-slot angle benchmark uses 33 training and 33
interlaced hold-out angles.  Its committed result reports phase-flux relative
error below 2e-15, Maxwell-vs-ROM torque relative RMSE about 1.90e-3, and
direct virtual-work-vs-ROM relative RMSE below 4e-9.  The C ABI/Python
1000-step lock reports current differences near 1e-14 A and power-balance
residuals below 1.1e-8 W.  These are regression locks for this benchmark, not a
claim that every motor needs the same parent order or angle sample count.
""",
    "limits": """\
# Honest Scope Boundary

The reusable motor ROM, curved mesh gate, C ABI, skew/end/temperature law,
functional hysteresis port, and angle-table validation are implemented.  A
specific production machine still needs its own 3D Cubit mesh, material and
winding data, measured B-H/hysteresis law, positive-real SIBC realization,
and operating envelope.  End-region and skew corrections must be checked on
that geometry; the 2D 24-slot benchmark does not replace a final 3D design
qualification.
""",
}


def get_angle_periodic_rom_knowledge(topic: str = "architecture") -> str:
    key = (topic or "architecture").strip().lower()
    if key in ("all", "*"):
        return "\n\n---\n\n".join(SECTIONS[name] for name in SECTIONS)
    if key not in SECTIONS:
        return "Unknown topic %r. Valid topics: %s (or 'all')." % (
            topic,
            ", ".join(sorted(SECTIONS)),
        )
    return SECTIONS[key]


__all__ = ["SECTIONS", "get_angle_periodic_rom_knowledge"]
