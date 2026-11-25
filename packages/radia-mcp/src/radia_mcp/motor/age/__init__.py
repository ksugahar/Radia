"""Air-gap element (AGE) method notes served by radia-mcp.motor ``motor_age_quality``.

This package must not import ``radia``; the AGE no-radia-import exception
applies here. The text was moved unchanged from the former ngsolve-mcp.age server.
"""

AGE_TOPICS = {
    'harmonic_coupling': "AGE replaces the physical air-gap mesh with an analytical Fourier-series\n"
        "harmonic coupling between stator and rotor cylindrical interfaces.",
    'torque': "Electromagnetic torque is computed directly from the AGE Fourier coefficients\n"
        "without remeshing across arbitrary rotor angles.",
    'periodic_motion': "Supports fractional-slot / anti-periodic sector symmetry and smooth\n"
        "angular sweeps for cogging torque and back-EMF.",
}
