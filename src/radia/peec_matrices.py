"""
PEEC Matrix Construction - C++ pybind11 Module Required

This module provides Loop-Star PEEC matrix construction for quasi-static
electromagnetic analysis using Darwin approximation.

IMPORTANT: This module requires the C++ pybind11 module (peec_matrices.pyd).
No Python fallback is provided. Build Radia with Build.ps1 to generate
the required .pyd file.

Usage:
    from radia.peec_matrices import PEECBuilder, create_wire_peec, create_loop_peec

    # Method 1: Using PEECBuilder class
    builder = PEECBuilder()
    builder.create_wire([0,0,0], [0.1,0,0], 1e-3, 1e-3, 10)
    L, R, P, M_LS = builder.build()

    # Port impedance at 1 MHz
    Z = builder.compute_impedance(1e6, [1.0] + [0.0]*9)

    # Method 2: Using convenience functions
    L, R, P, M_LS = create_wire_peec([0,0,0], [0.1,0,0], 1e-3, 1e-3, 10)

Reference:
    A. Ruehli, "Equivalent Circuit Models for Three-Dimensional Multiconductor
    Systems," IEEE Trans. MTT, 1974.
"""

# Import from C++ pybind11 module (REQUIRED - no fallback)
try:
    from radia.peec_matrices import PEECBuilder, create_wire_peec, create_loop_peec
except ImportError as e:
    raise ImportError(
        "PEEC C++ module (peec_matrices.pyd) not found. "
        "Build Radia with Build.ps1 to generate the required module. "
        f"Original error: {e}"
    ) from e

# Export public API
__all__ = ['PEECBuilder', 'create_wire_peec', 'create_loop_peec']
