#!/usr/bin/env python
"""
Unified Benchmark Conditions for Linear Material Tests

This file defines the common parameters for comparing tetrahedral and hexahedral
element benchmarks with linear materials (matching nonlinear benchmark structure).

Author: Radia Development Team
Date: 2025-12-12
"""

import numpy as np

# =============================================================================
# Physical Constants
# =============================================================================
MU_0 = 4 * np.pi * 1e-7  # Vacuum permeability [T/(A/m)]

# =============================================================================
# Geometry Parameters
# =============================================================================
CUBE_SIZE = 1.0      # Cube edge length [m]
CUBE_HALF = 0.5      # Half of cube edge [m]

# =============================================================================
# Material Parameters (Linear)
# =============================================================================
MU_R = 1000          # Relative permeability (industry standard input)
# Note: As of v1.3.14, rad.MatLin() accepts mu_r directly (not chi)
# Internal conversion: chi = mu_r - 1
CHI = MU_R - 1       # For reference/calculations only

# =============================================================================
# External Field
# =============================================================================
H_EXT = 200000.0     # External H field magnitude [A/m] (same as nonlinear benchmark)
B_EXT = MU_0 * H_EXT  # External B field magnitude [T]

# =============================================================================
# Solver Parameters
# =============================================================================
SOLVER_TOLERANCE = 0.0001   # Relative tolerance for convergence
MAX_ITERATIONS = 1000       # Maximum iterations

# =============================================================================
# Mesh Parameters
# =============================================================================
# Hexahedron: n_div subdivisions per edge -> n_div^3 elements
# Tetrahedron: maxh = CUBE_SIZE / n_div (approximate correspondence)

# Standard mesh sizes for comparison
STANDARD_N_DIVS = [5, 10, 15]
STANDARD_MAXH = [0.35, 0.25, 0.2, 0.15, 0.1]

def n_div_to_maxh(n_div):
    """Convert hexahedron n_div to approximate tetrahedral maxh."""
    return CUBE_SIZE / n_div

def maxh_to_n_div(maxh):
    """Convert tetrahedral maxh to approximate hexahedron n_div."""
    return int(round(CUBE_SIZE / maxh))

# =============================================================================
# Analytical Solution for Validation
# =============================================================================
def M_analytical_z():
    """
    Analytical solution for uniform magnetization in a cube with linear material.

    For a linear material (chi) in external field H_ext:
    M = chi * H_internal

    For a cube, the demagnetizing factor N ~= 1/3 (approximate).
    H_internal = H_ext - N * M
    M = chi * (H_ext - N * M)
    M = chi * H_ext / (1 + chi * N)

    For high chi and N=1/3: M ~= H_ext / N = 3 * H_ext
    """
    N_demag = 1.0 / 3.0  # Approximate demagnetizing factor for cube
    M_z = CHI * H_EXT / (1 + CHI * N_demag)
    return M_z

# Reference value for validation
M_ANALYTICAL_Z = M_analytical_z()

# =============================================================================
# Summary of Conditions
# =============================================================================
def print_conditions():
    """Print summary of benchmark conditions."""
    print('=' * 70)
    print('LINEAR BENCHMARK CONDITIONS (Unified)')
    print('=' * 70)
    print()
    print('Geometry:')
    print('  Cube size:     %.1f m x %.1f m x %.1f m' % (CUBE_SIZE, CUBE_SIZE, CUBE_SIZE))
    print('  Center:        (0, 0, 0)')
    print()
    print('Material (Linear):')
    print('  mu_r:          %d' % MU_R)
    print('  chi:           %d' % CHI)
    print()
    print('External Field:')
    print('  H_ext:         %.0f A/m (z-direction)' % H_EXT)
    print('  B_ext:         %.6f T' % B_EXT)
    print()
    print('Analytical Solution (N_demag = 1/3):')
    print('  M_z:           %.0f A/m' % M_ANALYTICAL_Z)
    print()
    print('Solver:')
    print('  Tolerance:     %.5f' % SOLVER_TOLERANCE)
    print('  Max iter:      %d' % MAX_ITERATIONS)
    print()
    print('Mesh Correspondence (n_div <-> maxh):')
    print('  %8s  %10s  %12s  %12s' % ('n_div', 'maxh (m)', 'Hexa elems', 'Tetra elems'))
    print('  ' + '-' * 46)
    for n in STANDARD_N_DIVS:
        maxh = n_div_to_maxh(n)
        n_hexa = n**3
        # Approximate tetra count (empirical: ~5-6x hexahedra for same maxh)
        n_tetra_approx = int(n_hexa * 5.5)
        print('  %8d  %10.3f  %12d  %12s' % (n, maxh, n_hexa, '~%d' % n_tetra_approx))
    print()
    print('=' * 70)


if __name__ == '__main__':
    print_conditions()
