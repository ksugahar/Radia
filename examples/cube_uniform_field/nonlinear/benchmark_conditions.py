#!/usr/bin/env python
"""
Unified Benchmark Conditions for Nonlinear Material Tests

This file defines the common parameters for comparing tetrahedral and hexahedral
element benchmarks with nonlinear materials.

Author: Radia Development Team
Date: 2025-12-06
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
# External Field
# =============================================================================
H_EXT = 50000.0      # External H field magnitude [A/m]
B_EXT = MU_0 * H_EXT  # External B field magnitude [T]

# =============================================================================
# B-H Curve Data (Nonlinear Material)
# =============================================================================
# Format: [H (A/m), B (T)]
# This represents a typical soft iron B-H curve
BH_DATA = [
    [0.0, 0.0],
    [100.0, 0.1],
    [200.0, 0.3],
    [500.0, 0.8],
    [1000.0, 1.2],
    [2000.0, 1.5],
    [5000.0, 1.7],
    [10000.0, 1.8],
    [50000.0, 2.0],
    [100000.0, 2.1],
]

# Convert B-H to H-M format for Radia's MatSatIsoTab
# M = B/mu_0 - H
HM_DATA = [[h, b/MU_0 - h] for h, b in BH_DATA]

# =============================================================================
# Solver Parameters
# =============================================================================
SOLVER_TOLERANCE = 0.001    # Relative tolerance for convergence
MAX_ITERATIONS = 1000       # Maximum nonlinear iterations

# =============================================================================
# Mesh Parameters
# =============================================================================
# Hexahedron: n_div subdivisions per edge -> n_div^3 elements
# Tetrahedron: maxh = CUBE_SIZE / n_div (approximate correspondence)

# Standard mesh sizes for comparison
# n_div: Number of divisions per cube edge
# For tetra: maxh ~= 1.0/n_div
STANDARD_N_DIVS = [3, 4, 5, 6, 7, 8]

def n_div_to_maxh(n_div):
    """Convert hexahedron n_div to approximate tetrahedral maxh."""
    return CUBE_SIZE / n_div

def maxh_to_n_div(maxh):
    """Convert tetrahedral maxh to approximate hexahedron n_div."""
    return int(round(CUBE_SIZE / maxh))

# =============================================================================
# Summary of Conditions
# =============================================================================
def print_conditions():
    """Print summary of benchmark conditions."""
    print('=' * 70)
    print('BENCHMARK CONDITIONS (Unified)')
    print('=' * 70)
    print()
    print('Geometry:')
    print('  Cube size:     %.1f m x %.1f m x %.1f m' % (CUBE_SIZE, CUBE_SIZE, CUBE_SIZE))
    print('  Center:        (0, 0, 0)')
    print()
    print('External Field:')
    print('  H_ext:         %.0f A/m (z-direction)' % H_EXT)
    print('  B_ext:         %.6f T' % B_EXT)
    print()
    print('Material (B-H curve):')
    print('  %10s  %10s  %12s' % ('H (A/m)', 'B (T)', 'mu_r'))
    print('  ' + '-' * 36)
    for h, b in BH_DATA:
        if h > 0:
            mu_r = b / (MU_0 * h)
            print('  %10.0f  %10.3f  %12.1f' % (h, b, mu_r))
        else:
            print('  %10.0f  %10.3f  %12s' % (h, b, '-'))
    print()
    print('Solver:')
    print('  Tolerance:     %.4f' % SOLVER_TOLERANCE)
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
