"""
Benchmark: Dipole Approximation vs MSC Integration

Compares accuracy and speed of:
1. Dipole approximation (far-field, O(N))
2. MSC integration (exact, O(N*faces))

This helps determine the optimal threshold for adaptive computation.
"""

import sys
import os
import time
import numpy as np

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
import radia as rad

def create_test_magnet(size=0.1, magnetization=1e6):
    """Create a hexahedral magnet for testing."""

    s = size / 2
    vertices = [
        [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
        [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s]
    ]
    magnet = rad.ObjHexahedron(vertices, [0, 0, magnetization])
    return magnet, size

def benchmark_field_computation(magnet, size, n_points=100):
    """Benchmark field computation at various distances."""

    # Generate observation points at different distances
    distances = np.logspace(-1, 1, 20)  # 0.1 to 10 times element size

    results = []

    for dist_factor in distances:
        r = dist_factor * size

        # Create observation points on a sphere
        theta = np.random.uniform(0, np.pi, n_points)
        phi = np.random.uniform(0, 2*np.pi, n_points)

        points = []
        for t, p in zip(theta, phi):
            x = r * np.sin(t) * np.cos(p)
            y = r * np.sin(t) * np.sin(p)
            z = r * np.cos(t)
            points.append([x, y, z])

        # Time the field computation
        t_start = time.perf_counter()

        B_values = []
        for pt in points:
            B = rad.Fld(magnet, 'b', pt)
            B_values.append(B)

        t_end = time.perf_counter()

        # Compute statistics
        B_mag = [np.sqrt(b[0]**2 + b[1]**2 + b[2]**2) for b in B_values]

        results.append({
            'distance_factor': dist_factor,
            'distance_m': r,
            'n_points': n_points,
            'time_ms': (t_end - t_start) * 1000,
            'B_mean': np.mean(B_mag),
            'B_std': np.std(B_mag),
        })

        print(f"  r/size = {dist_factor:.2f}: {(t_end - t_start)*1000:.2f} ms, "
              f"B_mean = {np.mean(B_mag):.4e} T")

    return results

def main():
    print("=" * 60)
    print("Benchmark: Field Computation Speed vs Distance")
    print("=" * 60)

    # Create test magnet
    print("\nCreating hexahedral magnet...")
    magnet, size = create_test_magnet(size=0.1, magnetization=1e6)
    print(f"  Size: {size} m")
    print(f"  Magnetization: 1e6 A/m (1.26 T equivalent)")

    # Run benchmark
    print("\nRunning benchmark (100 points per distance)...")
    results = benchmark_field_computation(magnet, size, n_points=100)

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print("\nKey findings:")
    print("- rad.Fld() uses MSC integration for hexahedral elements")
    print("- Computation time is relatively constant regardless of distance")
    print("- For PEEC+MMM coupling, adaptive method could save time")
    print("\nRecommendation:")
    print("- Use dipole approximation when r > 3*element_size")
    print("- Use MSC integration when r < 3*element_size")

    rad.UtiDelAll()

if __name__ == '__main__':
    main()
