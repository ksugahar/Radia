"""
NGSolve 2D Adaptive Mesh Refinement with Laplacian Smoothing
2次元アダプティブメッシュ + 節点移動（Laplacian Smoothing）

This example demonstrates:
  1. Creating 2D geometry (L-shaped domain)
  2. Adaptive mesh refinement with ZZ error estimator
  3. Laplacian mesh smoothing after each refinement
  4. VTU output at each step
"""
import os
from numpy import pi, sqrt, cos, sin, linspace, zeros, nan, isnan, meshgrid, array, log
from ngsolve import *
from netgen.geom2d import SplineGeometry

print("=" * 60)
print("NGSolve 2D Adaptive Mesh + Laplacian Smoothing")
print("=" * 60)

# ============================================================
# Parameters
# ============================================================
L = 1.0                 # Domain size [m]
maxh_initial = 0.25     # Initial mesh size [m] (coarser)
order = 2               # Polynomial order
max_iterations = 8      # Number of refinement iterations
theta = 0.5             # Dörfler marking parameter

# Smoothing parameters
smoothing_iterations = 3    # Number of Laplacian smoothing iterations
smoothing_weight = 0.4      # Relaxation weight (0 < w < 1)

print(f"\nParameters:")
print(f"  Domain size: {L} m x {L} m")
print(f"  Initial mesh size: {maxh_initial} m")
print(f"  Polynomial order: {order}")
print(f"  Max iterations: {max_iterations}")
print(f"  Doerfler theta: {theta}")
print(f"  Smoothing iterations: {smoothing_iterations}")
print(f"  Smoothing weight: {smoothing_weight}")


# ============================================================
# Geometry: L-shaped domain using SplineGeometry
# ============================================================
print("\nCreating L-shaped geometry...")

geo = SplineGeometry()

# L-shaped domain vertices
# (0,0) -> (1,0) -> (1,0.5) -> (0.5,0.5) -> (0.5,1) -> (0,1) -> (0,0)
pts = [(0, 0), (L, 0), (L, L/2), (L/2, L/2), (L/2, L), (0, L)]
p = [geo.AddPoint(*pt) for pt in pts]

# Add lines (counter-clockwise)
geo.Append(["line", p[0], p[1]], bc="bottom")
geo.Append(["line", p[1], p[2]], bc="right")
geo.Append(["line", p[2], p[3]], bc="inner")
geo.Append(["line", p[3], p[4]], bc="inner")
geo.Append(["line", p[4], p[5]], bc="top")
geo.Append(["line", p[5], p[0]], bc="left")

print(f"  L-shaped domain created")


# ============================================================
# Solve Poisson equation
# ============================================================
def solve_poisson(mesh, order):
    """Solve Poisson equation on given mesh"""
    fes = H1(mesh, order=order, dirichlet=".*")
    u, v = fes.TnT()

    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx
    a.Assemble()

    f = LinearForm(fes)
    f += 1 * v * dx
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

    return fes, gfu


# ============================================================
# ZZ Error Estimator
# ============================================================
def compute_error_estimator(mesh, fes, gfu):
    """Compute ZZ-type error estimator using flux recovery"""
    flux = grad(gfu)

    fes_flux = HDiv(mesh, order=fes.globalorder - 1)
    gf_flux = GridFunction(fes_flux)
    gf_flux.Set(flux)

    err = (flux - gf_flux) * (flux - gf_flux)
    element_errors = Integrate(err, mesh, element_wise=True)

    return element_errors


# ============================================================
# Dörfler Marking
# ============================================================
def mark_elements_doerfler(mesh, element_errors, theta):
    """Dörfler marking strategy"""
    sorted_indices = sorted(range(len(element_errors)),
                           key=lambda i: element_errors[i], reverse=True)

    total_error = sum(element_errors)
    marked_error = 0
    marked = []

    for idx in sorted_indices:
        marked.append(idx)
        marked_error += element_errors[idx]
        if marked_error >= theta * total_error:
            break

    return marked


# ============================================================
# Laplacian Smoothing
# ============================================================
def laplacian_smoothing(mesh, n_iter, weight):
    """
    Laplacian mesh smoothing with boundary constraints.

    For each interior vertex, move towards neighbor centroid:
      x_new = x_old + weight * (centroid - x_old)

    Boundary vertices: slide along boundary edges only.
    """
    print(f"    Applying Laplacian smoothing ({n_iter} iterations)...")

    ngmesh = mesh.ngmesh
    pnts = ngmesh.Points()
    n_vertices = len(pnts)

    # Build vertex adjacency
    neighbors = {i: set() for i in range(1, n_vertices + 1)}
    for el in ngmesh.Elements2D():
        verts = [v.nr for v in el.vertices]
        for i, v in enumerate(verts):
            for j, w in enumerate(verts):
                if i != j:
                    neighbors[v].add(w)

    # Identify boundary vertices
    boundary_vertices = set()
    boundary_edges = {}  # vertex -> list of boundary neighbor vertices

    for seg in ngmesh.Elements1D():
        verts = [v.nr for v in seg.vertices]
        for v in verts:
            boundary_vertices.add(v)
            if v not in boundary_edges:
                boundary_edges[v] = []
        # Add boundary edge connections
        if len(verts) == 2:
            boundary_edges[verts[0]].append(verts[1])
            boundary_edges[verts[1]].append(verts[0])

    # Corner vertices (have exactly 2 boundary neighbors)
    corner_vertices = set()
    for v in boundary_vertices:
        if len(boundary_edges.get(v, [])) != 2:
            corner_vertices.add(v)

    # Smoothing iterations
    for iteration in range(n_iter):
        new_positions = {}

        for vid in range(1, n_vertices + 1):
            px = pnts[vid][0]
            py = pnts[vid][1]

            if vid in corner_vertices:
                # Corner vertex - don't move
                new_positions[vid] = (px, py)

            elif vid in boundary_vertices:
                # Boundary vertex - keep fixed (simpler and safer)
                new_positions[vid] = (px, py)

            else:
                # Interior vertex - standard Laplacian smoothing
                if len(neighbors[vid]) > 0:
                    avg_x = sum(pnts[n][0] for n in neighbors[vid]) / len(neighbors[vid])
                    avg_y = sum(pnts[n][1] for n in neighbors[vid]) / len(neighbors[vid])

                    new_x = px + weight * (avg_x - px)
                    new_y = py + weight * (avg_y - py)

                    new_positions[vid] = (new_x, new_y)
                else:
                    new_positions[vid] = (px, py)

        # Update positions
        for vid in range(1, n_vertices + 1):
            pnts[vid][0] = new_positions[vid][0]
            pnts[vid][1] = new_positions[vid][1]

    ngmesh.Update()
    print(f"    Smoothing completed")


# ============================================================
# VTK Output
# ============================================================
def output_vtk(mesh, gfu, iteration):
    """Output mesh and solution to VTK file"""
    vtk_filename = f"mesh_iter_{iteration:02d}"
    vtk = VTKOutput(mesh, coefs=[gfu], names=["solution"],
                    filename=vtk_filename, subdivision=0)
    vtk.Do()
    return vtk_filename + ".vtu"


# ============================================================
# Main: Adaptive Refinement Loop
# ============================================================
print("\n" + "=" * 60)
print("Starting Adaptive Mesh Refinement")
print("=" * 60)

# Generate initial mesh
mesh = Mesh(geo.GenerateMesh(maxh=maxh_initial))

print(f"\nInitial mesh:")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")

# History tracking
history = {
    'ndof': [],
    'elements': [],
    'error': []
}

for iteration in range(max_iterations):
    print(f"\n{'=' * 60}")
    print(f"Iteration {iteration + 1}/{max_iterations}")
    print("=" * 60)

    # Solve
    fes, gfu = solve_poisson(mesh, order)

    # Compute error
    element_errors = compute_error_estimator(mesh, fes, gfu)
    total_error = sqrt(sum(element_errors))

    # Record history
    history['ndof'].append(fes.ndof)
    history['elements'].append(mesh.ne)
    history['error'].append(total_error)

    print(f"  Elements: {mesh.ne}")
    print(f"  Vertices: {mesh.nv}")
    print(f"  DOFs: {fes.ndof}")
    print(f"  Error estimator: {total_error:.6e}")

    # Output VTK
    vtk_file = output_vtk(mesh, gfu, iteration)
    print(f"  VTK saved: {vtk_file}")\

    # Mark and refine (except last iteration)
    if iteration < max_iterations - 1:
        marked = mark_elements_doerfler(mesh, element_errors, theta)
        print(f"  Marked elements: {len(marked)}")

        for el in mesh.Elements():
            mesh.SetRefinementFlag(el, False)
        for el_nr in marked:
            mesh.SetRefinementFlag(ElementId(VOL, el_nr), True)

        mesh.Refine()

        # Apply Laplacian smoothing after refinement
        laplacian_smoothing(mesh, smoothing_iterations, smoothing_weight)

# ============================================================
# Final Statistics
# ============================================================
print("\n" + "=" * 60)
print("Final Statistics")
print("=" * 60)

print(f"\n{'Iter':<6} {'Elements':<10} {'DOFs':<10} {'Error':<12}")
print("-" * 38)
for i in range(len(history['ndof'])):
    print(f"{i+1:<6} {history['elements'][i]:<10} {history['ndof'][i]:<10} "
          f"{history['error'][i]:<12.4e}")

print(f"\nInitial -> Final:")
print(f"  Elements: {history['elements'][0]} -> {history['elements'][-1]}")
print(f"  DOFs: {history['ndof'][0]} -> {history['ndof'][-1]}")
print(f"  Error: {history['error'][0]:.4e} -> {history['error'][-1]:.4e} "
      f"({history['error'][0]/history['error'][-1]:.1f}x reduction)")

# ============================================================
# Visualization
# ============================================================
print("\nGenerating plots...")

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('mathtext', **{'rm': 'serif', 'it': 'serif:italic',
                              'bf': 'serif:bold', 'fontset': 'cm'})

fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)

# Plot 1: Final mesh
ax1 = axes[0, 0]
for el in mesh.Elements():
    vertices = [mesh[v].point for v in el.vertices]
    xs = [v[0] for v in vertices] + [vertices[0][0]]
    ys = [v[1] for v in vertices] + [vertices[0][1]]
    ax1.plot(xs, ys, 'b-', linewidth=0.3)

ax1.set_xlabel('$x$ (m)')
ax1.set_ylabel('$y$ (m)')
ax1.set_title(f'Final Mesh ({mesh.ne} elements, with smoothing)')
ax1.set_aspect('equal')
ax1.set_xlim(-0.05, L + 0.05)
ax1.set_ylim(-0.05, L + 0.05)

# Plot 2: Solution contour
ax2 = axes[0, 1]
n_sample = 50
x_sample = linspace(0.01, L - 0.01, n_sample)
y_sample = linspace(0.01, L - 0.01, n_sample)
[xx, yy] = meshgrid(x_sample, y_sample)
uu = zeros(xx.shape)

for iy in range(n_sample):
    for ix in range(n_sample):
        x_val, y_val = x_sample[ix], y_sample[iy]
        if x_val > L/2 + 0.01 and y_val > L/2 + 0.01:
            uu[iy, ix] = nan
        else:
            try:
                uu[iy, ix] = gfu(mesh(x_val, y_val))
            except:
                uu[iy, ix] = nan

contour = ax2.contourf(xx, yy, uu, levels=20, cmap='viridis')
plt.colorbar(contour, ax=ax2, label='$u$')
ax2.set_xlabel('$x$ (m)')
ax2.set_ylabel('$y$ (m)')
ax2.set_title('Solution $u(x,y)$')
ax2.set_aspect('equal')

# Plot 3: Convergence
ax3 = axes[1, 0]
ax3.loglog(history['ndof'], history['error'], 'ko-', linewidth=2, markersize=6, label='Adaptive')
ndof_ref = array(history['ndof'])
err_ref = history['error'][0] * (ndof_ref[0] / ndof_ref) ** (order / 2)
ax3.loglog(ndof_ref, err_ref, 'r--', linewidth=1.5, label=f'$O(N^{{-{order}/2}})$')
ax3.set_xlabel('DOFs')
ax3.set_ylabel('Error Estimator')
ax3.set_title('Convergence')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xticks([10**2, 10**3])
ax3.set_xticklabels(['$10^2$', '$10^3$'])

# Plot 4: Error vs Iteration
ax4 = axes[1, 1]
iterations = range(1, len(history['error']) + 1)
ax4.semilogy(iterations, history['error'], 'ko-', linewidth=2, markersize=6)
ax4.set_xlabel('Iteration')
ax4.set_ylabel('Error Estimator')
ax4.set_title('Error vs Iteration')
ax4.grid(True, alpha=0.3)
ax4.set_xticks(list(iterations))

plt.tight_layout()

png_file = "convergence.png"
plt.savefig(png_file, dpi=150, bbox_inches='tight')
print(f"  Plot saved to: {png_file}")

os.startfile(png_file)

print("\n" + "=" * 60)
print("Computation completed successfully")
print("=" * 60)
