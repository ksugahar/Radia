"""3D cuboid in VACUUM with uniform applied B_z — true infinite domain problem.

This is the PHYSICAL Case B (different from earlier closed-PEC version):
- Conductor: 5 x 2 x 1 mm Cu cuboid (sigma in OMG_cond)
- Surrounding: vacuum (sigma=0 everywhere else)
- Boundary: A_ind -> 0 at infinity
- Applied: A_ext = (B_0/2)(-y, x, 0) (uniform B_z everywhere)

V0 implementation: AIR BOX TRUNCATION
- Build a large air box around the conductor (10x conductor scale)
- Set A_ind x n = 0 on the far air box boundary (PEC truncation)
- Truncation error scales as (R_cond/R_box)^3 (magnetic dipole far field)

Future versions:
- v1: Kelvin transform for rigorous infinite-domain treatment
- v2: BEM-FEM coupling via ngbem

Tanimoto-Kameari adapted for "impressed A_ext" rather than "impressed J":
  Stage 0: J_0 = sigma * A_ext (impressed source, only in conductor)
           Solve curl-curl A_0 / mu_0 = J_0 on full domain (cond + air),
           A_ind x n = 0 on far box.
  R_0 = 1 / <J_0, J_0/sigma>_conductor = 1 / (sigma V (a^2+b^2)/12)
  L_0 = R_0^2 * <J_0, A_0>_conductor
  Recursion: J_{n+1} = J_n - sigma * Apot_n / L_n  (only in conductor)
"""
from netgen.occ import Box, Pnt, OCCGeometry, Glue
from ngsolve import (
    Mesh, HCurl, BilinearForm, LinearForm, GridFunction,
    CoefficientFunction, curl, dx, x, y, z, BitArray, BND,
    Integrate, TaskManager, ngsglobals,
)
from collections import deque
from math import pi
import json
from pathlib import Path

mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7

# Conductor (centered at origin)
ax, ay, az = 5e-3, 2e-3, 1e-3
V_cond = ax * ay * az

# Air box scale (k * conductor)
AIR_SCALE = 5  # box is k times conductor size, centered at origin
H_COND = 0.30e-3   # conductor mesh size
H_AIR = 1.5e-3     # air mesh size (coarser)
ORDER = 2          # tradeoff: order vs total ndof
N_STAGES = 5       # Kameari stages

# Reference R_0 (Parseval, conductor centered at origin):
#   R_0 = 1 / <J,J/sigma>_cond, J = sigma*A_ext = (sigma/2)(-y,x,0)
#   <J,J/sigma> = (sigma/4) integral (x^2+y^2) dV over [-a/2,a/2]^3
#               = sigma * V * (a^2+b^2) / 48
# So R_0_centered = 48 / (sigma V (a^2+b^2))
R0_anal = 48 / (sigma_Cu * V_cond * (ax**2 + ay**2))


def classify_vertices(mesh):
    """Boundary vertices = those on outer_box OR conductor_surface.
       For the vacuum problem, only outer_box vertices have Dirichlet BC.
       Tree-cotree mask should exclude only outer_box-vertex Whitney edges."""
    boundary_v = set()
    for el in mesh.Elements(BND):
        # Only outer_box has Dirichlet
        if "outer_box" in mesh.GetMaterial(el):
            for v in el.vertices:
                boundary_v.add(v.nr)
    # Actually mesh.GetMaterial doesn't apply here for BND; use bcname
    # Simpler: check via fes.FreeDofs after dirichlet="outer_box"
    return boundary_v


def build_spanning_tree_interior(mesh, free_dofs_check, fes):
    """BFS spanning tree restricted to vertices whose Whitney edge DoF is free.

    More robust than vertex-based classification: just check the FE space's
    FreeDofs to determine which edges contribute to kernel.
    """
    nv = mesh.nv
    # Mark vertices as "boundary" if all their incident edges have non-free
    # lowest-order DoF (i.e., they are kept in the kernel as constants).
    # Simpler: BFS over all vertices, mask first DoF of each tree edge if free.
    visited = [False] * nv
    tree_edges = []
    adj = [[] for _ in range(nv)]
    for ed in mesh.edges:
        v0, v1 = ed.vertices[0].nr, ed.vertices[1].nr
        adj[v0].append((v1, ed.nr))
        adj[v1].append((v0, ed.nr))
    visited[0] = True
    queue = deque([0])
    while queue:
        v = queue.popleft()
        for vn, edn in adj[v]:
            if not visited[vn]:
                visited[vn] = True
                tree_edges.append(edn)
                queue.append(vn)
    return tree_edges


def build_geometry():
    # Conductor centered at origin (so symmetry around origin for Kelvin v1)
    cuboid = Box(Pnt(-ax/2, -ay/2, -az/2), Pnt(ax/2, ay/2, az/2))
    cuboid.mat("conductor").bc("conductor_surface")
    cuboid.maxh = H_COND

    # Air box (also centered)
    Bx, By, Bz = AIR_SCALE * ax / 2, AIR_SCALE * ay / 2, AIR_SCALE * az / 2
    air = Box(Pnt(-Bx, -By, -Bz), Pnt(Bx, By, Bz))
    air.mat("air").bc("outer_box")
    air.maxh = H_AIR

    full = Glue([air - cuboid, cuboid])
    return OCCGeometry(full)


def kameari_vacuum(mesh, n_stages):
    # HCurl on full domain (cond + air); far box PEC; nograds for higher-order kernel
    fes = HCurl(mesh, order=ORDER, dirichlet="outer_box", nograds=True)
    print(f"  HCurl ndof = {fes.ndof}")
    print(f"  Materials = {mesh.GetMaterials()}")
    print(f"  Boundaries = {set(mesh.GetBoundaries())}")

    # Tree-cotree mask in-place into fes.FreeDofs() to remove
    # lowest-order Whitney gradient kernel.
    # BFS over all vertices, mask first DoF of each tree edge if currently free.
    tree_edges = build_spanning_tree_interior(mesh, None, fes)
    fd = fes.FreeDofs()
    masked = 0
    for edge_nr in tree_edges:
        edge = mesh.edges[edge_nr]
        dofs = fes.GetDofNrs(edge)
        if dofs and fd[dofs[0]]:
            fd[dofs[0]] = False
            masked += 1
    print(f"  Tree-cotree masked {masked} Whitney edge DoFs (in-place)")
    print(f"  Active DoFs after gauge: {sum(fd)}")

    # Material function: sigma in conductor, 0 in air (mu=mu0 everywhere)
    sigma_cf = mesh.MaterialCF({"conductor": sigma_Cu, "air": 0.0})

    u, v = fes.TnT()

    # External potential A_ext = (B_0/2)(-y, x, 0) with B_0 = 1 (normalized)
    A_ext = CoefficientFunction((-y, x, 0)) * 0.5

    # Initial impressed J_0 = sigma * A_ext (in conductor only)
    J = sigma_cf * A_ext

    # Stage 0: solve curl-curl A / mu_0 = J on full domain
    a_form = BilinearForm(fes)
    a_form += (1.0/mu0) * curl(u) * curl(v) * dx

    print("  Assembling curl-curl...")
    with TaskManager():
        a_form.Assemble()
        # Use sparsecholesky for clean direct solve
        # sigma in conductor only -> RHS only nonzero there
        # But the A_ind exists everywhere
        # The system is invertible because of full Dirichlet on outer_box
        # plus nograds=True kills higher-order gradient kernel
        # For lower-order Whitney kernel: PEC outer_box might or might not kill it
        # If sparsecholesky fails (singular), need to add tree-cotree mask later
        try:
            inv = a_form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        except Exception as e:
            print(f"  Sparsecholesky failed: {e}")
            print("  Need tree-cotree mask. Aborting v0.")
            return [], []

    sigma_inv_cf = mesh.MaterialCF({"conductor": 1.0/sigma_Cu, "air": 0.0})
    R_list, L_list = [], []
    Apot = None

    for n in range(n_stages):
        print(f"\n[Stage {n}]")
        f = LinearForm(fes)
        f += J * v * dx("conductor")  # source only in conductor
        with TaskManager():
            f.Assemble()
        gfA = GridFunction(fes, name=f"A_{n}")
        with TaskManager():
            gfA.vec.data = inv * f.vec

        # Compute R_n, L_n (Tanimoto accumulated)
        # R_n = 1 / <J_n, J_n / sigma>_cond
        R_inv_int = Integrate(J * J * sigma_inv_cf * dx("conductor"), mesh)
        if abs(R_inv_int) < 1e-25:
            print(f"  R_inv too small ({R_inv_int}), stop")
            break
        Rn = 1.0 / float(R_inv_int)
        if Apot is None:
            Apot = Rn * gfA
        else:
            Apot = Apot + Rn * gfA
        Ln_int = Integrate(J * Apot * dx("conductor"), mesh)
        Ln = Rn * float(Ln_int)
        R_list.append(Rn)
        L_list.append(Ln)
        tau_us = Ln / Rn * 1e6 if Rn > 0 and Ln > 0 else 0
        print(f"  R_{n} = {Rn:.6e}, L_{n} = {Ln:.6e} ({Ln*1e9:.4e} nH), tau = {tau_us:.4f} us")
        if Ln <= 0:
            print(f"  L_{n} <= 0, stop")
            break

        # Update J for next stage (only in conductor)
        J = J - sigma_cf * Apot / Ln

    return R_list, L_list


def main():
    ngsglobals.msg_level = 0
    print(f"3D cuboid in VACUUM v0 (air box truncation)")
    print(f"  Conductor: {ax*1000}x{ay*1000}x{az*1000} mm Cu")
    print(f"  Air box: {AIR_SCALE}x conductor (centered at origin)")
    print(f"  h_cond = {H_COND*1000} mm, h_air = {H_AIR*1000} mm, order = {ORDER}")
    print(f"  N_stages = {N_STAGES}")
    print()

    geo = build_geometry()
    print("Generating mesh...")
    mesh = Mesh(geo.GenerateMesh(maxh=H_AIR))
    print(f"  ne = {mesh.ne}, nv = {mesh.nv}")

    print(f"\nReference R_0 = 12/(sigma V (a^2+b^2)) = {R0_anal:.6e}")
    print(f"  (Parseval, exact, comes from <J,J/sigma> = sigma V (a^2+b^2)/12)")
    print()

    R, L = kameari_vacuum(mesh, N_STAGES)

    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"R_0 / R_0_analytic = {R[0]/R0_anal:.6f}  (should be 1.0 exactly)")

    out = {
        "method": "Cuboid in vacuum (air box truncation v0)",
        "geometry_mm": [ax*1000, ay*1000, az*1000],
        "air_scale": AIR_SCALE,
        "h_cond_mm": H_COND*1000,
        "h_air_mm": H_AIR*1000,
        "order": ORDER,
        "ne": mesh.ne,
        "R0_analytic": R0_anal,
        "R": R,
        "L": L,
    }
    out_path = Path(__file__).parent / "cuboid_521_vacuum_caseB_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
