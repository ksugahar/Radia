"""Debug: find why VectorEddyCurrentFEMBEM system is singular."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4e-7 * np.pi

from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh, BND
from ngbem_eddy import VectorEddyCurrentFEMBEM, create_conductor_mesh

width, height, depth = 0.02, 0.02, 0.01
sigma_al = 3.7e7
freq = 100e3
maxh = 0.01

mesh = create_conductor_mesh(width, height, depth, maxh=maxh)
n_vol = sum(1 for _ in mesh.Elements())
n_bnd = sum(1 for _ in mesh.Elements(BND))
print(f"Mesh: {n_vol} vol, {n_bnd} bnd elements")

solver = VectorEddyCurrentFEMBEM(mesh, sigma_al, mu_r=1.0, order=1)
solver.assemble(kappa=0.01, intorder=4)

n1 = solver._n_hcurl
n2 = solver._n_hdiv
n3 = solver._n_l2
print(f"\nDOFs: H(curl)={n1}, HDivSurface={n2}, SurfaceL2={n3}")

# Check rank of B
B = solver._B_np
U_b, s_b, Vt_b = np.linalg.svd(B, full_matrices=False)
rank_B = np.sum(s_b > 1e-10 * s_b[0])
print(f"\nB matrix ({n2}x{n1}):")
print(f"  Singular values: min={s_b[-1]:.6e}, max={s_b[0]:.6e}")
print(f"  Rank: {rank_B} (of {min(n2, n1)})")
print(f"  Nonzero entries: {np.count_nonzero(np.abs(B) > 1e-15)}")

# Check rank of S_bem
S = solver._S_bem_np
U_s, s_s, Vt_s = np.linalg.svd(S, full_matrices=False)
rank_S = np.sum(s_s > 1e-10 * s_s[0])
print(f"\nS_bem matrix ({n2+n3}x{n2+n3}):")
print(f"  Singular values: min={s_s[-1]:.6e}, max={s_s[0]:.6e}")
print(f"  Rank: {rank_S} (of {n2+n3})")

# Check rank of individual BEM blocks
A_k = solver._A_k_np
s_Ak = np.linalg.svd(A_k, compute_uv=False)
rank_Ak = np.sum(s_Ak > 1e-10 * s_Ak[0])
print(f"\nA_k matrix ({n2}x{n2}):")
print(f"  Singular values: min={s_Ak[-1]:.6e}, max={s_Ak[0]:.6e}")
print(f"  Rank: {rank_Ak}")

V_k = solver._V_k_np
s_Vk = np.linalg.svd(V_k, compute_uv=False)
rank_Vk = np.sum(s_Vk > 1e-10 * s_Vk[0])
print(f"\nV_k matrix ({n3}x{n3}):")
print(f"  Singular values: min={s_Vk[-1]:.6e}, max={s_Vk[0]:.6e}")
print(f"  Rank: {rank_Vk}")

Q_k = solver._Q_k_np
s_Qk = np.linalg.svd(Q_k, compute_uv=False)
rank_Qk = np.sum(s_Qk > 1e-10 * s_Qk[0])
print(f"\nQ_k matrix ({n3}x{n2}):")
print(f"  Singular values: min={s_Qk[-1]:.6e}, max={s_Qk[0]:.6e}")
print(f"  Rank: {rank_Qk}")

# Check rank of a_FEM
omega = 2 * np.pi * freq
a_FEM = solver._a_curl_np + 1j * omega * sigma_al * solver._a_mass_np
s_aFEM = np.linalg.svd(a_FEM, compute_uv=False)
rank_aFEM = np.sum(s_aFEM > 1e-10 * s_aFEM[0])
print(f"\na_FEM matrix ({n1}x{n1}):")
print(f"  Singular values: min={s_aFEM[-1]:.6e}, max={s_aFEM[0]:.6e}")
print(f"  Rank: {rank_aFEM}")

# Build and check the full system matrix
N = n1 + n2 + n3
A_sys = np.zeros((N, N), dtype=complex)
A_sys[:n1, :n1] = a_FEM
A_sys[:n1, n1:n1+n2] = solver._B_np.T  # +B^T
A_sys[n1:n1+n2, :n1] = solver._B_np    # +B
A_sys[n1:, n1:] = -solver._S_bem_np    # -S_bem

s_full = np.linalg.svd(A_sys, compute_uv=False)
rank_full = np.sum(s_full > 1e-10 * s_full[0])
print(f"\nFull system matrix ({N}x{N}):")
print(f"  Singular values: min={s_full[-1]:.6e}, max={s_full[0]:.6e}")
print(f"  Rank: {rank_full}")
print(f"  Null space dimension: {N - rank_full}")
print(f"  Condition number: {s_full[0]/max(s_full[-1], 1e-300):.2e}")

# Also try with -B^T (original sign)
A_sys2 = np.zeros((N, N), dtype=complex)
A_sys2[:n1, :n1] = a_FEM
A_sys2[:n1, n1:n1+n2] = -solver._B_np.T  # -B^T
A_sys2[n1:n1+n2, :n1] = solver._B_np     # +B
A_sys2[n1:, n1:] = -solver._S_bem_np     # -S_bem

s_full2 = np.linalg.svd(A_sys2, compute_uv=False)
rank_full2 = np.sum(s_full2 > 1e-10 * s_full2[0])
print(f"\nSystem with -B^T ({N}x{N}):")
print(f"  Rank: {rank_full2}")
print(f"  Smallest SV: {s_full2[-1]:.6e}")

# Try with +S_bem (positive BEM)
A_sys3 = np.zeros((N, N), dtype=complex)
A_sys3[:n1, :n1] = a_FEM
A_sys3[:n1, n1:n1+n2] = solver._B_np.T
A_sys3[n1:n1+n2, :n1] = solver._B_np
A_sys3[n1:, n1:] = +solver._S_bem_np  # +S_bem

s_full3 = np.linalg.svd(A_sys3, compute_uv=False)
rank_full3 = np.sum(s_full3 > 1e-10 * s_full3[0])
print(f"\nSystem with +B^T, +S_bem ({N}x{N}):")
print(f"  Rank: {rank_full3}")
print(f"  Smallest SV: {s_full3[-1]:.6e}")

# Check null vector
if rank_full < N:
    U, s, Vt = np.linalg.svd(A_sys)
    null_vec = Vt[-1, :]  # last row = smallest singular value direction
    print(f"\n  Null vector structure:")
    print(f"    |A part| = {np.linalg.norm(null_vec[:n1]):.6e}")
    print(f"    |j part| = {np.linalg.norm(null_vec[n1:n1+n2]):.6e}")
    print(f"    |rho part| = {np.linalg.norm(null_vec[n1+n2:]):.6e}")

# Check FreeDofs
print("\n--- FreeDofs check ---")
fd_hc = solver._fes_hcurl.FreeDofs()
fd_hd = solver._fes_hdiv.FreeDofs()
fd_l2 = solver._fes_l2.FreeDofs()
n_free_hc = sum(1 for i in range(n1) if fd_hc[i])
n_free_hd = sum(1 for i in range(n2) if fd_hd[i])
n_free_l2 = sum(1 for i in range(n3) if fd_l2[i])
print(f"  H(curl): {n_free_hc} free of {n1}")
print(f"  HDivSurface: {n_free_hd} free of {n2}")
print(f"  SurfaceL2: {n_free_l2} free of {n3}")

# List non-free DOFs
nonfree_hd = [i for i in range(n2) if not fd_hd[i]]
nonfree_l2 = [i for i in range(n3) if not fd_l2[i]]
if nonfree_hd:
    print(f"  Non-free HDivSurface DOFs: {nonfree_hd}")
if nonfree_l2:
    print(f"  Non-free SurfaceL2 DOFs: {nonfree_l2}")

# Check zero rows/cols in A_k
print("\n--- Zero row/col analysis ---")
Ak_row_norms = np.linalg.norm(A_k, axis=1)
Ak_col_norms = np.linalg.norm(A_k, axis=0)
zero_rows_Ak = np.where(Ak_row_norms < 1e-15)[0]
zero_cols_Ak = np.where(Ak_col_norms < 1e-15)[0]
print(f"  A_k: {len(zero_rows_Ak)} zero rows, {len(zero_cols_Ak)} zero cols")
if len(zero_rows_Ak) > 0:
    print(f"    Zero row indices: {zero_rows_Ak.tolist()}")

B_row_norms = np.linalg.norm(B, axis=1)
B_col_norms = np.linalg.norm(B, axis=0)
zero_rows_B = np.where(B_row_norms < 1e-15)[0]
print(f"  B: {len(zero_rows_B)} zero rows")
if len(zero_rows_B) > 0:
    print(f"    Zero row indices: {zero_rows_B.tolist()}")

# Common zero rows (DOFs dead in both operators)
common_zero = set(zero_rows_Ak.tolist()) & set(zero_rows_B.tolist())
print(f"  Common zero rows (dead DOFs): {len(common_zero)}")
if common_zero:
    print(f"    Indices: {sorted(common_zero)}")

# Check S_bem zero rows
S_row_norms = np.linalg.norm(S, axis=1)
S_col_norms = np.linalg.norm(S, axis=0)
zero_rows_S = np.where(S_row_norms < 1e-15)[0]
zero_cols_S = np.where(S_col_norms < 1e-15)[0]
print(f"  S_bem: {len(zero_rows_S)} zero rows, {len(zero_cols_S)} zero cols")
if len(zero_rows_S) > 0:
    print(f"    Zero row indices: {zero_rows_S.tolist()}")

# Smallest 10 singular values of A_k
print(f"\n  A_k smallest 10 SVs: {s_Ak[-10:]}")
# Smallest 10 singular values of S_bem
print(f"  S_bem smallest 10 SVs: {s_s[-10:]}")
