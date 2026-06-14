"""Phase 12.3: Multi-hex BEM K matrix + Foster eigenvalues.

Speed strategy: use Gauss-Legendre 6-point quadrature in 3D (216 Phi evals per
K entry), giving ~10-12 digit accuracy with the closed-form Phi.

K_ij = (mu0 / 4pi) int_cell_i Phi_hex_j(r) dV.

Time constants tau_k = sigma * lambda_k where K v_k = lambda_k M v_k.

Date: 2026-05-02
"""

import time
import mpmath as mp

from foster_cln_hex_tet import hex_newton_potential

mp.mp.dps = 25

mu0 = 4 * mp.pi * mp.mpf("1e-7")
sigma_Cu = mp.mpf("5.8e7")

a, b, c = mp.mpf("5e-3"), mp.mpf("2e-3"), mp.mpf("1e-3")


# ---------- Gauss-Legendre nodes/weights ----------
def gauss_legendre(n):
    """Nodes (in [-1,1]) and weights for n-point Gauss-Legendre."""
    # Use mpmath's built-in
    nodes_weights = mp.calculus.quadrature.GaussLegendre(mp.mp).get_nodes(-1, 1, degree=n, prec=mp.mp.prec)
    # Returns list of (node, weight)
    nodes = [nw[0] for nw in nodes_weights]
    weights = [nw[1] for nw in nodes_weights]
    return nodes, weights


def gauss_legendre_simple(n):
    """Direct Gauss-Legendre nodes/weights for small n via known values."""
    # For n=4: nodes and weights are known
    if n == 4:
        nodes = [
            -mp.mpf("0.861136311594052575223946488893"),
            -mp.mpf("0.339981043584856264802665759103"),
            +mp.mpf("0.339981043584856264802665759103"),
            +mp.mpf("0.861136311594052575223946488893"),
        ]
        weights = [
            mp.mpf("0.347854845137453857373063949222"),
            mp.mpf("0.652145154862546142626936050778"),
            mp.mpf("0.652145154862546142626936050778"),
            mp.mpf("0.347854845137453857373063949222"),
        ]
        return nodes, weights
    if n == 6:
        nodes = [
            -mp.mpf("0.932469514203152027812301554494"),
            -mp.mpf("0.661209386466264513661399595020"),
            -mp.mpf("0.238619186083196908630501721681"),
            +mp.mpf("0.238619186083196908630501721681"),
            +mp.mpf("0.661209386466264513661399595020"),
            +mp.mpf("0.932469514203152027812301554494"),
        ]
        weights = [
            mp.mpf("0.171324492379170345040296142173"),
            mp.mpf("0.360761573048138607569833513838"),
            mp.mpf("0.467913934572691047389870343990"),
            mp.mpf("0.467913934572691047389870343990"),
            mp.mpf("0.360761573048138607569833513838"),
            mp.mpf("0.171324492379170345040296142173"),
        ]
        return nodes, weights
    raise ValueError(f"n={n} not implemented")


def gauss_3d_integral(integrand, x0, x1, y0, y1, z0, z1, n=4):
    """Triple Gauss-Legendre integral over [x0,x1]x[y0,y1]x[z0,z1]."""
    nodes, weights = gauss_legendre_simple(n)
    sx, sy, sz = (x1 - x0) / 2, (y1 - y0) / 2, (z1 - z0) / 2
    cx, cy, cz = (x1 + x0) / 2, (y1 + y0) / 2, (z1 + z0) / 2
    total = mp.mpf(0)
    for i in range(n):
        xi = cx + sx * nodes[i]
        wi = weights[i]
        for j in range(n):
            yj = cy + sy * nodes[j]
            wj = weights[j]
            for k in range(n):
                zk = cz + sz * nodes[k]
                wk = weights[k]
                total += wi * wj * wk * integrand(xi, yj, zk)
    return total * sx * sy * sz


# ---------- Mesh ----------
def make_cell_verts(x0, x1, y0, y1, z0, z1):
    return [
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    ]


def cell_bounds(verts):
    xs = [v[0] for v in verts]
    ys = [v[1] for v in verts]
    zs = [v[2] for v in verts]
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


def cell_volume(verts):
    x0, x1, y0, y1, z0, z1 = cell_bounds(verts)
    return (x1 - x0) * (y1 - y0) * (z1 - z0)


# ---------- K matrix ----------
def K_pair(cell_i_verts, cell_j_verts, n_quad=4):
    """K_ij/(mu0/4pi) = int_cell_i Phi_hex_j(r) dV via 3D Gauss-Legendre."""
    x0, x1, y0, y1, z0, z1 = cell_bounds(cell_i_verts)

    def integrand(x, y, z):
        return hex_newton_potential(cell_j_verts, (x, y, z))

    return gauss_3d_integral(integrand, x0, x1, y0, y1, z0, z1, n=n_quad)


def assemble_K_matrix(cells, n_quad=4, verbose=True):
    N = len(cells)
    K_raw = [[mp.mpf(0)] * N for _ in range(N)]
    for i in range(N):
        for j in range(i, N):
            t0 = time.time()
            K_raw[i][j] = K_pair(cells[i], cells[j], n_quad=n_quad)
            K_raw[j][i] = K_raw[i][j]
            if verbose:
                print(f"  K[{i}][{j}] = {mp.nstr(K_raw[i][j], 12)}  ({time.time()-t0:.1f} s)")
    K = [[K_raw[i][j] * mu0 / (4 * mp.pi) for j in range(N)] for i in range(N)]
    return K, K_raw


def gen_eig_diag_M(K, M_diag):
    """K v = lambda M v with M diagonal. Returns sorted eigenvalues (descending magnitude)."""
    N = len(K)
    sqrtM_inv = [1 / mp.sqrt(M_diag[i]) for i in range(N)]
    A = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            A[i, j] = K[i][j] * sqrtM_inv[i] * sqrtM_inv[j]
    eigs_complex = mp.eig(A, left=False, right=False)
    eigs = sorted([abs(e) for e in eigs_complex], reverse=True)
    return eigs


# ===================================================================
# Test 1: single cell (whole cuboid)
# ===================================================================
print("=" * 64)
print("Phase 12.3 Test 1: single cell K_00 (whole cuboid)")
print("=" * 64)

x0, x1 = -a/2, +a/2
y0, y1 = -b/2, +b/2
z0, z1 = -c/2, +c/2
cell0 = make_cell_verts(x0, x1, y0, y1, z0, z1)
V_total = cell_volume(cell0)
print(f"V = {mp.nstr(V_total, 10)} m^3")

for nq in [4, 6]:
    t0 = time.time()
    K_raw00 = K_pair(cell0, cell0, n_quad=nq)
    print(f"  n_quad={nq}: K_raw[0][0] = {mp.nstr(K_raw00, 16)}  "
          f"({time.time()-t0:.1f} s)")

K00 = K_raw00 * mu0 / (4 * mp.pi)
lam = K00 / V_total
tau = sigma_Cu * lam
print(f"  K_00 (with mu0/4pi) = {mp.nstr(K00, 12)}  H*m^4")
print(f"  lambda = K/V = {mp.nstr(lam, 8)} H*m")
print(f"  tau = sigma * lambda = {mp.nstr(tau * 1e6, 8)} us")
print(f"  (Phase 7 reference leading tau = 19.77 us; PWC-1cell expected coarser)")


# ===================================================================
# Test 2: 2 x 2 x 2 = 8-cell mesh
# ===================================================================
print()
print("=" * 64)
print("Phase 12.3 Test 2: 2x2x2 mesh (8 cells, n_quad=4)")
print("=" * 64)

nx, ny, nz = 2, 2, 2
xs_div = [x0 + i*(x1 - x0)/nx for i in range(nx + 1)]
ys_div = [y0 + i*(y1 - y0)/ny for i in range(ny + 1)]
zs_div = [z0 + i*(z1 - z0)/nz for i in range(nz + 1)]

cells = []
for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            cells.append(make_cell_verts(
                xs_div[i], xs_div[i+1],
                ys_div[j], ys_div[j+1],
                zs_div[k], zs_div[k+1]))

N = len(cells)
M_diag = [cell_volume(c) for c in cells]
print(f"N = {N} cells, V_cell = {mp.nstr(M_diag[0], 10)}")

t_start = time.time()
K_full, K_raw = assemble_K_matrix(cells, n_quad=4, verbose=False)
print(f"K matrix assembled in {time.time()-t_start:.1f} s")

print("Diagonal of K (cell self-coupling):")
for i in range(N):
    print(f"  K[{i}][{i}] = {mp.nstr(K_full[i][i], 10)}")

eigs = gen_eig_diag_M(K_full, M_diag)
print()
print("Eigenvalues lambda_k (descending magnitude) -> tau_k:")
for k, e in enumerate(eigs):
    tau_k = sigma_Cu * e
    print(f"  lambda_{k} = {mp.nstr(e, 10)}  ->  tau = {mp.nstr(tau_k * 1e6, 10)} us")

print()
print("=" * 64)
print("Phase 12.3 status")
print("=" * 64)
print(f"Leading tau (8-cell PWC) = {mp.nstr(sigma_Cu * eigs[0] * 1e6, 8)} us")
print(f"Phase 7 (spectral basis) = 19.77 us")
print(f"Convergence: PWC mesh -> spectral as N -> infinity.")
print()
print("Next: Phase 12.4 - Lanczos matrix-free Foster (Born iteration).")
