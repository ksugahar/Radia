"""
hdiv_vim_structured.py -- HAND-BUILT structured-hex RT0 design spec for the C++ VIM-in-core port.

This mirrors EXACTLY what the C++ `rad_hdiv_vim` will do: enumerate faces/cells of a structured
nx*ny*nz hex grid, build the charge map B (M -> (rho=-div M per cell, sigma=M.n per bnd face)) with
deterministic signs (global face normals +x/+y/+z), build a dense Coulomb Gram G, and form the
symmetric demag operator N = B^T G B.  NO NGSolve -- the topology is enumerated by hand, the way
C++ must.  Validates the structural foundation the C++ first-increment golden test will lock:
  ndof (faces), n_loop = ndof - rank(B), ||N - N^T||, max ||N . loop||/||N||  (loops = ker B).
Cross-checks against the NGSolve prototype golden (hdiv_demag_quad_self.json): regular 3x3x3 ->
ndof=108, n_loop=28, asym~3.9e-16, loop_res~4.8e-16.

Loop-nullity is AUTOMATIC by construction (B.loop=0 => N.loop = B^T G (B loop) = 0 for ANY G), so
this increment validates (a) the topology/charge-map B is correct (-> the right n_loop), (b) G is
assembled symmetric.  The physics of G (spectrum) is validated later with the accurate self-energy.
"""
import json, os
import numpy as np
from math import pi

HERE = os.path.dirname(os.path.abspath(__file__))
c_cube = 1.88231  # unit-cube mean reciprocal distance (self-energy constant); see quad_self

def build_structured_rt0(nx, ny, nz, h=1.0):
    """Enumerate a structured hex grid's RT0 faces. Returns face table + cell table.
    Global face normals are +x / +y / +z.  Face index layout: all x-faces, then y, then z."""
    nodes = lambda i, j, k: np.array([i*h, j*h, k*h])
    faces = []   # each: dict(normal_axis, area, centroid(3), cell_lo, cell_hi, is_bnd)
    def cell_id(i, j, k): return (i*ny + j)*nz + k       # cell (i,j,k), i in[0,nx) etc.
    # x-faces: i in [0,nx], j in [0,ny), k in [0,nz); normal +x; area h*h; separates (i-1,*) | (i,*)
    fx0 = 0
    for i in range(nx+1):
        for j in range(ny):
            for k in range(nz):
                c = nodes(i, j, k) + np.array([0, 0.5*h, 0.5*h])
                lo = cell_id(i-1, j, k) if i > 0 else -1
                hi = cell_id(i,   j, k) if i < nx else -1
                faces.append(dict(ax=0, area=h*h, c=c, lo=lo, hi=hi, bnd=(i == 0 or i == nx)))
    fy0 = len(faces)
    for i in range(nx):
        for j in range(ny+1):
            for k in range(nz):
                c = nodes(i, j, k) + np.array([0.5*h, 0, 0.5*h])
                lo = cell_id(i, j-1, k) if j > 0 else -1
                hi = cell_id(i, j,   k) if j < ny else -1
                faces.append(dict(ax=1, area=h*h, c=c, lo=lo, hi=hi, bnd=(j == 0 or j == ny)))
    fz0 = len(faces)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz+1):
                c = nodes(i, j, k) + np.array([0.5*h, 0.5*h, 0])
                lo = cell_id(i, j, k-1) if k > 0 else -1
                hi = cell_id(i, j, k)   if k < nz else -1
                faces.append(dict(ax=2, area=h*h, c=c, lo=lo, hi=hi, bnd=(k == 0 or k == nz)))
    n_cell = nx*ny*nz
    cell_c = np.zeros((n_cell, 3)); cell_V = np.full(n_cell, h**3)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                cell_c[cell_id(i, j, k)] = nodes(i, j, k) + 0.5*h
    return faces, n_cell, cell_c, cell_V

def charge_map(faces, n_cell, cell_V):
    """B: charge DOFs (n_cell rho + n_bnd sigma) x n_face.  Matches the NGSolve prototype
    normalization Bv/V (rho per unit vol), Bb/area (sigma per unit area).
    rho_cell = -div M = -(sum over the cell's faces of  s * flux),  s=+1 if the global face normal
    points OUT of the cell (it's the cell's HI face), -1 if IN (cell's LO face).  flux dof = M_f."""
    nf = len(faces)
    bnd_faces = [f for f in range(nf) if faces[f]["bnd"]]
    bnd_row = {f: n_cell + r for r, f in enumerate(bnd_faces)}   # charge row for each bnd face
    n_charge = n_cell + len(bnd_faces)
    B = np.zeros((n_charge, nf))
    for f, fc in enumerate(faces):
        # cell on the LO side: this face is that cell's HI face -> global normal points OUT -> +1 div
        if fc["lo"] >= 0:
            B[fc["lo"], f] += -(+1.0) / cell_V[fc["lo"]]     # rho = -div, contribution -(+1)
        # cell on the HI side: this face is that cell's LO face -> global normal points IN -> -1 div
        if fc["hi"] >= 0:
            B[fc["hi"], f] += -(-1.0) / cell_V[fc["hi"]]     # -(-1) = +1
        if fc["bnd"]:
            # sigma = M . n_OUTWARD.  Global face normal is +axis; outward (out of the domain) is
            # +global if the cell sits on the LO side (this is the domain's HIGH boundary) and
            # -global if the cell sits on the HI side (domain's LOW boundary).  Using the global
            # normal for all boundary faces (the earlier bug) flips sigma on the low boundary ->
            # a spurious monopole surface charge -> unphysical demag factors (>1).
            out_sign = 1.0 if fc["lo"] >= 0 else -1.0
            B[bnd_row[f], f] += out_sign / fc["area"]
    return B, n_charge, len(bnd_faces)

def coulomb_gram(faces, n_cell, cell_c, cell_V):
    """dense charge-charge Coulomb Gram G (symmetric).  Charge cells = volume cells then bnd faces.
    off-diag: centroid monopole meas_a meas_b/(4pi r);  diag: cube/square self-energy."""
    bnd_faces = [f for f in range(len(faces)) if faces[f]["bnd"]]
    cent = np.vstack([cell_c, np.array([faces[f]["c"] for f in bnd_faces])]) if bnd_faces else cell_c
    meas = np.concatenate([cell_V, np.array([faces[f]["area"] for f in bnd_faces])])
    D = np.linalg.norm(cent[:, None, :] - cent[None, :, :], axis=2); np.fill_diagonal(D, np.inf)
    G = (meas[:, None]*meas[None, :]) / (4*pi*D)
    diag = np.empty(len(meas))
    diag[:n_cell] = c_cube * cell_V**(5/3.) / (4*pi)
    diag[n_cell:] = 2.97321 * np.array([faces[f]["area"] for f in bnd_faces])**1.5 / (4*pi)
    np.fill_diagonal(G, diag)
    return G

def cell_faces(faces, n_cell):
    """per cell -> {axis: [lo_face, hi_face]} (the two faces normal to each axis)."""
    cf = [{0: [None, None], 1: [None, None], 2: [None, None]} for _ in range(n_cell)]
    for fi, f in enumerate(faces):
        a = f["ax"]
        if f["lo"] >= 0: cf[f["lo"]][a][1] = fi   # cell is on the LO side -> face is its HI face
        if f["hi"] >= 0: cf[f["hi"]][a][0] = fi   # cell is on the HI side -> face is its LO face
    return cf

def assemble_mass(faces, n_cell, h=1.0):
    """Lowest-order RT0 HDiv mass matrix (regular hex grid), dense (nf x nf).  Unit-flux basis ->
    per-cell per-axis 2x2 block (1/h)[[1/3,1/6],[1/6,1/3]] on (lo_face, hi_face); shared faces
    accumulate from both adjacent cells (interior diag 2/(3h), boundary 1/(3h))."""
    nf = len(faces); M = np.zeros((nf, nf))
    cf = cell_faces(faces, n_cell)
    blk = (1.0/h) * np.array([[1/3., 1/6.], [1/6., 1/3.]])
    for c in range(n_cell):
        for a in (0, 1, 2):
            idx = cf[c][a]
            for p in range(2):
                for q in range(2):
                    M[idx[p], idx[q]] += blk[p, q]
    return M

def run(nx, ny, nz, tag):
    faces, n_cell, cell_c, cell_V = build_structured_rt0(nx, ny, nz)
    ndof = len(faces)
    B, n_charge, n_bnd = charge_map(faces, n_cell, cell_V)
    G = coulomb_gram(faces, n_cell, cell_c, cell_V)
    N = B.T @ G @ B
    # loops = ker(B); rank via SVD
    U, S, Vt = np.linalg.svd(B)
    rankB = int(np.sum(S > 1e-9*S.max())); n_loop = ndof - rankB
    loops = Vt[rankB:, :]
    Nn = np.linalg.norm(N, 2)
    asym = np.linalg.norm(N - N.T) / Nn
    loop_res = max(np.linalg.norm(N @ loops[k]) for k in range(n_loop)) / Nn if n_loop else 0.0
    print(f"[{tag}] {nx}x{ny}x{nz}: ndof={ndof} (x{nx+1 if False else ''}) n_cell={n_cell} "
          f"n_bnd={n_bnd} n_charge={n_charge} rankB={rankB} n_loop={n_loop}")
    print(f"    ||N-N^T||/||N|| = {asym:.2e}   max||N.loop||/||N|| = {loop_res:.2e}")
    return dict(tag=tag, ndof=ndof, n_cell=n_cell, n_bnd=n_bnd, n_charge=n_charge,
                rankB=rankB, n_loop=n_loop, asym=float(asym), loop_res=float(loop_res))

if __name__ == "__main__":
    print("=== hand-built structured-hex RT0 (C++ design spec) -- structural validation ===")
    res = [run(1, 1, 1, "1x1x1"), run(2, 2, 2, "2x2x2"), run(3, 3, 3, "3x3x3")]
    print("\nGOLDEN cross-check (NGSolve prototype, hdiv_demag_quad_self.json): 3x3x3 -> ndof=108, "
          "n_loop=28, asym~3.9e-16, loop_res~4.8e-16")
    r3 = res[-1]
    ok = (r3["ndof"] == 108 and r3["n_loop"] == 28 and r3["asym"] < 1e-12 and r3["loop_res"] < 1e-10)
    print(f"MATCH: ndof {r3['ndof']}==108 {r3['ndof']==108}, n_loop {r3['n_loop']}==28 {r3['n_loop']==28}, "
          f"asym<1e-12 {r3['asym']<1e-12}, loop_res<1e-10 {r3['loop_res']<1e-10}  => {'PASS' if ok else 'FAIL'}")
    with open(os.path.join(HERE, "hdiv_vim_structured.json"), "w") as f:
        json.dump({"cases": res, "golden_match_3x3x3": bool(ok)}, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_vim_structured.json"))
