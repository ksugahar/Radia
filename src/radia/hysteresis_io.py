"""Hysteresis material data I/O and shape function construction.

Reads various input formats (.hys, .mat, B-H loops) and constructs
B-input Play model shape functions. Converts to energy-based parameters
for EnergyBasedHysteresis.

Functions:
    load_hys(filepath)              - JMAG .hys file reader
    load_mat(filepath)              - MATLAB .mat file reader (Potter-Schmulian)
    build_shape_functions(loops, K) - ShapeFunction.m port
    play_hysteron(f_k, Bx, By, eta, px, py) - PlayHysteron.m port (verification)
    convert_play_to_energy(eta, f_k_tables)  - B-input Play -> energy params
    hys_to_radia(filepath, K)       - One-step .hys -> Radia params
"""

import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid


# =====================================================================
# File readers
# =====================================================================

def load_hys(filepath):
    """Read JMAG .hys file.

    Format (multi-loop JMAG, as written by JilesAtherton.m):
        Line 1: model type ("Jiles Atherton")
        Line 2: type  maxpts  nloops
        Then, repeated nloops times (one self-contained block per loop):
            header: loopID  npts  flag
            units:  "A/m;tesla" or "tesla;A/m"
            npts rows: col1  col2

    Parameters
    ----------
    filepath : str or Path
        Path to .hys file.

    Returns
    -------
    loops : list of dict
        Each: {'B': ndarray, 'H': ndarray, 'Bmax': float}
        Descending branch (Bmax -> -Bmax).
    """
    filepath = Path(filepath)
    lines = filepath.read_text().splitlines()
    nloops = int(lines[1].split()[2])              # line 2: type  maxpts  nloops

    loops = []
    idx = 2
    for _ in range(nloops):
        # Each loop is a self-contained block: per-loop header + units + npts
        # data rows. JMAG (e.g. JilesAtherton.m) writes a header AND a units
        # line before EVERY loop; the previous single-header + np.loadtxt reader
        # only worked for nloops == 1 and crashed on multi-loop files
        # ("number of columns changed from 2 to 3" at the second loop header).
        hdr = lines[idx].split(); idx += 1
        npts = int(hdr[1])
        units = lines[idx].strip(); idx += 1
        b_first = units.lower().split(';')[0].startswith('t')
        block = np.array([[float(v) for v in lines[idx + j].split()[:2]]
                          for j in range(npts)])
        idx += npts
        if b_first:
            B, H = block[:, 0], block[:, 1]
        else:
            H, B = block[:, 0], block[:, 1]
        # Enforce the DESCENDING branch convention (Bmax -> -Bmax) promised by
        # this function's docstring and required by build_shape_functions; JMAG
        # stores each branch B-ascending, which otherwise yields garbage shape
        # functions (~1500 A/m reconstruction error, near-zero material).
        if B[0] < B[-1]:
            B, H = B[::-1], H[::-1]
        loops.append({
            'B': B.copy(),
            'H': H.copy(),
            'Bmax': float(np.max(np.abs(B)))
        })

    loops.sort(key=lambda x: x['Bmax'], reverse=True)
    return loops


def load_mat(filepath, bh_var='BH', db_var='dB', bmax_var='BMax'):
    """Read MATLAB .mat file (Potter-Schmulian B_input.mat format).

    Parameters
    ----------
    filepath : str or Path
        Path to .mat file.

    Returns
    -------
    loops : list of dict
        Each: {'B': ndarray, 'H': ndarray, 'Bmax': float}
    dB : float
        B step size.
    BMax : ndarray
        Array of Bmax values.
    """
    from scipy.io import loadmat
    data = loadmat(filepath, squeeze_me=True)

    BH_struct = data[bh_var]
    dB = float(data[db_var])
    BMax = np.asarray(data[bmax_var])

    loops = []
    for item in BH_struct:
        # Bug fix 2026-05-02: previous version did `item['B'].flat[0]` which
        # returns only the FIRST scalar element of the flat iterator, not the
        # whole branch. This produced HdataLUT with a single non-zero row per
        # column, breaking the JMAG identification at the input stage.
        B = np.asarray(item['B']).flatten()
        H = np.asarray(item['H']).flatten()
        loops.append({
            'B': B,
            'H': H,
            'Bmax': float(np.max(np.abs(B)))
        })

    loops.sort(key=lambda x: x['Bmax'], reverse=True)
    return loops, dB, BMax


# =====================================================================
# ShapeFunction.m port
# =====================================================================

def build_shape_functions(loops, dB=None):
    """Construct B-input Play shape functions from B-H loop data.

    Port of MATLAB ShapeFunction.m.

    Parameters
    ----------
    loops : list of dict
        B-H loops sorted by Bmax (largest first).
        Each: {'B': ndarray, 'H': ndarray, 'Bmax': float}
    dB : float, optional
        B step size. If None, estimated from data.

    Returns
    -------
    eta : ndarray, shape (Np,)
        B-space thresholds for play operators.
    f_k_tables : list of (r, f) tuples
        Shape function tables. r = |p| values, f = H contribution.
    Bplay : ndarray
        B-axis grid points for shape functions.
    """
    # Sort loops by Bmax descending (largest first).
    #
    # The ShapeFunction.m algorithm requires column 0 = LARGEST Bmax (longest
    # branch). The "last point correction" uses last_col = column1 - 2*ii - 1,
    # which decreases with ii — only consistent when ii=0 corresponds to the
    # longest branch.
    #
    # Bug fix 2026-05-02: prior version sorted ASCENDING, producing shape
    # functions that diverge from the input branches by ~1500 A/m in
    # reconstruction. Verified against bqm_model.compute_shape_functions
    # (independent port) and the MATLAB +bie/computeShapeFunctions reference.
    loops_desc = sorted(loops, key=lambda x: x['Bmax'], reverse=True)

    if dB is None:
        if len(loops_desc) >= 2:
            # loops_desc is sorted DESCENDING by Bmax, so [-1] is the smallest
            # and [-2] the next-smallest: take abs so the step is positive.
            # (The bare subtraction was negative -> np.zeros((negative, ...))
            #  crashed for any caller that did not pass dB explicitly.)
            dB = abs(loops_desc[-2]['Bmax'] - loops_desc[-1]['Bmax'])
        else:
            dB = loops_desc[0]['Bmax'] / 10.0

    n_loops = len(loops_desc)  # = column0 = row1
    Bmax = loops_desc[0]['Bmax']

    # Build HdataLUT: column 0 = LARGEST Bmax (longest branch), column k = k-th
    # largest. Padded with trailing zeros for shorter branches.
    max_n_B = int(round(2 * Bmax / dB)) + 1
    HdataLUT = np.zeros((max_n_B, n_loops))

    for col, loop in enumerate(loops_desc):
        n_pts = len(loop['B'])
        HdataLUT[:n_pts, col] = loop['H'][:n_pts]

    # --- ShapeFunction.m algorithm ---
    row0, column0 = HdataLUT.shape  # (max_n_B, n_loops)
    row1 = column0       # n_loops
    column1 = row0 - 1   # max_n_B - 1

    mu = np.zeros((row1, column1))

    # First difference
    for ii in range(row1):
        mu[ii, 0] = HdataLUT[0, ii] - HdataLUT[1, ii]

    # Interior differences (difference matrix D)
    for ii in range(row1 - 1):
        for jj in range(column1 - 2 * (ii + 1)):
            mu[ii, jj + 1] = (HdataLUT[jj + 1, ii] - HdataLUT[jj + 2, ii]
                              - (HdataLUT[jj, ii + 1] - HdataLUT[jj + 1, ii + 1]))

    # Row sum for last-point correction
    mu_temp = np.sum(mu, axis=1)

    # Last point correction
    for ii in range(row1):
        last_col = column1 - 2 * ii - 1  # 0-indexed: column1-2*(ii+1-1)-1 = column1-2*ii-1
        if last_col >= 0:
            mu[ii, last_col] = (HdataLUT[last_col, ii] - HdataLUT[last_col + 1, ii]
                                - mu_temp[ii])

    # Mirror and halve
    mu2 = np.zeros((row1 * 2, column1))
    ii0 = 0
    for jj in range(column1):
        for ii in range(row1 * 2 - ii0):
            mu2[ii + ii0, jj] = mu[(ii) // 2, jj] / 2  # ceil(ii/2) in 1-indexed = (ii)//2 in 0-indexed with +1 offset
        ii0 += 1

    # Flip upside-down
    mu2 = mu2[::-1, :]

    # Cumulative sum -> descending branch
    Hdown = -np.cumsum(mu2, axis=0)

    # Ascending branch (flip and negate)
    Hup = -Hdown[::-1, :]

    # Full symmetric shape function: [Hup; 0; Hdown]
    Hfunc = np.vstack([Hup, np.zeros((1, column1)), Hdown])

    # Build Bplay grid
    Bplay = np.arange(Bmax, -Bmax - dB / 4, -dB / 2)
    # Ensure length matches Hfunc rows
    n_func_rows = Hfunc.shape[0]  # should be 2*row1*2 + 1 = 4*n_loops + 1

    if len(Bplay) != n_func_rows:
        Bplay = np.linspace(Bmax, -Bmax, n_func_rows)

    # Create shape function interpolants
    f_k_tables = []
    Bplay_asc = Bplay[::-1]  # ascending for interpolation
    # Bug fix 2026-05-02: tolerance on the mask to keep the r=0 point even
    # when Bplay's middle element is something like -5e-17 due to arange
    # rounding. Without this, the (0, f_k(0)=0) anchor was dropped, causing
    # ~30 A/m residual reconstruction error from extrapolation near r=0.
    bp_tol = (dB / 2) * 1e-6
    for jj in range(column1):
        Hfunc_asc = Hfunc[::-1, jj]
        mask = Bplay_asc >= -bp_tol
        r_pos = np.maximum(Bplay_asc[mask], 0.0)   # snap tiny-negative to 0
        f_pos = Hfunc_asc[mask]
        # Force f_k(0) = 0 exactly (anti-symmetry); guard against accumulated
        # rounding in cumsum when r=0 is the first index.
        if r_pos[0] == 0:
            f_pos = f_pos.copy()
            f_pos[0] = 0.0
        f_k_tables.append((r_pos.copy(), f_pos.copy()))

    # Thresholds
    eta = np.arange(0, Bmax - dB / 4, dB / 2)
    Np = len(eta)

    # Trim to Np shape functions (eta and f_k should match)
    if len(f_k_tables) > Np:
        f_k_tables = f_k_tables[:Np]
    elif Np > len(f_k_tables):
        eta = eta[:len(f_k_tables)]

    return eta, f_k_tables, Bplay


# =====================================================================
# PlayHysteron.m port (for verification)
# =====================================================================

def play_hysteron(f_k_interps, Bx, By, eta, px, py):
    """B-input Play operator (vector). Port of PlayHysteron.m.

    Parameters
    ----------
    f_k_interps : list of callable
        Shape function interpolants f_k(|p|) -> H contribution.
    Bx, By : float
        Current B-field components.
    eta : ndarray, shape (Np,)
        Play thresholds.
    px, py : ndarray, shape (Np,)
        Current play operator states.

    Returns
    -------
    Hx, Hy : float
        Resulting H-field components.
    px, py : ndarray
        Updated play operator states.
    """
    Np = len(eta)
    px = px.copy()
    py = py.copy()

    # Update play states
    B_p = np.sqrt((Bx - px) ** 2 + (By - py) ** 2)
    denom = np.maximum(eta, B_p)
    denom[denom < 1e-30] = 1e-30
    px = Bx - eta * (Bx - px) / denom
    py = By - eta * (By - py) / denom

    # First operator (eta=0) always tracks B exactly
    px[0] = Bx
    py[0] = By

    # Evaluate shape functions
    p_mag = np.sqrt(px ** 2 + py ** 2)
    Hx = 0.0
    Hy = 0.0
    for k in range(Np):
        if p_mag[k] < 1e-30:
            continue
        fk_val = f_k_interps[k](p_mag[k])
        Hx += fk_val * px[k] / p_mag[k]
        Hy += fk_val * py[k] / p_mag[k]

    return Hx, Hy, px, py


# =====================================================================
# Play-to-Energy conversion
# =====================================================================

def convert_play_to_energy(eta, f_k_tables):
    """Convert B-input Play parameters to energy-based model parameters.

    Refactor 2026-05-XX: After consolidating Type 5 (energy hysteresis)
    as a thin subclass of Type 6 (Play model), the "chi" parameter of
    MatEnergyHysteresis is numerically identical to the Play threshold
    eta_k. The conversion is therefore a no-op pass-through; the
    function is retained for backward-compatible call sites.

    Parameters
    ----------
    eta : ndarray, shape (K,)
        B-input Play thresholds.
    f_k_tables : list of (r, f) tuples
        Shape function tables from build_shape_functions().

    Returns
    -------
    params : dict
        Keys: 'K', 'chi', 'f_k_tables'
        Ready for rad.MatEnergyHysteresis(**params).
    """
    K = len(f_k_tables)
    chi = eta.copy()  # numerically identical to eta in the unified Type 5/6 formulation

    return {
        'K': K,
        'chi': chi,
        'f_k_tables': f_k_tables,
    }


# =====================================================================
# Convenience: one-step file -> Radia params
# =====================================================================

def hys_to_radia(filepath, K=None, eps=1e-8):
    """Convert .hys file to energy-based hysteresis parameters.

    Returns a dict of parameters ready for rad.MatEnergyHysteresis().

    Parameters
    ----------
    filepath : str or Path
        Path to .hys file.
    K : int, optional
        Number of play operators. If None, determined from data.
    eps : float
        Regularization parameter.

    Returns
    -------
    params : dict
        Keys: 'K', 'chi', 'f_k_tables', 'eps'.
        Usage: rad.MatEnergyHysteresis(**params)
    """
    loops = load_hys(filepath)
    eta, f_k_tables, Bplay = build_shape_functions(loops)

    if K is not None and K < len(f_k_tables):
        indices = np.linspace(0, len(f_k_tables) - 1, K, dtype=int)
        f_k_tables = [f_k_tables[i] for i in indices]
        eta = eta[indices]

    p = convert_play_to_energy(eta, f_k_tables)
    return {'K': p['K'], 'chi': p['chi'],
            'f_k_tables': p['f_k_tables'], 'eps': eps}


def mat_to_radia(filepath, eps=1e-8):
    """Convert MATLAB .mat file (Potter-Schmulian) to energy-based params.

    Returns a dict of parameters ready for rad.MatEnergyHysteresis().

    Parameters
    ----------
    filepath : str or Path
        Path to B_input.mat.

    Returns
    -------
    params : dict
        Keys: 'K', 'chi', 'f_k_tables', 'eps'.
        Usage: rad.MatEnergyHysteresis(**params)
    """
    loops, dB, BMax = load_mat(filepath)
    eta, f_k_tables, Bplay = build_shape_functions(loops, dB=dB)

    p = convert_play_to_energy(eta, f_k_tables)
    return {'K': p['K'], 'chi': p['chi'],
            'f_k_tables': p['f_k_tables'], 'eps': eps}


# =====================================================================
# Convenience: one-step file -> Play model params
# =====================================================================

def hys_to_play_radia(filepath, K=None):
    """Convert .hys file to Play model parameters.

    Returns (K, eta, f_k_tables) ready for rad.MatPlayHysteresis().

    Parameters
    ----------
    filepath : str or Path
        Path to .hys file.
    K : int, optional
        Number of play operators. If None, determined from data.

    Returns
    -------
    K : int
        Number of play operators.
    eta : ndarray
        Play thresholds in Tesla.
    f_k_tables : list of (r_array, f_array) tuples
        Shape function tables.
    """
    loops = load_hys(filepath)
    eta, f_k_tables, Bplay = build_shape_functions(loops)

    if K is not None and K < len(f_k_tables):
        indices = np.linspace(0, len(f_k_tables) - 1, K, dtype=int)
        f_k_tables = [f_k_tables[i] for i in indices]
        eta = eta[indices]

    return len(f_k_tables), eta, f_k_tables


def mat_to_play_radia(filepath):
    """Convert MATLAB .mat file (Potter-Schmulian) to Play model parameters.

    Returns (K, eta, f_k_tables) ready for rad.MatPlayHysteresis().

    Parameters
    ----------
    filepath : str or Path
        Path to B_input.mat.

    Returns
    -------
    K : int
        Number of play operators.
    eta : ndarray
        Play thresholds in Tesla.
    f_k_tables : list of (r_array, f_array) tuples
        Shape function tables.
    """
    loops, dB, BMax = load_mat(filepath)
    eta, f_k_tables, Bplay = build_shape_functions(loops, dB=dB)
    return len(f_k_tables), eta, f_k_tables
