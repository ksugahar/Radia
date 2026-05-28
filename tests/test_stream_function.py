"""Tests for radia.stream_function -- generic (ACA+)+TSVD least-norm solver.

The C++ core (src/core/rad_stream_function.cpp) is kernel-agnostic: it does
(ACA+)+TSVD of an M x N matrix whose entries A(i,j) come from a caller-supplied
callback.  ACA+ itself is HACApK's cHACApK_acaplus (single source of truth);
only the TSVD recompression (manuscript Method 2/3, IEEJ SA-25-020) lives in
Radia.

Two kernels are exercised:
  1. a coil Biot-Savart kernel (mirrored rectangular current loops) implemented
     in numpy HERE, to cross-check bit-for-bit against the validated Fortran
     reference coil_solver.f90 (method_aca_tsvd_1/2);
  2. Radia's OWN field computation (radia.Fld over permanent-magnet objects via
     radia_field_kernel), to prove the same machinery serves magnetic materials,
     not just coils.

Method mapping (radia <-> f90):
  radia method=3 (E = diag(Sc) Vc' D, 2 SVDs)  <->  f90 method_aca_tsvd_2
  radia method=2 (SVD(C)+SVD(D)+Middle)         <->  f90 method_aca_tsvd_1
"""
import os
import sys

import numpy as np
import pytest

from radia.stream_function import (
    aca_tsvd, pseudo_inverse_solve, radia_field_kernel,
)

# f2py reference (LAB-only network drive).
REF_DIR = r"W:\04_卒論論文関係\2025年度\046_伊藤海人\2026_01_06_f2py_matlab比較\f2py"


# --------------------------------------------------------------------------
# Coil Biot-Savart kernel (mirrored rectangular loops) -- numpy reference that
# replicates coil_solver.f90's biot_savart_entry exactly, used to drive the
# f90 cross-check.  This lives in the TEST only; the production module embeds
# no field kernel.
# --------------------------------------------------------------------------
def _seg_Hz(O, P1, P2):
    eps = 1.0e-15
    d = P2 - P1
    R21 = float(d @ d)
    o1 = O - P1
    o2 = O - P2
    OR1 = np.sqrt(o1 @ o1)
    OR2 = np.sqrt(o2 @ o2)
    O121 = float(o1 @ d)
    O221 = float(o2 @ d)
    Rc12 = O121 / OR1
    Rc21 = -O221 / OR2
    L = o1 - O121 * d / R21
    L1 = float(L @ L) + eps
    f = (Rc12 + Rc21) / (4.0 * np.pi * L1 * R21)
    return (d[0] * L[1] - d[1] * L[0]) * f  # Hz


def _loop_Hz(obs, center, offsets):
    obs = np.asarray(obs, float)
    C = np.asarray(center, float) + np.asarray(offsets, float)  # (4,3) upper corners
    hz = 0.0
    for c in range(4):
        n = (c + 1) % 4
        hz += _seg_Hz(obs, C[c], C[n])
    Cl = C.copy()
    Cl[:, 2] = -Cl[:, 2]  # mirror through z
    for c in range(4):
        n = (c + 1) % 4
        hz += _seg_Hz(obs, Cl[c], Cl[n])
    return hz


def _square_loop(a):
    return np.array([[-a / 2, -a / 2, 0.0],
                     [ a / 2, -a / 2, 0.0],
                     [ a / 2,  a / 2, 0.0],
                     [-a / 2,  a / 2, 0.0]])


def _coil_geometry(obs_z, n_side=10, m_side=5, z0=0.005, a=0.02,
                   span=0.1, obs_span=0.08):
    xs = np.linspace(-span, span, n_side)
    ys = np.linspace(-span, span, n_side)
    cx, cy = np.meshgrid(xs, ys, indexing="ij")
    centers = np.column_stack([cx.ravel(), cy.ravel(), np.full(cx.size, z0)])
    n = centers.shape[0]
    offsets = np.broadcast_to(_square_loop(a), (n, 4, 3)).copy()

    oxs = np.linspace(-obs_span, obs_span, m_side)
    oys = np.linspace(-obs_span, obs_span, m_side)
    ox, oy = np.meshgrid(oxs, oys, indexing="ij")
    obs = np.column_stack([ox.ravel(), oy.ravel(), np.full(ox.size, obs_z)])
    return obs, centers, offsets


def _coil_dense_A(obs, centers, offsets):
    M, N = obs.shape[0], centers.shape[0]
    A = np.empty((M, N))
    for i in range(M):
        for j in range(N):
            A[i, j] = _loop_Hz(obs[i], centers[j], offsets[j])
    return A


def _recon_err(A, res, k=None):
    if k is None:
        k = res.k_aca
    Alr = res.U[:, :k] @ np.diag(res.S[:k]) @ res.V[:, :k].T
    return np.linalg.norm(A - Alr) / np.linalg.norm(A)


# --------------------------------------------------------------------------
# Reconstruction vs true dense A (coil kernel)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("method", [3, 2])
def test_reconstruction_near_field_full_rank(method):
    """Near field is effectively full rank (k_aca = min(M,N)); the
    (ACA+)+TSVD factorization reconstructs A to machine precision."""
    obs, centers, offsets = _coil_geometry(obs_z=0.05)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    entry = lambda i, j: A[i, j]
    res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                   aca_eps=1.0e-8, method=method)
    assert res.k_aca == min(M, N)
    assert _recon_err(A, res) < 1.0e-10
    assert np.all(res.S >= 0.0)
    assert np.all(np.diff(res.S[:res.k_aca]) <= 1.0e-12)


def test_low_rank_far_field_compresses():
    """Far field is smooth -> ACA+ compresses (k_aca << min(M,N))."""
    obs, centers, offsets = _coil_geometry(obs_z=0.5)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    entry = lambda i, j: A[i, j]
    res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                   aca_eps=1.0e-8, method=3)
    assert res.k_aca < min(M, N)        # compression happened
    assert res.k_aca >= 1
    assert _recon_err(A, res) < 1.0e-4


def test_truncation_monotonic():
    """Reconstruction error decreases (weakly) as more TSVD modes are kept."""
    obs, centers, offsets = _coil_geometry(obs_z=0.05)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    entry = lambda i, j: A[i, j]
    res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                   aca_eps=1.0e-10, method=3)
    errs = [_recon_err(A, res, k=k) for k in (5, 10, 15, res.k_aca)]
    for a, b in zip(errs, errs[1:]):
        assert b <= a + 1.0e-12


def test_methods_2_and_3_agree():
    """Method 2 and Method 3 yield the same singular spectrum."""
    obs, centers, offsets = _coil_geometry(obs_z=0.1)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    entry = lambda i, j: A[i, j]
    r3 = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                  aca_eps=1.0e-8, method=3)
    r2 = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                  aca_eps=1.0e-8, method=2)
    assert r2.k_aca == r3.k_aca
    ns = r3.k_aca
    rel = np.linalg.norm(r2.S[:ns] - r3.S[:ns]) / np.linalg.norm(r3.S[:ns])
    assert rel < 1.0e-9


# --------------------------------------------------------------------------
# Pseudo-inverse solve
# --------------------------------------------------------------------------
def test_pseudo_inverse_solve_recovers_range():
    """For a target B = A @ phi0 in range(A), the least-norm solve reproduces
    B (near field is full row rank, so A phi == B exactly)."""
    obs, centers, offsets = _coil_geometry(obs_z=0.05)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    entry = lambda i, j: A[i, j]
    res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                   aca_eps=1.0e-10, method=3)
    rng = np.random.default_rng(0)
    phi0 = rng.standard_normal(N)
    B = A @ phi0
    phi = pseudo_inverse_solve(res, B, k_mode=res.k_aca)
    assert phi.shape == (N,)
    assert np.linalg.norm(A @ phi - B) / np.linalg.norm(B) < 1.0e-8
    assert np.linalg.norm(phi) <= np.linalg.norm(phi0) + 1.0e-9


def test_pseudo_inverse_solve_validates_B_length():
    obs, centers, offsets = _coil_geometry(obs_z=0.05)
    M, N = obs.shape[0], centers.shape[0]
    A = _coil_dense_A(obs, centers, offsets)
    res = aca_tsvd(M, N, lambda i, j: A[i, j], modes=min(M, N),
                   kmax=min(M, N), aca_eps=1.0e-8, method=3)
    with pytest.raises(ValueError):
        pseudo_inverse_solve(res, np.zeros(M + 1))


def test_aca_tsvd_validates_args():
    with pytest.raises(ValueError):
        aca_tsvd(0, 5, lambda i, j: 0.0, modes=1, kmax=1)
    with pytest.raises(TypeError):
        aca_tsvd(5, 5, 123, modes=1, kmax=1)  # not callable


# --------------------------------------------------------------------------
# Generic path through Radia's OWN field computation (magnetic materials)
# --------------------------------------------------------------------------
def test_radia_field_kernel_magnets():
    """The same (ACA+)+TSVD machinery works with Radia's MMM/MSC field over
    permanent-magnet objects (not coils), via radia_field_kernel + radia.Fld.
    Proves the solver is kernel-agnostic and reuses Radia's existing kernels."""
    import radia as rad
    rad.UtiDelAll()
    try:
        # N small permanent-magnet cubes on a plane, M observation points above.
        xs = np.linspace(-0.04, 0.04, 6)
        ys = np.linspace(-0.04, 0.04, 6)
        gx, gy = np.meshgrid(xs, ys, indexing="ij")
        centers = np.column_stack([gx.ravel(), gy.ravel(), np.zeros(gx.size)])
        N = centers.shape[0]
        sources = [rad.ObjRecMag(c.tolist(), [0.01, 0.01, 0.01], [0.0, 0.0, 954930.0])
                   for c in centers]

        oxs = np.linspace(-0.03, 0.03, 4)
        oys = np.linspace(-0.03, 0.03, 4)
        ox, oy = np.meshgrid(oxs, oys, indexing="ij")
        obs = np.column_stack([ox.ravel(), oy.ravel(), np.full(ox.size, 0.06)])
        M = obs.shape[0]

        entry = radia_field_kernel(obs, sources, component=2, field="b")
        A = np.array([[entry(i, j) for j in range(N)] for i in range(M)])

        res = aca_tsvd(M, N, entry, modes=min(M, N), kmax=min(M, N),
                       aca_eps=1.0e-8, method=3)
        assert 1 <= res.k_aca <= min(M, N)
        assert _recon_err(A, res) < 1.0e-5
    finally:
        rad.UtiDelAll()


# --------------------------------------------------------------------------
# f90 cross-check (LAB only) -- runs in a fresh subprocess to dodge conftest's
# DLL-search pollution (the f2py module bundles its own Intel/MKL DLLs and
# cannot be imported alongside conftest's preloaded ngsolve/cubit).
# Verified bit-exact: k_aca identical, ||S_f90 - S_radia||/||S_f90|| ~ 1e-15.
# --------------------------------------------------------------------------
_F90_SUBPROCESS = r"""
import os, sys
import numpy as np
import radia  # noqa: F401  (loads radia + ngsolve MKL first, as in production)
from radia.stream_function import aca_tsvd

ref = sys.argv[1]
try:
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(ref)
    if ref not in sys.path:
        sys.path.insert(0, ref)
    import coil_solver
except Exception as exc:  # noqa: BLE001
    print("SKIP", type(exc).__name__)
    sys.exit(0)

def seg_Hz(O, P1, P2):
    eps = 1.0e-15
    d = P2 - P1; R21 = float(d @ d)
    o1 = O - P1; o2 = O - P2
    OR1 = np.sqrt(o1 @ o1); OR2 = np.sqrt(o2 @ o2)
    O121 = float(o1 @ d); O221 = float(o2 @ d)
    Rc12 = O121 / OR1; Rc21 = -O221 / OR2
    L = o1 - O121 * d / R21
    L1 = float(L @ L) + eps
    f = (Rc12 + Rc21) / (4.0 * np.pi * L1 * R21)
    return (d[0] * L[1] - d[1] * L[0]) * f

def loop_Hz(obs, center, offsets):
    C = center + offsets
    hz = 0.0
    for c in range(4):
        n = (c + 1) % 4; hz += seg_Hz(obs, C[c], C[n])
    Cl = C.copy(); Cl[:, 2] = -Cl[:, 2]
    for c in range(4):
        n = (c + 1) % 4; hz += seg_Hz(obs, Cl[c], Cl[n])
    return hz

n, m, z0, a, span, ospan, obs_z = 10, 5, 0.005, 0.02, 0.1, 0.08, 0.1
xs = np.linspace(-span, span, n); ys = np.linspace(-span, span, n)
cx, cy = np.meshgrid(xs, ys, indexing="ij")
centers = np.column_stack([cx.ravel(), cy.ravel(), np.full(cx.size, z0)])
N = centers.shape[0]
loop = np.array([[-a/2,-a/2,0.],[a/2,-a/2,0.],[a/2,a/2,0.],[-a/2,a/2,0.]])
off = np.broadcast_to(loop, (N, 4, 3)).copy()
oxs = np.linspace(-ospan, ospan, m); oys = np.linspace(-ospan, ospan, m)
ox, oy = np.meshgrid(oxs, oys, indexing="ij")
obs = np.column_stack([ox.ravel(), oy.ravel(), np.full(ox.size, obs_z)])
M = obs.shape[0]; kmax = min(M, N); modes = kmax; eps = 1.0e-8

A = np.empty((M, N))
for i in range(M):
    for j in range(N):
        A[i, j] = loop_Hz(obs[i], centers[j], off[j])
entry = lambda i, j: A[i, j]

fails = []
for f90name, method in (("method_aca_tsvd_2", 3), ("method_aca_tsvd_1", 2)):
    r = aca_tsvd(M, N, entry, modes=modes, kmax=kmax, aca_eps=eps, method=method)
    fn = getattr(coil_solver, f90name)
    Uf, Sf, Vf, kf, _t1, _t2, _pk = fn(
        np.ascontiguousarray(obs[:, 0]), np.ascontiguousarray(obs[:, 1]),
        np.ascontiguousarray(obs[:, 2]), np.ascontiguousarray(centers[:, 0]),
        np.ascontiguousarray(centers[:, 1]), np.ascontiguousarray(centers[:, 2]),
        np.ascontiguousarray(off[:, :, 0]), np.ascontiguousarray(off[:, :, 1]),
        np.ascontiguousarray(off[:, :, 2]), modes, kmax, eps)
    ns = min(int(kf), r.k_aca)
    rel = np.linalg.norm(Sf[:ns] - r.S[:ns]) / np.linalg.norm(Sf[:ns])
    if int(kf) != r.k_aca or rel > 1.0e-9:
        fails.append(f"{f90name}: kf={kf} k_aca={r.k_aca} relS={rel:.2e}")

if fails:
    print("FAIL", "; ".join(fails)); sys.exit(1)
print("PASS"); sys.exit(0)
"""


@pytest.mark.skipif(not os.path.isdir(REF_DIR),
                    reason="f2py coil_solver reference (W: drive) not available")
def test_matches_f90_reference():
    """radia C++ port matches the f2py coil_solver.f90 reference exactly
    (same k_aca, same singular values) for BOTH method pairs, using a coil
    Biot-Savart kernel.  Runs in a fresh subprocess to dodge conftest's
    DLL-search pollution."""
    import subprocess
    proc = subprocess.run([sys.executable, "-c", _F90_SUBPROCESS, REF_DIR],
                          capture_output=True, text=True, timeout=300)
    out = (proc.stdout + "\n" + proc.stderr).strip()
    if proc.stdout.startswith("SKIP"):
        pytest.skip(f"f2py coil_solver could not be imported: {proc.stdout.strip()}")
    assert "PASS" in proc.stdout, f"f90 cross-check failed (rc={proc.returncode}):\n{out}"
