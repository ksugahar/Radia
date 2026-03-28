"""
PMCHWT-SIBC: corrected EFIE using M = -Z_s*(n x J) magnetic current.

The standard EFIE: Z_s*J + jw*mu0*SL(J) = -jw*A_inc   (WRONG for finite Z_s)

The corrected EFIE adds the magnetic current contribution to A_scat:
  A_scat = mu0*SL(J) + (Z_s/jw)*curl_SL(n x J)

On the surface (exterior limit):
  curl_SL(F)_t = (1/2)*F_rotated + K_rotated(F)
  where F_rotated means tangential projection after curl.

Corrected system:
  [Z_s*M + jw*mu0*SL + Z_s*C, D^T] [J]   [-jw*A_inc]
  [D,                          0  ] [p] = [0        ]

where C = correction from magnetic current.
"""

import math
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi
R = 0.01
B0 = 0.001
freq = 1000
omega = 2 * math.pi * freq
H0 = B0 / MU_0
jw_mu0_R = omega * MU_0 * R

from ngsolve import (Mesh, HDivSurface, SurfaceL2, BilinearForm,
                     LinearForm, GridFunction, Integrate, CF, ds, div,
                     BND, TaskManager, x, y, specialcf, Cross)
from ngsolve.bem import LaplaceDL, LaplaceSL
from netgen.occ import Sphere, Pnt, OCCGeometry, Glue

# --- Mesh ---
sph = Sphere(Pnt(0, 0, 0), R)
for f in sph.faces:
    f.name = "sphere"
geo = OCCGeometry(Glue(sph.faces))
mesh = Mesh(geo.GenerateMesh(maxh=R / 3))
mesh.Curve(2)

fes = HDivSurface(mesh, order=0)
fes_l2 = SurfaceL2(mesh, order=0)
u, v = fes.TnT()
ndof = fes.ndof
n_p = fes_l2.ndof
ne = mesh.GetNE(BND)
print(f"Mesh: {ne} elems, {ndof} J-DOFs, {n_p} p-DOFs")

# --- Assemble operators ---
with TaskManager():
    DL_op = LaplaceDL(u.Trace() * ds) * v.Trace() * ds
    SL_op = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds

DL = np.zeros((ndof, ndof))
SL = np.zeros((ndof, ndof), dtype=complex)
for i in range(ndof):
    ei = GridFunction(fes)
    ei.vec[:] = 0
    ei.vec[i] = 1.0
    r1 = ei.vec.CreateVector()
    r1.data = DL_op.mat * ei.vec
    DL[:, i] = r1.FV().NumPy().copy()
    r2 = ei.vec.CreateVector()
    r2.data = SL_op.mat * ei.vec
    SL[:, i] = r2.FV().NumPy().copy()

# Mass
mass_bf = BilinearForm(fes)
mass_bf += u.Trace() * v.Trace() * ds
mass_bf.Assemble()
M = np.zeros((ndof, ndof))
rows, cols, vals = mass_bf.mat.COO()
for r, c, val in zip(rows, cols, vals):
    M[int(r), int(c)] = val

# Rotation R: <v, n x u>
n_cf = specialcf.normal(3)
rot_bf = BilinearForm(fes)
rot_bf += Cross(n_cf, u.Trace()) * v.Trace() * ds
rot_bf.Assemble()
R_rot = np.zeros((ndof, ndof))
rows, cols, vals = rot_bf.mat.COO()
for r, c, val in zip(rows, cols, vals):
    R_rot[int(r), int(c)] = val

# Verify: R^2 = -M (rotation by 180 = negation for tangential fields)
check = np.linalg.solve(M, R_rot @ np.linalg.solve(M, R_rot))
print(f"R*M^-1*R*M^-1 + I ~ 0? max = {np.max(np.abs(check + np.eye(ndof))):.4e}")

# Divergence
q = fes_l2.TestFunction()
u_J = fes.TrialFunction()
bf_D = BilinearForm(trialspace=fes, testspace=fes_l2)
bf_D += div(u_J.Trace()) * q * ds
bf_D.Assemble()
D_mat = np.zeros((n_p, ndof))
try:
    rows, cols, vals = bf_D.mat.COO()
    for r, c, val in zip(rows, cols, vals):
        D_mat[int(r), int(c)] = float(val)
except Exception:
    pass
rank_D = np.linalg.matrix_rank(D_mat)
D_red = D_mat[:rank_D, :]

# --- RHS ---
lf_A = LinearForm(fes)
lf_A += CF((-(B0 / 2) * y, (B0 / 2) * x, 0)) * v.Trace() * ds
lf_A.Assemble()
rhs_A_base = lf_A.vec.FV().NumPy().copy().astype(complex)

lf_H = LinearForm(fes)
lf_H += CF((y * H0 / R, -x * H0 / R, 0)) * v.Trace() * ds
lf_H.Assemble()
rhs_H = lf_H.vec.FV().NumPy().copy()

M_inv = np.linalg.inv(M)


def H_rms(J_vec):
    gf_re = GridFunction(fes)
    gf_im = GridFunction(fes)
    gf_re.vec.FV().NumPy()[:] = J_vec.real
    gf_im.vec.FV().NumPy()[:] = J_vec.imag
    ea = Integrate(CF(1), mesh, BND, element_wise=True)
    jr = Integrate(gf_re * gf_re, mesh, BND, element_wise=True)
    ji = Integrate(gf_im * gf_im, mesh, BND, element_wise=True)
    at = 0.0
    jt = 0.0
    for el in mesh.Elements(BND):
        at += abs(ea[el.nr])
        jt += abs(jr[el.nr]) + abs(ji[el.nr])
    return math.sqrt(jt / at)


# --- Correction matrix ---
# A_scat_t = mu0*SL(J)_t + (Z_s/jw)*curl_SL(nxJ)_t
# curl_SL(F)_t on exterior surface:
#   n x curl SL(F)|_ext = (1/2)*F + K(F)  [MFIE jump relation]
#   For tangential F: K(F) = -DL(F)
#   So n x curl SL(F)|_ext = (1/2)*F - DL(F)
#   curl SL(F)_t = [n x curl SL(F)] x n = ... rotation back
#
# In weak form with tangential test v:
#   <v, curl_SL(F)_t> = <v, [n x curl SL(F)] x n>
#   Using v . (a x n) = -(n x v) . a = -<n x v, a>:
#   Wait, actually for tangential v and a = n x curl SL:
#   v . (a x n) where a is tangential: a x n = -n x a, so v . (-n x a) = -(n x v) . a
#   Hmm this is not right.
#
# Simpler: tangential projection: F_t = F - n(n.F)
# For n x curl SL(F) which is tangential, rotating back:
# [n x curl SL(F)] x n gives the tangential field rotated 90deg back = curl SL(F)_t
# But [a_tang x n] is also tangential, rotated by -90deg.
# So: curl_SL(F)_t = [(1/2)*F - DL(F)] x n = n x [(1/2)*F - DL(F)] with sign...
# Actually a_tang x n = -n x a_tang. So:
# curl_SL(F)_t = -n x [(1/2)*F - DL(F)]
#
# In weak form: <v, curl_SL(F)_t> = -<v, n x [(1/2)*F - DL(F)]>
#              = -R_rot @ [(1/2)*I - M_inv @ DL] @ F_coeff
#              Wait, need to be careful with weak/strong.
#
# Let me compute directly: for coefficient vector c representing F = n x J:
# Step 1: F_coeff = M_inv @ R_rot @ c_J (coefficients of n x J in HDivSurface)
# Step 2: DL_F = DL @ F_coeff (weak form of DL applied to F)
# Step 3: half_minus_DL_F = 0.5*M @ F_coeff - DL_F (weak form of (1/2)F - DL(F))
# Step 4: curl_SL_F_t_weak = -R_rot @ M_inv @ half_minus_DL_F
#         = -R_rot @ M_inv @ (0.5*M - DL) @ F_coeff
#         = -R_rot @ (0.5*I - M_inv @ DL) @ M_inv @ R_rot @ c_J

# Correction matrix C (weak form): <v, (Z_s/jw)*curl_SL(n x J)_t>
# C_weak = (Z_s/jw) * [-R_rot @ (0.5*I - M_inv @ DL) @ M_inv @ R_rot]
# Note: the Z_s/jw factor is applied per ratio below.

# Pre-compute the Z_s-independent part:
half_minus_K = 0.5 * np.eye(ndof) - M_inv @ DL  # (1/2)*I - M^{-1}*DL = M^{-1}*(1/2*M - DL)
C_base = -R_rot @ half_minus_K @ M_inv @ R_rot  # correction operator (weak form, without Z_s/jw)

print()
print("PMCHWT-SIBC vs EFIE vs MFIE vs Analytical")
print(f"{'ratio':>8s} {'Ana':>10s} {'EFIE':>10s} {'PMCHWT':>10s} {'MFIE':>10s} "
      f"{'E/A':>6s} {'P/A':>6s}")
print("-" * 68)

jw = 1j * omega

for ratio in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    Zs = ratio * jw_mu0_R * (1 + 1j) / math.sqrt(2)
    beta = 1j * omega * MU_0 * R
    H_ana = abs((3.0 / 2.0) * jw * B0 * R / (beta + Zs)) * math.sqrt(2.0 / 3.0)

    rhs_A = -jw * rhs_A_base

    # Original EFIE
    n_t = ndof + rank_D
    A_efie = np.zeros((n_t, n_t), dtype=complex)
    A_efie[:ndof, :ndof] = Zs * M + jw * MU_0 * SL
    A_efie[:ndof, ndof:] = D_red.T
    A_efie[ndof:, :ndof] = D_red
    rhs_e = np.zeros(n_t, dtype=complex)
    rhs_e[:ndof] = rhs_A
    J_efie = np.linalg.solve(A_efie, rhs_e)[:ndof]

    # Corrected EFIE (PMCHWT-SIBC)
    # -jw*A_scat = -jw*mu0*SL*J + (-Z_s)*curl_SL(nxJ)_t
    # SIBC: Z_s*M*J + jw*mu0*SL*J + Z_s*C_base*J = -jw*A_inc (RHS same)
    # Wait: the sign. A_scat = mu0*SL + (Z_s/jw)*curl_SL(nxJ)
    # -jw*A_scat = -jw*mu0*SL - Z_s*curl_SL(nxJ)
    # SIBC: Z_s*J = -jw*A_inc - jw*A_scat = -jw*A_inc - jw*mu0*SL - Z_s*curl_SL(nxJ)
    # Z_s*M*J + jw*mu0*SL*J + Z_s*curl_SL_weak*J = -jw*A_inc
    # curl_SL_weak = C_base (precomputed without Z_s/jw, but including the sign)
    # Actually C_base already includes the negative sign from the derivation.
    # The additional term is: Z_s * C_base * J (not Z_s/jw * C_base)
    # Wait, let me re-derive:
    # A_M_t = (Z_s/jw)*curl_SL(nxJ)_t
    # -jw*A_M_t = -Z_s*curl_SL(nxJ)_t
    # In the EFIE, this adds: -Z_s*C_base (with C_base = curl_SL_weak without Z_s/jw)
    # Hmm, C_base has the -R*(1/2 - M^{-1}*DL)*M^{-1}*R factor
    # So the additional term in the system is: -Z_s*C_base or +Z_s*C_base?
    #
    # Let me verify with the sign: for the EFIE,
    # Z_s*M + jw*mu0*SL acts on J, equals -jw*A_inc
    # Adding -Z_s*curl_SL_weak to the LHS:
    # Z_s*M + jw*mu0*SL - Z_s*C_base = -jw*A_inc
    # Or adding +Z_s*C_base:
    # Z_s*M + jw*mu0*SL + Z_s*C_base = -jw*A_inc
    # Let me try both.

    A_p1 = np.zeros((n_t, n_t), dtype=complex)
    A_p1[:ndof, :ndof] = Zs * M + jw * MU_0 * SL + Zs * C_base
    A_p1[:ndof, ndof:] = D_red.T
    A_p1[ndof:, :ndof] = D_red
    rhs_p = np.zeros(n_t, dtype=complex)
    rhs_p[:ndof] = rhs_A
    try:
        J_p1 = np.linalg.solve(A_p1, rhs_p)[:ndof]
        H_p1 = H_rms(J_p1)
    except Exception:
        H_p1 = float('nan')

    A_p2 = np.zeros((n_t, n_t), dtype=complex)
    A_p2[:ndof, :ndof] = Zs * M + jw * MU_0 * SL - Zs * C_base
    A_p2[:ndof, ndof:] = D_red.T
    A_p2[ndof:, :ndof] = D_red
    try:
        J_p2 = np.linalg.solve(A_p2, rhs_p)[:ndof]
        H_p2 = H_rms(J_p2)
    except Exception:
        H_p2 = float('nan')

    re = H_rms(J_efie) / H_ana
    rp1 = H_p1 / H_ana
    rp2 = H_p2 / H_ana
    print(f"{ratio:8.3f} {H_ana:10.2f} {H_rms(J_efie):10.2f} "
          f"{H_p1:10.2f} {H_p2:10.2f} {re:6.3f} {rp1:6.3f} {rp2:6.3f}")

print()
print("Column 'PMCHWT' = +C_base, last column = -C_base")
