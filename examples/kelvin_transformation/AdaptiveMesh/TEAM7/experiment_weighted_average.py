"""
TEAM7: A法とOmega法の重み付き平均の収束性検討

渦電流問題（複素数）での重み付き平均の挙動を調べる。
真値が不明なので、メッシュ refinement での収束性を観察。

複素数の場合、誤差ノルムの定義に注意が必要：
  ||f||^2 = Re(∫ f* · f dV) = ∫ |f|^2 dV
"""
import os
import sys

os.environ['NETGENDIR'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin'
os.environ['PATH'] = r'S:\NGSolve\01_GitHub\install_ngsolve\bin;' + os.environ.get('PATH', '')
sys.path.insert(0, r'S:\NGSolve\01_GitHub\install_ngsolve\Lib\site-packages')

from numpy import pi, sqrt
import numpy as np
from ngsolve import *
from ngsolve.krylovspace import CGSolver

from team7_geometry import (
    create_team7_geometry,
    get_coil_current_density,
    get_material_properties,
    mu0, sigma_al
)

print("=" * 70)
print("TEAM7: Weighted Average Convergence Study")
print("=" * 70)

freq = 50.0  # Hz
omega = 2 * pi * freq

# Test measurement point (inside the plate)
test_x = 72e-3
test_y = 72e-3
test_z = 9.5e-3  # Middle of plate (z = 0 to 19mm)

print(f"\nFrequency: {freq} Hz")
print(f"Test point: ({test_x*1000:.0f}, {test_y*1000:.0f}, {test_z*1000:.0f}) mm")

# Test different mesh sizes
maxh_list = [0.1, 0.08, 0.06]
order = 1

results = []

for maxh in maxh_list:
    print(f"\n{'='*70}")
    print(f"maxh = {maxh*1000:.0f} mm")
    print("=" * 70)

    mesh = create_team7_geometry(maxh=maxh)
    mesh.Curve(order)

    mu, sigma, inv_mu = get_material_properties(mesh)
    J0 = get_coil_current_density(mesh)

    # ============================================================
    # A method (時間調和渦電流問題)
    # curl(1/mu * curl(A)) + j*omega*sigma*A = J0
    # ============================================================
    print("\n--- A method ---")

    fes_A = HCurl(mesh, order=order, dirichlet="outer", complex=True)
    A, v_A = fes_A.TnT()
    print(f"  DOFs: {fes_A.ndof}")

    a_A = BilinearForm(fes_A)
    a_A += inv_mu * InnerProduct(curl(A), curl(v_A)) * dx
    a_A += 1j * omega * sigma * InnerProduct(A, v_A) * dx
    a_A += 1e-12 * InnerProduct(A, v_A) * dx  # regularization
    a_A.Assemble()

    f_A = LinearForm(fes_A)
    f_A += InnerProduct(J0, v_A) * dx
    f_A.Assemble()

    gf_A = GridFunction(fes_A)

    try:
        gf_A.vec.data = a_A.mat.Inverse(fes_A.FreeDofs(), inverse="umfpack") * f_A.vec
        print("  Solved with UMFPACK")
    except:
        from ngsolve.krylovspace import GMResSolver
        pre = a_A.mat.CreateSmoother(fes_A.FreeDofs())
        solver = GMResSolver(a_A.mat, pre, maxiter=3000, tol=1e-8, printrates=False)
        gf_A.vec.data = solver * f_A.vec
        print("  Solved with GMRes")

    B_A = curl(gf_A)
    H_A = inv_mu * B_A

    # ============================================================
    # Omega method (時間調和、導体領域のみ)
    # 静磁場Omega法は渦電流問題には直接使えないが、
    # 比較のためにB_Aから求めたHとの差を見る
    # ============================================================
    # 渦電流問題ではT-Omega法を使う必要があるが、
    # ここでは単純にA法のB_AからH_AとH1投影を比較

    print("\n--- H1 Recovery from A method ---")

    # H_AをH1空間に投影（curl-free成分を抽出）
    fes_rec = H1(mesh, order=order+1, dirichlet="outer", complex=True)
    Omega_rec, psi_rec = fes_rec.TnT()

    # <mu * grad(Omega), grad(psi)> = <B_A, grad(psi)>
    a_rec = BilinearForm(fes_rec)
    a_rec += mu * InnerProduct(grad(Omega_rec), grad(psi_rec)) * dx
    a_rec.Assemble()

    f_rec = LinearForm(fes_rec)
    f_rec += InnerProduct(B_A, grad(psi_rec)) * dx
    f_rec.Assemble()

    gf_rec = GridFunction(fes_rec)
    gf_rec.vec.data = a_rec.mat.Inverse(fes_rec.FreeDofs(), inverse="umfpack") * f_rec.vec

    H_rec = grad(gf_rec)

    # ============================================================
    # Compare at test point
    # ============================================================
    print("\n--- Comparison at test point ---")

    try:
        mp = mesh(test_x, test_y, test_z)
        H_A_val = H_A(mp)
        H_rec_val = H_rec(mp)

        print(f"  H_A   = ({H_A_val[0]:.4f}, {H_A_val[1]:.4f}, {H_A_val[2]:.4f})")
        print(f"  H_rec = ({H_rec_val[0]:.4f}, {H_rec_val[1]:.4f}, {H_rec_val[2]:.4f})")

        # 差分
        H_diff_val = [H_A_val[i] - H_rec_val[i] for i in range(3)]
        print(f"  Diff  = ({H_diff_val[0]:.4f}, {H_diff_val[1]:.4f}, {H_diff_val[2]:.4f})")
    except Exception as e:
        print(f"  Test point evaluation failed: {e}")
        H_A_val = [0, 0, 0]
        H_rec_val = [0, 0, 0]

    # ============================================================
    # Global norms
    # ============================================================
    print("\n--- Global norms ---")

    # ||H_A||
    H_A_norm = sqrt(Integrate(InnerProduct(H_A, Conj(H_A)), mesh).real)
    print(f"  ||H_A|| = {H_A_norm:.6e}")

    # ||H_rec||
    H_rec_norm = sqrt(Integrate(InnerProduct(H_rec, Conj(H_rec)), mesh).real)
    print(f"  ||H_rec|| = {H_rec_norm:.6e}")

    # ||H_A - H_rec||
    H_diff = H_A - H_rec
    diff_norm = sqrt(Integrate(InnerProduct(H_diff, Conj(H_diff)), mesh).real)
    print(f"  ||H_A - H_rec|| = {diff_norm:.6e}")
    print(f"  Relative diff = {diff_norm/H_A_norm*100:.2f}%")

    # Energy norm (mu-weighted)
    H_diff_mu = sqrt(Integrate(mu * InnerProduct(H_diff, Conj(H_diff)), mesh).real)
    print(f"  ||H_A - H_rec||_mu = {H_diff_mu:.6e}")

    results.append({
        'maxh': maxh,
        'ne': mesh.ne,
        'ndof_A': fes_A.ndof,
        'H_A_norm': H_A_norm,
        'H_rec_norm': H_rec_norm,
        'diff_norm': diff_norm,
        'rel_diff': diff_norm/H_A_norm*100,
        'diff_mu': H_diff_mu,
    })

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\n{'maxh[mm]':>10s} {'Elements':>10s} {'DOFs':>10s} {'||H_A-H_rec||':>14s} {'Rel.Diff %':>12s}")
print("-" * 60)

for r in results:
    print(f"{r['maxh']*1000:10.0f} {r['ne']:10d} {r['ndof_A']:10d} {r['diff_norm']:14.6e} {r['rel_diff']:12.2f}")

print("""
Note:
- H_rec is the H1-projected (curl-free) field from B_A
- In eddy current problems, curl(H) = J (not zero)
- The difference ||H_A - H_rec|| shows the rotational component of H
- For equilibrated error estimation in eddy current problems,
  a different approach (T-Omega formulation) is needed.
""")
