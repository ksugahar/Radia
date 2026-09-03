#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A-Φ法による渦電流解析（改訂版）
A-Phi Method for Eddy Current Analysis with Rotating Magnet (Revised)

A-Φ法 (ベクトルポテンシャル-スカラーポテンシャル法):
- A_r: 渦電流による反応ベクトルポテンシャル (HCurl空間、全領域、求める未知数)
- A_ext: 外部磁場のソースポテンシャル (Radiaから取得、既知)
- A_total = A_ext + A_r (全ポテンシャル)
- Φ: 電気スカラーポテンシャル (H1空間、全領域)
- B = curl(A_total) = curl(A_ext) + curl(A_r) (磁束密度)
- E = -∂A_total/∂t - grad(Φ) (電場)
- J = sigmaE = -sigma(∂A_total/∂t + grad(Φ)) (渦電流密度)

支配方程式（導体内、A_rに関する方程式）:
式1: ∇x(1/mu ∇xA_r) + sigma(∂A_r/∂t + ∇φ) = -sigma∂A_ext/∂t  (アンペール＋ファラデー)
式2: ∇·[sigma(∂A_r/∂t + ∇φ)] = -∇·[sigma∂A_ext/∂t]  (電流連続性)

弱形式（時間離散化後、Backward Euler法）:
左辺:
  (1/mu)∫curl(A_r)·curl(W) + (sigma/Δt)∫A_r·W + (sigma/Δt)∫∇φ·W  (式1)
  (sigma/Δt)∫A_r·∇ψ + (sigma/Δt)∫∇φ·∇ψ  (式2)
右辺:
  (sigma/Δt)∫A_r,old·W  (式1: A_r履歴項)
  (sigma/Δt)∫∇φ_old·W  (式1: Φ履歴項)
  -(sigma/Δt)∫(A_ext - A_ext,old)·W  (式1: 外部ポテンシャル変化)
  (sigma/Δt)∫A_r,old·∇ψ  (式2: A_r履歴項)
  (sigma/Δt)∫∇φ_old·∇ψ  (式2: Φ履歴項)
  -(sigma/Δt)∫(A_ext - A_ext,old)·∇ψ  (式2: 外部ポテンシャル変化)

重要: A_extの正しい取得方法
- rad.RadiaField(magnet, 'a') でA_extを直接取得可能
- これによりA-Phi法の正しい実装が実現される
- 従来のB_extを使う方法は次元が不一致で誤り（README.md参照）

ゲージ条件:
nograds=True により tree-cotree gauge が自動適用される
明示的な div(A) = 0 項は不要（両方を混在させると問題が発生）

T-Omega法との違い:
- T-Omega法: 電流ポテンシャルT、J = curl(T)
- A-Φ法: 磁気ポテンシャルA、J = -sigma(∂A/∂t + grad(Φ))

更新履歴:
- 2025-12-31: Radiaから'a'フィールドでA_extを取得するように修正
  - README.mdの知見を反映し、A-Phi法の正しい実装を実現
  - B_extではなくA_extを使用することで次元の一致を確保
  - 外部磁場の時間微分をGridFunctionで正しく計算（参照渡し問題の回避）
- 2025-12-26 (第2回): ユーザー指摘により符号と係数を完全修正
  - 電流連続性の符号を正に修正（部分積分なし）
  - 渦電流密度でΦ項を復活（J = -sigma(dA/dt + ∇Φ)）
  - Φ項は「渦を巻く」ために不可欠
- 2025-12-26 (第1回): 支配方程式に基づいて定式化を完全に見直し
  式1（アンペール＋ファラデー）と式2（電流連続性）を正しく実装
- 2025-12-25: Phi結合項と電流保存則方程式を追加（標準的なA-Phi法に準拠）
  参考: https://github.com/kamearia/EMPY_Analysis および学術文献

参考文献:
- Fully discrete potential-based FEM for transient eddy current problems
- A-Phi formulation with tree-cotree gauge for eddy currents
- README.md「技術的議論: A-Phi法における外部磁場の取り扱い」
"""

import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import gc
import ngsolve
from ngsolve import *
from netgen.occ import *
import psutil

# Radiaのインポート
try:
    import radia as rad
    print(f"Radia version: {rad.UtiVer()}")
except ImportError as e:
    print(f"Error importing Radia: {e}")
    sys.exit(1)

# メモリ使用量を取得
def print_memory(label=""):
    process = psutil.Process()
    mem_mb = process.memory_info().rss / 1024 / 1024
    if label:
        print(f"  [メモリ] {label}: {mem_mb:.1f} MB")
    return mem_mb

print("=" * 80)
print("A-Φ法による渦電流解析")
print("A-Phi Method Eddy Current Analysis")
print("=" * 80)

mem_start = print_memory("プログラム起動")

# ========================================================================
# パラメータ設定
# ========================================================================

print("\n[設定] ELFMAGICモデル用パラメータ...")

# 磁石パラメータ (ELFMAGICのSample3.maiより)
magnet_size = 0.001  # m - 磁石サイズ (1mm立方体)
Br = 0.2  # T - 残留磁束密度
Hc = 127323  # A/m - 保磁力
mu0 = 4 * np.pi * 1e-7  # H/m
M = Br / mu0  # A/m - 磁化

# 移動パラメータ
x_start = -0.006  # m
x_end = 0.004  # m
y_fixed = 0.002  # m
z_fixed = 0.0  # m

# 回転パラメータ (ELFMAGICモデルに合わせる)
# 【試行5】時間刻みを細かくして過渡状態の精度を改善
time_steps = 180  # 180ステップで十分（計算時間の制約のためステップ増加は禁止）
rotation_per_step = 4.0  # deg/step - 反時計回り
total_rotation = 2 * 360  # 2回転
reset_interval = 90  # ステップごとに角度リセット
total_simulation_time = 2.0  # s
dt = total_simulation_time / time_steps  # s

# 銅板パラメータ (Sample3.maiと同じ: Y方向に配置)
copper_x_min = -0.006  # m
copper_x_max = 0.006   # m
copper_y_min = -0.0005  # m
copper_y_max = 0.0     # m
copper_z_min = -0.006  # m
copper_z_max = 0.006   # m
copper_thickness = 0.0005  # m
copper_resistivity = 1.7241e-8  # Ohm*m
copper_conductivity = 1.0 / copper_resistivity  # S/m

# メッシュパラメータ
mesh_domain = 0.006  # m - [-6, 6]^3 mm^3
mesh_maxh = 0.0005  # m (0.5mm - 高精度)

# 磁場測定グリッドのパラメータ（12mm x 12mm XZ平面、Y=0）
measurement_plane_y = 0.0  # m - 測定平面のY座標
measurement_grid_range = 0.006  # m - 測定範囲 [-6, 6] mm
measurement_grid_points = 31  # 各軸のグリッド点数 (31x31 = 961点, 0.4mm間隔)

# 回転周波数の計算（渦電流計算用）
# 2回転を2秒で実行 -> 1 Hz
frequency = total_rotation / 360.0 / total_simulation_time  # Hz

print(f"\n[渦電流計算パラメータ]")
print(f"  総シミュレーション時間: {total_simulation_time} s")
print(f"  タイムステップ数: {time_steps}")
print(f"  時間刻み dt: {dt:.6f} s ({dt*1000:.3f} ms)")
print(f"  銅の導電率: {copper_conductivity:.3e} S/m")
print(f"  銅の抵抗率: {copper_resistivity:.3e} Ohm*m")
print(f"  移動範囲: X={x_start*1000:.1f}mm -> {x_end*1000:.1f}mm")
print(f"  回転: {rotation_per_step}°/step")

# 出力設定
# スクリプトと同じディレクトリに出力フォルダを作成
script_dir = Path(__file__).parent
output_dir = script_dir / "output_comparison_A_method_order2"
output_dir.mkdir(exist_ok=True)
print(f"  出力ディレクトリ: {output_dir}")
print(f"  総ステップ数: {time_steps + 1} (0～{time_steps})")

# ========================================================================
# Step 1: メッシュ生成
# ========================================================================

print(f"\n[Step 1] NGSolve 3Dメッシュの作成（銅板含む）...")

# 銅板（導体領域）- Y方向に配置
copper_plate = Box(
    (copper_x_min, copper_y_min, copper_z_min),
    (copper_x_max, copper_y_max, copper_z_max)
)
copper_plate.mat("copper")
copper_plate.faces.name = "copper_boundary"

# 空気領域（全体から銅板を除く）
air_box = Box(
    (-mesh_domain, -mesh_domain, -mesh_domain),
    (mesh_domain, mesh_domain, mesh_domain)
)
air_box.faces.name = "outer_boundary"

# 空気 = 全体 - 銅板
air = air_box - copper_plate
air.mat("air")

# Φのポイント接地用頂点（解の一意性のため）
# 空気領域の角（外部境界上）に配置
gnd_vertex = Vertex(Pnt(mesh_domain, mesh_domain, mesh_domain))
gnd_vertex.name = "GND"

# ジオメトリ統合（GND頂点を含む）
geometry = Glue([copper_plate, air, gnd_vertex])

# メッシュ生成
geo = OCCGeometry(geometry)
ngmesh = geo.GenerateMesh(maxh=mesh_maxh)
mesh = Mesh(ngmesh)

print(f"  メッシュ要素数: {mesh.ne} 要素")
print(f"  領域: 全領域 [-6.0, 6.0]^3 mm^3")
print(f"  銅板: X=[{copper_x_min*1000:.1f}, {copper_x_max*1000:.1f}]mm, "
      f"Y=[{copper_y_min*1000:.1f}, {copper_y_max*1000:.1f}]mm, "
      f"Z=[{copper_z_min*1000:.1f}, {copper_z_max*1000:.1f}]mm")
print(f"  銅板厚: {copper_thickness*1000:.2f}mm")
print(f"  導電率: {copper_conductivity:.2e} S/m")

mem_after_mesh = print_memory("メッシュ生成後")

# ========================================================================
# 位置と回転角を計算する関数
# ========================================================================
def get_magnet_position_and_rotation(step):
    """
    指定ステップでの磁石の位置と回転角を計算

    ELFMAGICのMOV1コマンドに基づく:
    - X位置: 線形補間 (-6mm -> 4mm)
    - 回転角: 90ステップごとにリセット (0°->360°を2回繰り返し)

    ステップ0-90: 0° -> 360° (4°/step)
    ステップ91-180: 0° -> 360° (4°/step)
    """
    # X位置の線形補間
    x_pos = x_start + (x_end - x_start) * step / time_steps

    # 回転角の計算（T-Omega法と同じく連続的に増加）
    rotation_angle = step * rotation_per_step

    return x_pos, rotation_angle

# ========================================================================
# Step 1.5: 磁場測定グリッドの作成（XZ平面）
# ========================================================================
print(f"\n[Step 1.5] 磁場測定グリッドの作成（XZ平面）...")

grid_range = measurement_grid_range
grid_points_per_axis = measurement_grid_points

# 1D配列の作成
x_1d = np.linspace(-grid_range, grid_range, grid_points_per_axis)
z_1d = np.linspace(-grid_range, grid_range, grid_points_per_axis)

# XZ平面 (Y=measurement_plane_y)
x_grid_xz, z_grid_xz = np.meshgrid(x_1d, z_1d)
grid_points_xz = []
for i in range(grid_points_per_axis):
    for j in range(grid_points_per_axis):
        grid_points_xz.append([x_grid_xz[i, j], measurement_plane_y, z_grid_xz[i, j]])
grid_points_xz = np.array(grid_points_xz)

grid_spacing = 2 * grid_range / (grid_points_per_axis - 1)
print(f"  XZ平面グリッド: {grid_points_per_axis}x{grid_points_per_axis} = {len(grid_points_xz)}点 (Y={measurement_plane_y})")
print(f"  範囲: [{-grid_range*1000:.1f}, {grid_range*1000:.1f}] mm")
print(f"  グリッド間隔: {grid_spacing*1000:.2f} mm")

# データ保存用のリスト
measurement_data_xz = []

# ========================================================================
# Step 2: A-Φ法の有限要素空間設定
# ========================================================================

print(f"\n[Step 2] A-Φ法の有限要素空間設定...")

# A: ベクトルポテンシャル (HCurl空間、銅板のみで定義)
# 2026-01-05修正: 定義領域を銅板のみに制限（T-Omega法と同じ方針）
# 理由: 空気中に不要な自由度を持たせない（物理的に妥当）
# 境界条件: 全境界でAxn = 0 (接線成分ゼロ)
fesA = HCurl(mesh, order=2, nograds=True, definedon=mesh.Materials("copper"),
             dirichlet=".*", complex=False)

# Φ: スカラーポテンシャル (H1空間、全領域)
# ポイント接地: 解の一意性のため、1点（GND頂点）でΦ=0を固定
# 3Dメッシュでは dirichlet_bbbnd を使用（bbbnd = boundary of boundary of boundary = 頂点）
# 参照: Axisymmetric_sphere_with_Kelvin.py
fesPhi = H1(mesh, order=2, dirichlet_bbbnd="GND", complex=False)

# 複合空間
fesAPhi = fesA * fesPhi

# ベクトル場の投影用空間（2種類）
# 2026-01-05修正: curl演算用と高精度投影用を分離
fesVec_curl = HCurl(mesh, order=2, complex=False)  # curl演算用
fesVec_h1 = VectorH1(mesh, order=2, complex=False)  # 高精度投影用（T-Omega法と同じ）

print(f"  A空間（HCurl, 銅板のみ）の自由度数: {fesA.ndof}")
print(f"  Φ空間（H1, 全領域）の自由度数: {fesPhi.ndof}")
print(f"  A-Φ複合空間の自由度数: {fesAPhi.ndof}")
print(f"  ベクトル場空間（curl用）の自由度数: {fesVec_curl.ndof}")
print(f"  ベクトル場空間（高精度投影用）の自由度数: {fesVec_h1.ndof}")

# 試行関数と検査関数
(A, phi), (W, psi) = fesAPhi.TnT()

# 前ステップの解を保持
gfAPhi_old = GridFunction(fesAPhi)
gfAPhi_old.vec[:] = 0
gfA_old, gfPhi_old = gfAPhi_old.components

# 材料係数
sigma_dict = {"copper": copper_conductivity, "air": 0}
mu_dict = {"copper": mu0, "air": mu0}

Sigma = CF([sigma_dict.get(mat, 0) for mat in mesh.GetMaterials()])
Mu = CF([mu_dict.get(mat, mu0) for mat in mesh.GetMaterials()])

print(f"  導電率 (銅板): {copper_conductivity:.2e} S/m")
print(f"  透磁率 (真空): {mu0:.3e} H/m")
print(f"  時間刻み: {dt} s")

# ========================================================================
# 行列の事前準備（試行関数と検査関数の定義）
# ========================================================================
print("\n[Step 2.5] 行列の準備（毎ステップ再構築方式）...")

# 試行関数と検査関数
(A, phi) = fesAPhi.TrialFunction()
(W, psi) = fesAPhi.TestFunction()

print("  注意: 行列は時間ループ内で毎ステップ再構築されます")
print("  （回転・並進磁石による非一様磁場に対応）")

# 前ステップの外部A_ext保持用
Av_old = CF((0, 0, 0))  # 初期値: ゼロベクトル

# ========================================================================
# Step 3: 時間ステップループ
# ========================================================================

print(f"\n[Step 3] 時間ステップループ開始（A-Φ法による渦電流計算）...")
print(f"  方式: RadiaField + A-Φ法")
print("-" * 80)

# CSV保存用
eddy_current_data = []

for step in range(time_steps + 1):

    print(f"\n[Time Step {step}/{time_steps}]", end=" ")

    # ------------------------------------------------------------------
    # 3.1: 磁石位置・回転の計算
    # ------------------------------------------------------------------
    x_pos, rotation_angle = get_magnet_position_and_rotation(step)
    print(f"X = {x_pos*1000:.4f} mm, 回転角 = {rotation_angle:.2f} deg")

    mem_step_start = print_memory("ステップ開始")

    # ------------------------------------------------------------------
    # 3.2: Radiaで磁石を作成（T-Omega法と同じパターン）
    # ------------------------------------------------------------------

    # Delete all previous Radia objects (T-Omega法のline 420と同じ)
    rad.UtiDelAll()

    # 回転角度（ラジアン）
    theta = np.deg2rad(rotation_angle)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    # 磁石を原点で作成（初期磁化方向: Y方向）
    # ObjHexahedronを使用（ベクトルポテンシャルA対応）
    half_size = magnet_size / 2
    vertices_origin = [
        [-half_size, -half_size, -half_size],  # vertex 1
        [ half_size, -half_size, -half_size],  # vertex 2
        [ half_size,  half_size, -half_size],  # vertex 3
        [-half_size,  half_size, -half_size],  # vertex 4
        [-half_size, -half_size,  half_size],  # vertex 5
        [ half_size, -half_size,  half_size],  # vertex 6
        [ half_size,  half_size,  half_size],  # vertex 7
        [-half_size,  half_size,  half_size],  # vertex 8
    ]
    magnet = rad.ObjHexahedron(vertices_origin, [0, M, 0])

    # RadiaFieldの座標変換パラメータ（T-Omega法のline 449-452と同じ）
    origin = [x_pos, y_fixed, z_fixed]  # 並進
    u_axis = [cos_theta, sin_theta, 0]   # 回転後のX軸
    v_axis = [-sin_theta, cos_theta, 0]  # 回転後のY軸
    w_axis = [0, 0, 1]                    # 回転後のZ軸（変化なし）

    # ------------------------------------------------------------------
    # 3.3: RadiaFieldで外部ベクトルポテンシャルA_extを取得
    # ------------------------------------------------------------------

    # T-Omega法のline 455-459と同じパターンでRadiaField作成
    A_cf = rad.RadiaField(magnet, 'a',
                                    origin=origin,
                                    u_axis=u_axis,
                                    v_axis=v_axis,
                                    w_axis=w_axis)

    # GridFunctionに投影（curl演算用）
    gfA_ext = GridFunction(fesVec_curl)
    gfA_ext.Set(A_cf)

    # CoefficientFunctionも保存（T-Omega法のline 471と同じパターン）
    # 重要: A-Phi法で直接使用するため、CoefficientFunctionを保持
    A_ext_cf = A_cf

    # ★修正: 外部磁場BもRadiaFieldから直接取得（T-Omega法と同様）
    # curl(A_ext)ではなく、Radiaから直接Bを取得することで数値微分誤差を回避
    B_ext_cf = rad.RadiaField(magnet, 'b',
                                        origin=origin,
                                        u_axis=u_axis,
                                        v_axis=v_axis,
                                        w_axis=w_axis)
    gfB_ext = GridFunction(fesVec_curl)
    gfB_ext.Set(B_ext_cf)

    # デバッグ: A_extの値とdA/dtを確認、さらにcurl(A_ext) = B_extを検証
    if step <= 3:
        try:
            # 磁石に近い点でテスト（磁石中心: x_pos, y_fixed=-0.0001m, z_fixed=0.006m）
            test_point_world = [x_pos, y_fixed + 0.002, z_fixed - 0.002]  # 磁石の少し上方 (m)
            mip_center = mesh(test_point_world[0], test_point_world[1], test_point_world[2])

            # rad.RadiaField経由で取得（変換を適用した座標系）
            A_ext_val = gfA_ext(mip_center)
            print(f"  [NGSolve] A_ext: ({A_ext_val[0]:.6e}, {A_ext_val[1]:.6e}, {A_ext_val[2]:.6e})")

            # Maxwell関係式 curl(A_ext) = B_ext/mu_0 の検証
            MU_0 = 4 * np.pi * 1e-7  # H/m
            B_from_curl_A = curl(gfA_ext)(mip_center)

            # B_directも同じ変換を使用
            B_cf = rad.RadiaField(magnet, 'b',
                                            origin=origin,
                                            u_axis=u_axis,
                                            v_axis=v_axis,
                                            w_axis=w_axis)
            B_direct = B_cf(mip_center)
            print(f"  [検証] curl(A_ext) (NGSolve): ({B_from_curl_A[0]:.6e}, {B_from_curl_A[1]:.6e}, {B_from_curl_A[2]:.6e})")
            print(f"  [検証] B_direct (rad.RadiaField): ({B_direct[0]:.6e}, {B_direct[1]:.6e}, {B_direct[2]:.6e})")

            # curl(A)/B比は 1/mu_0 になるはず
            B_ratio = [B_from_curl_A[i] / B_direct[i] if abs(B_direct[i]) > 1e-12 else 0
                       for i in range(3)]
            expected_ratio = 1.0 / MU_0
            print(f"  [検証] curl(A)/B比: ({B_ratio[0]:.6e}, {B_ratio[1]:.6e}, {B_ratio[2]:.6e})")
            print(f"  [検証] 期待値 1/MU_0: {expected_ratio:.6e}")
        except Exception as e:
            print(f"  [デバッグ] 評価失敗: {e}")
            import traceback
            traceback.print_exc()

    # ------------------------------------------------------------------
    # 3.2.5: curl(A_ext) と B_ext の比較（Step 0 のみ）
    # ------------------------------------------------------------------
    # Maxwell関係式: B = curl(A) の検証
    # curl(A_ext)とRadiaから直接取得したB_extを比較
    if step == 0:
        print(f"\n  [検証] curl(A_ext) と B_ext の比較...")
        curl_B_comparison = []

        # 銅板内のグリッド点で比較
        # 銅板範囲: x: -5mm to 5mm, y: -2.5mm to 2.5mm, z: 0 to 2mm
        test_points_x = np.linspace(-0.004, 0.004, 9)  # 9点
        test_points_y = np.linspace(-0.002, 0.002, 5)  # 5点
        test_points_z = np.linspace(0.0005, 0.0015, 3)  # 3点

        total_points = 0
        valid_points = 0
        rel_errors = []

        B_cf = rad.RadiaField(magnet, 'b',
                                        origin=origin,
                                        u_axis=u_axis,
                                        v_axis=v_axis,
                                        w_axis=w_axis)

        for px in test_points_x:
            for py in test_points_y:
                for pz in test_points_z:
                    try:
                        mip = mesh(px, py, pz)

                        # curl(A_ext)を計算
                        curl_A = curl(gfA_ext)(mip)

                        # B_extをRadiaから直接取得
                        B_direct = B_cf(mip)

                        # 各成分の相対誤差を計算
                        for i in range(3):
                            if abs(B_direct[i]) > 1e-10:
                                rel_err = abs(curl_A[i] - B_direct[i]) / abs(B_direct[i]) * 100
                                rel_errors.append(rel_err)

                        # 大きさの比較
                        curl_A_mag = np.sqrt(curl_A[0]**2 + curl_A[1]**2 + curl_A[2]**2)
                        B_direct_mag = np.sqrt(B_direct[0]**2 + B_direct[1]**2 + B_direct[2]**2)

                        if B_direct_mag > 1e-10:
                            mag_rel_err = abs(curl_A_mag - B_direct_mag) / B_direct_mag * 100
                            curl_B_comparison.append({
                                'x': px, 'y': py, 'z': pz,
                                'curl_Ax': curl_A[0], 'curl_Ay': curl_A[1], 'curl_Az': curl_A[2],
                                'Bx': B_direct[0], 'By': B_direct[1], 'Bz': B_direct[2],
                                'curl_A_mag': curl_A_mag, 'B_mag': B_direct_mag,
                                'rel_error_%': mag_rel_err
                            })
                            valid_points += 1

                        total_points += 1
                    except Exception as e:
                        total_points += 1
                        continue

        if len(rel_errors) > 0:
            mean_rel_error = np.mean(rel_errors)
            max_rel_error = np.max(rel_errors)
            print(f"  [検証結果] curl(A_ext) vs B_ext:")
            print(f"    有効比較点: {valid_points}/{total_points}")
            print(f"    成分別平均相対誤差: {mean_rel_error:.2f}%")
            print(f"    成分別最大相対誤差: {max_rel_error:.2f}%")

            if len(curl_B_comparison) > 0:
                mag_errors = [d['rel_error_%'] for d in curl_B_comparison]
                print(f"    大きさ平均相対誤差: {np.mean(mag_errors):.2f}%")
                print(f"    大きさ最大相対誤差: {np.max(mag_errors):.2f}%")

        # CSVに出力
        if len(curl_B_comparison) > 0:
            curl_B_df = pd.DataFrame(curl_B_comparison)
            curl_B_csv = output_dir / "curl_A_vs_B_comparison.csv"
            curl_B_df.to_csv(curl_B_csv, index=False)
            print(f"  [出力] curl(A) vs B 比較CSV: {curl_B_csv}")

    mem_after_cf = print_memory("RadiaField作成後")

    # ------------------------------------------------------------------
    # 3.3: A-Φ法による渦電流計算（一様磁場コード準拠）
    # ------------------------------------------------------------------

    print(f"  [A-Φ法] 渦電流計算開始...")

    # GridFunction for current step
    gfAPhi = GridFunction(fesAPhi)

    # ========================================================================
    # 2026-01-05修正: 行列を毎ステップ再構築
    # 理由: 回転・並進する磁石による非一様時間変化磁場に対応
    # ========================================================================

    # 外部磁場の変化: dAv = Av_new - Av_old
    Av_new = A_ext_cf
    dAv = Av_new - Av_old

    # 左辺行列: (M + dt*K) の形
    # M: sigma*(Ar + grad(phi))*(W + grad(psi))
    # K: (1/mu)*curl(Ar)*curl(W)
    a = BilinearForm(fesAPhi, check_unused=False)

    # M の項（質量行列的）
    a += Sigma * A * W * dx("copper")
    a += Sigma * grad(phi) * W * dx("copper")
    a += Sigma * A * grad(psi) * dx("copper")
    a += Sigma * grad(phi) * grad(psi) * dx("copper")

    # dt*K の項（剛性行列的）
    a += dt * (1/Mu) * curl(A) * curl(W) * dx

    # 正則化項（BDDC前処理の安定化のため）
    penalty = 1e-6
    a += penalty * dt * (1/Mu) * A * W * dx
    a += penalty * dt * (1/Mu) * grad(phi) * grad(psi) * dx

    # 右辺 LinearForm
    f = LinearForm(fesAPhi)

    # 前ステップの履歴項: M * u_old
    # M の各項を展開
    f += Sigma * gfAPhi_old.components[0] * W * dx("copper")
    f += Sigma * grad(gfAPhi_old.components[1]) * W * dx("copper")
    f += Sigma * gfAPhi_old.components[0] * grad(psi) * dx("copper")
    f += Sigma * grad(gfAPhi_old.components[1]) * grad(psi) * dx("copper")

    # 外部磁場による誘導項: -sigma*(Av^{n+1} - Av^n)
    f += -Sigma * dAv * W * dx("copper")
    f += -Sigma * dAv * grad(psi) * dx("copper")

    # Step 0では初期条件として解をゼロに設定（T-Omega法と同様）
    if step == 0:
        print(f"  [A-Φ法] Step 0: 初期条件として渦電流をゼロに設定")
        gfAPhi.vec[:] = 0.0
    else:
        # 解く（BDDC前処理付きCG法）
        try:
            from ngsolve import Preconditioner, solvers

            # BDDC前処理器を構築（Assemble前に作成）
            pre = Preconditioner(a, "bddc")

            # アセンブル（TaskManagerを使用してマルチスレッド化）
            with TaskManager():
                a.Assemble()
                f.Assemble()
                pre.Update()
                solvers.CG(sol=gfAPhi.vec, rhs=f.vec, mat=a.mat, pre=pre.mat,
                           tol=1e-8, maxsteps=500, printrates=False)
            

            if step <= 5:
                print(f"  [A-Φ法ソルバー] BDDC+CG法で求解完了")
        except Exception as e:
            # フォールバック: 直接法
            if step <= 5:
                print(f"  [A-Φ法ソルバー] CG法失敗、直接法を使用: {e}")

            # 直接法用に再アセンブル（前処理器なしの場合）
            if not hasattr(a, 'mat') or a.mat is None:
                with TaskManager():
                    a.Assemble()
                    f.Assemble()

            gfAPhi.vec.data = a.mat.Inverse(fesAPhi.FreeDofs()) * f.vec

    # GridFunctionのコンポーネント取得
    gfA, gfPhi = gfAPhi.components

    # デバッグ: 解いたA, Phiの値を確認
    if step <= 3:
        try:
            mip_center = mesh(0, -0.00025, 0)
            A_val = gfA(mip_center)
            Phi_val = gfPhi(mip_center)
            print(f"  [デバッグ] 解いたA: ({A_val[0]:.6e}, {A_val[1]:.6e}, {A_val[2]:.6e})")
            print(f"  [デバッグ] 解いたPhi: {Phi_val:.6e}")
        except Exception as e:
            print(f"  [デバッグ] 解の評価失敗: {e}")

    # ------------------------------------------------------------------
    # 3.4: 渦電流密度の計算（時間微分ベース）
    # ------------------------------------------------------------------

    # 正しい定式化: J = -sigma * d(A + grad(Phi))/dt
    #
    # 重要: 一様磁場コードは定常正弦波を前提として J = sigma * omega * |At + grad(Phi)| を使用
    #       しかし、我々のコードは過渡解析（磁石の回転・並進）なので、
    #       時間微分から直接計算する必要があります
    #
    # At_cf = Av_new + gfA (全ベクトルポテンシャル)
    # Phi = gfPhi (スカラーポテンシャル)
    #
    # d(A + grad(Phi))/dt = (A^{n+1} + grad(Phi^{n+1}) - A^n - grad(Phi^n)) / dt

    if step > 0:
        # 現在ステップの全ポテンシャル
        At_cf_new = Av_new + gfA
        At_plus_gradPhi_new = At_cf_new + grad(gfPhi)

        # 前ステップの全ポテンシャル
        At_cf_old = Av_old + gfAPhi_old.components[0]
        At_plus_gradPhi_old = At_cf_old + grad(gfAPhi_old.components[1])

        # 時間微分
        d_At_gradPhi_dt = (At_plus_gradPhi_new - At_plus_gradPhi_old) / dt

        # 渦電流: J = -sigma * d(A + grad(Phi))/dt
        J_eddy = -Sigma * d_At_gradPhi_dt
        J_sq = J_eddy * J_eddy
        J_magnitude = sqrt(J_sq)
    else:
        # ステップ0では渦電流はゼロ
        J_eddy = CF((0, 0, 0))
        J_sq = J_eddy * J_eddy
        J_magnitude = sqrt(J_sq)

    # 磁束密度の計算
    # ★重要な修正: B = curl(A) のみ！
    # Φは電気スカラーポテンシャルであり、磁場には寄与しない
    # grad(Φ)を含めていたのは根本的な誤り
    B_internal = curl(gfA)  # <- grad(gfPhi)を削除
    # ★修正: 外部磁場はRadiaFieldから直接取得（curl(A_ext)ではなく）
    # これによりT-Omega法と同じ外部磁場を使用し、数値微分誤差を回避
    B_ext = gfB_ext  # Radiaから直接取得した外部磁場
    # 全磁束密度: B_total = B_internal + B_ext
    B_total = B_internal + B_ext

    # 渦電流RMS計算（一様磁場コード line 227準拠）
    try:
        copper_volume = Integrate(1.0, mesh, definedon=mesh.Materials("copper"))
        # J_sq を積分して平均の平方根を取る
        J_rms_val = np.sqrt(Integrate(J_sq, mesh, definedon=mesh.Materials("copper")) / copper_volume)

        print(f"  [A-Φ法] 渦電流RMS: {J_rms_val:.2e} A/m^2")

        # 最初の5ステップで詳細出力
        if step <= 5:
            try:
                mip_copper = mesh(0, -0.00025, 0)
                J_at_copper = J_eddy(mip_copper)
                J_mag_at_copper = J_magnitude(mip_copper)
                print(f"  [A-Φ法] 銅板中心の渦電流: ({J_at_copper[0]:.3e}, {J_at_copper[1]:.3e}, {J_at_copper[2]:.3e}) A/m^2")
                print(f"  [A-Φ法] 銅板中心の渦電流強度: {J_mag_at_copper:.3e} A/m^2")
            except Exception as e:
                print(f"  [警告] 銅板中心での評価失敗: {e}")

    except Exception as e:
        J_rms_val = 0.0
        print(f"  [警告] 渦電流RMS計算失敗: {e}")

    # ------------------------------------------------------------------
    # 3.4.5: ジュール発熱量の計算
    # ------------------------------------------------------------------
    # P [W] = ∫_Omega (1/sigma) |J|^2 dOmega
    # ρ = 1/sigma: 抵抗率 [Omega·m]
    # 導体領域（copper）のみで積分
    try:
        # 抵抗率 ρ = 1/sigma
        rho_copper = copper_resistivity  # Omega·m (= 1.0 / copper_conductivity)

        # 発熱密度: q = ρ |J|^2 = (1/sigma) |J|^2 [W/m^3]
        # T-Omega法と同じスカラー定数を使用（CoefficientFunctionの1/Sigmaは空気領域で不安定）
        q_density = rho_copper * J_sq

        # 導体領域での積分
        P_joule = Integrate(q_density, mesh, definedon=mesh.Materials("copper"))  # [W]

        if step <= 5 or step % 10 == 0:
            print(f"  [A-Φ法] ジュール発熱量: P = {P_joule:.6e} W")
    except Exception as e:
        P_joule = 0.0
        print(f"  [警告] 発熱量計算失敗: {e}")

    # ------------------------------------------------------------------
    # 3.4.5.5: 磁束密度BのRMS値を計算（導体領域のみ）
    # ------------------------------------------------------------------
    try:
        # B_totalの各成分の二乗和を積分
        B_squared_integral = Integrate(B_total*B_total, mesh, definedon=mesh.Materials("copper"))
        B_rms_val = np.sqrt(B_squared_integral / copper_volume)

        if step <= 5 or step % 10 == 0:
            print(f"  [A-Φ法] 磁束密度RMS: {B_rms_val:.6e} T")
    except Exception as e:
        B_rms_val = 0.0
        print(f"  [警告] 磁束密度RMS計算失敗: {e}")

    # ------------------------------------------------------------------
    # 3.4.6: ローレンツ力の計算
    # ------------------------------------------------------------------
    # ローレンツ力密度: f = J x B
    # J_eddy: 渦電流密度 [A/m^2]
    # ★試行3修正: ローレンツ力計算にはB_ext（外部磁場）のみを使用
    # 理由: B_internalは渦電流自身が作る磁場であり、T-Omega法との比較で
    #       B_internal計算方法の違いによる誤差を排除するため
    try:
        # Cross積でローレンツ力密度を計算（B_extのみ使用）
        from ngsolve import Cross
        f_lorentz_density = Cross(J_eddy, gfB_ext)  # B_total -> gfB_ext に変更

        # 銅板領域での全ローレンツ力を積分（T-Omega法と同じorder=5を使用）
        integration_order = 5
        force_x = Integrate(f_lorentz_density[0], mesh, definedon=mesh.Materials("copper"), order=integration_order)
        force_y = Integrate(f_lorentz_density[1], mesh, definedon=mesh.Materials("copper"), order=integration_order)
        force_z = Integrate(f_lorentz_density[2], mesh, definedon=mesh.Materials("copper"), order=integration_order)

        # 力の大きさを計算
        force_magnitude = np.sqrt(force_x**2 + force_y**2 + force_z**2)

        if step <= 5 or step % 10 == 0:
            print(f"  [A-Φ法] ローレンツ力: Fx={force_x:.6e} N, Fy={force_y:.6e} N, Fz={force_z:.6e} N")
            print(f"  [A-Φ法] ローレンツ力の大きさ: |F|={force_magnitude:.6e} N")
    except Exception as e:
        force_x = force_y = force_z = force_magnitude = 0.0
        print(f"  [警告] ローレンツ力計算失敗: {e}")

    # CSV用データ保存
    eddy_current_data.append({
        'step': step,
        'time': step * dt,
        'x_pos': x_pos,
        'rotation_angle': rotation_angle,
        'J_rms': J_rms_val,
        'P_joule': P_joule,
        'B_rms': B_rms_val,
        'force_x': force_x,
        'force_y': force_y,
        'force_z': force_z,
        'force_magnitude': force_magnitude
    })

    # ------------------------------------------------------------------
    # 3.4.7: A-Φ法とT-Omega法の比較用データ出力（B_internal）
    # ------------------------------------------------------------------
    # A-Φ法: B_internal = curl(A_r) （渦電流による反応磁束密度）
    # T-Omega法: B_internal = mu(T + ∇Omega)
    # この比較により、両手法の内部整合性を検証
    if step == 0 or step == 30 or step == 60 or step == 90 or step == 120 or step == 150 or step == 180:
        print(f"\n  [比較用データ] B_internal (curl(A_r)) をグリッド点で出力...")
        B_internal_data = []

        # 銅板内のグリッド点
        # 銅板は Y=[-0.5, 0.0]mm にあるので、その範囲内でサンプリング
        test_x = np.linspace(-0.003, 0.003, 7)
        test_y = np.linspace(-0.0004, -0.0001, 4)  # Y=[-0.4, -0.1]mm（銅板内部）
        test_z = np.linspace(-0.003, 0.003, 7)

        for px in test_x:
            for py in test_y:
                for pz in test_z:
                    try:
                        mip = mesh(px, py, pz)
                        B_int = B_internal(mip)
                        B_int_mag = np.sqrt(B_int[0]**2 + B_int[1]**2 + B_int[2]**2)
                        B_internal_data.append({
                            'step': step,
                            'x': px, 'y': py, 'z': pz,
                            'Bx_internal': B_int[0],
                            'By_internal': B_int[1],
                            'Bz_internal': B_int[2],
                            'B_internal_mag': B_int_mag
                        })
                    except:
                        continue

        # CSVに追記
        B_internal_csv = output_dir / "B_internal_A_method.csv"
        if step == 0:
            B_internal_df = pd.DataFrame(B_internal_data)
            B_internal_df.to_csv(B_internal_csv, index=False)
        else:
            B_internal_df = pd.DataFrame(B_internal_data)
            B_internal_df.to_csv(B_internal_csv, mode='a', header=False, index=False)

        print(f"    出力点数: {len(B_internal_data)}")

    # ------------------------------------------------------------------
    # 3.6: グリッド点での磁場測定（XZ平面）
    # ------------------------------------------------------------------
    print(f"  [測定] XZ平面で磁場を測定中...")

    # B_ext (外部磁場 = curl(A_ext))を使用して測定
    step_measurements_xz = []
    for idx, point in enumerate(grid_points_xz):
        try:
            # メッシュポイントを作成してB_extで評価
            mesh_point = mesh(point[0], point[1], point[2])
            B_vec = B_ext(mesh_point)
            B_magnitude = np.sqrt(B_vec[0]**2 + B_vec[1]**2 + B_vec[2]**2)
            step_measurements_xz.append({
                'point_id': idx + 1,
                'Gx': point[0],
                'Gy': point[1],
                'Gz': point[2],
                'Bx': B_vec[0],
                'By': B_vec[1],
                'Bz': B_vec[2],
                'B': B_magnitude
            })
        except Exception as e:
            step_measurements_xz.append({
                'point_id': idx + 1,
                'Gx': point[0],
                'Gy': point[1],
                'Gz': point[2],
                'Bx': 0.0,
                'By': 0.0,
                'Bz': 0.0,
                'B': 0.0
            })

    measurement_data_xz.append({
        'step': step,
        'x_pos': x_pos,
        'rotation_angle': rotation_angle,
        'measurements': step_measurements_xz
    })

    print(f"  [測定完了] XZ: {len(step_measurements_xz)}点")

    # ------------------------------------------------------------------
    # 3.7: 前ステップの状態を更新（次ステップのための準備）
    # ------------------------------------------------------------------
    # 一様磁場コード line 203-204準拠
    gfAPhi_old.vec.data = gfAPhi.vec
    Av_old = Av_new

    # メモリクリア
    try:
        del gfAPhi, f
    except:
        pass
    gc.collect()

    mem_step_end = print_memory("ステップ終了")
    mem_delta = mem_step_end - mem_step_start
    print(f"  [メモリ増加] このステップ: {mem_delta:+.1f} MB")

print("\n" + "=" * 80)
print("A-Φ法による渦電流解析完了")
print("=" * 80)

# ========================================================================
# Step 5: CSV出力
# ========================================================================

print(f"\n[Step 5] 渦電流データをCSVに保存...")

import csv
csv_file = output_dir / "eddy_current_data.csv"
with open(csv_file, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=['step', 'time', 'x_pos', 'rotation_angle', 'J_rms', 'P_joule', 'B_rms', 'force_x', 'force_y', 'force_z', 'force_magnitude'])
    writer.writeheader()
    writer.writerows(eddy_current_data)

print(f"  CSV保存完了: {csv_file}")
print(f"  データ行数: {len(eddy_current_data)}")

# ========================================================================
# Step 5.5: 磁場測定データをCSVファイルに書き込み
# ========================================================================
print(f"\n[Step 5.5] 磁場測定データをCSVファイルに書き込み中...")

csv_measurement_file = output_dir / "magnetic_field_measurements.csv"
with open(csv_measurement_file, 'w', newline='', encoding='utf-8') as f:
    csvwriter = csv.writer(f)

    # ヘッダー行
    csvwriter.writerow(["TimeStep", "近傍点番号", "Bx(T)", "By(T)", "Bz(T)", "B(T)", "Gx(m)", "Gy(m)", "Gz(m)"])

    # データ行
    for step_data in measurement_data_xz:
        step = step_data['step']
        for meas in step_data['measurements']:
            csvwriter.writerow([
                step,
                meas['point_id'],
                f"{meas['Bx']:.12e}",
                f"{meas['By']:.12e}",
                f"{meas['Bz']:.12e}",
                f"{meas['B']:.12e}",
                f"{meas['Gx']:.6e}",
                f"{meas['Gy']:.6e}",
                f"{meas['Gz']:.6e}"
            ])

print(f"磁場測定CSVファイル出力完了: {csv_measurement_file}")
print(f"  総データ行数: {len(measurement_data_xz) * len(grid_points_xz)}")
print(f"  タイムステップ数: {len(measurement_data_xz)}")
print(f"  各ステップの測定点数: {len(grid_points_xz)}")

# 最初の10ステップを表示
print(f"\n[結果] 最初の10ステップの渦電流RMS:")
print("  Step   J_rms (A/m^2)   Angle (deg)")
print("  " + "-" * 40)
for i in range(min(10, len(eddy_current_data))):
    data = eddy_current_data[i]
    print(f"  {data['step']:4d}   {data['J_rms']:12.2f}   {data['rotation_angle']:6.1f}")

print(f"\n完了: A-Φ法による渦電流解析")
print(f"出力ディレクトリ: {output_dir.absolute()}")
