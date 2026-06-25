r"""
回転磁石と銅板の渦電流シミュレーション（過渡解析版 - T-Omega法 - CF形式アプローチ）
Rotating Magnet with Copper Plate Eddy Current (Transient Analysis - T-Omega Method - CF Approach)

【修正前版】
- 定式化: 分離形式（T*W + grad(omega)*grad(psi)）
- ソルバー: BDDC+CG法（A-Phi法と同様）
- 出力フォルダ: output_magnet1_with_eddy_current_CF_approach

モデルパラメータ（private rotating-magnet reference model より）:
- 磁石サイズ: 1mm x 1mm x 1mm
- 移動範囲: X方向 -6mm -> 4mm (合計10mm)
- Y位置: 2mm (固定)
- Z位置: 0mm (固定)
- 回転: Z軸周りに反時計回り(CCW), 4°/step
- 磁化方向: Y軸方向（磁石と同期回転）
- タイムステップ数: 180ステップ
- 総回転数: 2回転 (720°, 90ステップごとに0°リセット)

銅板パラメータ（Sample3.mai より）:
- サイズ: X=[-6, 6]mm, Y=[-0.5, 0]mm, Z=[-6, 6]mm
- 厚さ: 0.5mm (THIN 2 0.0005)
- 電気抵抗率: 1.7241x10⁻⁸ Omega·m (OHM2 2 0 1.7241e-8)
- 導電率: sigma = 5.8x10⁷ S/m

実装方式:
- RadiaFieldの座標変換機能（origin, u_axis, v_axis, w_axis）で回転・並進
- T-Omega法による渦電流計算
- 過渡解析: 各タイムステップで外部磁場の変化により渦電流が誘導される
- 【修正前版】分離形式の定式化とBDDC+CGソルバー

Author: Generated with Claude Code
Date: 2025-12-10
Modified: CF approach version (before EMPY_Analysis reference modification)
"""

import sys
import os
from pathlib import Path
import traceback

# Import radia and ngsolve
import radia as rad
import ngsolve
from ngsolve import Mesh, H1, HDiv, HCurl, L2, GridFunction, Integrate
from ngsolve import SetHeapSize  # ヒープサイズ設定用
from ngsolve import grad, curl, CF, sqrt, IfPos, x, y, z, Cross
from ngsolve import BilinearForm, LinearForm, dx
from netgen.occ import Box, Glue, OCCGeometry

import numpy as np
import pandas as pd
import gc
import psutil

print("=" * 80)
print("ELFMAGICモデル 磁石1 シミュレーション（CF形式アプローチ版）")
print("Magnet 1 Simulation from ELFMAGIC Model (CF Approach)")
print("=" * 80)

# メモリ監視用の関数
def get_memory_usage():
	"""現在のプロセスのメモリ使用量を取得（MB単位）"""
	process = psutil.Process(os.getpid())
	mem_info = process.memory_info()
	return mem_info.rss / 1024 / 1024

def print_memory(label):
	"""メモリ使用量を表示"""
	mem_mb = get_memory_usage()
	print(f"  [メモリ] {label}: {mem_mb:.1f} MB")
	return mem_mb

# 初期メモリ使用量
print("\n[メモリ監視] 初期状態")
initial_memory = print_memory("プログラム起動時")

# ========================================================================
# シミュレーション設定 (ELFMAGICモデルに基づく)
# ========================================================================
print("\n[設定] ELFMAGICモデルパラメータ...")

# Set Radia to use meters (standard for Radia and NGSolve)

# 磁石のパラメータ (Sample0.maiより)
magnet_size = 0.001  # m - 立方体の一辺の長さ (1mm x 1mm x 1mm)
Br = 0.2  # T - 残留磁束密度 (ELFMAGICのHBCU + HBSCから: 2.0 kG)
Hc = 127323  # A/m - 保磁力 (HBCU + HBSCから計算)

# 磁化Mを計算 (Radiaでは磁化をA/m単位で指定)
mu0 = 4 * np.pi * 1e-7  # H/m - 真空の透磁率
M = Br / mu0  # A/m - 磁化 (約159,155 A/m)

# 移動パラメータ (MOV1コマンドより) - 単位をメートルに統一
x_start = -0.006  # m - 初期X位置（メッシュ境界から余裕を持たせる）
x_end = 0.004     # m - 最終X位置
y_fixed = 0.002   # m - Y位置(固定) - 磁石をY=2mmに配置
z_fixed = 0.0     # m - Z位置(固定)

# 回転パラメータ (MOV1コマンドより)
# 【試行5】時間刻みを細かくして過渡状態の精度を改善
time_steps = 180  # 180ステップで十分（計算時間の制約のためステップ増加は禁止）
rotation_per_step = 4.0  # deg/step - 反時計回り
total_rotation = 2 * 360  # 2回転
reset_interval = 90  # ステップごとに角度リセット（未使用）

# ELFMAGICのシミュレーション時間設定から推定
total_simulation_time = 2.0  # s - 2回転に相当する時間

# メッシュパラメータ - 単位をメートルに統一
mesh_domain = 0.006  # m - メッシュの範囲を12mmの立方体に変更 ([-6, 6]^3 mm^3)
mesh_maxh = 0.0005  # m - 空気領域のメッシュサイズ (0.5mm - 高精度)
mesh_maxh_magnet = 0.0004  # m - 磁石領域のメッシュサイズ（速度重視のため2倍: 0.2mm -> 0.4mm）

# 磁場測定グリッドのパラメータ
measurement_plane_y = 0.0  # m - 測定平面のY座標 (Y=0mm, 銅板上面, A-Phi法と一致)
measurement_grid_range = 0.006  # m - 測定範囲 [-6, 6] mm
measurement_grid_points = 31  # 各軸のグリッド点数 (31x31 = 961点, 0.4mm間隔)

# 銅板パラメータ (Sample3.maiより)
copper_x_min = -0.006  # m
copper_x_max = 0.006   # m
copper_y_min = -0.0005 # m (厚さ0.5mm)
copper_y_max = 0.0     # m
copper_z_min = -0.006  # m
copper_z_max = 0.006   # m
copper_thickness = 0.0005  # m (0.5mm)
copper_resistivity = 1.7241e-8  # Omega·m (Sample3.mai: OHM2 2 0 1.7241e-8)
copper_conductivity = 1.0 / copper_resistivity  # S/m (約 5.8e7 S/m)

# 渦電流計算を有効化
calculate_eddy_currents = True

# T-Omega法の時間積分パラメータ
dt = total_simulation_time / time_steps  # s - タイムステップ幅（約0.0111秒）
mu_copper = mu0  # H/m - 銅の透磁率（非磁性体なのでmu_0）

print(f"\n[渦電流計算パラメータ]")
print(f"  総シミュレーション時間: {total_simulation_time} s")
print(f"  タイムステップ数: {time_steps}")
print(f"  時間刻み dt: {dt:.6f} s ({dt*1000:.3f} ms)")
print(f"  銅の導電率: {copper_conductivity:.3e} S/m")
print(f"  銅の抵抗率: {copper_resistivity:.3e} Ohm*m")
print(f"  【CF形式アプローチ】BDDC+CGソルバー使用")

# 出力設定 - 別のフォルダ名に変更
output_dir = Path(__file__).parent / "output_comparison_T_Omega_method_order2"
output_dir.mkdir(exist_ok=True)

print(f"  磁石サイズ: {magnet_size*1000:.1f}x{magnet_size*1000:.1f}x{magnet_size*1000:.1f} mm (立方体)")
print(f"  残留磁束密度 Br: {Br} T")
print(f"  保磁力 Hc: {Hc} A/m")
print(f"  移動範囲: X={x_start*1000:.1f}mm -> {x_end*1000:.1f}mm (Y={y_fixed*1000:.1f}mm, Z={z_fixed*1000:.1f}mm固定)")
print(f"  回転: {rotation_per_step}°/step, 反時計回り")
print(f"  総回転: {total_rotation/360:.0f}回転")
print(f"  タイムステップ数: {time_steps}")
print(f"  出力ディレクトリ: {output_dir}")

# ========================================================================
# Step 1: NGSolveメッシュの作成（銅板を含む）
# ========================================================================
print("\n[Step 1] NGSolve 3次元メッシュの作成（銅板を含む）...")

# 3次元立方体の空気領域を作成（磁石を囲む）
domain_m = mesh_domain
magnet_half = magnet_size / 2.0

# 銅板領域を作成（Y=-0.5mm から Y=0mm）
copper_plate = Box((copper_x_min, copper_y_min, copper_z_min),
                   (copper_x_max, copper_y_max, copper_z_max)).mat("copper")
copper_plate.maxh = 0.0005  # 銅板メッシュサイズ（0.5mm）

# 全体の空気領域（銅板を含む領域全体）
air_region = Box((-domain_m, -domain_m, -domain_m),(domain_m, domain_m, domain_m)).mat("air")
air_region.maxh = mesh_maxh

# Omegaのポイント接地用頂点（解の一意性のため）
from netgen.occ import Vertex, Pnt
gnd_vertex = Vertex(Pnt(domain_m, domain_m, domain_m))
gnd_vertex.name = "GND"

# ジオメトリを結合
geo = OCCGeometry(Glue([air_region, copper_plate, gnd_vertex]))

# メッシュ生成
mesh = Mesh(geo.GenerateMesh(maxh=mesh_maxh))

print(f"  メッシュ作成完了: {mesh.ne} 要素")
print(f"  領域: 立方体 [-{mesh_domain*1000:.1f}, {mesh_domain*1000:.1f}]^3 mm^3")
print(f"  銅板: X=[{copper_x_min*1000:.1f}, {copper_x_max*1000:.1f}]mm, "
      f"Y=[{copper_y_min*1000:.1f}, {copper_y_max*1000:.1f}]mm, "
      f"Z=[{copper_z_min*1000:.1f}, {copper_z_max*1000:.1f}]mm")
print(f"  銅板厚さ: {copper_thickness*1000:.2f}mm")
print(f"  導電率: {copper_conductivity:.2e} S/m")

# T-Omega法の有限要素空間を設定
fesT = HCurl(mesh, order=2, nograds=True, definedon=mesh.Materials("copper"),
             dirichlet=".*", complex=False)

fesOmega = H1(mesh, order=2, dirichlet_bbbnd="GND", complex=False)

# 複合空間
fesTOmega = fesT * fesOmega

print(f"  T空間（HCurl, 銅板のみ）の自由度数: {fesT.ndof}")
print(f"  Omega空間（H1, 空気領域のみ）の自由度数: {fesOmega.ndof}")
print(f"  T-Omega複合空間の自由度数: {fesTOmega.ndof}")
print(f"  磁石サイズ: {magnet_size*1000:.1f}mm立方体")
print(f"  メッシュ次元: {mesh.dim}D")
print_memory("メッシュ生成後")

# ========================================================================
# 位置と回転角を計算する関数
# ========================================================================
def get_magnet_position_and_rotation(step):
	"""指定ステップでの磁石の位置と回転角を計算"""
	x_pos = x_start + (x_end - x_start) * step / time_steps
	rotation_angle = step * rotation_per_step
	return x_pos, rotation_angle

# ========================================================================
# Step 1.6: 磁場測定グリッドの作成（3つの平面：XY, XZ, YZ）
# ========================================================================
print("\n[Step 1.6] 磁場測定グリッドの作成（XY, XZ, YZ平面）...")

# グリッド範囲
grid_range = measurement_grid_range
grid_points_per_axis = measurement_grid_points

# 1D配列の作成
x_1d = np.linspace(-grid_range, grid_range, grid_points_per_axis)
y_1d = np.linspace(-grid_range, grid_range, grid_points_per_axis)
z_1d = np.linspace(-grid_range, grid_range, grid_points_per_axis)

# XZ平面 (Y=0)
x_grid_xz, z_grid_xz = np.meshgrid(x_1d, z_1d)
grid_points_xz = []
for i in range(grid_points_per_axis):
	for j in range(grid_points_per_axis):
		grid_points_xz.append([x_grid_xz[i, j], measurement_plane_y, z_grid_xz[i, j]])
grid_points_xz = np.array(grid_points_xz)

# XY平面 (Z=0)
x_grid_xy, y_grid_xy = np.meshgrid(x_1d, y_1d)
grid_points_xy = []
for i in range(grid_points_per_axis):
	for j in range(grid_points_per_axis):
		grid_points_xy.append([x_grid_xy[i, j], y_grid_xy[i, j], 0.0])
grid_points_xy = np.array(grid_points_xy)

# YZ平面 (X=0)
y_grid_yz, z_grid_yz = np.meshgrid(y_1d, z_1d)
grid_points_yz = []
for i in range(grid_points_per_axis):
	for j in range(grid_points_per_axis):
		grid_points_yz.append([0.0, y_grid_yz[i, j], z_grid_yz[i, j]])
grid_points_yz = np.array(grid_points_yz)

# XZ平面をメインのCSV出力用として保持
grid_points = grid_points_xz
n_grid_points = len(grid_points)

grid_spacing = 2 * grid_range / (grid_points_per_axis - 1)
print(f"  XZ平面グリッド: {grid_points_per_axis}x{grid_points_per_axis} = {len(grid_points_xz)}点 (Y={measurement_plane_y})")
print(f"  XY平面グリッド: {grid_points_per_axis}x{grid_points_per_axis} = {len(grid_points_xy)}点 (Z=0)")
print(f"  YZ平面グリッド: {grid_points_per_axis}x{grid_points_per_axis} = {len(grid_points_yz)}点 (X=0)")
print(f"  範囲: [{-grid_range}, {grid_range}] mm")
print(f"  グリッド間隔: {grid_spacing:.2f} mm (ELFMAGICと一致)")

# VTK出力ディレクトリの作成
vtk_xy_dir = output_dir / "vtk_xy_plane"
vtk_xz_dir = output_dir / "vtk_xz_plane"
vtk_yz_dir = output_dir / "vtk_yz_plane"
vtk_xy_dir.mkdir(exist_ok=True)
vtk_xz_dir.mkdir(exist_ok=True)
vtk_yz_dir.mkdir(exist_ok=True)

# CSV出力ファイルの準備（XZ平面のデータ）
csv_output_file = output_dir / "magnetic_field_measurements.csv"
print(f"  CSV出力ファイル（XZ平面）: {csv_output_file}")

# ========================================================================
# Step 1.7: T-Omega法の初期化
# ========================================================================
print("\n[Step 1.7] T-Omega法の初期化...")

# 前ステップの解を保持するGridFunction（複合空間）
gfTOmega_old = GridFunction(fesTOmega)
gfTOmega_old.vec[:] = 0

# 前ステップの外部磁場を保持するGridFunction
fesVec = HDiv(mesh, order=2)
gfB_ext_old = GridFunction(fesVec)
gfB_ext_old.vec[:] = 0

# 前ステップの外部磁場CoefficientFunctionを保持
B_ext_cf_old = CF((0, 0, 0))  # 初期値: ゼロベクトル

# 材料係数
sigma_dict = {
    "copper": copper_conductivity,
    "air": 0,
}
mu_dict = {
    "copper": mu0,
    "air": mu0,
}

# CoefficientFunctionとして定義
Sigma = CF([sigma_dict.get(mat, 0) for mat in mesh.GetMaterials()])
Mu = CF([mu_dict.get(mat, mu0) for mat in mesh.GetMaterials()])

sigma = Sigma
mu = Mu

# T-Omega法の試行関数と試験関数
(T, omega), (W, psi) = fesTOmega.TnT()

print(f"  導電率 (銅板): {copper_conductivity:.2e} S/m")
print(f"  透磁率 (真空): {mu0:.3e} H/m")
print(f"  時間刻み: {dt} s")

# ========================================================================
# T-Omega法の行列を事前に構築（CF形式アプローチ - 分離形式 + BDDC+CG）
# ========================================================================
print("\n[Step 1.8] T-Omega法の行列を事前に構築...")
print("  定式化: 分離形式（修正前版）")
print("  空気領域: Mu*grad(omega)*grad(psi)")
print("  導体領域: Mu*T*W + Mu*grad(omega)*grad(psi) + (1/sigma)*curl(T)*curl(W)")
print("  ソルバー: BDDC+CG法")

from ngsolve import Preconditioner, solvers, TaskManager

# 左辺行列 A: 分離形式の定式化
a_mat = BilinearForm(fesTOmega, check_unused=False)
# 空気領域（Omegaのみ）
a_mat += (mu0/dt) * grad(omega) * grad(psi) * dx("air")
# 導体領域（分離形式: T*W と grad(omega)*grad(psi) を分離）
a_mat += (mu0/dt) * T * W * dx("copper")
a_mat += (mu0/dt) * grad(omega) * grad(psi) * dx("copper")
a_mat += (1.0/sigma) * curl(T) * curl(W) * dx("copper")

# BDDC前処理器を構築（Assemble前に作成）
pre_a = Preconditioner(a_mat, "bddc")
a_mat.Assemble()
pre_a.Update()

# 質量行列 M（右辺の前ステップ項用）- 分離形式
m_mat = BilinearForm(fesTOmega, check_unused=False)
m_mat += (mu0/dt) * grad(omega) * grad(psi) * dx("air")
m_mat += (mu0/dt) * T * W * dx("copper")
m_mat += (mu0/dt) * grad(omega) * grad(psi) * dx("copper")
m_mat.Assemble()

print(f"  行列アセンブル完了、BDDC前処理器構築完了")

# ========================================================================
# Step 2: 時間ステップループ（T-Omega法による渦電流計算）
# ========================================================================
print("\n[Step 2] 時間ステップループ開始（T-Omega法による渦電流計算）...")
print(f"  評価領域: 3次元立方体空気領域 + 銅板")
print(f"  方式: RadiaField座標変換 + T-Omega法（BDDC+CG法）")
print(f"  - 毎ステップで磁石を削除・再作成して変換を適用")
print("-" * 80)

# 結果保存用
position_data = []
rotation_data = []
integral_data = []
lorentz_force_data = []  # ローレンツ力データ
eddy_current_data = []  # 渦電流データ（CSV用）
measurement_data_xz = []  # XZ平面グリッド測定データ（CSV用）
measurement_data_xy = []  # XY平面グリッド測定データ（VTK用）
measurement_data_yz = []  # YZ平面グリッド測定データ（VTK用）

for step in range(0, time_steps + 1):  # 0から180まで
	x_pos, rotation_angle = get_magnet_position_and_rotation(step)
	theta = rotation_angle * np.pi / 180.0  # radians

	print(f"\n[Time {step/10:.1f} msec] X = {x_pos*1000:.4f} mm, 回転角 = {rotation_angle:.2f} deg")
	mem_step_start = print_memory("ステップ開始")

	# ------------------------------------------------------------------
	# 2.1: Radia変換APIを使用して磁石を作成・配置
	# ------------------------------------------------------------------
	rad.UtiDelAll()

	# 磁石を原点で作成（初期磁化方向: Y方向）
	magnet = rad.ObjRecMag([0, 0, 0], [magnet_size, magnet_size, magnet_size], [0, M, 0])

	print(f"    初期磁化ベクトル: [0, {M:.0f}, 0] A/m (Br相当: [0, {Br:.6f}, 0] T)")
	print(f"    最終位置: ({x_pos*1000:.2f}, {y_fixed*1000:.2f}, {z_fixed*1000:.2f}) mm")
	print(f"    回転角: {rotation_angle:.2f}° (Z軸周り)")

	cos_theta = np.cos(theta)
	sin_theta = np.sin(theta)

	origin = [x_pos, y_fixed, z_fixed]  # 並進
	u_axis = [cos_theta, sin_theta, 0]   # 回転後のX軸
	v_axis = [-sin_theta, cos_theta, 0]  # 回転後のY軸
	w_axis = [0, 0, 1]                    # 回転後のZ軸（変化なし）

	# RadiaField CoefficientFunctionを作成（座標変換を指定）
	B_cf = rad.RadiaField(magnet, 'b',
	                                origin=origin,
	                                u_axis=u_axis,
	                                v_axis=v_axis,
	                                w_axis=w_axis)

	mem_after_cf = print_memory("RadiaField作成後")

	# ------------------------------------------------------------------
	# 2.2: 外部磁場をGridFunctionに投影（T-Omega法用）
	# ------------------------------------------------------------------
	gfB_ext = GridFunction(fesVec)
	gfB_ext.Set(B_cf)
	B_ext_cf = B_cf  # CoefficientFunctionも保存

	mem_after_set = print_memory("外部磁場GridFunction.Set()後")

	# デバッグ: 磁石中心と銅板中心での磁場確認
	try:
		mesh_center = mesh(x_pos, y_fixed, z_fixed)
		B_at_magnet = gfB_ext(mesh_center)
		print(f"  磁石中心 @ ({x_pos*1000:.2f},{y_fixed*1000:.2f},{z_fixed*1000:.2f})mm: "
		      f"Bx={B_at_magnet[0]:+.5f}, By={B_at_magnet[1]:+.5f}, Bz={B_at_magnet[2]:+.5f} T")

		copper_center_point = mesh(0, -0.00025, 0)
		B_at_copper = gfB_ext(copper_center_point)
		print(f"  銅板中心 @ (0.00,-0.25,0.00)mm: "
		      f"Bx={B_at_copper[0]:+.5f}, By={B_at_copper[1]:+.5f}, Bz={B_at_copper[2]:+.5f} T")
	except Exception as e:
		print(f"  磁場評価エラー: {e}")

	# 外部磁場の時間変化（step>0のみ）
	if step > 0:
		try:
			copper_center_point = mesh(0, -0.00025, 0)
			B_current = gfB_ext(copper_center_point)
			B_old_val = gfB_ext_old(copper_center_point)
			dB = [(B_current[i] - B_old_val[i])/dt for i in range(3)]
			dB_norm = np.sqrt(dB[0]**2 + dB[1]**2 + dB[2]**2)
			print(f"  銅板中心での dB_ext/dt: ({dB[0]:.3e}, {dB[1]:.3e}, {dB[2]:.3e}) T/s")
			print(f"  |dB_ext/dt|: {dB_norm:.3e} T/s")
		except Exception as e:
			print(f"  dB_ext/dt評価エラー: {e}")

	# ------------------------------------------------------------------
	# 2.3: T-Omega法による渦電流計算（BDDC+CG法）
	# ------------------------------------------------------------------
	print(f"  [T-Omega] 渦電流計算開始（BDDC+CG法）...")

	# GridFunction for current step
	gfTOmega = GridFunction(fesTOmega)

	# Step 0では初期条件として解をゼロに設定
	if step == 0:
		print(f"  [T-Omega] Step 0: 初期条件として渦電流をゼロに設定")
		gfTOmega.vec[:] = 0.0

		# ------------------------------------------------------------------
		# curl(A_ext) vs B_ext の検証（Radiaから両方取得して比較）
		# ------------------------------------------------------------------
		print(f"\n  [検証] curl(A_ext) vs B_ext の比較（Radia内部検証）...")

		# RadiaからA_ext（ベクトルポテンシャル）を取得
		A_cf = rad.RadiaField(magnet, 'a',
		                                origin=origin,
		                                u_axis=u_axis,
		                                v_axis=v_axis,
		                                w_axis=w_axis)

		# A_extをHCurl空間のGridFunctionに投影してcurlを計算
		fesHCurl = HCurl(mesh, order=2)
		gfA_ext = GridFunction(fesHCurl)
		gfA_ext.Set(A_cf)

		# 銅板内のグリッド点で比較
		test_points_x = np.linspace(-0.004, 0.004, 9)  # 9点
		test_points_y = np.linspace(-0.002, 0.002, 5)  # 5点
		test_points_z = np.linspace(0.0005, 0.0015, 3)  # 3点

		curl_B_comparison = []
		rel_errors = []
		valid_points = 0
		total_points = 0

		for px in test_points_x:
			for py in test_points_y:
				for pz in test_points_z:
					try:
						mip = mesh(px, py, pz)

						# curl(A_ext)を計算
						curl_A = curl(gfA_ext)(mip)

						# B_extをRadiaから直接取得（既に取得済みのB_cfを使用）
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
			print(f"  [検証結果] curl(A_ext) vs B_ext (Radia内部):")
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
			print(f"  [検証] curl(A)=B比較結果をCSVに保存: {curl_B_csv}")

		# メモリ解放
		del gfA_ext, fesHCurl

	else:
		# dB/dt項: (B_ext^n - B_ext^{n-1})/dt
		dBdt = (B_ext_cf - B_ext_cf_old) / dt

		# 右辺ベクトルを構築
		f = LinearForm(fesTOmega)
		# 導体領域での誘導項: -dB_ext/dt * W（ファラデーの法則）
		f += -dBdt * W * dx("copper")
		# 外部磁場のgrad(psi)項（静磁場ポテンシャル）
		f += -B_ext_cf * grad(psi) * dx
		f.Assemble()

		# 右辺: M * u_old + f
		rhs = m_mat.mat * gfTOmega_old.vec + f.vec

		# BDDC+CG法で解く
		try:
			with TaskManager():
				solvers.CG(sol=gfTOmega.vec, rhs=rhs, mat=a_mat.mat, pre=pre_a.mat,
				           tol=1e-8, maxsteps=500, printrates=False)

			if step <= 5:
				print(f"  [ソルバー] BDDC+CG法で求解完了")
		except Exception as e:
			if step <= 5:
				print(f"  [ソルバー] CG法失敗、直接法を使用: {e}")
			# フォールバック: 直接法
			gfTOmega.vec.data = a_mat.mat.Inverse(fesTOmega.FreeDofs()) * rhs

	gfT, gfOmega = gfTOmega.components

	# 渦電流密度: J = curl(T)
	J_eddy_cf = curl(gfT)
	J_eddy = J_eddy_cf
	J_magnitude = sqrt(J_eddy * J_eddy)

	# 内部磁束密度: B_internal = muH_r = mu(T - ∇Omega)
	# H_r = T - ∇Omega（論文の定式化に従う）
	B_internal = mu * (gfT - grad(gfOmega))

	# 全磁束密度: B_total = B_internal + B_ext
	B_total_cf = B_internal + gfB_ext
	B_total = B_total_cf

	# 渦電流RMS（体積平均）を計算
	try:
		copper_volume = Integrate(1.0, mesh, definedon=mesh.Materials("copper"))
		J_squared_integral = Integrate(J_eddy*J_eddy, mesh, definedon=mesh.Materials("copper"))
		J_rms_val = np.sqrt(J_squared_integral / copper_volume)

		print(f"  [T-Omega] 渦電流RMS: {J_rms_val:.2e} A/m^2")

		if step == 0:
			print(f"  [T-Omega] 銅板の積分体積: {copper_volume:.6e} m^3 ({copper_volume*1e9:.3f} mm^3)")

		if step <= 5:
			try:
				mip_copper = mesh(0, -0.00025, 0)
				J_at_copper = J_eddy(mip_copper)
				J_mag_at_copper = J_magnitude(mip_copper)
				print(f"  [T-Omega] 銅板中心の渦電流: ({J_at_copper[0]:.3e}, {J_at_copper[1]:.3e}, {J_at_copper[2]:.3e}) A/m^2")
				print(f"  [T-Omega] 銅板中心の渦電流強度: {J_mag_at_copper:.3e} A/m^2")
			except Exception as e:
				print(f"  [T-Omega] 渦電流評価エラー: {e}")
	except Exception as e:
		J_rms_val = 0.0
		print(f"  [T-Omega] 渦電流RMS計算エラー: {e}")
		traceback.print_exc()

	# ------------------------------------------------------------------
	# 2.3.4.5: ジュール発熱量の計算
	# ------------------------------------------------------------------
	try:
		rho_copper = copper_resistivity
		q_density = rho_copper * (J_eddy * J_eddy)
		P_joule = Integrate(q_density, mesh, definedon=mesh.Materials("copper"))

		if step <= 5 or step % 10 == 0:
			print(f"  [T-Omega] ジュール発熱量: P = {P_joule:.6e} W")
	except Exception as e:
		P_joule = 0.0
		print(f"  [警告] 発熱量計算失敗: {e}")
		traceback.print_exc()

	# ------------------------------------------------------------------
	# 2.3.4.6: 磁束密度BのRMS値を計算（導体領域のみ）
	# ------------------------------------------------------------------
	try:
		B_squared_integral = Integrate(B_total*B_total, mesh, definedon=mesh.Materials("copper"))
		B_rms_val = np.sqrt(B_squared_integral / copper_volume)

		if step <= 5:
			print(f"  [T-Omega] 磁束密度RMS: {B_rms_val:.6e} T")
	except Exception as e:
		B_rms_val = 0.0
		print(f"  [警告] 磁束密度RMS計算失敗: {e}")

	# 渦電流データを保存
	eddy_current_data.append({
		'step': step,
		'time': step * dt,
		'x_pos': x_pos,
		'rotation_angle': rotation_angle,
		'J_rms': J_rms_val,
		'P_joule': P_joule,
		'B_rms': B_rms_val
	})

	# ------------------------------------------------------------------
	# 2.3.4.7: A-Φ法との比較用データ出力（B_reaction: 渦電流による反応磁場）
	# ------------------------------------------------------------------
	# ★重要な修正（2026-01-22）:
	# T-Omega法: B_total = mu(T + ∇Omega) は「全磁束密度」
	#        B_reaction = B_total - B_ext = mu(T + ∇Omega) - B_ext が「渦電流反応磁場」
	# A-Φ法: B_reaction = curl(A_r) が「渦電流反応磁場」
	# 比較すべきは両方の「反応磁場」部分
	if step == 0 or step == 30 or step == 60 or step == 90 or step == 120 or step == 150 or step == 180:
		print(f"\n  [比較用データ] B_internal (mu(T+∇Omega)) をグリッド点で出力...")
		B_internal_data = []

		# 銅板内のグリッド点
		# 銅板は Y=[-0.5, 0.0]mm にあるので、その範囲内でサンプリング
		test_x = np.linspace(-0.003, 0.003, 7)
		test_y = np.linspace(-0.0004, -0.0001, 4)  # Y=[-0.4, -0.1]mm（銅板内部）
		test_z = np.linspace(-0.003, 0.003, 7)

		# 渦電流反応磁場 = B_total - B_ext = mu(T+∇Omega) - B_ext
		B_reaction = B_internal - gfB_ext

		for px in test_x:
			for py in test_y:
				for pz in test_z:
					try:
						mip = mesh(px, py, pz)
						# 反応磁場（渦電流による磁場のみ）を出力
						B_react = B_reaction(mip)
						B_react_mag = np.sqrt(B_react[0]**2 + B_react[1]**2 + B_react[2]**2)
						B_internal_data.append({
							'step': step,
							'x': px, 'y': py, 'z': pz,
							'Bx_internal': B_react[0],
							'By_internal': B_react[1],
							'Bz_internal': B_react[2],
							'B_internal_mag': B_react_mag
						})
					except:
						continue

		# CSVに追記
		B_internal_csv = output_dir / "B_internal_T_Omega_method.csv"
		if step == 0:
			B_internal_df = pd.DataFrame(B_internal_data)
			B_internal_df.to_csv(B_internal_csv, index=False)
		else:
			B_internal_df = pd.DataFrame(B_internal_data)
			B_internal_df.to_csv(B_internal_csv, mode='a', header=False, index=False)

		print(f"    出力点数: {len(B_internal_data)}")

	# ------------------------------------------------------------------
	# 2.3.5: ローレンツ力の計算
	# ------------------------------------------------------------------
	# ★試行3修正: ローレンツ力計算にはB_ext（外部磁場）のみを使用
	# 理由: B_internalは渦電流自身が作る磁場であり、A-Φ法との比較で
	#       B_internal計算方法の違いによる誤差を排除するため
	try:
		f_lorentz = Cross(J_eddy, gfB_ext)  # B_total -> gfB_ext に変更

		integration_order = 5
		force_x = Integrate(f_lorentz[0], mesh, definedon=mesh.Materials("copper"), order=integration_order)
		force_y = Integrate(f_lorentz[1], mesh, definedon=mesh.Materials("copper"), order=integration_order)
		force_z = Integrate(f_lorentz[2], mesh, definedon=mesh.Materials("copper"), order=integration_order)

		force_magnitude = np.sqrt(force_x**2 + force_y**2 + force_z**2)

		if step <= 5:
			print(f"  [ローレンツ力] F = ({force_x:.3e}, {force_y:.3e}, {force_z:.3e}) N")
			print(f"  [ローレンツ力] |F| = {force_magnitude:.3e} N")
			if copper_volume > 0:
				f_density_avg = force_magnitude / copper_volume
				print(f"  [ローレンツ力] 平均力密度: {f_density_avg:.3e} N/m^3")

		lorentz_force_data.append({
			'step': step,
			'time': step * dt,
			'x_pos': x_pos,
			'rotation_angle': rotation_angle,
			'force_x': force_x,
			'force_y': force_y,
			'force_z': force_z,
			'force_magnitude': force_magnitude
		})

	except Exception as e:
		print(f"  [ローレンツ力] 計算エラー: {e}")
		import traceback
		traceback.print_exc()

	# 前ステップの解を更新
	gfTOmega_old.vec.data = gfTOmega.vec
	gfB_ext_old.vec.data = gfB_ext.vec
	B_ext_cf_old = B_ext_cf

	# ------------------------------------------------------------------
	# 2.4: 磁場積分の計算
	# ------------------------------------------------------------------
	try:
		B_magnitude_squared = gfB_ext[0]**2 + gfB_ext[1]**2 + gfB_ext[2]**2
		integral_B2 = Integrate(B_magnitude_squared, mesh)

		integral_Bx = Integrate(gfB_ext[0], mesh)
		integral_By = Integrate(gfB_ext[1], mesh)
		integral_Bz = Integrate(gfB_ext[2], mesh)

		print(f"  Int|B_ext|^2 dV = {integral_B2:.6e} T^2*m^3")
		print(f"  IntBx dV = {integral_Bx:+.6e} T*m^3")
		print(f"  IntBy dV = {integral_By:+.6e} T*m^3")
		print(f"  IntBz dV = {integral_Bz:+.6e} T*m^3")
	except Exception as e:
		print(f"  積分計算でエラー: {e}")
		import traceback
		traceback.print_exc()
		integral_B2 = 0.0
		integral_Bx = 0.0
		integral_By = 0.0
		integral_Bz = 0.0

	# ------------------------------------------------------------------
	# 2.6: グリッド点での磁場測定（3つの平面）
	# ------------------------------------------------------------------
	print(f"  [測定] 3つの平面で磁場を測定中...")

	# XZ平面の測定（CSV出力用）
	step_measurements_xz = []
	for idx, point in enumerate(grid_points_xz):
		try:
			mesh_point = mesh(point[0], point[1], point[2])
			B_vec = B_ext_cf(mesh_point)
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

	# XY平面の測定（VTK出力用）
	step_measurements_xy = []
	for idx, point in enumerate(grid_points_xy):
		try:
			mesh_point = mesh(point[0], point[1], point[2])
			B_vec = B_ext_cf(mesh_point)
			B_magnitude = np.sqrt(B_vec[0]**2 + B_vec[1]**2 + B_vec[2]**2)
			step_measurements_xy.append({
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
			step_measurements_xy.append({
				'point_id': idx + 1,
				'Gx': point[0],
				'Gy': point[1],
				'Gz': point[2],
				'Bx': 0.0,
				'By': 0.0,
				'Bz': 0.0,
				'B': 0.0
			})

	# YZ平面の測定（VTK出力用）
	step_measurements_yz = []
	for idx, point in enumerate(grid_points_yz):
		try:
			mesh_point = mesh(point[0], point[1], point[2])
			B_vec = B_ext_cf(mesh_point)
			B_magnitude = np.sqrt(B_vec[0]**2 + B_vec[1]**2 + B_vec[2]**2)
			step_measurements_yz.append({
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
			step_measurements_yz.append({
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

	measurement_data_xy.append({
		'step': step,
		'x_pos': x_pos,
		'rotation_angle': rotation_angle,
		'measurements': step_measurements_xy
	})

	measurement_data_yz.append({
		'step': step,
		'x_pos': x_pos,
		'rotation_angle': rotation_angle,
		'measurements': step_measurements_yz
	})

	print(f"  [測定完了] XZ: {len(step_measurements_xz)}点, XY: {len(step_measurements_xy)}点, YZ: {len(step_measurements_yz)}点")

	# データ保存
	position_data.append(x_pos)
	rotation_data.append(rotation_angle)
	integral_data.append({
		'B2': integral_B2,
		'Bx': integral_Bx,
		'By': integral_By,
		'Bz': integral_Bz
	})

	# メモリ解放
	try:
		del B_cf, gfB_ext
		del integral_B2, integral_Bx, integral_By, integral_Bz, B_magnitude_squared
		del magnet, origin, u_axis, v_axis, w_axis, cos_theta, sin_theta
		del gfTOmega, gfT, gfOmega
		del J_eddy, J_magnitude, B_internal, B_total
	except:
		pass

	mem_after_del = print_memory("ステップ終了")

	collected = gc.collect()
	mem_after_gc = print_memory("GC実行後")

	mem_delta = mem_after_gc - mem_step_start
	print(f"  [メモリ増加] このステップ: {mem_delta:+.1f} MB, 初期値からの増加: {mem_after_gc - initial_memory:+.1f} MB")

# ========================================================================
# Step 3: 結果の可視化とサマリー
# ========================================================================
print("\n" + "=" * 80)
print("[Step 3] シミュレーション完了 - 結果サマリー")
print("=" * 80)

# 時間変化のサマリーをファイルに保存
summary_file = output_dir / "simulation_summary.txt"
with open(summary_file, 'w', encoding='utf-8') as f:
	f.write("ELFMAGICモデル 磁石1 シミュレーション結果サマリー（CF形式アプローチ版）\n")
	f.write("=" * 80 + "\n\n")
	f.write(f"磁石形状: 立方体\n")
	f.write(f"磁石サイズ: {magnet_size*1000:.1f}x{magnet_size*1000:.1f}x{magnet_size*1000:.1f} mm\n")
	f.write(f"残留磁束密度 Br: {Br} T\n")
	f.write(f"保磁力 Hc: {Hc} A/m\n")
	f.write(f"磁化方向: Y方向（磁石と同期回転）\n")
	f.write(f"移動範囲: X={x_start*1000:.1f}mm -> {x_end*1000:.1f}mm (Y={y_fixed*1000:.1f}mm, Z={z_fixed*1000:.1f}mm固定)\n")
	f.write(f"回転: Z軸周り {rotation_per_step}°/step, 反時計回り, {total_rotation/360:.0f}回転\n")
	f.write(f"シミュレーション時間: {total_simulation_time} s\n")
	f.write(f"ソルバー: BDDC+CG法（CF形式アプローチ）\n\n")
	f.write("-" * 80 + "\n")
	f.write("時間ステップごとのデータ\n")
	f.write("-" * 80 + "\n")
	f.write(f"{'Time[s]':>8} {'X_pos[mm]':>10} {'Angle[deg]':>10} {'Int|B_ext|^2dV':>15}\n")
	f.write("-" * 80 + "\n")

	for step, (xp, angle, id_val) in enumerate(zip(position_data, rotation_data, integral_data)):
		time_s = step * dt
		f.write(f"{time_s:8.4f} {xp*1000:10.4f} {angle:10.2f} {id_val['B2']:15.6e}\n")

print(f"\nサマリーファイル: {summary_file}")
print(f"出力ディレクトリ: {output_dir}")

# ========================================================================
# Step 3.5: 磁場測定データをCSVファイルに書き込み
# ========================================================================
print("\n" + "=" * 80)
print("[Step 3.5] 磁場測定データをCSVファイルに書き込み中...")
print("=" * 80)

import csv

with open(csv_output_file, 'w', newline='', encoding='shift-jis') as csvfile:
	csvwriter = csv.writer(csvfile)

	csvwriter.writerow([str(csv_output_file.absolute())])
	csvwriter.writerow(["[M3GB] 空間磁束密度 (T)"])
	csvwriter.writerow([])
	csvwriter.writerow(["Time[s]", "近傍点番号", "Bx(T)", "By(T)", "Bz(T)", "B(T)", "Gx(m)", "Gy(m)", "Gz(m)"])

	for step_data in measurement_data_xz:
		step = step_data['step']
		time_s = step * dt
		for meas in step_data['measurements']:
			csvwriter.writerow([
				f"{time_s:.4f}",
				meas['point_id'],
				f"{meas['Bx']: .12E}",
				f"{meas['By']: .12E}",
				f"{meas['Bz']: .12E}",
				f"{meas['B']: .12E}",
				f"{meas['Gx']: .12E}",
				f"{meas['Gy']: .12E}",
				f"{meas['Gz']: .12E}"
			])

print(f"CSVファイル出力完了: {csv_output_file}")
print(f"  総データ行数: {len(measurement_data_xz) * n_grid_points}")
print(f"  タイムステップ数: {len(measurement_data_xz)}")
print(f"  各ステップの測定点数: {n_grid_points}")

# ========================================================================
# Step 3.6: 3つの平面の磁場データをParaView用VTKファイルに出力
# ========================================================================
print("\n" + "=" * 80)
print("[Step 3.6] 3つの平面の磁場データをVTKファイルに出力中...")
print("=" * 80)

# XZ平面（Y=0）
print("  XZ平面を出力中...")
for step_data in measurement_data_xz:
	step = step_data['step']
	measurements = step_data['measurements']
	vtk_file = vtk_xz_dir / f"magnetic_field_xz_step_{step:04d}.vtk"

	with open(vtk_file, 'w') as f:
		f.write("# vtk DataFile Version 3.0\n")
		f.write(f"Magnetic field on XZ plane at step {step}\n")
		f.write("ASCII\n")
		f.write("DATASET STRUCTURED_GRID\n")
		f.write(f"DIMENSIONS {grid_points_per_axis} 1 {grid_points_per_axis}\n")
		f.write(f"POINTS {len(measurements)} float\n")

		for meas in measurements:
			f.write(f"{meas['Gx']:.6f} {meas['Gy']:.6f} {meas['Gz']:.6f}\n")

		f.write(f"\nPOINT_DATA {len(measurements)}\n")
		f.write("VECTORS B_field float\n")
		for meas in measurements:
			f.write(f"{meas['Bx']:.12e} {meas['By']:.12e} {meas['Bz']:.12e}\n")

		f.write("\nSCALARS B_magnitude float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['B']:.12e}\n")

		f.write("\nSCALARS Bx float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['Bx']:.12e}\n")

		f.write("\nSCALARS By float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['By']:.12e}\n")

		f.write("\nSCALARS Bz float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['Bz']:.12e}\n")

# XY平面（Z=0）
print("  XY平面を出力中...")
for step_data in measurement_data_xy:
	step = step_data['step']
	measurements = step_data['measurements']
	vtk_file = vtk_xy_dir / f"magnetic_field_xy_step_{step:04d}.vtk"

	with open(vtk_file, 'w') as f:
		f.write("# vtk DataFile Version 3.0\n")
		f.write(f"Magnetic field on XY plane at step {step}\n")
		f.write("ASCII\n")
		f.write("DATASET STRUCTURED_GRID\n")
		f.write(f"DIMENSIONS {grid_points_per_axis} {grid_points_per_axis} 1\n")
		f.write(f"POINTS {len(measurements)} float\n")

		for meas in measurements:
			f.write(f"{meas['Gx']:.6f} {meas['Gy']:.6f} {meas['Gz']:.6f}\n")

		f.write(f"\nPOINT_DATA {len(measurements)}\n")
		f.write("VECTORS B_field float\n")
		for meas in measurements:
			f.write(f"{meas['Bx']:.12e} {meas['By']:.12e} {meas['Bz']:.12e}\n")

		f.write("\nSCALARS B_magnitude float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['B']:.12e}\n")

		f.write("\nSCALARS Bx float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['Bx']:.12e}\n")

		f.write("\nSCALARS By float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['By']:.12e}\n")

		f.write("\nSCALARS Bz float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['Bz']:.12e}\n")

# YZ平面（X=0）
print("  YZ平面を出力中...")
for step_data in measurement_data_yz:
	step = step_data['step']
	measurements = step_data['measurements']
	vtk_file = vtk_yz_dir / f"magnetic_field_yz_step_{step:04d}.vtk"

	with open(vtk_file, 'w') as f:
		f.write("# vtk DataFile Version 3.0\n")
		f.write(f"Magnetic field on YZ plane at step {step}\n")
		f.write("ASCII\n")
		f.write("DATASET STRUCTURED_GRID\n")
		f.write(f"DIMENSIONS 1 {grid_points_per_axis} {grid_points_per_axis}\n")
		f.write(f"POINTS {len(measurements)} float\n")

		for meas in measurements:
			f.write(f"{meas['Gx']:.6f} {meas['Gy']:.6f} {meas['Gz']:.6f}\n")

		f.write(f"\nPOINT_DATA {len(measurements)}\n")
		f.write("VECTORS B_field float\n")
		for meas in measurements:
			f.write(f"{meas['Bx']:.12e} {meas['By']:.12e} {meas['Bz']:.12e}\n")

		f.write("\nSCALARS B_magnitude float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['B']:.12e}\n")

		f.write("\nSCALARS Bx float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['Bx']:.12e}\n")

		f.write("\nSCALARS By float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['By']:.12e}\n")

		f.write("\nSCALARS Bz float 1\n")
		f.write("LOOKUP_TABLE default\n")
		for meas in measurements:
			f.write(f"{meas['Bz']:.12e}\n")

print(f"VTK平面データ出力完了")
print(f"  XZ平面: {len(measurement_data_xz)}ファイル -> {vtk_xz_dir}")
print(f"  XY平面: {len(measurement_data_xy)}ファイル -> {vtk_xy_dir}")
print(f"  YZ平面: {len(measurement_data_yz)}ファイル -> {vtk_yz_dir}")

# ========================================================================
# Step 4: 渦電流データのCSV出力
# ========================================================================
print("\n[Step 4] 渦電流データのCSV出力...")

eddy_csv_file = output_dir / "eddy_current_data.csv"

import pandas as pd
df_eddy = pd.DataFrame(eddy_current_data)
df_eddy.to_csv(eddy_csv_file, index=False)

print(f"  渦電流CSV出力完了: {eddy_csv_file}")
print(f"  データ数: {len(eddy_current_data)}ステップ")
print(f"  列: step, time, x_pos, rotation_angle, J_rms")

if len(df_eddy) > 0:
	print(f"\n  渦電流統計:")
	print(f"    J_rms 平均: {df_eddy['J_rms'].mean():.3e} A/m^2")
	print(f"    J_rms 最大: {df_eddy['J_rms'].max():.3e} A/m^2")
	print(f"    J_rms 最小: {df_eddy['J_rms'].min():.3e} A/m^2")

# ========================================================================
# ローレンツ力データのCSV出力
# ========================================================================
print("\n[Step 4.5] ローレンツ力データのCSV出力...")

lorentz_csv_file = output_dir / "lorentz_force_data.csv"

df_lorentz = pd.DataFrame(lorentz_force_data)
df_lorentz.to_csv(lorentz_csv_file, index=False)

print(f"  ローレンツ力CSV出力完了: {lorentz_csv_file}")
print(f"  データ数: {len(lorentz_force_data)}ステップ")
print(f"  列: step, time, x_pos, rotation_angle, force_x, force_y, force_z, force_magnitude")

if len(df_lorentz) > 0:
	print(f"\n  ローレンツ力統計:")
	print(f"    |F| 平均: {df_lorentz['force_magnitude'].mean():.3e} N")
	print(f"    |F| 最大: {df_lorentz['force_magnitude'].max():.3e} N")
	print(f"    |F| 最小: {df_lorentz['force_magnitude'].min():.3e} N")
	print(f"    Fx 平均: {df_lorentz['force_x'].mean():.3e} N")
	print(f"    Fy 平均: {df_lorentz['force_y'].mean():.3e} N")
	print(f"    Fz 平均: {df_lorentz['force_z'].mean():.3e} N")

print("\n" + "=" * 80)
print("SUCCESS: シミュレーション完了!（CF形式アプローチ版）")
print("=" * 80)

# ========================================================================
# Cleanup: Radia objects
# ========================================================================
print("\n[Cleanup] Radiaオブジェクトのクリーンアップ...")
rad.UtiDelAll()
print("  完了")
print("=" * 80)
