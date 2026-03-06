"""
TEAM Problem 7: コイル電流密度の正確な定義

レーストラックコイルの電流密度をCoefficientFunctionで定義する方法

コイル形状:
    +y方向
    ^
    |   ┌─────────────┐
    |   │  ┌─────┐   │
    |   │  │     │   │  <- 電流は反時計回り
    |   │  └─────┘   │
    |   └─────────────┘
    +-----------------> +x方向

電流の流れ:
- 上側直線部: -x方向
- 下側直線部: +x方向
- 右側直線部: +y方向
- 左側直線部: -y方向
- 角の円弧部: 接線方向（反時計回り）
"""
import os
import sys

os.environ['NETGENDIR'] = r'S:\NGSolve\01_GitHub\install_ksugahar\bin'
os.environ['PATH'] = r'S:\NGSolve\01_GitHub\install_ksugahar\bin;' + os.environ.get('PATH', '')
sys.path.insert(0, r'S:\NGSolve\01_GitHub\install_ksugahar\Lib\site-packages')

import math
from numpy import pi, sqrt
import numpy as np
from ngsolve import *
from netgen.occ import *

print("=" * 70)
print("TEAM Problem 7: Coil Current Density Definition")
print("=" * 70)

# ============================================================
# Coil parameters (TEAM Problem 7)
# ============================================================
current = 2742  # [A] Total coil current
n_turns = 1     # Number of turns (simplification)

# Coil cross-section
coil_a = 12.5e-3  # radial half-width [m]
coil_b = 50e-3    # axial half-height [m]
coil_area = (2 * coil_a) * (2 * coil_b)  # Cross-section area

# Coil path dimensions
coil_rc = 25e-3 + coil_a  # Corner radius (to center of coil cross-section)
coil_w1 = 50e-3           # Half-length of straight sections

# Coil center position
coil_x0 = 294e-3 - (coil_w1 + coil_rc + coil_a)
coil_y0 = coil_w1 + coil_rc + coil_a
coil_z0 = 19e-3 + 30e-3 + coil_b  # Center height above conductor

print(f"\nCoil parameters:")
print(f"  Current: {current} A")
print(f"  Cross-section: {2*coil_a*1000:.1f} x {2*coil_b*1000:.1f} mm")
print(f"  Area: {coil_area*1e6:.2f} mm^2")
print(f"  Corner radius: {coil_rc*1000:.1f} mm")
print(f"  Straight section length: {2*coil_w1*1000:.1f} mm")
print(f"  Center: ({coil_x0*1000:.1f}, {coil_y0*1000:.1f}, {coil_z0*1000:.1f}) mm")

# Current density magnitude
J0_mag = current * n_turns / coil_area
print(f"  Current density: {J0_mag:.2e} A/m^2")

# ============================================================
# Method 1: Segment-based current density
# ============================================================
print("\n1. Method 1: Segment-based current density")

def make_coil_current_segments():
    """
    コイルを8つのセグメントに分けて電流方向を定義
    - 4つの直線部分
    - 4つの円弧部分
    """
    # 各セグメントの境界
    # 直線部分の端点（コイル中心からの相対座標）
    inner_x = coil_w1       # 内側直線部分のx範囲
    outer_x = coil_w1       # 外側も同じ
    inner_y = coil_w1
    outer_y = coil_w1

    # コイル中心からの相対座標
    dx = x - coil_x0
    dy = y - coil_y0
    dz = z - coil_z0

    # コイル断面内かどうかをチェック
    # (これは後でジオメトリのマテリアルで制御する)

    # 各セグメントを判定
    # セグメント1: 上側直線 (y > w1, |x| < w1) -> Jx = -J0
    # セグメント2: 下側直線 (y < -w1, |x| < w1) -> Jx = +J0
    # セグメント3: 右側直線 (x > w1, |y| < w1) -> Jy = +J0
    # セグメント4: 左側直線 (x < -w1, |y| < w1) -> Jy = -J0
    # セグメント5-8: 四隅の円弧 -> 接線方向

    # 各領域を判定する関数
    # IfPos(a, b, c) = b if a > 0 else c

    # 上側直線部 (dy > w1 and |dx| < w1)
    in_top = IfPos(dy - coil_w1, 1, 0) * IfPos(coil_w1 - Norm(dx), 1, 0)

    # 下側直線部 (dy < -w1 and |dx| < w1)
    in_bottom = IfPos(-coil_w1 - dy, 1, 0) * IfPos(coil_w1 - Norm(dx), 1, 0)

    # 右側直線部 (dx > w1 and |dy| < w1)
    in_right = IfPos(dx - coil_w1, 1, 0) * IfPos(coil_w1 - Norm(dy), 1, 0)

    # 左側直線部 (dx < -w1 and |dy| < w1)
    in_left = IfPos(-coil_w1 - dx, 1, 0) * IfPos(coil_w1 - Norm(dy), 1, 0)

    # 直線部分の電流
    Jx_straight = -J0_mag * in_top + J0_mag * in_bottom
    Jy_straight = J0_mag * in_right - J0_mag * in_left

    # 円弧部分（四隅）
    # 右上隅: center at (w1, w1), 角度 0-90°
    # 左上隅: center at (-w1, w1), 角度 90-180°
    # 左下隅: center at (-w1, -w1), 角度 180-270°
    # 右下隅: center at (w1, -w1), 角度 270-360°

    # 各隅からの相対位置
    corner_ru = (dx - coil_w1, dy - coil_w1)  # 右上
    corner_lu = (dx + coil_w1, dy - coil_w1)  # 左上
    corner_ld = (dx + coil_w1, dy + coil_w1)  # 左下
    corner_rd = (dx - coil_w1, dy + coil_w1)  # 右下

    # 各隅での半径と角度
    r_ru = sqrt(corner_ru[0]**2 + corner_ru[1]**2 + 1e-20)
    r_lu = sqrt(corner_lu[0]**2 + corner_lu[1]**2 + 1e-20)
    r_ld = sqrt(corner_ld[0]**2 + corner_ld[1]**2 + 1e-20)
    r_rd = sqrt(corner_rd[0]**2 + corner_rd[1]**2 + 1e-20)

    # 円弧内かどうか（半径がcoil_rc ± coil_aの範囲内）
    in_arc_ru = IfPos(dx - coil_w1, 1, 0) * IfPos(dy - coil_w1, 1, 0)
    in_arc_lu = IfPos(-coil_w1 - dx, 1, 0) * IfPos(dy - coil_w1, 1, 0)
    in_arc_ld = IfPos(-coil_w1 - dx, 1, 0) * IfPos(-coil_w1 - dy, 1, 0)
    in_arc_rd = IfPos(dx - coil_w1, 1, 0) * IfPos(-coil_w1 - dy, 1, 0)

    # 接線方向（反時計回り）: (-y/r, x/r)
    # 右上: center at (w1, w1) -> tangent = (-(dy-w1)/r, (dx-w1)/r)
    Jx_arc_ru = -J0_mag * corner_ru[1] / r_ru * in_arc_ru
    Jy_arc_ru = J0_mag * corner_ru[0] / r_ru * in_arc_ru

    Jx_arc_lu = -J0_mag * corner_lu[1] / r_lu * in_arc_lu
    Jy_arc_lu = J0_mag * corner_lu[0] / r_lu * in_arc_lu

    Jx_arc_ld = -J0_mag * corner_ld[1] / r_ld * in_arc_ld
    Jy_arc_ld = J0_mag * corner_ld[0] / r_ld * in_arc_ld

    Jx_arc_rd = -J0_mag * corner_rd[1] / r_rd * in_arc_rd
    Jy_arc_rd = J0_mag * corner_rd[0] / r_rd * in_arc_rd

    # 全体の電流密度
    Jx = Jx_straight + Jx_arc_ru + Jx_arc_lu + Jx_arc_ld + Jx_arc_rd
    Jy = Jy_straight + Jy_arc_ru + Jy_arc_lu + Jy_arc_ld + Jy_arc_rd
    Jz = CF(0)

    return CF((Jx, Jy, Jz))

J0_segments = make_coil_current_segments()
print("  Segment-based current density defined")

# ============================================================
# Method 2: Simplified tangential current (for testing)
# ============================================================
print("\n2. Method 2: Simplified tangential current")

def make_coil_current_simple():
    """
    簡略化: コイル中心からの接線方向に電流が流れると仮定
    これは長方形コイルではなく円形コイルの近似
    """
    dx = x - coil_x0
    dy = y - coil_y0

    r = sqrt(dx**2 + dy**2 + 1e-20)

    # 接線方向（反時計回り）
    Jx = -J0_mag * dy / r
    Jy = J0_mag * dx / r
    Jz = CF(0)

    return CF((Jx, Jy, Jz))

J0_simple = make_coil_current_simple()
print("  Simple tangential current defined")

# ============================================================
# Test with simple geometry
# ============================================================
print("\n3. Testing with simple geometry...")

# Create simple test geometry
coil_outer = Box(
    gp_Pnt(coil_x0 - coil_w1 - coil_rc - coil_a,
           coil_y0 - coil_w1 - coil_rc - coil_a,
           coil_z0 - coil_b),
    gp_Pnt(coil_x0 + coil_w1 + coil_rc + coil_a,
           coil_y0 + coil_w1 + coil_rc + coil_a,
           coil_z0 + coil_b)
)
coil_inner = Box(
    gp_Pnt(coil_x0 - coil_w1 - coil_rc + coil_a,
           coil_y0 - coil_w1 - coil_rc + coil_a,
           coil_z0 - coil_b - 0.001),
    gp_Pnt(coil_x0 + coil_w1 + coil_rc - coil_a,
           coil_y0 + coil_w1 + coil_rc - coil_a,
           coil_z0 + coil_b + 0.001)
)

coil_geo = coil_outer - coil_inner
coil_geo.mat("coil")

air_box = Box(
    gp_Pnt(-0.1, -0.1, -0.1),
    gp_Pnt(0.5, 0.5, 0.3)
)
air_geo = air_box - coil_geo
air_geo.mat("air")

geo = Glue([coil_geo, air_geo])

mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=0.03))
print(f"  Mesh elements: {mesh.ne}")

# Test current density
mu0 = 4 * pi * 1e-7

# Apply current only in coil region
J0_coil = CoefficientFunction([
    J0_segments,  # coil
    CF((0, 0, 0)) # air
])

# Solve magnetostatic problem
fes = HCurl(mesh, order=2, dirichlet="")
A, v = fes.TnT()

a = BilinearForm(fes)
a += (1/mu0) * InnerProduct(curl(A), curl(v)) * dx
a += 1e-10 * InnerProduct(A, v) * dx
a.Assemble()

f = LinearForm(fes)
f += InnerProduct(J0_coil, v) * dx
f.Assemble()

gfA = GridFunction(fes)
gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec

B = curl(gfA)

# Calculate magnetic energy
W = 0.5 * Integrate((1/mu0) * InnerProduct(B, B), mesh)
print(f"  Magnetic energy: {W:.6e} J")

# Check B field at coil center
try:
    B_center = B(mesh(coil_x0, coil_y0, coil_z0 - 0.05))
    print(f"  B at (center, below coil): ({B_center[0]*1000:.4f}, {B_center[1]*1000:.4f}, {B_center[2]*1000:.4f}) mT")
except:
    print("  Could not evaluate B at center")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("Summary: Coil Current Definition")
print("=" * 70)
print(f"""
コイル電流密度の定義方法:

1. セグメント法 (Method 1):
   - コイルを直線部と円弧部に分割
   - 各セグメントで電流方向を定義
   - 最も正確だが、複雑

2. 簡略法 (Method 2):
   - コイル中心からの接線方向
   - 円形コイルの近似
   - 簡単だが、不正確

実際のコイル電流:
  J = J0 * (電流方向の単位ベクトル)
  J0 = I / A = {current} / {coil_area*1e6:.2f} mm^2 = {J0_mag:.2e} A/m^2

注意:
  - コイル領域のマテリアルでJ0を適用
  - A法: curl(1/mu * curl(A)) = J0
  - div(J0) = 0 を満たす必要あり（連続の式）
""")
