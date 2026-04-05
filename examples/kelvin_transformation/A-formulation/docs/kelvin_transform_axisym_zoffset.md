# 軸対称A-formulation用Kelvin変換（z-offset方式）

## 概要

本文書は、軸対称静磁場問題におけるKelvin変換の実装方法を説明する。特に、外部領域をz方向（WorkPlaneのy方向）にオフセットする方式について記述する。

## z-offset方式の利点

軸対称問題で外部領域を配置する際、2つの選択肢がある：

1. **x-offset（r方向オフセット）**: 積分因子`r = x`が領域間で異なり、複雑
2. **z-offset（z方向オフセット）**: 積分因子`r = x`が両領域で同一

z-offset方式では、軸対称の回転軸からの距離`r = x`が内部・外部両領域で変わらないため、積分因子`1/r`がシンプルに保たれる。

## 座標系

- **WorkPlane座標**: `(x, y)` = `(r, z)`
- **内部領域**: 原点`(0, 0)`を中心とした半円
- **外部領域**: `(0, z_offset)`を中心とした半円

Kelvin変換における球面動径距離：
- 内部: `ρ = sqrt(x² + y²)`
- 外部: `ρ' = sqrt(x² + (y - z_offset)²)`

## A-formulationの弱形式

軸対称A-formulation（`u = r·A_θ`）の弱形式：

**内部領域**:
$$\int_{\Omega_{in}} \frac{\nu}{r} \nabla u \cdot \nabla v \, dA = \int_{\Omega_{in}} f(v) \, dA$$

**外部領域（Kelvin変換後）**:

3D球面Kelvin変換では、双線形形式に`(R/r')²`係数がかかる：
$$\int_{\Omega_{out}} \frac{\nu'}{r} \nabla u \cdot \nabla v \, dA$$

ここで：
$$\nu' = \nu_0 \cdot \left(\frac{\rho'}{a}\right)^2$$

**重要**: `(a/ρ')⁴`ではなく`(ρ'/a)²`が正しい係数！

## NGSolve実装

### ジオメトリ作成

```python
from netgen.occ import *

# パラメータ
a = 1.0           # Kelvin境界半径
z_offset = 2 * a  # 外部領域のz方向オフセット

# 内部領域（フル円を作成後、左半分をカット）
wp1 = WorkPlane()
inner_full = wp1.Circle(a).Face()

# 外部領域（オフセット位置にフル円を作成）
wp2 = WorkPlane(Axes((0, z_offset, 0), n=Z, h=X))
outer_full = wp2.Circle(a).Face()

# 左半分をカット
cutter = MoveTo(-a-0.1, -a-0.1).Rectangle(a+0.1, 2*a+0.2).Face()
inner_half = inner_full - cutter
outer_half = outer_full - MoveTo(-a-0.1, z_offset-a-0.1).Rectangle(a+0.1, 2*a+0.2).Face()
```

### エッジの命名とPeriodic BC

```python
# 頂点距離を使用してKelvin境界エッジを識別
for edge in inner_air.edges:
    cx = edge.center.x
    try:
        v0, v1 = edge.vertices
        d0 = sqrt(v0.p.x**2 + v0.p.y**2)
        d1 = sqrt(v1.p.x**2 + v1.p.y**2)
        is_kelvin = abs(d0 - a) < 0.01 and abs(d1 - a) < 0.01 and cx > 0.01
    except:
        is_kelvin = False
    if cx < 0.01:
        edge.name = 'axis'
    elif is_kelvin:
        edge.name = 'kelvin_int'

# 外部領域も同様（y座標にz_offsetを考慮）
for edge in outer_half.edges:
    cx = edge.center.x
    try:
        v0, v1 = edge.vertices
        d0 = sqrt(v0.p.x**2 + (v0.p.y - z_offset)**2)
        d1 = sqrt(v1.p.x**2 + (v1.p.y - z_offset)**2)
        is_kelvin = abs(d0 - a) < 0.01 and abs(d1 - a) < 0.01 and cx > 0.01
    except:
        is_kelvin = False
    if cx < 0.01:
        edge.name = 'axis_ext'
    elif is_kelvin:
        edge.name = 'kelvin_ext'

# Glue後にy座標の符号でマッチング
shape = Glue([inner_air, sphere_half, outer_half, gnd_point])

kelvin_int_edges = [e for e in shape.edges if e.name == 'kelvin_int']
kelvin_ext_edges = [e for e in shape.edges if e.name == 'kelvin_ext']

for int_edge in kelvin_int_edges:
    int_y = int_edge.center.y
    for ext_edge in kelvin_ext_edges:
        ext_y_local = ext_edge.center.y - z_offset
        if (int_y > 0 and ext_y_local > 0) or (int_y < 0 and ext_y_local < 0):
            int_edge.Identify(ext_edge, 'periodic', IdentificationType.PERIODIC)
            break
```

### 材料特性（Kelvin係数）

```python
# 座標変数
y_local = y - z_offset
rho_prime = sqrt(x**2 + y_local**2)
rho_prime_safe = IfPos(rho_prime - 1e-10, rho_prime, 1e-10)

# 正しいKelvin係数: (ρ'/a)² for ν
kelvin_factor = (rho_prime_safe / a)**2

# 材料リスト
materials = mesh.GetMaterials()
nu_list = []
for mat in materials:
    if 'sphere' in mat:
        nu_list.append(nu0 / mu_r)
    elif 'outer' in mat.lower():
        nu_list.append(nu0 * kelvin_factor)  # 外部領域にKelvin係数
    else:
        nu_list.append(nu0)

nu_cf = CoefficientFunction(nu_list)
```

### 弱形式

```python
# FES（Periodic BC適用）
fes = H1(mesh, order=3, dirichlet='axis|axis_ext', dirichlet_bbnd='GND')
fes = Periodic(fes)

u = fes.TrialFunction()
v = fes.TestFunction()

# 軸対称積分因子
r_weight = IfPos(x - 1e-10, x, 1e-10)

# 双線形形式
a_form = BilinearForm(fes)
a_form += nu_cf / r_weight * grad(u) * grad(v) * dx
```

## 検証結果

磁性体球（R=0.5m）、解析領域（a=1.0m）、様々な透磁率での検証：

| μ_r | 解析解 (Bz) | 数値解 (Bz) | 誤差 (%) |
|-----|------------|------------|---------|
| 2 | 1.500000e-03 | 1.500004e-03 | 0.000 |
| 10 | 2.500000e-03 | 2.500021e-03 | 0.001 |
| 100 | 2.941176e-03 | 2.941210e-03 | 0.001 |
| 1000 | 2.994012e-03 | 2.994047e-03 | 0.001 |

解析解: `Bz = 3·μ_r/(μ_r+2)·B_0`

誤差がμ_rに依存せず、0.001%以下であることがKelvin変換の正確性を示す。

## 2次元と3次元の係数比較

| 次元 | 双線形形式のKelvin係数 (on μ) | A-formulation (on ν) |
|------|------------------------------|---------------------|
| 2D | 1（変化なし） | 1（変化なし） |
| 3D | (R/r')² | (r'/R)² |

## 参考ファイル

- `Sphere_A_formulation_simple.py`: ベースライン（Kelvin変換なし）
- `Sphere_A_formulation_zoffset_Kelvin.py`: z-offset Kelvin変換
- `Coil_A_formulation_spherical_Kelvin.py`: コイル問題への応用
