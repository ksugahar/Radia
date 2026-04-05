# Reduced A-formulation（還元ベクトルポテンシャル法）

本ドキュメントでは、静磁場解析における還元ベクトルポテンシャル法（Reduced A-formulation）について解説する。

## 1. 概要

### 1.1 基本的な考え方

静磁場問題において、全磁場を**ソース磁場**と**還元磁場**に分割する：

$$\mathbf{A} = \mathbf{A}_s + \mathbf{A}_r$$

ここで：
- $\mathbf{A}_s$：ソースベクトルポテンシャル（コイル電流等による既知の場）
- $\mathbf{A}_r$：還元ベクトルポテンシャル（磁性体等の影響を考慮した未知の場）

対応する磁束密度の分割は：

$$\mathbf{B} = \mathbf{B}_s + \mathbf{B}_r = \nabla \times \mathbf{A}_s + \nabla \times \mathbf{A}_r$$

### 1.2 利点

1. **計算領域の限定**：$\mathbf{A}_r$ は磁性体近傍でのみ有意な値を持つため、計算領域を限定できる
2. **境界条件の簡略化**：遠方境界で $\mathbf{A}_r = 0$ とできる
3. **精度向上**：磁性体による摂動場のみを計算するため、数値誤差が減少

## 2. 3次元定式化

### 2.1 支配方程式

マクスウェル方程式から出発する：

$$\nabla \times \mathbf{H} = \mathbf{J}$$
$$\nabla \cdot \mathbf{B} = 0$$
$$\mathbf{B} = \mu \mathbf{H}$$

ベクトルポテンシャルを導入：$\mathbf{B} = \nabla \times \mathbf{A}$

クーロンゲージ $\nabla \cdot \mathbf{A} = 0$ を課すと：

$$\nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A} \right) = \mathbf{J}$$

### 2.2 還元場の方程式

$\mathbf{A} = \mathbf{A}_s + \mathbf{A}_r$ を代入：

$$\nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A}_r \right) = \mathbf{J} - \nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A}_s \right)$$

ソース場 $\mathbf{A}_s$ が空気中（$\mu = \mu_0$）で $\nabla \times (\frac{1}{\mu_0} \nabla \times \mathbf{A}_s) = \mathbf{J}$ を満たす場合：

$$\nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A}_r \right) = \nabla \times \left( \frac{1}{\mu_0} - \frac{1}{\mu} \right) \nabla \times \mathbf{A}_s$$

### 2.3 弱形式

有限要素空間（HCurl空間）でテスト関数 $\mathbf{N}$ を用いて弱形式を導出：

$$\int_\Omega \frac{1}{\mu} (\nabla \times \mathbf{A}_r) \cdot (\nabla \times \mathbf{N}) \, d\Omega = \int_\Omega \left( \frac{1}{\mu_0} - \frac{1}{\mu} \right) (\nabla \times \mathbf{A}_s) \cdot (\nabla \times \mathbf{N}) \, d\Omega + \int_{\partial\Omega} (\mathbf{N} \times \mathbf{H}_s) \cdot \mathbf{n} \, d\Gamma$$

### 2.4 境界条件

**外部境界 $\Gamma_\infty$**：
- Bn=0 条件：$\mathbf{n} \times \mathbf{A}_r = 0$（Dirichlet）
- Ht=0 条件：$\mathbf{n} \times \mathbf{H}_r = 0$（Neumann）

**境界項**：
$$\int_{\partial\Omega} (\mathbf{N} \times \mathbf{H}_s) \cdot \mathbf{n} \, d\Gamma$$

この項はソース磁場の境界での接線成分を考慮する。

## 3. 軸対称A定式化（Z-offset Kelvin変換付き）

### 3.1 座標系と変数

軸対称問題では、円筒座標系 $(r, \theta, z)$ を用い、$\theta$ 方向の対称性を仮定する：

$$\mathbf{A} = A_\theta(r, z) \, \mathbf{e}_\theta$$

**変数変換**：$u = r A_\theta$ を導入すると、磁束密度は：

$$B_r = -\frac{1}{r} \frac{\partial u}{\partial z}, \quad B_z = \frac{1}{r} \frac{\partial u}{\partial r}$$

### 3.2 支配方程式

変数 $u$ に対する支配方程式：

$$-\nabla \cdot \left( \frac{\nu}{r} \nabla u \right) = J_\theta$$

ここで $\nu = 1/\mu$ は磁気抵抗率。

### 3.3 還元形式

$u = u_s + u_r$ と分割：

$$-\nabla \cdot \left( \frac{\nu}{r} \nabla u_r \right) = J_\theta + \nabla \cdot \left( \frac{\nu}{r} \nabla u_s \right)$$

空気中で $u_s$ が $J_\theta$ を満たす場合（電流源なし）：

$$-\nabla \cdot \left( \frac{\nu}{r} \nabla u_r \right) = \nabla \cdot \left[ \left( \frac{\nu_0 - \nu}{r} \right) \nabla u_s \right]$$

### 3.4 弱形式の導出

テスト関数 $v$ を乗じて積分：

$$\int_\Omega \frac{\nu}{r} \nabla u_r \cdot \nabla v \, r \, dr \, dz = -\int_\Omega \frac{\nu_0 - \nu}{r} \nabla u_s \cdot \nabla v \, r \, dr \, dz + \text{境界項}$$

$r$ 因子を約分すると：

$$\int_\Omega \nu \nabla u_r \cdot \nabla v \, dr \, dz = -\int_\Omega (\nu_0 - \nu) \nabla u_s \cdot \nabla v \, dr \, dz + \text{境界項}$$

**重要**：$1/r$ 重みが消えるのは、軸対称問題の特徴である。

### 3.5 ソースポテンシャル $u_s$

**一様磁束密度 $B_z = B_0$ を作る場合**：

$$B_z = \frac{1}{r} \frac{\partial u_s}{\partial r} = B_0$$

積分すると：

$$u_s = \frac{B_0 r^2}{2}$$

したがって：

$$\nabla u_s = B_0 r \, \mathbf{e}_r = B_0 \cdot (r, 0)$$

**注意**：この $u_s$ は一様 $B_z$ を与えるが、一様 $H_z$ を与えるわけではない。
$H_z = \nu B_z$ なので、$\nu$ が領域によって異なる場合、$H_z$ も異なる。

### 3.6 Z-offset Kelvin変換

軸対称問題では、外部領域を z 方向にオフセットする特殊なKelvin変換を使用する：

- **内部領域**：原点中心の半円（r-z平面）
- **外部領域**：$(0, z_{\text{offset}})$ 中心の半円

**利点**：$r$ 座標が内部と外部で同じ値を保持するため、変換後も軸対称性が維持される。

**座標変換**：
```
内部: (r, z)  →  外部: (r', z') = (r, z + z_offset)
```

**Kelvinファクター**：
$$\rho' = \sqrt{r'^2 + (z' - z_{\text{offset}})^2} = \sqrt{r^2 + z^2}$$
$$f_K = \left(\frac{\rho'}{a}\right)^2$$

ここで $a$ はKelvin境界の半径。

**透磁率の変換**：
$$\nu_{\text{outer}} = \nu_0 \cdot f_K$$

### 3.7 NGSolveによる実装

```python
# パラメータ
R_sphere = 0.5      # 磁性球の半径
a = 1.0             # Kelvin境界半径
z_offset = 3.0      # Z方向オフセット
mu_r = 100          # 比透磁率
mu0 = 4*pi*1e-7
nu0 = 1/mu0
B0 = mu0 * H0       # ソース磁束密度

# 有限要素空間（H1空間、周期境界条件付き）
fes = H1(mesh, order=3, dirichlet_bbnd='GND')
fes = Periodic(fes)
u, v = fes.TnT()

# 材料係数（内部領域）
nu_inner = nu0  # 空気
nu_magnetic = nu0 / mu_r  # 磁性体

# 材料係数（外部領域：Kelvin変換）
rho_prime = sqrt(x**2 + (y - z_offset)**2)
kelvin_factor = (rho_prime / a)**2
nu_outer = nu0 * kelvin_factor

# 双線形形式
a = BilinearForm(fes)
a += nu_inner * grad(u) * grad(v) * dx('air_inner')
a += nu_magnetic * grad(u) * grad(v) * dx('magnetic')
a += nu_outer * grad(u) * grad(v) * dx('air_outer')

# ソース項
# grad(u_s) = B0 * r * e_r = B0 * (x, 0) in (r,z) plane
grad_us_inner = B0 * CoefficientFunction((x, 0))
grad_us_outer = B0 * CoefficientFunction((x, 0))  # r座標は同じ

f = LinearForm(fes)
# 内部領域：磁性体でのみ非零（nu0 - nu_magnetic ≠ 0）
f += -(nu0 - nu_magnetic) * InnerProduct(grad_us_inner, grad(v)) * dx('magnetic')
# 外部領域：Kelvin変換によりnu0 → nu_outerだが、空気中なのでnu0 - nu0 = 0
```

## 4. Kelvin変換における境界項の消滅

### 4.1 通常の有限領域における境界項

還元ポテンシャル法の弱形式では、部分積分により境界項が発生する：

$$\int_{\partial\Omega} \nu \frac{\partial u_r}{\partial n} v \, d\Gamma$$

有限領域では、この境界項を適切に処理する必要がある（通常はDirichlet条件 $u_r = 0$ で消去）。

### 4.2 Kelvin変換による境界項の自動消滅

**周期境界条件の効果**：

Kelvin変換では、内部領域の境界と外部領域の境界が周期境界条件で結合される：

```python
inner_edge.Identify(outer_edge, "kelvin", IdentificationType.PERIODIC)
```

この結果、Kelvin境界上では：
1. 内部側からの境界項と外部側からの境界項が**符号反転**して加算される
2. 周期境界条件により $u_r$ が連続であれば、境界項は自動的に相殺される

**数学的説明**：

内部領域 $\Omega_{\text{int}}$ と外部領域 $\Omega_{\text{ext}}$ の界面 $\Gamma_K$ で：

$$\int_{\Gamma_K^-} \nu_{\text{int}} \frac{\partial u_r}{\partial n^-} v \, d\Gamma + \int_{\Gamma_K^+} \nu_{\text{ext}} \frac{\partial u_r}{\partial n^+} v \, d\Gamma = 0$$

ここで $n^- = -n^+$（法線方向が反対）かつ周期境界条件により $u_r$ と $v$ が連続。

### 4.3 境界項が消える条件

1. **周期境界条件**：Kelvin境界で $u_r$ の値が一致
2. **フラックスの連続性**：$\nu \frac{\partial u_r}{\partial n}$ が連続

Kelvin変換では透磁率が変換されるため、フラックス連続性は自動的に満たされる。

### 4.4 実装上の注意

```python
# 境界項を明示的に追加する必要がない
# 周期境界条件により自動的に処理される

a = BilinearForm(fes)
a += nu_inner * grad(u) * grad(v) * dx('air_inner')
a += nu_magnetic * grad(u) * grad(v) * dx('magnetic')
a += nu_outer * grad(u) * grad(v) * dx('air_outer')
# 境界項なし！

f = LinearForm(fes)
f += -(nu0 - nu_magnetic) * InnerProduct(grad_us, grad(v)) * dx('magnetic')
# 境界項なし！
```

## 5. 周期境界条件のためのジオメトリ構築

Kelvin変換で周期境界条件を使用する場合、NetGen/OCCでのエッジ検出が重要となる。

### 5.1 問題点

円弧を直接描画すると単一のエッジになり、周期境界条件のマッチングが困難：

```python
# NG: 単一エッジになる
wp.MoveTo(0, -a).Arc(a, 180).LineTo(0, -a).Close().Face()
```

### 5.2 解決策：完全円→カッターパターン

完全な円を作成し、矩形で切断することで上下2つのエッジを生成：

```python
# OK: 2つのエッジが生成される
inner_full = wp.Circle(a).Face()
cutter = MoveTo(-a-0.1, -a-0.1).Rectangle(a+0.1, 2*a+0.2).Face()
inner_half = inner_full - cutter
```

### 5.3 エッジ検出とマッチング

頂点距離を使用してKelvin境界のエッジを検出：

```python
kelvin_edges = []
for edge in shape.edges:
    try:
        v0, v1 = edge.vertices
        d0 = sqrt(v0.p.x**2 + v0.p.y**2)
        d1 = sqrt(v1.p.x**2 + v1.p.y**2)
        is_kelvin_arc = abs(d0 - a) < 0.01 and abs(d1 - a) < 0.01
        if is_kelvin_arc and edge.center.x > 0.01:
            edge.name = "kelvin"
            kelvin_edges.append(edge)
    except:
        pass

# y座標の符号でマッチング
for int_edge in inner_edges:
    for ext_edge in outer_edges:
        if (int_edge.center.y > 0) == (ext_edge.center.y - z_offset > 0):
            int_edge.Identify(ext_edge, "kelvin", IdentificationType.PERIODIC)
```

## 6. 検証結果

軸対称磁性球問題（$\mu_r = 100$、一様磁場中）で検証：

### 6.1 解析解

**内部** ($r < R$)：
$$H_z = \frac{3}{\mu_r + 2} H_0$$

**外部** ($r > R$)：
$$H_z = H_0 + \frac{\mu_r - 1}{\mu_r + 2} R^3 H_0 \frac{2z^2 - r^2}{(r^2 + z^2)^{5/2}}$$

### 6.2 数値結果

| 定式化 | 内部RMS誤差 | 外部RMS誤差 |
|--------|-------------|-------------|
| A-formulation + Kelvin | 0.00% | 0.02% |
| H-formulation + Kelvin | 0.00% | 0.02% |

両定式化で同等の精度が得られている。

## 7. まとめ

| 項目 | 3次元 | 軸対称 |
|------|-------|--------|
| 変数 | $\mathbf{A}$ (ベクトル) | $u = r A_\theta$ (スカラー) |
| 空間 | HCurl | H1 |
| 演算子 | $\nabla \times$ | $\nabla$ |
| 重み | なし | なし（$1/r$が約分） |
| ソース | $\nabla \times \mathbf{A}_s$ | $\nabla u_s$ |
| Kelvin変換 | X-offset | Z-offset |
| 境界項 | 周期BCで消滅 | 周期BCで消滅 |

## 参考文献

- EMPY_Analysis/Static/A_ReducedA.py (GitHub: kamearia/EMPY_Analysis)
- NGSolve documentation: https://docu.ngsolve.org/
- 実装コード: A_formulation_sphere_with_Kelvin.py
