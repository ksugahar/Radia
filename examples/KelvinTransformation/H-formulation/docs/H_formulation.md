# H-formulation による静磁界解析

## 1. 基本定式化

### 1.1 還元ポテンシャル法

全磁界 $\mathbf{H}$ を源磁界 $\mathbf{H}_s$ と還元磁界 $\mathbf{H}_r$ に分解する：

$$
\mathbf{H} = \mathbf{H}_s + \mathbf{H}_r
$$

ここで：
- $\mathbf{H}_s$：電流源や永久磁石による既知の磁界（$\nabla \times \mathbf{H}_s = \mathbf{J}$）
- $\mathbf{H}_r = -\nabla u$：スカラーポテンシャル $u$ から導かれる還元磁界

### 1.2 支配方程式

磁束密度の発散がゼロという条件 $\nabla \cdot \mathbf{B} = 0$ から：

$$
\nabla \cdot (\mu \mathbf{H}) = 0
$$

$$
\nabla \cdot (\mu (\mathbf{H}_s - \nabla u)) = 0
$$

$$
\nabla \cdot (\mu \nabla u) = \nabla \cdot (\mu \mathbf{H}_s)
$$

### 1.3 弱形式

テスト関数 $v$ を用いて：

$$
\int_\Omega \mu \nabla u \cdot \nabla v \, d\Omega = \int_\Omega \mu \mathbf{H}_s \cdot \nabla v \, d\Omega - \oint_{\partial\Omega} \mu v (\mathbf{H}_s \cdot \mathbf{n}) \, dS
$$

## 2. Kelvin変換との組み合わせ

### 2.1 境界項の自動消滅

Kelvin変換を用いると、周期境界条件により**境界積分項が自動的に消滅**する。

通常の有限領域解析では、外部境界で以下の境界項が必要：

$$
\oint_{\partial\Omega} \mu v (\mathbf{H}_s \cdot \mathbf{n}) \, dS
$$

しかし、Kelvin変換では：

1. 内部領域の外側境界と外部領域の外側境界（無限遠点に対応）が周期境界条件で結合される
2. 内部側からの境界項と外部側からの境界項が符号反転して加算される
3. $u$ が連続であれば、これらは完全に相殺される

**結果として、弱形式は単純化される：**

$$
\int_\Omega \mu \nabla u \cdot \nabla v \, d\Omega = \int_\Omega \mu \mathbf{H}_s \cdot \nabla v \, d\Omega
$$

### 2.2 実装上の利点

- 境界条件の明示的な指定が不要
- 無限遠での減衰条件が自動的に満たされる
- コードが簡潔になる

```python
# Kelvin変換を用いたH-formulation（境界項なし）
a = BilinearForm(fes)
a += mu * grad(u) * grad(v) * dx

f = LinearForm(fes)
f += mu * InnerProduct(grad(v), Hs) * dx
```

## 3. 次元別の定式化

### 3.1 2次元（平面問題）

2次元では $\mathbf{H}_s = (0, H_s)$ または $(H_s, 0)$ の形式：

$$
\int_\Omega \mu \nabla u \cdot \nabla v \, dA = \int_\Omega \mu \mathbf{H}_s \cdot \nabla v \, dA
$$

Kelvin変換後の外部領域では透磁率が位置依存となる：

$$
\mu'(\rho') = \mu_0 \left(\frac{R}{\rho'}\right)^2
$$

### 3.2 軸対称問題

軸対称問題では、積分に $r$ の重みがつく：

$$
\int_\Omega \mu \nabla u \cdot \nabla v \cdot r \, dr\,dz = \int_\Omega \mu \mathbf{H}_s \cdot \nabla v \cdot r \, dr\,dz
$$

Z-offset Kelvin変換を使用する場合、$r$ の重みは変換後も適切に変換される。

### 3.3 3次元問題

3次元では $\mathbf{H}_s = (H_{sx}, H_{sy}, H_{sz})$：

$$
\int_\Omega \mu \nabla u \cdot \nabla v \, dV = \int_\Omega \mu \mathbf{H}_s \cdot \nabla v \, dV
$$

Kelvin変換後の外部領域では：

$$
\mu'(\rho') = \mu_0 \left(\frac{R}{\rho'}\right)^4
$$

## 4. 源磁界 $\mathbf{H}_s$ の変換則

### 4.1 2次元での変換

一様磁界 $\mathbf{H}_s = H_0 \hat{y}$ に対して、Kelvin変換後：

$$
\mathbf{H}'_s = -H_0 \left(\frac{\rho'}{R}\right)^2 \hat{y}'
$$

双極子場（dipole）の場合：$\phi_s = y$ に対応
四極子場（quadrupole）の場合：$\phi_s = xy$ に対応

### 4.2 軸対称での変換

Z方向一様磁界 $\mathbf{H}_s = H_0 \hat{z}$ に対して：

- 内部領域：$\mathbf{H}_s = (0, H_0)$（r-z座標系）
- 外部領域：Z-offset Kelvin変換により適切に変換

### 4.3 3次元での変換

一様磁界 $\mathbf{H}_s = H_0 \hat{z}$ に対して、球面Kelvin変換後：

$$
\mathbf{H}'_s = -H_0 \left(\frac{\rho'}{R}\right)^3 \hat{z}'
$$

双極子場（dipole）の場合：$\phi_s = z$ に対応
四極子場（quadrupole）の場合：$\phi_s = xz$ に対応

## 5. 解析例と精度

### 5.1 磁性球の磁化問題

一様磁界中の透磁率 $\mu_r$ の球に対する解析解：

**球内部の磁界：**
$$
H_{\text{in}} = H_0 \frac{3}{\mu_r + 2}
$$

**球外部の磁界（$r > a$）：**
$$
H_r = H_0 \left(1 + \frac{2(\mu_r - 1)}{(\mu_r + 2)} \frac{a^3}{r^3}\right) \cos\theta
$$

### 5.2 数値計算結果

| 問題 | Kelvin変換なし | Kelvin変換あり |
|------|---------------|---------------|
| 2D双極子 | 境界依存 | 0.01%以下 |
| 2D四極子 | 境界依存 | 0.01%以下 |
| 軸対称双極子 | 境界依存 | 0.02%以下 |
| 3D双極子 | 境界依存 | 0.6%程度 |
| 3D四極子 | 境界依存 | 2-4%程度 |

## 6. 異方性透磁率への対応

H-formulationの2次形式では、透磁率テンソルを直接扱える：

$$
\int_\Omega \boldsymbol{\mu} \nabla u \cdot \nabla v \, d\Omega
$$

これにより、円筒異方性Kelvin変換で必要となる異方性透磁率にも対応可能。

```python
# 異方性透磁率の例
mu_tensor = CoefficientFunction((
    (mu_rr, 0),
    (0, mu_zz)
), dims=(2,2))

a += InnerProduct(mu_tensor * grad(u), grad(v)) * dx
```
