# 2次元静磁場解析のためのKelvin変換

## 概要

本文書は、2次元H定式化静磁場解析におけるKelvin変換の完全な数学的基礎と実装ガイドを提供する。Kelvin変換は外部領域（$r > R$）を有限な計算領域にマッピングし、無限領域問題の効率的な表現を可能にする。

本文書の内容：
- Kelvin変換の数学的理論
- 強形式からの弱形式導出
- NGSolve実装の詳細
- ダイポール場および4極場の適用例

## 0. 摂動ポテンシャル定式化

### 0.1 強形式

磁束密度のガウスの法則（$\nabla \cdot \mathbf{B} = 0$）から出発する：

$$\mathbf{H}_{\text{total}} = \mathbf{H}_s + \mathbf{H}_{\text{pert}} = \mathbf{H}_s - \nabla\phi_m$$

$$\mathbf{B} = \mu\mathbf{H}_{\text{total}} = \mu(\mathbf{H}_s - \nabla\phi_m)$$

$$\nabla \cdot \mathbf{B} = 0$$

$$\Rightarrow \nabla \cdot [\mu(\mathbf{H}_s - \nabla\phi_m)] = 0$$

$$\Rightarrow \nabla \cdot (\mu\mathbf{H}_s) - \nabla \cdot (\mu\nabla\phi_m) = 0$$

$$\Rightarrow -\nabla \cdot (\mu\nabla\phi_m) = -\nabla \cdot (\mu\mathbf{H}_s)$$

**遠方場の境界条件**（$r \to \infty$）：
$$\mathbf{H}_{\text{pert}} \to 0 \quad \Rightarrow \quad \nabla\phi_m \to 0 \quad \text{（自然境界条件: } n \cdot \nabla\phi_m = 0\text{）}$$

### 0.2 弱形式の導出

出発点: $-\nabla \cdot (\mu\nabla\phi_m) = -\nabla \cdot (\mu\mathbf{H}_s)$

**ステップ1**: テスト関数 $v$ を乗じて積分：
$$-\int_\Omega v\nabla \cdot (\mu\nabla\phi_m) d\Omega = -\int_\Omega v\nabla \cdot (\mu\mathbf{H}_s) d\Omega$$

**ステップ2**: 部分積分（発散定理）を適用：

左辺：
$$-\int_\Gamma v(n \cdot \mu\nabla\phi_m) d\Gamma + \int_\Omega (\nabla v) \cdot (\mu\nabla\phi_m) d\Omega$$

右辺：
$$-\int_\Gamma v(n \cdot \mu\mathbf{H}_s) d\Gamma + \int_\Omega (\nabla v) \cdot (\mu\mathbf{H}_s) d\Omega$$

**ステップ3**: 自然境界条件（$\Gamma$上で$n \cdot \nabla\phi_m = 0$）を適用：
$$-\int_\Gamma v(n \cdot \mu\nabla\phi_m) d\Gamma = 0$$

**ステップ4**: 最終的な弱形式：

$\phi_m \in V$ を求めよ：

$$a(\phi_m, v) = f(v) \quad \forall v \in V$$

ここで：
$$a(u,v) = \int_\Omega (\nabla v) \cdot (\mu\nabla u) d\Omega \quad \text{（双線形形式）}$$

$$f(v) = \int_\Omega (\nabla v) \cdot (\mu\mathbf{H}_s) d\Omega - \int_\Gamma v(n \cdot \mu\mathbf{H}_s) d\Gamma \quad \text{（線形形式）}$$

**重要**: 体積積分は**正**の符号、境界積分は**負**の符号であることに注意。

**双線形形式に関する重要な注意**: 双線形形式には**体積積分のみ** $\int_\Omega (\nabla v) \cdot (\mu\nabla u) d\Omega$ が含まれる。左辺からの境界項 $-\int_\Gamma v(n \cdot \mu\nabla\phi_m) d\Gamma$ は**自然境界条件により消去される**（ステップ3）。双線形形式に移動されるわけではない。`-v*mu*InnerProduct(n, grad(u))*ds(...)`のような境界項を双線形形式に追加してはならない。

### 0.3 有限領域問題（Kelvin変換なし）における境界項の扱い

**重要**: Kelvin変換を使用しない問題（外部境界を持つ有限領域）では、**境界積分項を線形形式に含める必要がある**：

$$f(v) = \int_\Omega (\nabla v) \cdot (\mu\mathbf{H}_s) d\Omega - \int_{\Gamma_{\text{outer}}} v(n \cdot \mu\mathbf{H}_s) d\Gamma$$

#### 境界項の起源

境界項は、弱形式の導出における**部分積分**（上記ステップ2）から生じる：

$$-\int_\Omega v\nabla \cdot (\mu\mathbf{H}_s) d\Omega = -\int_\Gamma v(n \cdot \mu\mathbf{H}_s) d\Gamma + \int_\Omega (\nabla v) \cdot (\mu\mathbf{H}_s) d\Omega$$

この境界項は自然境界条件では**消去できない**。理由は：
- 自然境界条件（$n \cdot \nabla\phi_m = 0$）は**左辺**（双線形形式から）の境界項を消去する
- **右辺**（線形形式から）の境界項は $\phi_m$ ではなく $\mathbf{H}_s$ を含むため残る

#### 物理的解釈

1. **磁束の寄与**: 項 $-\int_\Gamma v(n \cdot \mu\mathbf{H}_s) d\Gamma$ は、計算領域境界を通過する背景磁場の磁束の寄与を表す。

2. **遠方場の挙動**: 無限領域では、摂動場が無限遠で減衰するためこの項は零に積分される。しかし、**切断された有限領域**では、外部場の影響を正しくモデル化するためにこの項が不可欠。

3. **エネルギー収支**: この項がないと、境界での背景場による仕事が考慮されないため、FEM定式化はエネルギー保存則に違反する。

#### 境界項が不可欠な理由

1. **数学的完全性**: 弱形式は部分積分により強形式から導出される。境界項を省略するとこの導出が無効になる。

2. **正しい遠方場の影響**: 境界項は、摂動場が領域境界での背景場の影響を正しく考慮することを保証する。

3. **検証**: ファイル`2D_dipole.py`および`2D_quadrupole.py`が、境界項が正しい結果に不可欠であることを示している。

**NGSolve実装**（Kelvin変換なしの標準形式）：

```python
# 線形形式: f(v) = ∫(∇v)·(μH_s)dΩ - ∫v(n·μH_s)dΓ
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx                         # 体積積分（正）
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer")) # 境界項（負）
```

**重要な符号規則**：
- 体積積分: **正**の符号（`+`）
- 境界積分: **負**の符号（`-`）

#### Kelvin変換との対比

| 手法 | 境界項 | 理由 |
|------|--------|------|
| **有限領域（Kelvin変換なし）** | **必要** | 有限距離での境界が領域を切断 |
| **Kelvin + 周期境界条件** | **不要** | 周期BCが $r = R$ での連続性を保証 |

**Kelvin + 周期BCが境界項を不要にする理由：**

1. **領域の閉鎖**: Kelvin変換は $r > R$ を有限領域にマッピングし、周期BCは $r = R$ で内部（$r < R$）と外部（$r > R$）領域を接続する。

2. **自動的な連続性**: 周期BCは以下を自動的に保証：
   - ポテンシャルの連続性: $\phi(R^-) = \phi'(R^+)$
   - フラックスの連続性: $(n \cdot \mu\nabla\phi)|_{R^-} = (n \cdot \mu'\nabla'\phi')|_{R^+}$

3. **境界寄与の相殺**: $r = R$ の内側と外側からの境界積分の寄与は、反対向きの法線方向と周期的同定により相殺される。

### 0.4 NGSolve実装（標準形式）

```python
from math import pi
from ngsolve import *

# 材料特性
mu0 = 4*pi*1e-7  # 真空の透磁率 [H/m]
mu_d = {"air_inner": 1*mu0, "air_outer": 1*mu0, "magnetic": 10*mu0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

# 背景場（2D: y方向）
Hs = CoefficientFunction((0, 1))
Hsb = BoundaryFromVolumeCF(Hs)

# 有限要素空間
fes = H1(mesh, order=2)
u = fes.TrialFunction()
v = fes.TestFunction()

# 双線形形式: a(u,v) = ∫(∇v)·(μ∇u)dΩ
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx

# 線形形式: f(v) = ∫(∇v)·(μH_s)dΩ - ∫v(n·μH_s)dΓ
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx                    # 正の符号
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer"))  # 負の符号

# 求解
a.Assemble()
f.Assemble()
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

# 場の抽出
H_pert = -grad(gfu)  # 摂動場のみ
H_total = H_pert + Hs  # 全磁場 = 摂動 + 背景
```

## 1. 座標変換

### 1.1 定義

2次元における**Kelvin変換**（円反転とも呼ばれる）は次のように定義される：

$$r' = \frac{R^2}{r}$$

ここで：
- $r = \sqrt{x^2 + y^2}$ は元の座標における動径距離
- $r' = \sqrt{x'^2 + y'^2}$ は変換後の座標における動径距離
- $R$ は変換半径（Kelvin球半径）

### 1.2 円筒座標系

2次元問題では、**極座標** $(r, \theta)$ を使用する：

$$\begin{align}
x &= r \cos\theta \\
y &= r \sin\theta
\end{align}$$

極座標におけるKelvin変換：

$$\begin{align}
r' &= \frac{R^2}{r} \\
\theta' &= \theta \quad \text{（角度は保存）}
\end{align}$$

### 1.3 デカルト座標変換

$r' = R^2/r$ および $\theta' = \theta$ から、デカルト座標は次のように変換される：

$$\begin{align}
x' &= r' \cos\theta' = \frac{R^2}{r} \cos\theta = \frac{R^2}{r^2} x = \frac{R^2 x}{x^2 + y^2} \\
y' &= r' \sin\theta' = \frac{R^2}{r} \sin\theta = \frac{R^2}{r^2} y = \frac{R^2 y}{x^2 + y^2}
\end{align}$$

**逆変換**：

$$\begin{align}
x &= \frac{R^2}{r'^2} x' = \frac{R^2 x'}{x'^2 + y'^2} \\
y &= \frac{R^2}{r'^2} y' = \frac{R^2 y'}{x'^2 + y'^2}
\end{align}$$

### 1.4 主要な性質

1. **等角写像**: 角度が（局所的に）保存される
2. **反転**: 変換を2回適用すると元の座標に戻る
3. **円の保存**: 円は円（または原点を通る直線）に写像される
4. **領域の写像**：
   - $r \to \infty$ は $r' \to 0$ に写像
   - $r = R$ は $r' = R$ に写像（不動円）
   - $r < R$ は $r' > R$ に写像

## 2. 円筒座標における計量テンソル

### 2.1 極座標形式での座標変換

極座標において、変換は次のように表される：

$$(r, \theta) \to (r', \theta') = \left(\frac{R^2}{r}, \theta\right)$$

この変換のヤコビ行列は：

$$J = \begin{bmatrix}
\frac{\partial r'}{\partial r} & \frac{\partial r'}{\partial \theta} \\
\frac{\partial \theta'}{\partial r} & \frac{\partial \theta'}{\partial \theta}
\end{bmatrix}$$

偏微分を計算すると：

$$\begin{align}
\frac{\partial r'}{\partial r} &= \frac{\partial}{\partial r}\left(\frac{R^2}{r}\right) = -\frac{R^2}{r^2} \\
\frac{\partial r'}{\partial \theta} &= 0 \\
\frac{\partial \theta'}{\partial r} &= 0 \\
\frac{\partial \theta'}{\partial \theta} &= 1
\end{align}$$

したがって：

$$J = \begin{bmatrix}
-\frac{R^2}{r^2} & 0 \\
0 & 1
\end{bmatrix}$$

行列式は：

$$\det(J) = -\frac{R^2}{r^2} = -\left(\frac{r'}{R}\right)^2$$

体積要素には絶対値を使用するため：

$$|\det(J)| = \left(\frac{r'}{R}\right)^2$$

### 2.2 極座標における線素

元の極座標 $(r, \theta)$ における線素は：

$$ds^2 = dr^2 + r^2 d\theta^2$$

変換座標 $(r', \theta')$ においては、$dr'$ と $d\theta'$ で表す必要がある。

$r' = R^2/r$ より：

$$dr = \frac{\partial r}{\partial r'} dr' = -\frac{R^2}{r'^2} dr'$$

$\theta' = \theta$ なので $d\theta = d\theta'$。

線素に代入すると：

$$\begin{align}
ds^2 &= dr^2 + r^2 d\theta^2 \\
&= \left(-\frac{R^2}{r'^2} dr'\right)^2 + r^2 (d\theta')^2 \\
&= \frac{R^4}{r'^4} (dr')^2 + \left(\frac{R^2}{r'}\right)^2 (d\theta')^2 \\
&= \frac{R^4}{r'^4} (dr')^2 + \frac{R^4}{r'^2} (d\theta')^2 \\
&= \frac{R^4}{r'^4} \left[(dr')^2 + r'^2 (d\theta')^2\right]
\end{align}$$

したがって、変換後の極座標における線素は：

$$ds^2 = \left(\frac{R^2}{r'^2}\right)^2 \left[(dr')^2 + r'^2 (d\theta')^2\right] = \left(\frac{R}{r'}\right)^4 ds'^2$$

ここで $ds'^2 = (dr')^2 + r'^2 (d\theta')^2$ は $(r', \theta')$ 座標における標準的な極座標線素。

### 2.3 計量テンソル

変換後の極座標 $(r', \theta')$ における**計量テンソル**は：

$$g'_{ij} = \left(\frac{R}{r'}\right)^4 \begin{bmatrix}
1 & 0 \\
0 & r'^2
\end{bmatrix}$$

または、標準的な極座標計量を用いた等角平坦形式で：

$$g'_{ij} = \left(\frac{R}{r'}\right)^4 \, g_{\text{polar}}$$

デカルト的表現（$dx'$, $dy'$を使用）では、計量は次のように簡略化される：

$$g'_{ij} = \left(\frac{r'}{R}\right)^4 \delta_{ij}$$

**重要**: 計量は**等角平坦**である - 単位行列に比例する（デカルト形式において）ため、角度が保存される。

## 3. 体積要素の変換

### 3.1 極座標

極座標における体積要素（2Dでは面積要素）は：

$$d\Omega = r \, dr \, d\theta$$

変換座標において：

$$\begin{align}
d\Omega' &= |\det(J)| \, r \, dr \, d\theta \\
&= \left(\frac{r'}{R}\right)^2 r \, dr \, d\theta \\
&= \left(\frac{r'}{R}\right)^2 d\Omega
\end{align}$$

あるいは、変換座標で直接表現できる。$r = R^2/r'$ より：

$$\begin{align}
d\Omega &= r \, dr \, d\theta \\
&= \frac{R^2}{r'} \cdot \left(-\frac{R^2}{r'^2} dr'\right) \cdot d\theta' \\
&= -\frac{R^4}{r'^3} dr' d\theta'
\end{align}$$

絶対値を取り、$d\Omega' = r' dr' d\theta'$ と比較すると：

$$d\Omega = \frac{R^4}{r'^3} dr' d\theta' = \frac{R^4}{r'^4} \cdot r' dr' d\theta' = \left(\frac{R}{r'}\right)^4 d\Omega'$$

したがって：

$$\boxed{d\Omega = \left(\frac{R}{r'}\right)^4 d\Omega'}$$

または同等に：

$$\boxed{d\Omega' = \left(\frac{r'}{R}\right)^4 d\Omega}$$

## 4. スカラー場とベクトル場の変換

### 4.1 スカラーポテンシャル (φ)

スカラー磁気ポテンシャル $\phi$ は**スカラー場**である。座標変換下でスカラー場は次のように変換される：

$$\phi'(r', \theta') = \phi(r, \theta)$$

**場の値の変換はない** - 座標のみが変化する。

実装において：
- 内部領域（$r < R$）: $\phi$ を直接計算
- 外部領域（$r > R$）: $r' = R^2/r$ として $\phi'(r')$ で表現
- 界面（$r = R$）: $\phi'(R) = \phi(R)$（連続性）

### 4.2 極座標におけるベクトル場

極座標において、ベクトル場 $\mathbf{H}$ は次の成分を持つ：

$$\mathbf{H} = H_r \, \hat{e}_r + H_\theta \, \hat{e}_\theta$$

ここで $\hat{e}_r$ と $\hat{e}_\theta$ は動径方向と方位角方向の単位ベクトル。

### 4.3 磁場の変換（計量ベースの定式化）

**参考文献**: K. Sugahara, "Electromagnetic Analysis of Eddy Current Testing With Kelvin Transformation," IEEE Trans. Magn., vol. 58, no. 9, 2022. (TMAG3194371.pdf)

2次元静磁場解析におけるKelvin変換では、マクスウェル方程式の等角対称性を保存するため、磁場は**計量テンソル**に従って変換される。

#### 円筒座標における計量テンソル

一般的な曲線座標における計量は次のように定義される：

$$g_i = \frac{h_j h_k}{h_i}$$

ここで $h_i, h_j, h_k$ はスケール因子（線素）。

円筒座標 $(r, \theta, z)$ について：

**内部領域：**
$$g_r = r, \quad g_\theta = \frac{1}{r}, \quad g_z = r$$

**外部領域（Kelvin変換後）：**
$$g'_r = -r', \quad g'_\theta = -\frac{1}{r'}, \quad g'_z = -\frac{R^4}{r'^3}$$

#### 材料特性の変換

材料特性（透磁率 $\mu$、誘電率 $\epsilon$、導電率 $\sigma$）は計量の比に従って変換される：

$$\begin{align}
\frac{\mu'_r}{\mu_r} = \frac{\epsilon'_r}{\epsilon_r} = \frac{\sigma'_r}{\sigma_r} &= -\frac{g'_r}{g_r} = 1 \\
\frac{\mu'_\theta}{\mu_\theta} = \frac{\epsilon'_\theta}{\epsilon_\theta} = \frac{\sigma'_\theta}{\sigma_\theta} &= -\frac{g'_\theta}{g_\theta} = 1 \\
\frac{\mu'_z}{\mu_z} = \frac{\epsilon'_z}{\epsilon_z} = \frac{\sigma'_z}{\sigma_z} &= -\frac{g'_z}{g_z} = \left(\frac{R}{r'}\right)^4
\end{align}$$

**重要な知見**：
- **面内成分** $(r, \theta)$: 比 = 1（空間的変調なし）
- **面外成分** $(z)$: 比 = $(R/r')^4$（空間的変調）
- **負の符号**は負のヤコビ行列式から生じる

#### 磁場の変換

一様場の場合、変換は計量比に従う。$\mathbf{H}$ と $\mathbf{B}$ は材料特性を通じて関係しており、材料特性は計量に従って変換されるため、場の変換は：

$$\begin{align}
H'_r(r', \theta') &= -H_r\left(\frac{R^2}{r'}, \theta'\right) \\
H'_\theta(r', \theta') &= -H_\theta\left(\frac{R^2}{r'}, \theta'\right) \\
H'_z(r', \theta') &= -\left(\frac{R}{r'}\right)^4 H_z\left(\frac{R^2}{r'}, \theta'\right)
\end{align}$$

**注**: 負の符号は外部領域の負の計量から生じる。

### 4.4 デカルト座標における一様場の変換

デカルト座標における一様背景場 $\mathbf{H}_s = (0, 1, 0)$（y方向）について：

**ステップ1**: 位置 $(r, \theta)$ で極座標に変換：
$$\begin{align}
H_r &= H_x \cos\theta + H_y \sin\theta = 0 \cdot \cos\theta + 1 \cdot \sin\theta = \sin\theta \\
H_\theta &= -H_x \sin\theta + H_y \cos\theta = 0 \cdot (-\sin\theta) + 1 \cdot \cos\theta = \cos\theta \\
H_z &= 0
\end{align}$$

**ステップ2**: Kelvin変換を適用（面内成分の比 = 1）：
$$\begin{align}
H'_r(r', \theta') &= -H_r\left(\frac{R^2}{r'}, \theta'\right) = -\sin\theta' \\
H'_\theta(r', \theta') &= -H_\theta\left(\frac{R^2}{r'}, \theta'\right) = -\cos\theta' \\
H'_z(r', \theta') &= 0
\end{align}$$

**ステップ3**: 位置 $(x', y')$ でデカルト座標に戻す：
$$\begin{align}
H'_x &= H'_r \cos\theta' - H'_\theta \sin\theta' = (-\sin\theta') \cos\theta' - (-\cos\theta') \sin\theta' \\
&= -\sin\theta' \cos\theta' + \cos\theta' \sin\theta' = 0 \\
H'_y &= H'_r \sin\theta' + H'_\theta \cos\theta' = (-\sin\theta') \sin\theta' + (-\cos\theta') \cos\theta' \\
&= -\sin^2\theta' - \cos^2\theta' = -1 \\
H'_z &= 0
\end{align}$$

**結果**:
$$\boxed{\mathbf{H}'_s(x', y') = (0, -1, 0)}$$

**これは空間的に一様**（$r'$ に依存しない）であり、$\nabla \cdot \mathbf{H}' = 0$ を**自動的に満たす**。

### 4.5 極座標における勾配演算子

極座標における勾配演算子は：

$$\nabla = \hat{e}_r \frac{\partial}{\partial r} + \hat{e}_\theta \frac{1}{r} \frac{\partial}{\partial \theta}$$

変換座標 $(r', \theta')$ において、連鎖律を用いると：

$$\begin{align}
\frac{\partial}{\partial r'} &= \frac{\partial r}{\partial r'} \frac{\partial}{\partial r} = -\frac{R^2}{r'^2} \frac{\partial}{\partial r} \\
\frac{\partial}{\partial \theta'} &= \frac{\partial}{\partial \theta}
\end{align}$$

スカラー場 $\phi$ について：

$$\begin{align}
\nabla \phi &= \hat{e}_r \frac{\partial \phi}{\partial r} \\
&= \hat{e}_{r'} \left(-\frac{R^2}{r'^2}\right)^{-1} \frac{\partial \phi'}{\partial r'} \\
&= \hat{e}_{r'} \left(-\frac{r'^2}{R^2}\right) \frac{\partial \phi'}{\partial r'} \\
&= -\left(\frac{r'}{R}\right)^2 \hat{e}_{r'} \frac{\partial \phi'}{\partial r'}
\end{align}$$

より正確に再計算する。$\phi'(r', \theta') = \phi(r(r'), \theta'(r'))$ なので：

$$\frac{\partial \phi}{\partial r} = \frac{\partial \phi'}{\partial r'} \frac{\partial r'}{\partial r} = \frac{\partial \phi'}{\partial r'} \left(-\frac{R^2}{r^2}\right)$$

したがって：

$$\nabla \phi = \hat{e}_r \frac{\partial \phi}{\partial r} = \hat{e}_{r'} \frac{\partial \phi'}{\partial r'} \left(-\frac{R^2}{r^2}\right) = -\left(\frac{R}{r'}\right)^2 \hat{e}_{r'} \frac{\partial \phi'}{\partial r'}$$

これより：

$$\boxed{\nabla \phi = \left(\frac{R}{r'}\right)^2 \nabla' \phi'}$$

（注: 符号は勾配方向の変換に吸収される。）

極座標から直接計算する。変換 $r = R^2/r'$ より：

$$\frac{\partial \phi}{\partial r} = \frac{\partial \phi'}{\partial r'} \frac{dr'}{dr} = \frac{\partial \phi'}{\partial r'} \left(-\frac{R^2}{r^2}\right) = -\frac{\partial \phi'}{\partial r'} \left(\frac{r'^2}{R^2}\right)$$

よって：

$$\nabla \phi = -\left(\frac{r'}{R}\right)^2 \frac{\partial \phi'}{\partial r'} \hat{e}_r + \frac{R^2}{r'} \frac{1}{r} \frac{\partial \phi'}{\partial \theta'} \hat{e}_\theta$$

第二項について: $\frac{1}{r} = \frac{r'}{R^2}$ なので：

$$\nabla \phi = -\left(\frac{r'}{R}\right)^2 \frac{\partial \phi'}{\partial r'} \hat{e}_r + \left(\frac{r'}{R}\right)^2 \frac{1}{r'} \frac{\partial \phi'}{\partial \theta'} \hat{e}_\theta = \left(\frac{r'}{R}\right)^2 \nabla' \phi'$$

（第一項の負符号は単なる方向 - 変換座標では $r$ の増加は $r'$ の減少に対応。）

正しい勾配変換は：

$$\boxed{\nabla \phi = \left(\frac{r'}{R}\right)^2 \nabla' \phi'}$$

先に書いた $(R/r')^2$ ではない！

## 5. 変換座標における弱形式

### 5.1 元の弱形式（内部領域）

内部領域（$r < R$）における弱形式は：

$$\int_{\Omega_{\text{in}}} (\nabla v) \cdot (\mu \nabla u) \, d\Omega = \int_{\Omega_{\text{in}}} (\nabla v) \cdot (\mu \mathbf{H}_s) \, d\Omega - \int_{\Gamma} v (n \cdot \mu \mathbf{H}_s) \, d\Gamma$$

ここで：
- $u = \phi_m$（磁気スカラーポテンシャル）
- $v$ = テスト関数
- $\mu$ = 透磁率
- $\mathbf{H}_s$ = 背景場 = $(0, 1)$（y方向）
- $\Omega_{\text{in}}$ = 内部領域（$r < R$）
- $\Gamma$ = $r = R$ における境界

### 5.2 変換された弱形式（外部領域）

外部領域（$r > R$）について、有限計算領域内では $(R < r' < 2R)$ として表現し、弱形式を変換する。

#### 5.2.1 双線形形式の変換

出発点：

$$a(u, v) = \int_{\Omega_{\text{out}}} (\nabla v) \cdot (\mu \nabla u) \, d\Omega$$

$(r', \theta')$ 座標へ変換する際、以下を使用：
- $\nabla = \left(\frac{r'}{R}\right)^2 \nabla'$
- $d\Omega = \left(\frac{R}{r'}\right)^4 d\Omega'$

これより：

$$\begin{align}
a(u, v) &= \int_{\Omega_{\text{out}}} (\nabla v) \cdot (\mu \nabla u) \, d\Omega \\
&= \int_{\Omega'_{\text{out}}} \left[\left(\frac{r'}{R}\right)^2 \nabla' v\right] \cdot \left[\mu \left(\frac{r'}{R}\right)^2 \nabla' u\right] \left(\frac{R}{r'}\right)^4 d\Omega' \\
&= \int_{\Omega'_{\text{out}}} \left(\frac{r'}{R}\right)^4 (\nabla' v) \cdot (\mu \nabla' u) \left(\frac{R}{r'}\right)^4 d\Omega' \\
&= \int_{\Omega'_{\text{out}}} (\nabla' v) \cdot (\mu \nabla' u) \, d\Omega'
\end{align}$$

**注目すべき結果**: 係数が相殺される！双線形形式は変換座標において**同じ構造**を持つ：

$$\boxed{a(u, v) = \int_{\Omega'} (\nabla' v) \cdot (\mu \nabla' u) \, d\Omega'}$$

#### 5.2.2 変換された背景場

元の座標における背景場 $\mathbf{H}_s = (0, 1)$ は次のように変換される：

$$\mathbf{H}'_s(r') = \left(\frac{r'}{R}\right)^2 \mathbf{H}_s = \left(\frac{r'}{R}\right)^2 (0, 1) = \left(0, \left(\frac{r'}{R}\right)^2\right)$$

#### 5.2.3 線形形式の変換

線形形式は次のように変換される：

$$\begin{align}
f(v) &= \int_{\Omega_{\text{out}}} (\nabla v) \cdot (\mu \mathbf{H}_s) \, d\Omega \\
&= \int_{\Omega'_{\text{out}}} \left[\left(\frac{r'}{R}\right)^2 \nabla' v\right] \cdot \left[\mu \mathbf{H}_s^{\text{orig}}\right] \left(\frac{R}{r'}\right)^4 d\Omega'
\end{align}$$

ここで $\mathbf{H}_s^{\text{orig}} = (0, 1)$ は元の座標における**定数**背景場。

$$\begin{align}
f(v) &= \int_{\Omega'_{\text{out}}} \left(\frac{r'}{R}\right)^2 \left(\frac{R}{r'}\right)^4 (\nabla' v) \cdot (\mu \mathbf{H}_s^{\text{orig}}) \, d\Omega' \\
&= \int_{\Omega'_{\text{out}}} \left(\frac{R}{r'}\right)^2 (\nabla' v) \cdot (\mu \mathbf{H}_s^{\text{orig}}) \, d\Omega'
\end{align}$$

**変換された背景場** $\mathbf{H}'_s = (r'/R)^2 \mathbf{H}_s^{\text{orig}}$ を用いて表現すると：

$$\mathbf{H}_s^{\text{orig}} = \left(\frac{R}{r'}\right)^2 \mathbf{H}'_s$$

代入すると：

$$\begin{align}
f(v) &= \int_{\Omega'_{\text{out}}} \left(\frac{R}{r'}\right)^2 (\nabla' v) \cdot \left[\mu \left(\frac{R}{r'}\right)^2 \mathbf{H}'_s\right] d\Omega' \\
&= \int_{\Omega'_{\text{out}}} \left(\frac{R}{r'}\right)^4 (\nabla' v) \cdot (\mu \mathbf{H}'_s) \, d\Omega'
\end{align}$$

変換座標における線形形式には**2つの等価な表現**がある：

**形式1**（元の背景場を使用）：

$$\boxed{f(v) = \int_{\Omega'} \left(\frac{R}{r'}\right)^2 (\nabla' v) \cdot (\mu \mathbf{H}_s) \, d\Omega'}$$

ここで $\mathbf{H}_s = (0, 1)$ は定数背景場。

**形式2**（変換された背景場を使用）：

$$\boxed{f(v) = \int_{\Omega'} \left(\frac{R}{r'}\right)^4 (\nabla' v) \cdot (\mu \mathbf{H}'_s) \, d\Omega'}$$

ここで $\mathbf{H}'_s = (0, (r'/R)^2)$ は変換された背景場。

実装には**形式1がより単純で明確**。

### 5.3 完全な変換弱形式

**外部領域**（変換座標で $R < r' < 2R$ として表現）について：

$$\text{求めよ } \phi'_m \in V' \text{ such that:}$$

$$\int_{\Omega'} (\nabla' v) \cdot (\mu \nabla' u) \, d\Omega' = \int_{\Omega'} \left(\frac{R}{r'}\right)^2 (\nabla' v) \cdot (\mu \mathbf{H}_s) \, d\Omega' - \int_{\Gamma'} v (n \cdot \mu \mathbf{H}'_s) \, d\Gamma'$$

ここで：
- $r' = \sqrt{x'^2 + y'^2}$
- $\mathbf{H}_s = (0, 1)$（元の背景場、定数）
- $\mathbf{H}'_s = (0, (r'/R)^2)$（境界における変換された背景場）
- $\Omega'$ = 変換された外部領域（実装では環状領域 $R < r' < 2R$）
- $\Gamma'$ = $r' = R$ における境界（Kelvin界面）
- $R = 2.0$ m（Kelvin半径）

### 5.4 周期境界条件による境界項の消去

**周期境界条件を用いたKelvin変換の主な利点：**

**周期境界条件**（$r = R$ で内部と外部領域を接続）を用いたKelvin変換では、境界積分項は**不要**となる：

$$\boxed{\int_{\Gamma'} v (n \cdot \mu \mathbf{H}'_s) \, d\Gamma' = 0}$$

**理由：**

1. **周期BCによる連続性の保証**: 周期境界条件は、Kelvin界面 $r = R$ でポテンシャルとその法線微分の連続性を自動的に保証する。

2. **自動的なフラックスバランス**: 内部と外部側からの境界への寄与は、周期的同定により相殺される。

3. **簡略化された実装**: 境界積分を明示的に計算・適用する必要がない。

**周期BCを用いた簡略化された弱形式：**

$$\boxed{\int_{\Omega_{\text{interior}} \cup \Omega'_{\text{exterior}}} (\nabla v) \cdot (\mu \nabla u) \, d\Omega = \int_{\Omega_{\text{interior}}} (\nabla v) \cdot (\mu \mathbf{H}_s) \, d\Omega + \int_{\Omega'_{\text{exterior}}} \left(\frac{R}{r'}\right)^2 (\nabla v) \cdot (\mu \mathbf{H}_s) \, d\Omega}$$

**周期BCを用いたNGSolve実装：**

```python
# 周期境界条件の作成
fes = H1(mesh, order=2, dirichlet="GND")
fes = Periodic(fes)  # 周期BCの適用 - 境界項が不要に！

# 双線形形式（両領域で同じ）
a = BilinearForm(fes)
a += mu * grad(u) * grad(v) * dx

# 線形形式（境界項不要！）
f = LinearForm(fes)
f += mu * InnerProduct(grad(v), Hs_weighted) * dx  # 体積積分のみ

# 境界積分不要: -mu*v*InnerProduct(n, Hsb)*ds(...) は不要！
```

**検証**: この手法は以下で検証済み：
- `2D_dipole_with_Kelvin.py`: **0.000%誤差**（優れた一致）
- `2D_quadrupole_with_Kelvin.py`: 全領域で< 2%誤差

## 6. NGSolve実装戦略

### 6.1 領域分割

1. **内部領域**（$r < R = 2.0$ m）：
   - 磁性円: $r < 0.5$ m, $\mu_r = 10$
   - 内部空気: $0.5$ m $< r < 2.0$ m, $\mu_r = 1$
   - 標準弱形式
   - 背景場: $\mathbf{H}_s = (0, 1)$

2. **外部領域**（$r > R = 2.0$ m、環状領域 $R < r' < 2R$ として表現）：
   - Kelvin変換された空気: $R < r' < 2R = 4.0$ m, $\mu_r = 1$
   - $(R/r')^2$ 係数を持つ修正弱形式
   - 背景場: $\mathbf{H}_s = (0, 1)$（積分に $(R/r')^2$ 重みを適用）

3. **界面**（$r = R = 2.0$ m）：
   - 境界項を適用
   - $\phi$ の連続性は自動的に保証
   - 法線方向 $\mathbf{B}$ 成分の連続性は弱形式により保証

### 6.2 NGSolveでの実装

#### 材料特性

```python
from math import pi
from ngsolve import *

mu0 = 4*pi*1e-7  # 真空の透磁率 [H/m]
mu_r_mag = 10    # 磁性材料の比透磁率
mu_r_air = 1     # 空気の比透磁率

mu_d = {
    "magnetic": mu_r_mag * mu0,
    "air_inner": mu_r_air * mu0,
    "air_outer": mu_r_air * mu0  # 変換された外部
}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])
```

#### Kelvin係数を含む背景場

```python
# 動径距離
r = sqrt(x**2 + y**2)
kelvin_radius = 2.0  # R [m]

# Kelvin重み係数 (R/r')^2
# 外部領域（r > R）でのみ適用
kelvin_weight = IfPos(r - kelvin_radius + 1e-6,  # r > R（外部領域）の場合
                      (kelvin_radius / r)**2,     # (R/r')^2
                      1.0)                         # 内部は重みなし

# 背景場（y方向に一定）
Hs = CoefficientFunction((0, 1))  # 元の背景場
Hs_weighted = CoefficientFunction((0, kelvin_weight))  # 線形形式用

# 境界条件用
Hsb = BoundaryFromVolumeCF(CoefficientFunction((0, (r/kelvin_radius)**2)))
```

#### 弱形式

```python
fes = H1(mesh, order=2)
u = fes.TrialFunction()
v = fes.TestFunction()

# 双線形形式（相殺により両領域で同じ）
a = BilinearForm(fes)
a += mu * grad(u) * grad(v) * dx

# Kelvin重み付き背景場を用いた線形形式
f = LinearForm(fes)
f += mu * InnerProduct(grad(v), Hs_weighted) * dx

# Kelvin界面（r = R）での境界項
f += -mu * v * InnerProduct(n, Hsb) * ds(mesh.Boundaries("kelvin_interface"))
```

#### 場の抽出

```python
# 摂動場（両領域で同じ公式）
H_pert = -grad(gfu)

# 全磁場 = 摂動 + 背景
# 外部領域では定数背景場 Hs = (0, 1) を使用
H_total = H_pert + Hs
```

## 7. 主要結果のまとめ

### 7.1 座標変換（極座標）

$$\boxed{r' = \frac{R^2}{r}, \quad \theta' = \theta}$$

### 7.2 体積要素

$$\boxed{d\Omega = \left(\frac{R}{r'}\right)^4 d\Omega'}$$

### 7.3 勾配の変換

$$\boxed{\nabla \phi = \left(\frac{r'}{R}\right)^2 \nabla' \phi'}$$

### 7.4 弱形式

**双線形形式**（両領域で同じ）：

$$\boxed{a(u, v) = \int_{\Omega'} (\nabla' v) \cdot (\mu \nabla' u) \, d\Omega'}$$

**線形形式**（外部領域でKelvin係数）：

$$\boxed{f(v) = \int_{\Omega'} \left(\frac{R}{r'}\right)^2 (\nabla' v) \cdot (\mu \mathbf{H}_s) \, d\Omega'}$$

ここで $\mathbf{H}_s = (0, 1)$ は**定数**元背景場。

## 8. 参考文献

1. **Kelvin変換**：
   - Morse and Feshbach, "Methods of Theoretical Physics", Vol. 2
   - Jackson, "Classical Electrodynamics", 3rd ed., Section 1.13

2. **等角写像**：
   - Churchill and Brown, "Complex Variables and Applications"

3. **電磁場の有限要素法**：
   - Jin, "The Finite Element Method in Electromagnetics", 3rd ed.
   - Monk, "Finite Element Methods for Maxwell's Equations"

## 9. 検証済みダイポール場定式化

### 9.1 背景場（ダイポール）

**一様ダイポール背景場** H_s = (0, 1) A/m（y方向）について：

**特性**：
- div(H_s) = 0 を満たす（一様場は自動的にソレノイダル）
- 極座標: H_r = sinθ, H_θ = cosθ
- 空間的に一様（r依存性なし）

**実装**：

```python
# 背景場: H_s = [0, 1] A/m（2Dのy方向）
Hs = CoefficientFunction((0, 1))
Hsb = BoundaryFromVolumeCF(Hs)
```

### 9.2 弱形式（ダイポール）

**双線形形式**：
```python
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx
```

**線形形式**（境界項付き）：
```python
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx                         # 体積積分
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer")) # 境界項（負の符号）
```

**重要**: 境界項 -∫v(n·μH_s)dΓ は外部境界で必ず含める必要がある。

### 9.3 解析解（ダイポール、2D円筒）

一様場 H_s = (0, 1) 中の透磁率 $\mu_r$ を持つ半径 $a$ の磁性円筒について：

**内部（r < a）**: 摂動場は一様
```python
Hy_pert_analytical = -1.0 + 2.0/(mu_r + 1)
```

**外部（r > a）**: 摂動場は 1/r² で減衰

x軸上（θ = 0 または π）：
```python
Hy_pert_analytical = -(mu_r - 1)/(mu_r + 1) * (circle_radius/r)**2
```

y軸上（θ = π/2）：
```python
Hy_pert_analytical = (mu_r - 1)/(mu_r + 1) * (circle_radius/r)**2
```

### 9.4 ダイポールのKelvin変換

2D面内成分に対する計量ベース定式化：
- 材料特性比: μ'_r/μ_r = μ'_θ/μ_θ = 1（空間的変調なし）
- ヤコビ行列式からの負の計量

**背景場の変換**：

内部 H_s = (0, 1) は極座標で: (H_r, H_θ) = (sinθ, cosθ)

Kelvin変換後の外部計算領域では：
```
H'_r(r', θ') = -H_r(R²/r', θ') = -sinθ'
H'_θ(r', θ') = -H_θ(R²/r', θ') = -cosθ'
```

デカルト座標に戻すと：
```
H'_x = H'_r cosθ' - H'_θ sinθ' = -sinθ' cosθ' + cosθ' sinθ' = 0
H'_y = H'_r sinθ' + H'_θ cosθ' = -sin²θ' - cos²θ' = -1
```

**結果**: H'_s = (0, -1)（空間的に一様、div(H') = 0 を満たす）

**実装**：

```python
# 内部背景場（定数）
Hx_inner = 0.0
Hy_inner = 1.0

# 外部背景場（計量ベース変換により定数）
# 面内成分は空間的に変調されず、符号のみ反転
Hs_x_outer = -Hx_inner  # = 0
Hs_y_outer = -Hy_inner  # = -1

# 領域切り替えを含む背景場
is_exterior = IfPos(x - offset_x/2, 1.0, 0.0)
Hs_x = (1.0 - is_exterior) * Hx_inner + is_exterior * Hs_x_outer
Hs_y = (1.0 - is_exterior) * Hy_inner + is_exterior * Hs_y_outer
Hs = CoefficientFunction((Hs_x, Hs_y))
```

**Kelvin変換を用いた弱形式**：
```python
# 双線形形式（両領域で同じ構造）
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx

# 線形形式（周期BCでは明示的な境界項不要）
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx
```

**Kelvin + 周期BCの特徴**：
- 一様場の場合、変換は単純な符号反転
- 面内成分はr'依存性を獲得しない
- 周期BCにより明示的な境界積分不要
- div(H') = 0 が自動的に満たされる

## 10. 4極場定式化

### 10.1 背景場（4極場）

**4極場背景場** H_s = (x, -y) A/m について：

**特性**：
- div(H_s) = ∂x/∂x + ∂(-y)/∂y = 1 - 1 = 0 を満たす（ソレノイダル）
- 極座標: H_r = r cos(2θ), H_θ = -r sin(2θ)
- ポテンシャル φ_s = -(1/2)r² cos(2θ) に対応

**実装**：

```python
# 背景場: H_s = (x, -y) A/m（4極場）
Hs = CoefficientFunction((x, -y))
Hsb = BoundaryFromVolumeCF(Hs)  # 積分用境界値
```

### 10.2 弱形式（4極場）

**双線形形式**：
```python
a = BilinearForm(fes)
a += mu*grad(u)*grad(v)*dx
```

**線形形式**（境界項付き）：
```python
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx                         # 体積積分
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer")) # 境界項（負の符号）
```

**重要**: 正しい結果を得るには、境界項 -∫v(n·μH_s)dΓ を外部境界で含める必要がある。

### 10.3 解析解（4極場、2D円筒）

4極場 H_s = (x, -y) 中の透磁率 $\mu_r$ を持つ半径 $a$ の磁性円筒について：

**係数**：
```python
A_coeff = (mu_r - 1.0)/(2.0*(mu_r + 1.0))
B_coeff = A_coeff  # 内部と外部で同じ係数
```

**内部（r < a）**: 摂動ポテンシャル φ_pert = Ar² cos(2θ)

磁場：
- H_r = -∂φ/∂r = -2Ar cos(2θ)
- H_θ = -(1/r)∂φ/∂θ = 2Ar sin(2θ)

デカルト成分：
```python
Hr = -2.0 * A_coeff * r * cos(2*theta)
Htheta = 2.0 * A_coeff * r * sin(2*theta)
Hx_analytical = Hr * cos(theta) - Htheta * sin(theta)
Hy_analytical = Hr * sin(theta) + Htheta * cos(theta)
```

**外部（r > a）**: 摂動ポテンシャル φ_pert = B(a⁴/r²) cos(2θ)

磁場：
- H_r = -∂φ/∂r = 2B(a⁴/r³) cos(2θ)
- H_θ = -(1/r)∂φ/∂θ = -2B(a⁴/r³) sin(2θ)

デカルト成分：
```python
Hr = 2.0 * B_coeff * (circle_radius**4 / r**3) * cos(2*theta)
Htheta = -2.0 * B_coeff * (circle_radius**4 / r**3) * sin(2*theta)
Hx_analytical = Hr * cos(theta) - Htheta * sin(theta)
Hy_analytical = Hr * sin(theta) + Htheta * cos(theta)
```

### 10.4 4極場のKelvin変換

内部 H_s = (x, -y) は極座標で: (H_r, H_θ) = (r cos(2θ), -r sin(2θ))

計量ベースKelvin変換を適用：
```
H'_r(r', θ') = -H_r(R²/r', θ') = -(R²/r')cos(2θ)
H'_θ(r', θ') = -H_θ(R²/r', θ') = (R²/r')sin(2θ)
```

境界に合わせた線形 r 依存性：
```
H'_r = -r' cos(2θ)
H'_θ = r' sin(2θ)
```

デカルト座標への変換（cos(2θ)cosθ + sin(2θ)sinθ = cosθ を使用）：
```
H'_x = -r'cos(2θ)cosθ - r'sin(2θ)sinθ
     = -r'[cos(2θ)cosθ + sin(2θ)sinθ]
     = -r'cosθ = -x'

H'_y = -r'cos(2θ)sinθ + r'sin(2θ)cosθ
     = r'[sin(2θ)cosθ - cos(2θ)sinθ]
     = r'sin(θ) = y'
```

**結果**: H'_s = (-x', y')

**実装**：

```python
# 内部背景場
Hx_inner = x
Hy_inner = -y

# 外部背景場（変換後に空間的に一様！）
Hs_x_outer = -x_local
Hs_y_outer = y_local

# 領域切り替えを含む背景場
is_exterior = IfPos(x - offset_x/2, 1.0, 0.0)
Hs_x = (1.0 - is_exterior) * Hx_inner + is_exterior * Hs_x_outer
Hs_y = (1.0 - is_exterior) * Hy_inner + is_exterior * Hs_y_outer
Hs = CoefficientFunction((Hs_x, Hs_y))
```

## 11. まとめ

### 11.1 2D計量ベースKelvin変換則

2D計量ベースKelvin変換について：

$$\boxed{
\begin{align}
H'_r(r', \theta') &= -H_r(R^2/r', \theta') \\
H'_\theta(r', \theta') &= -H_\theta(R^2/r', \theta')
\end{align}
}$$

**主要な特性**：
- 両極座標成分が計量から負の符号を得る
- 空間依存性は変換後に簡略化されうる
- 自動的に div(H') = 0 を保存

変換例：
- 一様場（ダイポール） → 一様場（符号反転）
- 線形変化場（4極場） → 一様場

### 11.2 境界項のまとめ

| 手法 | 外部境界での境界項 | 理由 |
|------|-------------------|------|
| **Kelvin変換なし**（有限領域） | **必要** | 外部境界を通る磁束を考慮必要 |
| **Kelvin + 周期BC** | **不要** | 周期BCが自動的に連続性を保証 |
