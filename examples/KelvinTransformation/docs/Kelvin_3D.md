# 3次元静磁場解析のためのKelvin変換

## 概要

本文書は、3次元H-formulation静磁場解析におけるKelvin変換の完全な数学的基礎と実装ガイドを提供する。Kelvin変換は外部領域$(r > R)$を有限計算領域に写像し、無限領域問題の効率的な表現を可能にする。

本文書の内容：
- Kelvin変換の数学理論
- 強形式から弱形式への導出
- NGSolve実装の詳細
- ダイポール場および4極場の適用例

## 0. 摂動ポテンシャル定式化

### 0.1 強形式

磁束保存則（$\nabla \cdot \mathbf{B} = 0$）から出発する：

$$\mathbf{H}_{\text{total}} = \mathbf{H}_s + \mathbf{H}_{\text{pert}} = \mathbf{H}_s - \nabla\phi_m$$

$$\mathbf{B} = \mu\mathbf{H}_{\text{total}} = \mu(\mathbf{H}_s - \nabla\phi_m)$$

$$\nabla \cdot \mathbf{B} = 0$$

$$\Rightarrow \nabla \cdot [\mu(\mathbf{H}_s - \nabla\phi_m)] = 0$$

$$\Rightarrow \nabla \cdot (\mu\mathbf{H}_s) - \nabla \cdot (\mu\nabla\phi_m) = 0$$

$$\Rightarrow -\nabla \cdot (\mu\nabla\phi_m) = -\nabla \cdot (\mu\mathbf{H}_s)$$

**境界条件**：遠方界（$r \to \infty$）において：
$$\mathbf{H}_{\text{pert}} \to 0 \quad \Rightarrow \quad \nabla\phi_m \to 0 \quad \text{（自然境界条件: } n \cdot \nabla\phi_m = 0\text{）}$$

### 0.2 弱形式の導出

出発点：$-\nabla \cdot (\mu\nabla\phi_m) = -\nabla \cdot (\mu\mathbf{H}_s)$

**ステップ1**：試験関数$v$を掛けて積分：
$$-\int_\Omega v\nabla \cdot (\mu\nabla\phi_m) d\Omega = -\int_\Omega v\nabla \cdot (\mu\mathbf{H}_s) d\Omega$$

**ステップ2**：部分積分（発散定理）を適用：

左辺：
$$-\int_\Gamma v(n \cdot \mu\nabla\phi_m) d\Gamma + \int_\Omega (\nabla v) \cdot (\mu\nabla\phi_m) d\Omega$$

右辺：
$$-\int_\Gamma v(n \cdot \mu\mathbf{H}_s) d\Gamma + \int_\Omega (\nabla v) \cdot (\mu\mathbf{H}_s) d\Omega$$

**ステップ3**：自然境界条件（$n \cdot \nabla\phi_m = 0$、$\Gamma$上）を適用：
$$-\int_\Gamma v(n \cdot \mu\nabla\phi_m) d\Gamma = 0$$

**ステップ4**：最終弱形式：

$\phi_m \in V$を求めよ：

$$a(\phi_m, v) = f(v) \quad \forall v \in V$$

ここで：
$$a(u,v) = \int_\Omega (\nabla v) \cdot (\mu\nabla u) d\Omega \quad \text{（双線形形式）}$$

$$f(v) = \int_\Omega (\nabla v) \cdot (\mu\mathbf{H}_s) d\Omega - \int_\Gamma v(n \cdot \mu\mathbf{H}_s) d\Gamma \quad \text{（線形形式）}$$

**重要**：体積積分は**正**の符号、境界積分は**負**の符号であることに注意。

### 0.3 有限領域問題の境界項の取り扱い（Kelvin変換なし）

**重要**：Kelvin変換を使用しない問題（外部境界を持つ有限領域）では、**境界積分項を線形形式に含める必要がある**：

$$f(v) = \int_\Omega (\nabla v) \cdot (\mu\mathbf{H}_s) d\Omega - \int_{\Gamma_{\text{outer}}} v(n \cdot \mu\mathbf{H}_s) d\Gamma$$

#### 境界項の起源

境界項は、弱形式の導出における**部分積分**から生じる（上記ステップ2）：

$$-\int_\Omega v\nabla \cdot (\mu\mathbf{H}_s) d\Omega = -\int_\Gamma v(n \cdot \mu\mathbf{H}_s) d\Gamma + \int_\Omega (\nabla v) \cdot (\mu\mathbf{H}_s) d\Omega$$

この境界項は自然境界条件では**消去できない**。理由：
- 自然境界条件（$n \cdot \nabla\phi_m = 0$）は**左辺**（双線形形式から）の境界項を消去する
- **右辺**（線形形式から）の境界項は$\phi_m$ではなく$\mathbf{H}_s$を含むため残る

#### 物理的解釈

1. **磁束の寄与**：項$-\int_\Gamma v(n \cdot \mu\mathbf{H}_s) d\Gamma$は、計算領域境界を通過する背景磁束の寄与を表す。

2. **遠方界の挙動**：無限領域では、摂動場が無限遠で減衰するため、この項はゼロに積分される。しかし、**切断された有限領域**では、外部場の影響を正しくモデル化するためにこの項が不可欠である。

3. **エネルギー収支**：この項がなければ、境界での背景場による仕事が考慮されないため、FEM定式化はエネルギー保存則に違反する。

#### 境界項が不可欠な理由

1. **数学的完全性**：弱形式は部分積分を介して強形式から導出される。境界項を省略するとこの導出が無効になる。

2. **正しい遠方界の影響**：境界項は、摂動場が領域境界での背景場の影響を正しく考慮することを保証する。

3. **検証**：ファイル`3D_dipole.py`および`3D_quadrupole.py`は、正しい結果を得るために境界項が不可欠であることを実証している。

**NGSolve実装**（Kelvin変換なしの標準形式）：

```python
# 線形形式: f(v) = ∫(∇v)·(μH_s)dΩ - ∫v(n·μH_s)dΓ
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx                         # 体積積分（正）
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer")) # 境界項（負）
```

**重要な符号の規則**：
- 体積積分：**正**の符号（`+`）
- 境界積分：**負**の符号（`-`）

#### Kelvin変換との対比

| 手法 | 境界項 | 理由 |
|------|--------|------|
| **有限領域（Kelvin変換なし）** | **必要** | 有限距離での境界が領域を切断 |
| **Kelvin + 周期境界条件** | **不要** | 周期BCが$r = R$での連続性を保証 |

**Kelvin + 周期BCが境界項を消去する理由：**

1. **領域の閉包**：Kelvin変換は$r > R$を有限領域に写像し、周期BCは内部（$r < R$）と外部（$r > R$）領域を$r = R$で接続する。

2. **自動的な連続性**：周期BCは自動的に以下を強制する：
   - ポテンシャル連続性：$\phi(R^-) = \phi'(R^+)$
   - フラックス連続性：$(n \cdot \mu\nabla\phi)|_{R^-} = (n \cdot \mu'\nabla'\phi')|_{R^+}$

3. **境界寄与の相殺**：$r = R$の内側と外側からの境界積分寄与は、法線方向が逆であることと周期的な同一視により相殺する。

### 0.4 NGSolve実装（標準形式）

```python
from math import pi
from ngsolve import *

# 材料特性
mu0 = 4*pi*1e-7  # 真空の透磁率 [H/m]
mu_d = {"air_inner": 1*mu0, "air_outer": 1*mu0, "magnetic": 10*mu0}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])

# 背景磁場（3D: z方向）
Hs = CoefficientFunction((0, 0, 1))
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

3次元における**Kelvin変換**（球面反転とも呼ばれる）は以下で定義される：

$$r' = \frac{R^2}{r}$$

ここで：
- $r = \sqrt{x^2 + y^2 + z^2}$：元の座標における動径距離
- $r' = \sqrt{x'^2 + y'^2 + z'^2}$：変換後の座標における動径距離
- $R$：変換半径（Kelvin球半径）

### 1.2 球座標系

3次元問題では、**球座標**$(r, \theta, \phi)$を使用する：

$$\begin{align}
x &= r \sin\theta \cos\phi \\
y &= r \sin\theta \sin\phi \\
z &= r \cos\theta
\end{align}$$

ここで：
- $\theta$：極角（$0 \leq \theta \leq \pi$）
- $\phi$：方位角（$0 \leq \phi < 2\pi$）

球座標におけるKelvin変換は：

$$\begin{align}
r' &= \frac{R^2}{r} \\
\theta' &= \theta \quad \text{（極角は保存）} \\
\phi' &= \phi \quad \text{（方位角は保存）}
\end{align}$$

### 1.3 直交座標変換

$r' = R^2/r$、$\theta' = \theta$、$\phi' = \phi$より、直交座標は以下のように変換される：

$$\begin{align}
x' &= r' \sin\theta' \cos\phi' = \frac{R^2}{r} \sin\theta \cos\phi = \frac{R^2}{r^2} x = \frac{R^2 x}{x^2 + y^2 + z^2} \\
y' &= r' \sin\theta' \sin\phi' = \frac{R^2}{r} \sin\theta \sin\phi = \frac{R^2}{r^2} y = \frac{R^2 y}{x^2 + y^2 + z^2} \\
z' &= r' \cos\theta' = \frac{R^2}{r} \cos\theta = \frac{R^2}{r^2} z = \frac{R^2 z}{x^2 + y^2 + z^2}
\end{align}$$

**逆変換**：

$$\begin{align}
x &= \frac{R^2}{r'^2} x' = \frac{R^2 x'}{x'^2 + y'^2 + z'^2} \\
y &= \frac{R^2}{r'^2} y' = \frac{R^2 y'}{x'^2 + y'^2 + z'^2} \\
z &= \frac{R^2}{r'^2} z' = \frac{R^2 z'}{x'^2 + y'^2 + z'^2}
\end{align}$$

### 1.4 主要な性質

1. **等角写像**：角度が（局所的に）保存される
2. **対合性**：変換を2回適用すると元の座標に戻る
3. **球面の保存**：球面は球面に（または原点を通る平面に）写像される
4. **領域の写像**：
   - $r \to \infty$は$r' \to 0$に写像される
   - $r = R$は$r' = R$に写像される（不動球面）
   - $r < R$は$r' > R$に写像される

## 2. 球座標における計量テンソル

### 2.1 球座標形式での座標変換

球座標では、変換は：

$$(r, \theta, \phi) \to (r', \theta', \phi') = \left(\frac{R^2}{r}, \theta, \phi\right)$$

この変換のヤコビ行列は：

$$J = \begin{bmatrix}
\frac{\partial r'}{\partial r} & \frac{\partial r'}{\partial \theta} & \frac{\partial r'}{\partial \phi} \\
\frac{\partial \theta'}{\partial r} & \frac{\partial \theta'}{\partial \theta} & \frac{\partial \theta'}{\partial \phi} \\
\frac{\partial \phi'}{\partial r} & \frac{\partial \phi'}{\partial \theta} & \frac{\partial \phi'}{\partial \phi}
\end{bmatrix}$$

偏微分の計算：

$$\begin{align}
\frac{\partial r'}{\partial r} &= \frac{\partial}{\partial r}\left(\frac{R^2}{r}\right) = -\frac{R^2}{r^2} \\
\frac{\partial r'}{\partial \theta} &= 0, \quad \frac{\partial r'}{\partial \phi} = 0 \\
\frac{\partial \theta'}{\partial r} &= 0, \quad \frac{\partial \theta'}{\partial \theta} = 1, \quad \frac{\partial \theta'}{\partial \phi} = 0 \\
\frac{\partial \phi'}{\partial r} &= 0, \quad \frac{\partial \phi'}{\partial \theta} = 0, \quad \frac{\partial \phi'}{\partial \phi} = 1
\end{align}$$

したがって：

$$J = \begin{bmatrix}
-\frac{R^2}{r^2} & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix}$$

行列式は：

$$\det(J) = -\frac{R^2}{r^2} = -\left(\frac{r'}{R}\right)^2$$

体積要素では絶対値を考慮するため：

$$|\det(J)| = \left(\frac{r'}{R}\right)^2$$

### 2.2 球座標における線素

元の球座標$(r, \theta, \phi)$における線素は：

$$ds^2 = dr^2 + r^2 d\theta^2 + r^2 \sin^2\theta \, d\phi^2$$

変換後の座標$(r', \theta', \phi')$では、$dr'$、$d\theta'$、$d\phi'$で表す必要がある。

$r' = R^2/r$より：

$$dr = \frac{\partial r}{\partial r'} dr' = -\frac{R^2}{r'^2} dr'$$

$\theta' = \theta$、$\phi' = \phi$より$d\theta = d\theta'$、$d\phi = d\phi'$。

線素に代入すると：

$$\begin{align}
ds^2 &= dr^2 + r^2 d\theta^2 + r^2 \sin^2\theta \, d\phi^2 \\
&= \left(-\frac{R^2}{r'^2} dr'\right)^2 + r^2 (d\theta')^2 + r^2 \sin^2\theta' \, (d\phi')^2 \\
&= \frac{R^4}{r'^4} (dr')^2 + \left(\frac{R^2}{r'}\right)^2 (d\theta')^2 + \left(\frac{R^2}{r'}\right)^2 \sin^2\theta' \, (d\phi')^2 \\
&= \frac{R^4}{r'^4} (dr')^2 + \frac{R^4}{r'^2} (d\theta')^2 + \frac{R^4}{r'^2} \sin^2\theta' \, (d\phi')^2 \\
&= \frac{R^4}{r'^4} \left[(dr')^2 + r'^2 (d\theta')^2 + r'^2 \sin^2\theta' \, (d\phi')^2\right]
\end{align}$$

したがって、変換後の球座標における線素は：

$$ds^2 = \left(\frac{R^2}{r'^2}\right)^2 \left[(dr')^2 + r'^2 (d\theta')^2 + r'^2 \sin^2\theta' \, (d\phi')^2\right] = \left(\frac{R}{r'}\right)^4 ds'^2$$

ここで$ds'^2 = (dr')^2 + r'^2 (d\theta')^2 + r'^2 \sin^2\theta' \, (d\phi')^2$は$(r', \theta', \phi')$座標における標準的な球面線素。

### 2.3 計量テンソル

変換後の球座標$(r', \theta', \phi')$における**計量テンソル**は：

$$g'_{ij} = \left(\frac{R}{r'}\right)^4 \begin{bmatrix}
1 & 0 & 0 \\
0 & r'^2 & 0 \\
0 & 0 & r'^2 \sin^2\theta'
\end{bmatrix}$$

または、等価な共形平坦形式で：

$$g'_{ij} = \left(\frac{R}{r'}\right)^4 \, g_{\text{spherical}}$$

直交座標風の表現（$dx'$、$dy'$、$dz'$）では、計量は簡略化される：

$$g'_{ij} = \left(\frac{r'}{R}\right)^4 \delta_{ij}$$

**重要**：計量は**共形平坦**である - （直交座標形式で）単位行列に比例し、角度が保存されることを意味する。

## 3. 体積要素の変換

### 3.1 球座標

球座標における体積要素は：

$$d\Omega = r^2 \sin\theta \, dr \, d\theta \, d\phi$$

変換後の座標では：

$$\begin{align}
d\Omega' &= |\det(J)| \, r^2 \sin\theta \, dr \, d\theta \, d\phi \\
&= \left(\frac{r'}{R}\right)^2 r^2 \sin\theta \, dr \, d\theta \, d\phi \\
&= \left(\frac{r'}{R}\right)^2 d\Omega
\end{align}$$

あるいは、変換後の座標で直接表すこともできる。$r = R^2/r'$より：

$$\begin{align}
d\Omega &= r^2 \sin\theta \, dr \, d\theta \, d\phi \\
&= \left(\frac{R^2}{r'}\right)^2 \sin\theta' \cdot \left(-\frac{R^2}{r'^2} dr'\right) \cdot d\theta' \, d\phi' \\
&= -\frac{R^6}{r'^4} \sin\theta' \, dr' \, d\theta' \, d\phi'
\end{align}$$

絶対値を取り、$d\Omega' = r'^2 \sin\theta' \, dr' \, d\theta' \, d\phi'$と比較すると：

$$d\Omega = \frac{R^6}{r'^4} \sin\theta' \, dr' \, d\theta' \, d\phi' = \frac{R^6}{r'^6} \cdot r'^2 \sin\theta' \, dr' \, d\theta' \, d\phi' = \left(\frac{R}{r'}\right)^6 d\Omega'$$

したがって：

$$\boxed{d\Omega = \left(\frac{R}{r'}\right)^6 d\Omega'}$$

または同等に：

$$\boxed{d\Omega' = \left(\frac{r'}{R}\right)^6 d\Omega}$$

**注意**：3次元では、体積要素は$(R/r')^6$でスケールし、2次元の$(R/r')^4$とは異なる。

## 4. スカラー場とベクトル場の変換

### 4.1 スカラーポテンシャル (φ)

磁気スカラーポテンシャル$\phi$は**スカラー場**である。座標変換下でスカラー場は以下のように変換される：

$$\phi'(r', \theta', \phi') = \phi(r, \theta, \phi)$$

**場の値は変換されない** - 座標のみが変化する。

実装において：
- 内部領域$(r < R)$：$\phi$を直接計算
- 外部領域$(r > R)$：$\phi'(r')$として表現、ここで$r' = R^2/r$
- 界面$(r = R)$：$\phi'(R) = \phi(R)$（連続性）

### 4.2 球座標におけるベクトル場

球座標において、ベクトル場$\mathbf{H}$は成分を持つ：

$$\mathbf{H} = H_r \, \hat{e}_r + H_\theta \, \hat{e}_\theta + H_\phi \, \hat{e}_\phi$$

ここで$\hat{e}_r$、$\hat{e}_\theta$、$\hat{e}_\phi$は動径方向、極方向、方位方向の単位ベクトル。

### 4.3 磁場の変換（計量ベースの定式化）

**参考文献**：K. Sugahara, "Electromagnetic Analysis of Eddy Current Testing With Kelvin Transformation," IEEE Trans. Magn., vol. 58, no. 9, 2022. (TMAG3194371.pdf)

3次元静磁場のKelvin変換において、磁場はMaxwell方程式の共形対称性を保存するため、**計量テンソル**に従って変換される。

#### 球座標における計量テンソル

一般曲線座標における計量は以下で定義される：

$$g_i = \frac{h_j h_k}{h_i}$$

ここで$h_i, h_j, h_k$はスケール因子（線素）。

球座標$(r, \theta, \phi)$において：

**内部領域：**
$$g_r = r^2 \sin\theta, \quad g_\theta = \sin\theta, \quad g_\phi = \frac{1}{\sin\theta}$$

**外部領域（Kelvin変換後）：**
$$g'_r = -R^2 \sin\theta, \quad g'_\theta = -\frac{R^2}{r'^2} \sin\theta, \quad g'_\phi = -\frac{R^2}{r'^2 \sin\theta}$$

#### 材料特性の変換

材料特性（透磁率$\mu$、誘電率$\epsilon$、導電率$\sigma$）は計量比に従って変換される：

$$\begin{align}
\frac{\mu'_r}{\mu_r} = \frac{\epsilon'_r}{\epsilon_r} = \frac{\sigma'_r}{\sigma_r} &= -\frac{g'_r}{g_r} = \left(\frac{R}{r'}\right)^2 \\
\frac{\mu'_\theta}{\mu_\theta} = \frac{\epsilon'_\theta}{\epsilon_\theta} = \frac{\sigma'_\theta}{\sigma_\theta} &= -\frac{g'_\theta}{g_\theta} = \left(\frac{R}{r'}\right)^2 \\
\frac{\mu'_\phi}{\mu_\phi} = \frac{\epsilon'_\phi}{\epsilon_\phi} = \frac{\sigma'_\phi}{\sigma_\phi} &= -\frac{g'_\phi}{g_\phi} = \left(\frac{R}{r'}\right)^2
\end{align}$$

**3次元の重要な洞察**：
- **すべての成分**$(r, \theta, \phi)$：比 = $(R/r')^2$（空間変調）
- 負のヤコビアン行列式からの**負の符号**
- **2次元との違い**：2次元では面内成分の比は1だった；3次元ではすべての成分が一様にスケール

#### 磁場の変換

一様場の場合、変換は計量比に従う：

$$\begin{align}
H'_r(r', \theta', \phi') &= -\left(\frac{R}{r'}\right)^2 H_r\left(\frac{R^2}{r'}, \theta', \phi'\right) \\
H'_\theta(r', \theta', \phi') &= -\left(\frac{R}{r'}\right)^2 H_\theta\left(\frac{R^2}{r'}, \theta', \phi'\right) \\
H'_\phi(r', \theta', \phi') &= -\left(\frac{R}{r'}\right)^2 H_\phi\left(\frac{R^2}{r'}, \theta', \phi'\right)
\end{align}$$

**注意**：負の符号は外部領域における負の計量から生じる。

### 4.4 直交座標における一様場の変換

直交座標における一様背景場$\mathbf{H}_s = (0, 0, 1)$（z方向）の場合：

**ステップ1**：位置$(r, \theta, \phi)$で球座標に変換：
$$\begin{align}
H_r &= H_x \sin\theta \cos\phi + H_y \sin\theta \sin\phi + H_z \cos\theta \\
&= 0 \cdot \sin\theta \cos\phi + 0 \cdot \sin\theta \sin\phi + 1 \cdot \cos\theta = \cos\theta \\
H_\theta &= H_x \cos\theta \cos\phi + H_y \cos\theta \sin\phi - H_z \sin\theta \\
&= 0 \cdot \cos\theta \cos\phi + 0 \cdot \cos\theta \sin\phi - 1 \cdot \sin\theta = -\sin\theta \\
H_\phi &= -H_x \sin\phi + H_y \cos\phi = 0 \cdot (-\sin\phi) + 0 \cdot \cos\phi = 0
\end{align}$$

**ステップ2**：Kelvin変換を適用（すべての成分が$(R/r')^2$でスケール）：
$$\begin{align}
H'_r(r', \theta', \phi') &= -\left(\frac{R}{r'}\right)^2 H_r\left(\frac{R^2}{r'}, \theta', \phi'\right) = -\left(\frac{R}{r'}\right)^2 \cos\theta' \\
H'_\theta(r', \theta', \phi') &= -\left(\frac{R}{r'}\right)^2 H_\theta\left(\frac{R^2}{r'}, \theta', \phi'\right) = -\left(\frac{R}{r'}\right)^2 (-\sin\theta') = \left(\frac{R}{r'}\right)^2 \sin\theta' \\
H'_\phi(r', \theta', \phi') &= 0
\end{align}$$

**ステップ3**：位置$(x', y', z')$で直交座標に変換：
$$\begin{align}
H'_x &= H'_r \sin\theta' \cos\phi' + H'_\theta \cos\theta' \cos\phi' - H'_\phi \sin\phi' \\
&= -\left(\frac{R}{r'}\right)^2 \cos\theta' \sin\theta' \cos\phi' + \left(\frac{R}{r'}\right)^2 \sin\theta' \cos\theta' \cos\phi' - 0 \\
&= 0 \\
H'_y &= H'_r \sin\theta' \sin\phi' + H'_\theta \cos\theta' \sin\phi' + H'_\phi \cos\phi' \\
&= -\left(\frac{R}{r'}\right)^2 \cos\theta' \sin\theta' \sin\phi' + \left(\frac{R}{r'}\right)^2 \sin\theta' \cos\theta' \sin\phi' + 0 \\
&= 0 \\
H'_z &= H'_r \cos\theta' - H'_\theta \sin\theta' \\
&= -\left(\frac{R}{r'}\right)^2 \cos^2\theta' - \left(\frac{R}{r'}\right)^2 \sin^2\theta' \\
&= -\left(\frac{R}{r'}\right)^2 (\cos^2\theta' + \sin^2\theta') = -\left(\frac{R}{r'}\right)^2
\end{align}$$

**結果**：
$$\boxed{\mathbf{H}'_s(x', y', z') = \left(0, 0, -\left(\frac{R}{r'}\right)^2\right)}$$

**これは空間変動**を持ち$(R/r')^2$に比例するが、球対称性により$\nabla \cdot \mathbf{H}' = 0$を**満たす**。

**発散フリー条件の検証**：
球座標において：
$$\nabla \cdot \mathbf{H}' = \frac{1}{r'^2} \frac{\partial (r'^2 H'_r)}{\partial r'} + \frac{1}{r' \sin\theta'} \frac{\partial (\sin\theta' H'_\theta)}{\partial \theta'} + \frac{1}{r' \sin\theta'} \frac{\partial H'_\phi}{\partial \phi'}$$

$H'_r = -\frac{R^2}{r'^2} \cos\theta'$、$H'_\theta = \frac{R^2}{r'^2} \sin\theta'$、$H'_\phi = 0$で：
$$\nabla \cdot \mathbf{H}' = \frac{1}{r'^2} \frac{\partial (r'^2 \cdot (-\frac{R^2}{r'^2} \cos\theta'))}{\partial r'} + \frac{1}{r' \sin\theta'} \frac{\partial (\sin\theta' \cdot \frac{R^2}{r'^2} \sin\theta')}{\partial \theta'} = 0 + 0 = 0$$

### 4.5 球座標における勾配演算子

球座標における勾配演算子は：

$$\nabla = \hat{e}_r \frac{\partial}{\partial r} + \hat{e}_\theta \frac{1}{r} \frac{\partial}{\partial \theta} + \hat{e}_\phi \frac{1}{r\sin\theta} \frac{\partial}{\partial \phi}$$

変換後の座標$(r', \theta', \phi')$では、連鎖律を使用して：

$$\begin{align}
\frac{\partial}{\partial r'} &= \frac{\partial r}{\partial r'} \frac{\partial}{\partial r} = -\frac{R^2}{r'^2} \frac{\partial}{\partial r} \\
\frac{\partial}{\partial \theta'} &= \frac{\partial}{\partial \theta} \\
\frac{\partial}{\partial \phi'} &= \frac{\partial}{\partial \phi}
\end{align}$$

スカラー場$\phi$について、$\phi'(r', \theta', \phi') = \phi(r(r'), \theta'(r'), \phi'(r'))$であるから：

$$\frac{\partial \phi}{\partial r} = \frac{\partial \phi'}{\partial r'} \frac{\partial r'}{\partial r} = \frac{\partial \phi'}{\partial r'} \left(-\frac{R^2}{r^2}\right)$$

したがって、勾配の動径成分は：

$$\nabla \phi \cdot \hat{e}_r = \frac{\partial \phi}{\partial r} = \frac{\partial \phi'}{\partial r'} \left(-\frac{R^2}{r^2}\right) = -\left(\frac{R}{r'}\right)^2 \frac{\partial \phi'}{\partial r'}$$

角成分については：
$$\begin{align}
\frac{1}{r}\frac{\partial \phi}{\partial \theta} &= \frac{r'}{R^2} \frac{1}{r'}\frac{\partial \phi'}{\partial \theta'} = \frac{r'}{R^2} \frac{1}{r'}\frac{\partial \phi'}{\partial \theta'} \\
\frac{1}{r\sin\theta}\frac{\partial \phi}{\partial \phi} &= \frac{r'}{R^2} \frac{1}{r'\sin\theta'}\frac{\partial \phi'}{\partial \phi'}
\end{align}$$

すべての成分を組み合わせ（符号規約を考慮して）：

$$\boxed{\nabla \phi = \left(\frac{r'}{R}\right)^2 \nabla' \phi'}$$

これは2次元と同じスケーリング！

## 5. 変換後座標における弱形式

### 5.1 元の弱形式（内部領域）

内部領域$(r < R)$において、弱形式は：

$$\int_{\Omega_{\text{in}}} (\nabla v) \cdot (\mu \nabla u) \, d\Omega = \int_{\Omega_{\text{in}}} (\nabla v) \cdot (\mu \mathbf{H}_s) \, d\Omega - \int_{\Gamma} v (n \cdot \mu \mathbf{H}_s) \, d\Gamma$$

ここで：
- $u = \phi_m$（磁気スカラーポテンシャル）
- $v$ = 試験関数
- $\mu$ = 透磁率
- $\mathbf{H}_s$ = 背景磁場 = z方向の場合$(0, 0, 1)$
- $\Omega_{\text{in}}$ = 内部領域$(r < R)$
- $\Gamma$ = $r = R$での境界

### 5.2 変換後の弱形式（外部領域）

外部領域$(r > R)$について、有限計算において変換後領域$(R < r' < 2R)$として表現し、弱形式を変換する。

#### 5.2.1 双線形形式の変換

出発点：

$$a(u, v) = \int_{\Omega_{\text{out}}} (\nabla v) \cdot (\mu \nabla u) \, d\Omega$$

$(r', \theta', \phi')$座標に変換する。以下を使用：
- $\nabla = \left(\frac{r'}{R}\right)^2 \nabla'$
- $d\Omega = \left(\frac{R}{r'}\right)^6 d\Omega'$

次を得る：

$$\begin{align}
a(u, v) &= \int_{\Omega_{\text{out}}} (\nabla v) \cdot (\mu \nabla u) \, d\Omega \\
&= \int_{\Omega'_{\text{out}}} \left[\left(\frac{r'}{R}\right)^2 \nabla' v\right] \cdot \left[\mu \left(\frac{r'}{R}\right)^2 \nabla' u\right] \left(\frac{R}{r'}\right)^6 d\Omega' \\
&= \int_{\Omega'_{\text{out}}} \left(\frac{r'}{R}\right)^4 (\nabla' v) \cdot (\mu \nabla' u) \left(\frac{R}{r'}\right)^6 d\Omega' \\
&= \int_{\Omega'_{\text{out}}} \left(\frac{R}{r'}\right)^2 (\nabla' v) \cdot (\mu \nabla' u) \, d\Omega'
\end{align}$$

**2次元との重要な違い**：3次元では因子が**完全には相殺しない**！残余の$(R/r')^2$因子がある。

$$\boxed{a(u, v) = \int_{\Omega'} \left(\frac{R}{r'}\right)^2 (\nabla' v) \cdot (\mu \nabla' u) \, d\Omega'}$$

#### 5.2.2 変換後の背景磁場

元の座標における背景磁場$\mathbf{H}_s = (0, 0, 1)$は以下のように変換される：

$$\mathbf{H}'_s(r') = \left(\frac{r'}{R}\right)^2 \mathbf{H}_s = \left(\frac{r'}{R}\right)^2 (0, 0, 1) = \left(0, 0, \left(\frac{r'}{R}\right)^2\right)$$

#### 5.2.3 線形形式の変換

線形形式は以下のように変換される：

$$\begin{align}
f(v) &= \int_{\Omega_{\text{out}}} (\nabla v) \cdot (\mu \mathbf{H}_s) \, d\Omega \\
&= \int_{\Omega'_{\text{out}}} \left[\left(\frac{r'}{R}\right)^2 \nabla' v\right] \cdot \left[\mu \mathbf{H}_s^{\text{orig}}\right] \left(\frac{R}{r'}\right)^6 d\Omega'
\end{align}$$

ここで$\mathbf{H}_s^{\text{orig}} = (0, 0, 1)$は元の座標における**定数**背景磁場。

$$\begin{align}
f(v) &= \int_{\Omega'_{\text{out}}} \left(\frac{r'}{R}\right)^2 \left(\frac{R}{r'}\right)^6 (\nabla' v) \cdot (\mu \mathbf{H}_s^{\text{orig}}) \, d\Omega' \\
&= \int_{\Omega'_{\text{out}}} \left(\frac{R}{r'}\right)^4 (\nabla' v) \cdot (\mu \mathbf{H}_s^{\text{orig}}) \, d\Omega'
\end{align}$$

**線形形式**（外部領域にKelvin因子付き）：

$$\boxed{f(v) = \int_{\Omega'} \left(\frac{R}{r'}\right)^4 (\nabla' v) \cdot (\mu \mathbf{H}_s) \, d\Omega'}$$

ここで$\mathbf{H}_s = (0, 0, 1)$は**定数**の元の背景磁場。

または、変換後の背景磁場$\mathbf{H}'_s = (0, 0, (r'/R)^2)$を使用して：

$$\mathbf{H}_s^{\text{orig}} = \left(\frac{R}{r'}\right)^2 \mathbf{H}'_s$$

代入すると：

$$\boxed{f(v) = \int_{\Omega'} \left(\frac{R}{r'}\right)^6 (\nabla' v) \cdot (\mu \mathbf{H}'_s) \, d\Omega'}$$

実装においては、**$\mathbf{H}_s$を使用する最初の形式の方がシンプルで明確**。

### 5.3 完全な変換後弱形式

**外部領域**（変換後座標で$R < r' < 2R$として表現）について：

$$\text{$\phi'_m \in V'$を求めよ：}$$

$$\int_{\Omega'} \left(\frac{R}{r'}\right)^2 (\nabla' v) \cdot (\mu \nabla' u) \, d\Omega' = \int_{\Omega'} \left(\frac{R}{r'}\right)^4 (\nabla' v) \cdot (\mu \mathbf{H}_s) \, d\Omega' - \int_{\Gamma'} v (n \cdot \mu \mathbf{H}'_s) \, d\Gamma'$$

ここで：
- $r' = \sqrt{x'^2 + y'^2 + z'^2}$
- $\mathbf{H}_s = (0, 0, 1)$（元の背景磁場、定数）
- $\mathbf{H}'_s = (0, 0, (r'/R)^2)$（境界での変換後背景磁場）
- $\Omega'$ = 変換後外部領域（実装では球殻$R < r' < 2R$）
- $\Gamma'$ = $r' = R$での境界（Kelvin界面）
- $R = 2.0$ m（Kelvin半径）

**2次元との主な違い**：3次元では双線形形式と線形形式の両方にKelvin因子がある！

### 5.4 周期境界条件による境界項の消去

**周期境界条件を用いたKelvin変換の主な利点：**

**周期境界条件**（$r = R$で内部領域と外部領域を接続）を用いたKelvin変換を使用する場合、境界積分項は**不要**になる：

$$\boxed{\int_{\Gamma'} v (n \cdot \mu \mathbf{H}'_s) \, d\Gamma' = 0}$$

**理由：**

1. **周期BCが連続性を強制**：周期境界条件はKelvin界面$r = R$でのポテンシャルとその法線微分の連続性を自動的に強制する。

2. **自動的なフラックスバランス**：境界の内側と外側からの寄与は周期的な同一視により相殺する。

3. **簡略化された実装**：境界積分を明示的に計算または適用する必要がない。

**周期BCを用いた簡略化された弱形式：**

$$\boxed{\int_{\Omega_{\text{interior}} \cup \Omega'_{\text{exterior}}} \left(\frac{R}{r}\right)^2 (\nabla v) \cdot (\mu \nabla u) \, d\Omega = \int_{\Omega_{\text{interior}}} (\nabla v) \cdot (\mu \mathbf{H}_s) \, d\Omega + \int_{\Omega'_{\text{exterior}}} \left(\frac{R}{r'}\right)^4 (\nabla v) \cdot (\mu \mathbf{H}_s) \, d\Omega}$$

ここでKelvin重み因子は外部領域にのみ適用される。

**周期BCを用いたNGSolve実装：**

```python
# 周期境界条件を作成
fes = H1(mesh, order=2, dirichlet="GND")
fes = Periodic(fes)  # 周期BCを適用 - 境界項が不要に！

# 外部領域に(R/r')^2因子を持つ双線形形式
a = BilinearForm(fes)
a += kelvin_weight_bilinear * mu * grad(u) * grad(v) * dx

# 外部領域に(R/r')^4因子を持つ線形形式（境界項不要！）
f = LinearForm(fes)
f += mu * InnerProduct(grad(v), Hs_weighted) * dx  # 体積積分のみ

# 境界積分不要: -mu*v*InnerProduct(n, Hsb)*ds(...) は必要ない！
```


## 6. NGSolve実装戦略

### 6.1 領域分割

1. **内部領域**$(r < R = 2.0$ m$)$：
   - 磁性球：$r < 0.5$ m、$\mu_r = 10$
   - 内部空気：$0.5$ m $< r < 2.0$ m、$\mu_r = 1$
   - 標準弱形式
   - 背景磁場：$\mathbf{H}_s = (0, 0, 1)$

2. **外部領域**$(r > R = 2.0$ m、球殻$R < r' < 2R$として表現)$：
   - Kelvin変換された空気：$R < r' < 2R = 4.0$ m、$\mu_r = 1$
   - $(R/r')^2$と$(R/r')^4$因子を持つ修正弱形式
   - 背景磁場：$\mathbf{H}_s = (0, 0, 1)$（積分内で$(R/r')^4$の重みを持つ）

3. **界面**$(r = R = 2.0$ m$)$：
   - 境界項を適用
   - $\phi$の連続性は自動的に強制される
   - 法線$\mathbf{B}$成分の連続性は弱形式により強制される

### 6.2 NGSolveにおける実装

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
    "air_outer": mu_r_air * mu0  # 変換後外部
}
mu = CoefficientFunction([mu_d[mat] for mat in mesh.GetMaterials()])
```

#### Kelvin因子を持つ背景磁場

```python
# 動径距離
r = sqrt(x**2 + y**2 + z**2)
kelvin_radius = 2.0  # R [m]

# Kelvin重み因子
# 双線形形式用: 外部領域で(R/r')^2
kelvin_weight_bilinear = IfPos(r - kelvin_radius + 1e-6,  # r > R の場合（外部領域）
                               (kelvin_radius / r)**2,     # (R/r')^2
                               1.0)                         # 内部では重みなし

# 線形形式用: 外部領域で(R/r')^4
kelvin_weight_linear = IfPos(r - kelvin_radius + 1e-6,    # r > R の場合（外部領域）
                             (kelvin_radius / r)**4,       # (R/r')^4
                             1.0)                           # 内部では重みなし

# 背景磁場（z方向定数）
Hs = CoefficientFunction((0, 0, 1))  # 元の背景磁場
Hs_weighted = CoefficientFunction((0, 0, kelvin_weight_linear))  # 線形形式用

# 境界条件用
Hsb = BoundaryFromVolumeCF(CoefficientFunction((0, 0, (r/kelvin_radius)**2)))
```

#### 弱形式

```python
fes = H1(mesh, order=2)
u = fes.TrialFunction()
v = fes.TestFunction()

# 外部領域に(R/r')^2因子を持つ双線形形式
a = BilinearForm(fes)
a += kelvin_weight_bilinear * mu * grad(u) * grad(v) * dx

# 外部領域に(R/r')^4因子を持つ線形形式
f = LinearForm(fes)
f += mu * InnerProduct(grad(v), Hs_weighted) * dx

# Kelvin界面(r = R)での境界項
f += -mu * v * InnerProduct(n, Hsb) * ds(mesh.Boundaries("kelvin_interface"))
```

#### 場の抽出

```python
# 摂動場（両領域で同じ式）
H_pert = -grad(gfu)

# 全磁場 = 摂動 + 背景
# 外部領域では定数背景Hs = (0, 0, 1)を使用
H_total = H_pert + Hs
```

## 7. 2次元と3次元の比較

| 性質 | 2次元（極座標） | 3次元（球座標） |
|------|-----------------|-----------------|
| 座標変換 | $r' = R^2/r$、$\theta' = \theta$ | $r' = R^2/r$、$\theta' = \theta$、$\phi' = \phi$ |
| ヤコビアン行列式 | $\|\det(J)\| = (r'/R)^2$ | $\|\det(J)\| = (r'/R)^2$ |
| 体積要素 | $d\Omega = (R/r')^4 d\Omega'$ | $d\Omega = (R/r')^6 d\Omega'$ |
| 勾配スケーリング | $\nabla = (r'/R)^2 \nabla'$ | $\nabla = (r'/R)^2 \nabla'$ |
| 場スケーリング | $\mathbf{H}' = (r'/R)^2 \mathbf{H}$ | $\mathbf{H}' = (r'/R)^2 \mathbf{H}$ |
| **双線形形式** | $(r'/R)^4 \cdot (R/r')^4 = 1$ | $(r'/R)^4 \cdot (R/r')^6 = (R/r')^2$ |
| **線形形式** | $(r'/R)^2 \cdot (R/r')^4 = (R/r')^2$ | $(r'/R)^2 \cdot (R/r')^6 = (R/r')^4$ |

**重要な洞察**：3次元での異なる体積要素スケーリング（$(R/r')^6$対$(R/r')^4$）は、2次元では線形形式のみがKelvin重み因子を必要とするのに対し、3次元では双線形形式と線形形式の両方がKelvin重み因子を必要とすることを意味する。

## 8. 主要結果のまとめ

### 8.1 座標変換（球座標）

$$\boxed{r' = \frac{R^2}{r}, \quad \theta' = \theta, \quad \phi' = \phi}$$

### 8.2 体積要素

$$\boxed{d\Omega = \left(\frac{R}{r'}\right)^6 d\Omega'}$$

### 8.3 勾配変換

$$\boxed{\nabla \phi = \left(\frac{r'}{R}\right)^2 \nabla' \phi'}$$

### 8.4 弱形式

**双線形形式**（外部領域に$(R/r')^2$因子付き）：

$$\boxed{a(u, v) = \int_{\Omega'} \left(\frac{R}{r'}\right)^2 (\nabla' v) \cdot (\mu \nabla' u) \, d\Omega'}$$

**線形形式**（外部領域に$(R/r')^4$因子付き）：

$$\boxed{f(v) = \int_{\Omega'} \left(\frac{R}{r'}\right)^4 (\nabla' v) \cdot (\mu \mathbf{H}_s) \, d\Omega'}$$

ここで$\mathbf{H}_s = (0, 0, 1)$は**定数**の元の背景磁場。

## 9. 参考文献

1. **Kelvin変換**：
   - Morse and Feshbach, "Methods of Theoretical Physics", Vol. 2
   - Jackson, "Classical Electrodynamics", 3rd ed., Section 1.13

2. **等角写像**：
   - Churchill and Brown, "Complex Variables and Applications"

3. **電磁気学における有限要素法**：
   - Jin, "The Finite Element Method in Electromagnetics", 3rd ed.
   - Monk, "Finite Element Methods for Maxwell's Equations"

## 10. ダイポール場の実装と検証

### 10.1 Kelvin変換なしの実装（有限領域）

**ファイル**：`3D_dipole.py`

**背景磁場**：$\mathbf{H}_s = (0, 0, 1)$（z方向）

**ジオメトリ**：
- 磁性球：半径$a = 0.5$ m、比透磁率$\mu_r = 10$
- 空気領域：$0.5$ m $< r < 4.0$ m

**弱形式**：
```python
# 境界項付き線形形式（Kelvin変換なしでは必須）
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx                         # 体積積分
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer")) # 境界項
```

**注意**：外部境界での境界項は正しい結果を得るために不可欠。

### 10.2 Kelvin変換の実装

**ファイル**：`3D_dipole_with_Kelvin.py`

**背景磁場の変換**：
- 内部：$\mathbf{H}_s = (0, 0, 1)$
- 外部：$\mathbf{H}'_s = (0, 0, -(R/r')^2)$（符号反転と空間変調）

**注意**：Kelvin + 周期BCでは境界項不要。

## 11. 4極場の実装

### 11.1 Kelvin変換なしの実装（有限領域）

**ファイル**：`3D_quadrupole.py`

**背景磁場**：$\mathbf{H}_s = (x, 0, -z)$（X-Z平面における4極場）

この場はポテンシャルから導かれる：
$$\phi_s = xz$$

**ジオメトリ**：
- 磁性球：半径$a = 0.5$ m、比透磁率$\mu_r = 100$
- 空気領域：$0.5$ m $< r < 4.0$ m

**弱形式**：
```python
# 境界項付き線形形式（Kelvin変換なしでは必須）
f = LinearForm(fes)
f += mu*InnerProduct(grad(v), Hs)*dx                         # 体積積分
f += -mu*v*InnerProduct(n, Hsb)*ds(mesh.Boundaries("outer")) # 境界項
```

**注意**：外部境界での境界項は正しい結果を得るために不可欠。

### 11.2 問題設定（Kelvin変換）

[3D_quadrupole_with_Kelvin.py](3D_quadrupole_with_Kelvin.py)の実装は**4極場**構成を使用：

**背景磁場**：
$$\mathbf{H}_s = (-z, 0, -x) \quad \text{（X-Z平面における4極場）}$$

この場は**回転フリー**（カールフリー）であり、ポテンシャルから導かれる：
$$\phi_s = xz$$

**ジオメトリ**：
- 磁性球：半径$a = 0.5$ m、比透磁率$\mu_r = 100$
- 内部空気領域：$0.5$ m $< r < 2.0$ m
- Kelvin変換：$R = 2.0$ m
- 外部領域（変換後）：$0 < r' < 2.0$ m

**背景磁場のKelvin変換**：

外部領域では、背景磁場は以下のように変換される：
$$\mathbf{H}'_s = \left(\frac{r'}{R}\right)^2 (z', 0, x')$$

ここで$(x', y', z')$は外部領域中心を原点とする局所座標。

### 11.3 解析解

4極場中の磁性球について、摂動ポテンシャルは：

**内部**（$r < a$）：
$$\phi_{\text{pert}} = B \cdot xz$$

ここで係数$B$は境界条件から決定される。

**外部**（$r > a$）：
$$\phi_{\text{pert}} = A \cdot \frac{xz}{r^5}$$

ここで$A$は$A = B \cdot a^5$により$B$と関連する。

**境界条件**（$r = a$）：
1. ポテンシャル連続性：$\phi_{\text{in}} = \phi_{\text{out}}$
2. 法線磁束密度連続性：$\mu_r H_{r,\text{in}} = H_{r,\text{out}}$

4極場ポテンシャル$\phi = xz = r^2 \sin^2\theta \cos\theta \cos\phi$について、動径微分は$2r$因子を含む：
$$\frac{\partial \phi}{\partial r} = 2xz/r$$

境界条件を適用すると：

$$\boxed{B = -\frac{2(\mu_r - 1)}{2\mu_r + 5}}$$

$$\boxed{A = -\frac{2(\mu_r - 1)}{2\mu_r + 5} \cdot a^5}$$

$\mu_r = 100$、$a = 0.5$ mの場合：
- $B = -0.96618$
- $A = -0.03019$

**摂動場**：

内部（$r < a$）：
$$\mathbf{H}_{\text{pert}} = -\nabla\phi_{\text{pert}} = B(-z, 0, -x)$$

外部（$r > a$）：
$$\mathbf{H}_{\text{pert}} = -\nabla\left(\frac{A \cdot xz}{r^5}\right) = -A\left(\frac{z(r^2 - 5x^2)}{r^7}, 0, \frac{x(r^2 - 5z^2)}{r^7}\right)$$

## 12. まとめ

### 12.1 2次元との主要な違い

| 性質 | 2次元 | 3次元 |
|------|-------|-------|
| 体積要素スケーリング | $(R/r')^4$ | $(R/r')^6$ |
| 双線形形式因子 | 1（相殺） | $(R/r')^2$ |
| 線形形式因子 | $(R/r')^2$ | $(R/r')^4$ |
| 場変換 | H' = -H（符号のみ） | H' = -(R/r')² H（変調） |

### 12.2 境界項のまとめ

| 手法 | 外部での境界項 | 理由 |
|------|----------------|------|
| **Kelvin変換なし**（有限領域） | **必要** | 外部境界を通るフラックスを考慮する必要 |
| **Kelvin + 周期BC** | **不要** | 周期BCが自動的に連続性を強制 |
