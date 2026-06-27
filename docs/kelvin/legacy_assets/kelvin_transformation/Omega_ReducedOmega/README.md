# Ω-Reduced Ω 法（スカラー磁位法）

## 概要

Ω-Reduced Ω 法は、静磁場解析において全スカラー磁位（Total Scalar Potential）と還元スカラー磁位（Reduced Scalar Potential）を組み合わせた手法である。

## 符号規約

docs に昇格した classic 実装群（`docs/kelvin/kelvin_classic_demos.ipynb`
および `docs/kelvin/kelvin_classic_demos_results.json` に保存）では、
以下の符号規約を採用：

$$
\mathbf{H} = \nabla \Omega
$$

$$
\mathbf{B} = \mu \nabla \Omega
$$

**注意**: 一般的な教科書では $\mathbf{H} = -\nabla \Omega$ と定義されることが多いが、本実装では上記の符号規約を採用している。

## ソースポテンシャル

外部磁場 $\mathbf{H}_s = (0, 0, H_0)$ に対するソースポテンシャル：

$$
\Omega_s = H_0 \cdot z
$$

これにより $\nabla \Omega_s = (0, 0, H_0) = \mathbf{H}_s$ が成り立つ。

## 二領域法

### 領域定義

- **Total領域（全磁位領域）**: 磁性体内部。ソース項なしで解く。
- **Reduced領域（還元磁位領域）**: 空気領域。ソース項あり。

**重要**: Reduced領域には磁性体を置けない（μ = μ0 のみ許容）。

### 磁場の表現

**Total領域（磁性体）:**
$$
\mathbf{H}_t = \nabla \Omega_t
$$

**Reduced領域（空気）:**
$$
\mathbf{H}_r = \mathbf{H}_s + \nabla \Omega_r = \nabla \Omega_s + \nabla \Omega_r
$$

### 摂動場の計算

磁場の摂動成分（外部磁場からのずれ）を計算する場合：

**Total領域での摂動場:**
$$
\mathbf{B}_{pert,t} = \mu \nabla \Omega_t - \mathbf{B}_s = \mu \nabla \Omega_t - \mu_0 \mathbf{H}_s
$$

**Reduced領域での摂動場:**
$$
\mathbf{B}_{pert,r} = \mu_0 \nabla \Omega_r
$$

磁気エネルギーの計算には摂動場のみを使用（全場を使うとエネルギーが発散）。

## 弱形式

### 双一次形式

$$
a(\Omega, \psi) = \int_{\Omega_{total}} \mu (\nabla \Omega) \cdot (\nabla \psi) \, d\Omega + \int_{\Omega_{reduced}} \mu_0 (\nabla \Omega) \cdot (\nabla \psi) \, d\Omega
$$

### 線形形式

$$
f(\psi) = \int_{\Omega_{reduced}} \mu_0 (\nabla \Omega_s) \cdot (\nabla \psi) \, d\Omega + \int_{\Gamma} (\mathbf{n} \cdot \mathbf{B}_s) \psi \, ds
$$

### 境界条件

Total/Reduced領域の界面では：
- Dirichlet条件: $\Omega = \Omega_s$ を境界上に設定
- Neumann条件: $\mathbf{n} \cdot \mathbf{B}_s$ を界面項として追加

### 3D法線方向の注意

**重要**: NGSolveの3Dにおいて、`specialcf.normal(mesh.dim)` は球面境界上で**内向き**（空気から磁性体へ）を指す。Omega-Reduced Omega法では**外向き**法線が必要なため、3Dでは法線を反転する必要がある：

```python
# 3D: 法線を反転
normal = -specialcf.normal(mesh.dim)
```

2D軸対称では法線は既に外向きなので反転不要：

```python
# 2D axisymmetric: 反転不要
normal = specialcf.normal(mesh.dim)
```

## NGSolve実装

### ソースポテンシャルとソース場

```python
# コイル(外部ソース)からのスカラーポテンシャルと場
Ov = Ofield(coil)      # Ωs: ソーススカラーポテンシャル
Bv = Bfield(coil)      # Bs: ソース磁束密度
Hv = Bv / mu0          # Hs: ソース磁場

# 領域ごとのソース場（Reduced領域でのみ非ゼロ）
Bs_dic = {"iron": zero, "Omega_domain": Bv, ...}
Bs = CoefficientFunction([Bs_dic[mat] for mat in mesh.GetMaterials()])
```

**注意**: 一様磁場 $H_0$ in z方向の場合:
- $\Omega_s = H_0 \cdot z$
- $\mathbf{H}_s = (0, 0, H_0)$
- $\mathbf{B}_s = \mu_0 \mathbf{H}_s = (0, 0, \mu_0 H_0)$

### 双一次形式

```python
a = BilinearForm(fes)
a += Mu * (grad(omega) * grad(psi)) * dx(total_region)
a += Mu * (grad(omega) * grad(psi)) * dx(reduced_region)
```

### 線形形式

```python
# 1. 界面にソースポテンシャルを設定（Dirichlet条件のような形で）
gfOmega.Set(Ov, BND, mesh.Boundaries(total_boundary))

# 2. Reduced領域でのソース項（界面付近で寄与）
f = LinearForm(fes)
f += Mu * grad(gfOmega) * grad(psi) * dx(reduced_region)
f.Assemble()

# 3. FreeDofs処理（Dirichlet DOFを除去）
fcut = np.array(f.vec.FV())[fes.FreeDofs()]
np.array(f.vec.FV(), copy=False)[fes.FreeDofs()] = fcut

# 4. Neumann境界条件（界面でのBn連続性）
f += (normal * Bv) * psi * ds(total_boundary)
f.Assemble()
```

**重要**: `gfOmega.Set(Ov, BND, ...)` により、gfOmegaは界面境界上でのみΩsの値を持ち、内部では0である。したがって `grad(gfOmega)` は界面付近の要素でのみ非ゼロとなり、ソース項は界面を通じてのみ作用する。

### 磁束密度の計算

```python
# Total領域の磁束密度
Bt = grad(Ot) * Mu

# Reduced領域の摂動磁束密度
Br = (grad(Orr) - grad(Oxr)) * mu0

# ソース磁束密度（Reduced領域のみ）
Bs = mu0 * Hs

# 全磁束密度
BField = Bt + Br + Bs
```

**注意**: Total領域（磁性体）での摂動場は `Bt - Bs` で計算する必要がある。

### 磁性体内部の磁場

磁性体球の内部では、透磁率効果により磁場は一様で、外部磁場と同じ向きだが減少する：

$$
H_{in} = \frac{3}{2 + \mu_r} H_0
$$

μr = 100 の場合：$H_{in} \approx 0.0294 H_0$

磁束密度は逆に増加する：$B_{in} = \mu_r \mu_0 H_{in} = \frac{3\mu_r}{2 + \mu_r} \mu_0 H_0 \approx 2.94 \mu_0 H_0$

## 摂動場エネルギー

全場を使うと無限遠までの積分でエネルギーが発散するため、摂動場を使用してエネルギーを計算する。

### 摂動場の定義

**Total領域（磁性体）:**
$$
\mathbf{H}_{pert,t} = \nabla \Omega_t - \mathbf{H}_s = \nabla \Omega_t - (0, 0, H_0)
$$
$$
\mathbf{B}_{pert,t} = \mu_r \mu_0 \mathbf{H}_{pert,t}
$$

**Reduced領域（空気）:**
$$
\mathbf{H}_{pert,r} = \nabla \Omega_r
$$
$$
\mathbf{B}_{pert,r} = \mu_0 \mathbf{H}_{pert,r}
$$

### 摂動場エネルギーの計算

磁気エネルギーは以下で計算される：

$$
W_{pert} = \frac{1}{2} \int_V \mathbf{B}_{pert} \cdot \mathbf{H}_{pert} \, dV
$$

各領域での計算：

**Total領域:**
$$
W_t = \frac{1}{2} \int_{V_{magnetic}} \mu_r \mu_0 |\mathbf{H}_{pert,t}|^2 \, dV
$$

**Reduced領域:**
$$
W_r = \frac{1}{2} \int_{V_{air}} \mu_0 |\mathbf{H}_{pert,r}|^2 \, dV
$$

**Kelvin領域（変換後）:**
$$
W_k = \frac{1}{2} \int_{V_{kelvin}} \mu'(r') |\mathbf{H}'_{pert}|^2 \, dV'
$$

ここで $\mu'(r') = (R/r')^2 \mu_0$ はKelvin変換後の透磁率。

### 解析解との比較

#### 球内部エネルギー

球内部の摂動磁場は一様で、外部磁場と逆向き：

$$
H_{pert,in} = -\frac{\mu_r - 1}{\mu_r + 2} H_0
$$

内部エネルギー（$V_{sphere} = \frac{4\pi a^3}{3}$）：

$$
W_{in} = \frac{1}{2} \mu_r \mu_0 |H_{pert,in}|^2 V_{sphere} = \frac{2\pi}{3} \mu_r \mu_0 \left(\frac{\mu_r - 1}{\mu_r + 2}\right)^2 H_0^2 a^3
$$

#### 球外部エネルギー（双極子場）

磁性体球は磁気双極子として振る舞う。双極子モーメント：

$$
m = 4\pi a^3 \frac{\mu_r - 1}{\mu_r + 2} H_0
$$

双極子場のエネルギーは解析的に計算でき：

$$
W_{out} = \frac{\mu_0 m^2}{12 \pi a^3} = \frac{4\pi}{3} \mu_0 \left(\frac{\mu_r - 1}{\mu_r + 2}\right)^2 H_0^2 a^3
$$

#### 内部と外部のエネルギー比

$$
\frac{W_{in}}{W_{out}} = \frac{\mu_r}{2}
$$

高透磁率材料（$\mu_r \gg 1$）では内部エネルギーが支配的。

#### 総エネルギー

$$
W_{total} = W_{in} + W_{out} = \frac{2\pi}{3} \mu_0 \left(\frac{\mu_r - 1}{\mu_r + 2}\right)^2 H_0^2 a^3 \left(\mu_r + 2\right)
$$

#### 数値例

μr = 100、a = 0.5 m、H0 = 1 A/m の場合：
- $W_{in} = 3.099 \times 10^{-5}$ J
- $W_{out} = 6.198 \times 10^{-7}$ J
- $W_{in}/W_{out} = 50$（= μr/2）
- $W_{total} = 3.161 \times 10^{-5}$ J

### NGSolve実装

```python
# 摂動場の定義
H_pert_total = grad(gfu) - Hs  # Total領域
H_pert_reduced = grad(Orr) - grad(Oxr)  # Reduced領域

# エネルギー計算
energy_total = Integrate(0.5 * mu_r * mu0 * InnerProduct(H_pert_total, H_pert_total) * dx("magnetic"), mesh)
energy_reduced = Integrate(0.5 * mu0 * InnerProduct(H_pert_reduced, H_pert_reduced) * dx("air_inner"), mesh)
energy_kelvin = Integrate(0.5 * mu_kelvin * InnerProduct(H_pert_kelvin, H_pert_kelvin) * dx("air_outer"), mesh)
```

## Kelvin変換

無限遠境界条件を扱うため、Kelvin変換を適用：

```python
rs = model.rKelvin
xs = 2 * rs
r = sqrt((x - xs)**2 + y**2 + z**2)
fac = rs**2 / r**2

a += Mu * fac * (grad(omega) * grad(psi)) * dx("Kelvin")
```

変換後の透磁率: $\mu'(r') = (R/r')^2 \mu_0$

## 参考文献

- O. Bíró, K. Preis, "On the use of the magnetic vector potential in the finite element analysis of three-dimensional eddy currents," IEEE Trans. Magn., vol. 25, no. 4, pp. 3145-3159, 1989.
- J. Simkin, C.W. Trowbridge, "On the use of the total scalar potential on the numerical solution of fields problems in electromagnetics," Int. J. Numer. Methods Eng., vol. 14, pp. 423-440, 1979.
