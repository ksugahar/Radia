# 平衡誤差推定量 (Equilibrated Error Estimator)

有限要素法（FEM）における事後誤差評価（A Posteriori Error Estimation）と、計算電磁気学に特化した Edge Element（辺要素）に関する理論と実装をまとめます。
NGSolve を用いた電磁場解析、特に Kelvin 変換を用いた特異点を含む解析において、信頼性の高い適応型メッシュ制御（AMR）を実現するための理論的基盤となります。

---

## 目次

1. [事後誤差評価の基礎技術](#1-事後誤差評価の基礎技術)
2. [エッジ要素のための誤差評価（Braess-Schöberl）](#2-エッジ要素のための誤差評価)
3. [幾何学的定式化（Bossavit）](#3-幾何学的定式化)
4. [思想の合流点](#4-思想の合流点)
5. [NGSolve における実装](#5-ngsolve-における実装)
6. [境界条件の対称性](#6-境界条件の対称性)
7. [Kelvin変換の効率性](#7-kelvin変換の効率性)
8. [Ω法による平衡誤差推定](#8-ω法による平衡誤差推定)
   - [8.6 CG法による高速近似（CG-Smoother）](#86-cg法による高速近似cg-smoother)
   - [8.7 平衡誤差推定の双対性と実用上の制約](#87-平衡誤差推定の双対性と実用上の制約)

---

## 1. 事後誤差評価の基礎技術

**参考文献:** `A_posteriori_Error_Estimation_Techniques_in_Practical_Finite_Element_Analysis.pdf`

有限要素解析の信頼性を保証し、効率的な計算を行うための核心技術である「事後誤差評価」の体系的な解説です。

### 主なトピックとメカニズム
*   **AMR の指針**: 数値解そのものを用いて実際の誤差を推定し、メッシュ細分化（h-version）や次数上げ（p-version）を自動で行います。
*   **基本的なメカニズム（残差ベース）**:
    *   **Interior Residual**: 方程式 $Lu=f$ に対する要素内部での残差 $f - Lu_h$。
    *   **Jump Residual**: 要素境界でのフラックスの不連続性（ジャンプ）。
    これらを合計することで、最も直感的かつ低コストに誤差を推定します。
*   **評価指標**:
    *   **Reliability (信頼性)**: 誤差の上界を正しく評価（安全側）。
    *   **Efficiency (効率性)**: 推定が過大でないこと（無駄がない）。

---

## 2. エッジ要素のための誤差評価

**参考文献:** `Equilibrated Residual Error Estimator for Edge Elements.pdf`
**著者:** Dietrich Braess, Joachim Schöberl
**出版:** Mathematics of Computation, Vol. 77, No. 262, pp. 651–672, **April 2008**

Maxwell 方程式（H(curl) 問題）の **辺要素（Edge Elements）に特化** した、高度な誤差評価手法「Braess-Schöberl Estimator」についての詳細です。

> **歴史的背景**: 平衡誤差推定の基礎となる Prager-Synge の定理は **1947年** に発表されました。Braess-Schöberl (2008) は、この古典的定理を約60年後にエッジ要素（Maxwell方程式）へ拡張し、**汎用定数なしの厳密な誤差評価**を初めて数学的に証明しました。

### 2.1. 平衡化によるメカニズム

標準的な辺要素（A法）の数値解は、磁束密度 $B = \nabla \times A$ については厳密な適合性 (**B-strong**) を持ちますが、磁場 $H$ の保存則（Ampereの法則など）については弱形式の意味でしか満たされておらず、不連続性が残ります (**H-weak**)。

この手法の核心は、事後処理によって **「保存則を厳密に満たす（H-strong / Equilibrated）新しい場 $\sigma_{eq}$」** を構成し、元の数値解と比較することにあります。

*   **B-strong vs H-strong**:
    *   **FEM解 ($u_h$)**: B-strong, H-weak
    *   **Equilibrated Flux ($\sigma_{eq}$)**: H-strong (Conservation Exact)
    *   この「Strongな性質の違い」を利用して、真の解との距離（誤差）を挟み撃ちにするのがメカニズムの本質です。

*   **Prager-Synge の定理**: 保存則を満たす場 $\sigma_{eq}$ と数値解の勾配 $\nabla u_h$ の差は、常に**厳密な誤差の上界（Reliable Upper Bound）** となります。
*   **p-robustness**: 多項式次数 $p$ が変化しても、この上界の信頼性は揺らぎません。これは hp-FEM において極めて重要です。

### 2.2. 平衡誤差推定量が「最良」である数学的根拠

#### (1) Prager-Syngeの定理：等式評価の威力

**スカラー方程式（Poisson方程式 $-\Delta u = f$）の場合** (Theorem 1, p.653):

$\sigma \in H(\text{div})$ が平衡条件 $\text{div}\,\sigma + f = 0$ を満たすとき：

$$\|\nabla u - \nabla v\|^2 + \|\nabla u - \sigma\|^2 = \|\nabla v - \sigma\|^2$$

**Maxwell方程式（curl-curl方程式）の場合** (Theorem 10, p.665):

$\tilde{H} \in H(\text{curl})$ がAmpèreの法則 $\text{curl}\,\tilde{H} = j$ を満たすとき：

$$\|\mu^{-1/2} \text{curl}(u-v)\|^2 + \|\mu^{1/2}(H - \tilde{H})\|^2 = \|\mu^{-1/2}(\text{curl}\,v - \mu\tilde{H})\|^2$$

#### (2) なぜ「最良」なのか？— 直交性による等式評価

証明の核心は**直交性**にあります（p.665）：

$$\int_\Omega \text{curl}(u-v) \cdot (H - \tilde{H}) = \int_\Omega (u-v) \cdot \text{curl}(H - \tilde{H}) + \text{境界項} = \int_\Omega (u-v) \cdot (j - j) = 0$$

この直交性により、ピタゴラスの定理の形で誤差が**等式**として表現されます。

#### (3) 通常の誤差推定との決定的な違い

| 推定手法 | 評価形式 | 定数依存性 | 信頼性 |
|:---|:---|:---|:---|
| **残差ベース推定** | 不等式 $c_1 \eta \leq \Vert\text{error}\Vert \leq c_2 \eta$ | 汎用定数 $c_1, c_2$ に依存 | 定数の見積もりが困難 |
| **平衡誤差推定** | **等式** $\Vert\text{error}\Vert^2 + \cdots = \Vert\text{estimator}\Vert^2$ | **主要項に定数なし** | 数学的に厳密 |

#### (4) 最良性の具体的帰結

**Theorem 4 (p.656)** および **Theorem 13 (p.670)** より、平衡誤差推定量は以下を満たします：

$$c_0\|\sigma^\Delta\| - ch\|f - \bar{f}\| \leq \|\nabla(u-u_h)\| \leq \|\sigma^\Delta\| + ch\|f - \bar{f}\|$$

*   **右辺の主要項 $\|\sigma^\Delta\|$**: 汎用定数なしで誤差の上界
*   **データ振動項 $ch\|f - \bar{f}\|$**: 高次項のみに定数が現れる（無視できる）

#### (5) 技術的基盤：de Rham系列の完全性

平衡誤差推定量の構成には **de Rham系列の完全性（Exact Sequence）** が不可欠です（Section 3）：

$$\mathbb{R} \to H^1 \xrightarrow{\text{curl}} H(\text{div}) \xrightarrow{\text{div}} L_2 \to 0 \quad \text{（2D）}$$

$$0 \to H^1_0 \xrightarrow{\text{grad}} H_0(\text{curl}) \xrightarrow{\text{curl}} H_0(\text{div}) \xrightarrow{\text{div}} L_2 \to \mathbb{R} \to 0 \quad \text{（3D）}$$

この完全性により、発散がゼロの分布から局所的に平衡化されたフラックス $\sigma$ を構成できます（Lemma 3, p.654-655）。

#### (6) 結論：平衡誤差推定量が「最良」である理由のまとめ

1. **Prager-Syngeの定理による等式評価**: 不等式ではなく等式で誤差を捉える
2. **直交性の活用**: 真の誤差と推定量の間にピタゴラス的関係が成立
3. **汎用定数の排除**: 主要項に generic constant が入らない
4. **信頼性と効率性の同時達成**: 上界・下界ともに tight な評価
5. **de Rham複体との整合性**: 数学的構造に基づく厳密な構成

### 2.3. Cheap Construction（低コスト構成）

*   **局所パッチ（Local Patch）**: 全領域で連立方程式を解くのではなく、各頂点ごとの小さな領域（パッチ）で独立した問題を解き、それを貼り合わせる（Partition of Unity）ことで $\sigma_{eq}$ を構成します。
*   この局所化により、計算コストは主問題に比べて無視できるほど小さくなります（**Cheap**）。

### 2.4. 論文の適用範囲と PatchSolve による数値的拡張

#### (1) Braess-Schöberl 論文の適用範囲（限界）

Braess-Schöberl (2008) の論文で**厳密に証明されているのは最低次要素の場合のみ**です：

| 要素タイプ | 論文での扱い | 備考 |
|:---|:---|:---|
| **最低次 Nédélec 要素** | ✅ 厳密に証明 | Theorem 4, 13 で完全な証明 |
| **高次 Nédélec 要素 ($p > 1$)** | ❌ 直接は扱っていない | 拡張は自明ではない |
| **T-Ω法（2ポテンシャル）** | ❌ 対象外 | H-strong/B-weak の双対問題 |

論文の構成（式 (4.11), (4.12), (4.17), (4.18)）は、最低次要素の自由度構造に特化した**閉じた代数的公式**として導出されています。

#### (2) PatchSolve による数値的構成：高次・2ポテンシャルへの拡張

しかし、Prager-Synge の定理自体は**任意の次数・任意の定式化**で成立します。問題は「平衡条件を満たす場 $\sigma_{eq}$ をどう構成するか」です。

**PatchSolve アプローチ**は、閉じた代数公式ではなく、**局所的な最小化問題を数値的に解く**ことで $\sigma_{eq}$ を構成します：

```
各頂点パッチ ωV において：
    minimize  ‖σωV‖²
    subject to  div σωV = (残差の局所成分)
                σωV · n = 0  on ∂ωV
```

| 特性 | 閉じた公式（論文） | PatchSolve（数値的） |
|:---|:---|:---|
| **最低次要素** | ✅ 厳密解 | ✅ 同等の結果 |
| **高次要素 ($p > 1$)** | ❌ 公式なし | ✅ 数値的に構成可能 |
| **T-Ω法（2ポテンシャル）** | ❌ 対象外 | ✅ 適用可能 |
| **計算コスト** | O(1) per patch | O(n³) per patch（小規模） |

### 2.5. 補間の次数と効率性

$\sigma_{eq}$ を構成する際、どの多項式次数空間を使うべきか？

*   **A法 (B-strong / H-weak Formulation)**:
    *   **$HCurl(p-1)$ への補間が最適**（Efficiency Proof）。
    *   理由: A法では、磁束密度 $B = \nabla \times A$ が "Strong"（厳密な適合関係）として扱われます。この微分操作（curl）により、物理場の次数は自然に1つ下がります（$p \to p-1$）。

    | 補間次数の候補 | 計算コスト | 数学的十分性 | 判定 |
    | :--- | :--- | :--- | :--- |
    | **Order $p+1$** | **High** | **Excessive** | 無駄 |
    | **Order $p$** | Medium | Redundant | 無駄が多い |
    | **Order $p-1$** | **Lowest** | **Sufficient** | **最適** |

*   **T-Omega法 (H-strong / B-weak Formulation)**:
    *   **ターゲット空間**: **$H(div)(p)$** への補間が最適。

### 2.6. 空間の選択: なぜ H(curl) なのか？

| 空間 | 物理的連続性 | 判定 | 理由 |
| :--- | :--- | :--- | :--- |
| **$H(curl)$** | **Tangential Only** | **最適** | 磁場 $H$ は「接線成分のみ連続」であるべき物理量 |
| **$Vector H^1$** | Full | **不適** | 連続性が強すぎる |
| **$Vector L^2$** | None | **Not Reliable** | 連続性がない |
| **$H(div)$** | Normal Only | **Mismatch** | 演算子が適合しない |

---

## 3. 幾何学的定式化

**参考文献:** `Computational Electromagnetism/` (Bossavit 関連)

Alain Bossavit 氏による、電磁場を「ベクトル場」ではなく「微分形式」として捉える幾何学的アプローチです。
電磁場の物理量（電位、電場、磁束など）を、それぞれの性質に合わせて 0-form, 1-form, 2-form に割り当てることで、Maxwell 方程式の構造（div, curl の関係）を離散レベルでも保ちます。

---

## 4. 思想の合流点

Alain Bossavit（幾何学的定式化）と Joachim Schöberl（誤差評価とNGSolve）の仕事は、**「数学的構造（de Rham 複体）への深い洞察」** において本質的に合流しています。

### 4.1. de Rham 複体の活用
*   **Bossavit**: 離散化の **「出発点」** として、物理法則を $H^1 \xrightarrow{grad} H(curl) \xrightarrow{curl} H(div) \xrightarrow{div} L^2$ という完全列上に配置することを提唱。
*   **Joachim**: 誤差評価の **「道具」** としてこの系列を利用。

### 4.2. 双対性と保存則
*   **Bossavit**: 物理量は本来、Primal Grid だけでなく Dual Grid 上にも存在し、保存則は **Dual Grid 上で厳密に成り立つ**。
*   **Joachim**: 数値解（Primal）では破れがちなこの保存則を、事後処理（Equilibration）によって **局所的に回復**。

### 4.3. 物理学者と数学者の関係

| 観点 | Bossavit | Braess-Schöberl (2008) |
|:---|:---|:---|
| アプローチ | 物理的直観・微分形式 | 関数解析・変分法 |
| 主張の性質 | **設計原理** | **定理と証明** |
| Prager-Synge | 言及するが厳密証明なし | **Theorem 1, 10 として厳密に証明** |

**結論**: Bossavit の幾何学的直観が「正しい方向」を示し、Braess-Schöberl が約60年越しに Prager-Synge の定理を拡張して**厳密な数学的証明**を与えました。

---

## 5. NGSolve における実装

NGSolve は、上記の理論を実践するために設計された、**hp-FEM（次数混成）** に特化したソルバーです。

### 5.1. Variable Order

異なる次数の要素を混在させても、エッジやフェースごとの自由度管理により、数学的整合性が保たれます。

```python
from ngsolve import *
mesh = Mesh(unit_square.GenerateMesh(maxh=0.5))

# Default Order p=3
fes = HCurl(mesh, order=3)

# Per-Element Order Control
for el in mesh.Elements():
    if el.nr < 10:
        fes.SetOrder(NodeId(ELEMENT, el.nr), 1)

fes.Update()
print ("Total DoFs:", fes.ndof)
```

### 5.2. PatchwiseSolve

**ソースコード情報**（`ngsolve/comp/localsolve.cpp`）：
```cpp
/* File:   localsolve.cpp
 * Author: Joachim Schoeberl
 * Date:   12. May. 2020
 *
 * local solve - useful for equilibration
 */
void PatchwiseSolve(shared_ptr<SumOfIntegrals> bf,
                    shared_ptr<SumOfIntegrals> lf,
                    shared_ptr<GridFunction> gf,
                    LocalHeap & lh);
```

**Python API**：
```python
from ngsolve import PatchwiseSolve

lhs_integrals = grad(sigma) * grad(tau) * dx
rhs_integrals = residual * tau * dx
gf_sigma.vec.data = PatchwiseSolve(EA, fes, lhs_integrals, rhs_integrals)
```

---

## 6. 境界条件の対称性

Omega法（スカラーポテンシャル）とA法（ベクトルポテンシャル）では、**Dirichlet境界条件とNeumann境界条件が逆になります**。

### 6.1. 境界条件の対応関係

| 境界条件 | Omega法での意味 | A法での意味 |
|:---|:---|:---|
| **Dirichlet BC** | $\Omega = \text{const}$ → $\mathbf{n} \times \mathbf{H} = 0$ | $\mathbf{n} \times \mathbf{A} = 0$ → $\mathbf{n} \times \mathbf{B} = 0$ |
| **Neumann BC** | $\mathbf{n} \cdot \mathbf{B} = 0$ | $\mathbf{n} \times \mathbf{H} = 0$ |

### 6.2. 周期境界条件とNGSolveのバグ修正

**修正方法**（`libsrc/meshing/meshfunc.cpp`）:
```cpp
// 修正後（周期BC情報を保持）
// mesh.GetIdentifications().GetIdentifiedPoints().DeleteData();
```

**動作確認**:
```python
fes_before = H1(mesh, order=2, dirichlet="GND")
fes_after = Periodic(fes_before)
print(f"FreeDofs: {sum(fes_before.FreeDofs())} -> {sum(fes_after.FreeDofs())}")
# 正常: FreeDofs が減少している
```

---

## 7. Kelvin変換の効率性

### 7.1. 無限遠領域のメッシュ数

*   **事実**: Kelvin変換された「外部領域」は、変換後の計算空間においては特異点を含む球の内部にコンパクトにマッピングされます。
*   **3Dの幾何学的特性**: この外部領域を表現するために必要なメッシュ数は、内部領域に比べて **驚くほど少数** で済みます。

### 7.2. AMR適用時の振る舞い

*   **AMRの集中**: 誤差評価に基づくAMRを行うと、メッシュ細分化は物理的な変化の激しい「内部領域」や「エッジ付近」に集中します。
*   **外部領域の不変性**: Kelvin変換された外部領域では解が滑らかに減衰するため、AMRを繰り返しても **この領域の要素数はほとんど増えません。**

---

## 8. Ω法による平衡誤差推定

静磁場問題において、A法（ベクトルポテンシャル法）の解に対する平衡誤差推定量は、Ω法（スカラーポテンシャル法）を解くことで得られる。

### 8.1. 理論的背景

**A法**:
$$\text{curl}\left(\frac{1}{\mu} \text{curl}(\mathbf{A})\right) = 0, \quad \mathbf{H}_A = \frac{1}{\mu}\text{curl}(\mathbf{A})$$

**Ω法**:
$$\text{div}(\mu \nabla \Omega) = 0, \quad \mathbf{H}_\Omega = \nabla \Omega$$

### 8.2. Prager-Synge の定理

真の解を $\mathbf{H}$、A法の解を $\mathbf{H}_A$、curl-free な場を $\mathbf{H}_{eq}$ とすると：

$$\|\mathbf{H} - \mathbf{H}_A\|_\mu^2 + \|\mathbf{H} - \mathbf{H}_{eq}\|_\mu^2 = \|\mathbf{H}_A - \mathbf{H}_{eq}\|_\mu^2$$

この等式から直ちに上界が得られる：
$$\|\mathbf{H} - \mathbf{H}_A\|_\mu \leq \|\mathbf{H}_A - \mathbf{H}_{eq}\|_\mu$$

**Ω法の解が使える理由**:
- $\mathbf{H}_\Omega = \nabla \Omega$ は自動的に curl-free（$\text{curl}(\nabla \Omega) = 0$）
- 同じ境界条件を課せば、A法とΩ法は同じ真の解 $\mathbf{H}$ への近似
- したがって、平衡誤差推定量は：

$$\eta = \|\mathbf{H}_A - \mathbf{H}_\Omega\|_\mu$$

### 8.3. 計算コストの考察

AMR（適応型メッシュ細分化）のコンテキストでは：

| 問題 | 空間 | DOF | 備考 |
|------|------|-----|------|
| A法 | HCurl（ベクトル） | 大 | **主コスト** |
| Ω法 | H1（スカラー） | 小 | A法より軽い |

**結論**: Ω法を追加で解いても、A法1回分より明らかに軽い。ソルバーの選択（直接法 vs 反復法）は実装の詳細であり、本質ではない。

### 8.4. 実装例

```python
from ngsolve import *

# A法を解いた後、H_A = (1/mu) * curl(gf_A) が得られている

# Ω法を解く
fes_O = H1(mesh, order=order, dirichlet="outer")
Omega, psi = fes_O.TnT()

a = BilinearForm(fes_O)
a += mu_cf * InnerProduct(grad(Omega), grad(psi)) * dx
a.Assemble()

gf_O = GridFunction(fes_O)
gf_O.Set(H0 * z, BND)  # 外部磁場に対応する境界条件

res = gf_O.vec.CreateVector()
res.data = -a.mat * gf_O.vec

# 直接法でも反復法でも良い
inv = a.mat.Inverse(fes_O.FreeDofs(), inverse="sparsecholesky")
gf_O.vec.data += inv * res

# 平衡誤差推定量
H_Omega = grad(gf_O)
eta = sqrt(Integrate(mu_cf * InnerProduct(H_A - H_Omega, H_A - H_Omega) * dx, mesh))
```

### 8.5. 多項式次数の選択

A法で次数 $p$ の HCurl 空間を使用した場合、Ω法には同じ次数 $p$ の H1 空間を使用する。

```python
fes_A = HCurl(mesh, order=order, ...)  # A法
fes_O = H1(mesh, order=order, ...)     # Ω法（同じ次数）
```

### 8.6. CG法による高速近似（CG-Smoother）

Ω法を直接法で完全に解く代わりに、**CG法を少数回で打ち切る**ことで、計算コストを大幅に削減できる。

#### (1) 理論的背景

平衡誤差推定量 $\mathbf{H}_{eq} = \nabla \Omega$ を求める問題は：

$$\min_\Omega \|\nabla \Omega - \mathbf{H}_A\|^2$$

変分原理から、以下の H1 問題を解くことに帰着される：

$$(\nabla \Omega, \nabla \psi) = (\mathbf{H}_A, \nabla \psi) \quad \forall \psi \in H^1$$

この問題を CG 法で**少数回反復**して近似解を得る。

#### (1.5) CG打ち切りが使える数学的理由

**なぜCG打ち切りでも平衡誤差推定が有効なのか？**

**鍵となる性質**: $\mathbf{H}_{eq} = \nabla \Omega$ は、**$\Omega$ がどんな値であっても常に curl-free** である。

$$\text{curl}(\nabla \Omega) = 0 \quad \text{（ベクトル解析の恒等式）}$$

これは H1 問題を完全に解いても、途中で打ち切っても成り立つ。

**Prager-Synge の定理の適用条件**:

定理が成り立つために必要なのは、$\mathbf{H}_{eq}$ が curl-free であること**のみ**。最適解である必要はない。

$$\|\mathbf{H} - \mathbf{H}_A\|^2 + \|\mathbf{H} - \mathbf{H}_{eq}\|^2 = \|\mathbf{H}_A - \mathbf{H}_{eq}\|^2$$

この等式から導かれる**誤差の上界**：

$$\|\mathbf{H} - \mathbf{H}_A\| \leq \|\mathbf{H}_A - \mathbf{H}_{eq}\|$$

**CG打ち切りの影響**:

| 性質 | 完全解 | CG打ち切り |
|:-----|:-------|:-----------|
| curl-free | ✓ | ✓（恒等式により保証） |
| 誤差上界の保証 | ✓ | ✓ |
| 上界のタイトさ | 最適 | やや緩い |

**結論**: CG を打ち切っても、$\nabla \Omega_{approx}$ は curl-free なので、**誤差の上界としての性質は数学的に保証される**。ただし、上界が最適（タイト）でなくなる可能性がある。

**AMR への影響**: 誤差の**絶対値**は多少過大評価されるが、**相対的な分布**（どの要素の誤差が大きいか）は保たれる。これがAMRに十分使える理由である。

#### (2) 数値実験結果

軸対称磁気球問題（μr = 100）での検証結果：

| CG反復数 | $\|\mathbf{H}_A - \mathbf{H}_{eq}\|$ | 直接法との差 |
|:---------|:-------------------------------------|:-------------|
| 1 | 0.966 | 18.3% |
| 5 | 0.928 | 13.7% |
| 10 | 0.885 | 8.4% |
| **20** | **0.826** | **1.25%** |
| 50 | 0.816 | 0.00% |

#### (3) AMR への適用性

要素ごとの誤差分布を比較（CG 20回 vs 直接法）：

| 指標 | 値 | 意味 |
|:-----|:---|:-----|
| **相関係数** | **0.976** | 誤差分布がほぼ一致 |
| **細分化オーバーラップ (top 10%)** | **97.7%** | ほぼ同じ要素を選択 |

**結論**: AMR の目的（どの要素を細分化するか）では、CG 10〜20回で十分な精度が得られる。

#### (4) 実装例

```python
from ngsolve.krylovspace import CGSolver

# Ω法の線形システムを構築
fes_O = H1(mesh, order=order, dirichlet="outer")
Omega, psi = fes_O.TnT()

a_O = BilinearForm(fes_O)
a_O += grad(Omega) * grad(psi) * dx
a_O.Assemble()

f_O = LinearForm(fes_O)
f_O += H_A * grad(psi) * dx  # H_A はA法の解
f_O.Assemble()

# 前処理器（Jacobi スムーザー）
pre = a_O.mat.CreateSmoother(fes_O.FreeDofs())

# CG を 20 回で打ち切り
gf_O = GridFunction(fes_O)
inv_cg = CGSolver(a_O.mat, pre, maxiter=20, tol=1e-16, printrates=False)
gf_O.vec.data = inv_cg * f_O.vec

# 平衡誤差推定量
H_eq = grad(gf_O)
eta = sqrt(Integrate(InnerProduct(H_A - H_eq, H_A - H_eq) * dx, mesh))
```

#### (5) 計算コストの比較

| 方法 | 計算コスト | 精度 | AMR適用性 |
|:-----|:-----------|:-----|:----------|
| **直接法** | O(n^1.5〜2) | 厳密 | ◎ |
| **CG 収束まで** | O(n × iter) | 厳密 | ◎ |
| **CG 20回打ち切り** | O(n × 20) | 1〜2%誤差 | ◎ |
| **CG 10回打ち切り** | O(n × 10) | 8%誤差 | ○ |

大規模問題では、CG 打ち切りが直接法より高速になる場合がある。

### 8.7. 平衡誤差推定の双対性と実用上の制約

#### (1) 理論的な双対性

Prager-Synge の定理は、A法とΩ法の両方向で成り立つ：

| 主問題 | Strong 性質 | Weak 性質 | 平衡化問題 |
|:-------|:------------|:----------|:-----------|
| **A法** | ∇·B = 0 | curl H = J | Ω法（curl-free な H_eq） |
| **Ω法** | curl H = 0 | ∇·B = 0 | A法（div-free な B_eq） |

つまり、**逆方向（Ω法 → A法での平衡化）も理論的には可能**である。

#### (2) 逆方向が実用的でない理由

しかし、Ω法を主問題としてA法で平衡化するアプローチは、以下の理由で**実用的ではない**：

| 観点 | A法 → Ω法 平衡化 | Ω法 → A法 平衡化 |
|:-----|:-----------------|:-----------------|
| **平衡化問題の空間** | H1（スカラー） | HCurl（ベクトル） |
| **自由度** | 1 DOF/節点 | 3 DOF/辺（3D） |
| **計算コスト** | **軽い** | **重い** |
| **実用性** | ◎ | △ |

#### (3) 根本的な問題

平衡誤差推定の意義は、**主問題より軽い計算で誤差を推定する**ことにある。

- **A法が主問題の場合**: Ω法（H1）で平衡化 → 主問題より軽い → **有効**
- **Ω法が主問題の場合**: A法（HCurl）で平衡化 → 主問題より重い → **意味がない**

#### (4) 結論

| シナリオ | 推奨アプローチ |
|:---------|:---------------|
| A法で解いた | Ω法で平衡誤差推定（本章の手法） |
| Ω法で解いた | 残差ベース誤差推定、または別の手法 |

**Ω法を主問題とした場合の平衡誤差推定（A法での平衡化）は、計算コストの観点から実用的ではない。** 残差ベースの誤差推定や、ZZ（Zienkiewicz-Zhu）推定量など、他の手法を検討すべきである。

---

## まとめ

1. **平衡誤差推定量は数学的に最良**: Prager-Synge の定理による等式評価、汎用定数なし
2. **Ω法の解が平衡誤差推定量を与える**: $\mathbf{H}_\Omega = \nabla \Omega$ は自動的に curl-free
3. **計算コストは軽い**: H1空間はHCurlより小さく、AMRの中では追加コストは無視できる
4. **CG打ち切りで高速化可能**: 10〜20回のCG反復でAMRに十分な精度が得られる
5. **双対性の制約**: 逆方向（Ω法→A法平衡化）は計算コストの観点から実用的ではない

## 参考文献

- Braess, D., & Schöberl, J. (2008). Equilibrated residual error estimator for edge elements. *Mathematics of Computation*, 77(262), 651-672.
- Prager, W., & Synge, J. L. (1947). Approximations in elasticity based on the concept of function space. *Quarterly of Applied Mathematics*, 5(3), 241-269.
- Bossavit, A. (1998). *Computational Electromagnetism: Variational Formulations, Complementarity, Edge Elements*. Academic Press.
