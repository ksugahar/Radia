# アーキテクチャ全体像

## 4 層構造 (CLAUDE.md の "4-Layer Architecture" に準拠)

```
┌─────────────────────────────────────────────────────────────────┐
│  Layer 1: C++ Qt5 (.ccl)          radia_cubit.ccl                │
│  ─────────────────────────────────────────────────────────────  │
│  Export Mesh menu (GMSH/Nastran/VTK/Netgen Vol/FEMEEM/MEG)      │
│  Mesh Evaluation (_p1.vol ... _p5.vol)                           │
│  radia_export netgen/gmsh/nastran/vtk (APREPRO コマンド)          │
└─────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────┐
│  Layer 2: Cubit Python 3.10 + PySide6/Qt5                        │
│  ─────────────────────────────────────────────────────────────  │
│  register_toolbar.py -> Solve menu (IH panel 起動)               │
│  import cubit OK。import radia/ngsolve は 禁止 (DLL 競合)         │
└─────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────┐
│  Layer 3: PySide6 Panel (Python 3.12 + PySide6)                 │
│  ─────────────────────────────────────────────────────────────  │
│  radia_ih.py (IHWindow) — standalone PySide6 アプリ               │
│  Cubit とは別プロセス。import cubit は禁止。                       │
│  .vol ファイル (および STEP) をパラメータとして受け取る           │
└─────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────┐
│  Layer 4: Computation (Python 3.12, no GUI)                     │
│  ─────────────────────────────────────────────────────────────  │
│  calc_peec_bem.py (PEEC+BEM 1-way)                               │
│  calc_fem_coilmesh.py (FEM A-V)                                  │
│  import cubit 禁止。import PySide6 禁止。                         │
│  NGSolve + GMSH のみ。JSON stdout で結果返す。                    │
└─────────────────────────────────────────────────────────────────┘
```

**層間インターフェース**: `.vol` ファイル (text 形式、ABI 依存なし) + `.step` (STEP 3D CAD)。

## IH panel の 2 つの method

### (a) PEEC+BEM (1-way, P-focus)

**定式化**: PEEC (partial element equivalent circuit) フィラメントでコイルモデル化 + Scalar BIE SIBC で表面電位 φ を解く。

```
STEP coil file
 ├→ filaments_from_step → 9 本のフィラメント (nwinc=3, nhinc=3)
 ├→ PEEC MNA solve → I_fil (各フィラメント電流、複素数)
 │
.vol (wp hole + sibc sideset)
 ├→ _extract_bnd_only → sibc 表面メッシュ
 │   [重要] centroid-based winding 修正で normal 方向統一
 │
I_fil + sibc 表面
 ├→ compute_phi_inc_from_filaments (Biot-Savart 線積分)
 └→ BEM SIBC solve → φ on wp surface → H_t → P_wp
```

**成果物**: P_wp (主)、L_coil (空芯、参考値)  
**所要時間**: 3 分 (BEM assembly が O(N²)、wp mesh 595 nv)  
**利点**: コイル mesh 不要、wp も粗い surface mesh で OK  
**限界**: ワーク back-reaction (L) は未実装 (1-way = forward のみ)

### (b) FEM A-V (coil meshed + wp SIBC + Kelvin)

**定式化**: A-V compound FES (HCurl A + H1 phi on coil) で volumetric eddy を解く。

```
.vol (coil 体積 + source/sink port + wp hole/sibc + kelvin 外部)
 ├→ HCurl(A) x H1(phi on coil, Dirichlet source|sink)
 ├→ BF:
 │    nu * curl(A).curl(N) * dx                     (all domains)
 │    1e-6 * NU_0 * A.N * dx ('|'.join(非 kelvin))    (gauge)
 │    robin_wp * A.Trace().N.Trace() * ds(sibc)     (wp SIBC)
 │    s*sigma*(A+grad(phi))*(N+grad(psi))*dx('coil') (A-V)
 ├→ Dirichlet lift: phi = 1 on source, 0 on sink
 ├→ pardiso solve
 │
 ├→ 体積積分 current 抽出 (MCP 既知パターン):
 │    psi_n scalar H1 on coil, dirichlet=sink, Set(1) on source
 │    I_out = int J.grad(psi_n) dx('coil')    # Gauss consistent
 │    scale = I_target / I_out ⇒ gfu *= scale
 │
 ├→ L = 2 * W_m / I² + skin layer 項
 ├→ P_wp = 0.5 * Re(Z_s) * H_t_rms² * A_wp   (SIBC)
 └→ P_coil = 0.5 / sigma * int |J|² * dx('coil')   (volumetric, mesh-sensitive)
```

**成果物**: L (wp 込み)、P_wp、P_coil、I_out、R_total  
**所要時間**: 7-9 分 (ndof 874k、pardiso 54s + post 数秒)  
**利点**: 物理的に一番完全 (場も電流も体積で解く)  
**要求**: コイル体積 mesh (表面 ≤ δ/3、IH 用 sample は ~0.16mm)、source/sink port sideset

## 捨てた定式化 (歴史的経緯)

2026-04-18〜19 に複数 trial:

| Trial | 結果 | 退場理由 |
|---|---|---|
| T0 (Laplace scalar coil source) | 1/r cusp at gap corners | gapped geom で不安定 |
| Neumann-K coil surface + Robin SIBC | L collapse 10x | Z_s=(1+j)ρ/δ で K_response=-K 相殺 |
| Neumann-K no Robin | P_wp OK だが BEM と 24% 差 | uniform K が proximity 未反映 |
| Dowell 解析的 P_coil | 正確だが mesh 依存性を隠す | panel 降格、研究用 examples/ 残 |

**教訓**: 新しい formulation は必ず golden test + 2D ref との validation で lock。面白そうな idea でも、ladder (`tests/ → examples/ → panels/`) を経ずに panel に突っ込まない。

## 今の 2 methods を選んだ理由

- **User の primary metric = P_wp** (ワーク加熱量)。2 method で 1% 一致という強い相互 validation
- **Gapped torus の source/sink port は実機 IH の現実**。真の IH コイルは物理端子を持つ
- **MCP に既知の A-V + volume current extraction パターン** が効く → 再発明せず再利用

## 参考

- `memory/ih_panel_final_state_2026_04_19.md` — 定式化経緯の詳細
- `memory/feedback_panel_ladder.md` — 3段 ladder 原則
- `memory/feedback_demo_to_production_lock.md` — demo→production 移植ルール
- MCP `ih_knowledge.peec_bem_sibc` と `av_coil_sigma` — 物理詳細
