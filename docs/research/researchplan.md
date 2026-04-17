# Research Plan: LLM-Driven Coil Design with Full-OSS CAE Pipeline

Date: 2026-04-14 (initiated), updated through brainstorming session
Status: Planning / Brainstorming

## Core Thesis

**CAD 形状を先に作り、そこから解析モデルを自動生成する。**
LLM がこのループを自律的に回し、断面可変・曲率連続なコイルを自動設計する。

全工程 OSS。論文のコードを pip install だけで完全再現可能。

```
build123d (CAD, OSS)
    → CoilGeometry（CAD → 積分モデル自動変換）
    → Netgen + 高次カーブ要素 (mesh, OSS)
    → NGSolve + SIBC (solve, OSS)
    ↑ MCP servers (Radia, OSS) で LLM が全工程を駆動
```

---

## Background

### The real bottleneck in electromagnetic CAE

The research community optimizes for **solver speed and accuracy**, but the actual time breakdown in industrial electromagnetic design is:

- Model creation (CAD, defeaturing, healing): **days**
- Mesh generation and quality tuning: **hours**
- Solver execution: **minutes**

This mismatch is invisible to academics (students build the models), but obvious to anyone from industry (Sugahara: ex-Mitsubishi Electric).

### Existing work in the group

- **Mesh AI** (Shimizu, Compumag 2025): LLM-driven mesh generation
- **Radia MCP servers** (8 types, PyPI-distributed): LLM-driven CAE execution
- **3D printed coils** (GoTech collaboration): additive-manufactured accelerator magnet coils
- **SIBC + 高次カーブ要素** (比留間 → Joachim): 形状再現性の高い表面インピーダンス境界条件、ほぼ完成
- **ESIM-SIBC** (Hollaus, JSPS 2026): 積層鋼板の等価表面インピーダンス

### Three missing pieces — now all addressable

1. **CAD 自動生成**: build123d + mcp-server-build123d（2026-04-14 実装済み）
2. **CAD → 積分モデル自動変換**: CoilGeometry クラス（本計画の中核）
3. **自動設計ループ**: LLM が MCP 経由で全工程を駆動

---

## The Paradigm Shift: CAD First, Analysis Model Second

### 従来（CoilBuilder の設計思想）

```
arc + straight で電流パスを定義（積分モデル）
    → to_radia(): Biot-Savart 源
    → to_occ(): CAD 形状（後付け、おまけ）
```

**積分モデルが先、CAD が後。** 形状は近似的で、arc + straight に制限される。

### 新しい発想

```
build123d で物理形状を自由に作る（CAD が真実）
    → CoilGeometry: CAD から積分モデルを自動生成
    → 加速器: Biot-Savart ワイヤーセグメント
    → IH: 表面電流 SIBC モデル
    → 製造: STL/STEP をそのまま出力
```

**CAD が先、積分モデルは派生。**

### なぜ逆転すべきか

| | 積分モデル → CAD（従来） | CAD → 積分モデル（提案） |
|---|---|---|
| 形状自由度 | arc + straight に制限 | **任意形状**（スプライン、loft、断面可変） |
| 断面 | 矩形一定のみ | **可変、任意形状** |
| エンド部 | 円弧（曲率不連続） | **スプライン（C² 連続、物理的に正しい）** |
| CAD の品質 | 後付け → 計算モデルとずれうる | **CAD が唯一の真実** |
| 3D プリンタ | CAD を別途作り直す | **そのまま STL 出力** |
| LLM 連携 | arc/straight パラメータ列を生成 | **build123d スクリプトを生成** |

---

## CoilGeometry: Unified Coil Representation

### 設計

```python
class CoilGeometry:
    """CAD 形状から解析モデルを自動生成する統一コイル表現。
    
    入力: build123d Solid + 電流パラメータ
    出力: 用途に応じた解析モデル
    """
    
    # --- 入力 ---
    
    @staticmethod
    def from_build123d(solid, current, freq=0, sigma=5.8e7, centerline=None):
        """build123d Solid からコイルを生成（汎用）"""
        ...
    
    @staticmethod
    def from_coil_builder(coil_builder):
        """既存 CoilBuilder からも生成可能（後方互換）"""
        ...
    
    # --- 加速器用出力（DC） ---
    
    def to_wire_segments(self, n=100):
        """DC Biot-Savart 用ワイヤーセグメント
        中心線に沿って直線素片に分割"""
        ...
    
    def to_radia(self):
        """Radia Biot-Savart オブジェクト"""
        ...
    
    # --- IH 用出力（AC） ---
    
    def to_surface_current(self):
        """高周波近似: 導体表面の電流シート
        build123d Solid の外表面 → SIBC モデル"""
        ...
    
    def to_mesh_region(self, label="coil"):
        """コイルを FEM 解析領域として返す（メッシュ対象）
        SIBC or フル FEM 用"""
        self.solid.label = label
        return self.solid
    
    # --- 製造用出力 ---
    
    def export_step(self, path): ...
    def export_stl(self, path): ...
```

### CoilBuilder との関係

CoilBuilder は廃止しない。位置づけが変わる：

```
旧: CoilBuilder = コイルの唯一の定義方法
新: CoilBuilder = 古典的な記述方法（arc+straight、後方互換）
    build123d   = 汎用的な記述方法（任意形状）
    どちらも → CoilGeometry → 共通の解析モデル出力
```

---

## Coil Physics: Why Variable Cross-Section and Continuous Curvature

### 従来のコイル記述の問題

arc + straight は「CAD が手間だから」の近似。物理的根拠はない。

- **円弧のエンド部**: 曲率が不連続（直線 → 円弧の接続点で jump）
  → 応力集中点。3D プリンタ製コイルでは物理的に影響する
- **断面一定**: 製造制約（巻線機、CICC）であって物理的最適ではない
  → 電流密度、構造強度、冷却が位置によって異なるのに断面が一定

### build123d で解決

- **スプラインパス**: C² 連続 → 曲率連続 → 応力集中なし
- **loft による断面変化**: 位置の関数として断面形状を定義
- **パラメトリック**: Python で物理式から直接形状を生成

```python
# 物理式から断面を決定
for t in np.linspace(0, 1, 100):
    z = path.position_at(t)
    force = em_force_at(z)          # 電磁力
    j_target = j_max / safety(force) # 許容電流密度
    area = I / j_target              # 必要断面積
    sections.append(path.location_at(t) * Rectangle(w(area), h(area)))
coil = loft(sections)
```

### IH コイルの電流分布（AC の場合）

IH コイルは 1〜数ターンの太い導体。kHz〜MHz 帯では表皮深さ δ << 導体寸法。

電流分布を決める 3 つの効果（重要度順）：

| 効果 | 原因 | DC でも？ | 記述 |
|------|------|----------|------|
| **1. 経路長効果** | コイルの曲率 → 内側が短い → J ∝ 1/r | Yes | 幾何学的 |
| **2. 表皮効果** | 自己磁場 → 表面に電流集中 | No | Bessel / exp(-d/δ) |
| **3. 近接効果** | ワークピース・隣接ターンの磁場 | No | 連立 / SIBC |

**高周波では 3 つとも表面電流近似に帰着する。** 導体表面のどこに電流が集中するか、の問題。

### 解析モデルの階層

| Level | 近似 | 用途 | 速度 |
|-------|------|------|------|
| 0 | DC Biot-Savart（J 一様） | 初期設計 | 最速 |
| 1 | + 経路長効果（1/r） | DC 最適化 | 速い |
| 2 | + 表皮効果（表面電流） | AC 最適化 | 速い |
| 3 | 高次カーブ要素 + SIBC（近接効果込み） | **最終検証** | 遅い |

**最適化ループは Level 0-2 で回し、最終検証だけ Level 3。**

---

## SIBC + 高次カーブ要素: The Precision Layer

### 背景

比留間先生が XFEM で SIBC を研究 → **形状再現性が精度の鍵** と判明
→ Joachim（Netgen/NGSolve 開発者）に高次カーブ要素を依頼 → ほぼ完成

### なぜ高次カーブ要素が必要か

SIBC は導体表面での境界条件。表面の形状が不正確だと：
- 法線方向が狂う → 表面インピーダンスの方向が狂う
- 曲率が不正確 → 表面電流分布が不正確
- **低次要素の直線近似では不十分**

高次カーブ要素 → 曲面を正確に表現 → SIBC の精度が出る

### build123d との連携

```
build123d: コイル表面を正確に定義（OCCT BREP）
    → STEP/BREP → Netgen（高次カーブ要素でメッシュ）
    → NGSolve + SIBC で表面電流分布を計算
```

**build123d の CAD 品質 → Netgen の高次要素品質 → SIBC の解析精度。**
CAD が先であることの価値が、ここで最大化される。

### JSPS (Hollaus) との接続

- Hollaus の ESIM: 積層鋼板の等価表面インピーダンス
- build123d で cell-problem geometry を Cubit なしで生成可能
- Full OSS stack: build123d + Netgen + NGSolve + ESIM-SIBC + Radia

---

## Paper: LLM-Driven Coil Design with Variable Cross-Section

### Story

1. **Problem**: 3D printing removes manufacturing constraints on coil geometry.
   Cross-section can vary, path curvature can be continuous.
   But the design parameter space is too large for human exploration.

2. **Key insight**: CAD should come first, analysis model second.
   CoilGeometry class automatically derives Biot-Savart sources (DC)
   or surface current SIBC models (AC) from build123d CAD shapes.

3. **Method**: LLM orchestrates a fully-OSS pipeline via MCP servers:
   build123d (CAD) → Netgen + curved elements (mesh) → NGSolve + SIBC (solve)
   → evaluate → modify CAD → iterate.

4. **Demonstration**:
   - Accelerator dipole: variable cross-section + spline path → field harmonics optimization
   - IH coil: variable cross-section → heating pattern optimization with SIBC validation

5. **Evaluation**:
   - Field quality: multipole harmonics (accelerator) / power density distribution (IH)
   - 3D printability: minimum wall thickness, overhang angle
   - Convergence: iterations to specification, where LLM makes wrong judgments
   - Comparison: vs. constant-cross-section baseline, vs. expert design

### Why this paper is strong

- **Novel**: CAD-first coil design with automatic analysis model derivation
- **Timely**: 3D printed magnets are a hot topic (CERN, BNL, GoTech)
- **Fully reproducible**: entire pipeline is OSS (pip install で完全再現)
- **Two applications**: accelerator (DC) + IH (AC) from the same CoilGeometry
- **Precision story**: build123d → 高次カーブ要素 → SIBC の精度連鎖が明確
- **Interesting whether it succeeds or fails**: failure analysis reveals LLM limits

### All-OSS pipeline (critical differentiator)

```
build123d (Apache-2.0) → Netgen (LGPL) → NGSolve (LGPL) → Radia MCP (BSD-3)
    CAD                    mesh              solve             orchestration
```

Everything installs with `pip install`. Paper includes complete code.
**Reader reproduces the exact same design with zero license cost.**

Cubit (commercial) is optional: hex mesh for higher accuracy,
but the core pipeline does not depend on it.

### Target venues

- COMPUMAG 2027 (primary target — full paper)
- CEFC 2026 (if timing works — digest)
- IEEE Transactions on Magnetics (journal, after conference)
- IEEJ 静止器・回転機合同研究会 (Japanese domestic, preliminary results)

### Application scope

| Domain | Current practice | Variable cross-section + continuous curvature benefit |
|--------|-----------------|------------------------------------------------------|
| Accelerator (GoTech) | 3D printed, constant cross-section, arc+straight | Field quality + structural + spline path (no stress concentration) |
| IH | Wound, constant cross-section | Heating pattern control, surface current optimization via SIBC |
| Fusion (stretch) | CICC, constant cross-section | Weight reduction, integrated cooling, force management |

---

## Kubota's Theme

### テーマ

**LLM × MCP × build123d による電磁機器コイルの自動設計**

### Phase 1: Foundation（M1 前期）
- build123d でのコイル CAD 生成を習得（レーストラック、スプラインパス）
- CoilGeometry クラスの実装（from_build123d + to_wire_segments）
- MCP server 経由の build123d → Netgen → NGSolve パイプライン構築

### Phase 2: Variable cross-section + AC model（M1 後期）
- 断面可変コイルの build123d loft 生成
- CoilGeometry に AC 出力追加（表面電流、経路長効果）
- IH コイルで加熱パターン評価

### Phase 3: Automated design loop（M2 前期）
- LLM が設計仕様 → CAD → 解析 → 評価 → 修正のフルループを自律実行
- 加速器ダイポール + IH で実証
- 高次カーブ要素 + SIBC による最終検証
- 論文執筆

### Phase 4: Fusion application（M2 後期 / stretch）
- トカマク TF コイルへの適用検討
- 大型コイルでのスケーラビリティ

### Positioning

- Mesh AI (清水, Compumag 2025): **mesh generation** automation
- Coil Design AI (久保田, 2027-): **coil geometry + analysis model** automation
- 「Mesh AI がメッシュを自動化した。本研究はコイル設計そのものを自動化し、
  CAD から解析モデルの自動生成まで含む」

---

## MCP Server: mcp-server-build123d (2026-04-14, implemented)

### Background

Web 調査の結果、**build123d の API を tool として expose する MCP server は世界に存在しない**
ことが確認された（2026-04-14 時点）。存在するのは:
- **ocp-viewer-mcp** (dmilad): OCP CAD Viewer のスクリーンショット取得のみ（「目」）
- **CAD Agent** (Svetlana-DAO-LLC): Docker コンテナ内の build123d 実行環境

Note: Gemini が「xintlabs/build123d-mcp が存在する」と主張したが、
GitHub 404 確認済み。ハルシネーション。

### 実装済み: mcp-server-build123d

Radia MCP server と同一パターン（FastMCP + knowledge base + selftest）で実装。

**Tools (3つ):**

| Tool | 機能 |
|------|------|
| `build123d_usage(topic)` | CAE 向け API リファレンス（9 トピック） |
| `execute_build123d(script)` | スクリプト実行 + 品質検証（volume, area, min edge length, CAE warnings） |
| `inspect_geometry(file_path)` | 既存 STEP/BREP の CAE 品質診断 |

**Files:**
- `src/radia/mcp_server/build123d/server.py`
- `src/radia/mcp_server/build123d/build123d_knowledge.py`
- Entry point: `mcp-server-build123d` (packages/radia-mcp)

### CAE-Safe Subset

build123d は何でもできるが、CAE には「きれいな CAD」が必要。
Cubit が「汚い CAD を作りにくい」のは便利機能がないから（= 制約が品質を守る）。
build123d でも **プリミティブ + boolean に限定すれば Cubit と同等の品質**。

ただしコイル設計では loft, sweep, spline が本質的に必要。
CAE-safe の判定を execute_build123d の自動検証（min edge length 等）で担保する。

---

## build123d as the Lab's Standard CAD

### Rationale

Cubit Python でも build123d でも形状は書ける。差は：
- **Cubit**: 有償。卒業後使えない。論文の再現性がない
- **build123d**: OSS (Apache-2.0)。pip install。論文にコードを載せられる

### Cubit の役割（維持）

Cubit は hex メッシャー + 品質保証モデラーとして維持。
ただしコア・パイプラインは Cubit に依存しない。

### Migration strategy

- **build123d**: standard for all geometry creation
- **netgen.occ**: deprecated for geometry creation (OCCGeometry loader only)
- **Cubit**: optional hex mesher
- **Existing .jou / netgen.occ code**: LLM-driven reconstruction to build123d

---

## Open questions

- **Label round-trip**: build123d Face `label` → STEP → netgen.occ で生き残るか
- **Coil loft quality**: 断面数を増やしたとき OCCT がどこまで滑らかにつなぐか
- **Centerline extraction**: build123d Solid からコイル中心線を自動抽出する方法
- **Surface current accuracy**: 表面電流 Biot-Savart 近似と SIBC の差はどの程度か
- **LLM の物理理解**: 多極展開と鉄飽和（加速器）/ 加熱パターンと形状（IH）をどこまで理解できるか
- **Fusion scale**: 大型コイル（数m級）での build123d/OCCT のメモリ・精度

## Next actions

- [ ] レーストラックコイルの build123d 生成（スプラインエンド、断面一定）
- [ ] CoilGeometry クラス設計・プロトタイプ実装
- [ ] build123d → Netgen → NGSolve 一気通貫テスト（磁場解析）
- [ ] 断面可変コイルの build123d loft テスト
- [ ] 表面電流 Biot-Savart 近似の実装と SIBC 結果との比較
- [ ] LLM に「IH コイルを設計して」と投げるデモ
- [ ] Kubota テーマ定義書
- [ ] GoTech との接点整理（既存 3D プリンタコイルデータ）
