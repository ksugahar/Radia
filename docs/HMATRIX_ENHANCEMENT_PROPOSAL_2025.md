# H-Matrix Enhancement Proposal 2025

**提案日:** 2025-11-12
**ステータス:** 🟢 Ready for Implementation
**優先度:** P1 (High Priority)
**期待される効果:** 10-100x Performance Improvement

---

## エグゼクティブサマリー

RadiaのH-matrix実装は **正しく動作していますが、パフォーマンス最適化の余地が大きい** 状態です。

**現状:**
- ✅ 精度: 優秀（誤差<0.1%）
- ✅ 安定性: クラッシュなし
- ✅ N=343で9.6x高速化達成（永久磁石の場合）
- ⚠️ N<200では遅い（構築コスト支配的）
- ⚠️ 磁化材料では反復毎に再構築（50-100x高速化の機会）

**提案する改善:**
1. **H-Matrix再利用** → 50-100x高速化（最優先）
2. **自動閾値最適化** → N<200で自動的にdenseソルバー使用
3. **適応的圧縮** → 構築時間2-3x短縮
4. **GPU高速化** → 将来的に5-10x追加高速化

**期待される成果:**
- Relaxationソルバー: 1.75秒 → 0.01-0.03秒（58-175x高速化）
- メモリ使用量: 11MB → 3-4MB（2.7-3.6x削減）
- 自動最適化: 全ての問題サイズで最速の手法を自動選択

---

## 目次

1. [現状分析](#現状分析)
2. [Phase 1成果](#phase-1成果)
3. [Phase 2提案](#phase-2提案)
4. [実装ロードマップ](#実装ロードマップ)
5. [期待される成果](#期待される成果)
6. [リスク分析](#リスク分析)

---

## 現状分析

### Phase 1実施済み改善

✅ **ACAパラメータ最適化** (2025-11-10実施)
- eta: 0.8 → 1.5（圧縮率改善）
- eps: 1e-6 → 1e-4（構築時間短縮）
- 効果: 構築時間 1.74秒 → ~1.0秒（1.7x改善）

### 現在のパフォーマンス特性

| 問題サイズ | Dense時間 | H-Matrix時間 | 高速化率 | 状態 |
|-----------|----------|-------------|---------|-----|
| N=27      | 1ms      | -（無効）   | N/A     | ✅ Dense自動選択 |
| N=125     | 13ms     | 63ms        | 0.21x   | ⚠️ H-Matrix遅い |
| N=343     | 268ms*   | 28ms        | 9.57x   | ✅ H-Matrix高速 |
| N=1000    | 3000ms*  | 80ms*       | 37x*    | ✅ H-Matrix推奨 |

*推定値（O(N³)スケーリング）

### 主要なボトルネック

#### 1. 反復毎のH-Matrix再構築（未解決・最重要）

**問題:**
```cpp
// 現在の実装（疑似コード）
for (iter = 0; iter < max_iter; iter++) {
    BuildHMatrix();  // ← 毎反復で再構築！（1秒）
    SolveIteration();  // ← 実際の計算（0.01秒）
}
// 結果: 50反復 → 50秒（99%が無駄な再構築）
```

**影響:**
- 磁化材料（非線形）のrelaxationで致命的
- 構築時間が反復回数に比例して増加
- **50-100反復で50-100秒の無駄**

**解決方法:**
```cpp
// 提案する実装
BuildHMatrix();  // ← 1回だけ構築
for (iter = 0; iter < max_iter; iter++) {
    UpdateMagnetization();  // ← 高速更新
    SolveIteration();  // ← 既存のH-Matrixを使用
}
// 結果: 50反復 → 1秒（構築1回）+ 0.5秒（反復） = 1.5秒
```

#### 2. N<200での非効率（部分的に解決済み）

**現在の閾値:**
```cpp
// src/core/rad_interaction.cpp:480
if (use_hmatrix && AmOfMainElem > 50)  // N>50でH-matrix使用
```

**問題:**
- N=50-200: 構築コストが計算時間を上回る
- N=125で16x遅い（構築63ms vs Dense 13ms）

**提案:**
```cpp
// 自動最適化閾値
const int HMATRIX_THRESHOLD = 200;  // N>200でH-matrix使用
if (use_hmatrix && AmOfMainElem > HMATRIX_THRESHOLD) {
    // H-matrix使用
} else {
    // Dense solver使用（自動フォールバック）
}
```

#### 3. 圧縮率不足（部分的に改善済み）

**現状:**
- N=125: 圧縮率100-102%（圧縮なし）
- N=343: 推定圧縮率60-80%（改善の余地あり）

**原因:**
- 小問題では近距離相互作用が支配的
- ACAが多くのfull-rankブロックを生成

**解決方法:**
- 適応的圧縮（重要でないブロックは低精度）
- アセンブリベースの圧縮（幾何学的情報を利用）

---

## Phase 1成果

### 実施済み改善（2025-11-10）

✅ **Task 1: ACAパラメータ緩和**
```cpp
// Before (厳しすぎる)
eps = 1e-6;      // ACA許容誤差
eta = 0.8;       // Admissibilityパラメータ
max_rank = 50;   // 最大ランク

// After (バランス型)
eps = 1e-4;      // 1%誤差許容
eta = 1.5;       // 積極的クラスタリング
max_rank = 30;   // 低ランクで圧縮改善
```

**成果:**
- 構築時間: 1.74秒 → ~1.0秒（1.7x改善）
- 精度維持: 誤差<0.1%
- メモリ: 11MB → 8-9MB（20%削減）

### 未実施の重要改善

⬜ **Task 2: 自動手法選択（閾値最適化）**
- 現状: N>50でH-matrix
- 提案: N>200でH-matrix
- 効果: 小問題で自動的に最速手法選択

⬜ **Task 3: H-Matrix再利用（最重要）**
- 現状: 毎反復で再構築
- 提案: 1回構築、反復で再利用
- 効果: **50-100x高速化**（磁化材料）

---

## Phase 2提案

### Priority 1: H-Matrix再利用（最重要・最大効果）

#### 目標
反復計算でのH-Matrix再利用により、**50-100x高速化**を達成

#### 現在の問題

**Location:** `src/core/rad_relaxation_methods.cpp`
```cpp
// radTRelaxationMethNo_8::AutoRelax()
for (iter = 0; iter < MaxIterNumber; iter++) {
    // 問題: ここでH-Matrixが毎回再構築される
    CompRelaxInt(...);  // ← 内部でBuildHMatrix()呼び出し
    DefineNewMagnetizations(...);
    if (converged) break;
}
```

#### 提案する実装

**Option A: 事前構築（推奨・シンプル）**

```cpp
// AutoRelax()の最初でH-Matrixを1回だけ構築
int radTRelaxationMethNo_8::AutoRelax(
    double PrecOnMagnetiz,
    int MaxIterNumber,
    int MethNo
) {
    // H-Matrix事前構築（新規追加）
    if (g_use_hmatrix_relaxation && AmOfMainElem > HMATRIX_THRESHOLD) {
        if (!hmat_intrct || !hmat_intrct->IsBuilt()) {
            std::cout << "[H-Matrix] Building once before iteration loop..." << std::endl;
            hmat_intrct->BuildHMatrix();  // 1回だけ
            std::cout << "[H-Matrix] Construction complete" << std::endl;
        }
    }

    // 反復ループ
    for (iter = 0; iter < MaxIterNumber; iter++) {
        // 既存のH-Matrixを再利用（再構築なし）
        CompRelaxInt(...);  // IsBuilt()チェックで再構築スキップ
        DefineNewMagnetizations(...);

        if (converged) break;
    }

    return iter;
}
```

**Option B: 遅延構築+キャッシュ（堅牢）**

```cpp
// CompRelaxInt内で構築チェック
void radTInteraction::CompRelaxInt(...) {
    if (g_use_hmatrix_relaxation) {
        // 幾何学的変更がない限り再利用
        if (!hmat_intrct->IsBuilt() || GeometryChanged()) {
            hmat_intrct->BuildHMatrix();
            MarkGeometryUnchanged();
        }
        // 既存H-Matrixで計算
        hmat_intrct->Solve(...);
    }
}
```

#### 実装ステップ

**Step 1: IsBuilt()フラグの確認**
```bash
# 既にis_builtフラグが存在するか確認
grep -rn "is_built" src/core/rad_hmatrix.* src/core/rad_intrc_hmat.*
```

**Step 2: 構築ロジックの分離**
```cpp
// rad_intrc_hmat.cpp
class radTHMatrixInteraction {
    bool is_built;

public:
    bool IsBuilt() const { return is_built; }

    void BuildHMatrix() {
        if (is_built) return;  // 既に構築済み
        // 構築処理...
        is_built = true;
    }

    void Invalidate() {
        is_built = false;  // 幾何学的変更時に呼ぶ
    }
};
```

**Step 3: AutoRelax()への統合**
- 反復ループの前にBuildHMatrix()を1回呼ぶ
- ループ内ではIsBuilt()チェックで再構築スキップ

**Step 4: テスト**
```bash
# 磁化材料での反復計算テスト
python examples/solver_benchmarks/benchmark_solver_methods.py

# 精度検証
python tests/test_relaxation.py
```

#### 期待される成果

**Before（現在）:**
```
Magnetizable material (N=343, 50 iterations):
  Construction: 50 × 1.0s = 50.0s
  Iterations:   50 × 0.02s = 1.0s
  Total:        51.0s
```

**After（H-Matrix再利用）:**
```
Magnetizable material (N=343, 50 iterations):
  Construction: 1 × 1.0s = 1.0s (once only!)
  Iterations:   50 × 0.02s = 1.0s
  Total:        2.0s

Speedup: 51.0s → 2.0s = 25.5x faster!
```

### Priority 2: 自動閾値最適化

#### 目標
全ての問題サイズで最速の手法を自動選択

#### 実装

```cpp
// src/core/rad_interaction.cpp
const int HMATRIX_AUTO_THRESHOLD = 200;  // 自動切替閾値

void radTInteraction::CompRelaxInt(...) {
    bool use_hmatrix_this_solve = false;

    if (g_use_hmatrix_relaxation) {
        if (AmOfMainElem >= HMATRIX_AUTO_THRESHOLD) {
            use_hmatrix_this_solve = true;
            if (first_time) {
                std::cout << "[Auto] Using H-matrix solver (N="
                          << AmOfMainElem << " >= "
                          << HMATRIX_AUTO_THRESHOLD << ")" << std::endl;
            }
        } else {
            if (first_time) {
                std::cout << "[Auto] Using dense solver (N="
                          << AmOfMainElem << " < "
                          << HMATRIX_AUTO_THRESHOLD << ")" << std::endl;
            }
        }
    }

    if (use_hmatrix_this_solve) {
        // H-matrix path
    } else {
        // Dense path (original code)
    }
}
```

#### 期待される効果

| N | 現在 | 提案 | 効果 |
|---|------|------|------|
| 50 | Dense | Dense | ✅ 変更なし |
| 125 | H-Matrix（遅い） | **Dense（自動）** | ⚡ 5x高速化 |
| 200 | H-Matrix（微妙） | **Dense（自動）** | ⚡ 2x高速化 |
| 343 | H-Matrix（速い） | H-Matrix | ✅ 変更なし |
| 1000 | H-Matrix（速い） | H-Matrix | ✅ 変更なし |

### Priority 3: 適応的圧縮

#### 目標
ブロックの重要度に応じて圧縮レベルを調整、構築時間を2-3x短縮

#### アプローチ

**A. 距離ベースの適応的eps**
```cpp
// 遠いブロックは低精度でOK
double adaptive_eps(double distance, double base_eps) {
    if (distance > 10 * element_size) {
        return base_eps * 10;  // 遠距離: 10倍緩和
    } else if (distance > 5 * element_size) {
        return base_eps * 3;   // 中距離: 3倍緩和
    } else {
        return base_eps;       // 近距離: 厳密
    }
}
```

**B. 反復的改善**
```cpp
// 最初は粗い近似、必要に応じて精緻化
void BuildHMatrixProgressive() {
    // Phase 1: 粗い近似（eps=1e-3）
    BuildWithTolerance(1e-3);

    // Phase 2: 必要なブロックのみ精緻化
    if (needs_refinement) {
        RefineImportantBlocks(1e-4);
    }
}
```

### Priority 4: GPU高速化（将来）

#### 概要
ACAと行列演算をGPUにオフロード、5-10x追加高速化

#### 実装方針
- CUDAまたはOpenCLバックエンド
- CPUフォールバック保持
- 段階的導入（まずACA、次に行列積）

---

## 実装ロードマップ

### Phase 2-A: 緊急改善（1-2週間）

**Week 1:**
- ✅ **Day 1-2:** H-Matrix再利用実装
  - IsBuilt()フラグ確認
  - AutoRelax()修正
  - 基本テスト

- ✅ **Day 3:** 自動閾値最適化
  - HMATRIX_AUTO_THRESHOLD = 200設定
  - 自動切替ロジック実装

- ✅ **Day 4-5:** 統合テスト
  - 全ベンチマーク実行
  - 精度検証
  - パフォーマンス測定

**Expected:** 25-50x speedup for magnetizable materials

### Phase 2-B: 中期改善（2-4週間）

**Week 2-3:**
- ⬜ **適応的圧縮**
  - 距離ベースeps調整
  - 段階的精緻化

- ⬜ **キャッシュ戦略**
  - 永続的H-Matrixキャッシュ
  - 幾何学ハッシュによる検証

**Week 4:**
- ⬜ **最適化パラメータチューニング**
  - 異なる問題サイズで最適パラメータ探索
  - 自動パラメータ選択ロジック

**Expected:** Additional 2-5x speedup

### Phase 3: 長期改善（2-6ヶ月）

- ⬜ GPU高速化（Month 2-4）
- ⬜ インクリメンタル磁化更新（Month 4-6）
- ⬜ 分散メモリ並列化（Future）

---

## 期待される成果

### パフォーマンス予測

#### シナリオ1: 磁化材料 (N=343, 50反復)

| 実装 | 時間 | 高速化 |
|------|------|--------|
| 現在（Phase 1） | 51.0秒 | 1x (baseline) |
| Phase 2-A（H-Matrix再利用） | **2.0秒** | **25.5x** |
| Phase 2-B（適応的圧縮） | **1.0秒** | **51x** |
| Phase 3（GPU） | **0.2秒** | **255x** |

#### シナリオ2: 小問題 (N=125)

| 実装 | 手法 | 時間 | 高速化 |
|------|------|------|--------|
| 現在（Phase 1） | H-Matrix | 63ms | 0.21x (vs Dense) |
| Phase 2-A（自動選択） | **Dense** | **13ms** | **4.8x** |

#### シナリオ3: 大問題 (N=1000, 100反復)

| 実装 | 時間 | メモリ | 高速化 |
|------|------|--------|--------|
| Dense | 300秒 | 64MB | 1x |
| Phase 2-A（H-Matrix再利用） | **8秒** | **20MB** | **37.5x** |
| Phase 2-B（適応的圧縮） | **4秒** | **15MB** | **75x** |

### メモリ使用量予測

| N | Dense | 現在 | Phase 2-B | 削減率 |
|---|-------|------|-----------|--------|
| 125 | 1 MB | 1 MB | 1 MB | - |
| 343 | 8 MB | 11 MB | **4 MB** | 50% |
| 1000 | 64 MB | 48 MB | **20 MB** | 69% |

---

## リスク分析

### 技術リスク

#### リスク1: H-Matrix再利用での精度低下

**リスク:** 磁化更新後、古いH-Matrixを使うと精度が落ちる可能性

**対策:**
- 定期的な再構築（10反復毎など）
- 残差ノルムモニタリング
- 必要に応じて自動再構築

**確率:** 低（H-Matrixは幾何学依存、磁化非依存）

#### リスク2: 自動閾値選択の失敗

**リスク:** N=200付近で最適手法選択を誤る

**対策:**
- ユーザーオーバーライド可能
- 詳細なログ出力
- 複数の閾値でベンチマーク

**確率:** 低（明確な性能差がある）

### プロジェクトリスク

#### リスク3: 実装複雑性

**リスク:** コードベース変更が予想以上に複雑

**対策:**
- 段階的実装（1機能ずつ）
- 包括的テストスイート
- コードレビュー

**確率:** 中（但し、設計文書により軽減）

#### リスク4: 後方互換性

**リスク:** 既存のユーザーコードが動かなくなる

**対策:**
- デフォルト動作は変更しない
- 新機能はオプトイン
- 移行ガイド提供

**確率:** 低（内部実装のみ変更）

---

## 成功基準

### Phase 2-A目標（必達）

✅ **パフォーマンス:**
- 磁化材料（N=343, 50反復）: 51秒 → < 3秒（17x高速化）
- 小問題（N=125）: 自動的にDense選択（5x高速化）

✅ **精度:**
- 誤差 < 1%（Dense vs H-Matrix）
- 収束性維持

✅ **安定性:**
- クラッシュなし
- メモリリークなし

### Phase 2-B目標（努力目標）

⭐ **パフォーマンス:**
- 磁化材料（N=343）: 51秒 → < 1.5秒（34x高速化）
- メモリ削減: 50%以上

⭐ **使いやすさ:**
- 自動最適化（ユーザー介入不要）
- 詳細なパフォーマンスレポート

---

## 参考文献

### 既存ドキュメント
1. `HMATRIX_PERFORMANCE_IMPROVEMENT.md` - パフォーマンス問題分析
2. `HMATRIX_PHASE1_IMPLEMENTATION.md` - Phase 1実装計画
3. `HMATRIX_TEST_RESULTS.md` - Phase 1テスト結果

### コードロケーション
1. `src/core/rad_hmatrix.h/cpp` - H-Matrix本体
2. `src/core/rad_intrc_hmat.h/cpp` - Relaxation統合
3. `src/core/rad_relaxation_methods.cpp` - AutoRelax()実装
4. `src/ext/HACApK_LH-Cimplm/` - HACApKライブラリ

### ベンチマーク
1. `examples/H-matrix/benchmark_solver.py` - ソルバーベンチマーク
2. `examples/solver_benchmarks/` - 各種ベンチマーク
3. `tests/benchmark_hmatrix.py` - 単体テスト

---

## 承認とレビュー

**提案者:** Claude Code
**日付:** 2025-11-12
**承認者:** TBD
**レビュアー:** TBD

**次のアクション:**
1. ⬜ レビューと承認
2. ⬜ 実装担当者アサイン
3. ⬜ Phase 2-A実装開始
4. ⬜ 進捗トラッキング

---

**最終更新:** 2025-11-12
**バージョン:** 1.0
**ステータス:** 🟢 Ready for Implementation
