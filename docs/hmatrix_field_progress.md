# H-Matrix Field Evaluation Implementation Progress

**Project:** Radia - H-matrix accelerated field evaluation
**Target:** `rad.Fld()` batch evaluation with 50-200x speedup
**Date Started:** 2025-11-08
**Status:** Phase 1 Complete (Infrastructure)

---

## 目標

- **目的**: `rad.Fld()`のバッチ評価をH-matrixで高速化
- **目標性能**: 50-200x speedup for N≥100, M≥100
- **計算量**: O(M×N) → O((M+N)log(M+N))
- **対象**: 非線形反復問題、大規模磁場計算

---

## 完了したタスク ✅

### Phase 1: Core Infrastructure (完了)

#### 1. 設計文書作成 ✅
- **ファイル**: `docs/hmatrix_field_design.md`
- **内容**:
  - API設計（グローバル設定 + バッチ関数）
  - キャッシング戦略
  - 実装段階の定義
  - 性能見積もり

#### 2. HMatrixFieldEvaluatorクラス実装 ✅
- **ヘッダー**: `src/core/radhmat.h` (line 211-388)
- **実装**: `src/core/radhmat_field.cpp` (新規作成)
- **主要機能**:
  ```cpp
  class radTHMatrixFieldEvaluator {
      // Constructor/Destructor
      radTHMatrixFieldEvaluator(const radTHMatrixConfig&);
      ~radTHMatrixFieldEvaluator();

      // Core methods
      int Build(radTGroup* source_group);           // H-matrix構築
      int EvaluateField(...)                         // フィールド評価
      bool IsValid(radTGroup*);                      // キャッシュ検証
      void Clear();                                  // メモリ解放

      // Private implementation
      int ExtractSourceGeometry(...);                // 幾何抽出
      int EvaluateFieldDirect(...);                  // 直接計算（フォールバック）
  };
  ```

#### 3. 幾何抽出機能 ✅
- **実装箇所**: `radhmat_field.cpp::ExtractSourceGeometry()`
- **機能**:
  - radTGroupから磁化要素を抽出
  - 中心座標、磁気モーメントを取得
  - radTg3dRelaxへのdynamic_cast使用
  - 単位変換処理（mm → m, mm³ → m³）
- **処理**:
  ```cpp
  for each element in source_group:
      cast to radTg3dRelax
      extract Magn (magnetization, A/m)
      extract CentrPoint (center, mm)
      extract Volume() (volume, mm³)
      compute magnetic moment = M * V (A·m²)
  ```

#### 4. 直接計算フォールバック ✅
- **実装箇所**: `radhmat_field.cpp::EvaluateFieldDirect()`
- **機能**:
  - H-matrix未構築時のフォールバック
  - 磁気双極子カーネル実装
  - OpenMP並列化対応（M>100）
- **カーネル**:
  ```
  H(r) = (1/4π) * Σ [3(m·r̂)r̂ - m] / |r|³
  ```

#### 5. キャッシュ無効化機構 ✅
- **実装箇所**:
  - `ComputeGeometryHash()`: ジオメトリのハッシュ計算
  - `IsValid()`: キャッシュの有効性チェック
- **ハッシュアルゴリズム**:
  - 要素数
  - 最初の10要素の座標（高速化のため）
  - 幾何変更時に自動無効化

#### 6. ビルドシステム統合 ✅
- **変更ファイル**: `CMakeLists.txt` (line 63)
- **追加内容**:
  ```cmake
  ${CORE_DIR}/radhmat_field.cpp
  ```
- **ビルド確認**: ✅ 成功（警告のみ、エラーなし）

---

## 現在の実装状況

### ✅ 完了
- [x] 設計文書
- [x] HMatrixFieldEvaluatorクラス構造
- [x] 幾何抽出
- [x] 直接計算フォールバック
- [x] キャッシュ無効化
- [x] ビルド成功

### ⏳ 部分実装（TODO付き）
- [ ] BuildTargetClusterTree() - HACApK cluster tree構築
- [ ] BuildFieldHMatrix() - H-matrix構築
- [ ] FieldKernel() - H-matrix用カーネル関数
- [ ] EvaluateField() - H-matrix加速評価

### 📋 未実装
- [ ] Python API bindings
- [ ] グローバル設定機構
- [ ] rad.FldBatch() 関数
- [ ] 性能ベンチマーク
- [ ] 精度検証

---

## 次のステップ

### Phase 2: HACApK Integration (次の実装)

#### 1. 観測点クラスタリング
**実装**: `BuildTargetClusterTree()`
```cpp
int BuildTargetClusterTree(const std::vector<TVector3d>& obs_points) {
    // 1. Convert TVector3d to hacapk::Point3D
    // 2. Build cluster tree using HACApK
    // 3. Store in target_cluster_tree
}
```

#### 2. フィールドカーネル実装
**実装**: `FieldKernel()`
```cpp
static double FieldKernel(int i, int j, void* kernel_data) {
    // Compute field at target i due to source j
    // Kernel: H(r) = (3(m·r̂)r̂ - m) / (4π|r|³)
    // Return scalar kernel value for HACApK
}
```

#### 3. H-matrix構築
**実装**: `BuildFieldHMatrix()`
```cpp
int BuildFieldHMatrix() {
    // 1. Define admissibility criterion
    // 2. Build H-matrix using HACApK
    // 3. Store in hmatrix_data
    // 4. Compute memory usage
}
```

#### 4. H-matrix評価
**実装**: `EvaluateField()`
```cpp
int EvaluateField(...) {
    // 1. BuildTargetClusterTree()
    // 2. BuildFieldHMatrix()
    // 3. H-matrix * moment_vector
    // 4. Extract field components
}
```

### Phase 3: Python API Integration

#### 1. グローバル設定
```python
rad.SetHMatrixFieldEval(True, eps=1e-6)
```

#### 2. バッチ関数
```python
H = rad.FldBatch(obj, 'h', obs_points, use_hmatrix=True)
```

#### 3. キャッシュ管理
```python
rad.ClearHMatrixCache()
```

### Phase 4: Testing & Optimization

#### 1. 性能ベンチマーク
- 様々な問題サイズ（N=100, 1000, 10000）
- 観測点数（M=100, 1000, 10000）
- 直接計算との速度比較

#### 2. 精度検証
- 直接計算との誤差評価
- 許容誤差ε=1e-6での精度確認

#### 3. メモリ使用量測定
- H-matrix メモリ使用量
- キャッシュ効率の評価

---

## 技術的詳細

### データ構造

#### Source Geometry
```cpp
std::vector<double> source_positions;  // [x1,y1,z1, x2,y2,z2, ...]
std::vector<double> source_moments;    // [mx1,my1,mz1, mx2,my2,mz2, ...]
```

#### Cache
```cpp
size_t geometry_hash;  // キャッシュ無効化用
int num_evaluations;    // 再利用回数
```

### 単位系

| 量 | Radia内部 | 変換 | SI単位 |
|----|----------|------|--------|
| 座標 | mm | ×1e-3 | m |
| 体積 | mm³ | ×1e-9 | m³ |
| 磁化 | A/m | - | A/m |
| 磁気モーメント | - | M×V×1e-9 | A·m² |
| 磁場H | A/m | - | A/m |

### 性能目標

| N | M | Direct | H-matrix | 目標Speedup |
|---|---|--------|----------|-------------|
| 100 | 100 | 10k ops | ~1.3k ops | 8x |
| 1,000 | 1,000 | 1M ops | ~20k ops | **50x** |
| 10,000 | 10,000 | 100M ops | ~200k ops | **500x** |

---

## ファイル構成

```
src/core/
├── radhmat.h                  # HMatrixFieldEvaluator宣言
├── radhmat.cpp                # HMatrixFieldSource実装（既存）
└── radhmat_field.cpp          # HMatrixFieldEvaluator実装（新規）

docs/
├── hmatrix_field_design.md    # 設計文書
└── hmatrix_field_progress.md  # 進捗文書（本ファイル）
```

---

## ビルド状況

### 最新ビルド
- **日時**: 2025-11-08
- **結果**: ✅ **成功**
- **警告**: Unicode encoding (C4819) - 無害
- **エラー**: なし
- **出力**: `build/Release/radia.cp312-win_amd64.pyd`

### ビルドコマンド
```bash
cmake --build build --config Release --target radia
```

---

## 依存関係

### 既存ライブラリ
- **HACApK**: H-matrix構築・演算
- **OpenMP**: 並列化

### Radiaクラス
- `radTGroup`: 磁性要素グループ
- `radTg3d`: 3D幾何基底クラス
- `radTg3dRelax`: 磁化可能要素
- `radTField`: フィールド構造体

---

## 今後の課題

### 実装課題
1. **HACApK integration**: cluster tree, H-matrix construction
2. **Kernel definition**: 3-component vector field kernel
3. **Memory management**: HACApK data structure lifecycle
4. **Error handling**: construction failure, numerical issues

### 設計課題
1. **Vector field H-matrix**: 3×3 block structure for vector fields?
2. **Adaptive precision**: ε調整の自動化
3. **Cache persistence**: ファイルへの保存・読み込み
4. **Multiple field types**: B, A, M への対応

### 性能課題
1. **Construction overhead**: 一回目の構築コスト
2. **Crossover point**: 直接計算とH-matrixの閾値最適化
3. **OpenMP scaling**: スレッド数最適化
4. **Memory footprint**: 大規模問題でのメモリ使用量

---

## まとめ

### 達成状況
✅ **Phase 1 完了**
- 基本インフラ構築完了
- ビルド成功
- フォールバック機能実装

### 次のマイルストーン
⏳ **Phase 2 開始**
- HACApK integration
- Kernel implementation
- H-matrix construction

### 最終目標への進捗
📊 **進捗: ~30%**
- Infrastructure: 100%
- HACApK integration: 0%
- Python API: 0%
- Testing: 0%

---

**更新日**: 2025-11-08
**次回更新予定**: Phase 2 完了時

## Python API Implementation ✅

### Phase 1.5: Python API Integration (完了)

#### 1. C API Implementation ✅
- **ファイル**: `src/lib/radentry_hmat.h`, `src/lib/radentry_hmat.cpp`
- **実装関数**:
  - `RadFldBatch()` - バッチフィールド評価
  - `RadSetHMatrixFieldEval()` - H-matrix 有効化/無効化
  - `RadClearHMatrixCache()` - キャッシュクリア
  - `RadGetHMatrixStats()` - 統計情報取得

#### 2. Python Wrappers ✅
- **ファイル**: `src/python/radpy_hmat.h`
- **実装内容**:
  - `radia_FldBatch()` - Python wrapper
  - `radia_SetHMatrixFieldEval()` - Python wrapper
  - `radia_ClearHMatrixCache()` - Python wrapper
  - `radia_GetHMatrixStats()` - Python wrapper

#### 3. Python API ✅
**利用可能な関数**:
```python
# Batch field evaluation
field = rad.FldBatch(obj, 'b', points, use_hmatrix=0)

# Global settings
rad.SetHMatrixFieldEval(enabled, eps)  # enabled: 0/1, eps: 1e-6

# Cache management
rad.ClearHMatrixCache()

# Statistics
stats = rad.GetHMatrixStats()  # [is_enabled, num_cached, total_memory_MB]
```

#### 4. Testing ✅
- **テストスクリプト**: `test_hmatrix_python_api.py`
- **テスト結果**: ✅ 全テスト成功
- **検証項目**:
  - rad.FldBatch() が rad.Fld() と完全一致（誤差 0）
  - 全API関数が正常動作
  - 直接計算モードで動作確認
  - バッチ評価が高速（3.97 ms vs 6.00 ms for 100 points）

---

## 更新されたビルド状況

### 最新ビルド（Python API版）
- **日時**: 2025-11-08
- **結果**: ✅ **成功**
- **警告**: Unicode encoding (C4819), LIBCMT conflict - 無害
- **エラー**: なし
- **出力**: `build/Release/radia.cp312-win_amd64.pyd`

### 追加ソースファイル
```
src/lib/
├── radentry_hmat.h        # C API 宣言
└── radentry_hmat.cpp      # C API 実装

src/python/
└── radpy_hmat.h           # Python wrappers

src/python/radpy.cpp       # メソッドエントリ追加
CMakeLists.txt             # ビルド設定更新
```

---

## 更新された実装状況

### ✅ 完了（Phase 1 + Python API）
- [x] 設計文書
- [x] HMatrixFieldEvaluatorクラス構造
- [x] 幾何抽出
- [x] 直接計算フォールバック
- [x] キャッシュ無効化
- [x] **C API実装**
- [x] **Python API bindings**
- [x] **Python API テスト**
- [x] **ビルド成功**

### ⏳ 部分実装（TODO付き）
- [ ] BuildTargetClusterTree() - HACApK cluster tree構築
- [ ] BuildFieldHMatrix() - H-matrix構築
- [ ] FieldKernel() - H-matrix用カーネル関数
- [ ] EvaluateField() - H-matrix加速評価（現在は直接計算のみ）

### 📋 未実装
- [ ] HACApK integration (Phase 2)
- [ ] H-matrix加速計算（現在は直接計算フォールバック）
- [ ] 性能ベンチマーク（大規模問題）
- [ ] 精度検証（H-matrix vs 直接計算）

---

## まとめ（更新）

### 達成状況
✅ **Phase 1 + Python API 完了**
- 基本インフラ構築完了
- Python API 実装・テスト完了
- ビルド成功
- 直接計算モード動作確認

### 次のマイルストーン
⏳ **Phase 2: HACApK Integration**
- HACApK cluster tree構築
- H-matrix field kernel実装
- H-matrix構築と評価
- 性能ベンチマーク

### 最終目標への進捗
📊 **進捗: ~50%**
- Infrastructure: 100%
- Python API: 100%
- HACApK integration: 0%
- Testing (basic): 100%
- Testing (performance): 0%

---

**更新日**: 2025-11-08 (Python API追加)
**次回更新予定**: Phase 2 (HACApK integration) 完了時
