# Refactoring Plan for Legacy Files

## 概要

27個の古い命名規則ファイル（`rad_` で始まらず、`_` がないファイル）をリファクタリングします。

## 推奨アプローチ

### ⚠️ 重要な注意事項

大規模なリファクタリングは以下の理由により**段階的に**進めることを強く推奨します：

1. **破壊的変更のリスク**: 一度に多数のファイルを変更すると、予期しないバグが発生しやすい
2. **ビルドエラー**: CMakeLists.txt、include文、依存関係の更新が必要
3. **テストの必要性**: 各段階でビルドとテストを実行して安全性を確認
4. **レビュー可能性**: 小さな変更の方がレビューしやすい

---

## 提案: 3つのオプション

### オプション A: 段階的リファクタリング（推奨）

**利点**: 安全、テスト可能、段階的な進行
**欠点**: 時間がかかる

#### Phase 1: ヘッダーのみファイル（12ファイル）
低リスク、影響範囲が明確

1. radappl.h → rad_application.h
2. radauxst.h → rad_auxiliary_structures.h
3. radcnvrg.h → rad_convergence.h
4. radexcep.h → rad_exception.h
5. radg.h → rad_geometry_base.h
6. radg3da1.h → rad_geometry_3d_aux.h
7. radhandl.h → rad_handle.h
8. radmtra1.h → rad_material_aux.h
9. radopnam.h → rad_operation_names.h
10. radstlon.h → rad_string_long.h
11. radtrans.h → rad_transform_def.h
12. radyield.h → rad_yield.h

**作業手順**:
```bash
# 1. ファイル名変更
git mv radappl.h rad_application.h

# 2. 全ファイルでinclude文を更新
find . -name "*.cpp" -o -name "*.h" | xargs sed -i 's/#include "radappl.h"/#include "rad_application.h"/g'

# 3. CMakeLists.txt更新
# 4. ビルドとテスト
cmake --build build --config Release
pytest tests/

# 5. 次のファイルへ
```

#### Phase 2: コア実装ファイル（8ファイル）
中〜高リスク、単体試験追加が必須

1. **radcast.cpp/h** → rad_type_cast.cpp/h
2. **radg3d.cpp/h** → rad_geometry_3d.cpp/h
3. **radgroup.cpp/h** → rad_group.cpp/h
4. **radinter.cpp** → rad_interaction.cpp

**作業手順**:
```bash
# 各ファイルごとに:
# 1. 単体試験を追加（リファクタリング前）
# 2. ファイル名変更
# 3. include文更新
# 4. ビルドとテスト
# 5. 次のファイルへ
```

#### Phase 3: Material/Transform/Serialization（7ファイル）

1. radmater.cpp/h → rad_material_def.cpp/h
2. radmaterial.cpp → rad_material_impl.cpp
3. radtransform.cpp → rad_transform_impl.cpp
4. radsend.cpp/h → rad_serialization.cpp/h
5. radhmat.cpp/h → rad_hmatrix.cpp/h

---

### オプション B: 一括リファクタリング（非推奨）

**利点**: 速い
**欠点**: 高リスク、デバッグが困難

⚠️ **推奨しません** - バグの特定が困難になります

---

### オプション C: シンボリックリンク/エイリアス方式

旧ファイル名を残したまま、新ファイル名へのエイリアスを作成

**利点**: 互換性維持、段階的移行
**欠点**: 二重管理の手間

例:
```cpp
// radappl.h (旧ファイル、deprecated)
#pragma message("Warning: radappl.h is deprecated, use rad_application.h")
#include "rad_application.h"
```

---

## 単体試験の追加

### 優先度: High

以下のファイルは単体試験が不足している可能性が高い：

1. **radcast.cpp** - Type casting
   ```python
   # tests/test_type_cast.py
   def test_cast_to_g3d():
       # Test casting to radTg3d

   def test_cast_to_group():
       # Test casting to radTGroup
   ```

2. **radgroup.cpp** - Grouping operations
   ```python
   # tests/test_group.py
   def test_create_group():
       # Test group creation

   def test_group_operations():
       # Test add, remove, etc.
   ```

3. **radmaterial.cpp** - Material system
   ```python
   # tests/test_material_system.py
   def test_linear_material():
       # Test rad.MatLin

   def test_nonlinear_material():
       # Test rad.MatSatIsoTab
   ```

4. **radtransform.cpp** - Transformations
   ```python
   # tests/test_transformations.py
   def test_translation():
       # Test rad.TrfTrsl

   def test_rotation():
       # Test rad.TrfRot
   ```

---

## 推奨実行順序

### Step 1: 単体試験追加（1-2日）
- [ ] test_type_cast.py
- [ ] test_group.py
- [ ] test_material_system.py
- [ ] test_transformations.py
- [ ] test_serialization.py

### Step 2: Phase 1 リファクタリング（1日）
ヘッダーのみファイル12個を一括リネーム

### Step 3: Phase 2 リファクタリング（2-3日）
コア実装ファイルを1つずつ慎重にリネーム

### Step 4: Phase 3 リファクタリング（2日）
残りのファイルをリネーム

### Step 5: 最終検証（1日）
- 全テスト実行
- ビルド確認
- ドキュメント更新

**合計見積もり**: 7-9日

---

## 質問

リファクタリングを進める前に、以下を確認させてください：

1. **どのオプションを選びますか？**
   - A: 段階的リファクタリング（推奨）
   - B: 一括リファクタリング（非推奨）
   - C: シンボリックリンク方式

2. **優先度は？**
   - 単体試験追加を優先？
   - ファイル名変更を優先？
   - 両方同時進行？

3. **スコープ**
   - 全27ファイルを一度に？
   - まずは一部（例：ヘッダーのみ）から？

---

**推奨**: オプションA + 単体試験優先 + ヘッダーファイルから開始

ご希望を教えていただければ、具体的な作業を開始します！
