# H-matrix フィールド評価問題の根本原因

**Date**: 2025-11-08
**Issue**: H-matrixがGridFunction.Set()で効果を発揮しない理由の詳細分析

---

## 🔍 コード解析結果

### フィールド評価の実装パス

1. **Python**: `rad.Fld(obj, 'b', [x, y, z])`
2. **C API**: `RadFld()` → `FieldArbitraryPointsArray()` (radentry.cpp:1011, radinter.cpp:1645)
3. **Application**: `radTApplication::ComputeField()` (radapl3.cpp:145)
4. **実装**:

```cpp
// radapl3.cpp:145-187
void radTApplication::ComputeField(int ElemKey, char* FieldChar, double** Points, long Np)
{
    // ...
    #pragma omp parallel for if(Np > 100)
    for(long i=0; i<Np; i++)
    {
        double *t = Points[i];
        TVector3d v; v.x = *(t++); v.y = *(t++); v.z = *t;

        FieldArray[i] = radTField(FieldKey, CompCriterium, v, ZeroVect, ...);
        g3dPtr->B_genComp(&(FieldArray[i]));  // ← ここで各点のフィールド計算
    }
    // ...
}
```

### 重要な発見

**各評価点で`B_genComp()`を個別に呼び出している**

- OpenMP並列化はあり（Np > 100の場合）
- しかし、各点の計算は完全に独立
- **H-matrixは一切使われていない**

---

## 🚫 H-matrixが使われない理由

### H-matrixが使われる場所

H-matrixは以下の場所**のみ**で使用されます：

#### 1. **相互作用行列の構築** (radintrc_hmat.cpp)
```cpp
// 要素i-j間の相互作用行列 K[i][j] を階層的に構築
void radTInteraction::CreateInteractionMatrix()
{
    if (UseHMatrix) {
        BuildHMatrixApproximation();  // H-matrix approximation
    } else {
        BuildDenseMatrix();  // Dense matrix
    }
}
```

#### 2. **線形ソルバー** (Gauss-Seidel, LU分解)
```cpp
// M_new = K^-1 * b を解く
void radTRelax::SolveLinearSystem()
{
    if (UseHMatrix) {
        HMatrixVectorProduct();  // O(N log N) 演算
    } else {
        DenseMatrixVectorProduct();  // O(N^2) 演算
    }
}
```

### H-matrixが使われない場所

#### ❌ **フィールド評価** (radapl3.cpp, B_genComp)

```cpp
// 各要素が個別にフィールド貢献を計算
void radTg3dGroup::B_genComp(radTField* FieldPtr)
{
    // グループの各要素を順次処理
    for (auto& element : GroupElements) {
        element->B_genComp(FieldPtr);  // 各要素の貢献を加算
    }
}

void radTRecMag::B_genComp(radTField* FieldPtr)
{
    // 直接計算：Biot-Savart積分
    // H-matrixは使用されない！
    ComputeFieldFromMagnetization();
}
```

**なぜH-matrixが使われないか**:
- フィールド評価は「既知の磁化分布M」から「点rでのフィールドB(r)」を計算
- 相互作用行列Kは要素間の相互作用を表すもの（M_i → M_j）
- フィールド評価は「要素i → 任意点r」で、Kとは異なる演算

---

## 📊 ベンチマーク結果の説明

### 期待していた挙動

```
H-matrixを使えば：
- 要素数N増加 → H-matrixで効率化
- 評価点M増加 → M個の評価でもH-matrix利用で高速化
```

### 実際の挙動

```
| N_elements | M_vertices | Time(ms) | H-matrix? |
|-----------|-----------|----------|-----------|
| 27        | 135       | ~15      | No effect |
| 343       | 5034      | ~7300    | No effect |

計算量: O(M × N) （H-matrixの有無に関わらず）
```

**理由**: `B_genComp()`が各要素を直接ループで処理しているため

---

## 🔧 なぜOpenMP並列化も限定的か

### 現在の並列化

```cpp
#pragma omp parallel for if(Np > 100)
for(long i=0; i<Np; i++)
{
    g3dPtr->B_genComp(&(FieldArray[i]));
}
```

**問題点**:
1. **Python overhead**: 各点でPython ↔ C++の呼び出し
   - `rad.Fld()` が個別に呼ばれる（rad_ngsolve.cpp:220）
   - 並列化はC++側でのみ有効（複数点を一度に渡した場合のみ）

2. **NGSolveとの統合**:
   ```cpp
   // rad_ngsolve.cpp: RadiaFieldCF::Evaluate()
   py::object field_result = rad.attr("Fld")(radia_obj, field_type, coords);
   ```
   - 各NGSolveメッシュ頂点で個別にPython呼び出し
   - 100点バッチでも無い限り、OpenMP並列化されない

---

## 💡 解決策の検討

### Option A: H-matrixをフィールド評価に拡張

**アイデア**: 既存のH-matrix構造を流用

```cpp
// 擬似コード
B(r) = Σ G(r, r_i) * M_i
     ≈ H-matrix的に階層評価
```

**課題**:
- H-matrixは要素間の相互作用（K[i][j]）を表現
- フィールド評価は「要素 → 任意点」で、行列形式ではない
- 大規模な実装変更が必要

**効果**: O(M × N) → O(M × log(N)) 程度の改善可能

### Option B: Fast Multipole Method (FMM)

**原理**: 階層的なフィールド評価
- 遠方の要素群をまとめて多重極展開
- 近傍の要素のみ直接計算

**効果**: O(M × N) → O(M + N) に削減

**課題**: 大規模な新規実装

### Option C: ベクトル化API（最小の変更）

**現状の問題**:
```python
# NGSolve側で1点ずつ評価
for vertex in mesh.vertices:
    B = rad.Fld(obj, 'b', vertex.coords)  # Python call overhead
```

**改善案**:
```python
# 全頂点を一度に評価
vertices = mesh.get_all_vertex_coords()  # M×3 array
B_array = rad.FldVec(obj, 'b', vertices)  # 1回のPython call
```

**効果**:
- Python overhead削減: M回 → 1回
- OpenMP並列化が確実に有効
- 既存のComputeField()が既に複数点対応

**実装規模**: 小（rad_ngsolve.cppの修正のみ）

### Option D: rad_ngsolveでのキャッシング

**アイデア**: 粗いグリッドでフィールド計算 → 補間

```python
# 粗いグリッドで計算（例: h=0.1m）
coarse_grid = generate_grid(h=0.1)
B_coarse = [rad.Fld(obj, 'b', p) for p in coarse_grid]

# 細かいメッシュでは補間
for vertex in fine_mesh:
    B[vertex] = interpolate(B_coarse, vertex.coords)
```

**効果**: M回の評価 → M_coarse回の評価
**精度**: トレードオフあり

---

## 📌 結論

### 根本原因

1. **H-matrixはソルバー用**: 相互作用行列Kの操作のみ
2. **フィールド評価は別経路**: `B_genComp()`で直接計算
3. **NGSolve統合のオーバーヘッド**: 各点で個別にPython呼び出し

### 現実的な対応

**短期**（すぐ実装可能）:
- ✅ **Option C**: ベクトル化API `rad.FldVec()` の実装
  - 最小の変更でPython overhead削減
  - OpenMP並列化が確実に有効

**中期**（数週間の開発）:
- **Option D**: キャッシング/補間による高速化
  - 精度とのバランス調整

**長期**（大規模プロジェクト）:
- **Option B**: FMM実装
  - 根本的なO(N)化
  - 他の用途にも有益

### ユーザーへの推奨

**現時点**:
- H-matrixパラメータは`use_hmatrix`として利用可能
- ただし、GridFunction.Set()には効果なし
- ソルバー（rad.Solve）には効果あり

**ワークアラウンド**:
1. メッシュサイズを抑える（M を減らす）
2. Radia要素数を減らす（N を減らす）
3. 対称性を利用してフィールド評価領域を限定

---

## 📚 参考

- **H-matrix文献**: Hackbusch, Wolfgang. "A sparse matrix arithmetic based on H-matrices."
- **FMM文献**: Greengard & Rokhlin, "A fast algorithm for particle simulations"
- **実装例**: bempp-cl (BEM with FMM)

---

**次のアクション**: Option C (ベクトル化API) の実装を検討
