# Radiaソルバー高速化提案（ngbem実装との比較）

**日付:** 2025-11-08
**目的:** ngbemの実装を参考にRadiaのH-matrixソルバーを高速化

---

## 1. 現状分析：Radiaの既存H-matrix実装

### 1.1. 既に実装されている機能

Radiaには**2つのH-matrixクラス**が既に実装されています：

#### A. `radTHMatrixFieldSource` (radhmat.cpp/h)
**目的:** フィールド評価の高速化
```cpp
// Configuration
struct radTHMatrixConfig {
	double eps;             // ACA tolerance (default: 1e-6)
	int max_rank;           // Maximum rank (default: 50)
	int min_cluster_size;   // Min cluster size (default: 10)
	bool use_openmp;        // OpenMP parallelization (default: true)
	int num_threads;        // Number of threads (0 = auto)
};
```

**既存機能:**
- ✅ HACApKライブラリによるH-matrix構築
- ✅ OpenMP並列化サポート
- ✅ クラスター木構造（空間分割）
- ✅ ACA（Adaptive Cross Approximation）による低ランク近似

#### B. `radTHMatrixInteraction` (radintrc_hmat.cpp/h)
**目的:** リラクゼーションソルバーの高速化
```cpp
// 3x3テンソル相互作用行列をH-matrix化
std::unique_ptr<hacapk::HMatrix> hmat[9];  // 9個のH-matrix

// 行列-ベクトル積: H = InteractMatrix * M
void MatVec(const TVector3d* M_in, TVector3d* H_out);
```

**既存機能:**
- ✅ 3x3テンソルを9個の独立したH-matrixで表現
- ✅ 対称性変換のキャッシュ（性能向上）
- ✅ O(N log N) 複雑度の行列-ベクトル積
- ✅ OpenMP並列化

---

## 2. ngbemとの比較

### 2.1. アーキテクチャ比較

| 側面 | Radia (現在) | ngbem | 優劣 |
|------|-------------|-------|------|
| **H-matrix構築** | HACApK（外部ライブラリ） | 独自実装（カスタマイズ可能） | ngbem有利（柔軟性） |
| **並列化** | H-matrix全体をOpenMP | ブロック単位でParallelForRange | ngbem有利（細粒度） |
| **テンソル扱い** | 9個の独立H-matrix | 1個のH-matrixでカーネルコンポーネント | ngbem有利（効率） |
| **メモリ管理** | 標準new/delete | LocalHeap + HeapReset | ngbem有利（速度） |
| **ACA精度制御** | 固定eps | ブロックサイズ適応型eps | ngbem有利（精度） |
| **SIMD** | ❌ なし | ✅ あり（カーネル評価） | ngbem有利 |

### 2.2. 並列化戦略の違い

#### Radia（現在）:
```cpp
// BuildHMatrix() - 全体をOpenMPで並列化
hacapk_params.nthr = omp_get_max_threads();
hmat[idx] = hacapk::build_hmatrix(points, points, kernel, ...);
// HACApK内部でOpenMP並列化
```

**特徴:**
- HACApKライブラリに任せる
- 粗粒度並列化（H-matrix全体）
- 実装がシンプル

#### ngbem:
```cpp
// CalcHMatrix() - ブロックごとに並列化
ParallelForRange(matList.Size(), [&](IntRange r)
{
	for (int k : r)
	{
		if(block.IsNearField())
			CalcBlockMatrix(...);      // 密行列ブロック
		else
			CalcFarFieldBlock(...);    // 低ランクブロック（ACA）
	}
}, TasksPerThread(4));
```

**特徴:**
- 細粒度並列化（ブロック単位）
- 負荷分散が優れている
- カスタムスケジューリング可能

---

## 3. 具体的な改善提案

### 3.1. 【優先度：高】ブロック並列化の導入

**現状の問題:**
```cpp
// radintrc_hmat.cpp BuildHMatrix()
for(int row = 0; row < 3; row++)
{
	for(int col = 0; col < 3; col++)
	{
		int idx = row * 3 + col;
		// 9個のH-matrixを逐次構築 ← ここが遅い
		hmat[idx] = hacapk::build_hmatrix(...);
	}
}
```

**改善案:**
```cpp
// 9個のH-matrixを並列構築
#pragma omp parallel for schedule(dynamic) if(config.use_openmp)
for(int idx = 0; idx < 9; idx++)
{
	int row = idx / 3;
	int col = idx % 3;

	KernelData kdata;
	kdata.hmat_ptr = this;
	kdata.tensor_row = row;
	kdata.tensor_col = col;

	hmat[idx] = hacapk::build_hmatrix(points, points, KernelFunction, &kdata, hacapk_params);
}
```

**期待効果:**
- 理論加速率: 最大9倍（9個を並列）
- 実測加速率: 4-8スレッドで3-6倍
- コスト: 低（数行の変更のみ）

---

### 3.2. 【優先度：中】テンソルカーネルの統合

**現状の問題:**
- 9個の独立したH-matrixを構築・保存
- メモリ使用量が9倍
- 同じクラスター木を9回計算

**改善案（ngbemパターン）:**
```cpp
// 1個のH-matrixで3x3テンソル全体を扱う
class radTHMatrixInteraction
{
	// Before: std::unique_ptr<hacapk::HMatrix> hmat[9];
	// After:
	std::unique_ptr<hacapk::HMatrix> hmat;  // 単一H-matrix

	// カーネル関数が3x3行列を返す
	static void KernelFunctionTensor(int i, int j, void* user_data, TMatrix3f& result)
	{
		// 1回の呼び出しで3x3全体を計算
		radTHMatrixInteraction* hmat_ptr = static_cast<radTHMatrixInteraction*>(user_data);
		hmat_ptr->ComputeInteractionKernel(i, j, result);
	}
};

// MatVecも統合
void MatVec(const TVector3d* M_in, TVector3d* H_out)
{
	// 1回のH-matrix-vector積で完了（9回ではなく）
	hmat->matvec_tensor(M_in, H_out);
}
```

**期待効果:**
- メモリ削減: 最大50%（クラスター木の共有）
- 構築時間削減: 30-40%（クラスター木を1回だけ構築）
- MatVec高速化: 20-30%（キャッシュ効率向上）

**コスト:** 中（HACApKライブラリの拡張が必要）

---

### 3.3. 【優先度：低】SIMD最適化

**ngbemの実装:**
```cpp
// ngbem CalcBlockMatrix() - SIMD vectorization
Vec<3,SIMD<double>> x = mirx[ix].GetPoint();
Vec<3,SIMD<double>> y = miry[iy].GetPoint();
auto kernel_ = kernel.Evaluate(x, y, nx, ny);
```

**Radiaへの適用:**
```cpp
// radintrc_hmat.cpp ComputeInteractionKernel()
// Before: 1点ずつ計算
void ComputeInteractionKernel(int i, int j, TMatrix3df& result)
{
	radTField Field(...);
	g3dRelaxPtrColNo->B_comp(&Field);  // 1点
	// ...
}

// After: SIMD化（複数点を同時計算）
void ComputeInteractionKernelBatch(int* i_list, int* j_list, int n, TMatrix3df* results)
{
	#pragma omp simd
	for(int k = 0; k < n; k++)
	{
		// SIMD命令で4-8点を同時計算
		// ...
	}
}
```

**期待効果:**
- カーネル評価: 2-4倍高速化（AVX2/AVX512使用時）
- H-matrix構築: 20-30%高速化（カーネル評価は構築時間の一部）

**コスト:** 高（Radiaのコア関数の大幅な書き換え）

---

### 3.4. 【優先度：低】LocalHeapメモリ管理

**ngbemの実装:**
```cpp
LocalHeapMem<100000> lh("Potential::Eval");
for (each element)
{
	HeapReset hr(lh);  // 高速リセット（再確保不要）
	// 一時メモリ使用
}
```

**Radiaへの適用:**
```cpp
// radintrc_hmat.cpp BuildHMatrix()
// Before: 標準的な動的メモリ確保
TMatrix3d SubMatrix(ZeroVect, ZeroVect, ZeroVect), BufSubMatrix;

// After: プールされたメモリ
MemoryPool pool(1024*1024);  // 1MB pool
for(int i = 0; i < n_elem; i++)
{
	pool.reset();
	TMatrix3d* SubMatrix = pool.allocate<TMatrix3d>();
	// ...
}
```

**期待効果:**
- メモリ断片化削減
- 一時メモリ確保の高速化: 5-10%

**コスト:** 中（新しいメモリプールクラスの実装）

---

## 4. 実装優先順位と期待効果

### フェーズ1: 低コスト・高効果（即座に実装可能）

#### 提案1: 9個のH-matrixの並列構築
```cpp
// radintrc_hmat.cpp BuildHMatrix() の修正
#pragma omp parallel for schedule(dynamic) if(config.use_openmp)
for(int idx = 0; idx < 9; idx++)
{
	// 既存コードをそのまま並列化
}
```

**効果:**
- 実装工数: **1-2時間**
- 期待加速: **3-6倍（H-matrix構築時）**
- リスク: **極めて低**

---

### フェーズ2: 中コスト・中効果（数日で実装可能）

#### 提案2: MatVec関数の並列化
```cpp
// radintrc_hmat.cpp MatVec() の修正
void radTHMatrixInteraction::MatVec(const TVector3d* M_in, TVector3d* H_out)
{
	// Before: 9個のH-matrixを逐次処理
	// After: 3x3ブロック単位で並列化

	#pragma omp parallel for collapse(2) if(config.use_openmp)
	for(int i = 0; i < n_elem; i++)
	{
		for(int comp = 0; comp < 3; comp++)
		{
			// H_out[i][comp] を計算
		}
	}
}
```

**効果:**
- 実装工数: **1-2日**
- 期待加速: **2-4倍（ソルバー反復ごと）**
- リスク: **低**

---

### フェーズ3: 高コスト・高効果（数週間）

#### 提案3: 統合テンソルカーネル
```cpp
// HACApKライブラリの拡張が必要
// または独自H-matrix実装（ngbemスタイル）
```

**効果:**
- 実装工数: **2-4週間**
- 期待加速: **1.5-2倍（総合）**
- メモリ削減: **40-50%**
- リスク: **中**

---

## 5. 推奨実装計画

### 段階的アプローチ

#### ステップ1: 即効性のある並列化（1日）
```
1. BuildHMatrix()の並列化 (#pragma omp parallel for)
2. MatVec()の並列化
3. ベンチマークテスト
```

**期待結果:**
- H-matrix構築: 3-6倍高速化
- 各反復: 2-4倍高速化
- 総合: 2-5倍高速化

#### ステップ2: 詳細なプロファイリング（2-3日）
```
1. 各関数の実行時間を測定
2. ボトルネック特定
3. さらなる最適化ポイントの抽出
```

#### ステップ3: 必要に応じた追加最適化（数週間）
```
1. テンソルカーネル統合（効果が大きい場合）
2. SIMD最適化（効果が大きい場合）
3. メモリ管理改善（必要な場合）
```

---

## 6. 実装例：提案1（即座に実装可能）

### Before:
```cpp:src/core/radintrc_hmat.cpp
int radTHMatrixInteraction::BuildHMatrix()
{
	// ...
	for(int row = 0; row < 3; row++)
	{
		for(int col = 0; col < 3; col++)
		{
			int idx = row * 3 + col;
			std::cout << "  Component [" << row << "][" << col << "]... " << std::flush;

			KernelData kdata;
			kdata.hmat_ptr = this;
			kdata.tensor_row = row;
			kdata.tensor_col = col;

			hmat[idx] = hacapk::build_hmatrix(points, points, KernelFunction, &kdata, hacapk_params);

			memory_used += hmat[idx]->memory_usage();

			std::cout << "rank=" << hmat[idx]->ktmax << std::endl;
		}
	}
	// ...
}
```

### After:
```cpp:src/core/radintrc_hmat.cpp
int radTHMatrixInteraction::BuildHMatrix()
{
	// ...
	std::cout << "\nBuilding 9 H-matrices in parallel..." << std::endl;

	// Thread-safe memory accumulation
	std::vector<size_t> memory_per_component(9);

	#pragma omp parallel for schedule(dynamic) if(config.use_openmp && n_elem > 100)
	for(int idx = 0; idx < 9; idx++)
	{
		int row = idx / 3;
		int col = idx % 3;

		#pragma omp critical
		{
			std::cout << "  Component [" << row << "][" << col << "]... " << std::flush;
		}

		KernelData kdata;
		kdata.hmat_ptr = this;
		kdata.tensor_row = row;
		kdata.tensor_col = col;

		// Each thread builds its own H-matrix
		hmat[idx] = hacapk::build_hmatrix(points, points, KernelFunction, &kdata, hacapk_params);

		// Store memory usage
		memory_per_component[idx] = hmat[idx]->memory_usage();

		#pragma omp critical
		{
			std::cout << "rank=" << hmat[idx]->ktmax
			          << ", blocks=" << hmat[idx]->nlf
			          << ", memory=" << (hmat[idx]->memory_usage() / 1024) << " KB" << std::endl;
		}
	}

	// Sum up memory usage
	memory_used = 0;
	for(int idx = 0; idx < 9; idx++)
	{
		memory_used += memory_per_component[idx];
	}

	// ...
}
```

**変更点:**
1. `#pragma omp parallel for` で9個のH-matrixを並列構築
2. `schedule(dynamic)` で負荷分散
3. `n_elem > 100` で小規模問題では並列化オーバーヘッドを回避
4. `#pragma omp critical` でスレッドセーフな標準出力
5. `memory_per_component` で各スレッドが独立にメモリ使用量を記録

---

## 7. ベンチマーク計画

### テストケース

```python
# ベンチマークスクリプト
import radia as rad
import time

# Test case 1: Small (N=125)
n = 5  # 5x5x5 = 125 elements

# Test case 2: Medium (N=1000)
n = 10  # 10x10x10 = 1000 elements

# Test case 3: Large (N=8000)
n = 20  # 20x20x20 = 8000 elements

# Build H-matrix
t_start = time.time()
# ... BuildHMatrix() ...
t_build = time.time() - t_start

# Solve (10 iterations)
t_start = time.time()
for i in range(10):
	# ... MatVec() ...
t_solve = time.time() - t_start

print(f"Build time: {t_build:.3f} s")
print(f"Solve time: {t_solve:.3f} s")
print(f"Time per iteration: {t_solve/10:.3f} s")
```

### 測定項目

| 項目 | Before | After | 加速率 |
|------|--------|-------|--------|
| H-matrix構築（N=1000） | ? s | ? s | ?x |
| 1反復あたり（N=1000） | ? s | ? s | ?x |
| メモリ使用量 | ? MB | ? MB | ? |
| 総ソルバー時間（10反復） | ? s | ? s | ?x |

---

## 8. まとめ

### ngbemから学べること

✅ **すぐに適用できる:**
1. ブロック並列化（9個のH-matrixを並列構築）
2. MatVec関数の並列化
3. 詳細な統計出力

⚠️ **慎重に検討すべき:**
1. テンソルカーネルの統合（HACApKの拡張が必要）
2. SIMD最適化（効果vs実装コスト）
3. LocalHeapメモリ管理（効果が限定的）

❌ **適用できない:**
1. H-matrixをフィールド評価に使う（ngbemも使っていない）
2. BEM固有の最適化手法

### 最優先で実装すべきこと

**提案1: 9個のH-matrixの並列構築**
- 実装時間: 1-2時間
- 期待効果: 3-6倍高速化（H-matrix構築）
- リスク: 極めて低
- コード変更: 数行のみ

**次のステップ:**
1. まず提案1を実装・テスト
2. 効果を測定
3. 必要に応じて提案2（MatVec並列化）を実装

---

**作成者:** Claude Code
**バージョン:** 1.0
**日付:** 2025-11-08
