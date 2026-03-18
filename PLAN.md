# PLAN: OpenMP → NGSolve TaskManager 移行

## 背景

### 問題
NGSolve の `import` 時に `mkl_set_num_threads(1)` が呼ばれ、Radia の MKL ベース LU ソルバー (`dgesv_`) がシングルスレッドで動作する。i7-9700K (8コア) で CPU 利用率が約 25% に低下。

### 原因チェーン
```
import radia
  -> esim_cell_problem.py
    -> from ngsolve import *
      -> ngslib C++ init
        -> mkl_set_num_threads(1)   <- ここで MKL が 1 スレッドに制限される
```

Radia は `mkl_service.h` を include しているが、`mkl_set_num_threads()` を一度も呼ばない。

### 方針
本 Fork は NGSolve を import することが前提。
NGSolve は TaskManager ベースの並列化を採用し、OpenMP を意図的に無効化している。
Radia も NGSolve と同じ TaskManager に移行することで、スレッディング競合を根本的に解決する。

**OpenMP フォールバックは不要** — NGSolve が常に利用可能。

---

## ブランチ戦略

- **作業ブランチ**: `main` （直接作業）
- `peec-dev` にも同じ変更を反映する（マージまたはチェリーピック）

---

## 対象ファイル一覧

### 1. C++ ソースコード

#### Radia Core (`src/core/`) — 50 箇所

| ファイル | #pragma omp 数 | 用途 | 優先度 |
|----------|---------------|------|--------|
| `rad_relaxation_methods.cpp` | 23 | BiCGSTAB, GMRES, Jacobi, dgesv_, 残差 | **最高** |
| `rad_interaction.cpp` | 6 | 相互作用行列構築 | **最高** |
| `rad_field_unified.cpp` | 4 | 磁場計算 (B, H, A) | 高 |
| `rad_peec_matrices.cpp` | 5 | PEEC 行列 (L, P, R) | 中 |
| `rad_transform_impl.cpp` | 3 | 対称性変換 | 中 |
| `rad_mmm_matrices.cpp` | 2 | MMM 行列構築 | 高 |
| `rad_hacapk.cpp` | 2 | H行列ブロック充填 | **最高** |
| `rad_exafmm.cpp` | 2 | ExaFMM direct 計算 | 中 |
| `rad_material_impl.cpp` | 1+1 | 材料場計算 + omp_get_max_threads | 低 |
| `rad_poly_analytical.cpp` | 1 | 解析的多面体場 | 低 |
| `rad_point_classify.cpp` | 1 | 内外判定 | 低 |
| `rad_cln.cpp` | 1 | CLN 周波数スイープ | 低 |

#### HACApK (`src/ext/HACApK/`) — 14 箇所

| ファイル | #pragma omp 数 | 用途 | 優先度 |
|----------|---------------|------|--------|
| `cHACApK_base.c` | 8 | H行列構築 (ACA+) | **最高** |
| `cHACApK_cpp_impl.c` | 6 + omp_* 3 | H行列 matvec | **最高** |

#### ExaFMM-t (`external/exafmm-t/`) — ~67 箇所

| ファイル | #pragma omp 数 | 用途 | 優先度 |
|----------|---------------|------|--------|
| `include/fmm.h` | ~25 | P2M, M2L, L2P 等 FMM 演算子 | 高 |
| `include/fmm_scale_invariant.h` | ~20 | FMM scale-invariant 版 | 高 |
| `include/fmm_base.h` | ~8 | P2P, M2M, L2L 基底 | 高 |
| `include/test.h` | ~6 | テスト用 | 低 |
| `python/exafmm.cpp` | ~9 | Python バインディング | 中 |
| 他ヘッダ | ~5 | geometry, build_list, timer | 低 |

### 2. ビルドシステム

| ファイル | OpenMP 参照数 | 内容 |
|----------|-------------|------|
| `CMakeLists.txt` | ~60 | 全ターゲットの OpenMP リンク設定 |
| `CMakeLists_ngsolve_msvc.txt` | ~10 | NGSolve MSVC ビルド用 |
| `Build.ps1` | 4 | `-NoOpenMP` フラグ |
| `external/exafmm-t/CMakeLists.txt` | 3 | ExaFMM OpenMP 設定 |

### 3. ドキュメント

| ファイル | OpenMP 参照数 | 更新内容 |
|----------|-------------|----------|
| `CLAUDE.md` | 5 | ポリシー: "OpenMP: Intel OpenMP Only" -> "TaskManager" |
| `BUILD.md` | ~15 | ビルド手順、依存関係、トラブルシューティング |
| `CHANGELOG.md` | ~10 | 歴史的記録（移行エントリ追加のみ、過去は変更しない） |
| `README.md` | 要確認 | 概要の並列化説明 |
| `pyproject.toml` | 1 | description: "OpenMP Parallelization" を変更 |
| `SECURITY.md` | 1 | libiomp5md.dll の記述 |
| `docs/api/API_REFERENCE.md` | 1 | FldVTS の OpenMP 記述 |
| `tests/README.md` | ~10 | テスト説明、OpenMP ベンチマーク |
| `examples/cube_uniform_field/README.md` | 1 | OMP_NUM_THREADS 設定 |
| `examples/fmm_field_evaluation/README.md` | 1 | OpenMP 並列化説明 |
| `examples/solver_benchmarks/README.md` | ~5 | OpenMP 並列構築説明 |
| `examples/solver_benchmarks/BENCHMARK_RESULTS.md` | 3 | OpenMP 有効記述 |
| `docs/peec/PEEC_PANEL_IMPLEMENTATION.md` | 1 | OpenMP 最適化記述 |

### 4. Python スクリプト（OMP_NUM_THREADS 設定）

| ファイル | 内容 |
|----------|------|
| `examples/cube_uniform_field/benchmark_common.py` | `os.environ['OMP_NUM_THREADS'] = '8'` |
| `examples/cube_uniform_field/benchmark_cube_block_jacobi.py` | 同上 |
| `examples/cube_uniform_field/tetrahedron/benchmark_tetra.py` | 同上 |
| `examples/cube_uniform_field/hexahedron/benchmark_hex.py` | 同上 |
| `tests/test_parallel_performance.py` | OMP_NUM_THREADS 動的設定 |
| `tests/benchmarks/benchmark_threads.py` | 同上 |
| `tests/benchmarks/benchmark_openmp.py` | OpenMP スケーリングテスト |
| `tests/benchmarks/benchmark_heavy.py` | 同上 |
| `tests/benchmarks/benchmark_correct.py` | 同上 |

### 5. テストファイル

| ファイル | 内容 |
|----------|------|
| `tests/test_fldbatch_openmp.py` | FldBatch OpenMP テスト -> TaskManager テストに変更 |
| `tests/test_parallel_performance.py` | 並列性能テスト -> TaskManager 版に変更 |
| `tests/benchmarks/benchmark_openmp.py` | OpenMP スケーリング -> TaskManager スケーリングに変更 |

---

## 変換パターン

### パターン 1: parallel for (42 箇所) -> `ParallelFor`

```cpp
// Before (OpenMP)
#pragma omp parallel for if(n > 100)
for (int i = 0; i < n; i++) {
    result[i] = compute(i);
}

// After (TaskManager)
ParallelFor(Range(n), [&](size_t i) {
    result[i] = compute(i);
});
```

### パターン 2: parallel for schedule(dynamic) (12 箇所) -> `ParallelFor`

```cpp
// Before
#pragma omp parallel for schedule(dynamic) if(n > 20)
for (int i = 0; i < n; i++) { ... }

// After -- TaskManager は動的スケジューリングがデフォルト
ParallelFor(Range(n), [&](size_t i) { ... });
```

### パターン 3: reduction (3 箇所) -> MKL BLAS or AtomicAdd

```cpp
// Before
double sum = 0;
#pragma omp parallel for reduction(+:sum) if(n > 100)
for (int i = 0; i < n; i++) {
    sum += x[i] * y[i];
}

// After -- MKL cblas_ddot（BLAS最適化、MKL内部で並列化）
double sum = cblas_ddot(n, x, 1, y, 1);

// または AtomicAdd (NGSolve提供)
double sum = 0;
ParallelFor(Range(n), [&](size_t i) {
    AtomicAdd(sum, x[i] * y[i]);
});
```

### パターン 4: critical (7 箇所) -> `std::mutex` or `AtomicAdd`

```cpp
// Before
#pragma omp critical
{ BufMaxModM = std::max(BufMaxModM, localMax); }

// After
static std::mutex mtx;
{
    std::lock_guard<std::mutex> lock(mtx);
    BufMaxModM = std::max(BufMaxModM, localMax);
}
```

### パターン 5: parallel + omp for (HACApK matvec) -> C ラッパー経由 ParallelFor

```cpp
// Before (cHACApK_cpp_impl.c)
int nthreads = omp_get_max_threads();
double *ytmp = calloc(nthreads * nd, sizeof(double));
#pragma omp parallel
{
    int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic, 32)
    for (int ip = 0; ip < nlf; ip++) { ... ytmp[tid*nd + j] += ...; }
}

// After -- C ラッパー経由で TaskManager を使用
int nthreads = hacapk_get_num_threads();
double *ytmp = calloc(nthreads * nd, sizeof(double));
// hacapk_parallel_for() 内部で ParallelFor を呼ぶ
```

### パターン 6: omp task (ExaFMM, M2M/L2L 木走査) -> TaskManager task

```cpp
// Before (MSVC では EXAFMM_NO_OMP_TASKS で既に無効化)
#pragma omp task untied
M2M(parent, child);
#pragma omp taskwait

// After -- NGSolve TaskManager の task API
// AddTask() + WaitForTasks() に置換
```

---

## MKL スレッディング

### 現状の問題
- `rad_relaxation_methods.cpp:48` で `#include "mkl_service.h"` はあるが呼び出しなし
- NGSolve が `mkl_set_num_threads(1)` を呼んだ後、Radia が復元しない

### 解決策
TaskManager に移行しても、`dgesv_()` 等の MKL LAPACK は内部で独自スレッディングを使用する。
NGSolve 方式に合わせて、MKL 呼び出し前後でスレッド数を制御する。

```cpp
// dgesv_ 呼び出し前
int saved_mkl_threads = mkl_get_max_threads();
mkl_set_num_threads(TaskManager::GetNumThreads());

dgesv_(&totalDOF, &nrhs, ...);

// 呼び出し後に復元
mkl_set_num_threads(saved_mkl_threads);
```

---

## HACApK (C 言語) の対応

HACApK は C で書かれており、C++ の TaskManager API を直接使えない。

### 方法: C ラッパー関数

```c
// hacapk_taskmanager.h -- C から呼べるラッパー
#ifdef __cplusplus
extern "C" {
#endif

void hacapk_parallel_for(int start, int end, void (*func)(int, void*), void* data);
int hacapk_get_num_threads(void);
int hacapk_get_thread_id(void);

#ifdef __cplusplus
}
#endif
```

```cpp
// hacapk_taskmanager.cpp -- C++ 実装
#include <ngsolve/taskmanager.hpp>

extern "C" void hacapk_parallel_for(int start, int end,
                                     void (*func)(int, void*), void* data) {
    ParallelFor(Range(start, end), [&](size_t i) {
        func((int)i, data);
    });
}

extern "C" int hacapk_get_num_threads(void) {
    return TaskManager::GetNumThreads();
}

extern "C" int hacapk_get_thread_id(void) {
    return TaskManager::GetThreadId();
}
```

HACApK の `cHACApK_base.c`, `cHACApK_cpp_impl.c` 内の `#pragma omp` を
`hacapk_parallel_for()` 呼び出しに置換する。

---

## ExaFMM-t の対応

### 課題
- 外部ライブラリ（Fork 管理）
- `omp task` を使用している（M2M/L2L 木走査）
- MSVC では `EXAFMM_NO_OMP_TASKS` で既に task 無効化済み
- GridFunction との整合性のため対応必須

### 方針
1. `#pragma omp parallel for` -> TaskManager `ParallelFor` に置換
2. `omp task` -> TaskManager の task API に置換
3. `omp_get_max_threads()`, `omp_set_num_threads()` -> TaskManager API に置換
4. `EXAFMM_NO_OMP_TASKS` マクロ -> 削除（TaskManager で統一）

---

## ビルドシステム変更

### CMakeLists.txt

```cmake
# Before: OpenMP を検索・リンク
option(RADIA_ENABLE_OPENMP "Enable OpenMP parallelization" ON)
find_package(OpenMP)
target_link_libraries(radia_ngsolve PRIVATE ${INTEL_OPENMP_LIB})

# After: NGSolve TaskManager を使用（OpenMP 検索は不要）
# NGSolve は既に find_package(NGSolve) で見つかっている
# TaskManager は ngstd ライブラリに含まれる
target_link_libraries(radia_ngsolve PRIVATE ngsolve::ngstd)
target_compile_definitions(radia_ngsolve PRIVATE RADIA_USE_TASKMANAGER)
```

変更対象ターゲット:
- `radia_ngsolve` (~20行)
- `peec_matrices` (~15行)
- `cln_core` (~15行)
- `mmm_core` (~15行)
- `_radia_pybind` (~15行)

### Build.ps1

- `-NoOpenMP` フラグ -> 削除または `-NoTaskManager` に変更
- `RADIA_ENABLE_OPENMP` -> `RADIA_USE_TASKMANAGER`

### libiomp5md.dll

OpenMP を使わなくなるため:
- `libiomp5md.dll` のコピー処理 -> 削除
- MKL は引き続き `mkl_rt.2.dll` 経由でリンク
- MKL 内部の OpenMP ランタイムは MKL 自身が管理

---

## ドキュメント更新詳細

### CLAUDE.md

| セクション | 変更内容 |
|-----------|----------|
| "OpenMP: Intel OpenMP Only" (L285-287) | -> "Parallelization: NGSolve TaskManager" |
| "Required MKL DLLs" (L283) | libiomp5md.dll の記述を更新 |
| ComputeFieldBatch "OpenMP parallelized" (L154) | -> "TaskManager parallelized" |
| "OpenMP parallelized batch" (L171) | -> "TaskManager parallelized batch" |

### BUILD.md

| セクション | 変更内容 |
|-----------|----------|
| `-NoOpenMP` オプション (L59) | 削除 |
| "Intel OpenMP" 記述 (L75, 93, 99, 528, 540, 546, 553) | -> TaskManager に更新 |
| "OpenMP issues" トラブルシューティング (L324-335) | -> TaskManager に更新 |
| "Intel OpenMP (`libiomp5md.dll`)" (L604) | 更新 |

### pyproject.toml

```toml
# Before
description = "Radia 3D Magnetostatics with NGSolve Integration and OpenMP Parallelization"
# After
description = "Radia 3D Magnetostatics with NGSolve Integration and TaskManager Parallelization"
```

### CHANGELOG.md

過去のエントリは変更しない。新しい移行エントリを追加:
```markdown
### 20XX-XX-XX
- **Parallelization: OpenMP -> NGSolve TaskManager**
  - Migrated all parallel loops from OpenMP to NGSolve TaskManager
  - Resolved MKL threading conflict (NGSolve mkl_set_num_threads(1))
  - HACApK: C wrapper for TaskManager integration
  - ExaFMM-t: Full TaskManager migration for GridFunction compatibility
  - Removed OpenMP dependency (libiomp5md.dll no longer required by Radia)
```

### テスト README・ベンチマーク

| ファイル | 変更内容 |
|----------|----------|
| `tests/README.md` | "OpenMP" -> "TaskManager"、テスト説明更新 |
| `tests/test_fldbatch_openmp.py` | リネーム: `test_fldbatch_parallel.py` |
| `tests/benchmarks/benchmark_openmp.py` | リネーム: `benchmark_parallel.py` |
| `examples/solver_benchmarks/README.md` | "OpenMP parallelization" -> "TaskManager" |
| `examples/cube_uniform_field/README.md` | `OMP_NUM_THREADS` 設定削除 |

### Python ベンチマークスクリプト

`os.environ['OMP_NUM_THREADS']` 設定 -> 削除（TaskManager が自動管理）

対象 (9 ファイル):
- `examples/cube_uniform_field/benchmark_common.py`
- `examples/cube_uniform_field/benchmark_cube_block_jacobi.py`
- `examples/cube_uniform_field/tetrahedron/benchmark_tetra.py`
- `examples/cube_uniform_field/hexahedron/benchmark_hex.py`
- `tests/test_parallel_performance.py`
- `tests/benchmarks/benchmark_threads.py`
- `tests/benchmarks/benchmark_openmp.py`
- `tests/benchmarks/benchmark_heavy.py`
- `tests/benchmarks/benchmark_correct.py`

---

## 実装フェーズ

### Phase 1: 基盤 + MKL 修正 (即効性あり)
1. `radia_parallel.h` 統一ヘッダ作成（TaskManager 直接使用、フォールバックなし）
2. `CMakeLists.txt` から OpenMP 関連設定を削除、TaskManager リンクに変更
3. `rad_relaxation_methods.cpp`: `dgesv_()` 前後で `mkl_set_num_threads()` 追加
4. ビルド確認・ベンチマーク比較

### Phase 2: Radia Core ソルバー (最高優先度)
5. `rad_relaxation_methods.cpp`: 23 箇所 -> `ParallelFor`
6. `rad_interaction.cpp`: 6 箇所
7. `rad_hacapk.cpp`: 2 箇所
8. ベンチマーク: LU, BiCGSTAB, HACApK の性能確認

### Phase 3: HACApK (C ラッパー)
9. `hacapk_taskmanager.h/cpp` C ラッパー作成
10. `cHACApK_base.c`: 8 箇所 -> `hacapk_parallel_for()`
11. `cHACApK_cpp_impl.c`: 6 箇所 -> `hacapk_parallel_for()`
12. H行列構築・matvec のベンチマーク

### Phase 4: 磁場計算・その他 Core
13. `rad_field_unified.cpp`: 4 箇所
14. `rad_mmm_matrices.cpp`: 2 箇所
15. `rad_transform_impl.cpp`: 3 箇所
16. 残り 6 ファイル (各 1-5 箇所)

### Phase 5: ExaFMM-t
17. `include/fmm.h`: ~25 箇所
18. `include/fmm_base.h`: ~8 箇所
19. `include/fmm_scale_invariant.h`: ~20 箇所
20. `python/exafmm.cpp`: ~9 箇所
21. その他ヘッダ: ~5 箇所
22. GridFunction との整合性テスト

### Phase 6: ビルドシステム・スクリプト
23. `CMakeLists.txt`: OpenMP 設定削除 (~60行)
24. `CMakeLists_ngsolve_msvc.txt`: 同上 (~10行)
25. `Build.ps1`: `-NoOpenMP` 削除
26. `external/exafmm-t/CMakeLists.txt`: OpenMP -> TaskManager
27. Python ベンチマーク: `OMP_NUM_THREADS` 設定削除 (9 ファイル)

### Phase 7: ドキュメント全更新
28. `CLAUDE.md`: OpenMP ポリシー -> TaskManager ポリシー
29. `BUILD.md`: ビルド手順更新
30. `CHANGELOG.md`: 移行エントリ追加
31. `pyproject.toml`: description 更新
32. `SECURITY.md`: libiomp5md.dll 記述更新
33. `docs/api/API_REFERENCE.md`: OpenMP 記述更新
34. `tests/README.md`: テスト説明更新
35. テストファイルリネーム (`test_fldbatch_openmp.py` 等)
36. `examples/*/README.md`: OpenMP 記述更新 (4 ファイル)
37. `examples/solver_benchmarks/BENCHMARK_RESULTS.md`: OpenMP 記述更新

### Phase 8: テスト・検証
38. 全ソルバーの回帰テスト（数値精度一致確認）
39. スレッドスケーリング測定（1, 2, 4, 8 スレッド）
40. cube_uniform_field + c_type_electromagnet ベンチマーク再実行
41. peec-dev ブランチへの反映・マージ

---

## リスク

| リスク | 影響 | 対策 |
|--------|------|------|
| TaskManager API の学習コスト | 中 | NGSolve ソースから実例を参照 |
| HACApK (C) との結合 | 高 | C ラッパーで分離 |
| ExaFMM-t の変更量 (~67 箇所) | 高 | Phase 5 で独立実施、段階的テスト |
| reduction パターンの変換 | 中 | MKL BLAS (cblas_ddot) 活用 |
| omp task -> TaskManager task | 中 | MSVC では既に無効化済み |
| 性能劣化 | 低 | Phase 毎にベンチマーク比較 |
| peec-dev とのコンフリクト | 中 | Phase 8 で対応 |

---

## 備考

- NGSolve の TaskManager は `with TaskManager(): ...` (Python) で起動する
- C++ では `RunWithTaskManager([&]() { ... })` で明示的に起動
- `ParallelFor` は TaskManager 起動中のみ並列化、それ以外はシーケンシャル実行
- MKL BLAS/LAPACK (dgesv_, cblas_ddot 等) は TaskManager と独立に動作する
  -> `mkl_set_num_threads()` で明示制御が必要

---

**作成日**: 2026-03-01
**対象リポジトリ**: S:\Radia\01_GitHub
**作業ブランチ**: main + peec-dev
