# 開発者ガイド

## アーキテクチャとソースコード構成

### ディレクトリ構成

```
ngsolve-sparsesolv/
├── include/sparsesolv/           # ヘッダオンリーライブラリ本体
│   ├── sparsesolv.hpp            # メインヘッダ（全コンポーネントをinclude）
│   ├── core/                     # 基盤コンポーネント
│   │   ├── types.hpp             # 型定義 (index_t, complex_t)
│   │   ├── constants.hpp         # 数値定数（許容誤差、閾値）
│   │   ├── solver_config.hpp     # SolverConfig構造体
│   │   ├── sparse_matrix_view.hpp # 非所有CSR行列ビュー（ゼロコピー）
│   │   ├── sparse_matrix_coo.hpp # COO形式（アセンブリ用）
│   │   ├── sparse_matrix_csr.hpp # CSR形式（格納用）
│   │   ├── dense_matrix.hpp      # 密行列 + LU逆行列
│   │   ├── preconditioner.hpp    # 前処理基底クラス（テンプレート）
│   │   ├── parallel.hpp          # 並列化抽象レイヤ (TaskManager/OpenMP/逐次)
│   │   ├── level_schedule.hpp    # レベルスケジューリング（三角求解の並列化）
│   │   ├── abmc_ordering.hpp     # ABMC順序付け（三角求解の並列化）
│   │   └── rcm_ordering.hpp      # RCMバンド幅縮小順序付け
│   ├── preconditioners/          # 前処理の実装
│   │   ├── ic_preconditioner.hpp # 不完全コレスキー (IC) 分解
│   │   ├── sgs_preconditioner.hpp # 対称ガウス・ザイデル (SGS)
│   │   ├── compact_amg.hpp     # CompactAMG（古典的AMG、ヘッダオンリー）
│   │   ├── compact_ams.hpp     # CompactAMS（Hiptmair-Xu補助空間前処理）
│   │   └── complex_compact_ams.hpp # ComplexCompactAMS（渦電流向け融合Re/Im）
│   ├── solvers/                  # 反復法ソルバー
│   │   ├── iterative_solver.hpp  # 反復法ソルバー基底クラス
│   │   ├── cg_solver.hpp         # 共役勾配法 (CG)
│   │   └── sgs_mrtr_solver.hpp   # SGS-MRTR（分割公式）
│   └── ngsolve/                  # NGSolve統合レイヤ
│       ├── sparsesolv_precond.hpp # BaseMatrixラッパー (IC, SGS, Compact AMS, Solver)
│       └── sparsesolv_python_export.hpp # pybind11バインディング + ファクトリ関数
├── ngsolve/
│   └── python_module.cpp         # pybind11モジュールエントリポイント
├── tests/
│   └── test_sparsesolv.py        # ソルバーおよび前処理のテスト
├── docs/                         # ドキュメント（このフォルダ）
├── CMakeLists.txt                # ビルド設定
└── LICENSE
```

### 設計方針

#### ヘッダオンリー

すべてのC++コードは`.hpp`ヘッダファイルで実装されている。
コンパイル時にすべてインライン展開されるため、リンクの問題がなく配布も容易である。

#### テンプレート設計

すべてのアルゴリズムクラスは`template<typename Scalar>`でパラメータ化されている:

```cpp
template<typename Scalar = double>
class ICPreconditioner : public Preconditioner<Scalar> { ... };
```

`Scalar`は`double`または`std::complex<double>`でインスタンス化される。
これにより実数問題と複素数問題（例: 渦電流）を単一のコードベースで扱うことができる。

#### 並列化抽象レイヤ

`core/parallel.hpp`がコンパイル時にバックエンドを切り替える:

| ビルド設定 | バックエンド | 用途 |
|---|---|---|
| `SPARSESOLV_USE_NGSOLVE_TASKMANAGER` | NGSolve TaskManager | NGSolve統合時 |
| `_OPENMP` | OpenMP | スタンドアロン |
| （いずれも未定義） | 逐次実行 | デバッグ |

主要API:

```cpp
sparsesolv::parallel_for(n, [&](index_t i) { ... });
sparsesolv::parallel_reduce(n, init, [&](index_t i) -> T { ... });
sparsesolv::get_num_threads();
```

### ヘッダ依存関係

```
sparsesolv.hpp (メインヘッダ)
├── core/types.hpp
├── core/constants.hpp
├── core/solver_config.hpp
├── core/sparse_matrix_view.hpp ← core/parallel.hpp
├── core/preconditioner.hpp ← core/sparse_matrix_view.hpp
├── core/abmc_ordering.hpp
├── preconditioners/ic_preconditioner.hpp
│   ← core/preconditioner.hpp, core/level_schedule.hpp,
│      core/abmc_ordering.hpp, core/rcm_ordering.hpp
├── preconditioners/sgs_preconditioner.hpp
│   ← core/preconditioner.hpp
├── solvers/iterative_solver.hpp
│   ← core/preconditioner.hpp, core/sparse_matrix_view.hpp
├── solvers/cg_solver.hpp
│   ← solvers/iterative_solver.hpp
└── solvers/sgs_mrtr_solver.hpp
    ← solvers/iterative_solver.hpp, core/level_schedule.hpp
```

NGSolve統合レイヤ（NGSolveとビルドする場合のみ）:

```
ngsolve/sparsesolv_precond.hpp   ← sparsesolv.hpp + NGSolveヘッダ
ngsolve/sparsesolv_python_export.hpp ← sparsesolv_precond.hpp + pybind11
```

### NGSolve統合レイヤ

#### SparseMatrixView（ゼロコピーラッパー）

NGSolveの`SparseMatrix<SCAL>`はCSR的な内部構造を持つ。
`SparseSolvPrecondBase::prepare_matrix_view()`がこれをSparseSolvの
`SparseMatrixView<SCAL>`にデータコピーなしで変換する。

FreeDofs処理: 拘束DOFの行は単位行（対角=1、非対角=0）に置換される。

```cpp
// sparsesolv_precond.hpp: SparseSolvPrecondBase::prepare_matrix_view()
if (!freedofs_->Test(i)) {
    modified_values_[k] = (j == i) ? SCAL(1) : SCAL(0);  // identity row
} else if (!freedofs_->Test(j)) {
    modified_values_[k] = SCAL(0);  // zero coupling to constrained DOF
}
```

#### SparseSolvPrecondBase（BaseMatrixラッパー）

すべての前処理の基底クラス。NGSolveの`BaseMatrix`を継承し、
NGSolveの`CGSolver`との互換性のために`Mult()` / `MultAdd()`を実装する。

```
SparseSolvPrecondBase<SCAL>  （抽象基底）
├── SparseSolvICPreconditioner<SCAL>    → ICPreconditioner<SCAL>
├── SparseSolvSGSPreconditioner<SCAL>   → SGSPreconditioner<SCAL>
└── SparseSolvSolver<SCAL>              → ICCG/SGSMRTR/CG

CompactAMS  （BaseMatrixを直接継承）
  → 実数HCurl AMS前処理（静磁場、Update()対応）
  → Python: CompactAMSPreconditionerImpl

ComplexCompactAMS  （BaseMatrixを直接継承）
  → 複素渦電流向け融合Re/Im AMS前処理（Update()対応）
  → Python: ComplexCompactAMSPreconditionerImpl
```

非線形ソルバーへの対応（Update()）:
- `Update()`: 現在の行列で前処理を再構築（幾何情報は保持）
- `Update(new_mat)`: 新しい行列で前処理を再構築
- 幾何情報（G行列、Pi行列、転置行列）は初回構築時のみ計算される

#### Pythonファクトリ関数（自動ディスパッチ）

`sparsesolv_python_export.hpp`のファクトリ関数は行列の型を自動判定し、
適切なテンプレートインスタンスを生成する:

```cpp
// ICPreconditioner factory
m.def("ICPreconditioner", [](shared_ptr<BaseMatrix> mat, ...) {
    if (mat->IsComplex()) {
        // → SparseSolvICPreconditioner<Complex>
    } else {
        // → SparseSolvICPreconditioner<double>
    }
});
```

#### 型登録 (AMS)

CompactAMS / ComplexCompactAMSは`py::class_`で具象型を登録し、
ファクトリ関数が具象型を返すことで、Pythonから`Update()`にアクセス可能になる:

```cpp
// 型登録
py::class_<CompactAMS, shared_ptr<CompactAMS>, BaseMatrix>
    (m, "CompactAMSPreconditionerImpl")
    .def("Update", py::overload_cast<>(&CompactAMS::Update))
    .def("Update", py::overload_cast<shared_ptr<SparseMatrix<double>>>(&CompactAMS::Update));

// ファクトリ（具象型を返す）
m.def("CompactAMSPreconditioner", [...] -> shared_ptr<CompactAMS> { ... });
```

アルゴリズムの詳細は[algorithms.md](algorithms.md)を参照。

---

## ビルド手順

### 前提条件

- NGSolve（ソースビルドまたはpip install）+ CMake設定ファイル
- CMake 3.16以上
- C++17対応コンパイラ（MSVC 2022、GCC 10+、Clang 10+）
- pybind11（CMakeが自動取得）

### ビルドコマンド

```bash
git clone https://github.com/ksugahar/ngsolve-sparsesolv.git
cd ngsolve-sparsesolv
mkdir build && cd build

cmake .. -DSPARSESOLV_BUILD_NGSOLVE=ON \
         -DNGSOLVE_INSTALL_DIR=/path/to/ngsolve/install/cmake

cmake --build . --config Release
```

**Intel MKL使用時**（CompactAMG/AMS内部で使用するBLAS）:

```bash
cmake .. -DSPARSESOLV_BUILD_NGSOLVE=ON \
         -DNGSOLVE_INSTALL_DIR=/path/to/ngsolve/install/cmake \
         -DCMAKE_PREFIX_PATH="C:/Program Files (x86)/Intel/oneAPI/mkl/latest"
```

### インストール

```bash
# NGSolveのsite-packagesにコピー
SITE_PACKAGES=$(python -c "import ngsolve, pathlib; print(pathlib.Path(ngsolve.__file__).parent.parent)")
mkdir -p "$SITE_PACKAGES/sparsesolv_ngsolve"
cp build/Release/sparsesolv_ngsolve*.pyd "$SITE_PACKAGES/sparsesolv_ngsolve/"
echo "from .sparsesolv_ngsolve import *" > "$SITE_PACKAGES/sparsesolv_ngsolve/__init__.py"
```

---

## テストの実行

```bash
python -m pytest tests/test_sparsesolv.py -v --tb=short
```

テスト構成:
- `test_sparsesolv.py`: ソルバーおよび前処理のテスト（46ケース）
  - ICCG、SGSMRTR、CG、IC、SGS（各種問題に対して）
  - 2D/3D、H1/VectorH1/HCurl、実数/複素数
  - ABMC順序付け、対角スケーリング、自動シフト

---

## 過去のバグと教訓

### 1. DenseMatrixのP^Tバグ

**問題**: 密行列LU逆行列において、置換行列がPではなくP^Tとして構築されていた。

```cpp
// 誤り (P^T):
inv(piv[j], j) = Scalar(1);

// 正解 (P):
inv(k, piv[k]) = Scalar(1);
```

**教訓**: PA = LU において、Pは「行kを行piv[k]に移動する」行列である。
正しい構築は P[k, piv[k]] = 1 となる。

### 2. 複素内積（非共役と共役）

**問題**: CGソルバーが共役内積 `std::conj(a[i]) * b[i]` を使用していた。

FEM行列は**複素対称**（A^T = A）であり、エルミート（A^H = A）ではない。
例: 渦電流方程式 curl-curl + i*sigma*mass。

```cpp
// 誤り（エルミート行列向け）:
sum += std::conj(a[i]) * b[i];

// 正解（複素対称行列向け）:
sum += a[i] * b[i];
```

**症状**: 渦電流問題で5000回反復しても発散。
修正後: 58回で収束。

**教訓**: 複素FEM問題では常に非共役内積を使用すること。
これはNGSolveの`CGSolver(conjugate=False)`に対応する。

### 3. SGS-MRTRの複素比較

**問題**: SGS-MRTRのzeta計算で`denom >= 0`を使用していたが、
`std::complex<double>`には`>=`演算子が存在しない。

```cpp
// 修正後:
denom = (std::real(denom) >= 0) ? ... : ...;
```

---

## NGSolve APIに関する注意事項

### AutoVector

NGSolveの`BaseVector::CreateVector()`は`AutoVector`を返す。
shared_ptrに似たセマンティクスを持つが、独自のムーブ/コピー規約がある。

### CalcElementMatrix

要素行列の計算には`LocalHeap`が必要:

```cpp
LocalHeap lh(10000000, "name", true);  // mult_by_threads=true

IterateElements(*fes, VOL, lh,
    [&](FESpace::Element el, LocalHeap& lh_thread) {
        auto& fe = el.GetFE();
        auto& trafo = el.GetTrafo();
        FlatMatrix<SCAL> elmat(ndof, ndof, lh_thread);
        integrator->CalcElementMatrix(fe, trafo, elmat, lh_thread);
    });
```

### IsRegularDof

要素のDOFリストには無効なDOF（< 0）が含まれる場合がある。
`IsRegularDof(dnum)`でフィルタする:

```cpp
auto dnums = el.GetDofs();
for (int i = 0; i < dnums.Size(); ++i) {
    if (IsRegularDof(dnums[i])) {
        // 有効なDOF
    }
}
```

---

## 参考文献

1. J. A. Meijerink, H. A. van der Vorst,
   "An Iterative Solution Method for Linear Systems of Which the
   Coefficient Matrix is a Symmetric M-Matrix",
   *Math. Comp.*, Vol. 31, No. 137, pp. 148-162, 1977.
   [DOI: 10.1090/S0025-5718-1977-0438681-4](https://doi.org/10.1090/S0025-5718-1977-0438681-4)

2. M. R. Hestenes, E. Stiefel,
   "Methods of Conjugate Gradients for Solving Linear Systems",
   *J. Research of the National Bureau of Standards*,
   Vol. 49, No. 6, pp. 409-436, 1952.
   [DOI: 10.6028/jres.049.044](https://doi.org/10.6028/jres.049.044)

3. T. Iwashita, H. Nakashima, Y. Takahashi,
   "Algebraic Block Multi-Color Ordering Method for Parallel Multi-Threaded
   Sparse Triangular Solver in ICCG Method",
   *Proc. IEEE IPDPS*, 2012.
   [DOI: 10.1109/IPDPS.2012.51](https://doi.org/10.1109/IPDPS.2012.51)

4. E. Cuthill, J. McKee,
   "Reducing the Bandwidth of Sparse Symmetric Matrices",
   *Proc. 24th Nat. Conf. ACM*, pp. 157-172, 1969.
   [DOI: 10.1145/800195.805928](https://doi.org/10.1145/800195.805928)

5. 圓谷友紀 (T. Tsuburaya), 三船泰 (Y. Mifune), 岩下武史 (T. Iwashita), 高橋英治 (E. Takahashi),
   "MRTR法に基づく前処理付き反復法の数値実験"
   (Numerical experiments on preconditioned iterative methods based on the MRTR method),
   *電気学会研究会資料* (IEEJ Technical Meeting Papers), SA-12-64, 2012.

6. T. Tsuburaya, Y. Okamoto, K. Fujiwara, S. Sato,
   "Improvement of the Preconditioned MRTR Method With Eisenstat's Technique
   in Real Symmetric Sparse Matrices",
   *IEEE Trans. Magnetics*, Vol. 49, No. 5, pp. 1641-1644, 2013.
   [DOI: 10.1109/TMAG.2013.2240283](https://doi.org/10.1109/TMAG.2013.2240283)

7. 圓谷友紀 (T. Tsuburaya),
   "大規模電磁界問題の有限要素解析のための反復法の開発"
   (Development of iterative methods for finite element analysis of large-scale electromagnetic field problems),
   Doctoral dissertation, 宇都宮大学 (Utsunomiya University), 2016.

8. JP-MARs/SparseSolv,
    https://github.com/JP-MARs/SparseSolv
