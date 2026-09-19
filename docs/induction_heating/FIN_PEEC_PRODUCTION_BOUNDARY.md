# フィン PEEC と生産バックエンドの境界設計（2026-09-19）

対象：`radia.fin_section` / `radia.fin_sweep` / `peec_fin_topology.assemble_experimental_fin_peec`
で構成される「フィン対応表面 PEEC」を、生産経路（`panels/calc_inductance.py` と
`simulink/ih_operator_assembly.py`）へどう接続するかの契約。数値精度の受け入れ
（`validation_test/induction_heating/FIN_PEEC_IMPLEMENTATION_2026-09-19.md`）とは
独立に、**インターフェイスだけ**を決める。

## 1. 接続点は 1 つ

IH アセンブラ `_solve_unit_current` は自前で電磁解を持たず、`calc_inductance.run_inductance`
を `--coil-solver {peec,bem-a}` で呼び、ワーク側 `qsurf.sol`（1 A での SIBC 表面損失、
H1 P1）だけを受け取る。したがってフィン PEEC は **`calc_inductance` の第 3 のコイル
ソルバ `fin-surface`** として実装すれば、IH アセンブラは CLI 引数の追加のみで到達する。

```
STEP ─ fin_sweep.fin_graph_from_step ─ assemble_experimental_fin_peec ─ 枝電流 I_b
        │                                                              │
        └ 断面解析（フィン・先端・探針）                                  ▼
                                   coil_data = {source_type:"filament", paths:[枝], I_fil:I_b, ...}
                                                                       │
                       既存: A_from_filaments / phi_inc / Biot-Savart / export_filaments_msh
                                                                       ▼
                                          ワーク BEM-SIBC（weak）→ qsurf.sol → IH 作用素
```

## 2. `coil_data` 契約（変更なし、拡張のみ）

| キー | series bundle (`peec`) | `fin-surface` |
|---|---|---|
| `source_type` | `"filament"` | `"filament"`（下流は 2 点ポリラインの列として消費できる） |
| `paths` | レーンごとのポリライン | **枝ごと**の 1 セグメントポリライン `[(p0, p1)]`（縦枝＋横枝） |
| `I_fil` | レーン電流 (K,) | 枝電流 (n_branch,)（横枝も含む；KCL は保存） |
| `L_coil`, `R_coil` | 端子量 | `Re(IᴴLI)`, `Re(IᴴRI)`（シートリアクタンス除外、BEM-A 規約） |
| `n_filaments` | K | n_branch |
| `cross_section_kind` | rect / circle / unknown | `"fin-surface"` |
| `internal_impedance_model` | dowell / bessel | `"leontovich-sheet"` |
| 追加 | — | `fin_sweep`（station ごとのフィン検出）、`fin_metrics`（分担・重心・探針 H） |

下流で「レーン = 全長で 1 電流」を前提にした箇所は `_export_coil_msh_viz` の
`export_filaments_msh` を含めて無い（すべて `paths × I_fil` を独立に積分する）。
`peec_coupled_bem_solver`（strong）は K×K bundle の `R_f, L_f` を前提にするので、
`fin-surface` は **weak のみ**。strong は明示的に拒否する（gate 4 まで）。

## 3. `fin-surface` の CLI

```
--coil-solver fin-surface --coil-step coil.step
--fin-n-lanes 64  --fin-n-stations 20  --fin-lane-grading {auto,uniform}
--fin-tip-lanes 8 --fin-route {auto,straight_prism,section_planes}
```

`--peec-*` は `fin-surface` に効かない（近接反復・Dowell/Bessel は使わない）。
`--coupling-mode strong` はエラー。

## 4. IH アセンブラ側

- `IHOperatorAssemblyOptions.coil_step_solver: "peec" | "fin-surface"`（既定 `"peec"`）。
  STEP 入力のときだけ意味を持つ。`.vol` は従来どおり `bem-a`。
- ネイティブ設定の `eddy_solver` は MEX 契約 `{fem, peec, bem-a, bim}` を守るため
  **`"peec"` のまま**。区別は `eddy_method`（`"fin-surface PEEC + BEM-SIBC unit-current
  response"`）と `geometry.coil_backend`（`"fin-surface"`）に書く。ランタイム
  （`EddyRuntime`、n_unknown=1、`heat_projection = 1 A 損失密度`）は無変更。
- `unit_current` 記録に `fin_metrics`（station ごと）を添える。

## 5. 受け入れ基準（境界の妥当性、精度ではない）

1. `--coil-solver fin-surface --coil-only` が直線 fixture で `fin_peec_from_step.py`
   と同じ R / L を返す（同じ関数を呼ぶので同値であること）。
2. 同じ fixture + ワーク `.vol` で weak 結合が完走し、`qsurf.sol` が有限・非負。
3. IH アセンブラで `coil_step_solver="fin-surface"` の `native_ih.json` が
   `validateIHNativeConfig` を通る（`eddy_solver="peec"` のまま）。
4. `peec` と `fin-surface` の `qsurf` を同一ワークで比べ、**隅領域の損失分担・損失重心**
   を記録する（差の大小は判定しない。フィン側の精度受け入れは別文書）。

## 6. 後段（この文書の範囲外）

- gate 4：ワーク反作用を **枝ごとの emf** で KVL に入れる。BEM 側は枝の線積分
  `∫ A_wp·dl` を返す関数が要る。`HybridSurfaceTopology.solve_reference(emf=...)` は受け口を
  既に持つ。
- 別ソリッドのフィン（CAD で unite）、フィレット根元、断面積 ±30% フィルタ。
