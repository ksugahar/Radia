# 久保田君 初週 onboarding checklist

## 環境構築 (Day 1)

### アカウント / アクセス
- [ ] GitHub `ksugahar/Radia` の contributor 権限取得 (菅原に依頼)
- [ ] Slack `#radia-dev` 参加
- [ ] LAB PC に座れる or 100号機 SSH access

### LAB PC (または自分の開発機)
- [ ] SMB share `//labserver/work/` マウント確認 (= LAB の `S:\`)
- [ ] 100号機 への SSH alias 設定 (`ssh 100`)
  ```
  # ~/.ssh/config
  Host 100
      HostName 192.168.11.100
      User Administrator
  ```
- [ ] Cubit 2025.3 インストール済みか確認
- [ ] Python 3.12 + pip install radia + cubit-mesh-export
- [ ] `cubit-plugin-install` 実行 (Cubit plugin 配置)
- [ ] Claude Code 初期設定 + `.mcp.json` 配置 (MCP 3 server 設定)

### 動作確認
- [ ] Cubit 起動 → Tools → Play Journal File で `src/radia/panels/samples/ih_peec_bem_coarse.jou` を再生
- [ ] Tools → Solve → Radia-NGSolve → IH window 起動
- [ ] PEEC+BEM 検証を走らせて P_wp ≈ 6.48e-5 W が出るか確認

## アーキテクチャ理解 (Day 2-3)

### 読む資料
- [ ] `docs/kubota/00_README.md` (この folder 全体)
- [ ] `docs/kubota/01_architecture.md`
- [ ] `CLAUDE.md` (project-wide policy)
- [ ] MCP `ih_knowledge("peec_bem_sibc")` と `("av_coil_sigma")` を Claude に聞いて要約してもらう

### 手を動かすタスク
- [ ] `src/radia/panels/calc_peec_bem.py` を読んで、どの順で処理してるか俯瞰
- [ ] `src/radia/panels/calc_fem_coilmesh.py` を同様に
- [ ] Golden test 2 本を走らせて、期待通り通ることを確認:
  ```bash
  pytest tests/panels/test_peec_bem_golden.py -v
  pytest tests/panels/test_fem_coilmesh_golden.py -v
  ```
- [ ] `ih_fem_kelvin_skin_fine.jou` を読んで、.vol がどう生成されるか追う

## 研究課題に入る前 (Day 4-5)

### 質問リスト作成
Validate 理解:
- [ ] PEEC+BEM の "1-way" が何を意味するか 30 秒で説明できる
- [ ] A-V 定式化でなぜ source/sink Dirichlet lift が必要か説明できる
- [ ] Neumann-K がなぜ retired されたか説明できる
- [ ] Kelvin transformation が open boundary を解く原理を理解
- [ ] SIBC (Surface Impedance Boundary Condition) と volumetric σ の使い分け

不明点は Slack で質問 or Claude に聞く (MCP を引用してくれる)。

### 小タスク (初週の仕上げ)
- [ ] `ih_peec_bem_coarse.jou` を改造して、major radius R=40mm の版を作る
- [ ] `panel` で走らせて P_wp がどう変わるか観察
- [ ] 結果を `memory/` に書くか、`examples/` のスクリプトで参考値として残すか決める

## 2 週目 以降

### 実際の研究課題に着手
研究内容に応じて:
- **wp back-reaction 付き PEEC+BEM (2-way)**: calc_peec_bem 拡張。行列連成 or 2 次摂動
- **Nonlinear SIBC (ESIM) の panel 統合**: calc_fem_coilmesh に ESIM Karl iteration 追加
- **Multi-coil / multi-wp**: 現 panel は 1 コイル 1 ワーク前提。拡張
- **HACApK-BEM 化**: BEM assembly O(N²) → O(N log N)

### 研究発表の流れ
1. `examples/` に prototype  
2. `tests/` に golden 付けて validate  
3. `docs/research/` に internal memo (gitignored) or `docs/papers/` に public draft  
4. panel 昇格を決めたら: method combo に追加 + MCP knowledge 更新 + golden lock

## よくある先輩アドバイス

### Do
- Claude に聞いて MCP 引用させる (知識の再発明防止)
- 小さい変更 → pytest → 確認 の cycle を短く
- .jou 変更したら必ず `.vol` を Cubit で再生成、古い .vol と混同しない
- `memory/` と `docs/kubota/` を更新してから push

### Don't
- `.vol` / `.step` / `.pyd` を git commit しない (bin は release で配る)
- 目視で "たぶん合ってる" → panel 昇格しない (必ず golden)
- PEEC と FEM を同じ .vol で動かそうとしない (topology 要件違う、別 .jou)
- Cubit 開いたまま deploy (`.ccm`/`.ccl` ロック)

## エスカレーション

- 物理式の疑問 → 菅原  
- パッケージング / deploy 周り → 菅原 or LAB (see `memory/` "Deploy Flow")
- MCP / LLM の活用法 → Claude Code directly, or Slack  
- 論文作成 → 菅原 + 共同研究者

## チェックインミーティング

- 週次 15 分程度、Slack DM or 対面
- 週末までに上記タスクの進捗を更新 (check 記入)
- 不明点は翌週のミーティング前までに Slack で共有
