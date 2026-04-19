# 久保田君 引き継ぎ資料 (2026-04-19)

Radia IH パネル開発の引き継ぎ資料一覧。順に読むことを推奨。

## 読む順序

| # | ファイル | 所要時間 | 目的 |
|---|---|---|---|
| 1 | [01_architecture.md](01_architecture.md) | 15 分 | アーキテクチャ全体像、2 panel methods の原理 |
| 2 | [02_method_selection.md](02_method_selection.md) | 10 分 | どの method をいつ使うか判断フロー |
| 3 | [03_jou_template.md](03_jou_template.md) | 20 分 | .jou ファイル書き方、annotated examples |
| 4 | [04_golden_tests.md](04_golden_tests.md) | 15 分 | Regression test の回し方、追加方法 |
| 5 | [05_troubleshooting.md](05_troubleshooting.md) | 参考 | 問題切り分け playbook |
| 6 | [06_mcp_quickstart.md](06_mcp_quickstart.md) | 10 分 | MCP server の活用 (LLM 連携) |
| 7 | [07_onboarding_checklist.md](07_onboarding_checklist.md) | 1 週間分 | 初週のタスクリスト |
| A | [IH_PANEL_GUI_TEST_CHECKLIST.md](IH_PANEL_GUI_TEST_CHECKLIST.md) | 参考 | GUI 検証手順 (本番前毎回) |

## 全体の文脈 (3 分サマリ)

Radia は誘導加熱 (IH) 解析の open-source tool。Cubit でメッシュ生成 → Radia の 2 panel methods で解析:

- **PEEC+BEM (1-way)**: ワーク加熱量 P_wp を速く (3 分) 出したい時。STEP ファイル + 粗い表面メッシュで OK。
- **FEM A-V (coil meshed + SIBC + Kelvin)**: 正確な L + P_wp + P_coil 欲しい時。コイル体積メッシュが必要 (7-9 分)。

2D 軸対称 reference と両 method で **P_wp 1% 一致**、L 6%、P_coil はメッシュ依存。

基本的な開発フロー (ladder 原則):
```
tests/ (数学) → examples/ (電磁場研究) → panels/ (工学応用、GUI)
```

新しい method を panel に追加するときは golden test 必須。mesh / physics / CLI の変更は `tests/panels/test_*_golden.py` で lock される。

## 連絡先

- Slack: #radia-dev
- GitHub: https://github.com/ksugahar/Radia
- 菅原 (ksugahar@ele.kindai.ac.jp)
