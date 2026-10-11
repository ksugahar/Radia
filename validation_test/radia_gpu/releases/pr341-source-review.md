# PR 341 independent source review

Date: 2026-10-11. Reviewer: independent Codex agent
`gpu_distribution_review`, operating read-only.

Reviewed range: `12c82487a..d67476ccd2a1b5300439fdb72f07a33c95ead442`
(30 changed files). The review covered the stable publisher and exact-payload
identity, CPU wheel ownership guard, legacy Simulink inventory and archive
compatibility, transfer-inclusive Poisson benchmark, installed machine/device
identity, numerical comparison acceptance, and corresponding contracts.

One P2 finding was resolved: the new numerical predicate tests imported SciPy
although the package test extra did not require it. Commit `d67476ccd` moves
`scipy.io.loadmat` into the actual MATLAB comparison execution function.
The author ran all 13 predicate tests and an import check with SciPy explicitly
unavailable. These are author checks, not reviewer execution.

Completed reviewer report:

> 独立レビュー完了。対象は PR #341 の `12c82487a..d67476ccd2a1b5300439fdb72f07a33c95ead442`、変更30ファイルです。
>
> 前回の P2 は解消しています。SciPy の import が `main()` 内へ移り、数値判定テストの収集時に SciPy／MATLAB が不要になりました。残る具体的な指摘はありません。
>
> 読み取り専用のソースレビューです。テスト・GPU実行は行っておらず、A100 の実機受入や正式版の数値・性能要件を認証するものではありません。

This record is not a stable release qualification. A final release wheel must
still satisfy its exact-source independent review and same-wheel A100/CPU
parity gates. No A100 acceptance or performance result is inferred here.
