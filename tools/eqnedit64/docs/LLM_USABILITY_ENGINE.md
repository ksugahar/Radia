# Eqnedit64 LLM操作感エンジン

## 1. 目的

この仕組みは、Eqnedit32で蓄積された熟練操作、Eqnedit64の意味状態、利用者の
違和感を、LLMが検討できる証拠へ変換する。LLMにGUIを勝手に操作・変更させる
仕組みではない。検出器は候補を示し、LLMは仮説を作り、人間が採否を決め、
採用した動作だけをバックグラウンド回帰試験へ固定する。

操作感の判断順序は次とする。

```text
Eqnedit32原本・ヘルプ・逆アセンブル・既存の選好台帳
                          +
Eqnedit64 operation log v2（操作＋意味状態＋時間差）
                          ↓
決定論的な摩擦候補抽出（欠陥判定はしない）
                          ↓
LLMが証拠列を引用して最小変更と再現試験を提案
                          ↓
利用者が採用／却下／保留を決定
                          ↓
実装＋非表示イベント再生試験＋選好台帳
```

## 2. 操作トレース

操作ログv2は各イベントに次を記録する。

- 単調時刻による開始後時間と直前イベントからの時間差
- イベント名と詳細
- 構造キャレット、選択、編集中のTeX、dirty状態
- canvas/sourceなどの焦点
- 入力スタイル、左／中央／右の表示配置、ズーム
- `equation` / `equation*`、raw TeX欄の表示状態
- 未完了のショートカット接頭辞

旧ログv1は時刻から時間差を復元して引き続き解析できる。

## 3. 摩擦候補

`tools/analyze_usability_trace.py` は現在、次を候補として抽出する。

| 検出器 | 意味 | 単独で欠陥とみなすか |
|---|---|---|
| `explicit_user_marker` | 利用者がF12で違和感を明示 | 強い証拠だが原因はLLMと人間が判断 |
| `invalid_shortcut` | 覚えた／推測した接頭辞キーが未登録 | 旧版キーとの照合が必要 |
| `ineffective_command` | コマンドが構造状態を変えなかった | 境界での正常なno-opの場合がある |
| `immediate_undo` | 構造変更を2.5秒以内にUndo | 意図的な試行の場合がある |
| `navigation_reversal` | 1.5秒以内に逆方向へ戻った | 視認のための通常移動の場合がある |
| `correction_burst` | 4秒以内に修正操作が3回以上 | 通常の文字修正と構造摩擦を区別する |

閾値は候補窓を狭めるための値であり、自然さの点数ではない。回数を減らすために
ショートカットやキャレット動作を自動変更してはならない。

## 4. LLMレビュー束の作成

操作ログを記録し、違和感の直後にF12を押した後、次を実行する。

```powershell
pwsh -NoProfile -File build\analyze_last_operation_log.ps1
```

既定では最新ログを読み、`C:\temp\Eqnedit64-usability-*.json` を作る。
`structure`プライバシーモードでは、TeX、選択、入力内容を長さ、ハッシュ、
TeXコマンド／環境名へ置換する。内容が原因判断に必要で、ローカルの許可された
LLMだけへ渡す場合は明示的に次を使う。

```powershell
pwsh -NoProfile -File build\analyze_last_operation_log.ps1 -Privacy full
```

ログも`full`分析束も研究内容を含み得る。自動送信しない。生成器はネットワークを
使わず、マウス、キーボード、前面ウィンドウを操作しない。

## 5. LLMの出力契約

LLMは候補ごとに次を返す。

- `candidate_id`
- `verdict`: `friction` / `expected` / `uncertain`
- `hypothesis`
- `evidence_seq`: 根拠となる操作番号
- `legacy_evidence`: Eqnedit32資料または既存選好
- `proposed_minimal_change`
- `regression_test`: 外部入力を送らない再現試験
- `human_question`: 人間にしか決められない一点

本文を伏せた分析束から数式内容を推測してはならない。F12のない低信頼候補だけで
動作を変更してはならない。複数の原因を一つの大きなGUI変更へまとめず、可逆な
最小差分を提案する。

## 6. 選好台帳と完了条件

採用した判断は [`USABILITY_PREFERENCES.jsonl`](USABILITY_PREFERENCES.jsonl) へ
1行1決定で追記する。却下した案も、同じ提案を繰り返さないため結論と理由を残す。

操作感変更の完了条件は次のすべてである。

1. 操作ログまたはEqnedit32資料の根拠がある。
2. LLM提案と人間の採否が区別されている。
3. `docs/GUI_SPEC.md` が更新されている。
4. 非表示のイベント再生試験が追加されている。
5. 選好台帳へ決定と試験名が記録されている。
6. 全編集、fuzz、GUI、外部貼り付け試験が合格する。

これにより「何となく自然になった」を、証拠、判断、実行可能仕様の連鎖として
研究室に残す。
