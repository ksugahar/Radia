# Golden test の運用

## Golden test とは

> Demo で動く実装を production に porting しても動くとは限らない。
> 期待値と tolerance を JSON で lock 、pytest で定期的に verify する。

2026-04-19 の教訓: 24% の P_wp 不一致、4% mesh 感度差、20x normal 方向反転、いずれも golden 無しでは発見に週単位。

## 既存 golden tests

| test | 対象 | 所要時間 | 検証内容 |
|---|---|---|---|
| `tests/panels/test_peec_bem_golden.py` | calc_peec_bem.py | ~3 分 | P_wp -2.2% vs 2D SIBC ref (tol 15%) |
| `tests/panels/test_fem_coilmesh_golden.py` | calc_fem_coilmesh.py | ~9 分 | P_wp -1.3% vs 2D SIBC (tol 5%), L -6.4% vs 2D (tol 10%), P_coil captured ±15% |

## 回し方

```bash
cd S:/Radia/01_GitHub
pytest tests/panels/test_peec_bem_golden.py -v
pytest tests/panels/test_fem_coilmesh_golden.py -v

# まとめて (約 12 分):
pytest tests/panels/test_peec_bem_golden.py tests/panels/test_fem_coilmesh_golden.py -v
```

### pytest オプション
- `-v`: verbose (どのテストが実行されたか表示)
- `-s`: stdout/stderr を表示 (進捗見たい時)
- `--tb=short`: traceback 短く
- `-x`: 最初の失敗で止める
- `-k 'fem'`: テスト名 filter

## Golden が失敗したら

1. **まず panel が壊れたのか physics が変わったのか切り分け**:
   ```bash
   # CLI 直接実行で再現するか
   python src/radia/panels/calc_peec_bem.py \
     --peec-step src/radia/panels/samples/ih_peec_bem_coarse_coil.step \
     --vol src/radia/panels/samples/ih_peec_bem_coarse.vol \
     --frequency 7000 --current 1.0 --coil-sigma 5.8e7 \
     --sigma 5.8e7 --half-thickness 0.0125 --mu-r 1.0 \
     --peec-nwinc 3 --peec-nhinc 3 --wp-label sibc
   ```
   
2. **JSON 出力を見て diff を特定**:
   - L_coil 変わった? → PEEC 側 (filament topology, MNA)
   - P_wp 変わった? → BEM 側 (extractor, winding, H_t)
   - 両方変わった? → 共通 source (radia DLL 、numpy 等)

3. **git log で最近の変更を確認**:
   ```bash
   git log --oneline src/radia/panels/calc_peec_bem.py
   git log --oneline src/radia/panels/calc_fem_coilmesh.py
   ```

4. **tolerance 範囲内の drift なら golden JSON 更新**:
   ```bash
   # tests/panels/golden/*.json を編集
   # expected.P_wp_W, tolerance_Pwp_pct を調整
   ```
   コミットメッセージに理由明記:
   ```
   test(golden): update P_wp tolerance from 5% to 8%
   
   After coil mesh refinement (size 0.0005→0.0003 in sample .jou),
   P_wp improved from -1.3% to -0.5% vs 2D ref.
   ```

## 新しい golden を追加したい時

新 method が panel に入った時の手順:

### 1. Golden JSON を書く (tests/panels/golden/)

```json
{
  "description": "One-line description",
  "rationale": "Why this test exists",
  "sample": {
    "vol": "path/relative/to/repo",
    "generator_jou": "path/for/regeneration"
  },
  "physics": {
    "frequency_Hz": 7000,
    "current_A": 1.0,
    "...": "..."
  },
  "reference_XXX": {
    "L_nH": 50.54,
    "note": "..."
  },
  "expected": {
    "P_wp_W": 6.5e-5,
    "tolerance_Pwp_pct": 5.0
  },
  "captured_at": "YYYY-MM-DD",
  "source_of_truth": "which commit + which script produced the captured values"
}
```

### 2. Test ファイル

```python
# tests/panels/test_XXX_golden.py
import json, os, subprocess, sys
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
GOLDEN = os.path.join(HERE, "golden", "my_method_golden.json")
CALC = os.path.join(ROOT, "src", "radia", "panels", "calc_XXX.py")


def _load_golden():
    with open(GOLDEN) as f:
        return json.load(f)


@pytest.mark.slow
def test_my_method_matches_ref():
    g = _load_golden()
    vol = os.path.join(ROOT, g["sample"]["vol"])
    if not os.path.isfile(vol):
        pytest.skip(f"Sample missing: {vol}. Regenerate via .jou.")
    
    phys = g["physics"]
    cmd = [sys.executable, CALC, ...]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, ...
    result = json.loads(proc.stdout)
    
    exp = g["expected"]
    err = abs(result["P_wp_W"] - exp["P_wp_W"]) / exp["P_wp_W"] * 100
    assert err < exp["tolerance_Pwp_pct"], \
        f"Drift: {err:.2f}% (tol {exp['tolerance_Pwp_pct']}%)"
```

### 3. サンプルファイル

- .jou/.py (再生成可能) を `src/radia/panels/samples/` に
- .vol/.step (生成済み) も同じ場所 (wheel に含める)

### 4. pyproject.toml の package-data 確認

`src/radia/panels/samples/*.jou/.vol/.step` が含まれているか:
```toml
[tool.setuptools.package-data]
radia = [
    ...
    "panels/samples/*.jou",
    "panels/samples/*.vol",
    ...
]
```

## CI で golden を回すには

現状は手動実行。将来的には GitHub Actions で golden を回したい:

```yaml
# .github/workflows/golden-tests.yml (未実装)
name: Golden regression tests
on: [push, pull_request]
jobs:
  golden:
    runs-on: [self-hosted, lab-gpu]  # NGSolve + MKL 環境
    steps:
      - uses: actions/checkout@v4
      - run: pip install -e . -e packages/cubit-mesh-export
      - run: pytest tests/panels/ -v -m slow --timeout=900
```

実装時の注意:
- NGSolve + MKL + pardiso が必要 → self-hosted runner
- ndof 数十万の pardiso solve はメモリ多め (16GB+)
- slow マーク付きテストは週次 or nightly、毎 PR は軽いテストのみ

## Golden policy

- **tolerance を緩めるときは必ず commit message で justify**
- **新 method を panel に入れる時は必ず golden 追加**
- **mesh refinement で golden 値が改善したら captured を更新 + tolerance 締める**
- **golden JSON の `source_of_truth` フィールドは必ず埋める** (再現可能性)
