# MCP server quickstart (LLM 連携)

## MCP 概要

Radia には 3 つの MCP server:
- `mcp-server-radia` (radia_ngsolve): NGSolve + BEM + HCurl/HDiv + sparsesolv + Kelvin
- `mcp-server-ih`: IH 解析の method 選択、SIBC, ESIM, PEEC+BEM, A-V
- `mcp-server-cubit`: Cubit scripting、export format、panel convention

Claude Code の `.mcp.json` で LAB では 3 つ全部、mdx では PyPI 経由で同じく 3 つ。

## 使い方 (Claude Code との対話例)

### Q: IH 解析の panel、PEEC+BEM と FEM どっち使うべき?

→ Claude は `mcp__radia_ih__ih_knowledge(topic="peec_bem_sibc")` と `("av_coil_sigma")` 呼んで引用して答える。

### Q: A-V 定式化で source/sink Dirichlet lift どう書く?

→ Claude は `ih_knowledge("av_coil_sigma")` を引いて、MCP 内にあるコードスニペット (`src/radia/mcp_server/ih/ih_knowledge.py:1513-1545` 相当) を表示。

### Q: NGSolve で HCurl SIBC Robin BC 書きたい

→ `mcp__radia_ngsolve__ngsolve_rules()` で既知 gotcha list + Robin BC 定式化を返す。

### Q: Cubit で imprint+merge が効かないんだけど

→ `mcp__radia_cubit__cubit_scripting("sideset_pitfalls")` 等で共通パターン。

## Topics リスト

### mcp-server-ih
- `ih_knowledge`: overview, gmsh_mesh, eddy_current, thermal, rotating, postprocess, pitfalls, esim_kelvin, **peec_bem_sibc (2026-04-19)**, **av_coil_sigma (2026-04-19)**, failed_approaches
- `ih_sibc`: all, peec_fem (historical), overview, esim, biot_savart, screening

### mcp-server-radia-ngsolve
- `ngsolve_rules`: HCurl/HDiv gotchas (nograds, Periodic, trace 等)
- `kelvin_knowledge`: Kelvin transform + 2-sphere implementation
- `ngsbem_inductance_knowledge`: BEM inductance (LaplaceSL saddle point EFIE)
- `radia_knowledge`: MMM/MSC physics
- `panel_gui_pitfalls`: combo state, mode switch, GMSH viz, subprocess args 等
- `panel_describer`: panel parameter auto-documentation
- `sparsesolv_knowledge`: Compact AMS/COCR preconditioner

### mcp-server-cubit
- `cubit_scripting_knowledge`: block, sideset, group, webcut pattern
- `export_knowledge`: GMSH v4.1, Nastran BDF, VTK, Netgen vol
- `netgen_workflow_knowledge`: .vol order, curvedelements
- `panel_conventions_knowledge`: ラベル命名規則

## MCP server files (コード的な位置)

```
src/radia/mcp_server/
├─ ih/
│  ├─ server.py                 # MCP tool entry (ih_knowledge, ih_sibc)
│  ├─ ih_knowledge.py           # Main topics (overview, peec_bem_sibc, av_coil_sigma)
│  └─ sibc_knowledge.py         # IH architecture detail
├─ radia_ngsolve/
│  ├─ server.py                 # NGSolve tools
│  ├─ radia_knowledge.py
│  ├─ ngsolve_knowledge.py
│  ├─ rules.py                  # Static rule checks (HCurl nograds 等)
│  ├─ kelvin_knowledge.py
│  ├─ ngsbem_inductance_knowledge.py
│  └─ panel_gui_pitfalls_knowledge.py
└─ cubit/
   ├─ server.py
   ├─ cubit_scripting_knowledge.py
   ├─ export_knowledge.py
   ├─ netgen_workflow_knowledge.py
   ├─ panel_conventions_knowledge.py
   └─ rules.py
```

## 知識を追加する時 (開発者向け)

1. トピック追加 (例: `ih_knowledge.py` の新セクション):
   ```python
   INDUCTION_HEATING_MY_NEW_TOPIC = """
   # My new physics concept
   
   ## Motivation
   ...
   
   ## Formulation
   ...
   
   ## Validation
   | Method | L | P |
   |---|---|---|
   | reference | ... | ... |
   | implementation | ... | ... |
   
   ## References
   - `src/radia/panels/calc_xxx.py` — implementation
   - `tests/panels/test_xxx_golden.py` — regression
   """
   ```

2. Topic dict に登録:
   ```python
   def get_induction_heating_documentation(topic: str = "all") -> str:
       topics = {
           ...
           "my_new_topic": INDUCTION_HEATING_MY_NEW_TOPIC,
       }
   ```

3. 動作確認:
   ```bash
   python -c "
   import sys; sys.path.insert(0, 'src/radia')
   from mcp_server.ih import ih_knowledge
   print(ih_knowledge.get_induction_heating_documentation('my_new_topic')[:500])
   "
   ```

4. Claude Code restart or /mcp reload で知識が読み込まれる

## LAB / 100号機 で MCP が動くメカニズム

- 両機とも `S:/Radia/01_GitHub` (LAB) = `W:/00_CAE/Radia/01_GitHub` (100号機) の **同一 SMB share**
- `.mcp.json` で server の python path + file path 指定:
  ```json
  {
    "mcpServers": {
      "radia-ih": {
        "command": "python",
        "args": ["C:/Program Files/Python312/Lib/site-packages/radia/mcp_server/ih/server.py"]
      }
    }
  }
  ```
- `src/radia/mcp_server/...` を編集 → 即座に両機で反映 (restart 必要)

mdx (一般 PyPI) では `pip install radia-mcp` で別 install。挙動は同じ。

## 使いこなし tips

- **Claude に panel を書かせる時**: 先に `ih_knowledge("av_coil_sigma")` 呼ばせると既知パターン再利用される
- **Retired な option を Claude に見せたくない時**: `failed_approaches` topic で明示
- **ML 的な behavior 修正**: `panel_gui_pitfalls` で既知 bug pattern 登録すると LLM が同じ間違いしにくい
- **PR review に MCP を活用**: "このコードで HCurl nograds=False になってるけどいい?" みたいな自動指摘
