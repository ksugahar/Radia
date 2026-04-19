# IH Panel GUI テスト チェックリスト (2026-04-19 refactored)

パネル refactor 後の GUI 検証手順。

## 事前確認 (claude 実施済み)

| 項目 | 結果 |
|---|---|
| Panel import + method combo (2 択) | ✔ |
| Widget 数 (19) + section header 構造 | ✔ |
| Material presets (Cu/Al/Brass/Steel100/Steel500/Stainless/Custom) × 2 (coil + wp 独立) | ✔ |
| Impedance model (Linear SIBC / ESIM WIP) + BH file/max_iter/tol | ✔ |
| Linear solver per method (PEEC: Dense LU/HACApK, FEM: pardiso/shifted_ams/BDDC/iccg) | ✔ |
| .vol label validator on load | ✔ |
| Physics sanity check (delta, R/delta warning) | ✔ |
| Formatted output (IH Summary section) | ✔ |
| `--help` on both calc scripts (WIP options plumbed) | ✔ |
| Golden tests PASS (705s) | ✔ |
| 100号機 site-packages sync | ✔ |

## 当日の GUI テスト (ユーザ実施)

### A. 起動 + status label
1. Cubit で `ih_peec_bem_coarse.jou` を再生 (.vol 生成)
2. Tools → Solve → Radia-NGSolve → IH window 起動
3. Status label が **緑で ".vol OK"** + **"WP delta = ... R/delta = ..."** 表示されるか

### B. PEEC+BEM (Fast) method 初期動作
1. Method = "Fast workpiece heating (PEEC+BEM, 1-way)" (default)
2. Tooltip (hover) で "±5%, 3 min, 1-way forward..." 出るか
3. 見える widgets:
   - Method, Drive (freq, current)
   - Coil material (preset: Copper, sigma 5.8e7, mu_r 1.0 disabled)
   - Coil geometry (PEEC+BEM only): STEP file
   - Workpiece material (preset + sigma/mu_r/half_thickness)
   - Workpiece impedance model (Linear SIBC)
   - Linear solver (Dense LU / HACApK)
   - Advanced: PEEC nwinc/nhinc
4. 見えない widgets:
   - FES order (FEM 専用)
   - BH file / ESIM max iter / tolerance (ESIM 選択時のみ)

### C. Material preset
1. Coil material を `Stainless 304` に変更 → sigma = 1.4e6, mu_r = 1 auto-fill、disabled
2. Coil material を `Custom` に変更 → sigma/mu_r lineedit 編集可能に
3. WP material = `Copper` に変更しても Coil は変わらないこと (独立性)

### D. FEM method に切替
1. Method = "Full simulation (FEM A-V + wp SIBC + Kelvin)"
2. Widget 可視性変化:
   - 消える: PEEC geometry (STEP file, nwinc, nhinc)
   - 出る: FES order
   - Solver 選択肢: pardiso/shifted_ams/BDDC/iccg
3. もし FEM サンプル (.vol) を load してなければ status label に "ERROR: Missing material 'coil'" 等の赤文字
4. ih_fem_kelvin_skin_fine.vol を load → 全緑

### E. ESIM UI
1. Impedance model = "Nonlinear ESIM (experimental, WIP)" 選択
2. 出る widgets: BH file browser, max_iter, tolerance
3. Run ボタン押す → `{"error":"ESIM is not implemented..."}` stderr
4. Impedance model を Linear SIBC に戻す → ESIM widgets 消える

### F. Physics sanity warning
1. WP material = `Copper`、周波数 = 100 Hz に変更
2. delta 計算される → status label 更新
3. R/delta < 3 なら amber warning 表示

### G. 実行 (PEEC+BEM)
1. Method: Fast workpiece heating
2. PEEC STEP: `ih_peec_bem_coarse_coil.step`
3. .vol: `ih_peec_bem_coarse.vol`
4. freq=7000, current=1.0, Cu coil, Cu wp, mu_r=1, half_thickness=0.0125
5. Run → 約 3 分後、"IH Summary" 出力:
   - L_coil (vacuum): 88.57 nH
   - P_workpiece: ~6.48e-5 W (0.065 mW)
   - Heating efficiency: -- (P_coil が無いので N/A)

### H. 実行 (FEM A-V)
1. Method: Full simulation
2. .vol: `ih_fem_kelvin_skin_fine.vol`
3. freq=7000, current=1.0, Cu coil, Cu wp, mu_r=1, half_thickness=0.0125, FES order=1
4. Run → 約 9 分後、IH Summary:
   - L_total: 47.32 nH
   - P_workpiece: 6.54e-5 W
   - P_coil: 1.48e-4 W (mesh-sensitive ±15%)
   - P_total: 2.13e-4 W
   - Heating efficiency: ~30%
   - coil delta/h_max を表示 (0.79/0.97 → WARN)

### I. Settings persistence
1. Method, Materials, freq, current を変更して window 閉じる
2. 再起動 → text-based save で同じ状態に復元

## 判定

| テスト | 合格基準 |
|---|---|
| A-F (UI 動作) | 可視性切替、material preset、validator、warning が想定通り |
| G (PEEC+BEM) | L_coil ±1%, P_wp ±15% vs captured |
| H (FEM A-V) | L ±10%, P_wp ±5%, P_coil ±15% vs captured |
| I | method + materials が復元 |

全合格 → ユーザへの引き継ぎ完了。

## サポート資料

- `docs/kubota/00_README.md` 以下全 7 本
- `memory/ih_panel_final_state_2026_04_19.md` — 定式化経緯
- `tests/panels/test_*_golden.py` — 自動テスト (CLI)
- MCP `ih_knowledge("peec_bem_sibc")` と `("av_coil_sigma")` — 物理詳細

## 不具合 playbook

1. **panel が method combo 空**: settings save 形式の互換問題 → settings ファイル削除
2. **.vol 読んでも status label 出ない**: `inspect_vol_labels` 失敗、NGSolve import error の可能性
3. **Run ボタン disabled**: .vol 必要 label 欠落 → status label の ERROR 文確認
4. **calc subprocess 落ちる**: stderr 見る、`tests/panels/test_*_golden.py` で CLI 再現
5. **ESIM 選んで error**: 想定どおり WIP、Linear SIBC に戻す
