# .jou テンプレート + annotated 例

## Cubit .jou (journal) の基本

IH 解析の .jou は以下を生成:
- コイル体積 + ワーク (hole or 体積) + 空気 + Kelvin 外部
- ブロック (材料ラベル): `coil`, `air`, `kelvin` (必要に応じて `workpiece`)
- サイドセット (境界ラベル): `source`, `sink` (port), `sibc` (wp 表面), `coil_surface` (必要なら)
- メッシュ
- STEP 出力 (コイルのみ、PEEC+BEM 用)
- .vol 出力 (Netgen 形式、計算用)

## Method 別の .jou 要件

### PEEC+BEM (1-way)

**最小構成 (参考: [ih_peec_bem_coarse.jou](../../src/radia/panels/samples/ih_peec_bem_coarse.jou))**:

```bash
reset
# 1. コイル (gapped torus; gap ≥ 5 deg 推奨)
create vertex 0.030000 0 0
create vertex 0.033000 0 0
create vertex 0.030000 0 0.003000
create curve arc center vertex 1 vertex 2 vertex 3 normal 0 1 0 full
create surface curve 1
sweep surface 1 axis 0 0 0 0 0 1 angle 355   # 5 deg gap
${coil_vid = Id("volume")}

# 2. ワーク (穴として subtract、BEM は surface のみ使用)
create cylinder radius 0.025000 height 0.025000
${wp_vid = Id("volume")}

# 3. 空気 (webcut で対称化、接続メッシュ用)
create sphere radius 0.060000
${air_vid = Id("volume")}
${id_before_webcut = Id("volume")}
webcut volume {air_vid} with plane zplane
${air_top = air_vid}
${air_bot = id_before_webcut + 1}

# 4. subtract + imprint + merge
subtract volume {coil_vid} {wp_vid} from volume {air_top} {air_bot} keep_tool
imprint volume {coil_vid} {wp_vid} {air_top} {air_bot}
merge volume {coil_vid} {wp_vid} {air_top} {air_bot}

# 5. STEP export (PEEC フィラメント用、コイルのみ)
export step "ih_my_sample_coil.step" volume {coil_vid} overwrite

# 6. メッシュ (coarse でも OK、BEM 精度は wp 表面密度のみ依存)
volume {coil_vid} scheme tetmesh; volume {coil_vid} size 0.003; mesh volume {coil_vid}
volume {air_top} {air_bot} scheme tetmesh; volume {air_top} {air_bot} size 0.015
mesh volume {air_top} {air_bot}

# 7. ブロック
block 1 add volume {coil_vid}
block 1 name "coil"
block 2 add volume {air_top} {air_bot}
block 2 name "air"

# 8. サイドセット: wp 表面を "sibc" にする (BEM が読む)
sideset 1 add surface in volume {wp_vid}
sideset 1 name "sibc"

# 9. (optional) gap 面を source/sink に (FEM method でも使う場合)
group "coil_gaps" add surface in volume {coil_vid} with area < 0.0001
sideset 2 add surface in coil_gaps with y_coord > -0.001
sideset 2 name "source"
sideset 3 add surface in coil_gaps with y_coord < -0.001
sideset 3 name "sink"

# 10. Netgen vol export (BEM で使うのは表面のみ、vol 経由で BEM surface 抽出)
radia_export netgen "ih_my_sample.vol" order 1 overwrite
```

**注**: PEEC+BEM では `source`/`sink` は必須ではない (wp の `sibc` のみ必要)。でも FEM と互換に作るなら同じ .jou で両 method カバーできて便利。

### FEM A-V (coil meshed + SIBC + Kelvin)

**必須要件** (参考: [ih_fem_kelvin_skin_fine.jou](../../src/radia/panels/samples/ih_fem_kelvin_skin_fine.jou)):

```bash
reset
# 1. コイル (gap + source/sink face 必須)
create vertex 0.030000 0 0
create vertex 0.033000 0 0
create vertex 0.030000 0 0.003000
create curve arc center vertex 1 vertex 2 vertex 3 normal 0 1 0 full
create surface curve 1
sweep surface 1 axis 0 0 0 0 0 1 angle 355   # 5 deg gap
${coil_vid = Id("volume")}

# 2-4. ワーク・空気・subtract/imprint/merge は PEEC+BEM と同じ

# 5. STEP (PEEC との比較用、optional)
export step "ih_my_fem_coil.step" volume {coil_vid} overwrite

# 6. メッシュ (CRITICAL: コイルは delta/3 を目標に fine 化)
# Cu @ 7 kHz: delta = 0.79 mm
#   cross-section circles (length < 0.05 m): 120 intervals -> 0.16 mm
#   toroidal arcs (length > 0.05 m):         240 intervals -> 0.78 mm
#   volume interior:                          size 0.0005 (0.5 mm)
curve in volume {coil_vid} with length < 0.05 interval 120
curve in volume {coil_vid} with length > 0.05 interval 240
volume {coil_vid} scheme tetmesh
volume {coil_vid} size 0.0005
mesh volume {coil_vid}

# 空気も相応に
volume {air_top} {air_bot} scheme tetmesh
volume {air_top} {air_bot} size 0.010
mesh volume {air_top} {air_bot}

# 7. ブロック: coil は必須 (A-V で material として参照)
block 1 add volume {coil_vid}
block 1 name "coil"
block 2 add volume {air_top} {air_bot}
block 2 name "air"

# 8. サイドセット: source, sink, sibc ALL 必須
group "coil_gaps" add surface in volume {coil_vid} with area < 0.0001
sideset 1 add surface in coil_gaps with y_coord > -0.001
sideset 1 name "source"
sideset 2 add surface in coil_gaps with y_coord < -0.001
sideset 2 name "sink"
sideset 3 add surface in volume {wp_vid}
sideset 3 name "sibc"

# 9. Kelvin 外部ドメイン (.py 側で追加; .jou 単独では生成不可)
# → .py 版使用推奨:
#    cubit.cmd("sideset...")
#    info = add_kelvin_cubit(R=0.060, symmetry=["z"])
#    cubit.cmd('radia_export netgen "..." order 1 overwrite')
```

**Kelvin 外部ドメイン**: `.jou` だけでは生成できない。`.py` に書いて `add_kelvin_cubit()` を呼ぶ。[ih_closed_torus.py](../../src/radia/panels/samples/ih_closed_torus.py) 参照。

## ラベル設計の原則

### 必ず使うべきラベル
| ラベル | 種類 | 意味 | どの method が読む |
|---|---|---|---|
| `coil` | material (block) | コイル体積 | FEM A-V |
| `air` | material (block) | 空気 | FEM A-V |
| `kelvin` | material (block) | Kelvin 外部 (periodic) | FEM A-V |
| `sibc` | sideset | ワーク表面 | PEEC+BEM + FEM A-V |
| `source` | sideset | コイル gap 入力面 | FEM A-V |
| `sink` | sideset | コイル gap 出力面 | FEM A-V |

### 名前は sidebar 管理
- 小文字 + アンダースコア (`coil_surface`, `wp_top` 等)  
- 日本語 (漢字、カタカナ) は使わない — NGSolve の ReadGmsh で問題になる場合あり

### よく間違えるパターン
- **wp 材料ブロック (`workpiece`) を作ってしまう**: BEM は wp を穴 (subtract) とし sibc で表現するのが推奨。ブロックにすると BEM の surface 抽出で面倒
- **sibc が 1 surface だけ**: 球への webcut で air が 2 分割されると sibc FD も複数。calc_peec_bem は name で複数 FD 束ねて処理する (2026-04-19 の fix)
- **gap = 0 (closed torus)**: FEM A-V 走らない (topological に source/sink が定義できない)。最低 3-5 度の gap が必要
- **source/sink y_coord filter の threshold が小さすぎ**: `y_coord > 0` だと gap 面両方とも片側判定されることがある。`y_coord > -0.001` 等で明示的に分ける

## メッシュ品質の目安

### PEEC+BEM (1-way)
- wp sibc 表面 tri: 500-2000 枚 (5% 精度)、5000+ は assembly が O(N²) で重い
- コイル体積: 粗くて OK (3 mm 前後)。BEM 1-way 自体は コイル体積 mesh 読まない

### FEM A-V
- コイル表面: ≤ δ/3 (Cu 7 kHz で 0.26 mm)、サンプル は 0.16 mm  
- コイル体積: ≤ δ (≤ 0.5 mm) 理想、現状サンプルは 0.5 mm で P_coil +75%
- ワーク sibc 表面: 粗くて OK (SIBC なので 1-2 mm)
- Kelvin 内/外: coil の 3 倍半径くらい

## 一度書いた .jou を他の条件で再利用

.jou は APREPRO 変数使えるので条件 sweep に便利:

```bash
#{R_coil = 0.030}    # major radius (m)
#{a_coil = 0.003}    # minor radius
create vertex {R_coil - a_coil} 0 0
create vertex {R_coil} 0 0
# ...
```

ただし **変数変更しても panel は気づかない** — .jou 再 play → .vol 再 export → panel で再読込が必要。

## 参考 sample ファイル

| File | 用途 |
|---|---|
| `src/radia/panels/samples/ih_peec_bem_coarse.{jou,vol,step}` | PEEC+BEM 最小 |
| `src/radia/panels/samples/ih_fem_kelvin_skin_fine.{jou,vol,step}` | FEM A-V 最小 (gapped + fine coil) |
| `src/radia/panels/samples/ih_closed_torus.{py,jou,vol,step}` | Reference (closed、現 panel 非対応だが研究用) |
| `src/radia/panels/samples/ih_fem_kelvin_sample.{jou,vol}` | 軽量 FEM 用サンプル |
