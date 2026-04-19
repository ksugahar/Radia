# Method 選択ガイド

ユーザ (研究室メンバ、学生) から「IH 解析したい」と言われた時のフロー。

## クイック判断表

| 欲しい物 | 解析タイプ | 推奨 method |
|---|---|---|
| ワーク加熱時間の概算 | 熱源としての P_wp のみ | **PEEC+BEM (1-way)** |
| ワーク加熱の設計検討 (複数条件 sweep) | P_wp を繰り返し評価 | **PEEC+BEM (1-way)** |
| コイルのインダクタンス L を knowledge paper に書く | L (wp 込み) と P の両方 | **FEM A-V** |
| 電源共振回路の設計 | L + R_total (Q 値) | **FEM A-V** |
| 鋼材ヒステリシス込み (ESIM) | 非線形 wp | **FEM A-V** (ESIM 版)、ただし現 panel は linear のみ 対応 |
| 時間ドメイン、過渡応答 | 時変解析 | 本 panel 対象外 (参考: examples/induction_heating/) |

## 判断フロー

```
┌─────────────────────────────────────────────────┐
│ ユーザ: IH 解析したい                            │
└─────────────┬───────────────────────────────────┘
              │
     ┌────────▼──────────────┐
     │ L (インダクタンス)     │
     │ が必要か?              │
     └──┬──────────┬──────────┘
        │ Yes      │ No
        │          │
        ▼          ▼
   ┌────────┐  ┌──────────────────┐
   │FEM A-V │  │PEEC+BEM 1-way    │
   │(7-9分) │  │(3分、fastest)    │
   └────────┘  └──────────────────┘
```

## ユーザ問診で聞くべきこと

IH 解析依頼を受けたら以下を確認:

1. **コイル形状**: 丸形 (torus / circular) か、racetrack か、その他?  
   → 今の panel は **円環状 (torus)** に最適化。他の形状は extend が必要。

2. **動作周波数**: kHz? MHz?  
   → Cu で f < 100 kHz なら SIBC 近似 OK。  
   → f > 1 MHz は表面集中が極端、BEM の scalar BIE が場合により適さない。

3. **ワーク材質**: Cu, Al (線形)? 鋼材 (非線形 BH)?  
   → Linear (Cu/Al) → どちらの method でも可  
   → Nonlinear (鋼材) → **FEM A-V + ESIM (Karl iteration)** 必須。panel は linear SIBC のみ現状対応

4. **知りたいのは熱源だけ? それとも電気回路パラメタ?**  
   → 熱源のみ → **PEEC+BEM** (fastest)  
   → 回路パラメタ込み → **FEM A-V**

5. **ワークは動く?** (コイルに対して相対運動 / 回転?)  
   → 動く → PEEC+BEM (wp mesh 小で、毎ステップ再 solve 可能)  
   → 静止 → どちらでも

6. **必要な精度**: engineering 5%? paper-publishable 1%?  
   → 5%: PEEC+BEM で十分  
   → 1%: FEM A-V + fine coil mesh

## Method の制約に注意

### PEEC+BEM (1-way) の前提
- コイルは **gapped torus (or similar)** — closed topology は scalar BIE で多価 → NG  
- 5% P 精度 OK なら wp mesh は粗くて良い (1-2 mm)  
- Back-reaction は 0 次近似 (wp の wp-induced 場が coil 電流分布を変えない仮定)

### FEM A-V の前提
- **gapped torus** + **source/sink port sideset** が必須  
- コイル体積 mesh (surface ≤ δ/3 理想、IH 用 sample は 0.16mm)  
- Kelvin 外部 (periodic-identified 2-sphere) で open boundary  
- Pardiso solve → メモリ 8GB 以上推奨 (ndof 874k で 15GB 程度)

## 周波数で変わる留意

| 周波数 | Cu skin depth δ | 注意 |
|---|---|---|
| 100 Hz | 6.6 mm | ほぼ DC、SIBC 近似厳しい。FEM 体積で直接解くほうが正確 |
| 1 kHz | 2.1 mm | SIBC 境界、mesh 注意 |
| 7 kHz (本 panel で validation 済) | 0.79 mm | SIBC OK、両 method でテスト済み |
| 50 kHz | 0.29 mm | SIBC よく効く、δ << 直径 |
| 500 kHz | 93 μm | SIBC ほぼ完璧、ただし mesh 細分必要 |

## 次の Method を追加する時の判断 (開発者向け)

新しい physics が必要になった時:

1. **まず examples/ で prototype** (数学・物理 validate)  
2. **tests/panels/ に golden test 追加** (2D ref or 解析解比較)  
3. **MCP `ih_knowledge` に知識追記** (formulation, gotchas)  
4. **panel に method entry 追加** (radia_ih.py の combo + `_build_XXX_command`)  
5. **golden が green のまま CI に載せる**

絶対やってはいけない:
- 目視で 5% 一致した → そのまま panel に入れる (必ず regression する)
- demo 動いた → golden なしで production 移植 (移植中に subtle bug 入る、今回の教訓)
