#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
テスト7: 水平キャパシタテスト
キャパシタを水平方向に配置したテスト回路
schemdraw PDF と LTSpice ASC を生成
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

import schemdraw
import schemdraw.elements as elm
from schemdraw_to_ltspice import SchemdrawToLTSpice

# 回路パラメータ
C1_VALUE = 1e-6    # 1uF
C2_VALUE = 100e-9  # 100nF

# ===== schemdraw で PDF 生成 =====
with schemdraw.Drawing(show=False) as d:
    d.config(unit=3, font='Times New Roman')

    # 電流源（下から上へ）
    I1 = d.add(elm.SourceI().up().label('I1\nAC 1', loc='left'))

    # キャパシタC1（右へ）水平配置
    C1 = d.add(elm.Capacitor().right().label('C1\n1u'))

    # 中間ノード
    d.add(elm.Dot())

    # キャパシタC2（右へ）水平配置
    C2 = d.add(elm.Capacitor().right().label('C2\n100n'))

    # 出力ノード
    d.add(elm.Dot())
    d.add(elm.Label().label('out', loc='right'))

    # 負荷抵抗（下へ）
    RL = d.add(elm.Resistor().down().label('RL\n50', loc='right'))

    # 下のライン：GNDへ戻る
    d.add(elm.Line().left().to(I1.start))
    d.add(elm.Ground())

    d.save('test_horizontal_cap.pdf')
    print('PDF saved: test_horizontal_cap.pdf')

# ===== LTSpice ASC 生成 =====
converter = SchemdrawToLTSpice()

# 回路構成（schemdraw座標系）
# 電流源: (0,0) -> (0,3)
# キャパシタC1: (0,3) -> (3,3) [水平]
# キャパシタC2: (3,3) -> (6,3) [水平]
# 負荷抵抗RL: (6,3) -> (6,0)
# ワイヤ: (6,0) -> (0,0)
# GND: (0,0)

converter.add_current_source('I1', 'AC 1', (0, 0), (0, 3))
converter.add_capacitor('C1', C1_VALUE, (0, 3), (3, 3))
converter.add_capacitor('C2', C2_VALUE, (3, 3), (6, 3))
converter.add_resistor('RL', 50, (6, 3), (6, 0))
converter.add_wire((6, 0), (0, 0))
converter.add_ground((0, 0))
converter.add_label((6, 3), 'out', is_output=True)
converter.add_spice_directive(0, 512, ".ac dec 100 1k 1Meg")

converter.save_asc('test_horizontal_cap.asc', sheet_width=1800, sheet_height=680)
print('ASC saved: test_horizontal_cap.asc')
