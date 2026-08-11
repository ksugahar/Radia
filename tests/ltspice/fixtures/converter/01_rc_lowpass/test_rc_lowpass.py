#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
テスト1: RCローパスフィルタ
schemdraw PDF と LTSpice ASC を生成
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

import schemdraw
import schemdraw.elements as elm
from schemdraw_to_ltspice import SchemdrawToLTSpice

# 回路パラメータ
R_VALUE = 1e3   # 1kΩ
C_VALUE = 1e-6  # 1uF

# ===== schemdraw で PDF 生成 =====
with schemdraw.Drawing(show=False) as d:
    d.config(unit=3, font='Times New Roman')

    # 電流源（下から上へ）
    I1 = d.add(elm.SourceI().up().label('I1\nAC 1', loc='left'))

    # 上のライン：抵抗
    R1 = d.add(elm.Resistor().right().label('R1\n1k'))

    # 出力ノード
    d.add(elm.Dot())
    d.add(elm.Label().label('out', loc='right'))

    # キャパシタ（下へ）
    C1 = d.add(elm.Capacitor().down().label('C1\n1u', loc='bottom'))

    # 下のライン：GNDへ戻る
    d.add(elm.Line().left().to(I1.start))
    d.add(elm.Ground())

    d.save('test_rc_lowpass.pdf')
    print('PDF saved: test_rc_lowpass.pdf')

# ===== LTSpice ASC 生成 =====
converter = SchemdrawToLTSpice()

# 回路構成（schemdraw座標系）
# 電流源: (0,0) -> (0,3)
# 抵抗: (0,3) -> (3,3)
# キャパシタ: (3,3) -> (3,0)
# ワイヤ: (3,0) -> (0,0)
# GND: (0,0)

converter.add_current_source('I1', 'AC 1', (0, 0), (0, 3))
converter.add_resistor('R1', R_VALUE, (0, 3), (3, 3))
converter.add_capacitor('C1', C_VALUE, (3, 3), (3, 0))
converter.add_wire((3, 0), (0, 0))
converter.add_ground((0, 0))
converter.add_label((3, 3), 'out', is_output=True)
converter.add_spice_directive(0, 512, ".ac dec 100 1 100k")

converter.save_asc('test_rc_lowpass.asc', sheet_height=680)
print('ASC saved: test_rc_lowpass.asc')
