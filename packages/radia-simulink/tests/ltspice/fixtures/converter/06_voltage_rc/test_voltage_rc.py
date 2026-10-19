#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
テスト4: 電圧源付きRC回路
schemdraw PDF と LTSpice ASC を生成
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

import schemdraw
import schemdraw.elements as elm
from schemdraw_to_ltspice import SchemdrawToLTSpice

# 回路パラメータ
R_VALUE = 10e3   # 10kΩ
C_VALUE = 10e-9  # 10nF

# ===== schemdraw で PDF 生成 =====
with schemdraw.Drawing(show=False) as d:
    d.config(unit=3, font='Times New Roman')

    # 電圧源（下から上へ）
    V1 = d.add(elm.SourceV().up().label('V1\nAC 1', loc='left'))

    # 上のライン：抵抗
    R1 = d.add(elm.Resistor().right().label('R1\n10k'))

    # 出力ノード
    d.add(elm.Dot())
    d.add(elm.Label().label('out', loc='right'))

    # キャパシタ（下へ）
    C1 = d.add(elm.Capacitor().down().label('C1\n10n', loc='bottom'))

    # 下のライン：GNDへ戻る
    d.add(elm.Line().left().to(V1.start))
    d.add(elm.Ground())

    d.save('test_voltage_rc.pdf')
    print('PDF saved: test_voltage_rc.pdf')

# ===== LTSpice ASC 生成 =====
converter = SchemdrawToLTSpice()

# 回路構成（schemdraw座標系）
converter.add_voltage_source('V1', 'AC 1', (0, 0), (0, 3))
converter.add_resistor('R1', R_VALUE, (0, 3), (3, 3))
converter.add_capacitor('C1', C_VALUE, (3, 3), (3, 0))
converter.add_wire((3, 0), (0, 0))
converter.add_ground((0, 0))
converter.add_label((3, 3), 'out', is_output=True)
converter.add_spice_directive(0, 512, ".ac dec 100 100 1Meg")

converter.save_asc('test_voltage_rc.asc', sheet_height=680)
print('ASC saved: test_voltage_rc.asc')
