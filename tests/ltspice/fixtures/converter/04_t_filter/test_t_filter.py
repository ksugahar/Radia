#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
テスト5: T型フィルタ（L-C-L）
schemdraw PDF と LTSpice ASC を生成
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

import schemdraw
import schemdraw.elements as elm
from schemdraw_to_ltspice import SchemdrawToLTSpice

# 回路パラメータ
L1_VALUE = 1e-3   # 1mH
C_VALUE = 100e-9  # 100nF
L2_VALUE = 1e-3   # 1mH

# ===== schemdraw で PDF 生成 =====
with schemdraw.Drawing(show=False) as d:
    d.config(unit=3, font='Times New Roman')

    # 電流源（下から上へ）
    I1 = d.add(elm.SourceI().up().label('I1\nAC 1', loc='left'))

    # 上のライン：インダクタ1
    L1 = d.add(elm.Inductor2().right().label('L1\n1m'))

    # 中間ノード
    mid = d.add(elm.Dot())

    # キャパシタ（下へ）
    d.push()
    C1 = d.add(elm.Capacitor().down().label('C1\n100n', loc='left'))
    d.pop()

    # インダクタ2（右へ）
    L2 = d.add(elm.Inductor2().right().label('L2\n1m'))

    # 出力ノード
    d.add(elm.Dot())
    d.add(elm.Label().label('out', loc='right'))

    # 負荷抵抗（下へ）
    RL = d.add(elm.Resistor().down().label('RL\n50', loc='right'))

    # 下のライン：GNDへ戻る
    d.add(elm.Line().left().to(C1.end))
    d.add(elm.Line().left().to(I1.start))
    d.add(elm.Ground())

    d.save('test_t_filter.pdf')
    print('PDF saved: test_t_filter.pdf')

# ===== LTSpice ASC 生成 =====
converter = SchemdrawToLTSpice()

# 回路構成（schemdraw座標系）
# 電流源: (0,0) -> (0,3)
# インダクタ1: (0,3) -> (3,3)
# キャパシタ: (3,3) -> (3,0)
# インダクタ2: (3,3) -> (6,3)
# 負荷抵抗: (6,3) -> (6,0)
# ワイヤ: (6,0) -> (0,0)
# GND: (0,0)

converter.add_current_source('I1', 'AC 1', (0, 0), (0, 3))
converter.add_inductor('L1', L1_VALUE, (0, 3), (3, 3))
converter.add_capacitor('C1', C_VALUE, (3, 3), (3, 0))
converter.add_inductor('L2', L2_VALUE, (3, 3), (6, 3))
converter.add_resistor('RL', 50, (6, 3), (6, 0))
converter.add_wire((6, 0), (3, 0))
converter.add_wire((3, 0), (0, 0))
converter.add_ground((0, 0))
converter.add_label((6, 3), 'out', is_output=True)
converter.add_spice_directive(0, 512, ".ac dec 100 1k 1Meg")

converter.save_asc('test_t_filter.asc', sheet_width=1800, sheet_height=680)
print('ASC saved: test_t_filter.asc')
