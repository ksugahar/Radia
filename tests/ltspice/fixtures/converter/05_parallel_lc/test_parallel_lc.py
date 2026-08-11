#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
テスト3: 並列LC回路（タンク回路）
schemdraw PDF と LTSpice ASC を生成
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

import schemdraw
import schemdraw.elements as elm
from schemdraw_to_ltspice import SchemdrawToLTSpice

# 回路パラメータ
L_VALUE = 1e-3   # 1mH
C_VALUE = 100e-9 # 100nF

# ===== schemdraw で PDF 生成 =====
with schemdraw.Drawing(show=False) as d:
    d.config(unit=3, font='Times New Roman')

    # 電流源（下から上へ）
    I1 = d.add(elm.SourceI().up().label('I1\nAC 1', loc='left'))
    start_pos = I1.start

    # 上のライン
    d.add(elm.Line().right().length(1.5))
    branch_top = d.add(elm.Dot())

    # 左側：インダクタ（下へ）
    d.push()
    L1 = d.add(elm.Inductor2().down().label('L1\n1m', loc='left'))
    branch_bottom_l = d.add(elm.Dot())
    d.pop()

    # 右側：キャパシタ（下へ）
    d.add(elm.Line().right().length(1.5))
    d.add(elm.Dot())
    d.add(elm.Label().label('out', loc='right'))
    C1 = d.add(elm.Capacitor().down().label('C1\n100n', loc='right'))

    # 下のライン
    d.add(elm.Line().left().to(branch_bottom_l.center))
    d.add(elm.Line().left().to(start_pos))
    d.add(elm.Ground())

    d.save('test_parallel_lc.pdf')
    print('PDF saved: test_parallel_lc.pdf')

# ===== LTSpice ASC 生成 =====
converter = SchemdrawToLTSpice()

# 回路構成（schemdraw座標系）
# 電流源: (0,0) -> (0,3)
# ワイヤ: (0,3) -> (1.5,3)
# インダクタ: (1.5,3) -> (1.5,0)
# ワイヤ: (1.5,3) -> (3,3)
# キャパシタ: (3,3) -> (3,0)
# ワイヤ: (3,0) -> (0,0)
# GND: (0,0)

converter.add_current_source('I1', 'AC 1', (0, 0), (0, 3))
converter.add_wire((0, 3), (1.5, 3))
converter.add_inductor('L1', L_VALUE, (1.5, 3), (1.5, 0))
converter.add_wire((1.5, 3), (3, 3))
converter.add_capacitor('C1', C_VALUE, (3, 3), (3, 0))
converter.add_wire((3, 0), (1.5, 0))
converter.add_wire((1.5, 0), (0, 0))
converter.add_ground((0, 0))
converter.add_label((3, 3), 'out', is_output=True)
converter.add_spice_directive(0, 512, ".ac dec 100 1k 100k")

converter.save_asc('test_parallel_lc.asc', sheet_height=680)
print('ASC saved: test_parallel_lc.asc')
