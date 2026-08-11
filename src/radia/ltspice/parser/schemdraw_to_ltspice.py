#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
schemdraw to LTSpice ASC converter

schemdrawで作成した受動素子回路をLTSpice ASCファイルに変換する汎用モジュール

LTSpiceシンボルのアンカーポイント（実測値）:
- R90抵抗: sym(x, y) → 右端子(x, y+16), 左端子(x-64, y+16)
- R0抵抗: sym(x, y) → 上端子(x+16, y), 下端子(x+16, y+64)
- R0インダクタ: sym(x, y) → 上端子(x+16, y-16), 下端子(x+16, y+96)
- R90インダクタ: sym(x, y) → 左端子(x, y+16), 右端子(x+112, y+16)
- R0キャパシタ: sym(x, y) → 上端子(x+16, y), 下端子(x+16, y+64)
- R90キャパシタ: sym(x, y) → 左端子(x, y+16), 右端子(x+64, y+16)
- R180電流源: sym(x, y) → 上端子(x, y-64), 下端子(x, y+64)
- R0電圧源: sym(x, y) → 上端子(x, y), 下端子(x, y+96)
"""

import schemdraw
import schemdraw.elements as elm
import numpy as np
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass


@dataclass
class LTSpiceSymbol:
    """LTSpiceシンボル情報"""
    symbol_type: str  # res, ind, cap, current, voltage
    x: int
    y: int
    rotation: str  # R0, R90, R180, R270
    name: str
    value: str
    terminals: Tuple[Tuple[int, int], Tuple[int, int]]  # (端子1, 端子2)


@dataclass
class LTSpiceWire:
    """LTSpiceワイヤ情報"""
    x1: int
    y1: int
    x2: int
    y2: int


class SchemdrawToLTSpice:
    """schemdrawからLTSpiceへの変換クラス"""

    # LTSpiceグリッドサイズ
    GRID_SIZE = 16

    # schemdrawからLTSpiceへのスケール係数
    SCALE = 64  # schemdraw 1unit = 64px in LTSpice

    # シンボルサイズ（LTSpice座標）
    SYMBOL_SIZES = {
        'res': {'R0': (32, 64), 'R90': (64, 32)},
        'ind': {'R0': (32, 112), 'R90': (112, 32)},
        'cap': {'R0': (32, 64), 'R90': (64, 32)},
        'current': {'R0': (32, 96), 'R180': (32, 96)},
        'voltage': {'R0': (32, 96), 'R180': (32, 96)},
    }

    # アンカーポイントオフセット（シンボル座標から各端子への相対位置）
    # 形式: (端子1オフセット, 端子2オフセット)
    ANCHOR_OFFSETS = {
        'res': {
            'R0': ((16, 0), (16, 64)),      # 上端子, 下端子
            'R90': ((-64, 16), (0, 16)),    # 左端子, 右端子（アンカーは右端子）
            'R180': ((16, 64), (16, 0)),
            'R270': ((64, 16), (0, 16)),
        },
        'ind': {
            'R0': ((16, -16), (16, 96)),    # 上端子, 下端子
            'R90': ((-96, 16), (0, 16)),    # 左端子, 右端子（アンカーは右端子）
            'R180': ((16, 112), (16, 0)),
            'R270': ((96, 16), (0, 16)),
        },
        'cap': {
            'R0': ((16, 0), (16, 64)),
            'R90': ((-64, 16), (0, 16)),    # 左端子, 右端子（アンカーは右端子）
            'R180': ((16, 64), (16, 0)),
            'R270': ((64, 16), (0, 16)),
        },
        'current': {
            'R0': ((16, 0), (16, 96)),
            'R180': ((0, -64), (0, 64)),    # 上端子, 下端子
        },
        'voltage': {
            'R0': ((0, 0), (0, 96)),        # 上端子(+), 下端子(-)
            'R180': ((0, 96), (0, 0)),
        },
    }

    def __init__(self, scale: float = 64, y_offset: int = 192):
        """
        Args:
            scale: schemdraw座標からLTSpice座標へのスケール係数
            y_offset: Y座標のオフセット（LTSpiceでは上が0）
        """
        self.scale = scale
        self.y_offset = y_offset
        self.symbols: List[LTSpiceSymbol] = []
        self.wires: List[LTSpiceWire] = []
        self.flags: List[Tuple[int, int, str]] = []
        self.iopins: List[Tuple[int, int, str]] = []
        self.texts: List[Tuple[int, int, str]] = []

    def snap_to_grid(self, x: float, y: float) -> Tuple[int, int]:
        """座標をLTSpiceグリッドにスナップ"""
        gx = round(x / self.GRID_SIZE) * self.GRID_SIZE
        gy = round(y / self.GRID_SIZE) * self.GRID_SIZE
        return int(gx), int(gy)

    def schemdraw_to_ltspice_coord(self, sx: float, sy: float) -> Tuple[int, int]:
        """schemdraw座標をLTSpice座標に変換

        schemdraw: Y軸上向き
        LTSpice: Y軸下向き
        """
        lx = sx * self.scale
        ly = self.y_offset - sy * self.scale + self.y_offset
        return self.snap_to_grid(lx, ly)

    def get_rotation_from_direction(self, start: Tuple[float, float],
                                    end: Tuple[float, float]) -> str:
        """schemdrawの方向からLTSpice回転を決定"""
        dx = end[0] - start[0]
        dy = end[1] - start[1]

        if abs(dx) > abs(dy):
            # 水平方向
            return 'R90' if dx > 0 else 'R270'
        else:
            # 垂直方向
            return 'R0' if dy < 0 else 'R180'

    def calculate_symbol_position(self, symbol_type: str, rotation: str,
                                  terminal1: Tuple[int, int],
                                  terminal2: Tuple[int, int]) -> Tuple[int, int]:
        """端子位置からシンボル配置位置を逆算"""
        offsets = self.ANCHOR_OFFSETS.get(symbol_type, {}).get(rotation)
        if offsets is None:
            return terminal1

        off1, off2 = offsets

        # R90の場合、アンカーは右端子（端子2）
        if rotation == 'R90':
            sym_x = terminal2[0] - off2[0]
            sym_y = terminal2[1] - off2[1]
        else:
            sym_x = terminal1[0] - off1[0]
            sym_y = terminal1[1] - off1[1]

        return self.snap_to_grid(sym_x, sym_y)

    def get_terminal_positions(self, symbol_type: str, rotation: str,
                               sym_x: int, sym_y: int) -> Tuple[Tuple[int, int], Tuple[int, int]]:
        """シンボル位置から端子位置を計算"""
        offsets = self.ANCHOR_OFFSETS.get(symbol_type, {}).get(rotation)
        if offsets is None:
            return ((sym_x, sym_y), (sym_x, sym_y))

        off1, off2 = offsets
        t1 = (sym_x + off1[0], sym_y + off1[1])
        t2 = (sym_x + off2[0], sym_y + off2[1])
        return t1, t2

    def format_value(self, value: float) -> str:
        """値をLTSpice形式にフォーマット"""
        if value == 0:
            return "0"
        elif value >= 1e6:
            return f"{value/1e6:.6g}Meg"
        elif value >= 1e3:
            return f"{value/1e3:.6g}k"
        elif value >= 1:
            return f"{value:.6g}"
        elif value >= 1e-3:
            return f"{value*1e3:.6g}m"
        elif value >= 1e-6:
            return f"{value*1e6:.6g}u"
        elif value >= 1e-9:
            return f"{value*1e9:.6g}n"
        elif value >= 1e-12:
            return f"{value*1e12:.6g}p"
        elif value >= 1e-15:
            return f"{value*1e15:.6g}f"
        else:
            return f"{value:.6g}"

    def add_resistor(self, name: str, value: float,
                     start: Tuple[float, float], end: Tuple[float, float]):
        """抵抗を追加"""
        t1 = self.schemdraw_to_ltspice_coord(start[0], start[1])
        t2 = self.schemdraw_to_ltspice_coord(end[0], end[1])
        rotation = self.get_rotation_from_direction(start, end)

        sym_pos = self.calculate_symbol_position('res', rotation, t1, t2)
        terminals = self.get_terminal_positions('res', rotation, sym_pos[0], sym_pos[1])

        self.symbols.append(LTSpiceSymbol(
            symbol_type='res',
            x=sym_pos[0],
            y=sym_pos[1],
            rotation=rotation,
            name=name,
            value=self.format_value(value),
            terminals=terminals
        ))

        # ワイヤを追加（端子間）
        self.wires.append(LTSpiceWire(t1[0], t1[1], terminals[0][0], terminals[0][1]))
        self.wires.append(LTSpiceWire(terminals[0][0], terminals[0][1],
                                       terminals[1][0], terminals[1][1]))
        self.wires.append(LTSpiceWire(terminals[1][0], terminals[1][1], t2[0], t2[1]))

    def add_inductor(self, name: str, value: float,
                     start: Tuple[float, float], end: Tuple[float, float]):
        """インダクタを追加"""
        t1 = self.schemdraw_to_ltspice_coord(start[0], start[1])
        t2 = self.schemdraw_to_ltspice_coord(end[0], end[1])
        rotation = self.get_rotation_from_direction(start, end)

        sym_pos = self.calculate_symbol_position('ind', rotation, t1, t2)
        terminals = self.get_terminal_positions('ind', rotation, sym_pos[0], sym_pos[1])

        self.symbols.append(LTSpiceSymbol(
            symbol_type='ind',
            x=sym_pos[0],
            y=sym_pos[1],
            rotation=rotation,
            name=name,
            value=self.format_value(value),
            terminals=terminals
        ))

        # ワイヤを追加
        self.wires.append(LTSpiceWire(t1[0], t1[1], terminals[0][0], terminals[0][1]))
        self.wires.append(LTSpiceWire(terminals[0][0], terminals[0][1],
                                       terminals[1][0], terminals[1][1]))
        self.wires.append(LTSpiceWire(terminals[1][0], terminals[1][1], t2[0], t2[1]))

    def add_capacitor(self, name: str, value: float,
                      start: Tuple[float, float], end: Tuple[float, float]):
        """キャパシタを追加"""
        t1 = self.schemdraw_to_ltspice_coord(start[0], start[1])
        t2 = self.schemdraw_to_ltspice_coord(end[0], end[1])
        rotation = self.get_rotation_from_direction(start, end)

        sym_pos = self.calculate_symbol_position('cap', rotation, t1, t2)
        terminals = self.get_terminal_positions('cap', rotation, sym_pos[0], sym_pos[1])

        self.symbols.append(LTSpiceSymbol(
            symbol_type='cap',
            x=sym_pos[0],
            y=sym_pos[1],
            rotation=rotation,
            name=name,
            value=self.format_value(value),
            terminals=terminals
        ))

        # ワイヤを追加
        self.wires.append(LTSpiceWire(t1[0], t1[1], terminals[0][0], terminals[0][1]))
        self.wires.append(LTSpiceWire(terminals[0][0], terminals[0][1],
                                       terminals[1][0], terminals[1][1]))
        self.wires.append(LTSpiceWire(terminals[1][0], terminals[1][1], t2[0], t2[1]))

    def add_current_source(self, name: str, value: str,
                           start: Tuple[float, float], end: Tuple[float, float]):
        """電流源を追加

        schemdrawのSourceI: startが-端子（下）、endが+端子（上）
        R180電流源: 上端子(t_top)、下端子(t_bottom)
        """
        t_bottom = self.schemdraw_to_ltspice_coord(start[0], start[1])
        t_top = self.schemdraw_to_ltspice_coord(end[0], end[1])

        # 電流源は通常R180で配置（矢印が上向き）
        rotation = 'R180'

        # R180: off1=上端子, off2=下端子
        # シンボル位置は下端子から計算（off2を使う）
        offsets = self.ANCHOR_OFFSETS['current']['R180']
        off2 = offsets[1]  # 下端子オフセット (0, 64)
        sym_x = t_bottom[0] - off2[0]
        sym_y = t_bottom[1] - off2[1]
        sym_pos = self.snap_to_grid(sym_x, sym_y)

        terminals = (t_top, t_bottom)

        self.symbols.append(LTSpiceSymbol(
            symbol_type='current',
            x=sym_pos[0],
            y=sym_pos[1],
            rotation=rotation,
            name=name,
            value=value,
            terminals=terminals
        ))

        # ワイヤを追加（上端子から下端子）
        self.wires.append(LTSpiceWire(t_top[0], t_top[1], t_bottom[0], t_bottom[1]))

    def add_voltage_source(self, name: str, value: str,
                           start: Tuple[float, float], end: Tuple[float, float]):
        """電圧源を追加

        schemdrawのSourceV: startが-端子（下）、endが+端子（上）
        R0電圧源: 上端子(t_top)、下端子(t_bottom)
        """
        t_bottom = self.schemdraw_to_ltspice_coord(start[0], start[1])
        t_top = self.schemdraw_to_ltspice_coord(end[0], end[1])

        rotation = 'R0'

        # R0: off1=上端子(16, 0), off2=下端子(16, 96)
        # シンボル位置は上端子から計算（off1を使う）
        offsets = self.ANCHOR_OFFSETS['voltage']['R0']
        off1 = offsets[0]  # 上端子オフセット (16, 0)
        sym_x = t_top[0] - off1[0]
        sym_y = t_top[1] - off1[1]
        sym_pos = self.snap_to_grid(sym_x, sym_y)

        terminals = (t_top, t_bottom)

        self.symbols.append(LTSpiceSymbol(
            symbol_type='voltage',
            x=sym_pos[0],
            y=sym_pos[1],
            rotation=rotation,
            name=name,
            value=value,
            terminals=terminals
        ))

        # ワイヤを追加（上端子から下端子）
        self.wires.append(LTSpiceWire(t_top[0], t_top[1], t_bottom[0], t_bottom[1]))

    def add_wire(self, start: Tuple[float, float], end: Tuple[float, float]):
        """ワイヤを追加"""
        t1 = self.schemdraw_to_ltspice_coord(start[0], start[1])
        t2 = self.schemdraw_to_ltspice_coord(end[0], end[1])
        self.wires.append(LTSpiceWire(t1[0], t1[1], t2[0], t2[1]))

    def add_ground(self, pos: Tuple[float, float]):
        """GNDフラグを追加"""
        lt_pos = self.schemdraw_to_ltspice_coord(pos[0], pos[1])
        self.flags.append((lt_pos[0], lt_pos[1], '0'))

    def add_label(self, pos: Tuple[float, float], name: str, is_output: bool = False):
        """ラベル/出力ピンを追加"""
        lt_pos = self.schemdraw_to_ltspice_coord(pos[0], pos[1])
        self.flags.append((lt_pos[0], lt_pos[1], name))
        if is_output:
            self.iopins.append((lt_pos[0], lt_pos[1], 'Out'))

    def add_spice_directive(self, x: int, y: int, directive: str):
        """SPICEディレクティブを追加"""
        self.texts.append((x, y, f"!{directive}"))

    def generate_asc(self, sheet_width: int = 1600, sheet_height: int = 400) -> str:
        """ASCファイル内容を生成"""
        lines = []
        lines.append("Version 4")
        lines.append(f"SHEET 1 {sheet_width} {sheet_height}")

        # ワイヤを追加（重複と長さゼロを除去）
        unique_wires = set()
        for w in self.wires:
            if w.x1 != w.x2 or w.y1 != w.y2:
                wire_tuple = (min(w.x1, w.x2), min(w.y1, w.y2),
                              max(w.x1, w.x2), max(w.y1, w.y2))
                if (w.x1, w.y1, w.x2, w.y2) not in unique_wires:
                    unique_wires.add((w.x1, w.y1, w.x2, w.y2))
                    lines.append(f"WIRE {w.x1} {w.y1} {w.x2} {w.y2}")

        # フラグを追加
        for x, y, name in self.flags:
            lines.append(f"FLAG {x} {y} {name}")

        # IOPINを追加
        for x, y, direction in self.iopins:
            lines.append(f"IOPIN {x} {y} {direction}")

        # シンボルを追加
        for sym in self.symbols:
            lines.append(f"SYMBOL {sym.symbol_type} {sym.x} {sym.y} {sym.rotation}")

            # シンボル固有のWINDOW設定
            if sym.symbol_type == 'res':
                if sym.rotation == 'R90':
                    lines.append("WINDOW 0 0 56 VBottom 2")
                    lines.append("WINDOW 3 32 56 VTop 2")
            elif sym.symbol_type == 'current':
                lines.append("WINDOW 0 24 80 Left 2")
                lines.append("WINDOW 3 24 0 Left 2")

            lines.append(f"SYMATTR InstName {sym.name}")
            lines.append(f"SYMATTR Value {sym.value}")

        # テキスト/ディレクティブを追加
        for x, y, text in self.texts:
            lines.append(f"TEXT {x} {y} Left 2 {text}")

        return "\n".join(lines)

    def save_asc(self, filename: str, sheet_width: int = 1600, sheet_height: int = 400):
        """ASCファイルを保存"""
        content = self.generate_asc(sheet_width, sheet_height)
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"ASCファイルを保存しました: {filename}")


def convert_schemdraw_element(converter: SchemdrawToLTSpice,
                               element: Any,
                               name: str,
                               value: Any):
    """schemdraw要素をLTSpiceに変換"""
    start = (element.absanchors['start'][0], element.absanchors['start'][1])
    end = (element.absanchors['end'][0], element.absanchors['end'][1])

    element_type = type(element).__name__

    if element_type == 'Resistor':
        converter.add_resistor(name, float(value), start, end)
    elif element_type in ['Inductor', 'Inductor2']:
        converter.add_inductor(name, float(value), start, end)
    elif element_type in ['Capacitor', 'Capacitor2']:
        converter.add_capacitor(name, float(value), start, end)
    elif element_type == 'SourceI':
        converter.add_current_source(name, str(value), start, end)
    elif element_type == 'SourceV':
        converter.add_voltage_source(name, str(value), start, end)
    elif element_type == 'Line':
        converter.add_wire(start, end)


# テスト用
if __name__ == '__main__':
    # 簡単なRC回路をテスト
    converter = SchemdrawToLTSpice()

    # 手動でテスト回路を追加
    converter.add_current_source('I1', 'AC 1', (0, 0), (0, 3))
    converter.add_resistor('R1', 100, (0, 3), (3, 3))
    converter.add_inductor('L1', 1e-9, (3, 3), (3, 0))
    converter.add_wire((0, 0), (3, 0))
    converter.add_ground((0, 0))
    converter.add_label((3, 3), 'out', is_output=True)
    converter.add_spice_directive(0, 512, ".ac dec 100 1 1Meg")

    # ASCファイルを生成
    converter.save_asc('LTSpice/test_converter.asc', sheet_height=680)
    print("\nGenerated ASC content:")
    print(converter.generate_asc())
