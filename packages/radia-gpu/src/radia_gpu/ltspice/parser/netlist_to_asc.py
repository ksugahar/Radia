#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPICE Netlist (.cir) to LTSpice ASC (.asc) Converter

SPICEネットリストを解析し、自動レイアウトを行い、
LTSpice .ascスキーマティックファイルを生成する。

参考: dominc8/netlist_converter (C++/OGDF)
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set
from enum import Enum
import re
import math


# =============================================================================
# 1. NetlistParser - ネットリスト解析
# =============================================================================

class ComponentType(Enum):
    """コンポーネント種別"""
    RESISTOR = 'R'
    CAPACITOR = 'C'
    INDUCTOR = 'L'
    VOLTAGE = 'V'
    CURRENT = 'I'
    DIODE = 'D'
    BJT = 'Q'
    MOSFET = 'M'
    JFET = 'J'
    SUBCIRCUIT = 'X'
    VCVS = 'E'
    CCCS = 'F'
    VCCS = 'G'
    CCVS = 'H'
    BEHAVIORAL = 'B'
    SWITCH = 'S'
    COUPLED = 'K'
    TLINE = 'T'
    OPAMP = 'U'
    GROUND = 'GND'
    NET_NODE = 'NET'


@dataclass
class Component:
    """パースされたコンポーネント"""
    name: str                    # R1, C1, V1 etc.
    comp_type: ComponentType
    node_pos: str                # 正端子ノード名 (C/D/out+)
    node_neg: str                # 負端子ノード名 (E/S/out-)
    value: str                   # 値（文字列のまま保持）
    node_ctrl: str = ''          # 制御端子+ (B/G/ctrl+) — 3/4端子素子用
    node_ctrl2: str = ''         # 制御端子- (ctrl-) — 4端子素子用
    node_out: str = ''           # 出力端子 — 5端子opamp用
    raw_line: str = ''           # 元のネットリスト行
    symbol_hint: str = ''        # 元の LTspice SYMBOL kind (ind2, schottky, pnp, ...)
    extra_nodes: List[str] = field(default_factory=list)  # SUBCIRCUIT 3+ pin の追加ノード


@dataclass
class SpiceDirective:
    """SPICEディレクティブ（.tran, .ac, .param 等）"""
    text: str


def _looks_like_pure_node_token(tok: str) -> bool:
    """Conservative substrate-vs-model disambiguation for BJT/JFET 4-term parse.

    LTspice's canonical netlist always emits the 4-terminal form for
    BJT/JFET, with substrate=0 (gnd) when the substrate is connected to
    ground:
        Q1 N001 N002 0 0 2N2222    → 4-term, substrate=0, model=2N2222
        Q1 N001 N002 0 2N2222      → 3-term, model=2N2222 (no substrate)

    Returning True means the token is *unambiguously* a node and not a
    model — so the parser should treat it as substrate when present
    between emitter/source and the model token.

    We only treat **pure-integer** tokens (e.g. `0`, `12`, `99`) as
    unambiguous nodes here. Named nodes like `N004` or `OUT` overlap
    too much with model-name conventions to safely auto-detect; for
    those we default to 3-terminal parsing. (Manually-written netlists
    almost never use a non-ground substrate; LTspice's canonical
    output uses substrate=0 for the cases where it matters.)
    """
    if not tok:
        return False
    try:
        int(tok)
        return True
    except ValueError:
        return False


class NetlistParser:
    """SPICEネットリスト (.cir) パーサー"""

    # コンポーネントの先頭文字 → ComponentType
    TYPE_MAP = {
        'R': ComponentType.RESISTOR,
        'C': ComponentType.CAPACITOR,
        'L': ComponentType.INDUCTOR,
        'V': ComponentType.VOLTAGE,
        'I': ComponentType.CURRENT,
        'D': ComponentType.DIODE,
        'Q': ComponentType.BJT,
        'M': ComponentType.MOSFET,
        'J': ComponentType.JFET,
        'X': ComponentType.SUBCIRCUIT,
        'E': ComponentType.VCVS,
        'F': ComponentType.CCCS,
        'G': ComponentType.VCCS,
        'H': ComponentType.CCVS,
        'B': ComponentType.BEHAVIORAL,
        'S': ComponentType.SWITCH,
        'K': ComponentType.COUPLED,
        'T': ComponentType.TLINE,
    }

    def __init__(self):
        self.components: List[Component] = []
        self.directives: List[SpiceDirective] = []
        self.title: str = ''

    def parse_file(self, filepath: str) -> 'NetlistParser':
        """ファイルからネットリストをパース"""
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        return self.parse_lines(lines)

    def parse_string(self, text: str) -> 'NetlistParser':
        """文字列からネットリストをパース"""
        lines = text.strip().split('\n')
        return self.parse_lines(lines)

    def parse_lines(self, lines: List[str]) -> 'NetlistParser':
        """行リストからネットリストをパース"""
        self.components = []
        self.directives = []
        self.title = ''
        # Lines the parser could not classify (unknown prefix or malformed).
        # Kept as (1-based line number, original text) for B4 lint reporting.
        self.unparsed_lines: List[Tuple[int, str]] = []

        # `* @sym=<kind>` hint preceding a component restores the original
        # LTspice SYMBOL kind (ind2, schottky, pnp, polcap, ...).
        pending_symbol_hint = ''

        # `.subckt ... .ends` body accumulator. While inside a subckt
        # block, EVERY line (components, models, params, comments) is
        # absorbed into a single multi-line directive so the regenerated
        # .asc can re-emit the block verbatim. This prevents subckt
        # members from leaking into top-level components.
        subckt_accumulator: Optional[List[str]] = None

        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue

            # 最初の非空行はタイトル
            if (i == 0 and subckt_accumulator is None
                    and not line.startswith('.') and not line.startswith('*')):
                first_char = line[0].upper()
                if first_char not in self.TYPE_MAP:
                    self.title = line
                    continue

            # Inside a .subckt block — absorb everything until .ends
            if subckt_accumulator is not None:
                subckt_accumulator.append(line)
                if line.lower().startswith('.ends'):
                    full_block = '\\n'.join(subckt_accumulator)
                    self.directives.append(SpiceDirective(text=full_block))
                    subckt_accumulator = None
                continue

            # Entering a .subckt block: start accumulating
            if line.lower().startswith('.subckt'):
                subckt_accumulator = [line]
                continue

            # コメント行 (@sym=... hint だけ拾う)
            if line.startswith('*'):
                if '@sym=' in line:
                    pending_symbol_hint = line.split('@sym=', 1)[1].strip()
                continue

            # インラインコメント除去
            if ';' in line:
                comment_pos = line.index(';')
                # 文字列リテラル内でないか簡易チェック
                line = line[:comment_pos].strip()
                if not line:
                    continue

            # ディレクティブ行
            if line.startswith('.'):
                directive = line.lower()
                if directive == '.end':
                    break
                self.directives.append(SpiceDirective(text=line))
                continue

            # K文（結合定数）/ A文（デジタル素子）はディレクティブとして扱う
            if line[0].upper() in ('K', 'A'):
                self.directives.append(SpiceDirective(text=line))
                continue

            # U文（opamp/サブサーキット）→ コンポーネントとしてパース
            if line[0].upper() == 'U':
                comp = self._parse_u_statement(line)
                if comp:
                    if pending_symbol_hint:
                        comp.symbol_hint = pending_symbol_hint
                    self.components.append(comp)
                else:
                    self.directives.append(SpiceDirective(text=line))
                pending_symbol_hint = ''
                continue

            # J文が3パーツ = ジャンパー（2端子短絡）→ ゼロ抵抗として扱う
            # (JFET は J name D G S model = 5パーツ以上)
            parts_check = line.split()
            if line[0].upper() == 'J' and len(parts_check) == 3:
                # J2 nodeA nodeB → R_J2 nodeA nodeB 0 (短絡)
                jname = parts_check[0]
                comp = Component(
                    name=f'R_{jname}', comp_type=ComponentType.RESISTOR,
                    node_pos=parts_check[1], node_neg=parts_check[2],
                    value='0', raw_line=line,
                )
                self.components.append(comp)
                pending_symbol_hint = ''
                continue
            # J文が4パーツ以下で3以外 = ディレクティブとして扱う
            if line[0].upper() == 'J' and len(parts_check) < 5:
                self.directives.append(SpiceDirective(text=line))
                pending_symbol_hint = ''
                continue

            # コンポーネント行をパース
            comp = self._parse_component_line(line)
            if comp:
                if pending_symbol_hint:
                    comp.symbol_hint = pending_symbol_hint
                self.components.append(comp)
            else:
                # Line did not match any known SPICE element class.
                # Remember it so the CLI can surface "what got skipped"
                # to the user instead of silently dropping it.
                self.unparsed_lines.append((i + 1, line))
            pending_symbol_hint = ''

        return self

    def _parse_component_line(self, line: str) -> Optional[Component]:
        """コンポーネント行をパース"""
        parts = line.split()
        if len(parts) < 3:
            return None

        name = parts[0]
        first_char = name[0].upper()

        if first_char not in self.TYPE_MAP:
            return None

        comp_type = self.TYPE_MAP[first_char]

        # 3端子素子: Q C B E model / M D G S B model / J D G S model
        node_ctrl = ''
        node_ctrl2 = ''
        extra_nodes: List[str] = []
        if comp_type in (ComponentType.BJT, ComponentType.JFET):
            if len(parts) >= 5:
                # SPICE BJT: Q<name> <C> <B> <E> [<sub>] <model> [<params>]
                # SPICE JFET: J<name> <D> <G> <S> [<sub>] <model> [<params>]
                # Substrate is optional. LTspice's canonical netlist always
                # includes it (`Q1 N001 N002 0 0 2N2222` has sub=0 before
                # model). Detect by token count: if there's an extra token
                # before what looks like a model name, treat it as substrate.
                node_pos = parts[1]  # Collector / Drain
                node_ctrl = parts[2] # Base / Gate
                node_neg = parts[3]  # Emitter / Source
                if len(parts) >= 6 and _looks_like_pure_node_token(parts[4]):
                    # 4-terminal form with substrate.
                    # parts[4] = substrate, parts[5..] = model + params
                    value = ' '.join(parts[5:])
                else:
                    # 3-terminal form (no substrate).
                    value = ' '.join(parts[4:])
            elif len(parts) == 4:
                # Modelless form: Q/J <name> <n1> <n2> <n3>. LTspice's npn/pnp/
                # njf/pjf symbols allow omitting SYMATTR Value, in which case the
                # extracted netlist has no trailing model token. Accept and keep
                # value empty; AscGenerator will round-trip without SYMATTR Value.
                node_pos = parts[1]
                node_ctrl = parts[2]
                node_neg = parts[3]
                value = ''
            else:
                return None
        elif comp_type == ComponentType.MOSFET:
            if len(parts) >= 6:
                # M name D G S B model
                node_pos = parts[1]  # Drain
                node_ctrl = parts[2]  # Gate
                node_neg = parts[3]  # Source
                value = parts[5] if len(parts) > 5 else ''
            elif len(parts) == 5:
                # Modelless form: M <name> <D> <G> <S> <B>. Same rationale as
                # the BJT/JFET modelless branch above.
                node_pos = parts[1]
                node_ctrl = parts[2]
                node_neg = parts[3]
                value = ''
            else:
                return None
        elif comp_type == ComponentType.COUPLED:
            # K文はコンポーネントではなくディレクティブとして扱う
            return None
        elif comp_type == ComponentType.SUBCIRCUIT:
            # X<name> <node>... <subckt_name> [param=value ...]
            # Inline parameters (LTspice opamps emit e.g.
            #   X1 in out V+ V- o level2 Avol=1Meg GBW=10Meg Rin=500Meg)
            # are NOT pins. The subckt name is the last bare token before
            # the first `key=value` token; everything before that is nodes,
            # everything from it on is parameters. The old code took the
            # very last token as the subckt name, mis-reading every
            # parameter as a node (inflating pin count and topology).
            body = parts[1:]
            pidx = next((i for i, t in enumerate(body) if '=' in t), len(body))
            pre = body[:pidx]          # nodes + subckt name (no params)
            if pidx == len(body):
                # No inline params — original semantics.
                if len(pre) >= 3:      # >=2 nodes + subckt name
                    node_pos = pre[0]
                    node_neg = pre[1]
                    extra_nodes = pre[2:-1]
                    value = pre[-1]
                elif len(pre) == 2:    # X1 N1 N2 (2 nodes, no model)
                    node_pos = pre[0]
                    node_neg = pre[1]
                    value = ''
                else:
                    return None
            else:
                # Inline params present — subckt name = last bare token.
                if len(pre) >= 2:      # >=1 node + subckt name
                    nodes = pre[:-1]
                    node_pos = nodes[0]
                    # Leave node_neg empty for a 1-pin subckt rather than
                    # fabricating a phantom ground ('0') connection.
                    node_neg = nodes[1] if len(nodes) >= 2 else ''
                    extra_nodes = nodes[2:]
                    value = pre[-1]    # subckt name (params dropped)
                else:
                    return None
        elif comp_type in (ComponentType.VCVS, ComponentType.VCCS):
            # E/G: 4端子 (out+ out- ctrl+ ctrl- gain/Laplace)
            node_pos = parts[1]
            node_neg = parts[2]
            rest = ' '.join(parts[3:])
            if 'laplace' in rest.lower() and len(parts) >= 6:
                # Laplace with control terminals: E out+ out- ctrl+ ctrl- Laplace=...
                node_ctrl = parts[3]
                node_ctrl2 = parts[4]
                value = ' '.join(parts[5:])
            elif 'laplace' in rest.lower():
                # Laplace without control: E out+ out- Laplace=...
                value = rest
            elif len(parts) >= 6:
                # 通常型: E name out+ out- ctrl+ ctrl- gain
                node_ctrl = parts[3]
                node_ctrl2 = parts[4]
                value = ' '.join(parts[5:])
            else:
                value = rest
        elif comp_type == ComponentType.SWITCH:
            # S name n+ n- ctrl+ ctrl- [model]
            if len(parts) >= 5:
                node_pos = parts[1]
                node_neg = parts[2]
                node_ctrl = parts[3]
                node_ctrl2 = parts[4]
                value = ' '.join(parts[5:]) if len(parts) > 5 else ''
            else:
                # 簡略形: S name n+ n- (制御なし)
                node_pos = parts[1]
                node_neg = parts[2]
                value = ' '.join(parts[3:])
        elif comp_type == ComponentType.TLINE:
            # T name port1+ port1- port2+ port2- params
            if len(parts) >= 5:
                node_pos = parts[1]
                node_neg = parts[2]
                node_ctrl = parts[3]
                node_ctrl2 = parts[4]
                value = ' '.join(parts[5:]) if len(parts) > 5 else ''
            else:
                node_pos = parts[1]
                node_neg = parts[2] if len(parts) > 2 else '0'
                value = ' '.join(parts[3:]) if len(parts) > 3 else ''
        else:
            # 2端子素子: name node+ node- value
            node_pos = parts[1]
            node_neg = parts[2]
            if len(parts) >= 4:
                value = ' '.join(parts[3:])
            else:
                value = ''

        return Component(
            name=name,
            comp_type=comp_type,
            node_pos=node_pos,
            node_neg=node_neg,
            value=value,
            node_ctrl=node_ctrl,
            node_ctrl2=node_ctrl2,
            raw_line=line,
            extra_nodes=extra_nodes,
        )

    def _parse_u_statement(self, line: str) -> Optional[Component]:
        """U文（opamp）をパース

        3-pin opamp.sub: U1 inv noninv out
        5-pin UniversalOpAmp: U1 in+ in- V+ V- out model
        """
        parts = line.split()
        if len(parts) < 4:
            return None

        name = parts[0]

        if len(parts) == 4:
            # 3-pin: U1 inv noninv out (no model name)
            return Component(
                name=name,
                comp_type=ComponentType.OPAMP,
                node_pos=parts[1],   # inverting input (SpiceOrder 1)
                node_neg=parts[2],   # non-inverting input (SpiceOrder 2)
                node_ctrl=parts[3],  # output (SpiceOrder 3)
                value='opamp',
                raw_line=line,
            )
        elif len(parts) >= 7:
            # 5-pin: U1 in+ in- V+ V- out model
            return Component(
                name=name,
                comp_type=ComponentType.OPAMP,
                node_pos=parts[1],    # In+ (SpiceOrder 1)
                node_neg=parts[2],    # In- (SpiceOrder 2)
                node_ctrl=parts[3],   # V+ (SpiceOrder 3)
                node_ctrl2=parts[4],  # V- (SpiceOrder 4)
                node_out=parts[5],    # OUT (SpiceOrder 5)
                value=parts[6],       # model name
                raw_line=line,
            )
        else:
            return None

    def get_all_nodes(self) -> Set[str]:
        """全ユニークノード名を取得"""
        nodes = set()
        for comp in self.components:
            nodes.add(comp.node_pos)
            nodes.add(comp.node_neg)
            if comp.node_ctrl:
                nodes.add(comp.node_ctrl)
            if comp.node_ctrl2:
                nodes.add(comp.node_ctrl2)
            if comp.node_out:
                nodes.add(comp.node_out)
        return nodes

    def get_ground_nodes(self) -> Set[str]:
        """グランドノード（'0' or 'gnd'）を取得"""
        nodes = self.get_all_nodes()
        return {n for n in nodes if n == '0' or n.lower() == 'gnd'}

    def get_signal_nodes(self) -> Set[str]:
        """信号ノード（グランド以外）を取得"""
        return self.get_all_nodes() - self.get_ground_nodes()


# =============================================================================
# 2. CircuitLayouter - グラフベースの自動レイアウト
# =============================================================================

@dataclass
class NodePosition:
    """ノードの座標"""
    x: int = 0
    y: int = 0


@dataclass
class PlacedComponent:
    """配置済みコンポーネント"""
    component: Component
    x: int = 0              # LTSpiceシンボル座標
    y: int = 0              # LTSpiceシンボル座標
    rotation: str = 'R0'    # R0, R90, R180, R270
    terminal1: Tuple[int, int] = (0, 0)  # 正端子座標 (C/D)
    terminal2: Tuple[int, int] = (0, 0)  # 負端子座標 (E/S)
    terminal3: Optional[Tuple[int, int]] = None  # 制御端子+ (B/G/ctrl+)
    terminal4: Optional[Tuple[int, int]] = None  # 制御端子- (ctrl-)
    terminal5: Optional[Tuple[int, int]] = None  # 出力端子 (opamp out)


# アンカーポイントオフセット（シンボル座標→端子座標）
# asc_parser.py の正規テーブルを使用（.asy実測値ベース）
from .asc_parser import TERMINAL_OFFSETS as _ASC_OFFSETS
from .asc_parser import TERMINAL_OFFSETS_3 as _ASC_OFFSETS_3
from .asc_parser import TERMINAL_OFFSETS_4 as _ASC_OFFSETS_4
from .asc_parser import TERMINAL_OFFSETS_5 as _ASC_OFFSETS_5

# ComponentType → asc_parser シンボル名のマッピング (2端子)
_TYPE_TO_SYM = {
    ComponentType.RESISTOR: 'res',
    ComponentType.CAPACITOR: 'cap',
    ComponentType.INDUCTOR: 'ind',
    ComponentType.VOLTAGE: 'voltage',
    ComponentType.CURRENT: 'current',
    ComponentType.DIODE: 'diode',
    ComponentType.VCVS: 'e',
    ComponentType.VCCS: 'g',
    ComponentType.CCCS: 'f',
    ComponentType.CCVS: 'h',
    ComponentType.BEHAVIORAL: 'bv',
    ComponentType.SWITCH: 'sw',
}

# 3端子素子のマッピング
_TYPE_TO_SYM_3 = {
    ComponentType.BJT: 'npn',     # デフォルト、PNPは配置時にモデル名で判定
    ComponentType.MOSFET: 'nmos',
    ComponentType.JFET: 'njf',
}

ANCHOR_OFFSETS = {}
for ct, sym_name in _TYPE_TO_SYM.items():
    if sym_name in _ASC_OFFSETS:
        ANCHOR_OFFSETS[ct] = _ASC_OFFSETS[sym_name]

# 3端子素子のオフセット（pin1=C/D, pin2=B/G, pin3=E/S）
ANCHOR_OFFSETS_3 = {}
for ct, sym_name in _TYPE_TO_SYM_3.items():
    if sym_name in _ASC_OFFSETS_3:
        ANCHOR_OFFSETS_3[ct] = _ASC_OFFSETS_3[sym_name]

# 端子間距離（pin1-pin2 の距離、R0方向）
COMPONENT_SPAN = {}
for ct, sym_name in _TYPE_TO_SYM.items():
    offs = _ASC_OFFSETS.get(sym_name, {}).get('R0')
    if offs:
        p1, p2 = offs
        COMPONENT_SPAN[ct] = abs(p2[1] - p1[1])
for ct, sym_name in _TYPE_TO_SYM_3.items():
    if sym_name in _ASC_OFFSETS_3:
        offs = _ASC_OFFSETS_3[sym_name].get('R0')
        if offs:
            COMPONENT_SPAN[ct] = abs(offs[2][1] - offs[0][1])  # C/D to E/S

GRID = 16  # LTSpiceグリッド


def snap(val: float) -> int:
    """グリッドスナップ"""
    return int(round(val / GRID) * GRID)


def calc_symbol_placement(comp_type: ComponentType, rotation: str,
                           t1: Tuple[int, int], t2: Tuple[int, int]):
    """端子座標からシンボル配置位置を逆算（2端子）"""
    offsets = ANCHOR_OFFSETS.get(comp_type, {}).get(rotation)
    if offsets is None:
        return (t1[0], t1[1], t1, t2)

    off1, off2 = offsets
    sym_x = snap(t1[0] - off1[0])
    sym_y = snap(t1[1] - off1[1])

    actual_t1 = (sym_x + off1[0], sym_y + off1[1])
    actual_t2 = (sym_x + off2[0], sym_y + off2[1])

    return (sym_x, sym_y, actual_t1, actual_t2)


def calc_symbol_placement_3t(comp_type: ComponentType, rotation: str,
                              t1: Tuple[int, int]):
    """3端子素子の配置位置を逆算

    t1 = node_pos (C/D) の目標座標。
    3端子オフセット: (off_cd, off_bg, off_es)
    Returns: (sym_x, sym_y, actual_cd, actual_bg, actual_es)
    """
    offsets = ANCHOR_OFFSETS_3.get(comp_type, {}).get(rotation)
    if offsets is None:
        return (t1[0], t1[1], t1, (t1[0], t1[1] + 48), (t1[0], t1[1] + 96))

    off_cd, off_bg, off_es = offsets
    sym_x = snap(t1[0] - off_cd[0])
    sym_y = snap(t1[1] - off_cd[1])

    actual_cd = (sym_x + off_cd[0], sym_y + off_cd[1])
    actual_bg = (sym_x + off_bg[0], sym_y + off_bg[1])
    actual_es = (sym_x + off_es[0], sym_y + off_es[1])

    return (sym_x, sym_y, actual_cd, actual_bg, actual_es)


def calc_symbol_placement_4t(comp_type: ComponentType, rotation: str,
                              t1: Tuple[int, int]):
    """4端子素子の配置位置を逆算

    t1 = 第1ピンの目標座標。
    Returns: (sym_x, sym_y, pin1, pin2, pin3, pin4)
    """
    sym_map = {
        ComponentType.VCVS: 'e', ComponentType.VCCS: 'g',
        ComponentType.SWITCH: 'sw', ComponentType.TLINE: 'tline',
    }
    sym_name = sym_map.get(comp_type)
    offsets = _ASC_OFFSETS_4.get(sym_name, {}).get(rotation) if sym_name else None
    if offsets is None:
        return (t1[0], t1[1], t1, (t1[0], t1[1]+96),
                (t1[0]-48, t1[1]+32), (t1[0]-48, t1[1]+80))

    off1, off2, off3, off4 = offsets
    sym_x = snap(t1[0] - off1[0])
    sym_y = snap(t1[1] - off1[1])

    return (sym_x, sym_y,
            (sym_x + off1[0], sym_y + off1[1]),
            (sym_x + off2[0], sym_y + off2[1]),
            (sym_x + off3[0], sym_y + off3[1]),
            (sym_x + off4[0], sym_y + off4[1]))


def calc_symbol_placement_5t(sym_name: str, rotation: str,
                              t1: Tuple[int, int]):
    """5端子素子（opamp）の配置位置を逆算

    t1 = 第1ピン(In+)の目標座標。
    Returns: (sym_x, sym_y, pin1, pin2, pin3, pin4, pin5)
    """
    offsets = _ASC_OFFSETS_5.get(sym_name, {}).get(rotation)
    if offsets is None:
        # フォールバック: UniversalOpAmp デフォルト
        offsets = _ASC_OFFSETS_5.get('universalopamp2', {}).get(rotation)
    if offsets is None:
        return (t1[0], t1[1], t1, (t1[0], t1[1]-32),
                (t1[0]+32, t1[1]-48), (t1[0]+32, t1[1]+16),
                (t1[0]+64, t1[1]-16))

    off1, off2, off3, off4, off5 = offsets
    sym_x = snap(t1[0] - off1[0])
    sym_y = snap(t1[1] - off1[1])

    return (sym_x, sym_y,
            (sym_x + off1[0], sym_y + off1[1]),
            (sym_x + off2[0], sym_y + off2[1]),
            (sym_x + off3[0], sym_y + off3[1]),
            (sym_x + off4[0], sym_y + off4[1]),
            (sym_x + off5[0], sym_y + off5[1]))


def calc_symbol_placement_3t_opamp(rotation: str, t1: Tuple[int, int]):
    """3端子opamp（opamp.sub）の配置位置を逆算

    t1 = 第1ピン(inv)の目標座標。
    Returns: (sym_x, sym_y, pin_inv, pin_noninv, pin_out)
    """
    offsets = _ASC_OFFSETS_3.get('opamp', {}).get(rotation)
    if offsets is None:
        return (t1[0], t1[1], t1, (t1[0], t1[1]+32), (t1[0]+64, t1[1]+16))

    off1, off2, off3 = offsets
    sym_x = snap(t1[0] - off1[0])
    sym_y = snap(t1[1] - off1[1])

    return (sym_x, sym_y,
            (sym_x + off1[0], sym_y + off1[1]),
            (sym_x + off2[0], sym_y + off2[1]),
            (sym_x + off3[0], sym_y + off3[1]))


class CircuitLayouter:
    """回路の自動レイアウト

    実際のLTSpice回路パターンに基づく配置:
    - ソース(V/I)は左側に縦配置（上が正端子）
    - 直列素子は上段を右へ水平配置
    - 分路素子は上段ノードから下のGNDへ縦配置
    - GNDは最下段
    """

    # ノード間隔
    H_SPACING = 192   # 水平ノード間隔
    V_SPACING = 192   # 垂直ノード間隔（上段→GND）
    TOP_Y = 192       # 上段のY座標
    GND_Y = 384       # GND段のY座標

    def __init__(self):
        self.node_positions: Dict[str, NodePosition] = {}
        self.placed_components: List[PlacedComponent] = []
        self._subckt_iso_idx = 0  # multi-pin SUBCIRCUIT isolation slots

    def layout(self, parser: NetlistParser) -> 'CircuitLayouter':
        """レイアウトを実行"""
        self.node_positions = {}
        self.placed_components = []

        ground_nodes = parser.get_ground_nodes()
        if not parser.components:
            return self

        # Step 1: コンポーネントを分類
        sources = []         # V/Iソース
        series_comps = []    # 直列素子（両端がシグナルノード）
        shunt_comps = []     # 分路素子（一端がGND）

        for comp in parser.components:
            is_source = comp.comp_type in (ComponentType.VOLTAGE,
                                            ComponentType.CURRENT)
            pos_is_gnd = comp.node_pos in ground_nodes
            neg_is_gnd = comp.node_neg in ground_nodes

            if is_source:
                sources.append(comp)
            elif pos_is_gnd or neg_is_gnd:
                shunt_comps.append(comp)
            else:
                series_comps.append(comp)

        # Step 2: ノード座標を決定
        self._assign_positions(parser, sources, series_comps, shunt_comps,
                               ground_nodes)

        # Step 3: 並列コンポーネントの水平オフセット計算
        # 同じ2ノード間に複数のコンポーネントがある場合、横にずらす
        parallel_offsets = self._calc_parallel_offsets(
            parser, sources, series_comps, shunt_comps, ground_nodes)

        # Step 4: コンポーネントを配置
        for comp in parser.components:
            offset = parallel_offsets.get(comp.name, 0)
            placed = self._place_component(comp, ground_nodes, offset)
            self.placed_components.append(placed)

        # Step 5: 重複解消 — 同一座標の部品を水平にずらす
        self._resolve_overlaps()

        return self

    def _resolve_overlaps(self):
        """部品重複と端子衝突を解消

        1. 同一座標のシンボルを水平にずらす
        2. 異なるノードの端子が同一座標にならないよう調整
        """
        for _ in range(20):
            # Phase 1: シンボル位置の重複
            occupied: Dict[Tuple[int, int], List[int]] = {}
            for i, pc in enumerate(self.placed_components):
                key = (pc.x, pc.y)
                occupied.setdefault(key, []).append(i)

            has_overlap = False
            for key, indices in occupied.items():
                if len(indices) <= 1:
                    continue
                has_overlap = True
                for rank, idx in enumerate(indices[1:], 1):
                    shift = rank * self.H_SPACING
                    self._shift_component(idx, shift)

            # Phase 2: 端子座標の衝突（異なるノードが同一座標）
            term_owner: Dict[Tuple[int, int], Tuple[str, int]] = {}  # coord -> (node, comp_idx)
            shift_set: Set[int] = set()
            for i, pc in enumerate(self.placed_components):
                if i in shift_set:
                    continue
                comp = pc.component
                term_list = [(pc.terminal1, comp.node_pos),
                             (pc.terminal2, comp.node_neg)]
                if pc.terminal3 is not None and comp.node_ctrl:
                    term_list.append((pc.terminal3, comp.node_ctrl))
                if pc.terminal4 is not None and comp.node_ctrl2:
                    term_list.append((pc.terminal4, comp.node_ctrl2))
                for term, node in term_list:
                    if term in term_owner:
                        existing_node, existing_idx = term_owner[term]
                        if existing_node != node:
                            has_overlap = True
                            shift_set.add(i)
                            break
                    else:
                        term_owner[term] = (node, i)

            # 衝突する部品を既存の全端子座標から離れた位置にシフト
            all_terms = set(term_owner.keys())
            for idx in sorted(shift_set):
                pc = self.placed_components[idx]
                # 右方向に空き位置を探す
                shift = self.H_SPACING
                while True:
                    new_t1 = (pc.terminal1[0] + shift, pc.terminal1[1])
                    new_t2 = (pc.terminal2[0] + shift, pc.terminal2[1])
                    ok = new_t1 not in all_terms and new_t2 not in all_terms
                    if ok and pc.terminal3 is not None:
                        ok = (pc.terminal3[0] + shift, pc.terminal3[1]) not in all_terms
                    if ok and pc.terminal4 is not None:
                        ok = (pc.terminal4[0] + shift, pc.terminal4[1]) not in all_terms
                    if ok:
                        break
                    shift += self.H_SPACING
                self._shift_component(idx, shift)
                # 新端子を登録
                pc2 = self.placed_components[idx]
                all_terms.add(pc2.terminal1)
                all_terms.add(pc2.terminal2)
                if pc2.terminal3 is not None:
                    all_terms.add(pc2.terminal3)
                if pc2.terminal4 is not None:
                    all_terms.add(pc2.terminal4)
                if pc2.terminal5 is not None:
                    all_terms.add(pc2.terminal5)

            if not has_overlap:
                break

    def _shift_component(self, idx: int, shift: int):
        """部品を水平にシフト"""
        pc = self.placed_components[idx]
        t3 = (pc.terminal3[0] + shift, pc.terminal3[1]) if pc.terminal3 else None
        t4 = (pc.terminal4[0] + shift, pc.terminal4[1]) if pc.terminal4 else None
        t5 = (pc.terminal5[0] + shift, pc.terminal5[1]) if pc.terminal5 else None
        self.placed_components[idx] = PlacedComponent(
            component=pc.component,
            x=pc.x + shift,
            y=pc.y,
            rotation=pc.rotation,
            terminal1=(pc.terminal1[0] + shift, pc.terminal1[1]),
            terminal2=(pc.terminal2[0] + shift, pc.terminal2[1]),
            terminal3=t3,
            terminal4=t4,
            terminal5=t5,
        )

    def _assign_positions(self, parser, sources, series_comps, shunt_comps,
                           ground_nodes):
        """ノード位置を割り当て"""
        # ソースの信号ノードを集める
        source_signal_nodes = set()
        for src in sources:
            if src.node_pos not in ground_nodes:
                source_signal_nodes.add(src.node_pos)
            if src.node_neg not in ground_nodes:
                source_signal_nodes.add(src.node_neg)

        # 直列チェーンを構築（ソースの信号ノードから右へ）
        # まず全信号ノードの順序を決める
        ordered_nodes = self._order_signal_nodes(
            parser, sources, series_comps, ground_nodes)

        # 信号ノードを上段に水平配置
        x = 0
        for node in ordered_nodes:
            self.node_positions[node] = NodePosition(x=x, y=self.TOP_Y)
            x += self.H_SPACING

        # GNDノード：接続されている信号ノードのX座標の最小値
        for gnd in ground_nodes:
            # GNDに直接繋がっている信号ノードのX座標を集める
            connected_x = []
            for comp in parser.components:
                other_node = None
                if comp.node_neg == gnd and comp.node_pos in self.node_positions:
                    other_node = comp.node_pos
                elif comp.node_pos == gnd and comp.node_neg in self.node_positions:
                    other_node = comp.node_neg

                if other_node:
                    connected_x.append(self.node_positions[other_node].x)

            if connected_x:
                gnd_x = min(connected_x)
            else:
                gnd_x = 0

            self.node_positions[gnd] = NodePosition(x=gnd_x, y=self.GND_Y)

    def _order_signal_nodes(self, parser, sources, series_comps,
                             ground_nodes) -> List[str]:
        """信号ノードを左から右の順序で並べる

        ソースの正端子から始めて、直列接続を辿る。
        """
        ordered = []
        visited = set()

        # ソースの信号側ノードを開始点にする
        start_nodes = []
        for src in sources:
            if src.node_pos not in ground_nodes:
                start_nodes.append(src.node_pos)
            elif src.node_neg not in ground_nodes:
                start_nodes.append(src.node_neg)

        if not start_nodes:
            # ソースがない場合、最初のコンポーネントのノードから
            if parser.components:
                comp = parser.components[0]
                if comp.node_pos not in ground_nodes:
                    start_nodes.append(comp.node_pos)
                if comp.node_neg not in ground_nodes:
                    start_nodes.append(comp.node_neg)

        # BFSで信号ノードを辿る
        queue = list(start_nodes)
        for node in queue:
            if node in visited or node in ground_nodes:
                continue
            visited.add(node)
            ordered.append(node)

            # このノードに接続された他の信号ノードを探す
            for comp in parser.components:
                neighbors = []
                all_nodes = [comp.node_pos, comp.node_neg]
                if comp.node_ctrl:
                    all_nodes.append(comp.node_ctrl)
                if comp.node_ctrl2:
                    all_nodes.append(comp.node_ctrl2)
                for n in all_nodes:
                    if n == node:
                        continue
                    if n not in ground_nodes and n not in visited:
                        neighbors.append(n)
                queue.extend(neighbors)

        # 残りの信号ノード（到達できなかったもの）
        all_signal = parser.get_signal_nodes()
        for node in sorted(all_signal):
            if node not in visited and node not in ground_nodes:
                ordered.append(node)

        return ordered

    def _calc_parallel_offsets(self, parser, sources, series_comps,
                                shunt_comps, ground_nodes) -> Dict[str, int]:
        """並列コンポーネントの水平オフセットを計算"""
        offsets: Dict[str, int] = {}

        # 同じノードペアを共有するコンポーネントをグループ化
        node_pair_groups: Dict[Tuple[str, str], List[Component]] = {}
        for comp in parser.components:
            # ソートしたノードペア（順序無関係に同じペアを検出）
            pair = tuple(sorted([comp.node_pos, comp.node_neg]))
            node_pair_groups.setdefault(pair, []).append(comp)

        for pair, comps in node_pair_groups.items():
            if len(comps) <= 1:
                continue

            # 並列: ソースを除いた素子にオフセットを付ける
            non_source = [c for c in comps
                         if c.comp_type not in (ComponentType.VOLTAGE,
                                                 ComponentType.CURRENT)]
            source = [c for c in comps
                     if c.comp_type in (ComponentType.VOLTAGE,
                                         ComponentType.CURRENT)]

            # ソースはオフセットなし
            all_to_offset = source + non_source
            n = len(all_to_offset)
            for i, comp in enumerate(all_to_offset):
                offsets[comp.name] = i * self.H_SPACING

        # Dense input rails such as boost/buck fronts often put a source,
        # a series inductor/switch path, and one or more input capacitors on
        # the same signal node.  Without a local escape offset the series
        # symbol can be placed almost on top of the source symbol.  Keep the
        # electrical connection via FLAGs, but move the series branch right.
        source_signal_nodes: Set[str] = set()
        for src in sources:
            if src.node_pos not in ground_nodes:
                source_signal_nodes.add(src.node_pos)
            if src.node_neg not in ground_nodes:
                source_signal_nodes.add(src.node_neg)

        shunt_nodes: Set[str] = set()
        for comp in shunt_comps:
            if comp.node_pos not in ground_nodes:
                shunt_nodes.add(comp.node_pos)
            if comp.node_neg not in ground_nodes:
                shunt_nodes.add(comp.node_neg)

        dense_source_nodes = source_signal_nodes & shunt_nodes
        if dense_source_nodes:
            escape = max(96, self.H_SPACING // 2)
            for comp in series_comps:
                if (
                    comp.node_pos in dense_source_nodes
                    or comp.node_neg in dense_source_nodes
                ):
                    offsets[comp.name] = max(offsets.get(comp.name, 0), escape)

        return offsets

    def _place_component(self, comp: Component,
                          ground_nodes: Set[str],
                          h_offset: int = 0) -> PlacedComponent:
        """コンポーネントを端子ノード間に配置"""
        # C2 (isolation zone): multi-pin SUBCIRCUITs go to a dedicated
        # vertical band far from regular components so the round-trip
        # `_estimate_terminals` only picks up our own FLAG positions
        # (and not adjacent components' WIRE endpoints).
        if (comp.comp_type == ComponentType.SUBCIRCUIT
                and len(comp.extra_nodes) > 0):
            idx = self._subckt_iso_idx
            self._subckt_iso_idx += 1
            ISO_X_BASE = 4096   # far right of any normal layout
            ISO_X_STEP = 384    # per-SUBCIRCUIT horizontal spacing
            ISO_Y = 0           # above the regular schematic
            iso_x = ISO_X_BASE + idx * ISO_X_STEP
            iso_y = ISO_Y
            # terminal1/terminal2 are nominal — AscGenerator._write_symbol
            # emits its own FLAG+WIRE grid around (iso_x, iso_y) for every
            # pin (including extras), and `_generate_wires` is suppressed
            # for this PlacedComponent (see flag below).
            pc = PlacedComponent(
                component=comp,
                x=iso_x, y=iso_y, rotation='R0',
                terminal1=(iso_x, iso_y),
                terminal2=(iso_x + 32, iso_y),
            )
            return pc

        pos_node = self.node_positions.get(comp.node_pos, NodePosition(0, 0))
        neg_node = self.node_positions.get(comp.node_neg, NodePosition(0, 0))

        is_source = comp.comp_type in (ComponentType.VOLTAGE,
                                        ComponentType.CURRENT)
        pos_is_gnd = comp.node_pos in ground_nodes
        neg_is_gnd = comp.node_neg in ground_nodes

        # 並列オフセット適用
        ox = h_offset

        # G(VCCS)素子はpin1(out+)が下(R0時)。通常のvoltage/resはpin1が上。
        # t1=node_pos→pin1, t2=node_neg→pin2 なので、
        # G素子でpos=GND → pin1(下)=GND → 自然にR0配置で正しい。
        # G素子でneg=GND → pin2(上)=GND → 信号が下、GNDが上 → t1/t2入替で対応。
        is_g_type = comp.comp_type in (ComponentType.VCCS,)

        # 分路素子: GND側の座標を信号側ノードの真下に調整
        if neg_is_gnd and not is_source:
            if is_g_type:
                # G neg=GND: pin2(上)=GND → swap t1/t2 so pin1(下)=signal
                t1 = (pos_node.x + ox, self.GND_Y)
                t2 = (pos_node.x + ox, pos_node.y)
            else:
                t1 = (pos_node.x + ox, pos_node.y)
                t2 = (pos_node.x + ox, self.GND_Y)
        elif pos_is_gnd and not is_source:
            if is_g_type:
                # G pos=GND(out+): pin1(下)=GND → R0 natural
                t1 = (neg_node.x + ox, self.GND_Y)  # pin1(out+)=GND=下
                t2 = (neg_node.x + ox, neg_node.y)   # pin2(out-)=signal=上
            else:
                t1 = (neg_node.x + ox, self.GND_Y)
                t2 = (neg_node.x + ox, neg_node.y)
        elif is_source:
            # ソース: 信号ノードを上、GNDを下に固定
            if neg_is_gnd:
                signal_node = pos_node
                t1 = (signal_node.x + ox, self.TOP_Y)
                t2 = (signal_node.x + ox, self.GND_Y)
            elif pos_is_gnd:
                signal_node = neg_node
                t1 = (signal_node.x + ox, self.GND_Y)
                t2 = (signal_node.x, self.TOP_Y)
            else:
                t1 = (pos_node.x, pos_node.y)
                t2 = (neg_node.x, neg_node.y)
        else:
            t1 = (pos_node.x + ox, pos_node.y)
            t2 = (neg_node.x + ox, neg_node.y)

        # 方向と回転を決定
        dx = t2[0] - t1[0]
        dy = t2[1] - t1[1]

        if is_source:
            # ソースは常に垂直配置
            if comp.comp_type == ComponentType.VOLTAGE:
                rotation = 'R0'
            else:
                rotation = 'R180'
        elif is_g_type:
            # G素子: pin1(out+)=下(R0), t1=pin1位置
            # t1が下(GND_Y)ならR0, t1が上ならR180
            rotation = 'R0' if t1[1] >= t2[1] else 'R180'
        elif abs(dx) > abs(dy):
            # 水平配置
            rotation = 'R90' if dx > 0 else 'R270'
        elif abs(dy) > 0:
            # 垂直配置
            rotation = 'R0' if dy > 0 else 'R180'
        else:
            rotation = 'R0'

        # 4端子素子の処理 (E, G, SW, T with control terminals)
        is_4t = (comp.comp_type in (ComponentType.VCVS, ComponentType.VCCS,
                                     ComponentType.SWITCH, ComponentType.TLINE)
                 and comp.node_ctrl != '')
        if is_4t:
            t1_pos = (pos_node.x + ox, pos_node.y)
            if pos_is_gnd:
                t1_pos = (neg_node.x + ox, self.GND_Y)
            rotation = 'R0'

            sym_x, sym_y, p1, p2, p3, p4 = \
                calc_symbol_placement_4t(comp.comp_type, rotation, t1_pos)

            return PlacedComponent(
                component=comp,
                x=sym_x, y=sym_y, rotation=rotation,
                terminal1=p1, terminal2=p2,
                terminal3=p3, terminal4=p4,
            )

        # Opamp素子の処理
        if comp.comp_type == ComponentType.OPAMP:
            rotation = 'R0'
            # In+の目標座標（node_pos）
            t1_pos = (pos_node.x + ox, pos_node.y)

            if comp.node_out:
                # 5-pin opamp: In+ In- V+ V- OUT model
                sym_name = comp.value.lower() if comp.value else 'universalopamp2'
                result = calc_symbol_placement_5t(sym_name, rotation, t1_pos)
                sym_x, sym_y, p1, p2, p3, p4, p5 = result
                return PlacedComponent(
                    component=comp,
                    x=sym_x, y=sym_y, rotation=rotation,
                    terminal1=p1, terminal2=p2,
                    terminal3=p3, terminal4=p4,
                    terminal5=p5,
                )
            else:
                # 3-pin opamp: inv noninv out
                result = calc_symbol_placement_3t_opamp(rotation, t1_pos)
                sym_x, sym_y, p_inv, p_noninv, p_out = result
                return PlacedComponent(
                    component=comp,
                    x=sym_x, y=sym_y, rotation=rotation,
                    terminal1=p_inv, terminal2=p_noninv,
                    terminal3=p_out,
                )

        # 3端子素子の処理
        is_3t = comp.comp_type in (ComponentType.BJT, ComponentType.MOSFET,
                                    ComponentType.JFET)
        if is_3t:
            # 3端子: C/D を上（pos_node）、E/S を下（GND or neg_node）に配置
            # B/G は左に出る（R0配置の場合）
            ctrl_node = self.node_positions.get(comp.node_ctrl, NodePosition(0, 0))

            # C/D 端子の目標座標
            t1_cd = (pos_node.x + ox, pos_node.y)
            if pos_is_gnd:
                t1_cd = (neg_node.x + ox, self.GND_Y)

            # R0 固定（縦配置）
            rotation = 'R0'

            sym_x, sym_y, actual_cd, actual_bg, actual_es = \
                calc_symbol_placement_3t(comp.comp_type, rotation, t1_cd)

            return PlacedComponent(
                component=comp,
                x=sym_x,
                y=sym_y,
                rotation=rotation,
                terminal1=actual_cd,
                terminal2=actual_es,
                terminal3=actual_bg,
            )

        # 2端子: シンボル位置を計算
        sym_x, sym_y, actual_t1, actual_t2 = calc_symbol_placement(
            comp.comp_type, rotation, t1, t2)

        return PlacedComponent(
            component=comp,
            x=sym_x,
            y=sym_y,
            rotation=rotation,
            terminal1=actual_t1,
            terminal2=actual_t2,
        )


# =============================================================================
# 3. AscGenerator - .ascファイル生成
# =============================================================================

class AscGenerator:
    """LTSpice .ascファイルジェネレータ"""

    # コンポーネント種別 → LTSpiceシンボル名
    SYMBOL_MAP = {
        ComponentType.RESISTOR: 'res',
        ComponentType.CAPACITOR: 'cap',
        ComponentType.INDUCTOR: 'ind',
        ComponentType.VOLTAGE: 'voltage',
        ComponentType.CURRENT: 'current',
        ComponentType.DIODE: 'diode',
        ComponentType.BJT: 'npn',
        ComponentType.MOSFET: 'nmos',
        ComponentType.JFET: 'njf',
        ComponentType.VCVS: 'e',
        ComponentType.VCCS: 'g',
        ComponentType.CCCS: 'f',
        ComponentType.CCVS: 'h',
        ComponentType.BEHAVIORAL: 'bv',
        ComponentType.SWITCH: 'sw',
        ComponentType.TLINE: 'tline',
        ComponentType.OPAMP: 'Opamps\\\\opamp',
    }

    # 回転別WINDOW設定
    WINDOW_CONFIGS = {
        'res': {
            'R90':  [('WINDOW 0 0 56 VBottom 2',), ('WINDOW 3 32 56 VTop 2',)],
            'R270': [('WINDOW 0 32 56 VTop 2',), ('WINDOW 3 0 56 VBottom 2',)],
        },
        'cap': {
            'R90':  [('WINDOW 0 0 32 VBottom 2',), ('WINDOW 3 32 32 VTop 2',)],
            'R270': [('WINDOW 0 0 32 VBottom 2',), ('WINDOW 3 32 32 VTop 2',)],
        },
        'ind': {
            'R270': [('WINDOW 0 32 56 VTop 2',), ('WINDOW 3 5 56 VBottom 2',)],
        },
        'current': {
            'R180': [('WINDOW 0 24 80 Left 2',), ('WINDOW 3 24 0 Left 2',)],
        },
    }

    def __init__(self):
        self.lines: List[str] = []
        # Optional .asy search dirs (set by NetlistToAsc.__init__).
        self.asy_search_dirs: List = []

    def generate(self, layouter: CircuitLayouter, parser: NetlistParser,
                 sheet_width: int = 0, sheet_height: int = 0) -> str:
        """ASCファイル内容を生成"""
        self.lines = []

        placed = layouter.placed_components
        node_positions = layouter.node_positions
        local_model_names = self._local_model_names(parser)
        self._raw_directive_component_ids = {
            id(pc.component)
            for pc in placed
            if self._should_emit_as_raw_directive(pc, local_model_names)
        }

        # シートサイズの自動計算
        if sheet_width == 0 or sheet_height == 0:
            sw, sh = self._calc_sheet_size(placed, node_positions)
            if sheet_width == 0:
                sheet_width = sw
            if sheet_height == 0:
                sheet_height = sh

        # ヘッダ
        self.lines.append('Version 4')
        self.lines.append(f'SHEET 1 {sheet_width} {sheet_height}')

        # ワイヤ生成
        wires = self._generate_wires(placed, node_positions, parser)
        for w in wires:
            self.lines.append(f'WIRE {w[0]} {w[1]} {w[2]} {w[3]}')

        # 全端子にスタブワイヤ + FLAGを配置
        # LTspiceはFLAGがワイヤ端点上にないとピンに接続しない
        # (C2 exclusion mirrors _generate_wires above.)
        node_terminal_map: Dict[str, List[Tuple[int, int]]] = {}
        for pc in placed:
            comp = pc.component
            if self._is_raw_directive_component(comp):
                continue
            if (comp.comp_type == ComponentType.SUBCIRCUIT
                    and len(comp.extra_nodes) > 0):
                continue
            node_terminal_map.setdefault(comp.node_pos, []).append(pc.terminal1)
            node_terminal_map.setdefault(comp.node_neg, []).append(pc.terminal2)
            if pc.terminal3 is not None and comp.node_ctrl:
                node_terminal_map.setdefault(comp.node_ctrl, []).append(pc.terminal3)
            if pc.terminal4 is not None and comp.node_ctrl2:
                node_terminal_map.setdefault(comp.node_ctrl2, []).append(pc.terminal4)
            if pc.terminal5 is not None and comp.node_out:
                node_terminal_map.setdefault(comp.node_out, []).append(pc.terminal5)

        ground_nodes = parser.get_ground_nodes()

        # スタブワイヤ端点にFLAGを配置
        flag_positions = getattr(self, '_flag_positions', {})
        placed_flags: Set[Tuple[int, int, str]] = set()

        for node_name in node_terminal_map:
            flag_name = '0' if node_name in ground_nodes else node_name
            pts = flag_positions.get(node_name, [])
            for pt in pts:
                key = (pt[0], pt[1], flag_name)
                if key not in placed_flags:
                    placed_flags.add(key)
                    self.lines.append(f'FLAG {pt[0]} {pt[1]} {flag_name}')

        # シンボル（コンポーネント）
        for pc in placed:
            self._write_symbol(pc)

        # ディレクティブ（.xxx, K文, A文）
        directive_y = sheet_height - 100
        for i, directive in enumerate(parser.directives):
            text = directive.text
            if text.startswith('.') or (text and text[0].upper() in ('K', 'A', 'J')):
                self.lines.append(
                    f'TEXT 0 {directive_y + i * 32} Left 2 !{text}'
                )

        # A文ディレクティブのノード名にFLAGを追加
        # (A文はコンポーネントとして配置されないため、ノードが未登録の場合がある)
        existing_flags = {f for _, _, f in placed_flags}
        flag_x = sheet_width - 200
        flag_y_start = directive_y - 100
        a_flag_count = 0
        for directive in parser.directives:
            text = directive.text
            if text and text[0].upper() == 'A':
                parts = text.split()
                if len(parts) >= 10:
                    # A-element: name pin1..pin8 model [params]
                    for node in parts[1:9]:
                        if (node != '0' and node.lower() != 'gnd'
                                and node not in existing_flags
                                and not node.startswith('NC_')):
                            y = flag_y_start - a_flag_count * 32
                            self.lines.append(f'FLAG {flag_x} {y} {node}')
                            existing_flags.add(node)
                            a_flag_count += 1

        return '\n'.join(self.lines)

    @staticmethod
    def _split_value_ac(comp: Component) -> tuple:
        """Split AC stimulus from value for voltage/current sources.

        LTspice .asc uses SYMATTR Value for DC/waveform and SYMATTR Value2
        for AC stimulus. In SPICE netlist they appear together:
            V1 n1 n2 SINE(0 1 1k) AC 1
        This splits them back:
            Value  = SINE(0 1 1k)   or  ""  (if DC=0)
            Value2 = AC 1
        """
        if comp.comp_type not in (ComponentType.VOLTAGE, ComponentType.CURRENT):
            return comp.value, ''

        val = comp.value
        if not val:
            return val, ''

        # Find " AC " (case-insensitive) not inside parentheses
        val_upper = val.upper()
        # Check for AC at various positions
        import re
        # Match AC followed by value, but not inside SINE/PULSE/etc parentheses
        m = re.search(r'\s+AC\s+', val_upper)
        if m:
            # Check we're not inside parentheses
            before = val[:m.start()]
            paren_depth = before.count('(') - before.count(')')
            if paren_depth == 0:
                dc_part = val[:m.start()].strip()
                ac_part = val[m.start():].strip()
                return dc_part, ac_part

        # Also handle standalone "AC" at the beginning: "AC 1"
        if val_upper.startswith('AC ') or val_upper == 'AC':
            return '', val

        return val, ''

    @staticmethod
    def _subckt_pins(comp: Component) -> List[str]:
        """Return positional pins for an X-subcircuit invocation."""
        return [
            p for p in [comp.node_pos, comp.node_neg] + list(comp.extra_nodes)
            if p
        ]

    def _subckt_terminal_offsets(self, comp: Component, rotation: str):
        """Look up .asy terminal offsets for a subcircuit symbol, if present."""
        if comp.comp_type != ComponentType.SUBCIRCUIT:
            return None
        sym_name = comp.value or ''
        pins = self._subckt_pins(comp)
        if not sym_name or not pins:
            return None
        from .asc_parser import AsyParser as _AsyParser  # noqa: E402
        return _AsyParser.get_terminal_offsets(
            sym_name,
            rotation,
            search_dirs=self.asy_search_dirs or None,
        )

    def _should_emit_subckt_as_raw_directive(self, pc: PlacedComponent) -> bool:
        """Prefer lossless raw X-lines when a matching .asy is unavailable.

        A generic unknown SYMBOL plus estimated pins can corrupt the X-line by
        inventing extra nodes. When we do not have the symbol definition, the
        simulation-correct representation is an LTspice SPICE directive.
        """
        comp = pc.component
        if comp.comp_type != ComponentType.SUBCIRCUIT or not comp.raw_line:
            return False
        pins = self._subckt_pins(comp)
        offsets = self._subckt_terminal_offsets(comp, pc.rotation)
        return not offsets or len(offsets) < len(pins)

    @staticmethod
    def _local_model_names(parser: NetlistParser) -> Set[str]:
        names: Set[str] = set()
        for directive in parser.directives:
            match = re.match(r"\s*\.model\s+(\S+)", directive.text, re.IGNORECASE)
            if match:
                names.add(match.group(1).lower())
        return names

    @staticmethod
    def _component_model_name(comp: Component) -> str:
        if not comp.value:
            return ''
        return comp.value.split()[0].lower()

    def _should_emit_as_raw_directive(
        self,
        pc: PlacedComponent,
        local_model_names: Set[str],
    ) -> bool:
        comp = pc.component
        if self._should_emit_subckt_as_raw_directive(pc):
            return True
        if not comp.raw_line:
            return False
        if comp.comp_type not in {
            ComponentType.BJT,
            ComponentType.MOSFET,
            ComponentType.JFET,
            ComponentType.DIODE,
        }:
            return False
        model_name = self._component_model_name(comp)
        return bool(model_name and model_name in local_model_names)

    def _is_raw_directive_component(self, comp: Component) -> bool:
        return id(comp) in getattr(self, '_raw_directive_component_ids', set())

    def _is_raw_directive_subckt(self, comp: Component) -> bool:
        return self._is_raw_directive_component(comp)

    # Opampモデル名 → LTSpiceシンボルパス
    OPAMP_SYMBOL_MAP = {
        'opamp':            'Opamps\\\\opamp',
        'universalopamp':   'Opamps\\\\UniversalOpamp',
        'universalopamp1':  'Opamps\\\\UniversalOpamp1',
        'universalopamp2':  'Opamps\\\\UniversalOpamp2',
        'universalopamp3':  'Opamps\\\\UniversalOpamp3',
        'universalopamp3a': 'Opamps\\\\UniversalOpamp3a',
        'universalopamp3b': 'Opamps\\\\UniversalOpamp3b',
        'universalopamp4':  'Opamps\\\\UniversalOpamp4',
        'universalopamp5':  'Opamps\\\\UniversalOpamp5',
        'lt1001':           'Opamps\\\\LT1001',
        'lt1001a':          'Opamps\\\\LT1001A',
        'lt1028':           'Opamps\\\\LT1028',
        'lt1028a':          'Opamps\\\\LT1028A',
    }

    def _write_symbol(self, pc: PlacedComponent):
        """コンポーネントシンボルを書き出す"""
        comp = pc.component
        if self._is_raw_directive_component(comp):
            self.lines.append(f'TEXT {pc.x} {pc.y} Left 2 !{comp.raw_line}')
            return
        sym_name = self.SYMBOL_MAP.get(comp.comp_type, 'res')

        # SUBCIRCUIT (X-class, including unknown vendor symbols): emit the
        # original symbol name so a round-trip preserves the SYMBOL line.
        # The netlist's X§<name> prefix is the asc_parser convention for
        # "this was a vendor symbol whose InstName did not start with X";
        # strip it so the InstName reads correctly in the new .asc.
        if comp.comp_type == ComponentType.SUBCIRCUIT:
            inst_name = comp.name
            if inst_name.startswith('X§'):
                inst_name = inst_name[2:]
            sym_name = comp.value or 'res'
            self.lines.append(f'SYMBOL {sym_name} {pc.x} {pc.y} {pc.rotation}')
            self.lines.append(f'SYMATTR InstName {inst_name}')
            if comp.value:
                # Preserve the model name for re-extraction.
                self.lines.append(f'SYMATTR SpiceModel {comp.value}')

            # Emit one FLAG per pin at the EXACT offset reported by
            # AsyParser.get_terminal_offsets(<sym>, <rotation>). asc_parser
            # uses the same lookup on re-extraction, so the pin-to-node
            # mapping survives intact — provided the .asy file is on the
            # LTspice library search path. Without a matching .asy, generate()
            # normally emits the original X-line as a raw SPICE directive.
            #
            # Imported here to avoid a top-level circular import; this
            # path is only hit for SUBCIRCUITs.
            all_pins = self._subckt_pins(comp)
            offsets = self._subckt_terminal_offsets(comp, pc.rotation)

            if offsets and len(offsets) >= len(all_pins):
                # .asy lookup succeeded — place each FLAG at the canonical
                # pin location. asc_parser will pick the same coordinates
                # on re-extraction, so node names round-trip exactly.
                for i, node in enumerate(all_pins):
                    ox, oy = offsets[i]
                    fx = pc.x + ox
                    fy = pc.y + oy
                    flag_name = '0' if node == '0' or node.lower() == 'gnd' else node
                    # 1-unit stub wire so the FLAG endpoint is a wire
                    # endpoint (LTspice convention: FLAGs only attach to
                    # wire endpoints).
                    self.lines.append(f'WIRE {fx} {fy} {fx+8} {fy}')
                    self.lines.append(f'FLAG {fx} {fy} {flag_name}')
            elif all_pins:
                # Defensive fallback for direct _write_symbol callers that did
                # not run generate()'s raw-directive classifier first.
                DX = 32
                DY = 16
                seen_flags: Set[Tuple[int, int, str]] = set()
                seen_wires: Set[Tuple[int, int, int, int]] = set()
                for i, node in enumerate(all_pins):
                    fx = pc.x + DX
                    fy = pc.y + DY * (i + 1)
                    flag_name = '0' if node == '0' or node.lower() == 'gnd' else node
                    # Long horizontal wire to give the FLAG a wire endpoint
                    # to attach to (LTspice convention).  Far end shifted by
                    # i so multiple stubs at the same y don't merge.
                    far_x = pc.x + DX + 200 + i * 4
                    far_y = fy
                    w = (fx, fy, far_x, far_y)
                    if w not in seen_wires:
                        seen_wires.add(w)
                        self.lines.append(f'WIRE {fx} {fy} {far_x} {far_y}')
                    k = (fx, fy, flag_name)
                    if k not in seen_flags:
                        seen_flags.add(k)
                        self.lines.append(f'FLAG {fx} {fy} {flag_name}')
            return

        # `* @sym=<kind>` hint takes priority for any class — restores
        # variants like ind2, schottky, polcap, pnp, npn3, etc.
        if comp.symbol_hint:
            sym_name = comp.symbol_hint

        # PNP/PMOS/PJF 判定: モデル名に "pnp"/"pmos"/"pjf" が含まれるか、
        # または名前のプレフィックスから推定
        # NOTE: when symbol_hint is set, trust it and skip the model-name
        # heuristic (which is lossy by definition).
        if comp.comp_type == ComponentType.BJT and not comp.symbol_hint:
            val_lower = (comp.value or '').lower()
            if 'pnp' in val_lower or comp.name.upper().startswith('QP'):
                sym_name = 'pnp'
        elif comp.comp_type == ComponentType.MOSFET and not comp.symbol_hint:
            val_lower = (comp.value or '').lower()
            if 'pmos' in val_lower or 'pch' in val_lower:
                sym_name = 'pmos'
        elif comp.comp_type == ComponentType.JFET and not comp.symbol_hint:
            val_lower = (comp.value or '').lower()
            if 'pjf' in val_lower:
                sym_name = 'pjf'
        elif comp.comp_type == ComponentType.OPAMP and not comp.symbol_hint:
            val_lower = (comp.value or '').lower()
            sym_name = self.OPAMP_SYMBOL_MAP.get(
                val_lower, 'Opamps\\\\UniversalOpamp2')
        elif (comp.comp_type in (ComponentType.VCVS, ComponentType.VCCS)
                and not comp.symbol_hint):
            # E/G without control terminals → use bv/bi (2-pin behavioral)
            # to avoid floating nc_ pins on the 4-pin g/e symbol
            if not comp.node_ctrl:
                if comp.comp_type == ComponentType.VCVS:
                    sym_name = 'bv'
                else:
                    sym_name = 'bi'

        self.lines.append(f'SYMBOL {sym_name} {pc.x} {pc.y} {pc.rotation}')

        # WINDOW設定
        win_config = self.WINDOW_CONFIGS.get(sym_name, {}).get(pc.rotation)
        if win_config:
            for win_lines in win_config:
                for wl in win_lines:
                    self.lines.append(wl)

        self.lines.append(f'SYMATTR InstName {comp.name}')
        if comp.comp_type == ComponentType.OPAMP:
            # Write SpiceModel for opamp/subcircuit so re-extraction preserves model name
            if comp.value:
                val_lower = comp.value.lower()
                if val_lower not in self.OPAMP_SYMBOL_MAP:
                    self.lines.append(f'SYMATTR SpiceModel {comp.value}')
        elif comp.value:
            value, value2 = self._split_value_ac(comp)
            if value:
                self.lines.append(f'SYMATTR Value {value}')
            if value2:
                if not value:
                    # AC-only source: promote Value2 to Value
                    self.lines.append(f'SYMATTR Value {value2}')
                else:
                    self.lines.append(f'SYMATTR Value2 {value2}')

    def _generate_wires(self, placed: List[PlacedComponent],
                         node_positions: Dict[str, NodePosition],
                         parser: NetlistParser) -> List[Tuple[int, int, int, int]]:
        """ワイヤ（配線）を生成 — FLAG + ピンスタブ方式

        FLAGはワイヤ端点上にないとLTspiceがピンに接続しない。
        各ピンから短いスタブワイヤを引き、その端点にFLAGを配置する。
        スタブは全て同じ方向（下向き16px）で、交差を最小限に。
        """
        wires = []
        wire_set = set()
        self._flag_positions = {}  # node_name -> list of (x, y) for FLAG

        # 全ピン位置を収集
        # C2: multi-pin SUBCIRCUITs are isolated (see _place_component) and
        # emit their own FLAG+WIRE grid in _write_symbol. Exclude them here
        # so the auto-router does not connect their nominal terminals back
        # to regular nodes (which would create spurious wire endpoints
        # near isolated SYMBOLs and defeat the isolation).
        node_pins: Dict[str, List[Tuple[int, int]]] = {}
        for pc in placed:
            comp = pc.component
            if self._is_raw_directive_subckt(comp):
                continue
            if (comp.comp_type == ComponentType.SUBCIRCUIT
                    and len(comp.extra_nodes) > 0):
                continue
            node_pins.setdefault(comp.node_pos, []).append(pc.terminal1)
            node_pins.setdefault(comp.node_neg, []).append(pc.terminal2)
            if pc.terminal3 is not None and comp.node_ctrl:
                node_pins.setdefault(comp.node_ctrl, []).append(pc.terminal3)
            if pc.terminal4 is not None and comp.node_ctrl2:
                node_pins.setdefault(comp.node_ctrl2, []).append(pc.terminal4)
            if pc.terminal5 is not None and comp.node_out:
                node_pins.setdefault(comp.node_out, []).append(pc.terminal5)

        # 全ピン座標と全ワイヤセグメントを追跡
        all_pin_coords = set()
        for pts in node_pins.values():
            all_pin_coords.update(pts)

        # スタブ端点とワイヤセグメントの追跡（T字接続回避）
        occupied_points = set(all_pin_coords)
        placed_segments: List[Tuple[int, int, int, int]] = []

        def point_on_segment(px, py, x1, y1, x2, y2):
            """点(px,py)がワイヤ(x1,y1)-(x2,y2)上にあるか（端点除く）"""
            if y1 == y2 == py:  # 水平
                return min(x1, x2) < px < max(x1, x2)
            if x1 == x2 == px:  # 垂直
                return min(y1, y2) < py < max(y1, y2)
            return False

        def is_safe_stub(pt, stub_end):
            """スタブ端点が既存ワイヤ上にないか確認"""
            if stub_end in occupied_points:
                return False
            # スタブ端点が既存ワイヤセグメントの中間にないか
            for seg in placed_segments:
                if point_on_segment(stub_end[0], stub_end[1], *seg):
                    return False
            # 新スタブワイヤが既存のスタブ端点を通過しないか
            for seg_end in occupied_points - all_pin_coords:
                if point_on_segment(seg_end[0], seg_end[1],
                                   pt[0], pt[1], stub_end[0], stub_end[1]):
                    return False
            return True

        for node_name, pins in node_pins.items():
            unique_pts = list(dict.fromkeys(pins))
            flag_pts = []

            for pt in unique_pts:
                stub_end = None
                for dx, dy in [(0, 16), (0, -16), (-16, 0), (16, 0),
                               (0, 32), (0, -32), (-32, 0), (32, 0)]:
                    candidate = (pt[0] + dx, pt[1] + dy)
                    if is_safe_stub(pt, candidate):
                        stub_end = candidate
                        break
                if stub_end is None:
                    stub_end = (pt[0], pt[1] + 16)  # fallback

                w = (pt[0], pt[1], stub_end[0], stub_end[1])
                if w not in wire_set:
                    wire_set.add(w)
                    wires.append(w)
                    placed_segments.append(w)
                occupied_points.add(stub_end)
                flag_pts.append(stub_end)

            self._flag_positions[node_name] = flag_pts

        return wires

    def _make_orthogonal_wires(self, x1: int, y1: int,
                                x2: int, y2: int
                                ) -> List[Tuple[int, int, int, int]]:
        """2点間の直交ワイヤを生成（L字ルーティング）"""
        wires = []

        if x1 == x2 and y1 == y2:
            return wires  # 長さ0

        if x1 == x2 or y1 == y2:
            # 既に直線
            wires.append((x1, y1, x2, y2))
        else:
            # L字ルーティング：まず水平、次に垂直
            wires.append((x1, y1, x2, y1))
            wires.append((x2, y1, x2, y2))

        return wires

    def _find_label_nodes(self, parser: NetlistParser) -> List[str]:
        """ラベルを付けるべきノードを特定

        - グランド以外で、名前が数字のみでないノード
        - 複数のコンポーネントが接続するノード
        """
        ground = parser.get_ground_nodes()
        node_count: Dict[str, int] = {}
        for comp in parser.components:
            node_count[comp.node_pos] = node_count.get(comp.node_pos, 0) + 1
            node_count[comp.node_neg] = node_count.get(comp.node_neg, 0) + 1

        labels = []
        for node_name in sorted(parser.get_signal_nodes()):
            if node_name in ground:
                continue
            # 名前が意味のある文字列（数字だけでない）の場合ラベル付け
            if not node_name.isdigit():
                labels.append(node_name)
            # または接続数が3以上（ジャンクション）の場合
            elif node_count.get(node_name, 0) >= 3:
                labels.append(node_name)

        return labels

    def _calc_sheet_size(self, placed: List[PlacedComponent],
                          node_positions: Dict[str, NodePosition]
                          ) -> Tuple[int, int]:
        """シートサイズを自動計算"""
        max_x = 400
        max_y = 400
        min_x = 0
        min_y = 0

        for pc in placed:
            max_x = max(max_x, pc.x + 200, pc.terminal1[0] + 100,
                       pc.terminal2[0] + 100)
            max_y = max(max_y, pc.y + 200, pc.terminal1[1] + 100,
                       pc.terminal2[1] + 100)
            min_x = min(min_x, pc.x, pc.terminal1[0], pc.terminal2[0])
            min_y = min(min_y, pc.y, pc.terminal1[1], pc.terminal2[1])

        for pos in node_positions.values():
            max_x = max(max_x, pos.x + 200)
            max_y = max(max_y, pos.y + 200)

        # ディレクティブ用スペース
        max_y += 150

        # 最小サイズ
        max_x = max(max_x, 880)
        max_y = max(max_y, 680)

        return (snap(max_x), snap(max_y))

    def _snap(self, val: float) -> int:
        """グリッドスナップ"""
        return snap(val)


# =============================================================================
# 4. NetlistToAsc - 統合クラス
# =============================================================================

class NetlistToAsc:
    """ネットリスト → ASC変換の統合インターフェース"""

    def __init__(self, asy_search_dirs: Optional[List[str]] = None):
        """Initialize converter.

        Args:
            asy_search_dirs: additional directories to scan for `.asy`
                files (used by multi-pin SUBCIRCUIT pin-offset lookup).
                Falls back to the ``LTSPICE_ASY_SEARCH_PATH`` env var
                (os.pathsep-separated) when None.
        """
        self.parser = NetlistParser()
        self.layouter = CircuitLayouter()
        self.generator = AscGenerator()
        # Resolve search dirs: explicit arg > env var > none.
        if asy_search_dirs is None:
            import os as _os
            env = _os.environ.get('LTSPICE_ASY_SEARCH_PATH', '')
            asy_search_dirs = [d for d in env.split(_os.pathsep) if d.strip()]
        self.generator.asy_search_dirs = [
            __import__('pathlib').Path(d) for d in (asy_search_dirs or [])
        ]

    def convert_file(self, input_path: str, output_path: str = None) -> str:
        """.cirファイルを.ascファイルに変換"""
        self.parser.parse_file(input_path)
        self.layouter.layout(self.parser)
        asc_content = self.generator.generate(self.layouter, self.parser)

        if output_path is None:
            output_path = input_path.rsplit('.', 1)[0] + '.asc'

        # LTspice reads .asc as latin-1/cp1252, not UTF-8. Writing UTF-8
        # corrupts characters like µ (0xB5 in latin-1) into the byte
        # sequence 0xC2 0xB5 which LTspice interprets as `Âµ` — making
        # `.tran 250µ` round-trip as `.tran 250Âµ` (interpreted as 250s
        # instead of 250µs, a 10⁶× factor in the .raw time axis).
        with open(output_path, 'w', encoding='latin-1', errors='replace') as f:
            f.write(asc_content)

        print(f'ASC saved: {output_path}')
        print(f'  Components: {len(self.parser.components)}')
        print(f'  Nodes: {len(self.parser.get_all_nodes())}')
        print(f'  Directives: {len(self.parser.directives)}')

        return asc_content

    def convert_string(self, netlist: str, output_path: str = None) -> str:
        """ネットリスト文字列を.ascに変換"""
        self.parser.parse_string(netlist)
        self.layouter.layout(self.parser)
        asc_content = self.generator.generate(self.layouter, self.parser)

        if output_path:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(asc_content)
            print(f'ASC saved: {output_path}')

        return asc_content


# =============================================================================
# テスト用
# =============================================================================

if __name__ == '__main__':
    # テスト1: 簡単なRC回路
    netlist_rc = """\
* RC Lowpass Filter
V1 in 0 AC 1
R1 in out 1k
C1 out 0 1u
.ac dec 100 1 100k
.end
"""

    print("=== Test 1: RC Lowpass ===")
    converter = NetlistToAsc()
    asc = converter.convert_string(netlist_rc)
    print(asc)
    print()

    # テスト2: RLC直列回路
    netlist_rlc = """\
* RLC Series Circuit
I1 0 in AC 1
R1 in mid 100
L1 mid out 10m
C1 out 0 1u
.ac dec 100 10 100k
.end
"""

    print("=== Test 2: RLC Series ===")
    asc2 = converter.convert_string(netlist_rlc)
    print(asc2)
    print()

    # テスト3: 並列LC回路
    netlist_plc = """\
* Parallel LC Tank
I1 0 top AC 1
L1 top 0 1m
C1 top 0 100n
.ac dec 100 1k 100k
.end
"""

    print("=== Test 3: Parallel LC ===")
    asc3 = converter.convert_string(netlist_plc)
    print(asc3)
