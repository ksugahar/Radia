#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
schemdraw Drawing to SPICE netlist (.cir) converter

schemdraw Drawingオブジェクトからトポロジーを抽出し、SPICEネットリストを生成する。

アルゴリズム:
1. drawing.elements を走査し、各要素のアンカー座標を取得
2. 同一座標（tolerance内）のアンカーを同一ノードとしてグループ化
3. Ground() の位置をノード "0" に割り当て
4. 各コンポーネント要素をSPICE行にマッピング
"""

import re
import sys
import math
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict


# =============================================================================
# schemdraw要素名 → SPICEプレフィクスのマッピング
# =============================================================================

SCHEMDRAW_TO_SPICE = {
    # Resistors (elm.Resistor() creates ResistorIEEE internally)
    'Resistor': 'R', 'ResistorIEEE': 'R', 'ResistorIEC': 'R',
    'Resistor2': 'R', 'ResistorVar': 'R',
    # Capacitors
    'Capacitor': 'C', 'Capacitor2': 'C', 'CapacitorVar': 'C',
    'CapacitorTrim': 'C',
    # Inductors
    'Inductor': 'L', 'Inductor2': 'L',
    # Sources
    'SourceV': 'V', 'BatteryCell': 'V', 'Battery': 'V',
    'SourceI': 'I',
    # Diodes
    'Diode': 'D', 'Schottky': 'D', 'Zener': 'D', 'LED': 'D',
    # BJTs
    'BjtNpn': 'Q', 'BjtPnp': 'Q',
    # MOSFETs
    'NFet': 'M', 'PFet': 'M',
    # JFETs
    'JFetN': 'J', 'JFetP': 'J',
    # Opamp
    'Opamp': 'X',
    # Controlled sources
    'SourceControlledV': 'E', 'SourceControlledI': 'G',
    # Switch
    'Switch': 'S',
    # Transmission line (cir_to_schemdraw emits ``elm.Coax()`` for T)
    'Coax': 'T',
}

# 3端子要素のアンカー名マッピング
THREE_TERMINAL_ANCHORS = {
    'BjtNpn': ('collector', 'base', 'emitter'),
    'BjtPnp': ('collector', 'base', 'emitter'),
    'NFet': ('drain', 'gate', 'source'),
    'PFet': ('drain', 'gate', 'source'),
    'JFetN': ('drain', 'gate', 'source'),
    'JFetP': ('drain', 'gate', 'source'),
}

# Opampのアンカー名
OPAMP_ANCHORS = ('in1', 'in2', 'out')

# 非コンポーネント要素（ノード接続やアノテーションのみ）
NON_COMPONENT_TYPES = {
    'Ground', 'GroundSignal', 'GroundChassis',
    'Dot', 'DotDotDot', 'Arrowhead', 'Gap',
    'Label', 'CurrentLabel', 'CurrentLabelInline',
    'Annotate', 'LoopCurrent', 'LoopArrow',
    'Vss', 'Vdd',
}

# Line系要素（ノード接続のみ、コンポーネントではない）
LINE_TYPES = {'Line', 'Arrow', 'LineDot', 'Wire'}


# =============================================================================
# ユーティリティ
# =============================================================================

def _point_key(pt, eps=0.05):
    """座標をグリッドに丸めてハッシュ可能なキーにする"""
    scale = 1.0 / eps
    x = round(float(pt[0]) * scale) / scale
    y = round(float(pt[1]) * scale) / scale
    return (x, y)


def _get_anchor_pos(element, anchor_name):
    """要素のアンカー位置を取得（absanchorsを使用）"""
    anchors = getattr(element, 'absanchors', None)
    if anchors is None:
        return None
    pos = anchors.get(anchor_name)
    if pos is None:
        return None
    try:
        return (float(pos[0]), float(pos[1]))
    except (TypeError, IndexError):
        return None


def _parse_label(element) -> Tuple[Optional[str], Optional[str]]:
    """要素のラベルからコンポーネント名と値を抽出する。

    schemdraw stores user labels in _userlabels list, each with a .label attribute.

    Returns:
        (name, value) - どちらもNoneの可能性あり
    """
    label_text = None

    # Primary method: _userlabels (schemdraw's actual storage)
    userlabels = getattr(element, '_userlabels', [])
    if userlabels:
        first = userlabels[0]
        lbl = getattr(first, 'label', None)
        if lbl:
            label_text = str(lbl)

    # Fallback: _userparams dict
    if not label_text:
        userparams = getattr(element, '_userparams', {})
        if isinstance(userparams, dict):
            lbl = userparams.get('label', None)
            if lbl:
                label_text = str(lbl)

    if not label_text:
        return None, None

    # "\n" で分割: "R1\n1k" → name="R1", value="1k"
    parts = label_text.split('\n')
    if len(parts) >= 2:
        return parts[0].strip(), '\n'.join(p.strip() for p in parts[1:])
    elif len(parts) == 1:
        text = parts[0].strip()
        tokens = text.split()
        if len(tokens) >= 2:
            etype = element.__class__.__name__
            expected_prefix = SCHEMDRAW_TO_SPICE.get(etype, '')
            first = tokens[0]
            if expected_prefix and first.upper().startswith(expected_prefix):
                return first, ' '.join(tokens[1:])
        # プレフィクス+数字の場合は名前のみ
        if re.match(r'^[A-Za-z]+\d+$', text):
            return text, None
        return None, text

    return None, None


# =============================================================================
# Union-Find for node merging
# =============================================================================

class UnionFind:
    """Union-Find data structure for merging coordinate-based nodes."""

    def __init__(self):
        self.parent: Dict[Tuple[float, float], Tuple[float, float]] = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]  # path compression
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


# =============================================================================
# メインコンバータ
# =============================================================================

class SchemdrawToCir:
    """schemdraw Drawing -> SPICE netlist converter"""

    def __init__(self, eps: float = 0.05):
        self.eps = eps
        self.uf = UnionFind()
        # グラウンド座標キーのセット
        self.ground_keys: Set[Tuple[float, float]] = set()
        # ルートキー → ノード名
        self.root_to_node: Dict[Tuple[float, float], str] = {}
        # ノードカウンタ
        self._next_node: int = 1
        # コンポーネントプレフィクスごとのカウンタ
        self._comp_counters: Dict[str, int] = defaultdict(int)
        # SPICE行リスト
        self.spice_lines: List[str] = []
        # Original component lines embedded by cir_to_schemdraw.  When
        # present, these are more faithful than geometry-derived netlists
        # for multi-pin / behavioral components.
        self.raw_spice_lines: List[str] = []
        # SPICEディレクティブ
        self.directives: List[str] = []
        # タイトル
        self.title: str = ''

    def _reset(self):
        self.uf = UnionFind()
        self.ground_keys = set()
        self.root_to_node = {}
        self._next_node = 1
        self._comp_counters = defaultdict(int)
        self.spice_lines = []
        self.raw_spice_lines = []
        self.directives = []
        self.title = ''

    def _key(self, pt):
        return _point_key(pt, self.eps)

    def _get_node(self, coord: Tuple[float, float]) -> str:
        """座標に対応するノード名を返す"""
        key = self._key(coord)
        root = self.uf.find(key)

        # グラウンドチェック
        if root in self.ground_keys:
            self.root_to_node[root] = '0'
            return '0'

        if root in self.root_to_node:
            return self.root_to_node[root]

        node_name = str(self._next_node)
        self._next_node += 1
        self.root_to_node[root] = node_name
        return node_name

    def _auto_name(self, prefix: str) -> str:
        """自動コンポーネント名を生成"""
        self._comp_counters[prefix] += 1
        return f"{prefix}{self._comp_counters[prefix]}"

    def _pass1_collect_all_connections(self, drawing):
        """Pass 1: 全要素のアンカー座標をUnion-Findに登録し、
        接続関係（Line, Dot, コンポーネント端子の位置一致）を構築。
        Ground位置も収集。"""

        all_anchor_keys = []

        for elem in drawing.elements:
            type_name = type(elem).__name__

            # Ground要素: startアンカーの位置をグラウンドとしてマーク
            if type_name in ('Ground', 'GroundSignal', 'GroundChassis'):
                pos = _get_anchor_pos(elem, 'start')
                if pos is not None:
                    key = self._key(pos)
                    self.uf.find(key)  # register
                    self.ground_keys.add(key)
                continue

            # Dot要素: 位置を登録（他の要素との接続点）
            if type_name == 'Dot':
                pos = _get_anchor_pos(elem, 'center') or _get_anchor_pos(elem, 'start')
                if pos is not None:
                    self.uf.find(self._key(pos))
                continue

            # Line系要素: start-endを同一ノードに統合
            if type_name in LINE_TYPES:
                start = _get_anchor_pos(elem, 'start')
                end = _get_anchor_pos(elem, 'end')
                if start is not None and end is not None:
                    self.uf.union(self._key(start), self._key(end))
                continue

            # コンポーネント要素: 各アンカーを登録
            if type_name in THREE_TERMINAL_ANCHORS:
                for aname in THREE_TERMINAL_ANCHORS[type_name]:
                    pos = _get_anchor_pos(elem, aname)
                    if pos is not None:
                        self.uf.find(self._key(pos))
            elif type_name == 'Opamp':
                for aname in OPAMP_ANCHORS:
                    pos = _get_anchor_pos(elem, aname)
                    if pos is not None:
                        self.uf.find(self._key(pos))
            else:
                for aname in ('start', 'end'):
                    pos = _get_anchor_pos(elem, aname)
                    if pos is not None:
                        self.uf.find(self._key(pos))

        # Ground キーをルートに伝搬
        new_ground_keys = set()
        for gk in self.ground_keys:
            root = self.uf.find(gk)
            new_ground_keys.add(root)
        self.ground_keys = new_ground_keys

    def _pass2_merge_coincident(self, drawing):
        """Pass 2: 同じ座標（tolerance内）にあるアンカーを統合。
        これはLine経由でなくても位置が一致する端子同士を接続する。

        全要素の全アンカー座標を集めて、同一キーのものをunionする。"""
        # 座標キー → 最初のキー（代表）
        key_to_representative: Dict[Tuple[float, float], Tuple[float, float]] = {}

        for elem in drawing.elements:
            type_name = type(elem).__name__
            if type_name in NON_COMPONENT_TYPES:
                # Ground等は既にpass1で処理済み
                # ただし位置は他の要素と統合する必要がある
                if type_name in ('Ground', 'GroundSignal', 'GroundChassis'):
                    pos = _get_anchor_pos(elem, 'start')
                    if pos is not None:
                        key = self._key(pos)
                        if key in key_to_representative:
                            self.uf.union(key, key_to_representative[key])
                        else:
                            key_to_representative[key] = key
                elif type_name == 'Dot':
                    pos = _get_anchor_pos(elem, 'center') or _get_anchor_pos(elem, 'start')
                    if pos is not None:
                        key = self._key(pos)
                        if key in key_to_representative:
                            self.uf.union(key, key_to_representative[key])
                        else:
                            key_to_representative[key] = key
                continue

            # 全アンカー位置を集める
            anchors_to_check = []
            if type_name in LINE_TYPES:
                anchors_to_check = ['start', 'end']
            elif type_name in THREE_TERMINAL_ANCHORS:
                anchors_to_check = list(THREE_TERMINAL_ANCHORS[type_name])
            elif type_name == 'Opamp':
                anchors_to_check = list(OPAMP_ANCHORS)
            else:
                anchors_to_check = ['start', 'end']

            for aname in anchors_to_check:
                pos = _get_anchor_pos(elem, aname)
                if pos is not None:
                    key = self._key(pos)
                    if key in key_to_representative:
                        self.uf.union(key, key_to_representative[key])
                    else:
                        key_to_representative[key] = key

        # Ground キーをルートに再伝搬
        new_ground_keys = set()
        for gk in self.ground_keys:
            root = self.uf.find(gk)
            new_ground_keys.add(root)
        # Also check all keys that share a root with any ground key
        for key in list(self.uf.parent.keys()):
            root = self.uf.find(key)
            if root in new_ground_keys:
                continue
            # Check if this key was originally a ground
            if key in self.ground_keys:
                new_ground_keys.add(root)
        self.ground_keys = new_ground_keys

    def _process_element(self, elem):
        """1つのschemdraw要素をSPICE行に変換"""
        type_name = type(elem).__name__

        # 非コンポーネント要素はスキップ
        if type_name in NON_COMPONENT_TYPES or type_name in LINE_TYPES:
            return

        # SPICEプレフィクスを取得
        prefix = SCHEMDRAW_TO_SPICE.get(type_name)
        if prefix is None:
            return

        # ラベルから名前と値を抽出
        label_name, label_value = _parse_label(elem)

        # コンポーネント名を決定
        comp_name = None
        if label_name:
            # ラベル名がSPICEプレフィクスで始まるか確認
            if label_name[0].upper() == prefix[0].upper():
                comp_name = label_name
            else:
                comp_name = prefix + label_name
        if not comp_name:
            comp_name = self._auto_name(prefix)

        # 値を決定
        value = label_value if label_value else ''

        # ノード接続を抽出
        if type_name in THREE_TERMINAL_ANCHORS:
            anchor_names = THREE_TERMINAL_ANCHORS[type_name]
            nodes = []
            for aname in anchor_names:
                pos = _get_anchor_pos(elem, aname)
                if pos is not None:
                    nodes.append(self._get_node(pos))
                else:
                    nodes.append('0')
            model = value if value else type_name
            if prefix == 'Q':
                self.spice_lines.append(
                    f"{comp_name} {nodes[0]} {nodes[1]} {nodes[2]} {model}")
            elif prefix == 'M':
                self.spice_lines.append(
                    f"{comp_name} {nodes[0]} {nodes[1]} {nodes[2]} {nodes[2]} {model}")
            elif prefix == 'J':
                self.spice_lines.append(
                    f"{comp_name} {nodes[0]} {nodes[1]} {nodes[2]} {model}")
        elif type_name == 'Opamp':
            nodes = []
            for aname in OPAMP_ANCHORS:
                pos = _get_anchor_pos(elem, aname)
                if pos is not None:
                    nodes.append(self._get_node(pos))
                else:
                    nodes.append('0')
            model = value if value else 'OPAMP'
            self.spice_lines.append(
                f"{comp_name} {nodes[0]} {nodes[1]} {nodes[2]} {model}")
        else:
            # 2端子要素
            start_pos = _get_anchor_pos(elem, 'start')
            end_pos = _get_anchor_pos(elem, 'end')
            if start_pos is None or end_pos is None:
                return

            node_start = self._get_node(start_pos)
            node_end = self._get_node(end_pos)

            if prefix in ('V', 'I'):
                # schemdraw convention: current flows from end(+) to start(-)
                node_pos = node_end
                node_neg = node_start
                val_str = value if value else 'DC 1'
                self.spice_lines.append(
                    f"{comp_name} {node_pos} {node_neg} {val_str}")
            elif prefix in ('E', 'G'):
                val_str = value if value else '1'
                self.spice_lines.append(
                    f"{comp_name} {node_start} {node_end} {val_str}")
            elif prefix == 'D':
                model = value if value else 'D'
                self.spice_lines.append(
                    f"{comp_name} {node_start} {node_end} {model}")
            elif prefix == 'T':
                # SPICE transmission line: T<name> A+ A- B+ B- [params]
                # schemdraw's Coax has 2 signal anchors (start, end); the
                # shield/return is implicit (and grounded in practice).
                # Map A+ = start, A- = 0, B+ = 0, B- = end; preserves
                # component count and signal-path semantics on round-trip.
                val_str = value if value else ''
                self.spice_lines.append(
                    f"{comp_name} {node_start} 0 0 {node_end} {val_str}".rstrip())
            else:
                # R, L, C
                val_str = value if value else '1'
                self.spice_lines.append(
                    f"{comp_name} {node_start} {node_end} {val_str}")

    def _collect_directives_and_node_names(self, drawing):
        """全要素のラベルからSPICEディレクティブとノード名マッピングを収集。

        cir_to_schemdraw が埋め込んだ Annotate ラベルを検出:
        - '.tran ...', '.ac ...' 等 → ディレクティブ
        - 'NODE:name' → ノード名マッピング（numeric→named に復元）
        """
        node_names_found = []  # (order_index, name) のリスト

        for elem in drawing.elements:
            label_name, label_value = _parse_label(elem)
            for text in (label_name, label_value):
                if not text:
                    continue
                text = text.strip()
                # NetlistParser classifies K (mutual inductance) and A
                # (digital primitive) lines as directives alongside the
                # `.tran` / `.ac` / `.model` / etc. families. cir_to_schemdraw
                # emits all four as Annotate labels, so we must accept all
                # four prefixes here too -- previously only `.`-prefix was
                # picked up and the K/A statements silently dropped on
                # round-trip.
                if not text:
                    continue
                first = text[0]
                if text.startswith('SPICELINE:'):
                    self.raw_spice_lines.append(text[len('SPICELINE:'):])
                elif first == '.' or first.upper() in ('K', 'A'):
                    self.directives.append(text)
                elif text.startswith('NODE:'):
                    node_names_found.append(text[5:])

        # ノード名の復元: numeric ノード → named ノード
        # spice_lines 中のノード番号を、発見順に名前付きノードに置換
        if node_names_found:
            # 現在使われている numeric ノード（0以外）を収集
            used_nodes = set()
            for line in self.spice_lines:
                parts = line.split()
                if len(parts) >= 3 and parts[0][0].isalpha():
                    for p in parts[1:]:
                        if p != '0' and re.match(r'^\d+$', p):
                            used_nodes.add(p)

            # numeric ノードを昇順で並べ、名前を割り当て
            sorted_nodes = sorted(used_nodes, key=lambda x: int(x))
            rename_map = {}
            for i, num_node in enumerate(sorted_nodes):
                if i < len(node_names_found):
                    rename_map[num_node] = node_names_found[i]

            # spice_lines のノード位置のみ置換（値部分は置換しない）
            if rename_map:
                new_lines = []
                for line in self.spice_lines:
                    parts = line.split()
                    if len(parts) >= 3 and parts[0][0].isalpha():
                        prefix = parts[0][0].upper()
                        # ノード位置を決定（コンポーネントタイプ依存）
                        if prefix == 'Q':      # Q name C B E model
                            node_indices = {1, 2, 3}
                        elif prefix == 'M':    # M name D G S B model
                            node_indices = {1, 2, 3, 4}
                        elif prefix == 'J':    # J name D G S model
                            node_indices = {1, 2, 3}
                        elif prefix == 'X':    # X name n1 n2 ... subckt
                            node_indices = set(range(1, len(parts) - 1))
                        else:                  # 2端子: name n+ n- value...
                            node_indices = {1, 2}
                        new_parts = []
                        for idx, p in enumerate(parts):
                            if idx in node_indices:
                                new_parts.append(rename_map.get(p, p))
                            else:
                                new_parts.append(p)
                        new_lines.append(' '.join(new_parts))
                    else:
                        new_lines.append(line)
                self.spice_lines = new_lines

    def convert(self, drawing, title: str = 'schemdraw circuit') -> str:
        """schemdraw Drawingオブジェクト -> SPICEネットリスト文字列"""
        self._reset()
        self.title = title

        # Pass 1: Ground収集 + Line統合
        self._pass1_collect_all_connections(drawing)

        # Pass 2: 同一座標のアンカーを統合
        self._pass2_merge_coincident(drawing)

        # Finalize ground roots: any key whose root shares a group with a ground key
        final_ground_roots = set()
        for gk in self.ground_keys:
            final_ground_roots.add(self.uf.find(gk))
        self.ground_keys = final_ground_roots

        # Pass 3: コンポーネントをSPICE行に変換
        for elem in drawing.elements:
            self._process_element(elem)

        # Pass 4: ディレクティブとノード名を収集・復元
        self._collect_directives_and_node_names(drawing)

        return self._format_netlist()

    def _format_netlist(self) -> str:
        """SPICEネットリスト文字列を生成"""
        lines = []
        lines.append(f"* {self.title}")
        lines.extend(self.raw_spice_lines if self.raw_spice_lines else self.spice_lines)
        for d in self.directives:
            # NetlistParser packs multi-line directives (e.g. .subckt ... .ends
            # bodies) into a single SpiceDirective.text with embedded real
            # newlines. cir_to_schemdraw stores them as Annotate labels with
            # the newlines escaped to literal "\n" so the Python string fits
            # on one source line. Restore real newlines on the way back out
            # so each internal component reappears on its own SPICE line and
            # any downstream tooling (lint, count-based round-trip checks)
            # sees the .subckt body as the multi-line block it really is.
            d_restored = d.replace('\\n', '\n')
            lines.append(d_restored)
        if not any(d.lower().startswith('.end') for d in self.directives):
            lines.append('.end')
        return '\n'.join(lines)


# =============================================================================
# 公開API関数
# =============================================================================

def schemdraw_to_cir(drawing, title: str = 'schemdraw circuit') -> str:
    """Extract SPICE netlist from a schemdraw Drawing object.

    Args:
        drawing: schemdraw.Drawing instance
        title: netlist title line

    Returns:
        SPICE netlist string
    """
    converter = SchemdrawToCir()
    return converter.convert(drawing, title)


def schemdraw_script_to_cir(
    script: str,
    title: str = 'schemdraw circuit',
    source_path: str | Path | None = None,
) -> str:
    """Execute a schemdraw script and extract the netlist.

    The script should create a Drawing and assign it to a variable 'd' or 'drawing',
    or use 'with schemdraw.Drawing() as d:' pattern.

    Args:
        script: Python source code string that creates a schemdraw Drawing
        title: netlist title line
        source_path: optional original script path.  When provided, it is
            exposed as ``__file__`` while the script runs so user scripts can
            resolve adjacent assets or helper modules.

    Returns:
        SPICE netlist string
    """
    import schemdraw
    import schemdraw.elements as elm
    script_path = Path(source_path) if source_path is not None else Path.cwd() / "<schemdraw-script>"

    namespace = {
        'schemdraw': schemdraw,
        'elm': elm,
        '__file__': str(script_path),
        '__name__': '__ltspice_converter_schemdraw__',
        '__builtins__': __builtins__,
    }

    # Capture Drawing objects by monkey-patching __init__
    original_init = schemdraw.Drawing.__init__
    captured_drawings = []

    def patched_init(self, *args, **kwargs):
        kwargs['show'] = False
        original_init(self, *args, **kwargs)
        captured_drawings.append(self)

    schemdraw.Drawing.__init__ = patched_init
    try:
        exec(script, namespace)
    finally:
        schemdraw.Drawing.__init__ = original_init

    # Find the Drawing object
    drawing = None
    if captured_drawings:
        drawing = captured_drawings[-1]
    else:
        for name in ('d', 'drawing', 'dwg', 'fig'):
            if name in namespace and isinstance(namespace[name], schemdraw.Drawing):
                drawing = namespace[name]
                break
        if drawing is None:
            for v in namespace.values():
                if isinstance(v, schemdraw.Drawing):
                    drawing = v
                    break

    if drawing is None:
        raise ValueError("No schemdraw Drawing found in script")

    return schemdraw_to_cir(drawing, title)


def schemdraw_file_to_cir(py_path: str, output_path: str = None) -> str:
    """Convert a schemdraw .py file to .cir.

    Args:
        py_path: path to the .py file containing schemdraw code
        output_path: path to write the .cir file (default: same name, .cir extension)

    Returns:
        SPICE netlist string
    """
    py_path = Path(py_path)
    script = py_path.read_text(encoding='utf-8')
    title = py_path.stem

    netlist = schemdraw_script_to_cir(script, title, source_path=py_path)

    if output_path is None:
        output_path = str(py_path.with_suffix('.cir'))

    Path(output_path).write_text(netlist, encoding='utf-8')
    print(f"Saved: {output_path}")
    return netlist


# =============================================================================
# テスト
# =============================================================================

if __name__ == '__main__':
    import schemdraw
    import schemdraw.elements as elm

    def run_test(name, drawing, expected_components):
        """Run a test and verify output."""
        netlist = schemdraw_to_cir(drawing, name)
        print(f"{'=' * 60}")
        print(f"Test: {name}")
        print(f"{'=' * 60}")
        print(netlist)

        # Basic validation
        lines = [l for l in netlist.split('\n') if l and not l.startswith('*') and not l.startswith('.')]
        ok = True
        for comp in expected_components:
            found = any(comp in line for line in lines)
            if not found:
                print(f"  [FAIL] Expected component '{comp}' not found!")
                ok = False
        if ok:
            print(f"  [PASS] All {len(expected_components)} expected components found")
        print()
        return ok

    results = []

    # --- Test 1: RC Lowpass Filter ---
    # V1(+)--R1--+--out
    #            |
    #            C1
    #            |
    # V1(-)-----+--GND
    with schemdraw.Drawing(show=False) as d:
        V1 = d.add(elm.SourceV().up().label('V1\nAC 1', loc='left'))
        d.add(elm.Line().right())
        R1 = d.add(elm.Resistor().right().label('R1\n1k'))
        d.add(elm.Dot())
        d.push()
        C1 = d.add(elm.Capacitor().down().label('C1\n1u', loc='left'))
        d.add(elm.Line().left().tox(V1.start))
        d.add(elm.Ground())

    results.append(run_test('RC Lowpass Filter', d, ['V1', 'R1', 'C1', '1k', '1u', 'AC 1']))

    # --- Test 2: Voltage Divider ---
    # V1(+)--R1--+--out
    #            |
    #            R2
    #            |
    # V1(-)-----+--GND
    with schemdraw.Drawing(show=False) as d:
        V1 = d.add(elm.SourceV().up().label('V1\nDC 10', loc='left'))
        d.add(elm.Line().right())
        R1 = d.add(elm.Resistor().right().label('R1\n1k'))
        d.add(elm.Dot())
        d.push()
        R2 = d.add(elm.Resistor().down().label('R2\n1k', loc='left'))
        d.add(elm.Line().left().tox(V1.start))
        d.add(elm.Ground())

    results.append(run_test('Voltage Divider', d, ['V1', 'R1', 'R2', 'DC 10']))

    # --- Test 3: RLC Circuit ---
    # V1(+)--R1--L1--+--out
    #                 |
    #                 C1
    #                 |
    # V1(-)----------+--GND
    with schemdraw.Drawing(show=False) as d:
        V1 = d.add(elm.SourceV().up().label('V1\nAC 1', loc='left'))
        d.add(elm.Line().right())
        R1 = d.add(elm.Resistor().right().label('R1\n100'))
        L1 = d.add(elm.Inductor2().right().label('L1\n10m'))
        d.add(elm.Dot())
        d.push()
        C1 = d.add(elm.Capacitor().down().label('C1\n1u', loc='left'))
        d.add(elm.Line().left().tox(V1.start))
        d.add(elm.Ground())

    results.append(run_test('RLC Circuit', d, ['V1', 'R1', 'L1', 'C1', '100', '10m', '1u']))

    # --- Test 4: Script-based conversion ---
    print("=" * 60)
    print("Test 4: Script-based conversion")
    print("=" * 60)
    test_script = """\
import schemdraw
import schemdraw.elements as elm

with schemdraw.Drawing(show=False) as d:
    V1 = d.add(elm.SourceV().up().label('V1\\nDC 5', loc='left'))
    R1 = d.add(elm.Resistor().right().label('R1\\n470'))
    d.add(elm.Line().down().toy(V1.start))
    d.add(elm.Line().left().tox(V1.start))
    d.add(elm.Ground())
"""
    netlist = schemdraw_script_to_cir(test_script, 'Script Test')
    print(netlist)
    ok = 'V1' in netlist and 'R1' in netlist and '470' in netlist and 'DC 5' in netlist
    print(f"  [{'PASS' if ok else 'FAIL'}] Script-based conversion")
    results.append(ok)
    print()

    # --- Summary ---
    print("=" * 60)
    print(f"Results: {sum(results)}/{len(results)} tests passed")
    print("=" * 60)
