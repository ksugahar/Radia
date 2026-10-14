#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPICE Netlist (.cir) to schemdraw Python script converter

SPICEネットリストからschemdraw回路図スクリプトを直接生成する。
LLMが.cirテキストを生成 → このコンバータで回路図化 → LTspiceで検証

アルゴリズム:
1. NetlistParserでネットリスト解析
2. ノードを信号フローで整列（BFS）
3. コンポーネントをseries/shunt/sourceに分類
4. 並列素子検出 (同一ノードペア)
5. 配置ベースでschemdrawコードを生成
"""

import re
import sys
from pathlib import Path
from typing import List, Dict, Set, Tuple, Optional
from dataclasses import dataclass, field
from collections import defaultdict, deque

from .netlist_to_asc import NetlistParser, Component, ComponentType, SpiceDirective


# =============================================================================
# SPICE → schemdraw 要素名マッピング
# =============================================================================

SPICE_TO_SCHEMDRAW = {
    ComponentType.RESISTOR: 'elm.Resistor()',
    ComponentType.CAPACITOR: 'elm.Capacitor()',
    ComponentType.INDUCTOR: 'elm.Inductor2()',
    ComponentType.VOLTAGE: 'elm.SourceV()',
    ComponentType.CURRENT: 'elm.SourceI()',
    ComponentType.DIODE: 'elm.Diode()',
    ComponentType.BJT: None,      # NPN/PNP で分岐
    ComponentType.MOSFET: None,   # NMOS/PMOS で分岐
    ComponentType.JFET: None,     # NJF/PJF で分岐
    ComponentType.VCVS: 'elm.SourceControlledV()',
    ComponentType.VCCS: 'elm.SourceControlledI()',
    ComponentType.CCVS: 'elm.SourceControlledV()',
    ComponentType.CCCS: 'elm.SourceControlledI()',
    ComponentType.BEHAVIORAL: 'elm.SourceV()',
    ComponentType.SWITCH: 'elm.Switch()',
    ComponentType.TLINE: 'elm.Coax()',
}


def _get_schemdraw_element(comp: Component) -> str:
    """コンポーネントからschemdraw要素コードを返す。

    Returns only the element constructor, never includes comments.
    """
    ct = comp.comp_type

    if ct == ComponentType.BJT:
        if 'pnp' in comp.value.lower():
            return 'elm.BjtPnp()'
        return 'elm.BjtNpn()'
    elif ct == ComponentType.MOSFET:
        if 'pmos' in comp.value.lower() or comp.name[0:2].upper() == 'MP':
            return 'elm.PFet()'
        return 'elm.NFet()'
    elif ct == ComponentType.JFET:
        if 'pjf' in comp.value.lower():
            return 'elm.JFetP()'
        return 'elm.JFetN()'
    elif ct == ComponentType.SUBCIRCUIT:
        # Check both the parsed value and the raw line for opamp identifiers
        search_text = (comp.value + ' ' + comp.raw_line).lower()
        opamp_hints = ('opamp', 'opa', 'lm741', 'lt10', 'ad8', 'ad7',
                       'tl07', 'tl08', 'lm358', 'ne5532',
                       'lt1001', 'lt1006', 'lt1007', 'lt1012', 'lt1013',
                       'lt1028', 'universalopamp')
        if any(x in search_text for x in opamp_hints):
            return 'elm.Opamp()'
        # Return a box element for unknown subcircuits (no comment!)
        return 'elm.RBox()'

    elem = SPICE_TO_SCHEMDRAW.get(ct)
    return elem or 'elm.RBox()'


def _parse_3terminal(comp: Component) -> Tuple[str, str, str]:
    """3端子素子のノード名を返す: (collector/drain, base/gate, emitter/source)"""
    parts = comp.raw_line.split()
    ct = comp.comp_type
    if ct == ComponentType.BJT and len(parts) >= 5:
        return parts[1], parts[2], parts[3]  # C, B, E
    elif ct == ComponentType.MOSFET and len(parts) >= 6:
        return parts[1], parts[2], parts[3]  # D, G, S
    elif ct == ComponentType.JFET and len(parts) >= 5:
        return parts[1], parts[2], parts[3]  # D, G, S
    return comp.node_pos, '', comp.node_neg


def _is_3terminal(comp: Component) -> bool:
    return comp.comp_type in (ComponentType.BJT, ComponentType.MOSFET, ComponentType.JFET)


def _is_source(comp: Component) -> bool:
    return comp.comp_type in (ComponentType.VOLTAGE, ComponentType.CURRENT)


def _make_var(name: str) -> str:
    """SPICE名 → Python変数名 (sanitize special chars like §, +, -)"""
    # Replace § separator (LTSpice .net format) with underscore
    v = name.replace('\xa7', '_').replace('§', '_')
    # Replace any other non-alphanumeric chars
    v = re.sub(r'[^a-zA-Z0-9_]', '_', v)
    if not v or v[0].isdigit():
        v = '_' + v
    if v in ('in', 'is', 'as', 'or', 'if', 'for', 'not', 'and', 'del'):
        v += '_'
    return v


def _sanitize_label(text: str) -> str:
    """Sanitize text for use inside a Python single-quoted string literal.

    - Escapes single quotes and backslashes
    - Replaces non-ASCII (e.g. µ) with ASCII equivalents
    - Truncates long values
    """
    # Replace common non-ASCII
    text = text.replace('\xb5', 'u')  # µ → u
    text = text.replace('\u00b5', 'u')
    text = text.replace('\u03bc', 'u')  # Greek mu
    text = text.replace('\u2126', 'Ohm')  # Ω
    # Strip any remaining non-ASCII
    text = text.encode('ascii', 'replace').decode('ascii')
    # Escape backslashes and single quotes
    text = text.replace('\\', '\\\\')
    text = text.replace("'", "\\'")
    return text


def _make_label(comp: Component) -> str:
    """コンポーネントのラベル (sanitized for Python string)"""
    name_clean = _sanitize_label(comp.name)
    parts = [name_clean]
    if comp.value:
        val = comp.value
        if len(val) > 20:
            val = val[:17] + '...'
        parts.append(_sanitize_label(val))
    return ' '.join(parts)


def _node_pair_key(comp: Component) -> Tuple[str, str]:
    """Return a canonical (sorted) node pair for detecting parallel components."""
    a, b = comp.node_pos, comp.node_neg
    return (min(a, b), max(a, b))


# =============================================================================
# メインコンバータ
# =============================================================================

class CirToSchemdraw:
    """SPICE netlist → schemdraw script converter"""

    def __init__(self):
        self.parser = NetlistParser()
        self.lines: List[str] = []
        self.gnd_nodes: Set[str] = set()
        self.signal_nodes: Set[str] = set()
        self.node_order: List[str] = []  # 信号フロー順
        self.node_to_comps: Dict[str, List[Component]] = defaultdict(list)

    def convert_string(self, netlist: str, name: str = 'circuit') -> str:
        """ネットリスト文字列 → schemdrawスクリプト"""
        self.parser.parse_string(netlist)
        return self._generate(name)

    def convert_file(self, cir_path: str, output_path: str = None) -> str:
        """ネットリストファイル → schemdrawスクリプト"""
        self._parse_file_robust(cir_path)
        name = Path(cir_path).stem
        script = self._generate(name)

        if output_path is not None:
            Path(output_path).write_text(script, encoding='utf-8')
        return script

    def _parse_file_robust(self, filepath: str):
        """Parse file with encoding fallback (utf-8 → latin-1 → cp1252)."""
        for enc in ('utf-8', 'latin-1', 'cp1252'):
            try:
                with open(filepath, 'r', encoding=enc) as f:
                    lines = f.readlines()
                self.parser.parse_lines(lines)
                return
            except (UnicodeDecodeError, UnicodeError):
                continue
        # Last resort: read with errors='replace'
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
        self.parser.parse_lines(lines)

    def _generate(self, name: str) -> str:
        self.lines = []
        self.gnd_nodes = self.parser.get_ground_nodes()
        self.signal_nodes = self.parser.get_signal_nodes()

        # ノード→コンポーネントの隣接マップ
        self.node_to_comps = defaultdict(list)
        for comp in self.parser.components:
            self.node_to_comps[comp.node_pos].append(comp)
            self.node_to_comps[comp.node_neg].append(comp)
            if _is_3terminal(comp):
                _, ctrl, _ = _parse_3terminal(comp)
                if ctrl:
                    self.node_to_comps[ctrl].append(comp)

        # ノード順序を決定（信号フロー）
        self.node_order = self._order_nodes()

        # コンポーネントを分類
        sources, series, shunts, three_terms, others = self._classify_components()

        # 並列素子を検出（同一ノードペアの素子グループ）
        parallel_groups = self._find_parallel_groups(series + shunts)

        # コード生成
        self._emit_header(name)
        self._emit_circuit(sources, series, shunts, three_terms, others, parallel_groups)
        self._emit_footer(name)

        return '\n'.join(self.lines)

    def _order_nodes(self) -> List[str]:
        """BFSでノードを信号フロー順に並べる"""
        source_comps = [c for c in self.parser.components if _is_source(c)]

        start_nodes = []
        for sc in source_comps:
            for n in [sc.node_pos, sc.node_neg]:
                if n not in self.gnd_nodes:
                    start_nodes.append(n)

        if not start_nodes and self.signal_nodes:
            start_nodes = [sorted(self.signal_nodes)[0]]

        visited = set()
        order = []
        queue = deque(start_nodes)

        while queue:
            node = queue.popleft()
            if node in visited or node in self.gnd_nodes:
                continue
            visited.add(node)
            order.append(node)

            for comp in self.node_to_comps[node]:
                for n in [comp.node_pos, comp.node_neg]:
                    if n not in visited and n not in self.gnd_nodes:
                        queue.append(n)

        # 未到達ノードを追加
        for n in sorted(self.signal_nodes - visited):
            order.append(n)

        return order

    def _classify_components(self):
        """コンポーネントをsource/series/shunt/3terminal/otherに分類"""
        sources = []
        series = []
        shunts = []
        three_terms = []
        others = []

        for comp in self.parser.components:
            if _is_source(comp):
                sources.append(comp)
            elif _is_3terminal(comp):
                three_terms.append(comp)
            elif comp.node_pos in self.gnd_nodes or comp.node_neg in self.gnd_nodes:
                shunts.append(comp)
            elif comp.node_pos in self.signal_nodes and comp.node_neg in self.signal_nodes:
                series.append(comp)
            else:
                others.append(comp)

        return sources, series, shunts, three_terms, others

    def _find_parallel_groups(self, components: List[Component]) -> Dict[Tuple[str, str], List[Component]]:
        """Find groups of parallel components (same node pair)."""
        pair_map: Dict[Tuple[str, str], List[Component]] = defaultdict(list)
        for comp in components:
            key = _node_pair_key(comp)
            pair_map[key].append(comp)
        # Only return groups with 2+ components
        return {k: v for k, v in pair_map.items() if len(v) >= 2}

    def _emit_header(self, name: str):
        # Sanitize the name for use in Python string literals
        safe_name = _sanitize_label(name)
        self.lines.extend([
            '#!/usr/bin/env python3',
            '# -*- coding: utf-8 -*-',
            f'"""Auto-generated schemdraw script from {safe_name}.cir"""',
            '',
            'import schemdraw',
            'import schemdraw.elements as elm',
            '',
            f'with schemdraw.Drawing(show=False) as d:',
            f"    d.config(unit=3, font='Times New Roman')",
            '',
        ])

    def _emit_footer(self, name: str):
        safe_name = _sanitize_label(name)

        # Preserve source component lines for lossless reverse conversion
        # when this script was generated from a SPICE netlist.  The visible
        # drawing remains best-effort, but schemdraw->cir can recover exact
        # multi-pin / behavioral / vendor component lines from this metadata.
        for comp in self.parser.components:
            safe_line = _sanitize_label(comp.raw_line)
            self.lines.append(f"    d.add(elm.Annotate().at((0, -4)).label('SPICELINE:{safe_line}').color('white'))")

        # ディレクティブを不可視Labelとして埋め込む（schemdraw→.cir で復元可能）
        # 同時にコメントとしても残す（可読性のため）
        for directive in self.parser.directives:
            safe_dir = _sanitize_label(directive.text)
            self.lines.append(f"    d.add(elm.Annotate().at((0, -2)).label('{safe_dir}').color('white'))")
            self.lines.append(f"    # SPICE: {safe_dir}")

        # ノード名マッピングを埋め込む（schemdraw→.cir でラベル復元用）
        for node in self.node_order:
            nl = node.lower()
            if node.isdigit() or re.match(r'^[Nn]\d+$', node):
                continue
            safe_node = _sanitize_label(node)
            self.lines.append(f"    d.add(elm.Annotate().at((0, -3)).label('NODE:{safe_node}').color('white'))")

        # If the circuit ended up with no drawn elements (.asc had only
        # directives / TEXT blocks, no SYMBOLs), schemdraw cannot compute
        # axis bounds and `d.save()` raises `ValueError: Axis limits cannot
        # be NaN or Inf`. Add a single invisible Line to give matplotlib a
        # finite bounding box.
        # Save: prefer .pdf, fall back to .svg if the backend (e.g. SVG-only
        # on headless Linux without matplotlib) does not support PDF.
        self.lines.extend([
            '',
            "    # Guard against empty drawings (no SYMBOLs in source .asc)",
            "    if not any(s.segments for s in d.elements):",
            "        d.add(elm.Line().right().length(d.unit).color('white'))",
            f"    _saved = False",
            f"    for _ext in ('.pdf', '.svg'):",
            f"        try:",
            f"            d.save(f'{safe_name}' + _ext)",
            f"            print(f'Saved: {safe_name}' + _ext)",
            f"            _saved = True",
            f"            break",
            f"        except (ValueError, NotImplementedError):",
            f"            continue",
            f"    if not _saved:",
            f"        print('Drawing skipped (no usable schemdraw backend)')",
        ])

    def _emit_circuit(self, sources, series, shunts, three_terms, others, parallel_groups):
        """配置ベースで回路コードを生成"""
        placed = set()
        first_source_var = None

        # Track which components are in parallel groups (will use push/pop)
        parallel_comp_names = set()
        for group in parallel_groups.values():
            for comp in group:
                parallel_comp_names.add(comp.name)

        # --- ソースを配置 (左端、上向き) ---
        for i, comp in enumerate(sources):
            var = _make_var(comp.name)
            elem = _get_schemdraw_element(comp)
            label = _make_label(comp)

            if i == 0:
                self.lines.append(f"    {var} = d.add({elem}.up().label('{label}', loc='top'))")
                first_source_var = var
            else:
                self.lines.append(f"    {var} = d.add({elem}.up().label('{label}', loc='top'))")
            placed.add(comp.name)

        if not first_source_var and series:
            first_source_var = '_start'
            self.lines.append(f"    _start = d.add(elm.Dot())")

        # --- ノード順に直列素子を配置 (右向き) ---
        node_idx = {n: i for i, n in enumerate(self.node_order)}

        def series_sort_key(c):
            i1 = node_idx.get(c.node_pos, 999)
            i2 = node_idx.get(c.node_neg, 999)
            return min(i1, i2)

        series_sorted = sorted(series, key=series_sort_key)

        last_var = first_source_var
        nodes_with_dot = set()

        for comp in series_sorted:
            if comp.name in placed:
                continue

            var = _make_var(comp.name)
            elem = _get_schemdraw_element(comp)
            label = _make_label(comp)

            # Determine input node
            node_in = comp.node_pos if node_idx.get(comp.node_pos, 999) < node_idx.get(comp.node_neg, 999) else comp.node_neg

            # Check if we need a branch point
            has_shunt = any(s for s in shunts if s.name not in placed and
                          (s.node_pos == node_in or s.node_neg == node_in))
            has_branch = any(s for s in series_sorted if s.name not in placed and s.name != comp.name and
                           (s.node_pos == node_in or s.node_neg == node_in))

            if (has_shunt or has_branch) and node_in not in nodes_with_dot:
                self.lines.append(f"    d.add(elm.Dot())")
                nodes_with_dot.add(node_in)

                # このノードのshunt素子を先にpush/popで配置
                self._emit_shunts_at_node(node_in, shunts, placed, parallel_groups)

            # Check if this series component is parallel with another
            pair_key = _node_pair_key(comp)
            if pair_key in parallel_groups and comp.name not in placed:
                group = parallel_groups[pair_key]
                first_in_group = True
                for pc in group:
                    if pc.name in placed:
                        continue
                    pvar = _make_var(pc.name)
                    pelem = _get_schemdraw_element(pc)
                    plabel = _make_label(pc)
                    if first_in_group:
                        self.lines.append(f"    d.push()")
                        self.lines.append(f"    {pvar} = d.add({pelem}.right().label('{plabel}'))")
                        first_in_group = False
                    else:
                        self.lines.append(f"    d.pop()")
                        self.lines.append(f"    d.push()")
                        self.lines.append(f"    d.add(elm.Line().down().length(d.unit/2))")
                        self.lines.append(f"    {pvar} = d.add({pelem}.right().label('{plabel}'))")
                        self.lines.append(f"    d.add(elm.Line().up().length(d.unit/2))")
                    placed.add(pc.name)
                    last_var = pvar
                self.lines.append(f"    d.pop()")
                # Move past the parallel group. `{last_var}.end` is only
                # defined on 2-terminal elements (Resistor, Capacitor, ...).
                # For multi-pin elements (Opamp, RBox, BJT, ...), there is
                # no `.end` anchor, so we fall back to `d.move()` which is
                # safe regardless of the element type.
                self.lines.append(
                    f"    _e = getattr({last_var}, 'end', None)\n"
                    f"    if _e is not None:\n"
                    f"        d.add(elm.Line().right().at(_e).length(0))\n"
                    f"    else:\n"
                    f"        d.move(d.unit, 0)"
                )
            else:
                self.lines.append(f"    {var} = d.add({elem}.right().label('{label}'))")
                placed.add(comp.name)
                last_var = var

        # --- 残りのshunt素子 (未配置) ---
        remaining_shunts = [comp for comp in shunts if comp.name not in placed]
        for i, comp in enumerate(remaining_shunts):
            var = _make_var(comp.name)
            elem = _get_schemdraw_element(comp)
            label = _make_label(comp)
            if len(remaining_shunts) > 1:
                # Dense same-node source/L/C cases otherwise stack every
                # shunt at the source output.  Walk along a small bus and
                # drop each shunt from its own branch point.
                self.lines.append(f"    d.add(elm.Line().right().length(d.unit*0.85))")
            self.lines.append(f"    d.add(elm.Dot())")
            self.lines.append(f"    d.push()")
            if len(remaining_shunts) == 1:
                self.lines.append(f"    d.add(elm.Line().right().length(d.unit/3))")
            self.lines.append(f"    {var} = d.add({elem}.down().label('{label}', loc='bottom'))")
            self.lines.append(f"    d.add(elm.Ground())")
            self.lines.append(f"    d.pop()")
            placed.add(comp.name)

        # --- 3端子素子 ---
        for comp in three_terms:
            if comp.name in placed:
                continue
            var = _make_var(comp.name)
            elem = _get_schemdraw_element(comp)
            label = _sanitize_label(comp.name)
            c_node, b_node, e_node = _parse_3terminal(comp)

            self.lines.append(f"    d.add(elm.Dot())")
            self.lines.append(f"    d.push()")
            self.lines.append(f"    {var} = d.add({elem}.right().label('{label}'))")
            self.lines.append(f"    d.pop()")
            placed.add(comp.name)

        # --- その他の素子 ---
        for comp in others:
            if comp.name in placed:
                continue
            var = _make_var(comp.name)
            elem = _get_schemdraw_element(comp)
            label = _make_label(comp)
            self.lines.append(f"    {var} = d.add({elem}.right().label('{label}'))")
            placed.add(comp.name)

        # --- GND ---
        # Place Ground at source's start (negative terminal)
        # Do NOT draw Lines to connect GND bus — that merges signal nodes
        # in the reverse (schemdraw→.cir) conversion.
        if first_source_var and first_source_var != '_start':
            self.lines.append(f"    d.add(elm.Ground().at({first_source_var}.start))")

        # --- ノードラベル (named nodes only) ---
        for node in self.node_order:
            nl = node.lower()
            # Skip purely numeric or generic node names
            if node.isdigit() or nl.startswith('n0') or nl.startswith('n00'):
                continue
            # Skip LTSpice auto-generated names like N001
            if re.match(r'^[Nn]\d+$', node):
                continue
            safe_node = _sanitize_label(node)
            self.lines.append(f"    # Node: {safe_node}")

    def _emit_shunts_at_node(self, node: str, shunts: List[Component],
                              placed: set, parallel_groups: dict):
        """Emit shunt components at a given node, handling parallel shunts with push/pop."""
        # Find shunts at this node
        node_shunts = [s for s in shunts if s.name not in placed and
                       (s.node_pos == node or s.node_neg == node)]
        if not node_shunts:
            return

        # Check for parallel shunt groups
        if len(node_shunts) >= 2:
            # Multiple shunts at same node - check if they share same node pair
            pair_key = _node_pair_key(node_shunts[0])
            same_pair = [s for s in node_shunts if _node_pair_key(s) == pair_key]

            if len(same_pair) >= 2:
                # Parallel shunts: use push/pop with horizontal offset
                for i, sc in enumerate(same_pair):
                    svar = _make_var(sc.name)
                    selem = _get_schemdraw_element(sc)
                    slabel = _make_label(sc)
                    self.lines.append(f"    d.push()")
                    # Give each shunt a short branch stub before dropping it.
                    # Without this, input shunts share the source x-coordinate
                    # and labels/components visually collide.
                    self.lines.append(f"    d.add(elm.Line().right().length(d.unit/3))")
                    if i > 0:
                        self.lines.append(f"    d.add(elm.Line().right().length(d.unit*0.85))")
                    self.lines.append(f"    {svar} = d.add({selem}.down().label('{slabel}', loc='bottom'))")
                    self.lines.append(f"    d.add(elm.Ground())")
                    self.lines.append(f"    d.pop()")
                    placed.add(sc.name)
                return

        # Single shunt or non-parallel shunts
        for sc in node_shunts:
            svar = _make_var(sc.name)
            selem = _get_schemdraw_element(sc)
            slabel = _make_label(sc)
            self.lines.append(f"    d.push()")
            self.lines.append(f"    d.add(elm.Line().right().length(d.unit/3))")
            self.lines.append(f"    {svar} = d.add({selem}.down().label('{slabel}', loc='bottom'))")
            self.lines.append(f"    d.add(elm.Ground())")
            self.lines.append(f"    d.pop()")
            placed.add(sc.name)


# =============================================================================
# 便利関数
# =============================================================================

def cir_to_schemdraw(cir_path: str, output_path: str = None) -> str:
    """CIRファイルをschemdrawスクリプトに変換"""
    converter = CirToSchemdraw()
    return converter.convert_file(cir_path, output_path)


def cir_string_to_schemdraw(netlist: str, name: str = 'circuit') -> str:
    """ネットリスト文字列をschemdrawスクリプトに変換"""
    converter = CirToSchemdraw()
    return converter.convert_string(netlist, name)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        script = cir_to_schemdraw(sys.argv[1])
        print(script)
    else:
        # テスト: RC ローパスフィルタ
        test_cir = """* RC Lowpass Filter
V1 in 0 AC 1
R1 in out 1k
C1 out 0 1u
.ac dec 20 1 100k
.end"""
        script = cir_string_to_schemdraw(test_cir, 'rc_lowpass')
        print(script)
