#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LTSpice ASC to schemdraw Python script converter

.ascファイルを解析し、同等の回路を描画するschemdrawスクリプトを生成する。
三角変換の最後のピース: .cir ↔ .asc ↔ .py

アルゴリズム:
1. ASCをパースしてコンポーネント・ノード接続を抽出
2. 回路グラフを構築
3. ソースから始めてグラフを辿り、描画順序を決定
4. 各ステップでschemdrawコードを生成
"""

from pathlib import Path
from typing import List, Dict, Set, Tuple, Optional
from dataclasses import dataclass

from .asc_parser import (AscParser, NetlistExtractor, AscSymbol,
                         SYMBOL_TO_SPICE, PASSIVE_SYMBOLS,
                         TWO_TERMINAL_SYMBOLS, THREE_TERMINAL_SYMBOLS)


# =============================================================================
# LTSpice シンボル → schemdraw 要素名マッピング
# =============================================================================

SYMBOL_TO_SCHEMDRAW = {
    # 受動素子
    'res': 'elm.Resistor()',
    'cap': 'elm.Capacitor()',
    'ind': 'elm.Inductor2()',
    'ind2': 'elm.Inductor2()',  # 個別のind2、Kの場合はTransformerに昇格
    'polcap': 'elm.Capacitor(polar=True)',
    # ソース
    'voltage': 'elm.SourceV()',
    'current': 'elm.SourceI()',
    'battery': 'elm.Battery()',
    # ダイオード
    'diode': 'elm.Diode()',
    'schottky': 'elm.Schottky()',
    'zener': 'elm.Zener()',
    'led': 'elm.LED()',
    'varactor': 'elm.Varactor()',
    # トランジスタ
    'npn': 'elm.BjtNpn()',
    'pnp': 'elm.BjtPnp()',
    'nmos': 'elm.NFet()',
    'pmos': 'elm.PFet()',
    'njf': 'elm.JFetN()',
    'pjf': 'elm.JFetP()',
    # 制御ソース
    'e': 'elm.SourceControlledV()',
    'g': 'elm.SourceControlledI()',
    'bv': 'elm.SourceV()',
    # スイッチ
    'sw': 'elm.Switch()',
    # その他
    'xtal': 'elm.Crystal()',
}

# LTSpice回転 → schemdraw方向（受動素子用）
ROTATION_TO_DIRECTION = {
    'R0': 'down',      # LTSpice R0: top→bottom = schemdraw down
    'R90': 'right',    # LTSpice R90: left→right = schemdraw right
    'R180': 'up',      # LTSpice R180: bottom→top = schemdraw up
    'R270': 'left',    # LTSpice R270: right→left = schemdraw left
    'M0': 'down',
    'M90': 'right',
    'M180': 'up',
    'M270': 'left',
}

# ソース素子の回転→方向
# LTSpice voltage R0: +端子が上 → schemdraw SourceV().up() (end=+)
# LTSpice current R180: 矢印上向き → schemdraw SourceI().up()
SOURCE_ROTATION_TO_DIRECTION = {
    'voltage': {
        'R0': 'up', 'R90': 'right', 'R180': 'down', 'R270': 'left',
        'M0': 'up', 'M90': 'right', 'M180': 'down', 'M270': 'left',
    },
    'current': {
        'R0': 'down', 'R90': 'left', 'R180': 'up', 'R270': 'right',
        'M0': 'down', 'M90': 'left', 'M180': 'up', 'M270': 'right',
    },
    'battery': {
        'R0': 'up', 'R90': 'right', 'R180': 'down', 'R270': 'left',
        'M0': 'up', 'M90': 'right', 'M180': 'down', 'M270': 'left',
    },
}


# =============================================================================
# schemdrawスクリプト生成器
# =============================================================================

@dataclass
class CircuitNode:
    """回路グラフのノード"""
    name: str
    x: int = 0
    y: int = 0
    is_ground: bool = False
    components: List[str] = None  # 接続されたコンポーネントのInstName

    def __post_init__(self):
        if self.components is None:
            self.components = []


@dataclass
class CircuitEdge:
    """回路グラフのエッジ（コンポーネント）"""
    inst_name: str
    symbol: AscSymbol
    node1: str   # 端子1のノード名
    node2: str   # 端子2のノード名


class AscToSchemdraw:
    """ASCファイルからschemdrawスクリプトを生成"""

    def __init__(self):
        self.asc: Optional[AscParser] = None
        self.nodes: Dict[str, CircuitNode] = {}
        self.edges: List[CircuitEdge] = []
        self.code_lines: List[str] = []

    def convert_file(self, asc_path: str, output_path: str = None) -> str:
        """ASCファイルをschemdrawスクリプトに変換"""
        self.asc = AscParser()
        self.asc.parse_file(asc_path)

        script = self._generate_script(Path(asc_path).stem)

        if output_path is None:
            output_path = str(Path(asc_path).with_suffix('.gen.py'))

        Path(output_path).write_text(script, encoding='utf-8')
        print(f'Script saved: {output_path}')
        return script

    def convert_string(self, asc_content: str, name: str = 'circuit') -> str:
        """ASC文字列をschemdrawスクリプトに変換"""
        self.asc = AscParser()
        self.asc.parse_string(asc_content)
        return self._generate_script(name)

    def _generate_script(self, name: str) -> str:
        """schemdrawスクリプトを生成"""
        # Step 1: ネット接続を解析
        extractor = NetlistExtractor(self.asc)
        extractor.extract()

        # Step 2: 回路グラフを構築
        self._build_circuit_graph(extractor)

        # Step 3: 描画順序を決定してコード生成
        self.code_lines = []
        self._emit_header(name)
        self._emit_circuit()
        self._emit_footer(name)

        return '\n'.join(self.code_lines)

    def _build_circuit_graph(self, extractor: NetlistExtractor):
        """回路グラフを構築"""
        self.nodes = {}
        self.edges = []

        # ノード情報を収集
        for flag in self.asc.flags:
            node_name = extractor._get_node_at((flag.x, flag.y))
            is_gnd = flag.name == '0'
            self.nodes[node_name] = CircuitNode(
                name=node_name, x=flag.x, y=flag.y, is_ground=is_gnd)

        # エッジ（コンポーネント）情報を収集
        for sym in self.asc.symbols:
            terms = self.asc.get_component_terminals(sym)
            if terms is None or len(terms) < 2:
                continue

            node1 = extractor._get_node_at(terms[0])
            node2 = extractor._get_node_at(terms[1])

            # ノードが未登録なら追加
            for nname, coord in [(node1, terms[0]), (node2, terms[1])]:
                if nname not in self.nodes:
                    self.nodes[nname] = CircuitNode(
                        name=nname, x=coord[0], y=coord[1])

            edge = CircuitEdge(
                inst_name=sym.inst_name,
                symbol=sym,
                node1=node1,
                node2=node2)
            self.edges.append(edge)

            self.nodes[node1].components.append(sym.inst_name)
            self.nodes[node2].components.append(sym.inst_name)

    def _emit_header(self, name: str):
        """ヘッダー部分を生成"""
        self.code_lines.extend([
            '#!/usr/bin/env python3',
            '# -*- coding: utf-8 -*-',
            f'"""Auto-generated schemdraw script from {name}.asc"""',
            '',
            'import schemdraw',
            'import schemdraw.elements as elm',
            '',
            f"with schemdraw.Drawing(show=False) as d:",
            f"    d.config(unit=3, font='Times New Roman')",
            '',
        ])

    def _emit_footer(self, name: str):
        """フッター部分を生成"""
        # ディレクティブをコメントとして追加
        for text in self.asc.texts:
            if text.is_directive:
                self.code_lines.append(
                    f"    # SPICE: {text.text}")

        self.code_lines.extend([
            '',
            f"    d.save('{name}.pdf')",
            f"    print('Saved: {name}.pdf')",
        ])

    def _emit_circuit(self):
        """回路のschemdrawコードを生成

        アルゴリズム:
        1. 全ソースを見つけて、それぞれからグラフを辿る
        2. ソースがない場合は最初のコンポーネントから
        3. 未到達のコンポーネントも救済して配置
        4. GNDへの帰線とGroundを最後に追加
        """
        # ノード→コンポーネントの隣接マップ
        node_to_edges: Dict[str, List[CircuitEdge]] = {}
        for node_name in self.nodes:
            node_to_edges[node_name] = []
        for edge in self.edges:
            node_to_edges.setdefault(edge.node1, []).append(edge)
            node_to_edges.setdefault(edge.node2, []).append(edge)

        placed: Set[str] = set()
        gnd_names = {n.name for n in self.nodes.values() if n.is_ground}
        first_source_var = None

        # Step 1: 全ソースを見つける
        source_edges = [e for e in self.edges
                        if e.symbol.symbol_type in ('voltage', 'current', 'battery')]

        if not source_edges:
            # ソースなし → 最初のコンポーネントから
            source_edges = self.edges[:1]

        # Step 2: 各ソースから辿る
        for source_edge in source_edges:
            if source_edge.inst_name in placed:
                continue

            source_var = self._emit_component(source_edge, None, None)
            placed.add(source_edge.inst_name)

            if first_source_var is None:
                first_source_var = source_var

            # ソースの信号側ノード（GNDでない方）
            if source_edge.node1 in gnd_names:
                start_node = source_edge.node2
            elif source_edge.node2 in gnd_names:
                start_node = source_edge.node1
            else:
                start_node = source_edge.node1

            self._traverse_from_node(start_node, source_edge, source_var,
                                      'end', placed, node_to_edges, gnd_names)

        # Step 3: 未到達コンポーネントを救済
        for edge in self.edges:
            if edge.inst_name not in placed:
                placed.add(edge.inst_name)
                var = self._emit_component(edge, None, None)
                # この素子の先もDFSで辿る
                for node_name in [edge.node1, edge.node2]:
                    if node_name not in gnd_names:
                        self._traverse_from_node(
                            node_name, edge, var, 'end',
                            placed, node_to_edges, gnd_names)

        # Step 4: GND帰線
        if first_source_var:
            self.code_lines.append(
                f"    d.add(elm.Line().left().to({first_source_var}.start))")
        self.code_lines.append(f"    d.add(elm.Ground())")

        # Step 5: 出力ラベルフラグ
        for flag in self.asc.flags:
            if flag.name != '0' and not flag.name.startswith('n'):
                self.code_lines.append(
                    f"    d.add(elm.Dot())")
                self.code_lines.append(
                    f"    d.add(elm.Label().label('{flag.name}', loc='right'))")

    def _traverse_from_node(self, node_name: str,
                             from_edge: CircuitEdge,
                             from_var: str,
                             from_anchor: str,
                             placed: Set[str],
                             node_to_edges: Dict[str, List[CircuitEdge]],
                             gnd_names: Set[str],
                             depth: int = 0):
        """ノードから出るコンポーネントを再帰的に辿る"""
        if depth > 50:
            return

        # このノードから出る未配置コンポーネントを収集
        available = []
        for edge in node_to_edges.get(node_name, []):
            if edge.inst_name not in placed:
                available.append(edge)

        if not available:
            return

        # 分類: 直列（GNDに行かない）と分路（GNDに行く）
        series = []
        shunt = []
        for edge in available:
            other_node = edge.node2 if edge.node1 == node_name else edge.node1
            if other_node in gnd_names:
                shunt.append(edge)
            else:
                series.append(edge)

        # メインパス: 直列素子を優先（なければ分路）
        # 複数の直列素子がある場合: 最初のものをメインパスに
        main_path = series[:1]
        branches = series[1:] + shunt

        # 分岐がある場合: Dotを追加
        if branches:
            self.code_lines.append(f"    d.add(elm.Dot())")

        # 分岐パスを先にpush/popで処理
        for branch_edge in branches:
            placed.add(branch_edge.inst_name)
            self.code_lines.append(f"    d.push()")

            branch_var = self._emit_component(
                branch_edge, from_var, from_anchor,
                node_name)
            placed.add(branch_edge.inst_name)

            # 分岐の先にさらに素子があれば辿る
            other_node = (branch_edge.node2
                         if branch_edge.node1 == node_name
                         else branch_edge.node1)
            if other_node not in gnd_names:
                self._traverse_from_node(
                    other_node, branch_edge, branch_var, 'end',
                    placed, node_to_edges, gnd_names, depth + 1)

            self.code_lines.append(f"    d.pop()")

        # メインパスを辿る
        for main_edge in main_path:
            placed.add(main_edge.inst_name)
            main_var = self._emit_component(
                main_edge, from_var, from_anchor,
                node_name)

            # メインパスの先のノードから再帰
            other_node = (main_edge.node2
                         if main_edge.node1 == node_name
                         else main_edge.node1)
            if other_node not in gnd_names:
                self._traverse_from_node(
                    other_node, main_edge, main_var, 'end',
                    placed, node_to_edges, gnd_names, depth + 1)

    def _emit_component(self, edge: CircuitEdge,
                         from_var: Optional[str],
                         from_anchor: Optional[str],
                         from_node: Optional[str] = None) -> str:
        """コンポーネント1個のschemdrawコードを出力し、変数名を返す"""
        sym = edge.symbol
        var = self._make_var_name(sym.inst_name)
        elem_code = self._symbol_to_schemdraw(sym)
        direction = self._get_direction(sym)
        label = self._get_label(sym)

        parts = [elem_code]

        # 位置指定
        if from_var and from_anchor:
            parts.append(f'.at({from_var}.{from_anchor})')

        # 方向
        parts.append(f'.{direction}()')

        # ラベル
        if label:
            # 方向に応じてラベル位置を調整
            if direction in ('down', 'up'):
                parts.append(f".label('{label}', loc='left')")
            else:
                parts.append(f".label('{label}')")

        code = ''.join(parts)
        self.code_lines.append(f"    {var} = d.add({code})")
        return var

    def _find_source(self) -> Optional[CircuitEdge]:
        """ソース素子を探す"""
        for edge in self.edges:
            if edge.symbol.symbol_type in ('voltage', 'current', 'battery'):
                return edge
        return None

    def _symbol_to_schemdraw(self, sym: AscSymbol) -> str:
        """LTSpiceシンボル → schemdraw要素コード"""
        return SYMBOL_TO_SCHEMDRAW.get(
            sym.symbol_type, f"elm.Line()  # unknown: {sym.symbol_type}")

    def _get_direction(self, sym: AscSymbol) -> str:
        """LTSpice回転 → schemdraw方向"""
        if sym.symbol_type in SOURCE_ROTATION_TO_DIRECTION:
            return SOURCE_ROTATION_TO_DIRECTION[sym.symbol_type].get(
                sym.rotation, 'up')
        return ROTATION_TO_DIRECTION.get(sym.rotation, 'right')

    def _get_label(self, sym: AscSymbol) -> str:
        """ラベル文字列を生成"""
        parts = []
        if sym.inst_name:
            parts.append(sym.inst_name)
        if sym.value:
            parts.append(sym.value)
        return ' '.join(parts) if parts else ''

    def _make_var_name(self, inst_name: str) -> str:
        """インスタンス名からPython変数名を生成"""
        # R1 → R1, C1 → C1 etc.
        name = inst_name.replace(' ', '_').replace('-', '_')
        # Python予約語チェック
        if name in ('in', 'is', 'as', 'or', 'if', 'for', 'not', 'and'):
            name = name + '_'
        return name


# =============================================================================
# 便利関数
# =============================================================================

def asc_to_schemdraw(asc_path: str, output_path: str = None) -> str:
    """ASCファイルをschemdrawスクリプトに変換"""
    converter = AscToSchemdraw()
    return converter.convert_file(asc_path, output_path)


# =============================================================================
# テスト
# =============================================================================

if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        asc_path = sys.argv[1]
        script = asc_to_schemdraw(asc_path)
        print(script)
    else:
        # テスト: 既存の変換例をschemdrawスクリプトに変換
        test_files = [
            r"examples\00_converter\01_rc_lowpass\test_rc_lowpass.asc",
            r"examples\00_converter\06_voltage_rc\test_voltage_rc.asc",
            r"examples\00_converter\05_parallel_lc\test_parallel_lc.asc",
        ]

        for f in test_files:
            full = str(Path(__file__).parent.parent / f)
            print(f"\n{'='*60}")
            print(f"Converting: {f}")
            print(f"{'='*60}")
            try:
                converter = AscToSchemdraw()
                script = converter.convert_file(full,
                    str(Path(full).with_suffix('.gen.py')))
                print(script)
            except Exception as e:
                print(f"ERROR: {e}")
                import traceback
                traceback.print_exc()
