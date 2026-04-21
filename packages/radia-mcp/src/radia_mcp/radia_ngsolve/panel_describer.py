"""
Static AST parser for Radia ``radia_*.py`` panel files.

Extracts the user-visible widget tree, mode-switch visibility
behaviour, and CLI command builders from the panel source WITHOUT
instantiating Qt or importing NGSolve. The result is fed to two
MCP tools:

  panel_describe_jp(panel_name)
      -- Japanese-language hierarchical description of the panel
         as it CURRENTLY exists in the source tree (not the cached
         panel_registry.json which can be stale).

  panel_widget_locations(panel_name, widget_key)
      -- File path + line numbers for every place a widget is
         defined, toggled visible/hidden, or read by the command
         builder. Lets an LLM make a precise edit without grep.

Why AST instead of runtime introspection: importing PySide6 +
instantiating IHPanel() pulls Qt + every panel module into the
MCP server process, which is heavy and side-effectful. AST
parsing is ~100x faster and has no runtime dependencies.
"""

from __future__ import annotations

import ast
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ============================================================
# Data classes
# ============================================================

@dataclass
class WidgetDef:
    """One add_line / add_combo / add_spin / add_browse call."""
    key: str
    label: str
    kind: str             # "line" | "combo" | "spin" | "browse"
    default: str = ""
    items: list = field(default_factory=list)  # combo items
    spin_range: tuple = ()                     # (min, max) for spin
    line_no: int = 0


@dataclass
class VisibilityRule:
    """One ``self._set_row_visible("key", expr)`` call inside a
    handler. ``expr`` is the literal source text of the second
    argument (e.g. ``"is_bem"``, ``"True"``, ``"not is_esim"``)."""
    key: str
    expr: str
    line_no: int


@dataclass
class CmdArg:
    """One CLI flag emitted by a _build_*_command method."""
    flag: str             # e.g. "--frequency"
    source: str           # e.g. 'self.val("freq")' or '"custom"'
    line_no: int


@dataclass
class HandlerInfo:
    """A method like _on_method_changed or _on_workpiece_changed."""
    name: str
    visibility_rules: list = field(default_factory=list)  # VisibilityRule
    branches: list = field(default_factory=list)          # ("is_bem", [rules])


@dataclass
class CommandBuilder:
    """A method like _build_bem_command / _build_fem_command."""
    name: str
    script: str = ""        # e.g. "calc_inductance.py"
    args: list = field(default_factory=list)  # CmdArg


@dataclass
class PanelInfo:
    panel_name: str
    panel_class: str
    file_path: str
    method_combo_items: list = field(default_factory=list)
    widgets: list = field(default_factory=list)         # WidgetDef
    handlers: list = field(default_factory=list)        # HandlerInfo
    builders: list = field(default_factory=list)        # CommandBuilder


# ============================================================
# AST helpers
# ============================================================

def _str_or_none(node: ast.expr) -> Optional[str]:
    """Return the string literal value if node is a Constant str."""
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def _list_of_strings(node: ast.expr) -> list:
    """Return [str, ...] if node is a list/tuple of string literals."""
    if not isinstance(node, (ast.List, ast.Tuple)):
        return []
    out = []
    for elt in node.elts:
        s = _str_or_none(elt)
        if s is not None:
            out.append(s)
    return out


def _unparse(node: ast.expr) -> str:
    try:
        return ast.unparse(node)
    except Exception:
        return "<expr>"


# ============================================================
# Visitors
# ============================================================

class _BuildUiVisitor(ast.NodeVisitor):
    """Walk the body of ``_build_ui`` and collect widget defs +
    the Method combo items."""

    def __init__(self):
        self.widgets: list = []
        self.method_items: list = []

    def visit_Call(self, node: ast.Call):
        # self.add_line(...) / self.add_combo(...) / etc.
        func = node.func
        if (isinstance(func, ast.Attribute)
                and isinstance(func.value, ast.Name)
                and func.value.id == "self"
                and func.attr in ("add_line", "add_combo",
                                   "add_spin", "add_browse")):
            kind = func.attr.split("_", 1)[1]  # "line", "combo", ...
            args = node.args
            if len(args) >= 2:
                key = _str_or_none(args[0])
                label = _str_or_none(args[1])
                if key and label:
                    w = WidgetDef(key=key, label=label, kind=kind,
                                  line_no=node.lineno)
                    if kind == "line" and len(args) >= 3:
                        w.default = _str_or_none(args[2]) or ""
                    elif kind == "combo" and len(args) >= 3:
                        w.items = _list_of_strings(args[2])
                        if key == "method":
                            self.method_items = list(w.items)
                    elif kind == "spin":
                        if len(args) >= 3 and isinstance(
                                args[2], ast.Constant):
                            w.default = repr(args[2].value)
                        if len(args) >= 5:
                            try:
                                w.spin_range = (
                                    args[3].value, args[4].value)
                            except AttributeError:
                                pass
                    self.widgets.append(w)
        self.generic_visit(node)


class _HandlerVisitor(ast.NodeVisitor):
    """Walk a single handler method body and collect:

      - direct ``self._set_row_visible(key, expr)`` calls (always
        applied), and
      - ``if`` branch conditions wrapped around them, so the caller
        can group rules by condition.

    Also picks up ``combo.addItems([...])`` calls so the dynamic
    Solver combo (BEM vs FEM items) is reported.
    """

    def __init__(self):
        self.rules_by_branch: dict = {"<top>": []}
        self.dynamic_combo_items: dict = {}  # widget_var -> [items]

    def _walk_block(self, body: list, branch_label: str):
        for stmt in body:
            if isinstance(stmt, ast.If):
                cond = _unparse(stmt.test)
                # Recurse into both branches with annotated labels
                if cond not in self.rules_by_branch:
                    self.rules_by_branch[cond] = []
                self._walk_block(stmt.body, cond)
                if stmt.orelse:
                    not_cond = f"not ({cond})"
                    if not_cond not in self.rules_by_branch:
                        self.rules_by_branch[not_cond] = []
                    self._walk_block(stmt.orelse, not_cond)
            elif isinstance(stmt, ast.Expr) and isinstance(
                    stmt.value, ast.Call):
                self._handle_call(stmt.value, branch_label)
            else:
                # Walk for nested calls (e.g. ``for ... self._set_row_visible(...)``)
                for child in ast.walk(stmt):
                    if isinstance(child, ast.Call):
                        self._handle_call(child, branch_label)

    def _handle_call(self, call: ast.Call, branch_label: str):
        f = call.func
        if (isinstance(f, ast.Attribute)
                and isinstance(f.value, ast.Name)
                and f.value.id == "self"
                and f.attr == "_set_row_visible"
                and len(call.args) >= 2):
            key = _str_or_none(call.args[0])
            if key:
                rule = VisibilityRule(
                    key=key, expr=_unparse(call.args[1]),
                    line_no=call.lineno)
                self.rules_by_branch.setdefault(
                    branch_label, []).append(rule)
        # combo.addItems([...]) inside a branch — captures the
        # dynamic Solver combo populated by _on_method_changed
        if (isinstance(f, ast.Attribute) and f.attr == "addItems"
                and len(call.args) == 1):
            items = _list_of_strings(call.args[0])
            if items:
                target = _unparse(f.value)
                self.dynamic_combo_items[
                    f"{target} ({branch_label})"] = items

    def visit_FunctionDef(self, node: ast.FunctionDef):
        self._walk_block(node.body, "<top>")
        # Don't recurse further — we're called per-handler.


class _BuilderVisitor(ast.NodeVisitor):
    """Walk a _build_*_command method and collect the CLI args."""

    def __init__(self):
        self.script: str = ""
        self.args: list = []

    def visit_Call(self, node: ast.Call):
        # calc_script("calc_FOO.py") -> capture FOO
        if (isinstance(node.func, ast.Name)
                and node.func.id == "calc_script"
                and len(node.args) == 1):
            s = _str_or_none(node.args[0])
            if s:
                self.script = s
        self.generic_visit(node)

    def collect_args(self, node: ast.FunctionDef):
        """Two patterns:

          cmd = [_PYTHON, calc_script("..."), "--vol", vol_path,
                 "--frequency", self.val("freq"), ...]
          cmd += ["--coil-sigma", coil_sigma]

        Collect every (--flag, source) pair from both."""
        for stmt in ast.walk(node):
            # cmd = [...] initial assignment
            if isinstance(stmt, ast.Assign):
                if (len(stmt.targets) == 1
                        and isinstance(stmt.targets[0], ast.Name)
                        and stmt.targets[0].id == "cmd"
                        and isinstance(stmt.value, ast.List)):
                    self._extract_pairs(stmt.value.elts, stmt.lineno)
            # cmd += [...] augmented assignment
            elif isinstance(stmt, ast.AugAssign):
                if (isinstance(stmt.target, ast.Name)
                        and stmt.target.id == "cmd"
                        and isinstance(stmt.op, ast.Add)
                        and isinstance(stmt.value, ast.List)):
                    self._extract_pairs(stmt.value.elts, stmt.lineno)

    def _extract_pairs(self, elts: list, lineno: int):
        i = 0
        while i < len(elts):
            elt = elts[i]
            flag = _str_or_none(elt)
            if flag and flag.startswith("--"):
                if i + 1 < len(elts):
                    src = _unparse(elts[i + 1])
                    self.args.append(CmdArg(flag=flag, source=src,
                                            line_no=lineno))
                    i += 2
                    continue
            i += 1


# ============================================================
# Public API
# ============================================================

def parse_panel_file(path: str | Path) -> PanelInfo:
    """Parse a ``radia_*.py`` panel file and return a PanelInfo."""
    path = Path(path)
    src = path.read_text(encoding="utf-8")
    tree = ast.parse(src, filename=str(path))

    panel_class = ""
    panel_name = path.stem.replace("radia_", "")
    info = PanelInfo(
        panel_name=panel_name, panel_class="", file_path=str(path))

    for node in ast.walk(tree):
        if (isinstance(node, ast.ClassDef)
                and node.name.endswith("Panel")):
            panel_class = node.name
            info.panel_class = panel_class

            for item in node.body:
                if not isinstance(item, ast.FunctionDef):
                    continue
                if item.name == "_build_ui":
                    bv = _BuildUiVisitor()
                    bv.visit(item)
                    info.widgets.extend(bv.widgets)
                    info.method_combo_items = bv.method_items
                elif item.name.startswith("_on_") and item.name.endswith(
                        "_changed"):
                    hv = _HandlerVisitor()
                    hv.visit_FunctionDef(item)
                    h = HandlerInfo(name=item.name)
                    for branch, rules in hv.rules_by_branch.items():
                        h.branches.append((branch, rules))
                    info.handlers.append(h)
                elif item.name.startswith("_build_") and item.name.endswith(
                        "_command"):
                    cv = _BuilderVisitor()
                    cv.visit(item)
                    cv.collect_args(item)
                    info.builders.append(CommandBuilder(
                        name=item.name, script=cv.script, args=cv.args))
            break  # only the first *Panel class

    return info


# ============================================================
# Japanese formatter
# ============================================================

# Japanese labels for common widget keys (extend as needed)
_JA_LABELS = {
    "method": "解析手法",
    "fes_order": "FE 次数",
    "solver": "ソルバ",
    "max_iter": "最大反復回数",
    "freq": "周波数 [Hz]",
    "current": "コイル電流 [A]",
    "coil_sigma": "コイル導電率 [S/m]",
    "workpiece_mode": "ワークピース",
    "wp_sigma": "WP 導電率 [S/m]",
    "half_thickness": "半厚さ [m]",
    "mu_r": "比透磁率 mu_r",
    "bh_file": "BH カーブファイル",
    "esim_geometry": "ESIM 形状座標系",
    "air_mode": "空気領域 B 場計算",
}


def _ja_label(key: str, fallback: str) -> str:
    return _JA_LABELS.get(key, fallback)


def describe_panel_jp(info: PanelInfo) -> str:
    """Format a PanelInfo as a Japanese hierarchical description."""
    lines = []
    lines.append(f"# パネル: {info.panel_class} ({info.panel_name})")
    lines.append(f"ファイル: `{info.file_path}`")
    lines.append("")

    # ---- Widget tree ----
    lines.append("## ウィジェット定義 (現在のソース)")
    lines.append("")
    if info.method_combo_items:
        lines.append(
            f"**Method combo**: {' / '.join(info.method_combo_items)}")
        lines.append("")
    lines.append("| Key | 日本語ラベル | 種類 | デフォルト | 詳細 |")
    lines.append("|-----|-------------|------|-----------|------|")
    for w in info.widgets:
        ja = _ja_label(w.key, w.label.rstrip(":"))
        detail = ""
        if w.kind == "combo" and w.items:
            detail = "/".join(w.items)
        elif w.kind == "spin" and w.spin_range:
            detail = f"範囲 {w.spin_range[0]}..{w.spin_range[1]}"
        elif w.kind == "browse":
            detail = "ファイル選択"
        default = w.default or "—"
        lines.append(f"| `{w.key}` | {ja} | {w.kind} | "
                     f"{default} | {detail} |")
    lines.append("")

    # ---- Handlers ----
    if info.handlers:
        lines.append("## モード切替・可視性ロジック")
        lines.append("")
        for h in info.handlers:
            jp_name = {"_on_method_changed": "Method 変更時",
                       "_on_workpiece_changed": "Workpiece 変更時",
                       "_on_impedance_changed": "Impedance 変更時",
                       "_on_validation_changed": "入力検証",
                       "_on_fem_material_changed": "FEM 材料変更時",
                       }.get(h.name, h.name)
            lines.append(f"### {jp_name} (`{h.name}`)")
            lines.append("")
            for branch, rules in h.branches:
                if not rules:
                    continue
                if branch == "<top>":
                    lines.append("**常時:**")
                else:
                    lines.append(f"**条件 `{branch}`:**")
                for r in rules:
                    ja = _ja_label(r.key, r.key)
                    lines.append(
                        f"- `{r.key}` ({ja}) → "
                        f"visible = `{r.expr}` (line {r.line_no})")
                lines.append("")

    # ---- Command builders ----
    if info.builders:
        lines.append("## サブプロセス起動コマンド")
        lines.append("")
        for b in info.builders:
            jp = b.name.replace("_build_", "").replace("_command", "")
            lines.append(f"### {jp.upper()} モード (`{b.name}`)")
            lines.append("")
            if b.script:
                lines.append(f"スクリプト: `panels/{b.script}`")
                lines.append("")
            if b.args:
                lines.append("CLI 引数:")
                lines.append("")
                lines.append("| フラグ | 値の生成元 |")
                lines.append("|--------|-----------|")
                for a in b.args:
                    lines.append(f"| `{a.flag}` | `{a.source}` |")
                lines.append("")

    return "\n".join(lines)


def widget_locations(info: PanelInfo, widget_key: str) -> dict:
    """Return file:line locations for everything that touches
    ``widget_key``: definition, visibility rules, command builder
    args. The result is a dict ready to JSON-serialize."""
    out = {
        "panel": info.panel_name,
        "file": info.file_path,
        "widget_key": widget_key,
        "definition": None,
        "visibility_rules": [],
        "command_builder_uses": [],
        "japanese_label": _JA_LABELS.get(widget_key, ""),
    }
    for w in info.widgets:
        if w.key == widget_key:
            out["definition"] = {
                "kind": w.kind,
                "label": w.label,
                "default": w.default,
                "items": w.items,
                "line": w.line_no,
            }
            break
    for h in info.handlers:
        for branch, rules in h.branches:
            for r in rules:
                if r.key == widget_key:
                    out["visibility_rules"].append({
                        "handler": h.name,
                        "branch": branch,
                        "expr": r.expr,
                        "line": r.line_no,
                    })
    for b in info.builders:
        for a in b.args:
            if widget_key in a.source or f'"{widget_key}"' in a.source:
                out["command_builder_uses"].append({
                    "builder": b.name,
                    "flag": a.flag,
                    "source": a.source,
                    "line": a.line_no,
                })
    return out


def find_panel_file(panel_name: str) -> Optional[str]:
    """Resolve "ih" / "em" / "pcb" -> radia_<name>.py path.

    panel_describer.py lives at
    ``src/radia/mcp_server/radia_ngsolve/panel_describer.py``,
    so parents[2] = ``src/radia`` (the panels directory). The
    ``radia_<name>.py`` files sit one level above this module's
    grandparent (mcp_server/), in src/radia/.
    """
    here = Path(__file__).resolve()
    # parents: [0]=radia_ngsolve [1]=mcp_server [2]=radia [3]=src
    radia_dir = here.parents[2]
    candidate = radia_dir / f"radia_{panel_name}.py"
    if candidate.is_file():
        return str(candidate)
    return None
