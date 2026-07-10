#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LTSpice ASC Parser & Netlist Extractor

.ascファイルを解析し、SPICEネットリスト(.cir)を生成する。
ラウンドトリップ試験（.cir → .asc → .cir）の逆変換側。
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Set, Optional
from pathlib import Path
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
import math


# =============================================================================
# LTspice executable detection
# =============================================================================

def _find_ltspice_exe() -> Optional[str]:
    """Locate LTspice.exe for use as a canonical netlist generator.

    Order:
    1. env LTSPICE_EXE (explicit override)
    2. PATH (LTspice.exe / XVIIx64.exe / scad3.exe)
    3. Standard install locations on Windows.
    """
    env = os.environ.get("LTSPICE_EXE")
    if env and Path(env).is_file():
        return env

    for cand in ("LTspice.exe", "XVIIx64.exe", "scad3.exe"):
        found = shutil.which(cand)
        if found:
            return found

    candidates = [
        Path(os.environ.get("PROGRAMFILES", "")) / "ADI" / "LTspice" / "LTspice.exe",
        Path(os.environ.get("LOCALAPPDATA", "")) / "Programs" / "ADI" / "LTspice" / "LTspice.exe",
        Path(os.environ.get("PROGRAMFILES", "")) / "LTC" / "LTspiceXVII" / "XVIIx64.exe",
    ]
    for p in candidates:
        try:
            if p.is_file():
                return str(p)
        except Exception:
            pass
    return None


def _dismiss_ltspice_dialogs(pid: int) -> int:
    """Dismiss visible modal dialogs owned by process ``pid``,
    non-destructively. Windows-only; a no-op (returns 0) elsewhere.

    ``LTspice.exe -netlist`` is normally headless, but some schematics
    make it pop a blocking message box -- most commonly *"Duplicate
    overlapping components found: C1 and C2 ..Delete?"* (Yes/No/Cancel).
    Left alone it blocks the whole call until timeout **and** flashes on
    the user's screen. We click **No** (keep the components, netlist the
    schematic exactly as drawn -- never let LTspice silently modify it);
    if there is no No button we fall back to Cancel, then to the
    standard IDNO / IDCANCEL dialog commands.

    Returns the number of dialogs dismissed.
    """
    import sys
    if sys.platform != "win32":
        return 0
    import ctypes
    from ctypes import wintypes

    u = ctypes.windll.user32
    WNDENUMPROC = ctypes.WINFUNCTYPE(wintypes.BOOL, wintypes.HWND,
                                     wintypes.LPARAM)
    BM_CLICK = 0x00F5
    WM_COMMAND = 0x0111
    IDNO, IDCANCEL = 7, 2

    def _text(hwnd) -> str:
        buf = ctypes.create_unicode_buffer(256)
        u.GetWindowTextW(hwnd, buf, 256)
        return buf.value

    def _cls(hwnd) -> str:
        buf = ctypes.create_unicode_buffer(64)
        u.GetClassNameW(hwnd, buf, 64)
        return buf.value

    dialogs: List[int] = []

    def _top(hwnd, _):
        if not u.IsWindowVisible(hwnd):
            return True
        if _cls(hwnd) != "#32770":          # standard dialog/message-box class
            return True
        wpid = wintypes.DWORD()
        u.GetWindowThreadProcessId(hwnd, ctypes.byref(wpid))
        if wpid.value == pid:
            dialogs.append(hwnd)
        return True

    u.EnumWindows(WNDENUMPROC(_top), 0)

    dismissed = 0
    for dlg in dialogs:
        buttons: List[tuple] = []

        def _child(hwnd, _):
            if _cls(hwnd).lower() == "button":
                buttons.append((hwnd, _text(hwnd)))
            return True

        u.EnumChildWindows(dlg, WNDENUMPROC(_child), 0)

        def _norm(s: str) -> str:
            return s.replace("&", "").strip().lower()

        clicked = False
        for wants in (("no", "いいえ"), ("cancel", "キャンセル")):
            for hwnd, txt in buttons:
                if _norm(txt) in wants:
                    u.SendMessageW(hwnd, BM_CLICK, 0, 0)
                    clicked = True
                    break
            if clicked:
                break
        if not clicked:
            # No recognisable safe button -- send the standard dialog
            # commands (No, then Cancel) so a stuck dialog still closes.
            u.SendMessageW(dlg, WM_COMMAND, IDNO, 0)
            u.SendMessageW(dlg, WM_COMMAND, IDCANCEL, 0)
        dismissed += 1
    return dismissed


def _kill_proc_tree(proc) -> None:
    """Terminate a process and its children (so no orphaned LTspice
    window/dialog is left behind on timeout)."""
    import sys
    try:
        if sys.platform == "win32":
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                           capture_output=True)
        else:
            proc.kill()
    except Exception:
        try:
            proc.kill()
        except Exception:
            pass


def _ltspice_emit_netlist(asc_path: str, ltspice_exe: str,
                          timeout: float = 20.0) -> str:
    """Use LTspice to convert .asc to canonical .net (LTspice's own naming).

    `LTspice.exe -netlist input.asc` writes `input.net` next to the .asc.
    We copy the .asc into a temp dir to avoid polluting the source tree.

    A background watchdog dismisses any LTspice modal dialog that pops up
    during the run (e.g. the "duplicate overlapping components" prompt),
    so the call stays headless and never blocks on -- or flashes -- a
    Yes/No box. On timeout the whole process tree is killed so nothing is
    left orphaned. Returns the .net text on success; raises RuntimeError
    on failure (the caller falls back to the pure-Python extractor).
    """
    import threading

    asc = Path(asc_path).resolve()
    if not asc.exists():
        raise FileNotFoundError(asc)
    with tempfile.TemporaryDirectory() as d:
        work = Path(d)
        asc_copy = work / asc.name
        # .asc files are commonly written in cp1252/latin-1; preserve byte-fidelity.
        asc_copy.write_bytes(asc.read_bytes())

        proc = subprocess.Popen(
            [ltspice_exe, "-netlist", str(asc_copy)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )

        stop = threading.Event()

        def _watch():
            # Poll for and dismiss LTspice modal dialogs while it runs.
            while not stop.is_set():
                try:
                    _dismiss_ltspice_dialogs(proc.pid)
                except Exception:
                    pass
                stop.wait(0.25)

        watcher = threading.Thread(target=_watch, daemon=True)
        watcher.start()
        try:
            proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            _kill_proc_tree(proc)
            try:
                proc.communicate(timeout=5)
            except Exception:
                pass
            raise RuntimeError(f"LTspice -netlist timed out after {timeout}s")
        finally:
            stop.set()
            watcher.join(timeout=1.0)

        net = asc_copy.with_suffix(".net")
        if not net.exists():
            # LTspice sometimes returns nonzero even on success, so the
            # .net file (not the exit code) is the real success signal.
            raise RuntimeError(
                f"LTspice -netlist did not produce {net.name} "
                f"(rc={proc.returncode})"
            )
        # Decode the .net robustly. LTspice 26 writes UTF-8 (`µ` as
        # 0xC2 0xB5); older versions wrote cp1252/latin-1 (a bare 0xB5),
        # and some emit UTF-16. The shared decoder handles all of them in
        # one place (a stray latin-1 read once corrupted `250µ` into
        # `250Âµ`, which LTspice then parsed as 250 *seconds*, dilating
        # the .raw time axis 10⁶×).
        from .._shared import decode_asc_bytes
        return decode_asc_bytes(net.read_bytes())


# =============================================================================
# ASCデータ構造
# =============================================================================

@dataclass
class AscWire:
    """WIRE要素"""
    x1: int
    y1: int
    x2: int
    y2: int


@dataclass
class AscFlag:
    """FLAG要素（ネットラベル・GND）"""
    x: int
    y: int
    name: str   # '0' = GND, otherwise net label


@dataclass
class AscSymbol:
    """SYMBOL要素（コンポーネント）"""
    symbol_type: str    # res, cap, ind, voltage, current, npn, etc.
    x: int
    y: int
    rotation: str       # R0, R90, R180, R270, M0, ...
    inst_name: str = ''
    value: str = ''
    value2: str = ''
    spice_model: str = ''
    spice_line: str = ''
    windows: List[str] = field(default_factory=list)
    symbol_path: str = ''  # full path from SYMBOL line (e.g. 'Opamps\\opamp')


@dataclass
class AscText:
    """TEXT要素（ディレクティブ・コメント）"""
    x: int
    y: int
    text: str
    is_directive: bool = False  # ! prefix


# =============================================================================
# 共通端子オフセットテンプレート
# =============================================================================
# .asy シンボルファイルの PIN 座標を基準に、回転/ミラー変換で全方向を生成。
# 変換行列 (実測検証済み):
#   R0:  (x, y)     R90: (-y, x)    R180: (-x, -y)   R270: (y, -x)
#   M0:  (-x, y)    M90: (y, x)     M180: (x, -y)    M270: (-y, -x)


def _make_offsets(pins):
    """R0のピン座標タプルから全8回転のオフセット辞書を生成"""
    tf = {
        'R0':   lambda x, y: (x, y),
        'R90':  lambda x, y: (-y, x),
        'R180': lambda x, y: (-x, -y),
        'R270': lambda x, y: (y, -x),
        'M0':   lambda x, y: (-x, y),
        'M90':  lambda x, y: (y, x),
        'M180': lambda x, y: (x, -y),
        'M270': lambda x, y: (-y, -x),
    }
    return {rot: tuple(f(px, py) for px, py in pins)
            for rot, f in tf.items()}


# .asy PIN 座標 (R0方向)
# res.asy:     PIN 16 16, PIN 16 96
# cap.asy:     PIN 16 0,  PIN 16 64
# ind.asy:     PIN 16 16, PIN 16 96
# ind2.asy:    PIN 16 16, PIN 16 96
# polcap.asy:  PIN 16 0,  PIN 16 64
# voltage.asy: PIN 0 16,  PIN 0 96
# current.asy: PIN 0 0,   PIN 0 80
# battery.asy: PIN 0 16,  PIN 0 96
# diode.asy:   PIN 16 0,  PIN 16 64
# npn.asy:     PIN 64 0,  PIN 0 48,  PIN 64 96  (C, B, E)
# pnp.asy:     PIN 64 0,  PIN 0 48,  PIN 64 96  (C, B, E)
# nmos.asy:    PIN 48 0,  PIN 0 80,  PIN 48 96  (D, G, S)
# pmos.asy:    PIN 48 0,  PIN 0 80,  PIN 48 96  (D, G, S)

_2T_RES = _make_offsets(((16, 16), (16, 96)))
_2T_CAP = _make_offsets(((16, 0), (16, 64)))
_2T_IND = _make_offsets(((16, 16), (16, 96)))
_2T_VOLTAGE = _make_offsets(((0, 16), (0, 96)))
_2T_CURRENT = _make_offsets(((0, 0), (0, 80)))
_2T_BATTERY = _make_offsets(((0, 16), (0, 96)))
_2T_DIODE = _make_offsets(((16, 0), (16, 64)))
_2T_G = _make_offsets(((0, 96), (0, 16)))  # VCCS: pin1(out+)=下, pin2(out-)=上

# opamp.asy: 3-pin (inv, noninv, out)
# PIN -32 48 (invin, SpiceOrder 1), PIN -32 80 (noninvin, SpiceOrder 2), PIN 32 64 (out, SpiceOrder 3)
_3T_OPAMP = _make_offsets(((-32, 48), (-32, 80), (32, 64)))

# UniversalOpAmp*.asy: 5-pin (In+, In-, V+, V-, OUT) — 全バリアント共通
# PIN -32 16 (In+), PIN -32 -16 (In-), PIN 0 -32 (V+), PIN 0 32 (V-), PIN 32 0 (OUT)
_5T_UOPAMP = _make_offsets(((-32, 16), (-32, -16), (0, -32), (0, 32), (32, 0)))

# LT1001.asy 等の5ピンopamp（オフセットが異なる）
# PIN -32 80 (In+), PIN -32 48 (In-), PIN 0 32 (V+), PIN 0 96 (V-), PIN 32 64 (OUT)
_5T_LT1001 = _make_offsets(((-32, 80), (-32, 48), (0, 32), (0, 96), (32, 64)))

_3T_NPN = _make_offsets(((64, 0), (0, 48), (64, 96)))
_3T_PNP = _make_offsets(((64, 0), (0, 48), (64, 96)))
_3T_NMOS = _make_offsets(((48, 0), (0, 80), (48, 96)))
_3T_PMOS = _make_offsets(((48, 0), (0, 80), (48, 96)))
_3T_NJF = _make_offsets(((48, 0), (0, 64), (48, 96)))  # njf.asy実測値
_3T_PJF = _make_offsets(((48, 0), (0, 64), (48, 96)))   # pjf.asy実測値

# アンカーポイントオフセット（シンボル座標 → 端子座標）
TERMINAL_OFFSETS = {
    # 受動素子 (2端子)
    'res':      _2T_RES,
    'cap':      _2T_CAP,
    'ind':      _2T_IND,
    'ind2':     _2T_IND,
    'polcap':   _2T_CAP,

    # ソース (2端子)
    'voltage':  _2T_VOLTAGE,
    'current':  _2T_CURRENT,
    'battery':  _2T_BATTERY,

    # ダイオード系 (2端子)
    'diode':    _2T_DIODE,
    'schottky': _2T_DIODE,
    'zener':    _2T_DIODE,
    'led':      _2T_DIODE,
    'varactor': _2T_DIODE,
    'tvs':      _2T_DIODE,

    # スイッチ (2端子 + 制御)
    'sw':       _2T_CAP,

    # 制御ソース (2端子出力 + 2端子制御) — 出力端子のみ
    'e':        _2T_VOLTAGE,    # VCVS
    'e2':       _2T_VOLTAGE,
    'g':        _2T_G,          # VCCS (pin order reversed vs voltage)
    'g2':       _2T_G,
    'f':        _2T_VOLTAGE,    # CCCS
    'h':        _2T_VOLTAGE,    # CCVS
    'bv':       _2T_VOLTAGE,    # Behavioral voltage
    'bi':       _2T_CURRENT,    # Behavioral current
    'bi2':      _2T_CURRENT,

    # クリスタル (2端子)
    'xtal':     _2T_CAP,

    # 伝送線路 (4端子→2端子近似)
    'tline':    _2T_CAP,
}

# 4端子素子のオフセット（出力2端子 + 制御2端子）
# e.asy:  PIN 0 16, PIN 0 96, PIN -48 32, PIN -48 80  (out+, out-, ctrl+, ctrl-)
# g.asy:  PIN 0 96, PIN 0 16, PIN -48 32, PIN -48 80  (out-, out+, ctrl+, ctrl-)
# sw.asy: PIN 0 16, PIN 0 96, PIN -48 80, PIN -48 32  (out+, out-, ctrl-, ctrl+)
# tline:  PIN -48 -16, PIN -48 16, PIN 48 -16, PIN 48 16
_4T_E = _make_offsets(((0, 16), (0, 96), (-48, 32), (-48, 80)))
_4T_E2 = _make_offsets(((0, 16), (0, 96), (-48, 80), (-48, 32)))  # e2: ctrl pins swapped
_4T_G = _make_offsets(((0, 96), (0, 16), (-48, 32), (-48, 80)))
_4T_G2 = _make_offsets(((0, 96), (0, 16), (-48, 80), (-48, 32)))  # g2: ctrl pins swapped
_4T_SW = _make_offsets(((0, 16), (0, 96), (-48, 80), (-48, 32)))
_4T_TLINE = _make_offsets(((-48, -16), (-48, 16), (48, -16), (48, 16)))

TERMINAL_OFFSETS_4 = {
    'e':     _4T_E,
    'e2':    _4T_E2,
    'g':     _4T_G,
    'g2':    _4T_G2,
    'sw':    _4T_SW,
    'tline': _4T_TLINE,
}

# 3端子素子のオフセット（別テーブル: C/D, B/G, E/S の3点）
TERMINAL_OFFSETS_3 = {
    'npn':   _3T_NPN,
    'pnp':   _3T_PNP,
    'nmos':  _3T_NMOS,
    'pmos':  _3T_PMOS,
    'njf':   _3T_NJF,
    'pjf':   _3T_PJF,
    'nigbt': _3T_NPN,
    'opamp': _3T_OPAMP,   # 3-pin opamp (opamp.sub)
}

# 5端子素子のオフセット（opamp: In+, In-, V+, V-, OUT）
TERMINAL_OFFSETS_5 = {
    'universalopamp':  _5T_UOPAMP,
    'universalopamp1': _5T_UOPAMP,
    'universalopamp2': _5T_UOPAMP,
    'universalopamp3': _5T_UOPAMP,
    'universalopamp3a': _5T_UOPAMP,
    'universalopamp3b': _5T_UOPAMP,
    'universalopamp4': _5T_UOPAMP,
    'universalopamp5': _5T_UOPAMP,
    'lt1001':  _5T_LT1001,
    'lt1001a': _5T_LT1001,
    'lt1028':  _5T_LT1001,
    'lt1028a': _5T_LT1001,
}

# シンボルタイプ → SPICEプレフィックス
SYMBOL_TO_SPICE = {
    'res': 'R', 'cap': 'C', 'ind': 'L', 'ind2': 'L',
    'polcap': 'C',
    'voltage': 'V', 'current': 'I', 'battery': 'V',
    'diode': 'D', 'schottky': 'D', 'zener': 'D', 'led': 'D',
    'varactor': 'D', 'tvs': 'D', 'smdiode': 'D',
    'npn': 'Q', 'pnp': 'Q',
    'npn3': 'Q', 'pnp3': 'Q', 'npn4': 'Q', 'pnp4': 'Q',  # substrate variants
    'nmos': 'M', 'pmos': 'M',
    'nmos4': 'M', 'pmos4': 'M',                          # 4-pin (with body)
    'njf': 'J', 'pjf': 'J',
    'njf4': 'J', 'pjf4': 'J',
    'nigbt': 'Z',
    'e': 'E', 'e2': 'E', 'g': 'G', 'g2': 'G',
    'f': 'F', 'h': 'H',
    'bv': 'B', 'bi': 'B', 'bi2': 'B',
    'sw': 'S',
    'xtal': 'X',
    'tline': 'T',
    'opamp': 'U',
    'universalopamp': 'U', 'universalopamp1': 'U', 'universalopamp2': 'U',
    'universalopamp3': 'U', 'universalopamp3a': 'U', 'universalopamp3b': 'U',
    'universalopamp4': 'U', 'universalopamp5': 'U',
    'lt1001': 'U', 'lt1001a': 'U', 'lt1028': 'U', 'lt1028a': 'U',
}

# 2端子素子シンボルタイプ（ラウンドトリップ対応可能なもの）
TWO_TERMINAL_SYMBOLS = {
    'res', 'cap', 'ind', 'ind2', 'polcap',
    'voltage', 'current', 'battery',
    'diode', 'schottky', 'zener', 'led', 'varactor', 'tvs',
    'sw', 'xtal', 'tline',
    'e', 'e2', 'g', 'g2', 'f', 'h', 'bv', 'bi', 'bi2',
}

# 3端子素子シンボルタイプ
THREE_TERMINAL_SYMBOLS = {
    'npn', 'pnp', 'nmos', 'pmos', 'njf', 'pjf', 'nigbt',
    'npn3', 'pnp3', 'npn4', 'pnp4',
    'nmos4', 'pmos4', 'njf4', 'pjf4',
}

# 受動素子のみ（旧定義 — 後方互換）
PASSIVE_SYMBOLS = {'res', 'cap', 'ind', 'ind2', 'voltage', 'current'}

# Zero-ohm "short" symbols: net-ties whose two pins LTspice merges into a
# single node, dropping the symbol from the netlist. The extractor unions
# their pins (see _build_net_groups) and emits no component for them.
_SHORT_SYMBOLS = {'jumper'}

# Single-letter SPICE element classes we accept from a symbol's
# ``SYMATTR Prefix`` when the symbol is not in the hardcoded SYMBOL_TO_SPICE
# map. Excludes A/K/U (special-function / coupling / opamp) which have their
# own dedicated handling.
_SPICE_CLASS_LETTERS = set('RCLVIDQMJEFGHBSTX')

# ラウンドトリップ可能な全シンボルタイプ
SUPPORTED_SYMBOLS = TWO_TERMINAL_SYMBOLS | THREE_TERMINAL_SYMBOLS


# =============================================================================
# LTspice標準ライブラリ検索パス
# =============================================================================

# LTspice lib.zip の場所候補
_LTSPICE_LIB_ZIP_CANDIDATES = [
    Path(r'C:/Program Files/ADI/LTspice/lib.zip'),
    Path(r'C:/Program Files (x86)/ADI/LTspice/lib.zip'),
]

def _find_ltspice_lib_zip() -> Optional[Path]:
    """LTspice lib.zip を探す"""
    import os
    for p in _LTSPICE_LIB_ZIP_CANDIDATES:
        if p.is_file():
            return p
    # LOCALAPPDATA
    local = os.environ.get('LOCALAPPDATA', '')
    if local:
        p = Path(local) / 'Programs' / 'ADI' / 'LTspice' / 'lib.zip'
        if p.is_file():
            return p
    return None


# =============================================================================
# .asy ファイルパーサー (Improvement A)
# =============================================================================

@dataclass
class AsyPin:
    """PIN definition from .asy file"""
    x: int
    y: int
    direction: str  # LEFT, RIGHT, TOP, BOTTOM, NONE
    name: str = ''
    spice_order: int = 0


class AsyParser:
    """Parse .asy files to extract pin definitions and ordering."""

    # Cache: symbol_path -> list of AsyPin (sorted by SpiceOrder)
    _cache: Dict[str, Optional[List[AsyPin]]] = {}
    _model_cache: Dict[str, str] = {}  # symbol_path -> SpiceModel name
    _attr_cache: Dict[str, str] = {}   # 'symbol_path|attr:Name|dirs' -> value
    _zip_file: Optional[zipfile.ZipFile] = None
    _zip_path: Optional[Path] = None

    @classmethod
    def _get_zip(cls) -> Optional[zipfile.ZipFile]:
        """Get (or open) the LTspice lib.zip"""
        zip_path = _find_ltspice_lib_zip()
        if zip_path is None:
            return None
        if cls._zip_file is None or cls._zip_path != zip_path:
            try:
                cls._zip_file = zipfile.ZipFile(str(zip_path), 'r')
                cls._zip_path = zip_path
            except Exception:
                return None
        return cls._zip_file

    @classmethod
    def parse_asy_text(cls, text: str,
                       _spice_model_out: Optional[List[str]] = None
                       ) -> Optional[List[AsyPin]]:
        """Parse .asy file content, return pins sorted by SpiceOrder.

        If _spice_model_out is provided (a list), appends the SpiceModel name.
        """
        pins: List[AsyPin] = []
        current_pin: Optional[AsyPin] = None

        for line in text.split('\n'):
            line = line.strip()
            if line.startswith('PIN '):
                if current_pin is not None:
                    pins.append(current_pin)
                parts = line.split()
                current_pin = AsyPin(
                    x=int(parts[1]) if len(parts) > 1 else 0,
                    y=int(parts[2]) if len(parts) > 2 else 0,
                    direction=parts[3] if len(parts) > 3 else 'NONE',
                )
            elif line.startswith('PINATTR ') and current_pin is not None:
                parts = line.split(' ', 2)
                if len(parts) >= 3:
                    if parts[1] == 'PinName':
                        current_pin.name = parts[2]
                    elif parts[1] == 'SpiceOrder':
                        try:
                            current_pin.spice_order = int(parts[2])
                        except ValueError:
                            pass
            elif line.startswith('SYMATTR SpiceModel '):
                model = line.split(' ', 2)[2] if len(line.split(' ', 2)) > 2 else ''
                if _spice_model_out is not None and model:
                    _spice_model_out.append(model)

        if current_pin is not None:
            pins.append(current_pin)

        if not pins:
            return None

        # Sort by SpiceOrder (0 means unspecified, put at end)
        pins.sort(key=lambda p: (p.spice_order == 0, p.spice_order))
        return pins

    @classmethod
    def find_and_parse(cls, symbol_path: str,
                       search_dirs: Optional[List[Path]] = None
                       ) -> Optional[List[AsyPin]]:
        """Find and parse .asy for a given symbol path.

        symbol_path: e.g. 'Opamps\\\\opamp', 'res', 'misc\\\\DIAC'
        search_dirs: additional directories to search (e.g. .asc file's dir).
                     If None or empty, falls back to the
                     ``LTSPICE_ASY_SEARCH_PATH`` env var
                     (os.pathsep-separated list).
        """
        # Env-var fallback: read ``LTSPICE_ASY_SEARCH_PATH`` so user-
        # supplied third-party symbol libraries propagate to every
        # call site without explicit plumbing.
        if not search_dirs:
            env = os.environ.get('LTSPICE_ASY_SEARCH_PATH', '')
            if env:
                search_dirs = [Path(d) for d in env.split(os.pathsep) if d.strip()]

        cache_key = symbol_path + '|' + str(search_dirs or [])
        if cache_key in cls._cache:
            return cls._cache[cache_key]

        # Normalize path separators
        norm_path = symbol_path.replace('\\', '/')

        # 1. Search in provided directories
        if search_dirs:
            for d in search_dirs:
                d = Path(d)
                # Try exact path
                asy_file = d / (norm_path + '.asy')
                if asy_file.is_file():
                    result = cls._parse_file(asy_file)
                    if result:
                        cls._cache[cache_key] = result
                        return result
                # Try just the basename
                base = norm_path.split('/')[-1]
                asy_file = d / (base + '.asy')
                if asy_file.is_file():
                    result = cls._parse_file(asy_file)
                    if result:
                        cls._cache[cache_key] = result
                        return result

        # 2. Search in lib.zip
        zf = cls._get_zip()
        if zf is not None:
            # Try lib/sym/<path>.asy
            zip_path = f'lib/sym/{norm_path}.asy'
            try:
                result = cls._parse_zip_entry(zf, zip_path)
                if result:
                    cls._cache[cache_key] = result
                    return result
            except KeyError:
                pass
            # Try case-insensitive search for the basename
            base = norm_path.split('/')[-1].lower()
            for name in zf.namelist():
                if name.lower().endswith(f'/{base}.asy') and name.startswith('lib/sym/'):
                    try:
                        result = cls._parse_zip_entry(zf, name)
                        if result:
                            cls._cache[cache_key] = result
                            return result
                    except Exception:
                        pass

        cls._cache[cache_key] = None
        return None

    @classmethod
    def _parse_zip_entry(cls, zf: zipfile.ZipFile,
                          entry_name: str) -> Optional[List[AsyPin]]:
        """Parse a .asy file from a zip entry, handling encoding."""
        raw = zf.read(entry_name)
        # Try encodings in order
        for enc in ('utf-8', 'utf-16-le', 'utf-16-be', 'ascii', 'latin-1'):
            try:
                text = raw.decode(enc)
                if text.startswith('\ufeff'):
                    text = text[1:]
                # Quick sanity check: should contain 'PIN' if valid
                if 'PIN ' in text:
                    return cls.parse_asy_text(text)
                elif 'Version' in text and 'PIN' not in text:
                    # Valid .asy but no PINs
                    return None
            except (UnicodeDecodeError, UnicodeError):
                continue
        return None

    @classmethod
    def _parse_file(cls, path: Path) -> Optional[List[AsyPin]]:
        """Parse a .asy file from disk."""
        for enc in ('utf-8', 'utf-16-le', 'ascii', 'latin-1'):
            try:
                text = path.read_text(encoding=enc)
                if text.startswith('\ufeff'):
                    text = text[1:]
                return cls.parse_asy_text(text)
            except (UnicodeDecodeError, UnicodeError):
                continue
        return None

    @classmethod
    def get_terminal_offsets(cls, symbol_path: str, rotation: str,
                             search_dirs: Optional[List[Path]] = None
                             ) -> Optional[Tuple[Tuple[int, int], ...]]:
        """Get terminal positions (offsets from symbol origin) in SpiceOrder.

        Returns absolute offsets after applying rotation/mirror transform.
        """
        pins = cls.find_and_parse(symbol_path, search_dirs)
        if pins is None:
            return None

        # Apply rotation/mirror to pin coordinates
        offsets = []
        for pin in pins:
            if pin.spice_order == 0:
                continue  # Skip pins without SpiceOrder
            ox, oy = cls._transform_point(pin.x, pin.y, rotation)
            offsets.append((ox, oy))

        if not offsets:
            # If no SpiceOrder defined, use all pins in order
            for pin in pins:
                ox, oy = cls._transform_point(pin.x, pin.y, rotation)
                offsets.append((ox, oy))

        return tuple(offsets) if offsets else None

    @classmethod
    def get_spice_model(cls, symbol_path: str,
                        search_dirs: Optional[List[Path]] = None) -> str:
        """Get SpiceModel name from .asy file (e.g., 'sample' → 'SAMPLEHOLD')."""
        cache_key = symbol_path + '|model|' + str(search_dirs or [])
        if cache_key in cls._model_cache:
            return cls._model_cache[cache_key]

        result = ''
        norm_path = symbol_path.replace('\\', '/')

        def _extract_model(text: str) -> str:
            for line in text.split('\n'):
                line = line.strip()
                if line.startswith('SYMATTR SpiceModel '):
                    parts = line.split(' ', 2)
                    return parts[2] if len(parts) > 2 else ''
            return ''

        # Search in provided directories
        if search_dirs:
            for d in search_dirs:
                d = Path(d)
                for candidate in [d / (norm_path + '.asy'),
                                  d / (norm_path.split('/')[-1] + '.asy')]:
                    if candidate.is_file():
                        for enc in ('utf-8', 'utf-16-le', 'latin-1'):
                            try:
                                text = candidate.read_text(encoding=enc)
                                result = _extract_model(text)
                                if result:
                                    cls._model_cache[cache_key] = result
                                    return result
                            except (UnicodeDecodeError, UnicodeError):
                                continue

        # Search in lib.zip
        zf = cls._get_zip()
        if zf is not None:
            base = norm_path.split('/')[-1].lower()
            for zname in zf.namelist():
                if zname.lower().endswith(f'/{base}.asy') and zname.startswith('lib/sym/'):
                    try:
                        raw = zf.read(zname)
                        for enc in ('utf-8', 'utf-16-le', 'latin-1'):
                            try:
                                text = raw.decode(enc)
                                result = _extract_model(text)
                                if result:
                                    cls._model_cache[cache_key] = result
                                    return result
                            except (UnicodeDecodeError, UnicodeError):
                                continue
                    except KeyError:
                        pass

        # Search on disk
        ltspice_sym = Path(os.environ.get('LOCALAPPDATA', '')) / 'LTspice' / 'lib' / 'sym'
        if ltspice_sym.exists():
            for candidate in [ltspice_sym / (norm_path + '.asy'),
                              ltspice_sym / (norm_path.split('/')[-1] + '.asy')]:
                if candidate.is_file():
                    for enc in ('utf-8', 'utf-16-le', 'latin-1'):
                        try:
                            text = candidate.read_text(encoding=enc)
                            result = _extract_model(text)
                            if result:
                                cls._model_cache[cache_key] = result
                                return result
                        except (UnicodeDecodeError, UnicodeError):
                            continue

        cls._model_cache[cache_key] = result
        return result

    @classmethod
    def get_symbol_attr(cls, symbol_path: str, attr: str,
                        search_dirs: Optional[List[Path]] = None) -> str:
        """Read a single ``SYMATTR <attr>`` value from the symbol's .asy.

        LTspice symbols declare their SPICE class and default model via
        ``SYMATTR Prefix`` / ``SYMATTR Value`` (e.g. ``xtal`` carries
        ``Prefix C``, ``TowTom2`` carries ``Prefix X`` + ``Value TowTom2``).
        Returns ``''`` when the attribute or the .asy is not found.
        """
        cache_key = f'{symbol_path}|attr:{attr}|{search_dirs or []}'
        if cache_key in cls._attr_cache:
            return cls._attr_cache[cache_key]

        norm_path = (symbol_path or '').replace('\\', '/')
        needle = f'SYMATTR {attr} '

        def _extract(text: str) -> str:
            for line in text.split('\n'):
                line = line.strip()
                if line.startswith(needle):
                    return line[len(needle):].strip()
            return ''

        result = ''
        # 1. provided dirs
        if search_dirs and not result:
            for d in search_dirs:
                d = Path(d)
                for cand in [d / (norm_path + '.asy'),
                             d / (norm_path.split('/')[-1] + '.asy')]:
                    if cand.is_file():
                        for enc in ('utf-8', 'utf-16-le', 'latin-1'):
                            try:
                                result = _extract(cand.read_text(encoding=enc))
                                break
                            except (UnicodeDecodeError, UnicodeError):
                                continue
                    if result:
                        break
                if result:
                    break
        # 2. lib.zip
        if not result:
            zf = cls._get_zip()
            if zf is not None:
                base = norm_path.split('/')[-1].lower()
                for zname in zf.namelist():
                    if (zname.lower().endswith(f'/{base}.asy')
                            and zname.startswith('lib/sym/')):
                        try:
                            raw = zf.read(zname)
                        except KeyError:
                            continue
                        for enc in ('utf-8', 'utf-16-le', 'latin-1'):
                            try:
                                result = _extract(raw.decode(enc))
                                break
                            except (UnicodeDecodeError, UnicodeError):
                                continue
                        if result:
                            break

        cls._attr_cache[cache_key] = result
        return result

    @staticmethod
    def _transform_point(x: int, y: int, rotation: str) -> Tuple[int, int]:
        """Transform a point by LTspice rotation/mirror.

        LTspice rotations: R0, R90, R180, R270, M0, M90, M180, M270
        In LTspice coordinate system (Y increases downward):
        - R90 = 90 degrees clockwise
        - M = mirror about Y axis (negate X), then apply rotation
        """
        mirror = rotation.startswith('M')
        angle = int(rotation[1:]) if len(rotation) > 1 else 0

        if mirror:
            x = -x

        if angle == 0:
            return (x, y)
        elif angle == 90:
            return (-y, x)
        elif angle == 180:
            return (-x, -y)
        elif angle == 270:
            return (y, -x)
        return (x, y)


# =============================================================================
# ASCパーサー
# =============================================================================

class AscParser:
    """LTSpice .ascファイルパーサー"""

    def __init__(self, asy_search_dirs: Optional[List[Path]] = None):
        self.wires: List[AscWire] = []
        self.flags: List[AscFlag] = []
        self.symbols: List[AscSymbol] = []
        self.texts: List[AscText] = []
        self.sheet_width: int = 880
        self.sheet_height: int = 680
        self.version: str = '4'
        self.asy_search_dirs: List[Path] = list(asy_search_dirs or [])
        self.source_dir: Optional[Path] = None  # directory of the .asc file

    def parse_file(self, filepath: str) -> 'AscParser':
        """ファイルからASCをパース"""
        path = Path(filepath)
        self.source_dir = path.parent
        # Add source directory to .asy search path
        if self.source_dir not in self.asy_search_dirs:
            self.asy_search_dirs.insert(0, self.source_dir)
        # LTSpice 17+ はUTF-16 LE、それ以前はUTF-8/ASCII
        # Robust decode (UTF-16 LE/BE w/ or w/o BOM, UTF-8, cp1252) via the
        # shared decoder -- single source of truth with the CLI reader. The
        # old per-encoding loop dropped BOM-less UTF-16 .asc (UniversalOpAmp*)
        # into a NUL-riddled string where 'Version'/'SHEET' never matched.
        from .._shared import decode_asc_bytes
        for encoding in ('robust',):
            try:
                text = decode_asc_bytes(path.read_bytes())
                if text.startswith('\ufeff'):
                    text = text[1:]  # BOM除去
                if 'Version' in text or 'SHEET' in text:
                    break
            except (UnicodeDecodeError, UnicodeError):
                continue
        else:
            raise ValueError(f"Cannot decode {filepath}")

        return self.parse_string(text)

    def parse_string(self, text: str) -> 'AscParser':
        """文字列からASCをパース"""
        self.wires = []
        self.flags = []
        self.symbols = []
        self.texts = []

        lines = text.split('\n')
        current_symbol: Optional[AscSymbol] = None

        for line in lines:
            line = line.rstrip()
            if not line:
                continue

            if line.startswith('Version '):
                self.version = line.split(' ', 1)[1].strip()

            elif line.startswith('SHEET '):
                parts = line.split()
                if len(parts) >= 4:
                    self.sheet_width = int(parts[2])
                    self.sheet_height = int(parts[3])

            elif line.startswith('WIRE '):
                parts = line.split()
                if len(parts) >= 5:
                    self.wires.append(AscWire(
                        x1=int(parts[1]), y1=int(parts[2]),
                        x2=int(parts[3]), y2=int(parts[4])
                    ))

            elif line.startswith('FLAG '):
                parts = line.split()
                if len(parts) >= 4:
                    self.flags.append(AscFlag(
                        x=int(parts[1]), y=int(parts[2]),
                        name=' '.join(parts[3:])
                    ))

            elif line.startswith('SYMBOL '):
                # 前のシンボルを確定
                if current_symbol:
                    self.symbols.append(current_symbol)

                parts = line.split()
                sym_type = parts[1] if len(parts) > 1 else ''
                # パス区切りを統一（backslash → 最後の部分のみ）
                sym_base = sym_type.split('\\')[-1].lower()

                current_symbol = AscSymbol(
                    symbol_type=sym_base,
                    x=int(parts[2]) if len(parts) > 2 else 0,
                    y=int(parts[3]) if len(parts) > 3 else 0,
                    rotation=parts[4] if len(parts) > 4 else 'R0',
                    symbol_path=sym_type,  # keep full path for .asy lookup
                )

            elif line.startswith('WINDOW ') and current_symbol:
                current_symbol.windows.append(line)

            elif line.startswith('SYMATTR ') and current_symbol:
                parts = line.split(' ', 2)
                attr_name = parts[1] if len(parts) > 1 else ''
                attr_value = parts[2] if len(parts) > 2 else ''

                if attr_name == 'InstName':
                    current_symbol.inst_name = attr_value
                elif attr_name == 'Value':
                    current_symbol.value = attr_value
                elif attr_name == 'Value2':
                    current_symbol.value2 = attr_value
                elif attr_name == 'SpiceModel':
                    current_symbol.spice_model = attr_value
                elif attr_name == 'SpiceLine':
                    current_symbol.spice_line = attr_value

            elif line.startswith('TEXT '):
                parts = line.split(' ', 5)
                if len(parts) >= 6:
                    raw_text = parts[5] if len(parts) > 5 else ''
                    is_dir = raw_text.startswith('!')
                    if is_dir:
                        raw_text = raw_text[1:]
                    self.texts.append(AscText(
                        x=int(parts[1]), y=int(parts[2]),
                        text=raw_text,
                        is_directive=is_dir,
                    ))
                elif len(parts) >= 5:
                    raw_text = parts[4] if len(parts) > 4 else ''
                    is_dir = raw_text.startswith('!')
                    if is_dir:
                        raw_text = raw_text[1:]
                    self.texts.append(AscText(
                        x=int(parts[1]), y=int(parts[2]),
                        text=raw_text,
                        is_directive=is_dir,
                    ))

        # 最後のシンボルを確定
        if current_symbol:
            self.symbols.append(current_symbol)

        return self

    def get_component_terminals(self, sym: AscSymbol) -> Optional[Tuple[Tuple[int,int], ...]]:
        """シンボルの端子座標を計算（4端子、2端子、3端子、5端子、.asy、または推定）"""
        # 4端子素子を先にチェック（E/G/SW/Tは2Tテーブルにも出力端子のみ登録されているため）
        offsets4 = TERMINAL_OFFSETS_4.get(sym.symbol_type, {}).get(sym.rotation)
        if offsets4 is not None:
            return tuple(
                (sym.x + off[0], sym.y + off[1]) for off in offsets4
            )

        # 2端子素子（ハードコードテーブル）
        offsets = TERMINAL_OFFSETS.get(sym.symbol_type, {}).get(sym.rotation)
        if offsets is not None:
            off1, off2 = offsets
            t1 = (sym.x + off1[0], sym.y + off1[1])
            t2 = (sym.x + off2[0], sym.y + off2[1])
            return (t1, t2)

        # 3端子素子（ハードコードテーブル）
        offsets3 = TERMINAL_OFFSETS_3.get(sym.symbol_type, {}).get(sym.rotation)
        if offsets3 is not None:
            return tuple(
                (sym.x + off[0], sym.y + off[1]) for off in offsets3
            )

        # 5端子素子（opamp: ハードコードテーブル）
        offsets5 = TERMINAL_OFFSETS_5.get(sym.symbol_type, {}).get(sym.rotation)
        if offsets5 is not None:
            return tuple(
                (sym.x + off[0], sym.y + off[1]) for off in offsets5
            )

        # .asy ファイルからピン定義を取得 (Improvement A)
        asy_offsets = AsyParser.get_terminal_offsets(
            sym.symbol_path or sym.symbol_type,
            sym.rotation,
            search_dirs=self.asy_search_dirs if self.asy_search_dirs else None
        )
        if asy_offsets is not None:
            return tuple(
                (sym.x + off[0], sym.y + off[1]) for off in asy_offsets
            )

        # 未知のシンボル → ワイヤ端点から端子位置を推定
        return self._estimate_terminals(sym)

    def _estimate_terminals(self, sym: AscSymbol,
                             search_radius: int = 160
                             ) -> Optional[Tuple[Tuple[int,int], ...]]:
        """ワイヤ端点からシンボルの端子位置を推定"""
        nearby = []
        for w in self.wires:
            for wx, wy in [(w.x1, w.y1), (w.x2, w.y2)]:
                dist = abs(wx - sym.x) + abs(wy - sym.y)
                if dist <= search_radius:
                    nearby.append((wx, wy, dist))

        if not nearby:
            return None

        # 距離でソートし、重複除去
        seen = set()
        unique = []
        for wx, wy, d in sorted(nearby, key=lambda x: x[2]):
            if (wx, wy) not in seen:
                seen.add((wx, wy))
                unique.append((wx, wy))

        if len(unique) >= 2:
            return tuple(unique[:max(len(unique), 8)])  # 最大8ピン
        elif len(unique) == 1:
            return (unique[0],)

        return None

    def get_passive_symbols(self) -> List[AscSymbol]:
        """受動素子（R,C,L,V,I）のシンボルのみ返す"""
        return [s for s in self.symbols if s.symbol_type in PASSIVE_SYMBOLS]

    def has_only_passives(self) -> bool:
        """受動素子のみで構成されているか"""
        for s in self.symbols:
            if s.symbol_type not in PASSIVE_SYMBOLS:
                return False
        return True

    def has_only_supported(self) -> bool:
        """サポートされた素子のみで構成されているか"""
        for s in self.symbols:
            if s.symbol_type not in SUPPORTED_SYMBOLS:
                return False
        return True

    def get_symbol_types(self) -> Set[str]:
        """全シンボルタイプを返す"""
        return {s.symbol_type for s in self.symbols}


# =============================================================================
# ネットリスト抽出器
# =============================================================================

class NetlistExtractor:
    """ASCからSPICEネットリストを抽出する

    ワイヤ接続を辿って各コンポーネントの端子がどのノードに
    繋がっているかを特定し、SPICEネットリストを生成する。
    """

    def __init__(self, asc: AscParser):
        self.asc = asc
        # 座標 → ノード名のマッピング
        self.coord_to_node: Dict[Tuple[int,int], str] = {}
        self._node_counter = 0
        # InstName → SPICE name rename map (for fixing K-statements etc.)
        self._name_remap: Dict[str, str] = {}

    def extract(self) -> str:
        """SPICEネットリストを生成"""
        # Step 1: ワイヤで接続された座標をグループ化（同一ネット）
        net_groups = self._build_net_groups()

        # Step 2: FLAGからノード名を割り当て
        self._assign_flag_names(net_groups)

        # Step 3: 全座標にノード名を割り当て
        self._assign_node_names(net_groups)

        # Step 4: コンポーネントのネットリスト行を生成
        lines = []
        lines.append('* Extracted from ASC')

        # Standard LTspice symbol_type per SPICE class. If the actual symbol
        # is something else (ind2, schottky, pnp, polcap, npn3, ...), emit a
        # hint comment so NetlistToAsc can restore the original SYMBOL kind.
        _DEFAULT_SYM_BY_PREFIX = {
            'R': {'res'}, 'C': {'cap'}, 'L': {'ind'},
            'V': {'voltage'}, 'I': {'current'},
            'D': {'diode'},
            'Q': {'npn'}, 'M': {'nmos'}, 'J': {'njf'},
            'E': {'e', 'e2'}, 'G': {'g', 'g2'}, 'F': {'f'}, 'H': {'h'},
            'B': {'bv', 'bi', 'bi2'}, 'S': {'sw'}, 'T': {'tline'},
        }
        for sym in self.asc.symbols:
            spice_line = self._symbol_to_spice(sym)
            if not spice_line:
                continue
            # Emit hint BEFORE the component line. NetlistParser's
            # ``pending_symbol_hint`` mechanism associates a ``* @sym=...``
            # comment with the NEXT component, so the comment must come
            # first.  Before v0.3.10 the hint was appended after the
            # component, which silently shifted every hint by one and
            # could orphan onto a SPICE class that doesn't match (e.g.
            # a polcap hint orphaning onto an Rload resistor → re-emitted
            # as a polcap, then dropped by the round-trip).
            first = spice_line.lstrip()[:1].upper()
            sym_type = (sym.symbol_type or '').lower()
            defaults = _DEFAULT_SYM_BY_PREFIX.get(first, set())
            # Do NOT emit a hint before an A-element line: NetlistParser
            # absorbs `A...` as a directive, so a preceding `* @sym=` would
            # have no component to bind to and would leak onto the *next*
            # real component (turning e.g. V1 into a bogus X§V1 ... sample).
            if (sym_type and sym_type not in defaults
                    and first not in ('X', 'A')):
                # Avoid leaking proprietary names: only emit the bare symbol
                # type string, no path or model.
                lines.append(f'* @sym={sym.symbol_type}')
            lines.append(spice_line)

        # Step 5: ディレクティブ (Improvement B & E)
        # TEXT directives in .asc use literal \n as line separator.
        # We must restore them as real newlines. Multi-line .subckt blocks,
        # .model, .param, .lib, .include, .ic, continuation lines (+) are
        # all preserved.
        for text in self.asc.texts:
            if text.is_directive:
                # Split on literal \n (which appears as \\n in the parsed text)
                sub_lines = text.text.split('\\n')
                for sub in sub_lines:
                    sub = sub.strip()
                    if sub:
                        # Apply name remap to K-statements and other directives
                        # that reference component names
                        sub = self._apply_name_remap(sub)
                        lines.append(sub)
            elif text.text.strip().startswith('.'):
                # Non-directive TEXT that looks like a SPICE directive
                # (sometimes users forget the ! prefix)
                sub_lines = text.text.split('\\n')
                for sub in sub_lines:
                    sub = sub.strip()
                    if sub and sub.startswith('.'):
                        sub = self._apply_name_remap(sub)
                        lines.append(sub)

        lines.append('.end')
        return '\n'.join(lines)

    def _apply_name_remap(self, line: str) -> str:
        """Apply InstName→SPICE name remap to directive lines (e.g. K-statements).

        The regex uses a negative lookbehind for `§` so an already-prefixed
        token (e.g. `X§Q1` inside a `.subckt` body line emitted earlier by
        this same extractor) is not re-prefixed into `X§X§Q1` on a second
        pass through the pipeline.
        """
        if not self._name_remap:
            return line
        import re
        # K-statement: K<name> L1 L2 [L3...] <coupling>
        # Replace component names that were remapped
        for old_name, new_name in self._name_remap.items():
            # `(?<!§)` blocks the match if the candidate is already
            # preceded by the rename separator. `\b` keeps the original
            # word-boundary safety.
            line = re.sub(r'(?<!§)\b' + re.escape(old_name) + r'\b',
                          new_name, line)
        return line

    def _build_net_groups(self) -> List[Set[Tuple[int,int]]]:
        """ワイヤで接続された座標をUnion-Findでグループ化

        LTSpiceのワイヤはコンポーネント端子を通り抜ける形で配線される。
        端子位置でワイヤを分割（切断）してからUnion-Findを実行する。
        """
        # Step 1: 全端子座標を集める
        terminal_coords: Set[Tuple[int,int]] = set()
        for sym in self.asc.symbols:
            terms = self.asc.get_component_terminals(sym)
            if terms:
                for t in terms:
                    terminal_coords.add(t)

        # Step 2: コンポーネントの端子ペア（内部接続）を集める
        # 隣接する端子間をペアとして登録
        internal_pairs: Set[Tuple[Tuple[int,int], Tuple[int,int]]] = set()
        for sym in self.asc.symbols:
            terms = self.asc.get_component_terminals(sym)
            if terms and len(terms) >= 2:
                for i in range(len(terms)):
                    for j in range(i + 1, len(terms)):
                        internal_pairs.add((terms[i], terms[j]))
                        internal_pairs.add((terms[j], terms[i]))

        # Step 3: ワイヤを端子位置で分割
        split_wires = []
        for w in self.asc.wires:
            segments = self._split_wire_at_terminals(w, terminal_coords)
            split_wires.extend(segments)

        # Step 4: コンポーネント内部を通るセグメントを除外
        # セグメントの両端がともに同一コンポーネントの端子ペアなら除外
        external_wires = []
        for seg in split_wires:
            p1 = (seg[0], seg[1])
            p2 = (seg[2], seg[3])
            if (p1, p2) in internal_pairs:
                continue  # コンポーネント内部のワイヤ
            external_wires.append(seg)
        split_wires = external_wires

        # (snap tolerance already applied in _build_net_groups)

        # Step 3: Union-Find
        parent: Dict[Tuple[int,int], Tuple[int,int]] = {}

        def find(c):
            if c not in parent:
                parent[c] = c
            while parent[c] != c:
                parent[c] = parent[parent[c]]
                c = parent[c]
            return c

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        # 全座標を集める
        all_coords: Set[Tuple[int,int]] = set()
        for w in split_wires:
            all_coords.add((w[0], w[1]))
            all_coords.add((w[2], w[3]))
        for f in self.asc.flags:
            all_coords.add((f.x, f.y))
        all_coords.update(terminal_coords)

        # 分割済みワイヤ端点を結合
        for w in split_wires:
            union((w[0], w[1]), (w[2], w[3]))

        # FLAG座標がワイヤ上にあるかチェック
        for f in self.asc.flags:
            fc = (f.x, f.y)
            for w in split_wires:
                wire = AscWire(w[0], w[1], w[2], w[3])
                if self._point_on_wire(f.x, f.y, wire):
                    union(fc, (w[0], w[1]))
                    break

        # Jumper / net-tie symbols are zero-ohm shorts, not components:
        # union their pins so the two nets collapse into one. This matches
        # LTspice's own netlister, which merges the nodes and drops the
        # jumper entirely (the component is suppressed in _symbol_to_spice).
        for sym in self.asc.symbols:
            if (sym.symbol_type or '').lower() in _SHORT_SYMBOLS:
                terms = self.asc.get_component_terminals(sym)
                if terms and len(terms) >= 2:
                    for t in terms[1:]:
                        union(terms[0], t)

        # Step 5: 端子座標をワイヤ端点にスナップ（近傍探索）
        # 3端子素子の端子オフセットは概算なので、余裕をもたせる
        SNAP_TOLERANCE = 64
        wire_endpoints = set()
        for seg in split_wires:
            wire_endpoints.add((seg[0], seg[1]))
            wire_endpoints.add((seg[2], seg[3]))
        for f in self.asc.flags:
            wire_endpoints.add((f.x, f.y))

        for coord in list(terminal_coords):
            if coord in wire_endpoints:
                continue
            best = None
            best_dist = SNAP_TOLERANCE + 1
            for wp in wire_endpoints:
                dist = abs(coord[0] - wp[0]) + abs(coord[1] - wp[1])
                if dist < best_dist:
                    best_dist = dist
                    best = wp
            if best is not None and best_dist <= SNAP_TOLERANCE:
                union(coord, best)

        # グループ化
        groups: Dict[Tuple[int,int], Set[Tuple[int,int]]] = {}
        for c in all_coords:
            root = find(c)
            groups.setdefault(root, set()).add(c)

        return list(groups.values())

    @staticmethod
    def _split_wire_at_terminals(w: AscWire,
                                  terminals: Set[Tuple[int,int]]
                                  ) -> List[Tuple[int,int,int,int]]:
        """ワイヤを端子位置で分割する

        ワイヤ上に端子がある場合、その位置で切断して
        複数のセグメントに分割する。
        """
        # ワイヤ上にある端子を探す
        on_wire = []
        for tx, ty in terminals:
            # 水平ワイヤ
            if w.y1 == w.y2 == ty:
                min_x = min(w.x1, w.x2)
                max_x = max(w.x1, w.x2)
                if min_x < tx < max_x:  # 端点は除外（厳密に中間のみ）
                    on_wire.append((tx, ty))
            # 垂直ワイヤ
            elif w.x1 == w.x2 == tx:
                min_y = min(w.y1, w.y2)
                max_y = max(w.y1, w.y2)
                if min_y < ty < max_y:
                    on_wire.append((tx, ty))

        if not on_wire:
            return [(w.x1, w.y1, w.x2, w.y2)]

        # 端点と中間点をソートしてセグメントに分割
        points = [(w.x1, w.y1)] + on_wire + [(w.x2, w.y2)]

        if w.y1 == w.y2:
            # 水平: X座標でソート
            points.sort(key=lambda p: p[0])
        else:
            # 垂直: Y座標でソート
            points.sort(key=lambda p: p[1])

        segments = []
        for i in range(len(points) - 1):
            p1 = points[i]
            p2 = points[i + 1]
            if p1 != p2:
                segments.append((p1[0], p1[1], p2[0], p2[1]))

        return segments

    @staticmethod
    def _point_on_wire(px: int, py: int, w: AscWire) -> bool:
        """点(px,py)がワイヤw上にあるかチェック（端点含む）"""
        # 水平ワイヤ
        if w.y1 == w.y2 == py:
            min_x = min(w.x1, w.x2)
            max_x = max(w.x1, w.x2)
            return min_x <= px <= max_x
        # 垂直ワイヤ
        if w.x1 == w.x2 == px:
            min_y = min(w.y1, w.y2)
            max_y = max(w.y1, w.y2)
            return min_y <= py <= max_y
        return False

    def _assign_flag_names(self, net_groups: List[Set[Tuple[int,int]]]):
        """FLAGの位置からノード名を割り当て"""
        for flag in self.asc.flags:
            coord = (flag.x, flag.y)
            # このFLAGが所属するグループを探す
            for group in net_groups:
                if coord in group:
                    name = flag.name
                    for c in group:
                        self.coord_to_node[c] = name
                    break

    def _assign_node_names(self, net_groups: List[Set[Tuple[int,int]]]):
        """名前未割り当てのグループに自動名を付ける"""
        for group in net_groups:
            # 既に名前があるか
            existing_name = None
            for c in group:
                if c in self.coord_to_node:
                    existing_name = self.coord_to_node[c]
                    break

            if existing_name is None:
                self._node_counter += 1
                # Match LTspice's canonical anonymous-node form: N001, N002,
                # ..., N999 (uppercase N, zero-padded width-3). LTspice uses
                # this whenever a wire net has no FLAG/label. Earlier we used
                # `n1`/`n2` lowercase which broke .raw trace-name matching
                # against LTspice's own output (V(N010) vs V(n10) miss).
                existing_name = f'N{self._node_counter:03d}'

            for c in group:
                if c not in self.coord_to_node:
                    self.coord_to_node[c] = existing_name

    def _get_node_at(self, coord: Tuple[int,int]) -> str:
        """座標のノード名を取得"""
        return self.coord_to_node.get(coord, f'?{coord}')

    def _symbol_to_spice(self, sym: AscSymbol) -> Optional[str]:
        """シンボルをSPICEコンポーネント行に変換"""
        # Jumper / net-tie symbols are shorts, not components — their pins
        # were already merged in _build_net_groups. Emit nothing (matches
        # LTspice, which drops the jumper from the netlist).
        if (sym.symbol_type or '').lower() in _SHORT_SYMBOLS:
            return None

        terms = self.asc.get_component_terminals(sym)
        if terms is None or len(terms) == 0:
            return None

        prefix = SYMBOL_TO_SPICE.get(sym.symbol_type)
        # For symbols that are unmapped, or only generically mapped to X
        # (subcircuit), consult the symbol's own .asy `SYMATTR Prefix` --
        # the class LTspice itself uses. e.g. `xtal` declares `Prefix C`,
        # so a crystal netlists as a capacitor (with Rser/Lser/Cpar), not a
        # generic subcircuit. The .asy is authoritative; a usable letter
        # overrides a hardcoded X.
        if prefix is None or prefix == 'X':
            asy_prefix = AsyParser.get_symbol_attr(
                sym.symbol_path or sym.symbol_type or '', 'Prefix',
                self.asc.asy_search_dirs or None,
            ).upper()
            if asy_prefix in _SPICE_CLASS_LETTERS:
                prefix = asy_prefix
            elif prefix is None:
                prefix = 'X'
            # else: keep the hardcoded 'X'

        name = sym.inst_name or f'{prefix}?'
        # Fix prefix mismatch: InstName's first letter disagrees with symbol's
        # SPICE prefix. e.g., ind2 symbol with InstName "T3a" → "L§T3a",
        # NJF symbol with InstName "Q1" → "J§Q1", nmos symbol with InstName
        # "Q1" → "M§Q1".
        # Use "§" (U+00A7) as the separator to match LTspice's own canonical
        # output (e.g. `M§Q1` appears verbatim in LTspice-generated .net files).
        # Earlier revisions used "_" which broke the round-trip with LTspice's
        # native simulator: probes/measurements that referenced the original
        # name lost track. The § form is what `LTspice.exe -netlist X.asc`
        # writes when InstName/symbol-prefix disagree.
        # Prefix-fix when InstName's first letter does not match the symbol's
        # SPICE prefix. Examples:
        #   ind2 symbol  + InstName "T3a"  → "L§T3a"  (T is another SPICE prefix)
        #   res  symbol  + InstName "NTC"  → "R§NTC"  (N is not a SPICE prefix
        #                                              for res, but LTspice still
        #                                              prefix-fixes user-given
        #                                              names that disagree)
        # The earlier "name[0] in _SPICE_PREFIXES" gate missed the second case
        # and silently dropped the component on the netlist->asc round-trip.
        if (prefix and prefix != 'X'
                and name and len(name) >= 2
                and name[0].isalpha()
                and name[0].upper() != prefix.upper()):
            new_name = f'{prefix}§{name}'
            self._name_remap[name] = new_name
            name = new_name
        # Unknown vendor symbol → X-prefix (subcircuit invocation). Preserves
        # the symbol type as the model name so the round-trip can re-emit it.
        # Exception: an InstName starting with 'A' is an LTspice special-
        # function device (sample&hold, varistor, phase detector, ...). It
        # must keep its 'A' prefix so the A-element path below emits a real
        # `A<n> ...` line — the same form LTspice's own netlister produces.
        # Prefixing it to `X§A1` mis-emits it as a generic subcircuit and
        # silently inflates the component count.
        elif (prefix == 'X'
                and name and name[0].isalpha()
                and name[0].upper() not in ('X', 'A')):
            new_name = f'X§{name}'
            self._name_remap[name] = new_name
            name = new_name
        value = sym.value or ''
        model = sym.spice_model or value

        # A-element (digital/special function): 8-slot pin format
        if name and name[0].upper() == 'A':
            asy_pins = AsyParser.find_and_parse(
                sym.symbol_path or sym.symbol_type,
                search_dirs=self.asc.asy_search_dirs if self.asc.asy_search_dirs else None
            )
            if asy_pins and terms:
                # Build 8-slot array, fill connected pins by SpiceOrder
                slots = ['0'] * 8
                for pin, term in zip(asy_pins, terms):
                    idx = pin.spice_order - 1  # SpiceOrder is 1-based
                    if 0 <= idx < 8:
                        slots[idx] = self._get_node_at(term)
                model_name = (sym.spice_model
                              or AsyParser.get_spice_model(
                                  sym.symbol_path or sym.symbol_type,
                                  self.asc.asy_search_dirs)
                              or sym.symbol_type.upper())
                params = f' {value}' if value else ''
                return f'{name} {" ".join(slots)} {model_name}{params}'.strip()

        if len(terms) > 3:
            nodes = [self._get_node_at(t) for t in terms]

            # E/G/SW/T: 4端子制御ソース → out+ out- ctrl+ ctrl- value
            if sym.symbol_type in ('e', 'e2', 'g', 'g2', 'sw', 'tline'):
                return f'{name} {" ".join(nodes)} {value}'.strip()

            # 5端子opamp → U name in+ in- V+ V- out model
            if sym.symbol_type in SYMBOL_TO_SPICE and SYMBOL_TO_SPICE[sym.symbol_type] == 'U':
                model = sym.spice_model or sym.symbol_type
                return f'{name} {" ".join(nodes)} {model}'.strip()

            # 多端子素子（サブサーキット/IC）
            model = sym.spice_model or sym.symbol_type
            return f'{name} {" ".join(nodes)} {model}'.strip()

        elif len(terms) == 3:
            # 3端子素子 (BJT, MOSFET, JFET)
            node1 = self._get_node_at(terms[0])  # C/D
            node2 = self._get_node_at(terms[1])  # B/G
            node3 = self._get_node_at(terms[2])  # E/S

            if sym.symbol_type in ('npn', 'pnp'):
                # Q C B E model
                return f'{name} {node1} {node2} {node3} {model}'.strip()
            elif sym.symbol_type in ('nmos', 'pmos'):
                # M D G S B model (B=S for simplicity)
                return f'{name} {node1} {node2} {node3} {node3} {model}'.strip()
            elif sym.symbol_type in ('njf', 'pjf'):
                # J D G S model
                return f'{name} {node1} {node2} {node3} {model}'.strip()
            elif sym.symbol_type in ('nigbt',):
                return f'{name} {node1} {node2} {node3} {model}'.strip()
            else:
                return f'{name} {node1} {node2} {node3} {value}'.strip()

        elif len(terms) == 1:
            # 1端子（電源ピン、constant source 等）
            node1 = self._get_node_at(terms[0])
            # X-prefixed → 1-pin subcircuit invocation (e.g. PowerSim CONST,
            # 1-pin constant source). Emit the symbol kind as the subckt name
            # so the round-trip preserves the component; without this, the
            # extractor produces a 2-token line like ``X10 N031`` which
            # NetlistParser drops as too-short.
            if name and name[0].upper() == 'X':
                sub_model = sym.spice_model or sym.symbol_type
                return f'{name} {node1} {sub_model}'.strip()
            return f'{name} {node1} {value}'.strip()

        elif len(terms) == 2:
            node1 = self._get_node_at(terms[0])
            node2 = self._get_node_at(terms[1])

            if sym.symbol_type in ('voltage', 'current', 'battery'):
                line = f'{name} {node1} {node2} {value}'
                if sym.value2:
                    line += f' {sym.value2}'
                return line.strip()
            elif sym.symbol_type in ('res', 'cap', 'ind', 'ind2', 'polcap'):
                line = f'{name} {node1} {node2} {value}'
                if sym.value2:
                    line += f' {sym.value2}'
                if sym.spice_line:
                    line += f' {sym.spice_line}'
                return line.strip()
            elif sym.symbol_type in ('diode', 'schottky', 'zener', 'led',
                                      'varactor', 'tvs'):
                return f'{name} {node1} {node2} {model}'.strip()
            elif sym.symbol_type in ('e', 'e2', 'g', 'g2'):
                # 制御ソース: 出力端子のみ。制御端子は別途。
                return f'{name} {node1} {node2} {value}'.strip()
            elif sym.symbol_type in ('bv', 'bi', 'bi2'):
                return f'{name} {node1} {node2} {value}'.strip()
            elif sym.symbol_type in ('sw',):
                return f'{name} {node1} {node2} {model}'.strip()
            else:
                # Unknown 2-terminal symbol: preserve symbol_type as model so
                # the round-trip back to .asc can emit the right SYMBOL line.
                # X-prefixed name already signals a subcircuit invocation
                # (see name-remap above).
                if name and name[0].upper() == 'X':
                    sub_model = sym.spice_model or sym.symbol_type
                    return f'{name} {node1} {node2} {sub_model}'.strip()
                return f'{name} {node1} {node2} {value}'.strip()

        return None


# =============================================================================
# 便利関数
# =============================================================================

def asc_to_netlist(filepath: str,
                   asy_search_dirs: Optional[List[Path]] = None,
                   prefer_ltspice: Optional[bool] = None,
                   ltspice_exe: Optional[str] = None,
                   ltspice_timeout: float = 20.0) -> str:
    """ASCファイルからネットリストを抽出。

    Two backends:
      1. LTspice -netlist (canonical, matches LTspice's own anonymous-node
         numbering exactly). Used when `prefer_ltspice=True` (default when
         LTspice.exe is available).
      2. Pure-Python NetlistExtractor (fallback when LTspice not available
         or `prefer_ltspice=False`).

    The LTspice backend is canonical because its anonymous-node numbering
    (N001, N002, …) is generated by the same algorithm LTspice itself
    uses when reading the .asc, so the resulting .cir round-trips through
    LTspice's simulator without trace-name drift in the .raw output.

    Set env `LTSPICE_NETLIST_PREFER=0` (or pass `prefer_ltspice=False`)
    to force the pure-Python path — useful for environments without
    LTspice or when reproducibility across machines is required.
    """
    if prefer_ltspice is None:
        env = os.environ.get("LTSPICE_NETLIST_PREFER", "1")
        prefer_ltspice = env not in ("0", "false", "False", "")

    if prefer_ltspice:
        exe = ltspice_exe or _find_ltspice_exe()
        if exe:
            try:
                return _ltspice_emit_netlist(filepath, exe, timeout=ltspice_timeout)
            except Exception:
                # Fall through to pure-Python on any LTspice failure
                # (timeout, license issue, malformed .asc, etc.).
                pass

    parser = AscParser(asy_search_dirs=asy_search_dirs)
    parser.parse_file(filepath)
    extractor = NetlistExtractor(parser)
    return extractor.extract()


def asc_to_cir(asc_path: str, cir_path: str = None,
               asy_search_dirs: Optional[List[Path]] = None) -> str:
    """ASCファイルを.cirファイルに変換"""
    netlist = asc_to_netlist(asc_path, asy_search_dirs=asy_search_dirs)
    if cir_path is None:
        cir_path = str(Path(asc_path).with_suffix('.cir'))
    with open(cir_path, 'w', encoding='utf-8') as f:
        f.write(netlist)
    return netlist


def classify_asc(filepath: str) -> dict:
    """ASCファイルを解析して分類情報を返す"""
    parser = AscParser()
    try:
        parser.parse_file(filepath)
    except Exception as e:
        return {'error': str(e), 'parseable': False}

    sym_types = parser.get_symbol_types()
    passive_only = parser.has_only_passives()

    supported_only = all(st in SUPPORTED_SYMBOLS for st in sym_types)

    return {
        'parseable': True,
        'num_symbols': len(parser.symbols),
        'num_wires': len(parser.wires),
        'num_flags': len(parser.flags),
        'num_texts': len(parser.texts),
        'symbol_types': sorted(sym_types),
        'passive_only': passive_only,
        'supported_only': supported_only,
        'has_directives': any(t.is_directive for t in parser.texts),
    }


# =============================================================================
# テスト
# =============================================================================

if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        filepath = sys.argv[1]
        print(f"Parsing: {filepath}")
        print(asc_to_netlist(filepath))
    else:
        print("usage: python -m ltspice_converter.parser.asc_parser <file.asc>",
              file=sys.stderr)
        sys.exit(2)
