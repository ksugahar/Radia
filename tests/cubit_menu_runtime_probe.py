"""Standalone PySide6 probe executed by Coreform Cubit's Python 3.10.

This file intentionally has no pytest or Radia-package dependency. The outer
pytest test launches it with ``bin/python3/python.exe`` from the installed
Cubit, so the regression covers the exact Qt binding that owns the production
menu. It models the cold-start sequence that caused the bug:

1. startup hook runs before the main window exists;
2. Cubit creates a stock menu bar;
3. Cubit clears/rebuilds that bar repeatedly;
4. the startup hook is replayed once more and must remain idempotent.
"""

from __future__ import annotations

import gc
import importlib.util
import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import PySide6
from PySide6.QtCore import qVersion
from PySide6.QtWidgets import QApplication, QMainWindow


EXPECTED_ACTIONS = [
    "Netgen Vol (.vol)...",
    "GMSH...",
    "Nastran BDF...",
    "VTK...",
    "FEMEEM...",
    "MEG (ELF/MAGIC)...",
]


def _require(condition, message):
    if not condition:
        raise AssertionError(message)


def _load_export_menu(path):
    spec = importlib.util.spec_from_file_location(
        "radia_export_menu_runtime_probe", path)
    module = importlib.util.module_from_spec(spec)
    _require(spec.loader is not None, "could not create module loader")
    spec.loader.exec_module(module)
    return module


def _radia_menus(module, main):
    return [
        action.menu() for action in main.menuBar().actions()
        if action.menu() is not None
        and action.menu().objectName() == module._INSTALLED_OBJECT_NAME
    ]


def _capture_state(module, main, stage):
    menus = _radia_menus(module, main)
    _require(len(menus) == 1,
             f"{stage}: expected one Radia menu, got {len(menus)}")
    actions = [action.text() for action in menus[0].actions()]
    _require(actions == EXPECTED_ACTIONS,
             f"{stage}: action contract mismatch: {actions!r}")
    return {
        "stage": stage,
        "menu_count": len(menus),
        "action_count": len(actions),
        "actions": actions,
    }


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    _require(
        len(argv) == 1,
        "usage: cubit_menu_runtime_probe.py <radia_export_menu.py>",
    )
    module_path = Path(argv[0]).resolve()
    _require(module_path.is_file(), f"module not found: {module_path}")

    app = QApplication.instance() or QApplication([])
    module = _load_export_menu(module_path)
    main_window = QMainWindow()
    main_window.setObjectName("claro")
    available = {"main": None}
    module.find_claro = lambda: available["main"]
    sys.modules["cubit"] = SimpleNamespace()

    try:
        _require(module.install_menu() is None,
                 "pre-window startup should defer menu creation")
        observer = module._menu_persistence_filter
        _require(observer is not None, "startup did not install the observer")
        parented = observer.parent() is app
        gc.collect()
        survived_gc = module._menu_persistence_filter is observer

        available["main"] = main_window
        rebuilds = []
        stock_sequences = [
            ["&File"],
            ["&File", "&Edit", "&Help"],
            ["&File", "&View", "&Tools", "&Help"],
            ["&File", "&Edit", "&View", "&Display", "&Tools", "&Help"],
        ]
        for index, stock_names in enumerate(stock_sequences):
            main_window.menuBar().clear()
            app.processEvents()
            for name in stock_names:
                main_window.menuBar().addMenu(name)
            app.processEvents()
            rebuilds.append(
                _capture_state(module, main_window, f"rebuild-{index + 1}")
            )

        _require(module.install_menu() is not None,
                 "startup replay after GUI creation did not return a menu")
        app.processEvents()
        replay_count = len(_radia_menus(module, main_window))
        _require(replay_count == 1,
                 f"startup replay created {replay_count} Radia menus")

        payload = {
            "ok": True,
            "runtime": {
                "python": sys.version.split()[0],
                "pyside6": PySide6.__version__,
                "qt": qVersion(),
            },
            "observer": {
                "installed_before_main_window": observer is not None,
                "survived_gc": survived_gc,
                "parented_to_qapplication": parented,
            },
            "rebuilds": rebuilds,
            "replay_menu_count": replay_count,
        }
        print(json.dumps(payload, ensure_ascii=True, sort_keys=True))
        return 0
    finally:
        observer = getattr(module, "_menu_persistence_filter", None)
        if observer is not None:
            app.removeEventFilter(observer)
        main_window.close()
        app.processEvents()


if __name__ == "__main__":
    raise SystemExit(main())
