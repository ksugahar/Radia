#!python
"""Explicit GUI test executed only inside a task-owned Cubit process."""
import json
import os
import sys
from pathlib import Path
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QApplication, QMainWindow

out = Path(os.environ['CME_GUI_TEST_OUTPUT'])
gui = Path(os.environ['CME_GUI_TEST_SOURCE'])
sys.path.insert(0, str(gui))
import cubit
baseline = os.environ.get('CME_GUI_TEST_BASELINE') == '1'
if not baseline:
    import radia_export_menu as menu


def probe():
    result = {'gui_started': True, 'passed': False, 'baseline': baseline,
              'menu_source': None if baseline else menu.__file__}
    try:
        if baseline:
            result['passed'] = True
            return
        assert 'radia' not in sys.modules
        assert 'netgen' not in sys.modules
        assert menu.install_menu(), 'Claro menu registration failed'
        app = QApplication.instance()
        main = next(w for w in app.topLevelWidgets()
                    if isinstance(w, QMainWindow) and 'Cubit' in w.windowTitle())
        assert main.isVisible()
        # Never wrap Claro-owned menus/actions in shiboken: ownership conflicts
        # can crash Cubit during shutdown. Observe registration through Claro.
        result['registered_actions'] = len(menu._claro_keepalive) - 1
        assert result['registered_actions'] == 6
        main.grab().save(str(out / 'window.png'))
        if os.environ.get('CME_GUI_TEST_MENU_ONLY') == '1':
            result['passed'] = True
            return
        for command in ('reset', 'create sphere radius 1', 'volume all scheme tetmesh',
                        'volume all size 0.3', 'mesh volume all', 'block 1 volume all',
                        'block 1 name "body"', 'sideset 1 surface all', 'sideset 1 name "outer"'):
            cubit.cmd(command)
        cubit.cmd('export netgen "' + (out / 'sphere.vol').as_posix() + '" order 2 overwrite')
        assert (out / 'sphere.vol').is_file()
        result['passed'] = True
    except Exception as exc:
        result['error'] = repr(exc)
    finally:
        (out / 'gui-result.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
        QTimer.singleShot(0, QApplication.instance().quit)


QTimer.singleShot(4000, probe)
