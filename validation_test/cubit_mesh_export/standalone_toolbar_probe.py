"""Asynchronous official toolbar test; return to Cubit before native actions.

Cubit must be free to invoke its own Python importer. Never synchronously click
Import while an outer Python callback is waiting in a nested editor event loop.
"""
import json
import os
from pathlib import Path
import sys
from PySide6.QtCore import QTimer, Qt, QPoint, QMetaObject, QCoreApplication, QEvent
from PySide6.QtGui import QContextMenuEvent, QKeyEvent
from PySide6.QtWidgets import QApplication, QAbstractButton, QDialog, QListWidget, QLineEdit, QWidget, QTextEdit, QPlainTextEdit


def start_toolbar_checks(out, mode, finish):
    app = QApplication.instance()
    destination = Path(os.environ['CME_GUI_TEST_TOOLBAR_DESTINATION'])
    destination.mkdir(parents=True, exist_ok=True)
    completed = []
    steps = []
    dispatches = []
    pending = []
    payload = {}

    def visible(classname):
        return next(w for w in app.allWidgets()
                    if w.isVisible() and w.metaObject().className() == classname)

    def click(button):
        assert button.isVisible() and button.isEnabled(), button.text()
        assert QMetaObject.invokeMethod(button, 'click', Qt.QueuedConnection)

    def later(callback):
        def call():
            if completed:
                return
            try:
                steps.append(callback.__name__)
                (out / 'steps.json').write_text(json.dumps(steps))
                callback()
            except Exception as exc:
                console = [w.toPlainText() for w in app.allWidgets()
                           if isinstance(w, (QTextEdit, QPlainTextEdit))]
                (out / 'failure-console.json').write_text(json.dumps(console, indent=2))
                for window in app.topLevelWidgets():
                    if window.isVisible() and window.objectName() == 'Claro':
                        window.grab().save(str(out / 'failure-window.png'))
                for window in app.topLevelWidgets():
                    if isinstance(window, QDialog) and window.isVisible():
                        window.grab().save(str(out / 'failure-dialog.png'))
                        controls = [{'class': w.metaObject().className(), 'name': w.objectName(),
                                     'text': w.text() if isinstance(w, QAbstractButton) else ''}
                                    for w in window.findChildren(QWidget) if w.isVisible()]
                        (out / 'failure-controls.json').write_text(json.dumps(controls, indent=2))
                completed.append(True)
                finish({'steps': steps}, repr(exc))
        QTimer.singleShot(400, call)

    def verify():
        scripts = list(destination.rglob('scripts/cubit_export_menu.py'))
        assert len(scripts) == 1, scripts
        source = Path(os.environ['CME_GUI_TEST_SOURCE']) / 'cubit_export_menu.py'
        assert scripts[0].read_bytes() == source.read_bytes()
        buttons = [w for w in app.allWidgets() if isinstance(w, QAbstractButton)
                   and destination.as_posix().lower() in w.toolTip().replace(chr(92), '/').lower()]
        assert len(buttons) == 6, [(w.text(), w.toolTip()) for w in buttons]
        assert 'radia' not in sys.modules and 'netgen' not in sys.modules
        buttons[0].window().showNormal()
        buttons[0].window().resize(1600, 1000)
        # Give the test toolbar its own row instead of inherited user overflow.
        # The outer controller restores the exact original docking preferences.
        buttons[0].window().addToolBar(Qt.BottomToolBarArea, buttons[0].parentWidget())
        payload.update({'mode': mode, 'buttons': [w.text() for w in buttons],
                        'menu_source': str(scripts[0]), 'steps': steps})
        states = []
        for button in buttons:
            parents = []
            widget = button
            while widget is not None:
                parents.append({'class': widget.metaObject().className(),
                                'name': widget.objectName(), 'visible': widget.isVisible()})
                widget = widget.parentWidget()
            states.append({'text': button.text(), 'enabled': button.isEnabled(), 'parents': parents})
        (out / 'button-states.json').write_text(json.dumps(states, indent=2))
        if mode == 'restart':
            import cubit
            cubit.cmd('create sphere radius 1')
            pending.extend(sorted(buttons, key=lambda w: w.text()))
            later(dispatch_next)
        else:
            later(complete)

    def complete():
        import runpy
        probe = Path(os.environ['CME_GUI_TEST_SOURCE']) / 'toolbar_probe.py'
        snapshot = runpy.run_path(str(probe))['_snapshot']()
        assert snapshot['ok'], snapshot
        payload['installed_release_probe'] = snapshot
        button = next(w for w in app.allWidgets() if isinstance(w, QAbstractButton)
                      and destination.as_posix().lower() in w.toolTip().replace(chr(92), '/').lower())
        button.window().grab().save(str(out / 'persisted-toolbar.png'))
        assert 'radia' not in sys.modules and 'netgen' not in sys.modules
        completed.append(True)
        payload.update({'passed': True, 'actual_button_dispatches': dispatches})
        finish(payload, None)

    def dispatch_next():
        if not pending:
            later(complete)
            return
        button = pending[0]
        click(button)
        later(cancel_export)

    def cancel_export():
        module = sys.modules['cubit_export_menu']
        assert Path(module.__file__).resolve() == Path(payload['menu_source']).resolve()
        dialog = next(w for w in app.topLevelWidgets()
                      if isinstance(w, module.ExportDialog) and w.isVisible())
        dialog.grab().save(str(out / ('toolbar-dialog-' + str(len(dispatches)) + '.png')))
        assert QMetaObject.invokeMethod(dialog, 'reject', Qt.QueuedConnection)
        dispatches.append(pending.pop(0).text())
        later(dispatch_next)

    def close_editor():
        dialog = visible('WTEditor')
        dialog.grab().save(str(out / 'imported-editor.png'))
        button = next(w for w in dialog.findChildren(QAbstractButton)
                      if w.isVisible() and w.text().replace('&', '') == 'OK')
        click(button)
        later(verify)

    def import_summary():
        wizard = visible('WImportWizard')
        wizard.grab().save(str(out / 'import-summary.png'))
        button = next(w for w in wizard.findChildren(QAbstractButton)
                      if w.isVisible() and w.isEnabled()
                      and w.text().replace('&', '') in ('Finish', 'OK'))
        click(button)
        later(close_editor)

    def import_file():
        wizard = visible('WImportWizard')
        for name, value in (('mFile', os.environ['CME_GUI_TEST_TOOLBAR_PACKAGE']),
                            ('mDirectory', str(destination))):
            wizard.findChild(QWidget, name).findChild(QLineEdit, 'leFilename').setText(value)
        wizard.grab().save(str(out / 'import-file.png'))
        button = next(w for w in wizard.findChildren(QAbstractButton)
                      if w.isVisible() and w.text().replace('&', '') == 'Import')
        click(button)
        later(import_summary)

    def context():
        popup = app.activePopupWidget()
        assert popup is not None
        for key in (Qt.Key_Down, Qt.Key_Return):
            QCoreApplication.postEvent(popup, QKeyEvent(QEvent.KeyPress, key, Qt.NoModifier))
            QCoreApplication.postEvent(popup, QKeyEvent(QEvent.KeyRelease, key, Qt.NoModifier))
        later(import_file)

    def editor():
        dialog = visible('WTEditor')
        items = dialog.findChild(QWidget, 'mToolbars').findChild(QListWidget, 'lwItemList')
        assert items.count() == 0, 'Real toolbars must not be loaded during import test'
        point = QPoint(30, 80)
        QCoreApplication.postEvent(items.viewport(), QContextMenuEvent(
            QContextMenuEvent.Mouse, point, items.viewport().mapToGlobal(point)))
        later(context)

    if mode == 'import':
        button = next(w for w in app.allWidgets() if isinstance(w, QAbstractButton)
                      and w.isVisible() and w.text() == 'Custom Toolbar Editor')
        click(button)
        later(editor)
    else:
        later(verify)
