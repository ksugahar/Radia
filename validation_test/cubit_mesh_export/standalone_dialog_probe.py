"""Exercise real exporter-owned Qt dialogs, never Claro-owned actions.

Only settings storage and the OS file viewer are redirected. Mesh export and
the successful check-vol subprocess use the installed wheel and real Cubit.
"""
import json
import os
from pathlib import Path

from PySide6.QtCore import QTimer, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QDialog, QDialogButtonBox, QMessageBox, QTableWidget


def run_dialog_checks(menu, out):
    app = QApplication.instance()
    settings = out / 'export_settings.json'
    original_settings = menu._settings_path
    original_open = menu._open_in_os
    original_python = os.environ['CUBIT_MESH_EXPORT_PYTHON']
    opened = []
    errors = []
    cases = []
    menu._settings_path = lambda: str(settings)
    menu._open_in_os = opened.append

    def interact(filename, accept):
        def callback():
            dialog = app.activeModalWidget()
            try:
                assert isinstance(dialog, menu.ExportDialog), type(dialog)
                dialog._dir.setText(out.as_posix())
                dialog._filename.setText(filename)
                dialog._order.setCurrentIndex(1)
                dialog.grab().save(str(out / (filename + '.dialog.png')))
                button = QDialogButtonBox.Ok if accept else QDialogButtonBox.Cancel
                QTest.mouseClick(dialog.findChild(QDialogButtonBox).button(button), Qt.LeftButton)
            except Exception as exc:
                errors.append(repr(exc))
                if dialog is not None:
                    dialog.reject()
        QTimer.singleShot(100, callback)

    try:
        for fmt in (menu.FMT_NETGEN, menu.FMT_GMSH, menu.FMT_NASTRAN,
                    menu.FMT_VTK, menu.FMT_FEMEEM, menu.FMT_MEG):
            interact('cancelled.vol', False)
            menu.launch_export(fmt)
            assert not errors, errors
            assert not settings.exists(), 'Cancel must not persist settings'
            assert not opened
            assert not (out / 'cancelled.vol').exists()
            cases.append({'case': fmt + '-cancel', 'passed': True})

        interact('dialog-sphere.vol', True)
        menu.launch_export(menu.FMT_NETGEN)
        assert not errors, errors
        report = json.loads((out / 'dialog-sphere.vol.check.json').read_text())
        assert report['passed'] and report['n_elements'] > 0
        assert report['materials'][0]['name'] == 'body'
        assert report['boundaries'][0]['name'] == 'outer'
        assert len(opened) == 1
        assert json.loads(settings.read_text())[menu.FMT_NETGEN]['order'] == 1
        results = [w for w in app.allWidgets() if isinstance(w, QDialog)
                   and w.isVisible() and w.windowTitle() == 'Netgen Vol Export - Order 2']
        assert len(results) == 1
        tables = results[0].findChildren(QTableWidget)
        assert tables[0].item(0, 0).text() == 'body'
        assert tables[1].item(0, 0).text() == 'outer'
        results[0].grab().save(str(out / 'verified-result.png'))
        results[0].close()
        cases.append({'case': 'netgen-accept-real-export-and-check', 'passed': True})

        os.environ['CUBIT_MESH_EXPORT_PYTHON'] = str(out / 'missing-python.exe')
        warnings = []
        timer = QTimer()
        def close_warning():
            dialog = app.activeModalWidget()
            if isinstance(dialog, QMessageBox):
                warnings.append(dialog.text())
                dialog.grab().save(str(out / 'verification-error.png'))
                dialog.accept()
        timer.timeout.connect(close_warning)
        timer.start(100)
        try:
            interact('missing-checker.vol', True)
            menu.launch_export(menu.FMT_NETGEN)
        finally:
            timer.stop()
        assert not errors, errors
        assert len(warnings) == 1 and 'Cannot start standalone check-vol' in warnings[0]
        assert len(opened) == 1, 'Unverified mesh must not open in an external viewer'
        cases.append({'case': 'missing-checker-visible-error', 'passed': True})
        return cases
    finally:
        menu._settings_path = original_settings
        menu._open_in_os = original_open
        os.environ['CUBIT_MESH_EXPORT_PYTHON'] = original_python
