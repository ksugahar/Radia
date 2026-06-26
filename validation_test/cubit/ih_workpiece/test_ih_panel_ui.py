"""Test IH panel UI mode switching and command generation.

No Cubit needed. Uses offscreen Qt platform.
"""

import os
import sys

# Offscreen Qt before any import
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                '..', '..', '..', 'src', 'radia'))


def row_shown(panel, key):
    """Check if a form row is shown (not hidden)."""
    from PySide6.QtWidgets import QFormLayout
    idx = panel._row_indices.get(key)
    if idx is None:
        return False
    label = panel._form.itemAt(idx, QFormLayout.LabelRole)
    if label and label.widget():
        return not label.widget().isHidden()
    return False


def test_bem_mode():
    """BEM mode: source/sink required, workpiece optional."""
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication(sys.argv)
    from radia_ih import IHPanel

    panel = IHPanel()
    panel._method_combo.setCurrentText('BEM')

    assert row_shown(panel, 'source'), 'source should be shown'
    assert row_shown(panel, 'sink'), 'sink should be shown'
    assert row_shown(panel, 'workpiece'), 'workpiece should be shown'
    assert row_shown(panel, 'coil_sigma'), 'coil_sigma should be shown'
    assert not row_shown(panel, 'wp_coil_radius'), 'wp_coil_radius hidden in BEM'
    assert not row_shown(panel, 'wp_material'), 'wp_material hidden in BEM'
    assert not row_shown(panel, 'current'), 'current hidden in BEM'

    # Not runnable without source/sink
    assert not panel.is_runnable()

    # Set source/sink -> runnable
    panel._widgets['source'].setText('source')
    panel._widgets['sink'].setText('sink')
    assert panel.is_runnable()

    # Command references calc_inductance.py
    cmd = panel.build_command('test.vol')
    assert 'calc_inductance.py' in cmd[1]
    assert '--source' in cmd
    assert '--sink' in cmd
    print('BEM mode: OK')


def test_bem_wp_mode():
    """BEM-SIBC (WP) mode: analytical coil, workpiece surface BEM."""
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication(sys.argv)
    from radia_ih import IHPanel

    panel = IHPanel()
    panel._method_combo.setCurrentText('BEM-SIBC (WP)')

    # WP-specific fields visible
    assert row_shown(panel, 'wp_coil_radius'), 'wp_coil_radius should be shown'
    assert row_shown(panel, 'wp_coil_current'), 'wp_coil_current should be shown'
    assert row_shown(panel, 'wp_gap_deg'), 'wp_gap_deg should be shown'
    assert row_shown(panel, 'wp_bnd_label'), 'wp_bnd_label should be shown'
    assert row_shown(panel, 'wp_material'), 'wp_material should be shown'
    assert row_shown(panel, 'wp_sigma'), 'wp_sigma should be shown'
    assert row_shown(panel, 'impedance'), 'impedance should be shown'

    # BEM coil fields hidden
    assert not row_shown(panel, 'source'), 'source hidden in WP mode'
    assert not row_shown(panel, 'sink'), 'sink hidden in WP mode'
    assert not row_shown(panel, 'workpiece'), 'workpiece hidden in WP mode'
    assert not row_shown(panel, 'coil_sigma'), 'coil_sigma hidden in WP mode'

    # FEM fields hidden
    assert not row_shown(panel, 'current'), 'current hidden in WP mode'
    assert not row_shown(panel, 'solver'), 'solver hidden in WP mode'

    # Always runnable (no source/sink required)
    assert panel.is_runnable()

    # Command references calc_heating_bem.py
    cmd = panel.build_command('model.vol')
    assert 'calc_heating_bem.py' in cmd[1]
    assert '--coil-radius' in cmd
    assert '--wp-label' in cmd
    assert '--material' in cmd
    assert '--msh-output' in cmd

    # Impedance model: esim shows BH file, dowell does not
    panel._widgets['impedance'].setCurrentText('esim')
    assert row_shown(panel, 'bh_file'), 'bh_file shown for esim'
    assert row_shown(panel, 'esim_geometry'), 'esim_geometry shown for esim'

    panel._widgets['impedance'].setCurrentText('dowell')
    assert not row_shown(panel, 'bh_file'), 'bh_file hidden for dowell'
    assert not row_shown(panel, 'esim_geometry'), 'esim_geometry hidden for dowell'

    print('BEM-SIBC (WP) mode: OK')


def test_fem_mode():
    """FEM mode: Kelvin + SIBC/ESIM fields."""
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication(sys.argv)
    from radia_ih import IHPanel

    panel = IHPanel()
    panel._method_combo.setCurrentText('FEM')

    assert row_shown(panel, 'current'), 'current should be shown'
    assert row_shown(panel, 'a_coil'), 'a_coil should be shown'
    assert row_shown(panel, 'solver'), 'solver should be shown'

    assert not row_shown(panel, 'source'), 'source hidden in FEM'
    assert not row_shown(panel, 'wp_coil_radius'), 'wp_coil_radius hidden in FEM'
    assert not row_shown(panel, 'wp_material'), 'wp_material hidden in FEM'

    assert panel.is_runnable()
    print('FEM mode: OK')


def test_mode_switching():
    """Switching modes preserves frequency value."""
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication(sys.argv)
    from radia_ih import IHPanel

    panel = IHPanel()
    panel._widgets['freq'].setText('12345')

    panel._method_combo.setCurrentText('BEM-SIBC (WP)')
    assert panel.val('freq') == '12345', 'freq lost on mode switch'

    panel._method_combo.setCurrentText('FEM')
    assert panel.val('freq') == '12345', 'freq lost on mode switch'

    panel._method_combo.setCurrentText('BEM')
    assert panel.val('freq') == '12345', 'freq lost on mode switch'
    print('Mode switching: OK')


def main():
    test_bem_mode()
    test_bem_wp_mode()
    test_fem_mode()
    test_mode_switching()
    print('\nAll IH panel UI tests passed!')
    return 0


if __name__ == '__main__':
    sys.exit(main())
