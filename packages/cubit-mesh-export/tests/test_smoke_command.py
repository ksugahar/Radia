"""Regression for Cubit's explicit command-plugin directory in batch mode."""

from cubit_mesh_export.smoke_test import _batch_command


def test_batch_command_passes_existing_plugin_directory(tmp_path):
    exe = tmp_path / "coreform_cubit.exe"
    plugin_dir = tmp_path / "plugins"
    plugin_dir.mkdir()
    driver = tmp_path / "driver.jou"
    assert _batch_command(exe, driver) == [
        str(exe), "-batch", "-nographics", "-nojournal",
        "-commandplugindir", str(plugin_dir), str(driver),
    ]


def test_batch_command_omits_missing_plugin_directory(tmp_path):
    exe = tmp_path / "cubit"
    driver = tmp_path / "driver.jou"
    assert _batch_command(exe, driver) == [
        str(exe), "-batch", "-nographics", "-nojournal", str(driver),
    ]
