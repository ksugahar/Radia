from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPORT_GMSH_CPP = ROOT / "src" / "cubit_plugin" / "ExportGmshCommand.cpp"


def test_export_gmsh_command_writes_geo_and_exact_opt_sidecars():
    source = EXPORT_GMSH_CPP.read_text(encoding="utf-8")

    assert 'replace_extension(msh_filename, ".geo")' in source
    assert 'geo + ".opt"' in source
    assert 'msh_filename + ".opt"' in source
    assert 'Merge \\"' in source
    assert "Open this .geo for normal review" in source
    assert "write_gmsh_launch_companions(filename)" in source
