"""Published PE CodeView records retain PDB filenames without private paths."""
from pathlib import Path, PureWindowsPath

import pefile
import pytest


PACKAGE = Path(__file__).resolve().parents[1] / "src" / "cubit_mesh_export"


@pytest.mark.parametrize("payload,pdb", [
    ("cubit_mesh_export.ccm", "cubit_mesh_export_ccm.pdb"),
    ("cubit_mesh_curver.pyd", "cubit_mesh_curver.pdb"),
])
def test_native_codeview_uses_only_the_pdb_filename(payload, pdb):
    path = PACKAGE / payload
    assert path.is_file(), f"Required rebuilt payload missing: {path}"
    with pefile.PE(str(path)) as pe:
        records = [pe.get_data(entry.struct.AddressOfRawData, entry.struct.SizeOfData)
                   for entry in pe.DIRECTORY_ENTRY_DEBUG if entry.struct.Type == 2]
    assert len(records) == 1
    record = records[0]
    assert record[:4] == b"RSDS"
    filename = record[24:].split(b"\0", 1)[0].decode("utf-8")
    assert filename == PureWindowsPath(filename).name, filename
    assert "/" not in filename and "\\" not in filename, filename
    assert filename == pdb
