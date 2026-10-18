"""Protect semantic text, encoded binaries and every wheel member."""
import importlib.util
from pathlib import Path
import zipfile

spec = importlib.util.spec_from_file_location("public_names", Path(__file__).with_name("check_public_names.py"))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_retired_product_semantic_text_and_encoded_binary():
    names = ["Equation" + " Editor" + " 3.0", "数式" + "エディタ" + "3.0", "eqn" + "edt32.exe"]
    for name in names:
        assert module.FORBIDDEN.search(name)
        for encoding in ("utf-8", "utf-16-le", "utf-16-be"):
            assert module.binary_hits(b"x" + name.encode(encoding) + b"\0\0")
    assert not module.binary_hits(b"Eqnedit64 version 3.1.5")


def test_wheel_metadata_native_and_member_names_are_scanned(tmp_path):
    forbidden = "Equation" + " Editor" + " 3.0"
    for member in ("package/tool.pyd", "package.dist-info/METADATA", "eqn" + "edt32.exe"):
        wheel = tmp_path / "candidate.whl"
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr(member, forbidden.encode("utf-16-le"))
        assert module.artifact_hits(wheel)
