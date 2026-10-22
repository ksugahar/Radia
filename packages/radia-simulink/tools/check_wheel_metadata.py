"""Check wheel metadata without building or certifying native GPU payloads."""
import argparse
from email.parser import BytesParser
from pathlib import Path
import shutil
import subprocess
import sys
import tomllib

from packaging.requirements import Requirement


def check(out):
    package = Path(__file__).resolve().parents[1]
    out = Path(out).resolve()
    source = out / "source"
    # The build backend may write egg-info. Keep every generated file outside
    # the checkout, and reject reuse of an old staging directory.
    shutil.copytree(package, source, ignore=shutil.ignore_patterns(
        "__pycache__", "*.egg-info", "build", "dist"))
    metadata = out / "metadata"
    metadata.mkdir()
    subprocess.run([sys.executable, "-c",
        "import setuptools.build_meta as b; b.prepare_metadata_for_build_wheel(" +
        repr(str(metadata)) + ")"], cwd=source, check=True)
    files = list(metadata.glob("*.dist-info/METADATA"))
    if len(files) != 1:
        raise ValueError("Expected one wheel metadata record")
    actual = BytesParser().parsebytes(files[0].read_bytes())
    project = tomllib.loads((package / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    if actual["Name"] != project["name"] or actual["Version"] != project["version"]:
        raise ValueError("Wheel name/version differs from package metadata")
    if actual["Requires-Python"] != project["requires-python"]:
        raise ValueError("Wheel Python requirement differs")
    expected = {str(Requirement(value)) for value in project.get('dependencies', [])}
    for extra, values in project['optional-dependencies'].items():
        for value in values:
            requirement = Requirement(value)
            marker = str(requirement.marker) if requirement.marker else None
            requirement.marker = None
            combined = str(requirement) + '; ' + (('(' + marker + ') and ') if marker else '') + 'extra == "' + extra + '"'
            expected.add(str(Requirement(combined)))
    actual_requirements = {str(Requirement(value)) for value in actual.get_all('Requires-Dist', [])}
    if actual_requirements != expected:
        raise ValueError('Wheel dependency metadata differs: ' + repr(actual_requirements ^ expected))
    return files[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    print(check(parser.parse_args().out))
