"""License metadata and external-runtime distribution boundary, without native builds."""
from email.parser import Parser
from io import StringIO
from pathlib import Path
import tomllib

import pytest
from packaging.licenses import canonicalize_license_expression

ROOT = Path(__file__).resolve().parents[1]
PROJECTS = [ROOT, *sorted((ROOT / "packages").glob("*/pyproject.toml"))]


@pytest.mark.parametrize("project", PROJECTS, ids=lambda p: p.parent.name if p.is_file() else "radia")
def test_distribution_license_metadata(project, monkeypatch):
    # Use the backend's metadata writer and sdist license-file selector only;
    # running setup.py would require native solver/MATLAB release artifacts.
    from setuptools.dist import Distribution
    from setuptools.config.pyprojecttoml import apply_configuration
    from setuptools.command.egg_info import FileList, manifest_maker

    directory = project.parent if project.is_file() else project
    data = tomllib.loads((directory / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    expression = data["license"]
    assert isinstance(expression, str)
    assert canonicalize_license_expression(expression) == expression
    assert "Apache-2.0" in expression and "LicenseRef-Radia" not in expression
    patterns = data["license-files"]
    assert {"LICENSE", "NOTICE", "THIRD_PARTY_NOTICES.md"} <= set(patterns)
    monkeypatch.chdir(directory)
    dist = Distribution({"name": data["name"], "version": data["version"]})
    apply_configuration(dist, directory / "pyproject.toml")
    expected = {str(p.relative_to(directory)).replace("\\", "/") for pattern in patterns for p in directory.glob(pattern) if p.is_file()}
    assert expected
    assert set(dist.metadata.license_files) == expected
    output = StringIO()
    dist.metadata.write_pkg_file(output)
    metadata = Parser().parsestr(output.getvalue())
    assert metadata["License-Expression"] == expression
    assert set(metadata.get_all("License-File")) == expected
    command = manifest_maker(dist)
    command.filelist = FileList()
    command.add_license_files()
    assert set(command.filelist.files) == expected


def test_no_external_runtime_dlls_in_package_sources():
    roots = [ROOT / "src", ROOT / "packages"]
    forbidden = ("mkl", "libiomp", "libng", "ngsolve", "netgen", "cudart", "cublas", "cusparse", "cusolver", "cufft", "curand", "nvrtc", "nvjitlink", "cudnn", "nvcuda")
    bad = [p.relative_to(ROOT).as_posix() for root in roots for p in root.rglob("*")
           if p.is_file() and p.suffix.lower() == ".dll" and p.name.lower().startswith(forbidden)]
    assert not bad, bad
