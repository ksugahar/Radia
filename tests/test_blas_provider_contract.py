"""Build-only provider admission contracts.

Guards: Architecture/provider mismatch rejection, explicit LP64 library closure,
missing-development-file rejection before configure, and keeping MKL Fortran
LAPACK declarations out of the NGSolve-facing provider boundary, while direct
Fortran callers explicitly include their provider's declaring header.
Assumes: Provider files are temporary fixtures; no native library is loaded.
"""
import importlib.util
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("radia_linux_provider", ROOT / "tools/build_linux.py")
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


@pytest.mark.parametrize("arch,provider", [("aarch64", "mkl"), ("x86_64", "nvpl"), ("riscv64", None)])
def test_wrong_provider_stops_before_library_discovery(arch, provider):
    with pytest.raises(RuntimeError, match="Provider mismatch|Unsupported Linux architecture"):
        builder.provider_configuration(arch, provider)


@pytest.mark.parametrize("arch,provider", [("x86_64", "mkl"), ("aarch64", "nvpl")])
def test_exact_provider_closure_and_missing_file(tmp_path, monkeypatch, arch, provider):
    include = tmp_path / "include"
    lib = tmp_path / "lib"
    include.mkdir()
    lib.mkdir()
    headers = ["mkl.h"] if provider == "mkl" else ["nvpl_blas_cblas.h", "nvpl_lapacke.h", "nvpl_fftw/fftw3.h"]
    libs = ["mkl_rt"] if provider == "mkl" else ["nvpl_lapack_lp64_seq", "nvpl_blas_lp64_seq", "nvpl_fftw"]
    for name in headers:
        p = include / name
        p.parent.mkdir(exist_ok=True)
        p.write_bytes(b"header fixture")
    for name in libs:
        (lib / ("lib" + name + ".so")).write_bytes(name.encode())
    monkeypatch.setenv("MKLROOT" if provider == "mkl" else "NVPL_ROOT", str(tmp_path))
    result = builder.provider_configuration(arch, provider)
    assert result["name"] == provider
    assert result["integer_abi"] == "lp64"
    assert len(result["libraries"]) == len(libs)
    assert all(len(v) == 64 for v in result["libraries"].values())
    (lib / ("lib" + libs[0] + ".so")).unlink()
    with pytest.raises(RuntimeError, match="Missing explicit"):
        builder.provider_configuration(arch, provider)

def test_shared_provider_header_selects_both_distributions():
    import json
    spec = importlib.util.spec_from_file_location("provider_ownership", ROOT / "tools/ci_distribution_contract.py")
    ownership = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ownership)
    manifest = json.loads((ROOT / "tools/ci_distribution_manifest.json").read_text(encoding="utf-8"))
    path = "src/core/rad_blas.h"
    assert set(ownership.owner_rule(manifest, path)["owners"]) == {"radia", "radia-simulink"}
    for lane in ("radia", "radia-simulink"):
        for event, projection in manifest["lanes"][lane]["events"].items():
            assert ownership.selected(path, projection), (lane, event, path)


def test_mkl_boundary_does_not_export_fortran_lapack_declarations():
    # NGSolve owns clapack.h. Loading the umbrella first causes C2116/C2733
    # for its complex Fortran signatures; only the C interfaces are needed.
    import re

    text = (ROOT / "src/core/rad_blas.h").read_text(encoding="utf-8")
    mkl_branch = text.split("#if defined(RADIA_BLAS_MKL)", 1)[1].split(
        "\n#else\n", 1)[0]
    includes = set(re.findall(r"^#include [<\"]([^>\"]+)[>\"]",
                              mkl_branch, re.MULTILINE))
    expected = {"mkl_cblas.h", "mkl_lapacke.h", "mkl_service.h"}
    assert includes == expected, (includes, expected)


def _fortran_lapack_units(root):
    """Scan source TUs and their quote includes without loading a provider."""
    import re

    units = {}
    paths = [p for p in (root / "src").rglob("*")
             if p.suffix in (".c", ".cc", ".cpp", ".cxx")]
    paths += list((root / "validation_test/numerics").glob("*.cpp"))
    for path in paths:
        text = path.read_text(encoding="utf-8")
        code = re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.S)
        names = set(re.findall(
            r"\b([sdcz](?:gesv|getrf|getrs|geqp3|orgqr|gesvd)_|RADIA_LAPACK_GETRS|RADIA_LAPACK_GESVD)\s*\(", code))
        if not names:
            continue
        relative = path.relative_to(root).as_posix()
        pattern = (r"#if defined\(RADIA_BLAS_MKL\)\s*"
                   r'#include <mkl_lapack\.h>.*?'
                   r"#elif defined\(RADIA_BLAS_NVPL\)\s*"
                   r'#include <nvpl_lapack\.h>.*?#endif')
        assert re.search(pattern, code, re.S), (
            relative, sorted(names), "missing provider Fortran declaring headers")
        assert not ({"dgetrs_", "zgetrs_", "dgesvd_"} & names), (
            relative, sorted(names), "string arguments require provider macro")
        if "RADIA_LAPACK_GETRS" in names:
            routine = "zgetrs" if "peec" in path.name else "dgetrs"
            assert "#define RADIA_LAPACK_GETRS " + routine + "_" in code, relative
            assert "#define RADIA_LAPACK_GETRS NVPL_LAPACK_" + routine in code, relative
        if "RADIA_LAPACK_GESVD" in names:
            assert "#define RADIA_LAPACK_GESVD dgesvd_" in code, relative
            assert "#define RADIA_LAPACK_GESVD NVPL_LAPACK_dgesvd" in code, relative
        seen = set()
        pending = [path]
        while pending:
            item = pending.pop()
            if item in seen:
                continue
            seen.add(item)
            content = item.read_text(encoding="utf-8")
            content = re.sub(r"/\*.*?\*/|//[^\n]*", "", content, flags=re.S)
            includes = re.findall(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]',
                                  content, re.M)
            assert not ({"la.hpp", "bla.hpp", "clapack.h"} & set(includes)), (
                relative, str(item), "Fortran declarations may reach NGSolve LAPACK")
            for include in re.findall(r'^\s*#\s*include\s*"([^"]+)"', content, re.M):
                for base in (item.parent, root / "src/core", root / "src/lib",
                             root / "src/ext/genmath"):
                    candidate = (base / include).resolve()
                    if candidate.is_relative_to(root.resolve()) and candidate.is_file():
                        pending.append(candidate)
                        break
        units[relative] = sorted(names)
    return units


def test_every_fortran_lapack_caller_has_private_provider_declarations():
    units = _fortran_lapack_units(ROOT)
    expected = {
        "src/core/rad_peec_matrices.cpp",
        "src/core/rad_relaxation_methods.cpp",
        "validation_test/numerics/test_provider_blas.cpp",
        "src/ext/HACApK/cHACApK_base.c",
    }
    assert set(units) == expected, (set(units), expected)
