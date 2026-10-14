"""Source-level contracts for the MATLAB/MEX runtime boundary."""

from __future__ import annotations

import tomllib
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _source(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def test_retired_cpu_gateway_cannot_be_loaded_by_matlab_setup():
    setup = _source("packages/radia-gpu/matlab/+radia/setup.m")
    boundary = _source("packages/radia-gpu/matlab/+radia/+internal/callMex.m")
    for source in (setup, boundary):
        assert "radia:mex:Retired" in source
        assert "radia_mex(" not in source
        assert "-MatlabMexOnly" not in source
    assert "MKL_THREADING_LAYER" not in setup


def test_matlab_release_uses_sequential_mkl_without_changing_python() -> None:
    build = _source("Build.ps1")
    cmake = _source("CMakeLists.txt")
    assert not (ROOT / "tools/package_simulink_release.py").exists()
    verifier = _source("tools/verify_simulink_release.py")

    assert build.index("$PythonLibrary") < build.index("$env:MKLROOT")
    assert 'MatlabMexOnly' not in build
    assert "mkl_intel_thread" in cmake
    assert "MKL_MATLAB_RUNTIME_DLLS" not in cmake
    assert "matlab/mkl_sequential.3.dll" in verifier
    assert "matlab/+radia/+internal/expectedNGSolveVersion.m" in verifier


def test_release_schema_preserves_historical_threaded_mkl_contracts() -> None:
    import importlib.util

    path = ROOT / "tools" / "verify_simulink_release.py"
    spec = importlib.util.spec_from_file_location("verify_runtime_contract", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    for historical in (
        module.LEGACY_FULL_REQUIRED_MEMBERS,
        module.PREVIEW_V2_REQUIRED_MEMBERS,
        module.FULL_REQUIRED_MEMBERS_V2,
        module.FULL_REQUIRED_MEMBERS_V3,
    ):
        assert "matlab/mkl_intel_thread.3.dll" in historical
        assert "matlab/mkl_sequential.3.dll" not in historical
        assert "matlab/+radia/+internal/expectedNGSolveVersion.m" not in historical
    for current in (module.REQUIRED_MEMBERS, module.FULL_REQUIRED_MEMBERS):
        assert "matlab/mkl_sequential.3.dll" in current
        assert "matlab/mkl_intel_thread.3.dll" not in current
        assert "matlab/+radia/+internal/expectedNGSolveVersion.m" in current


def test_verifier_accepts_each_historical_and_current_schema(tmp_path) -> None:
    import importlib.util

    path = ROOT / "tools" / "verify_simulink_release.py"
    spec = importlib.util.spec_from_file_location("verify_all_schemas", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    schemas = {
        "radia.simulink.ih-release-manifest.v1": module.LEGACY_REQUIRED_MEMBERS,
        "radia.simulink.ih-release-manifest.v2": module.PREVIEW_V2_REQUIRED_MEMBERS,
        "radia.simulink.ih-release-manifest.v3": module.REQUIRED_MEMBERS,
        "radia.simulink.library-release-manifest.v1": module.LEGACY_FULL_REQUIRED_MEMBERS,
        "radia.simulink.library-release-manifest.v2": module.FULL_REQUIRED_MEMBERS_V2,
        "radia.simulink.library-release-manifest.v3": module.FULL_REQUIRED_MEMBERS_V3,
        "radia.simulink.library-release-manifest.v4": module.FULL_REQUIRED_MEMBERS,
    }
    for schema, required in schemas.items():
        archive = tmp_path / f"{schema}.zip"
        payloads = {}
        for name in required - {"manifest.json"}:
            if name.endswith(".slx"):
                nested = io.BytesIO()
                with zipfile.ZipFile(nested, "w") as model:
                    model.writestr("simulink/blockdiagram.xml", "<Model/>")
                payloads[name] = nested.getvalue()
            else:
                payloads[name] = b"fixture"
        manifest = _manifest_for_schema(schema, payloads)
        with zipfile.ZipFile(archive, "w") as bundle:
            bundle.writestr("manifest.json", json.dumps(manifest))
            for name, payload in payloads.items():
                bundle.writestr(name, payload)
        assert module.verify_archive(archive)["schema"] == schema


def test_mex_handle_unlocks_happen_after_registry_critical_sections() -> None:
    gateway = _source("packages/radia-optuna/native/optuna_mex.cpp")
    helper = re.search(r"void DestroyRegistered\(.*?\n\}\n", gateway, flags=re.DOTALL)
    assert helper is not None
    helper_source = helper.group(0)
    assert helper_source.index("registry.erase") < helper_source.index("mexUnlock()")
    assert helper_source.index("if (!erased)") < helper_source.index("mexUnlock()")
    assert "while (lock_count > 0)" not in gateway



def _manifest_for_schema(schema: str, payloads: dict[str, bytes]) -> dict:
    is_preview = ".ih-release-" in schema
    version = int(schema.rsplit("v", 1)[1])
    manifest = {
        "schema": schema,
        "commit": "a" * 40,
        "matlab_release": "R2026a",
        "platform": "win64",
        "mex_extension": "mexw64",
        "required_matlab_products": ["MATLAB", "Simulink"],
        "release_channel": "preview" if is_preview else "production",
        "python_fallback": False,
        "python_per_step": False,
        "python_fallback_per_step": False,
        "operator_assembly": "preassembled",
        "files": [
            {"path": name, "size": len(payload),
             "sha256": hashlib.sha256(payload).hexdigest()}
            for name, payload in payloads.items()
        ],
    }
    if is_preview:
        manifest["backend"] = (
            "native-mex-sfunction" if version == 1
            else "matlab-level2+radia-mex-handles"
        )
        if version >= 2:
            _add_level2_manifest_fields(manifest, full=False, optuna=False)
    else:
        manifest.update({
            "package": "radia-simulink-library",
            "entry_model": "matlab/radia_simulink_library.slx",
            "backend": (
                "native-mex-sfunction-and-mex-handle"
                if version == 1 else "application-specific"
            ),
            "feature_toolbox_requirements": {
                "adjoint_topology_optimization": ["Optimization Toolbox"],
                "stream_function_topology_optimization": ["Optimization Toolbox"],
            },
        })
        if version == 1:
            manifest["required_mex"] = [
                "matlab/radia_mex.mexw64",
                "matlab/radia_ih_eddy_sfun.mexw64",
                "matlab/radia_ih_thermal_sfun.mexw64",
            ]
        else:
            manifest["feature_toolbox_requirements"][
                "electromagnet_topology_optimization"
            ] = ["Optimization Toolbox"]
            _add_level2_manifest_fields(
                manifest, full=True, optuna=version >= 3)
    return manifest


def _add_level2_manifest_fields(
        manifest: dict, *, full: bool, optuna: bool) -> None:
    required_mex = ["matlab/radia_mex.mexw64"]
    if optuna:
        required_mex.append("matlab/optuna_mex.mexw64")
    sfunctions = [
        "matlab/radia_ih_eddy_sfun.m",
        "matlab/radia_ih_thermal_sfun.m",
        "matlab/radia_ih_monitor_sfun.m",
        "matlab/+radia/+simulink/ihEddySFunction.m",
        "matlab/+radia/+simulink/ihThermalSFunction.m",
        "matlab/+radia/+simulink/ihMonitorSFunction.m",
    ]
    if full:
        sfunctions += [
            "matlab/radia_nonlinear_reactor_sfun.m",
            "matlab/+radia/+simulink/nonlinearReactorSFunction.m",
            "matlab/radia_streamfunction_optuna_sfun.m",
            "matlab/+radia/+simulink/streamFunctionOptunaSFunction.m",
        ]
        manifest["reactor_backend"] = "matlab-level2+radia-mex-handle"
        manifest["reactor_surrogate"] = False
    manifest.update({
        "ih_backend": "matlab-level2+radia-mex-handles",
        "required_mex": required_mex,
        "required_matlab_sfunctions": sfunctions,
        "standalone_mex_debug_api": True,
        "python_runtime_required_for_native_mex": True,
    })
