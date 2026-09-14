"""Run with a wheel-installed, isolated Python: python -I this_file --help.

No Radia/MCP import, editable install, user init file, or GUI is permitted.
Outputs remain in a new caller-selected directory, including failure evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import time
import zipfile


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--cubit-bin", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--use-installed-plugin", action="store_true",
                        help="Require deployed plugin/DLL hashes to match the wheel environment")
    parser.add_argument("--timeout", type=float, default=180)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    result = {"schema": "cubit-mesh-export.standalone-acceptance.v1",
              "passed": False, "gui_started": False,
              "host": socket.gethostname(), "python": sys.executable,
              "wheel": str(args.wheel.resolve())}
    started = time.monotonic()
    try:
        if not sys.flags.isolated or sys.prefix == sys.base_prefix:
            raise RuntimeError("Use an isolated venv interpreter with -I")
        for name in ("radia", "radia_mcp"):
            if importlib.util.find_spec(name) is not None:
                raise RuntimeError(f"Forbidden package is importable: {name}")
        result["radia_absent"] = result["radia_mcp_absent"] = True
        import cubit_mesh_export as cme
        from cubit_mesh_export import cubit_mesh_curver
        import netgen
        from cubit_mesh_export.smoke_test import _validate_exported_vol

        package = Path(cme.__file__).resolve().parent
        if not package.is_relative_to(Path(sys.prefix).resolve()):
            raise RuntimeError(f"Not installed in the isolated venv: {package}")
        result.update(version=cme.__version__, import_path=str(package),
                      wheel_sha256=sha256(args.wheel))
        result["native_import_path"] = str(Path(cubit_mesh_curver.__file__).resolve())
        result["dependencies"] = {name: metadata.version(name) for name in
                                  ("netgen-mesher", "ngsolve", "numpy")}
        # Check installed bytes against the exact wheel, not an editable tree.
        with zipfile.ZipFile(args.wheel) as archive:
            for entry in archive.namelist():
                if entry.startswith("cubit_mesh_export/") and not entry.endswith("/"):
                    installed = package / entry.removeprefix("cubit_mesh_export/")
                    if installed.read_bytes() != archive.read(entry):
                        raise RuntimeError(f"Installed wheel mismatch: {entry}")
        manifest = json.loads((package / "native_payloads.json").read_text())
        plugins = out / "plugins"
        plugins.mkdir()
        result["payloads"] = {}
        for name, expected in manifest["payloads"].items():
            source = package / name
            digest = sha256(source)
            if digest != expected["sha256"] or source.stat().st_size != expected["size"]:
                raise RuntimeError(f"Native manifest mismatch: {name}")
            shutil.copy2(source, plugins / name)
            result["payloads"][name] = digest
        for name in ("nglib.dll", "ngcore.dll"):
            shutil.copy2(Path(netgen.__file__).parent / name, plugins / name)
        console = args.cubit_bin.resolve() / "coreform_cubit.com"
        if not console.is_file():
            raise RuntimeError(f"Headless console missing: {console}")
        sample = Path(__file__).with_name("standalone_sphere.jou")
        shutil.copy2(sample, out / sample.name)
        result["journal_sha256"] = sha256(sample)
        vol = out / "sphere.vol"
        driver = out / "driver.jou"
        driver.write_text(f'play "{(out / sample.name).as_posix()}"\n'
                          f'export netgen "{vol.as_posix()}" order 2 overwrite\n'
                          'exit 0\n', encoding="utf-8")
        command = [str(console), "-batch", "-nographics", "-nojournal",
                   "-noinitfile"]
        deployed = console.parent / "plugins"
        if args.use_installed_plugin:
            for name in ("cubit_mesh_export.ccm", "nglib.dll", "ngcore.dll"):
                if sha256(deployed / name) != sha256(plugins / name):
                    raise RuntimeError(f"Deployed plugin/DLL differs from candidate: {name}")
            result["plugin_directory"] = str(deployed)
            command.extend(["-commandplugindir", str(deployed)])
        else:
            if (deployed / "cubit_mesh_export.ccm").exists():
                raise RuntimeError("Existing command plugin would be loaded twice; use "
                                   "--use-installed-plugin only for hash-identical deployment")
            command.extend(["-commandplugindir", str(plugins)])
            result["plugin_directory"] = str(plugins)
        command.append(str(driver))
        result["command"] = command
        env = os.environ.copy()
        for key in ("PYTHONPATH", "PYTHONHOME", "CUBIT_PLUGIN_DIR"):
            env.pop(key, None)
        env["PYTHONNOUSERSITE"] = "1"
        with (out / "cubit.log").open("w", encoding="utf-8") as log:
            proc = subprocess.run(command, cwd=out, env=env,
                                  stdin=subprocess.DEVNULL, stdout=log,
                                  stderr=subprocess.STDOUT, timeout=args.timeout)
        result["cubit_exit_code"] = proc.returncode
        if proc.returncode != 0:
            raise RuntimeError(f"Cubit failed with exit code {proc.returncode}")
        gate = _validate_exported_vol(vol, order=2, expect=["outer"],
                                     expect_materials=["body"], threshold=1.0)
        result["issues"] = gate["issues"]
        result["vol_sha256"] = sha256(vol)
        result["check_report"] = str(gate["report_path"])
        result["mesh"] = gate["report"]["mesh"]
        result["passed"] = gate["passed"]
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        result["elapsed_seconds"] = time.monotonic() - started
        (out / "result.json").write_text(json.dumps(result, indent=2) + "\n",
                                          encoding="utf-8")
        print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
