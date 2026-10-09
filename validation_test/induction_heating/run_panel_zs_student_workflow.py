"""Exercise existing panel-Zs wiring through the official MATLAB MCP/Toolkit.

Run on validation runtime with installed Radia native/MEX runtimes and job-local sources.
This is interface acceptance, not an independent electromagnetic accuracy test.
The model is generated in the job directory; no tracked SLX is changed.
"""
from __future__ import annotations

import argparse
import asyncio
from datetime import timedelta
import hashlib
import json
import os
from pathlib import Path
import platform
import socket
import subprocess
import sys


def activate_reviewed_sources(job: Path):
    """Use job-local Python adapters with a fixed installed native package.

    Preload native modules before putting source ahead of the wheel so a
    Python file with the same name cannot shadow its compiled extension.
    The bootstrap also applies to the assembler's Python subprocesses.
    """
    source = Path(__file__).resolve().parents[2] / "src/radia"
    if not source.is_dir():
        raise RuntimeError("Run from a staged source tree containing src/radia")
    bootstrap = job / "python-bootstrap"
    bootstrap.mkdir()
    code = (
        "from pathlib import Path\nimport importlib\nimport radia\n"
        "native = Path(radia.__file__).parent\n"
        "for extension in native.glob('*.pyd'):\n"
        "    importlib.import_module('radia.' + extension.stem)\n"
        f"radia.__path__.insert(0, {str(source)!r})\n"
    )
    (bootstrap / "sitecustomize.py").write_text(code, encoding="utf-8")
    os.environ["PYTHONPATH"] = str(bootstrap) + os.pathsep + os.environ.get("PYTHONPATH", "")
    exec(compile(code, str(bootstrap / "sitecustomize.py"), "exec"), {})


def prepare(job: Path):
    import ngsolve as ng
    import numpy as np
    from test_ih_operator_physics_golden import _assemble
    from radia.surface_impedance import main as template_cli

    ng.SetNumThreads(4)
    with ng.TaskManager():
        fixture = _assemble(job)
        em = json.loads((job / "run/electromagnetic_result.json").read_text())
        # Use the actual scalar solve's impedance, not a new material model.
        zr, zi = em["Z_s_wp_real"], em["Z_s_wp_imag"]
        template_cli([str(fixture["workpiece"]), "--frequency", "7000",
                      "--real-ohm", str(zr), "--imag-ohm", str(zi),
                      "--output", str(job / "zs-template.json")])
    table = json.loads((job / "zs-template.json").read_text())
    assert len(table["real_ohm"]) == len(table["centroids_m"])
    assert np.isfinite(zr) and zr > 0
    return fixture


async def exercise(args, job):
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    source = Path(__file__).resolve().parents[2]
    events = []
    # A new mode is intentional only when process/Engine admission is empty.
    inventory = subprocess.run(
        ["pwsh", "-NoProfile", "-Command",
         "@(Get-Process MATLAB -ErrorAction SilentlyContinue).Count"],
        check=True, capture_output=True, text=True)
    import matlab.engine
    if int(inventory.stdout.strip()) or matlab.engine.find_matlab():
        raise RuntimeError("An existing MATLAB must be explicitly reused; no new session started")
    parameters = StdioServerParameters(
        command=str(args.matlab_mcp),
        args=["--matlab-session-mode=new", "--matlab-display-mode=nodesktop",
              "--initialize-matlab-on-startup=false", "--disable-telemetry=true",
              f"--initial-working-folder={job}", f"--log-folder={job / 'mcp-logs'}",
              f"--extension-file={args.toolkit_tools}"],
        env=dict(os.environ))
    with (job / "mcp-stderr.log").open("w", encoding="utf-8") as errors:
        async with stdio_client(parameters, errlog=errors) as (read, write):
            async with ClientSession(read, write,
                                     read_timeout_seconds=timedelta(minutes=15)) as session:
                await session.initialize()

                async def call(name, arguments):
                    result = await session.call_tool(name, arguments)
                    payload = result.model_dump(mode="json")
                    events.append({"tool": name, "arguments": arguments, "result": payload})
                    (job / "tool-events.json").write_text(json.dumps(events, indent=2))
                    text = "\n".join(c.text for c in result.content if c.type == "text")
                    if result.isError or "Error using " in text or "Error in " in text:
                        raise RuntimeError(f"{name}: {text}")
                    return text

                async def evaluate(code):
                    # The server can return localized MATLAB errors with isError=false.
                    # Require an explicit completion marker rather than parsing errors.
                    marker = "RADIA_WORKFLOW_STEP_COMPLETE"
                    output = await call("evaluate_matlab_code", {
                        "code": code + f" fprintf('{marker}\\n');"})
                    if marker not in output:
                        raise RuntimeError(f"MATLAB step did not complete: {output}")
                    return output

                def quote(value):
                    return "'" + str(value).replace("'", "''") + "'"

                # Read library configuration before any model operation.
                await evaluate("gate=library.settingsLookup(); disp(jsonencode(gate)); "
                               "assert(~gate.found || gate.gatePass,'Library configuration needs setup');")
                await evaluate(
                    f"addpath({quote(args.mex_dir)}); addpath({quote(source / 'matlab')}); "
                    f"runtime=radia.setup(PythonExecutable={quote(sys.executable)},RequireMex=true); "
                    "disp(jsonencode(runtime)); "
                    f"job={quote(job)}; model='panel_zs_student'; "
                    "fprintf('HOST=%s PID=%d MATLAB=%s MEX=%s\\n',getenv('COMPUTERNAME'),feature('getpid'),version,which('radia_mex')); "
                    "table=jsondecode(fileread(fullfile(job,'zs-template.json'))); "
                    "zs=complex(table.real_ohm,table.imag_ohm); "
                    "topFace=abs(table.centroids_m(:,3)-.01)<1e-9; assert(any(topFace)); "
                    "zs(topFace)=2*zs(topFace); "
                    "radia.simulink.writeIHPanelImpedance(fullfile(job,'zs-template.json'),fullfile(job,'zs.json'),zs); "
                    "radia.simulink.buildIHNativeModel(ModelName=model,OutputDirectory=job,Open=true); "
                    "blocks=find_system(model,'Tag','RadiaIHGeometryUpdate'); "
                    "sid=Simulink.ID.getSID(blocks{1}); fprintf('GEOMETRY_SID=%s\\n',sid);")
                await call("model_read", {"model": "panel_zs_student", "scope": "root", "depth": "0"})
                sid_text = await evaluate("fprintf('SID_NUMBER=%s\\n',extractAfter(sid,':'));")
                import re
                match = re.search(r"SID_NUMBER=(\d+)", sid_text)
                if not match:
                    raise RuntimeError("Cannot resolve Geometry Update SID")
                target = "blk_" + match.group(1)

                async def configure(values):
                    output = await call("model_edit", {
                        "model": "panel_zs_student", "scope": "root",
                        "layout_mode": "incremental", "operations": json.dumps([
                            {"op": "configure", "target": target, "params": values}])})
                    if "status: ok" not in output:
                        raise RuntimeError(f"Toolkit configuration did not complete: {output}")

                await configure({"wp_vol": str(job / "workpiece.vol"),
                                 "coil_file": str(job / "gapped_torus_coil.step"),
                                 "config_file": str(job / "panel-native.json"),
                                 "python_executable": sys.executable,
                                 "panel_zs_file": (job / "zs.json").as_posix()})
                await evaluate(
                    "first=radia.simulink.updateIHGeometry(model); assert(first.rebuilt); "
                    "cfg=jsondecode(fileread(fullfile(job,'panel-native.json'))); "
                    "assert(string(cfg.surface_impedance.mode)=='specified-panel'); "
                    "assert(string(cfg.independent_fem_validation.status)=='not-performed'); "
                    "assert(cfg.unit_current.relative_power_error<=1e-8); "
                    "emPath=fullfile(job,'panel-native_artifacts','electromagnetic_result.json'); "
                    "panelEm=jsondecode(fileread(emPath)); assert(string(panelEm.impedance_model)=='specified-panel'); "
                    "usedZs=complex(panelEm.esim_per_panel_Z_s_real,panelEm.esim_per_panel_Z_s_imag); "
                    "assert(isequal(size(usedZs),size(zs)) && max(abs(usedZs-zs))<=1e-12*max(abs(zs))); "
                    "firstHash=string(cfg.surface_impedance.file_sha256); "
                    "fresh=radia.simulink.updateIHGeometry(model); "
                    "assert(~fresh.rebuilt && string(fresh.reason)=='up-to-date'); "
                    "workspace=get_param(model,'ModelWorkspace'); "
                    "workspace.assignin('radia_ih_geometry_revision',NaN); "
                    "reload=radia.simulink.updateIHGeometry(model); "
                    "assert(~reload.rebuilt && string(reload.reason)=='reloaded');")
                (job / "sub").mkdir(exist_ok=True)
                await configure({"panel_zs_file": str(job / "sub/../zs.json")})
                await evaluate(
                    "alias=radia.simulink.updateIHGeometry(model); assert(~alias.rebuilt); "
                    "save_system(model,fullfile(job,[model '.slx'])); close_system(model,0); "
                    "open_system(fullfile(job,[model '.slx'])); "
                    "reopened=radia.simulink.updateIHGeometry(model); assert(~reopened.rebuilt); "
                    "out=sim(model,'StopTime','1',ReturnWorkspaceOutputs='on'); "
                    "monitor=out.get('radia_ih_monitor'); assert(all(isfinite(monitor),'all')); "
                    "columns=strsplit(get_param([model '/IH Dashboard/Select Telemetry'],'OutputSignals'),','); "
                    "assert(size(monitor,2)==numel(columns),'IH monitor layout changed'); "
                    "meanIndex=find(strcmp(columns,'temperature_mean_K')); assert(isscalar(meanIndex)); "
                    "assert(monitor(end,meanIndex)>293.15); "
                    "radia.simulink.writeIHPanelImpedance(fullfile(job,'zs-template.json'),fullfile(job,'zs.json'),1.1*zs); "
                    "changed=radia.simulink.updateIHGeometry(model); assert(changed.rebuilt); "
                    "newcfg=jsondecode(fileread(fullfile(job,'panel-native.json'))); "
                    "assert(string(newcfg.surface_impedance.file_sha256)~=firstHash); "
                    "changedEm=jsondecode(fileread(emPath)); assert(string(changedEm.impedance_model)=='specified-panel'); "
                    "usedChangedZs=complex(changedEm.esim_per_panel_Z_s_real,changedEm.esim_per_panel_Z_s_imag); "
                    "assert(isequal(size(usedChangedZs),size(zs)) && max(abs(usedChangedZs-1.1*zs))<=1e-12*max(abs(zs)));")
                await configure({"panel_zs_file": ""})
                await evaluate(
                    "cleared=radia.simulink.updateIHGeometry(model); assert(cleared.rebuilt); "
                    "uniform=jsondecode(fileread(fullfile(job,'panel-native.json'))); "
                    "assert(string(uniform.surface_impedance.mode)=='uniform-linear'); "
                    "uniformEm=jsondecode(fileread(emPath)); assert(~uniformEm.esim_per_panel); "
                    "assert(abs(panelEm.P_wp_W-uniformEm.P_wp_W)>1e-10*max(panelEm.P_wp_W,uniformEm.P_wp_W)); "
                    "result=struct('status','passed','scope','existing panel-Zs interface wiring only', "
                    "'first',first,'fresh',fresh,'reload',reload,'alias',alias,'reopened',reopened, "
                    "'changed',changed,'cleared',cleared,'thermal_power_relative_error',cfg.unit_current.relative_power_error, "
                    "'independent_fem_validation','not-performed','native_monitor',monitor, "
                    "'monitor_columns',{columns},'panel_values_verified',true,'panel_count',numel(zs), "
                    "'panel_power_W',panelEm.P_wp_W,'changed_panel_power_W',changedEm.P_wp_W, "
                    "'uniform_power_W',uniformEm.P_wp_W); "
                    "fid=fopen(fullfile(job,'result.json'),'w'); cleanup=onCleanup(@()fclose(fid)); "
                    "fprintf(fid,'%s',jsonencode(result)); clear cleanup; close_system(model,0);")
                check = await call("model_check", {"model": str(job / "panel_zs_student.slx"),
                                                   "scope": "root", "checks": '["all"]'})
                if "status: healthy" not in check:
                    raise RuntimeError(f"Toolkit structural check failed: {check}")
                result = json.loads((job / "result.json").read_text())
                result["toolkit_structural_check"] = "healthy"
                (job / "result.json").write_text(json.dumps(result, indent=2))
    return events


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", type=Path, required=True)
    parser.add_argument("--matlab-mcp", type=Path, required=True)
    parser.add_argument("--toolkit-tools", type=Path, required=True)
    parser.add_argument("--mex-dir", type=Path, required=True)
    args = parser.parse_args()
    job = args.job.resolve()
    if not args.job.is_absolute() or job.exists():
        raise ValueError("Use a fresh absolute job directory")
    job.mkdir(parents=True)
    provenance = dict(platform_class=platform.system(), python=platform.python_version(),
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (job / "provenance.json").write_text(json.dumps(provenance, indent=2))
    try:
        activate_reviewed_sources(job)
        from importlib.metadata import version
        root = Path(__file__).resolve().parents[2]
        provenance["packages"] = {name: version(name) for name in
                                  ("radia", "ngsolve", "netgen-mesher", "mcp")}
        provenance["inputs"] = {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (args.matlab_mcp, args.toolkit_tools,
                         args.mex_dir / "radia_mex.mexw64",
                         args.mex_dir / "radia_mex.mexw64.build.json")
        }
        provenance["source_files"] = {
            path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for directory, suffix in ((root / "src/radia", "*.py"), (root / "matlab", "*.m"))
            for path in directory.rglob(suffix)
        }
        (job / "provenance.json").write_text(json.dumps(provenance, indent=2))
        prepare(job)
        asyncio.run(exercise(args, job))
        if not (job / "result.json").is_file():
            raise RuntimeError("Workflow did not produce its acceptance result")
    except Exception as exc:
        (job / "result.json").unlink(missing_ok=True)
        (job / "failure.json").write_text(json.dumps({"status": "failed", "error": str(exc)}))
        raise


if __name__ == "__main__":
    main()
