"""Hash-locked, isolated Windows compute runtimes and foreground jobs.

This is an operations helper, not a solver or release/publishing command.
Use the same wheelhouse/lock on every host. No package is installed globally.
"""
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import importlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import runpy
import shutil
import subprocess
import sys
import time
import venv
import zipfile
from email.parser import BytesParser

SCHEMA = "radia.compute-runtime.v1"
MODULES = ("radia", "radia._radia_pybind", "radia.bh_law",
           "radia.coil_builder", "radia.vim._solve", "radia.vim._nonlinear",
           "ngsolve", "netgen", "numpy", "scipy")
THREADS = {"OMP_NUM_THREADS": "8", "MKL_NUM_THREADS": "1",
           "OPENBLAS_NUM_THREADS": "1", "MKL_THREADING_LAYER": "TBB"}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def normalize(name: str) -> str:
    return name.lower().replace("_", "-").replace(".", "-")


def read_lock(path: Path) -> dict:
    lock = json.loads(path.read_text(encoding="utf-8"))
    payload = {k: v for k, v in lock.items() if k != "lock_sha256"}
    if lock.get("schema") != SCHEMA or hashlib.sha256(canonical(payload)).hexdigest() != lock.get("lock_sha256"):
        raise ValueError("Invalid or modified compute lock")
    return lock


def create_lock(wheelhouse: Path, output: Path, procedure_commit: str, runtime_id: str | None = None) -> dict:
    wheels, modules, records = [], {}, {}
    for path in sorted(wheelhouse.glob("*.whl")):
        with zipfile.ZipFile(path) as z:
            meta_names = [n for n in z.namelist() if n.endswith(".dist-info/METADATA")]
            if len(meta_names) != 1:
                raise ValueError(f"Ambiguous wheel metadata: {path}")
            meta = BytesParser().parsebytes(z.read(meta_names[0]))
            wheels.append(dict(file=path.name, sha256=digest(path),
                               name=normalize(meta["Name"]), version=meta["Version"]))
            record = meta_names[0].replace("METADATA", "RECORD")
            for name, value, _ in csv.reader(z.read(record).decode().splitlines()):
                if value.startswith("sha256="):
                    records.setdefault(name, set()).add(value.split("=", 1)[1])
            for module in MODULES:
                stem = module.replace(".", "/")
                candidates = [n for n in z.namelist() if n in (stem + ".py", stem + "/__init__.py")
                              or (n.startswith(stem + ".") and n.endswith(".pyd"))]
                if candidates:
                    if len(candidates) != 1 or module in modules:
                        raise ValueError(f"Ambiguous module {module}")
                    name = candidates[0]
                    modules[module] = dict(wheel=path.name, sha256=hashlib.sha256(z.read(name)).hexdigest())
    if set(MODULES) != set(modules):
        raise ValueError(f"Missing required wheel modules: {set(MODULES) - set(modules)}")
    names = [w["name"] for w in wheels]
    if len(set(names)) != len(names):
        raise ValueError("More than one wheel per distribution")
    lock = dict(schema=SCHEMA, python=platform.python_version(), architecture=platform.machine(),
                procedure_commit=procedure_commit, helper_sha256=digest(Path(__file__)),
                threads=THREADS, wheels=wheels, modules=modules)
    # Netgen and NGSolve both own share/__init__.py, with different delvewheel paths.
    # Accept only exact bytes present in the locked wheels and record the installed winner.
    lock["shared_record_hashes"] = {n: sorted(v) for n, v in records.items() if len(v) > 1}
    # Numerical environment identity is independent of operation-helper revisions.
    lock["runtime_id"] = runtime_id or hashlib.sha256(canonical(dict(
        wheels=wheels, python=lock["python"], architecture=lock["architecture"]))).hexdigest()[:16]
    if len(lock["runtime_id"]) != 16 or any(c not in "0123456789abcdef" for c in lock["runtime_id"]):
        raise ValueError("Runtime ID must contain 16 lowercase hexadecimal characters")
    lock["lock_sha256"] = hashlib.sha256(canonical(lock)).hexdigest()
    output.write_text(json.dumps(lock, indent=2), encoding="utf-8")
    return lock


def verify_wheels(lock: dict, wheelhouse: Path) -> None:
    expected = {row["file"] for row in lock["wheels"]}
    if {p.name for p in wheelhouse.glob("*.whl")} != expected:
        raise ValueError("Wheelhouse inventory differs from lock")
    for row in lock["wheels"]:
        if digest(wheelhouse / row["file"]) != row["sha256"]:
            raise ValueError(f"Wheel hash mismatch: {row['file']}")


def probe(lock: dict) -> dict:
    errors, loaded, packages, shared = [], {}, {}, {}
    prefix = Path(sys.prefix).resolve()
    os.environ.update(environment(lock, prefix))
    if sys.prefix == sys.base_prefix:
        errors.append("A dedicated virtual environment is required")
    if platform.python_version() != lock["python"] or platform.machine() != lock["architecture"]:
        errors.append("Python patch version or architecture differs")
    if digest(Path(__file__)) != lock["helper_sha256"]:
        errors.append("Execution helper differs from lock")
    allowed = {w["name"] for w in lock["wheels"]} | {"pip", "setuptools"}
    installed = {normalize(d.metadata["Name"]): d for d in metadata.distributions()}
    if set(installed) - allowed:
        errors.append("Unexpected installed distributions: " + str(sorted(set(installed) - allowed)))
    for row in lock["wheels"]:
        name = row["name"]
        try:
            dist = installed[name]
            packages[name] = dist.version
            if dist.version != row["version"]:
                errors.append(f"Version mismatch: {name}")
            # Check installed RECORD hashes, including numerical DLLs, not only version strings.
            if not dist.files:
                errors.append(f"Missing installed file records: {name}")
            for member in dist.files or []:
                if member.hash is None:
                    continue
                p = Path(dist.locate_file(member)).resolve()
                if not p.is_relative_to(prefix) or not p.is_file():
                    errors.append(f"Missing/outside-runtime file: {name}/{member}")
                    continue
                if member.hash.mode != "sha256":
                    errors.append(f"Unsupported RECORD hash: {name}/{member}")
                    continue
                actual = base64.urlsafe_b64encode(bytes.fromhex(digest(p))).decode().rstrip("=")
                known_shared = lock.get("shared_record_hashes", {}).get(str(member).replace("\\", "/"), [])
                if known_shared:
                    shared[str(member).replace("\\", "/")] = actual
                if actual != member.hash.value and actual not in known_shared:
                    errors.append(f"Installed file changed: {name}/{member}")
        except KeyError:
            errors.append(f"Missing distribution: {name}")
    for name, expected in lock["modules"].items():
        try:
            module = importlib.import_module(name)
            path = Path(module.__file__).resolve()
            loaded[name] = dict(path=str(path), sha256=digest(path))
            if not path.is_relative_to(prefix) or loaded[name]["sha256"] != expected["sha256"]:
                errors.append(f"Loaded source/binary mismatch: {name}")
        except Exception as exc:
            errors.append(f"Import failed: {name}: {exc!r}")
    import psutil
    numerical_dlls = {}
    for mapping in psutil.Process().memory_maps():
        path = Path(mapping.path)
        name = path.name.lower()
        if name.endswith(".dll") and any(s in name for s in ("mkl", "openblas", "ngcore", "nglib", "ngstd", "libngsolve", "tbb12", "libiomp")):
            path = path.resolve()
            numerical_dlls[name] = dict(path=str(path), sha256=digest(path))
            if not path.is_relative_to(prefix):
                errors.append(f"Numerical DLL loaded outside selected runtime: {path}")
    return dict(schema=SCHEMA, status="passed" if not errors else "failed", errors=errors,
                host=platform.node(), python=sys.executable, python_version=platform.python_version(),
                prefix=str(prefix), lock_sha256=lock["lock_sha256"], procedure_commit=lock["procedure_commit"],
                packages=packages, modules=loaded, threads=lock["threads"],
                shared_installed_files=shared, numerical_dlls=numerical_dlls,
                mklroot=os.environ["MKLROOT"])


def environment(lock: dict, runtime: Path | None = None) -> dict:
    env = os.environ.copy()
    for key in ("PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    # The supported MKLROOT override must be absolute and belong to this runtime.
    # Without it Radia's loader can mistake a cwd-relative bin directory for MKL.
    env.update(lock["threads"], MKLROOT=str((runtime or Path(sys.prefix)).resolve() / "Library"),
               PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
    return env


def deploy(lock_path: Path, wheelhouse: Path, runtime: Path) -> None:
    lock = read_lock(lock_path)
    verify_wheels(lock, wheelhouse)
    if platform.python_version() != lock["python"]:
        raise ValueError("Bootstrap Python patch version differs")
    # A deployed runtime is maintained software, never a case workspace.
    base = Path("C:/ProgramData/Radia/compute-runtimes").resolve()
    runtime = runtime.resolve()
    if runtime.parent != base or runtime.name != lock["runtime_id"]:
        raise ValueError("Use the dedicated runtime path recorded in the lock")
    if runtime.exists() and runtime.is_symlink():
        raise ValueError("Runtime must not be a link")
    if (runtime / "lock.json").exists():
        previous = read_lock(runtime / "lock.json")
        if any(previous[k] != lock[k] for k in ("wheels", "python", "architecture")):
            raise ValueError("An existing runtime cannot receive different numerical wheels")
    else:
        runtime.mkdir(parents=True, exist_ok=True)
        venv.EnvBuilder(with_pip=True).create(runtime)
    python = runtime / "Scripts/python.exe"
    req = runtime / "requirements.lock"
    req.write_text("\n".join(f"{r['name']}=={r['version']} --hash=sha256:{r['sha256']}" for r in lock["wheels"]) + "\n")
    subprocess.run([str(python), "-I", "-m", "pip", "install", "--no-index", "--no-deps",
                    "--only-binary=:all:", "--require-hashes", "--find-links", str(wheelhouse.resolve()),
                   "-r", str(req)], env=environment(lock, runtime), check=True)
    subprocess.run([str(python), "-I", "-m", "pip", "check"], env=environment(lock, runtime), check=True)
    shutil.copy2(lock_path, runtime / "lock.json")
    if Path(__file__).resolve() != runtime / "compute_runtime.py":
        shutil.copy2(__file__, runtime / "compute_runtime.py")
    subprocess.run([str(python), "-I", str(runtime / "compute_runtime.py"), "probe", "--lock",
                    str(runtime / "lock.json"), "--output", str(runtime / "acceptance.json")],
                   cwd=runtime, env=environment(lock, runtime), check=True)


def safe_member(root: Path, name: str) -> Path:
    path = (root / name).resolve()
    if not path.is_relative_to(root.resolve()) or Path(name).is_absolute():
        raise ValueError(f"Escaping bundle member: {name}")
    return path


def process_tree_private_bytes(process) -> int:
    """Windows venv launchers spawn the real interpreter; measure the whole tree."""
    import psutil
    total = 0
    for member in [process, *process.children(recursive=True)]:
        try:
            memory = member.memory_info()
            total += getattr(memory, "private", memory.rss)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return total


def verify_job(root: Path, job: dict, lock: dict) -> None:
    if job.get("schema") != "radia.compute-job.v1" or job.get("lock_sha256") != lock["lock_sha256"]:
        raise ValueError("Job/runtime lock mismatch")
    if not job.get("source_commit") or not job.get("physics_contract"):
        raise ValueError("Commit and physical assumptions must be recorded")
    for name, expected in job["files"].items():
        p = safe_member(root, name)
        if p.is_symlink() or not p.is_file() or digest(p) != expected:
            raise ValueError(f"Job file changed/missing: {name}")
    if job["entry"] not in job["files"]:
        raise ValueError("Entrypoint must be hashed")


def run_job(lock: dict, root: Path) -> int:
    report = probe(lock)
    if report["status"] != "passed":
        raise RuntimeError(json.dumps(report))
    root = root.resolve()
    if not root.is_relative_to(Path("C:/temp").resolve()) or root == Path("C:/temp").resolve():
        raise ValueError("Job must have its own C:/temp directory")
    job = json.loads((root / "job.json").read_text())
    verify_job(root, job, lock)
    (root / "runtime.json").write_text(json.dumps(report, indent=2))
    import psutil
    active = []
    for process in psutil.process_iter(["pid", "name", "cmdline"]):
        if process.pid == os.getpid():
            continue
        name = (process.info["name"] or "").lower()
        if name == "runner.worker.exe" or (job.get("workload", "analysis") != "smoke" and name in ("python.exe", "matlab.exe")):
            active.append(process.info)
    if active:
        raise RuntimeError("Compute admission refused: CI/other computation is active: " + json.dumps(active))
    free = psutil.virtual_memory().available
    if free < int(job["minimum_free_memory_bytes"]):
        raise RuntimeError("Insufficient free memory for declared job budget")
    # A foreground child keeps the SSH connection as the lifetime boundary.
    started = time.perf_counter()
    command = [sys.executable, "-I", "-u", str(Path(__file__).resolve()), "entry",
               "--lock", str(Path(sys.prefix) / "lock.json"), "--job-root", str(root)]
    with (root / "run.log").open("w", encoding="utf-8") as log:
        child = subprocess.Popen(command, cwd=root, env=environment(lock), stdout=log, stderr=subprocess.STDOUT)
        peak = 0
        try:
            while child.poll() is None:
                try:
                    peak = max(peak, process_tree_private_bytes(psutil.Process(child.pid)))
                except psutil.NoSuchProcess:
                    pass
                time.sleep(1)
        finally:
            if child.poll() is None:
                child.terminate()
                child.wait(timeout=30)
    result = dict(status="finished" if child.returncode == 0 else "failed", returncode=child.returncode,
                  host=platform.node(), command=command, lock_sha256=lock["lock_sha256"],
                  source_commit=job["source_commit"], wall_s=time.perf_counter()-started,
                  peak_process_private_bytes=peak)
    result["memory_scope"] = "Sum of child process tree private bytes, including Windows venv interpreter child; one-second samples"
    (root / "execution.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)
    return child.returncode


def execute_entry(lock: dict, root: Path) -> None:
    report = probe(lock)
    if report["status"] != "passed":
        raise RuntimeError(json.dumps(report))
    job = json.loads((root / "job.json").read_text())
    verify_job(root, job, lock)
    # Bundle helpers are visible after installed distributions, preventing source shadowing.
    sys.path.append(str(root.resolve()))
    sys.argv = [str(safe_member(root, job["entry"])), *job.get("args", [])]
    runpy.run_path(sys.argv[0], run_name="__main__")


def pack_job(spec_path: Path, lock: dict, output: Path) -> dict:
    """Export code from declared Git commits, never from uncommitted trial files."""
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    contents, provenance = {}, []
    for source in spec["sources"]:
        repo, commit = Path(source["repository"]), source["commit"]
        sha = subprocess.check_output(["git", "rev-parse", "--verify", commit + "^{commit}"], cwd=repo, text=True).strip()
        hashes = {}
        for target, tracked in source["files"].items():
            safe_member(Path("C:/temp/bundle-check"), target)
            if target in contents or target == "job.json":
                raise ValueError(f"Duplicate/reserved bundle member: {target}")
            contents[target] = subprocess.check_output(["git", "show", f"{sha}:{tracked}"], cwd=repo)
            hashes[tracked] = hashlib.sha256(contents[target]).hexdigest()
        provenance.append(dict(commit=sha, tracked_files=hashes))
    for target, filename in spec.get("inputs", {}).items():
        safe_member(Path("C:/temp/bundle-check"), target)
        if target in contents or target == "job.json":
            raise ValueError(f"Duplicate/reserved bundle member: {target}")
        contents[target] = Path(filename).read_bytes()
    job = dict(schema="radia.compute-job.v1", lock_sha256=lock["lock_sha256"],
               source_commit=provenance[0]["commit"], sources=provenance,
               entry=spec["entry"], args=spec.get("args", []),
               physics_contract=spec["physics_contract"],
               minimum_free_memory_bytes=spec["minimum_free_memory_bytes"],
               workload=spec.get("workload", "analysis"),
               files={n: hashlib.sha256(v).hexdigest() for n, v in contents.items()})
    if job["entry"] not in job["files"]:
        raise ValueError("Missing hashed entrypoint")
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as z:
        for name, value in contents.items():
            z.writestr(name, value)
        z.writestr("job.json", json.dumps(job, indent=2))
    return dict(archive=str(output), sha256=digest(output), job=job)


def recover(root: Path, output: Path) -> dict:
    root = root.resolve()
    if not root.is_relative_to(Path("C:/temp").resolve()) or root == Path("C:/temp").resolve():
        raise ValueError("Recovery requires an owned C:/temp child directory")
    files = {}
    for p in root.rglob("*"):
        if p.is_symlink() or getattr(p.lstat(), "st_file_attributes", 0) & 0x400:
            raise ValueError("Recovery refuses links/junctions")
        if p.is_file():
            files[p.relative_to(root).as_posix()] = digest(p)
    manifest = dict(root=str(root), host=platform.node(), files=files)
    if output.resolve().is_relative_to(root):
        raise ValueError("Recovery archive must be outside the job directory")
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as z:
        for name in files:
            z.write(root / name, name)
        z.writestr("recovery_manifest.json", json.dumps(manifest, indent=2))
    return dict(archive=str(output), archive_sha256=digest(output), root=str(root), host=platform.node())


def verify_recovery(archive: Path, destination: Path, expected: str) -> dict:
    if digest(archive) != expected:
        raise ValueError("Recovery archive hash mismatch")
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for member in z.infolist():
            safe_member(destination, member.filename)
            if (member.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError("Recovery archive contains a link")
        z.extractall(destination)
    manifest = json.loads((destination / "recovery_manifest.json").read_text())
    for name, expected_hash in manifest["files"].items():
        if digest(safe_member(destination, name)) != expected_hash:
            raise ValueError(f"Recovered file hash mismatch: {name}")
    return dict(status="verified", archive_sha256=expected, root=manifest["root"],
                host=manifest["host"], recovered_files=len(manifest["files"]), destination=str(destination))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="action", required=True)
    p = sub.add_parser("lock")
    p.add_argument("--wheelhouse", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--procedure-commit", required=True)
    p.add_argument("--runtime-id", help="Reuse a named runtime only when its numerical wheel identity is unchanged")
    p = sub.add_parser("deploy")
    p.add_argument("--lock", type=Path, required=True)
    p.add_argument("--wheelhouse", type=Path, required=True)
    p.add_argument("--runtime", type=Path, required=True)
    p = sub.add_parser("probe")
    p.add_argument("--lock", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    for name in ("run", "entry"):
        p = sub.add_parser(name)
        p.add_argument("--lock", type=Path, required=True)
        p.add_argument("--job-root", type=Path, required=True)
    p = sub.add_parser("pack")
    p.add_argument("--lock", type=Path, required=True)
    p.add_argument("--spec", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("recover")
    p.add_argument("--job-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("verify-recovery")
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--sha256", required=True)
    p.add_argument("--destination", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    a = ap.parse_args()
    if a.action == "lock":
        print(json.dumps(create_lock(a.wheelhouse, a.output, a.procedure_commit, a.runtime_id)))
    elif a.action == "deploy":
        deploy(a.lock, a.wheelhouse, a.runtime)
    elif a.action == "probe":
        result = probe(read_lock(a.lock))
        a.output.write_text(json.dumps(result, indent=2))
        print(json.dumps(dict(status=result["status"], host=result["host"], errors=result["errors"])))
        return 0 if result["status"] == "passed" else 2
    elif a.action == "entry":
        execute_entry(read_lock(a.lock), a.job_root)
    elif a.action == "pack":
        print(json.dumps(pack_job(a.spec, read_lock(a.lock), a.output)))
    elif a.action == "recover":
        print(json.dumps(recover(a.job_root, a.output)))
    elif a.action == "verify-recovery":
        report = verify_recovery(a.archive, a.destination, a.sha256)
        a.report.write_text(json.dumps(report, indent=2))
        print(json.dumps(report))
    else:
        return run_job(read_lock(a.lock), a.job_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
