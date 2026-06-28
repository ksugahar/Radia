"""Reusable ipywidgets command workbench for Radia notebook panels."""

from __future__ import annotations

from dataclasses import asdict, dataclass, is_dataclass, replace
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
from pprint import pformat
import queue
import shlex
import subprocess
import sys
import threading
import time
from typing import Any, Iterable


def shell_join(argv: list[str]) -> str:
    """Return a readable command line for display in notebooks."""

    return " ".join(shlex.quote(str(a)) for a in argv)


@dataclass(frozen=True, slots=True)
class NotebookFieldSpec:
    key: str
    label: str
    kind: str = "text"
    options: tuple[Any, ...] = ()
    section: str = "Inputs"
    width: str = "300px"


@dataclass(frozen=True, slots=True)
class CommandRunRecord:
    """Small durable summary for a notebook-launched local command."""

    command: list[str]
    run_dir: Path
    log_path: Path
    result_path: Path
    status: str
    returncode: int | None
    elapsed_s: float
    started_at_utc: str
    completed_at_utc: str


class CommandWorkbench:
    """Small browser-native command panel backed by a dataclass spec."""

    title = "Radia Workbench"
    field_specs: tuple[NotebookFieldSpec, ...] = ()
    section_order: tuple[str, ...] = ()

    def __init__(
        self,
        spec: Any,
        *,
        run_root: str | Path = "runs/radia_notebook",
        timeout_s: int = 3600,
    ):
        self.spec = spec
        self.run_root = Path(run_root)
        self.timeout_s = int(timeout_s)
        self._widgets = None
        self._output = None
        self._process: subprocess.Popen[str] | None = None
        self._cancel_requested = threading.Event()
        self._run_thread: threading.Thread | None = None
        self._is_running = False
        self.last_run: CommandRunRecord | None = None

    def display(self):
        """Display the workbench and return the top-level widget."""

        try:
            import ipywidgets as W
            from IPython.display import display
        except Exception as exc:  # pragma: no cover - environment dependent
            raise RuntimeError(
                f"{type(self).__name__} requires ipywidgets in Jupyter."
            ) from exc

        self._widgets = self._build_widgets(W)
        self._output = W.Output()
        body = self._build_layout(W)
        self._wire_widgets()
        self._refresh_visibility()
        display(body)
        return body

    def to_spec(self):
        if self._widgets is None:
            return self.spec
        values = {
            field.key: self._widgets[field.key].value
            for field in self.field_specs
            if field.key in self._widgets
        }
        return replace(self.spec, **values)

    def spec_kwargs(self) -> dict[str, Any]:
        """Return the current design state as ordinary Python keyword values."""

        spec = self.to_spec()
        if is_dataclass(spec):
            return asdict(spec)
        if hasattr(spec, "__dict__"):
            return dict(vars(spec))
        return {
            field.key: getattr(spec, field.key)
            for field in self.field_specs
            if hasattr(spec, field.key)
        }

    def spec_cell_source(
        self,
        *,
        variable_name: str = "spec",
        workbench_name: str = "workbench",
    ) -> str:
        """Return a notebook cell that recreates the current initial values."""

        spec = self.to_spec()
        spec_type = type(spec)
        workbench_type = type(self)
        kwargs = pformat(self.spec_kwargs(), width=88, sort_dicts=False)
        return (
            f"from {spec_type.__module__} import {spec_type.__qualname__}\n"
            f"from {workbench_type.__module__} import {workbench_type.__qualname__}\n\n"
            f"{variable_name} = {spec_type.__qualname__}(**{kwargs})\n"
            f"{workbench_name} = {workbench_type.__qualname__}({variable_name})\n"
            f"{workbench_name}.display()\n"
        )

    def build_command(self) -> list[str]:
        self.spec = self.to_spec()
        return self.spec.build_command()

    def run(self) -> subprocess.CompletedProcess[str]:
        cmd = self.build_command()
        return subprocess.run(
            cmd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=self.timeout_s if self.timeout_s > 0 else None,
        )

    def run_local(self, *, timeout_s: int | None = None) -> CommandRunRecord:
        """Run the current command locally and save notebook-friendly artifacts."""

        cmd = self.build_command()
        self._cancel_requested.clear()
        return self._execute_to_artifacts(cmd, timeout_s=timeout_s)

    def start_background_run(self, *, timeout_s: int | None = None) -> threading.Thread:
        """Start a local command without blocking the notebook kernel."""

        if self._run_thread and self._run_thread.is_alive():
            raise RuntimeError("A command is already running.")
        cmd = self.build_command()
        self._cancel_requested.clear()
        thread = threading.Thread(
            target=self._background_entry,
            args=(cmd, timeout_s),
            name=f"{type(self).__name__}-runner",
            daemon=True,
        )
        self._run_thread = thread
        thread.start()
        return thread

    def cancel(self) -> None:
        """Request cancellation of the active background process."""

        self._cancel_requested.set()
        proc = self._process
        if proc is not None and proc.poll() is None:
            proc.terminate()

    def _build_widgets(self, W):
        widgets = {}
        for field in self.field_specs:
            value = getattr(self.spec, field.key)
            widgets[field.key] = self._make_widget(W, field, value)
        btn_layout = W.Layout(width="150px", height="34px")
        widgets["build"] = W.Button(
            description="Build command", button_style="info",
            icon="wrench", layout=btn_layout,
        )
        widgets["run"] = W.Button(
            description="Run", button_style="success",
            icon="play", layout=btn_layout,
        )
        widgets["cancel"] = W.Button(
            description="Cancel", button_style="danger",
            icon="stop", disabled=True, layout=btn_layout,
        )
        widgets["spec_cell"] = W.Button(
            description="Show spec cell", icon="code", layout=btn_layout,
        )
        cfg_style = {"description_width": "80px"}
        widgets["timeout_s"] = W.IntText(
            description="timeout [s]",
            value=self.timeout_s,
            layout=W.Layout(width="200px"),
            style=cfg_style,
        )
        widgets["run_root"] = W.Text(
            description="run root",
            value=str(self.run_root),
            layout=W.Layout(width="460px"),
            style=cfg_style,
        )
        return widgets

    def _make_widget(self, W, field: NotebookFieldSpec, value):
        layout = W.Layout(width=field.width)
        # Wide description so field labels are not truncated to "Met...".
        style = {"description_width": "150px"}
        if field.kind == "dropdown":
            options = list(field.options)
            option_values = [
                opt[1] if isinstance(opt, tuple) and len(opt) == 2 else opt
                for opt in options
            ]
            if value not in option_values and option_values:
                value = option_values[0]
            return W.Dropdown(
                description=field.label,
                options=options,
                value=value,
                layout=layout,
                style=style,
            )
        if field.kind == "checkbox":
            return W.Checkbox(
                description=field.label, value=bool(value), layout=layout, style=style
            )
        if field.kind == "int":
            return W.IntText(
                description=field.label, value=int(value), layout=layout, style=style
            )
        if field.kind == "float":
            return W.FloatText(
                description=field.label, value=float(value), layout=layout, style=style
            )
        return W.Text(
            description=field.label, value=str(value), layout=layout, style=style
        )

    def _build_layout(self, W):
        by_section: dict[str, list[Any]] = {}
        for field in self.field_specs:
            by_section.setdefault(field.section, []).append(self._widgets[field.key])
        sections = [
            *self.section_order,
            *[name for name in by_section if name not in self.section_order],
        ]
        children = [
            W.VBox(by_section[name], layout=W.Layout(gap="4px"))
            for name in sections
            if name in by_section
        ]
        accordion = W.Accordion(children=children, selected_index=0)
        for i, name in enumerate([name for name in sections if name in by_section]):
            accordion.set_title(i, name)
        header = W.HTML(
            "<div style='font-size:17px;font-weight:600;color:#202124;"
            "padding:2px 2px 8px;border-bottom:2px solid #1a73e8;"
            "margin-bottom:10px;'>"
            f"{self.title}</div>"
        )
        action_row = W.HBox(
            [self._widgets["build"], self._widgets["run"],
             self._widgets["cancel"], self._widgets["spec_cell"]],
            layout=W.Layout(gap="8px", padding="8px 0 4px"),
        )
        config_row = W.HBox(
            [self._widgets["timeout_s"], self._widgets["run_root"]],
            layout=W.Layout(gap="12px", padding="2px 0 6px"),
        )
        if self._output is not None:
            self._output.layout = W.Layout(
                border="1px solid #e0e0e0", padding="8px",
                margin="8px 0 0", max_height="340px", overflow="auto",
            )
        return W.VBox(
            [header, accordion, action_row, config_row, self._output],
            layout=W.Layout(
                border="1px solid #d0d0d0", padding="16px",
                width="840px", gap="2px",
            ),
        )

    def _wire_widgets(self) -> None:
        for field in self.field_specs:
            self._widgets[field.key].observe(
                lambda _change: self._refresh_visibility(),
                names="value",
            )
        self._widgets["build"].on_click(lambda _button: self._show_command())
        self._widgets["spec_cell"].on_click(lambda _button: self._show_spec_cell())
        self._widgets["run"].on_click(lambda _button: self._run_command())
        self._widgets["cancel"].on_click(lambda _button: self._cancel_command())

    def _visible_fields(self) -> set[str]:
        self.spec = self.to_spec()
        if hasattr(self.spec, "visible_fields"):
            return set(self.spec.visible_fields())
        return {field.key for field in self.field_specs}

    def _refresh_visibility(self) -> None:
        if self._widgets is None:
            return
        visible = self._visible_fields()
        for field in self.field_specs:
            self._widgets[field.key].layout.display = (
                "" if field.key in visible else "none"
            )
        self._widgets["run"].disabled = self._is_running or not self._is_runnable()

    def _is_runnable(self) -> bool:
        self.spec = self.to_spec()
        if hasattr(self.spec, "is_runnable"):
            return bool(self.spec.is_runnable())
        if hasattr(self.spec, "missing_required_inputs"):
            return not self.spec.missing_required_inputs()
        return True

    def _show_command(self) -> None:
        with self._output:
            self._output.clear_output()
            try:
                print(shell_join(self.build_command()))
            except Exception as exc:
                print(f"Cannot build command: {exc}")

    def _run_command(self) -> None:
        with self._output:
            self._output.clear_output()
            spec = self.to_spec()
            if hasattr(spec, "missing_required_inputs"):
                missing = spec.missing_required_inputs()
                if missing:
                    print("Missing required inputs:")
                    for item in missing:
                        print(f"  - {item}")
                    return
            try:
                cmd = self.build_command()
            except Exception as exc:
                print(f"Cannot build command: {exc}")
                return
            print(shell_join(cmd))
            print("")
            if self._run_thread and self._run_thread.is_alive():
                print("A command is already running.")
                return
            self.timeout_s = int(self._widgets["timeout_s"].value)
            self.run_root = Path(self._widgets["run_root"].value)
            self._set_running(True)
            self._cancel_requested.clear()
            self._run_thread = threading.Thread(
                target=self._background_entry,
                args=(cmd, self.timeout_s),
                name=f"{type(self).__name__}-runner",
                daemon=True,
            )
            self._run_thread.start()

    def _cancel_command(self) -> None:
        with self._output:
            print("Cancel requested.")
        self.cancel()

    def _background_entry(self, cmd: list[str], timeout_s: int | None) -> None:
        try:
            record = self._execute_to_artifacts(cmd, timeout_s=timeout_s)
            self.last_run = record
            self._emit(
                f"\n[{record.status}; exit {record.returncode}; "
                f"{record.elapsed_s:.2f}s]\n"
                f"log: {record.log_path}\n"
                f"result: {record.result_path}\n"
            )
        except Exception as exc:  # pragma: no cover - defensive UI path
            self._emit(f"\nNotebook run failed before launch: {exc}\n")
        finally:
            self._set_running(False)

    def _set_running(self, running: bool) -> None:
        self._is_running = running
        if self._widgets is None:
            return
        self._widgets["run"].disabled = running or not self._is_runnable()
        self._widgets["cancel"].disabled = not running
        self._widgets["build"].disabled = running
        self._widgets["spec_cell"].disabled = running

    def _show_spec_cell(self) -> None:
        with self._output:
            self._output.clear_output()
            print(self.spec_cell_source())

    def _emit(self, text: str) -> None:
        if self._output is None:
            sys.stdout.write(text)
            return
        with self._output:
            sys.stdout.write(text)
            if text and not text.endswith("\n"):
                sys.stdout.flush()

    def _execute_to_artifacts(
        self,
        cmd: list[str],
        *,
        timeout_s: int | None = None,
    ) -> CommandRunRecord:
        timeout = self.timeout_s if timeout_s is None else int(timeout_s)
        run_dir = self._new_run_dir()
        run_dir.mkdir(parents=True, exist_ok=False)
        log_path = run_dir / "run.log"
        result_path = run_dir / "result.json"
        command_path = run_dir / "command.txt"
        command_path.write_text(shell_join(cmd) + "\n", encoding="utf-8")

        started_at_utc = datetime.now(timezone.utc).isoformat()
        started = time.monotonic()
        status = "failed"
        returncode: int | None = None
        self._process = None

        with log_path.open("w", encoding="utf-8") as log:
            log.write(shell_join(cmd) + "\n\n")
            log.flush()
            proc = subprocess.Popen(
                cmd,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=1,
            )
            self._process = proc
            assert proc.stdout is not None
            stdout_queue: queue.Queue[str | None] = queue.Queue()

            def read_stdout() -> None:
                assert proc.stdout is not None
                for line in proc.stdout:
                    stdout_queue.put(line)
                stdout_queue.put(None)

            reader = threading.Thread(
                target=read_stdout,
                name=f"{type(self).__name__}-stdout",
                daemon=True,
            )
            reader.start()
            reader_done = False
            while True:
                while True:
                    try:
                        line = stdout_queue.get_nowait()
                    except queue.Empty:
                        break
                    if line is None:
                        reader_done = True
                        continue
                    log.write(line)
                    log.flush()
                    self._emit(line)
                if self._cancel_requested.is_set():
                    status = "cancelled"
                    self._terminate_process(proc)
                    returncode = proc.wait(timeout=5)
                    reader.join(timeout=1)
                    break
                if timeout > 0 and time.monotonic() - started > timeout:
                    status = "timeout"
                    self._terminate_process(proc)
                    returncode = proc.wait(timeout=5)
                    reader.join(timeout=1)
                    break
                returncode = proc.poll()
                if returncode is not None and reader_done:
                    status = "passed" if returncode == 0 else "failed"
                    break
                time.sleep(0.05)
            reader.join(timeout=1)
            while True:
                try:
                    line = stdout_queue.get_nowait()
                except queue.Empty:
                    break
                if line:
                    log.write(line)
                    self._emit(line)
            proc.stdout.close()

        elapsed = time.monotonic() - started
        completed_at_utc = datetime.now(timezone.utc).isoformat()
        record = CommandRunRecord(
            command=cmd,
            run_dir=run_dir,
            log_path=log_path,
            result_path=result_path,
            status=status,
            returncode=returncode,
            elapsed_s=elapsed,
            started_at_utc=started_at_utc,
            completed_at_utc=completed_at_utc,
        )
        result_path.write_text(
            json.dumps(
                self._result_payload(record, timeout_s=timeout),
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )
        self._process = None
        return record

    def _new_run_dir(self) -> Path:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        root = Path(self.run_root)
        candidate = root / stamp
        suffix = 1
        while candidate.exists():
            suffix += 1
            candidate = root / f"{stamp}_{suffix:02d}"
        return candidate

    def _result_payload(self, record: CommandRunRecord, *, timeout_s: int) -> dict[str, Any]:
        try:
            import radia

            radia_version = getattr(radia, "__version__", "unknown")
        except Exception as exc:  # pragma: no cover - import environment dependent
            radia_version = f"not-importable: {exc}"
        timing = self._timing_payload(record)
        return {
            "radia_result": {
                "schema": "radia.notebook_panel_run.v2",
                "runtime_radia_version": radia_version,
                "executed_at_utc": record.started_at_utc,
                "completed_at_utc": record.completed_at_utc,
                "panel": self.title,
                "status": record.status,
                "returncode": record.returncode,
                "elapsed_s": round(record.elapsed_s, 6),
                "timing": timing,
                "timeout_s": timeout_s,
                "command": record.command,
                "command_line": shell_join(record.command),
                "run_dir": str(record.run_dir),
                "log": str(record.log_path),
                "runtime_python": sys.version.split()[0],
                "runtime_platform": platform.platform(),
            }
        }

    def _timing_payload(self, record: CommandRunRecord) -> dict[str, Any]:
        stages = self._collect_timing_stages(record.command, record.result_path)
        stages.sort(key=lambda item: item["elapsed_s"], reverse=True)
        return {
            "wall_elapsed_s": round(record.elapsed_s, 6),
            "top_stages": stages[:4],
            "source": (
                "numeric timing keys discovered in command JSON outputs; "
                "fallback is notebook wrapper wall time"
            ),
        }

    def _collect_timing_stages(
        self,
        command: list[str],
        result_path: Path,
    ) -> list[dict[str, Any]]:
        stages: list[dict[str, Any]] = []
        for path in self._command_json_candidates(command):
            if path == result_path or not path.is_file():
                continue
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except Exception:
                continue
            for key_path, seconds in self._walk_timing_values(payload):
                stages.append({
                    "name": key_path,
                    "elapsed_s": round(seconds, 6),
                    "source": str(path),
                })
        return stages

    @staticmethod
    def _command_json_candidates(command: list[str]) -> list[Path]:
        candidates: list[Path] = []
        flags_with_path = {
            "--output",
            "--json-output",
            "--result",
            "--results",
            "--summary",
            "--summary-json",
        }
        for i, arg in enumerate(command):
            text = str(arg)
            if text in flags_with_path and i + 1 < len(command):
                candidates.append(Path(command[i + 1]))
            elif text.lower().endswith(".json"):
                candidates.append(Path(text))
        unique: list[Path] = []
        seen: set[str] = set()
        for path in candidates:
            key = str(path.resolve()) if path.exists() else str(path)
            if key not in seen:
                seen.add(key)
                unique.append(path)
        return unique

    @classmethod
    def _walk_timing_values(
        cls,
        value: Any,
        *,
        prefix: str = "",
    ) -> Iterable[tuple[str, float]]:
        if isinstance(value, dict):
            for key, item in value.items():
                key_text = str(key)
                path = f"{prefix}.{key_text}" if prefix else key_text
                if cls._is_timing_key(key_text) and isinstance(item, (int, float)):
                    seconds = float(item)
                    if seconds >= 0:
                        yield path, seconds
                yield from cls._walk_timing_values(item, prefix=path)
        elif isinstance(value, list):
            for index, item in enumerate(value):
                yield from cls._walk_timing_values(item, prefix=f"{prefix}[{index}]")

    @staticmethod
    def _is_timing_key(key: str) -> bool:
        key_l = key.lower()
        if key_l in {"elapsed_s", "wall_elapsed_s", "total_elapsed_s"}:
            return True
        if key_l.startswith("t_") and key_l.endswith("_s"):
            return True
        if key_l.endswith("_time_s") or key_l.endswith("_elapsed_s"):
            return True
        return False

    @staticmethod
    def _terminate_process(proc: subprocess.Popen[str]) -> None:
        if proc.poll() is not None:
            return
        proc.terminate()
        try:
            proc.wait(timeout=3)
        except subprocess.TimeoutExpired:
            proc.kill()


def field_keys(fields: Iterable[NotebookFieldSpec]) -> tuple[str, ...]:
    return tuple(field.key for field in fields)
