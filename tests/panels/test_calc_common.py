"""Unit tests for panels/calc_common.py (no Cubit/NGSolve required).

Tests:
  - progress() writes to stderr in expected format
  - calc_main() handles success, error, --output file
"""

import argparse
import io
import json
import os
import sys
import tempfile

import pytest

# Ensure panels directory is on path
_panels_dir = os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia", "panels")
_panels_dir = os.path.abspath(_panels_dir)
if _panels_dir not in sys.path:
    sys.path.insert(0, _panels_dir)

from calc_common import progress, calc_main


class TestProgress:
    def test_writes_to_stderr(self, capsys):
        progress("TEST", "hello world")
        captured = capsys.readouterr()
        assert captured.out == ""
        assert "TEST:hello world\n" in captured.err

    def test_tag_and_message(self, capsys):
        progress("MESH_READY", "42 elements")
        captured = capsys.readouterr()
        assert "MESH_READY:42 elements\n" in captured.err

    def test_empty_message(self, capsys):
        progress("TAG", "")
        captured = capsys.readouterr()
        assert "TAG:\n" in captured.err


class TestCalcMain:
    def _run_calc_main(self, solve_func, extra_args=None, cli_args=None):
        """Helper: run calc_main capturing stdout/stderr."""
        parser = argparse.ArgumentParser()
        if extra_args:
            for arg_spec in extra_args:
                parser.add_argument(*arg_spec[0], **arg_spec[1])

        old_argv = sys.argv
        sys.argv = ["test"] + (cli_args or [])

        old_stdout = sys.stdout
        old_stderr = sys.stderr
        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()

        try:
            sys.stdout = stdout_buf
            sys.stderr = stderr_buf
            calc_main(solve_func, parser)
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
            sys.argv = old_argv

        stdout_text = stdout_buf.getvalue()
        stderr_text = stderr_buf.getvalue()
        return stdout_text, stderr_text

    def test_success_json_on_stdout(self):
        def solve(args):
            return {"value": 42}

        stdout, stderr = self._run_calc_main(solve)
        result = json.loads(stdout.strip())
        assert result["value"] == 42
        # calc_main stamps the radia version into every result dict
        # (auditability of the persisted .json).
        assert "radia_version" in result
        assert stderr == ""

    def test_error_returns_json(self):
        def solve(args):
            raise ValueError("test error")

        stdout, stderr = self._run_calc_main(solve)
        result = json.loads(stdout.strip())
        assert "error" in result
        assert "test error" in result["error"]
        # Traceback on stderr
        assert "ValueError" in stderr

    def test_output_file(self):
        def solve(args):
            return {"written": True}

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            tmp_path = f.name

        try:
            stdout, _ = self._run_calc_main(
                solve, cli_args=["--output", tmp_path])
            # JSON on stdout
            result_stdout = json.loads(stdout.strip())
            assert result_stdout["written"] is True
            # JSON in file
            with open(tmp_path) as f:
                result_file = json.load(f)
            assert result_file["written"] is True
        finally:
            os.unlink(tmp_path)

    def test_auto_adds_output_arg(self):
        """calc_main adds --output if parser doesn't have it."""
        parser = argparse.ArgumentParser()
        # No --output argument added

        def solve(args):
            return {"ok": True}

        old_argv = sys.argv
        sys.argv = ["test"]
        old_stdout = sys.stdout
        buf = io.StringIO()
        try:
            sys.stdout = buf
            calc_main(solve, parser)
        finally:
            sys.stdout = old_stdout
            sys.argv = old_argv

        result = json.loads(buf.getvalue().strip())
        assert result["ok"] is True

    def test_custom_args_passed(self):
        """Solver receives parsed args."""
        received = {}

        def solve(args):
            received["val"] = args.val
            return {"got": args.val}

        extra = [
            (["--val"], {"type": int, "default": 7}),
        ]
        stdout, _ = self._run_calc_main(
            solve, extra_args=extra, cli_args=["--val", "99"])
        assert received["val"] == 99
        result = json.loads(stdout.strip())
        assert result["got"] == 99

    def test_suppresses_print_from_solver(self):
        """print() inside solver is captured, not leaked to real stdout."""
        def solve(args):
            print("THIS SHOULD BE SUPPRESSED")
            return {"clean": True}

        stdout, _ = self._run_calc_main(solve)
        assert "SUPPRESSED" not in stdout
        result = json.loads(stdout.strip())
        assert result["clean"] is True
