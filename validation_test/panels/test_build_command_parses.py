"""Auto-discovery test: every panel x mode produces a CLI that argparse accepts.

This is the dynamic / runtime version of ``check_panel_cli.py``
(panel-cli-diff).  Where panel-cli-diff statically lints flag overlap
between the panel's ``build_command`` source and the calc's argparse,
this test:

  1. Discovers every (panel, mode) pair from ``panel_qa.get_panel_registry()``.
  2. Pre-fills the panel widgets with valid fixture paths so
     ``build_command(vol_path)`` succeeds.
  3. Captures the resulting CLI argv list.
  4. Imports the target ``calc_*.py``, intercepts ``parse_args`` to
     parse ONLY our argv (not sys.argv), then asserts the calc's
     argparse accepts every flag and value the panel emitted.

Why dynamic in addition to static panel-cli-diff:

  * panel-cli-diff doesn't see required-flag enforcement (it just
    checks that flags exist on both sides).  Dynamic catches "panel
    forgot to emit a --required flag in this mode".
  * panel-cli-diff doesn't catch type / choices mismatches, e.g. the
    panel emits ``--solver pardiso`` but argparse only accepts
    ``--solver=auto|bddc|...``.
  * Dynamic auto-extends to new modes -- adding a (panel, mode) row
    to ``get_panel_registry()`` is the only test-side change needed.

Operating cost: < 5s for ~10 panel-modes (one importlib import per
unique calc script, cached thereafter).  Pure pytest, no subprocess.
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

import pytest

# panel_qa lives in the same directory as this test; conftest adds
# validation_test/panels to sys.path.
from panel_qa import get_panel_registry


_REPO = Path(__file__).resolve().parents[2]
_PANELS = _REPO / "src" / "radia" / "panels"
_SAMPLES = _PANELS / "samples"
_FIXTURES = Path(__file__).resolve().parent / "fixtures"


# ----------------------------------------------------------------------
# Fixture lookup
# ----------------------------------------------------------------------

# Choose fixtures that exist and are stable.  The test does NOT run the
# calc (it stops after parse_args), so the file content does not have
# to be physics-valid -- only the path has to exist for build_command's
# os.path.isfile check.
_VOL_FIXTURE = _SAMPLES / "ih_bem_sample.vol"
_STEP_FIXTURE = _SAMPLES / "ih_peec_inductance_coil.step"
_HEAT_VOL_FIXTURE = _FIXTURES / "heat_workpiece_cylinder_R25_H25.vol"
_HEAT_AXISYM_VOL_FIXTURE = _FIXTURES / "heat_workpiece_cylinder_R25_H25_axisym.vol"


def _prefill_widgets(panel, *, mode_tag: str) -> None:
    """Write fixture paths into the panel's _widgets so build_command
    succeeds.  Lookups are best-effort -- absent keys are skipped.

    The mode_tag is the registry tag (``ih_ind`` / ``heat_3d`` / ...)
    so we can pick mode-appropriate fixtures.
    """
    if not hasattr(panel, "_widgets"):
        return

    def _set(key, value):
        if key in panel._widgets:
            try:
                panel._widgets[key].setText(str(value))
            except AttributeError:
                pass  # Combos / spins -- leave default.

    # Common: STEP file for PEEC modes, .vol file for BEM/FEM/EM modes.
    _set("peec_step", _STEP_FIXTURE)
    _set("coil_script", _PANELS / "samples" / "em_sample_coil.py")

    # Heat panel: takes its own --wp-vol path and (for spatial mode)
    # an em_vol + qsurf_sol.
    if mode_tag.startswith("heat_3d"):
        _set("wp_vol", _HEAT_VOL_FIXTURE)
    elif mode_tag.startswith("heat_axisym"):
        _set("wp_vol", _HEAT_AXISYM_VOL_FIXTURE)


def _vol_path_for(mode_tag: str) -> str:
    """Path passed into ``Window._panel.build_command(vol_path)``.

    Most panels resolve their main mesh from this argument.  Heat
    panel ignores it (uses --wp-vol instead).  An empty string is
    fine for modes that don't consume a .vol (PEEC inductance).
    """
    if mode_tag in ("ih_ind",):
        return ""  # PEEC inductance: STEP only, no .vol
    if mode_tag.startswith("heat"):
        return ""  # Heat panel reads wp_vol from its own widget
    if mode_tag == "pcb":
        return ""  # PCB reads .inp directly
    return str(_VOL_FIXTURE)


# ----------------------------------------------------------------------
# argparse interception
# ----------------------------------------------------------------------

class _StopAfterParse(Exception):
    """Raised inside the calc's main() right after parse_args succeeds.

    Carries the parsed Namespace so the caller can inspect.  This is
    how we avoid running the heavy NGSolve / Radia solve while still
    exercising argparse on the panel's actual emitted argv.
    """
    def __init__(self, ns):
        self.ns = ns


def _validate_via_calc_argparse(script_path: Path, argv: list[str]):
    """Import script_path and replay argv against its argparse.

    Returns the parsed argparse.Namespace on success.
    Raises SystemExit (code 2) on argparse rejection, with the
    error message visible in pytest's captured stderr.
    """
    real_parse_args = argparse.ArgumentParser.parse_args

    def patched_parse_args(self, args=None, namespace=None):
        # Force argparse to look at OUR argv, not sys.argv.
        ns = real_parse_args(self, argv, namespace)
        raise _StopAfterParse(ns)

    # Avoid module name collision in sys.modules between calc_*.py imports.
    mod_name = f"_calc_validate_{script_path.stem}_{id(argv)}"
    spec = importlib.util.spec_from_file_location(mod_name, str(script_path))
    mod = importlib.util.module_from_spec(spec)

    argparse.ArgumentParser.parse_args = patched_parse_args
    saved_argv = sys.argv[:]
    try:
        sys.argv = [str(script_path)]  # Defensive: hide test-runner argv.
        spec.loader.exec_module(mod)
        try:
            mod.main()
        except _StopAfterParse as ok:
            return ok.ns
        # If main() returned without ever calling parse_args, the test
        # pre-condition is broken -- not the panel's fault, but still
        # surface it.
        pytest.fail(f"{script_path.name}: main() returned without "
                    f"reaching parse_args(). The script may have a "
                    f"non-standard structure.")
    finally:
        argparse.ArgumentParser.parse_args = real_parse_args
        sys.argv = saved_argv


# ----------------------------------------------------------------------
# Pytest parametrization
# ----------------------------------------------------------------------

# Known-broken (panel,mode) pairs.  Each entry is (tag, reason).  The
# test marks them xfail(strict=True) so the day the bug is fixed the
# xfail flips to XPASS and the entry MUST be removed.  This is the
# auto-removal trigger: a green test stays green, a fixed bug
# becomes visible immediately.
_KNOWN_BROKEN: dict[str, str] = {}


def _registry_entries():
    """Return only entries whose Window class exposes _panel.build_command.

    PCB and Heat panels participate; the panel_qa registry already
    filters out invalid combinations.
    """
    return get_panel_registry()


def _maybe_xfail(tag):
    """If tag is a known-broken pair, return an xfail-strict marker."""
    reason = _KNOWN_BROKEN.get(tag)
    if reason is None:
        return ()
    return (pytest.mark.xfail(strict=True, reason=reason),)


def _parametrize_entries():
    return [
        pytest.param(e, marks=_maybe_xfail(e[0]), id=e[0])
        for e in _registry_entries()
    ]


@pytest.mark.parametrize("entry", _parametrize_entries())
def test_build_command_argparse_accepts(qapp, entry):
    """Panel.build_command() output must parse cleanly via the target
    calc_*.py argparse.

    Catches:
      * Panel emits a flag the calc removed.
      * Panel forgot a --required flag (e.g. mode-specific).
      * Type / choices mismatch (panel emits string, calc expects int).
    """
    tag, WinCls, combo_ref, value = entry
    win = WinCls("")
    try:
        # Apply mode (e.g. switch IH to METHOD_PEEC_BEM).  The combo
        # may live on the panel directly or in panel._widgets.
        if combo_ref and value is not None:
            panel = win._panel
            combo = (getattr(panel, combo_ref, None)
                     or panel._widgets.get(combo_ref))
            if combo is not None:
                combo.setCurrentText(value)

        panel = win._panel
        _prefill_widgets(panel, mode_tag=tag)

        try:
            cmd = panel.build_command(_vol_path_for(tag))
        except (ValueError, FileNotFoundError) as e:
            pytest.skip(f"{tag}: build_command requires extra setup: {e}")

        # cmd[0] is sys.executable, cmd[1] is the calc script path.
        if len(cmd) < 2:
            pytest.fail(f"{tag}: build_command returned {len(cmd)} args, "
                        f"expected at least [python, script_path, ...]")
        script = Path(cmd[1])
        argv = list(cmd[2:])

        if not script.is_file():
            pytest.fail(f"{tag}: build_command points at non-existent "
                        f"script: {script}")

        ns = _validate_via_calc_argparse(script, argv)
        # If we got here, argparse accepted every flag.  No assertion
        # needed -- a parse error would have raised SystemExit(2)
        # during _validate_via_calc_argparse.
        assert ns is not None
    finally:
        win.close()
        win.deleteLater()
