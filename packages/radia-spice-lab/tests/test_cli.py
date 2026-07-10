"""CLI smoke tests via subprocess.

These intentionally invoke the installed ``ltspice-convert`` console
script (not the function directly) so we exercise the same code path
end-users will hit, including exit-code propagation through the
setuptools-generated entry-point wrapper.
"""
from __future__ import annotations

import json
import math
import re
import subprocess
import sys
from pathlib import Path

import pytest


FIXTURES = Path(__file__).parent / "fixtures" / "bidirectional"
RC_ASC = FIXTURES / "00_converter_test_rc_lowpass.asc"
RC_CIR = FIXTURES / "00_converter_test_rc_lowpass.cir"


def run_cli(*args: str, check: bool = False) -> subprocess.CompletedProcess:
    """Invoke the CLI as `python -m ltspice_converter.cli ARGS`.

    Using `python -m` avoids depending on whether the console script
    is on PATH in the CI environment.
    """
    return subprocess.run(
        [sys.executable, "-m", "ltspice_converter.cli", *args],
        capture_output=True,
        text=True,
        check=check,
    )


# =============================================================================
# Basics
# =============================================================================

def test_version():
    r = run_cli("--version")
    assert r.returncode == 0
    assert r.stdout.strip().startswith("spice-circuit-lab ")


def test_help():
    r = run_cli("--help")
    assert r.returncode == 0
    assert "input file" in r.stdout.lower() or "input" in r.stdout.lower()


def test_error_on_unknown_extension(tmp_path):
    bad = tmp_path / "x.txt"
    bad.write_text("garbage")
    r = run_cli(str(bad))
    assert r.returncode != 0
    assert "unknown" in r.stderr.lower() or "extension" in r.stderr.lower()


def test_error_on_missing_file(tmp_path):
    r = run_cli(str(tmp_path / "nope.asc"))
    assert r.returncode != 0


# =============================================================================
# Convert mode
# =============================================================================

def test_convert_asc_to_cir(tmp_path):
    out = tmp_path / "rc.cir"
    r = run_cli(str(RC_ASC), "-o", str(out))
    assert r.returncode == 0, r.stderr
    assert out.is_file()
    text = out.read_text(encoding="utf-8")
    assert "R1" in text
    assert "C1" in text


def test_convert_cir_to_asc(tmp_path):
    out = tmp_path / "rc.asc"
    r = run_cli(str(RC_CIR), "-o", str(out))
    assert r.returncode == 0, r.stderr
    assert out.is_file()
    text = out.read_text(encoding="latin-1")
    assert "Version" in text or "SYMBOL" in text


def test_convert_auto_output_path(tmp_path):
    src = tmp_path / "rc.asc"
    src.write_text(RC_ASC.read_text(encoding="utf-8"))
    r = run_cli(str(src))   # no -o, expect rc.cir alongside
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "rc.cir").is_file()


def test_convert_with_to_flag(tmp_path):
    src = tmp_path / "rc.asc"
    src.write_text(RC_ASC.read_text(encoding="utf-8"))
    r = run_cli(str(src), "--to", "py")
    assert r.returncode == 0, r.stderr
    py_out = tmp_path / "rc.py"
    assert py_out.is_file()
    assert "schemdraw" in py_out.read_text(encoding="utf-8")


def test_convert_batch_to_dir(tmp_path):
    out_dir = tmp_path / "batch"
    inputs = list(FIXTURES.glob("00_converter_test_*.asc"))
    assert len(inputs) >= 5
    r = run_cli(*[str(p) for p in inputs], "-o", str(out_dir), "--to", "cir")
    assert r.returncode == 0, r.stderr
    cir_outputs = list(out_dir.glob("*.cir"))
    assert len(cir_outputs) == len(inputs)


# =============================================================================
# Check (lint) mode
# =============================================================================

def test_check_clean_file_pass():
    r = run_cli("--check", str(RC_ASC))
    assert r.returncode == 0, r.stderr
    assert "PASS" in r.stdout
    assert "[warn]" not in r.stdout


def test_check_strict_clean_still_pass():
    r = run_cli("--check", "--strict", str(RC_ASC))
    assert r.returncode == 0, r.stderr


def test_check_multiple_files():
    inputs = list(FIXTURES.glob("00_converter_test_*.asc"))[:3]
    r = run_cli("--check", *[str(p) for p in inputs])
    assert r.returncode == 0, r.stderr
    # Each file should produce its own section
    assert r.stdout.count("==") >= len(inputs) * 2  # header has == on both sides


# =============================================================================
# Info mode
# =============================================================================

def test_info_text():
    r = run_cli("--info", str(RC_ASC))
    assert r.returncode == 0, r.stderr
    assert "component_count" in r.stdout
    assert "component_types" in r.stdout


def test_info_json():
    r = run_cli("--info", "--json", str(RC_ASC))
    assert r.returncode == 0, r.stderr
    data = json.loads(r.stdout)
    assert isinstance(data, list)
    assert len(data) == 1
    item = data[0]
    assert item["format"] == "asc"
    assert item["component_count"] == 3
    assert item["component_types"] == {"I": 1, "R": 1, "C": 1}


# =============================================================================
# .asy search dir
# =============================================================================

def test_asy_dir_flag_accepted(tmp_path):
    # Flag should be parsed even if dir does not exist (just no-op there)
    r = run_cli("--asy-dir", str(tmp_path),
                "--info", str(RC_ASC))
    assert r.returncode == 0, r.stderr


# =============================================================================
# C5: static lint checks
# =============================================================================

def _write_cir(tmp_path: Path, name: str, body: str) -> Path:
    p = tmp_path / name
    p.write_text(body, encoding="utf-8")
    return p


def test_check_detects_duplicate_instance(tmp_path):
    cir = _write_cir(tmp_path, "dup.cir",
                     "* dup\nR1 a b 1k\nR1 b c 2k\n.end\n")
    r = run_cli("--check", str(cir))
    assert "duplicate instance name" in r.stdout
    assert "R1" in r.stdout


def test_check_detects_floating_node(tmp_path):
    cir = _write_cir(tmp_path, "float.cir",
                     "* float\nV1 in 0 5\nR1 in mid 1k\n"
                     "R2 mid out 1k\nR3 dangle 0 1k\n.end\n")
    r = run_cli("--check", str(cir))
    assert "floating node" in r.stdout
    assert "dangle" in r.stdout or "out" in r.stdout


def test_check_detects_orphan_model(tmp_path):
    cir = _write_cir(tmp_path, "orphan.cir",
                     "* orphan\nV1 in 0 5\nR1 in out 1k\n"
                     ".model UNUSED D\n.end\n")
    r = run_cli("--check", str(cir))
    assert "never used" in r.stdout
    assert "UNUSED" in r.stdout


def test_check_detects_undefined_model(tmp_path):
    cir = _write_cir(tmp_path, "undef_model.cir",
                     "* undef\nV1 in 0 5\nD1 in 0 MY_CUSTOM_DIODE\n.end\n")
    r = run_cli("--check", str(cir))
    # Should mention the undefined model since it's not a standard one
    assert "MY_CUSTOM_DIODE" in r.stdout


def test_check_detects_undefined_param(tmp_path):
    cir = _write_cir(tmp_path, "undef_param.cir",
                     "* undef\nR1 in out {missingparam}\n"
                     "R2 in out {Rval}\n.param Rval=1k\n.end\n")
    r = run_cli("--check", str(cir))
    assert "missingparam" in r.stdout
    # The defined param should NOT be flagged
    assert "Rval" not in r.stdout.split("parameter(s)")[1].split("\n")[0]


def test_check_clean_netlist_no_warnings(tmp_path):
    cir = _write_cir(tmp_path, "clean.cir",
                     "* clean\nV1 in 0 5\nR1 in out 1k\nC1 out 0 1u\n"
                     ".ac dec 100 1 100k\n.end\n")
    r = run_cli("--check", str(cir))
    assert r.returncode == 0
    assert "PASS" in r.stdout
    # Should have no [warn] lines from C5 checks
    static_warns = [
        "duplicate", "floating", "never used", "referenced but not defined",
        "parameter(s) referenced",
    ]
    for w in static_warns:
        assert w not in r.stdout, f"unexpected warning containing {w!r}: {r.stdout}"


# =============================================================================
# B4: unparsed-line surfacing + suggestions
# =============================================================================

def test_check_surfaces_unparsed_line(tmp_path):
    cir = _write_cir(tmp_path, "unparsed.cir",
                     "* test\nV1 in 0 5\nZorg foo bar baz\n.end\n")
    r = run_cli("--check", str(cir))
    assert "line 3" in r.stdout
    assert "Zorg" in r.stdout
    assert "unrecognised element" in r.stdout.lower() or \
           "SPICE elements start" in r.stdout


def test_check_typo_suggestion(tmp_path):
    # 'Resistor' starts with R and gets accepted by the parser (oddly named),
    # but a token like 'Diode' (starts with D, valid prefix) also passes.
    # Use an actually-unrecognised token so we hit the hint path.
    cir = _write_cir(tmp_path, "typo.cir",
                     "* typo test\nR1 in out 1k\nZorg foo bar\n.end\n")
    r = run_cli("--check", str(cir))
    assert "Zorg" in r.stdout


# =============================================================================
# C5 false-positive fixes (regression tests)
# =============================================================================

def test_check_x_subckt_invocation_not_treated_as_model_ref(tmp_path):
    """C5-fp-1: X<name> ... <subckt_name> must look up `.subckt`, not `.model`.

    Before v0.3.7 the lint pooled X-class refs into `referenced_models`
    along with D/Q/M/J, so a perfectly valid X1 invocation produced a
    spurious 'model INV referenced but not defined inline' warning.
    """
    cir = _write_cir(
        tmp_path, "xsubckt.cir",
        "* X invocation must NOT trigger a model-undefined warning\n"
        ".SUBCKT INV 1 2 3\n"
        "R1 1 2 1k\n"
        "R2 2 3 1k\n"
        ".ENDS INV\n"
        "V1 3 0 5\n"
        "X1 OUT 2 3 INV\n"
        "RL OUT 0 100k\n"
        ".op\n"
        ".end\n"
    )
    r = run_cli("--check", str(cir))
    # No warning about INV being an undefined model.
    assert "INV" not in r.stdout.split("model")[-1].split("\n")[0] \
        if "model(s) referenced" in r.stdout else True, \
        f"X-subckt INV was wrongly flagged: {r.stdout}"
    # And we should not see the bare 'INV referenced' warning either.
    assert "INV" not in [
        w.strip()
        for line in r.stdout.split("\n")
        for w in (line.split(":", 1)[-1].split(",") if "referenced" in line else [])
    ], f"INV mentioned as undefined: {r.stdout}"


def test_check_model_inside_subckt_not_flagged_as_orphan(tmp_path):
    """C5-fp-2: .model defined inside a .subckt body has local scope.

    Before v0.3.7 the orphan-model check scanned the WHOLE netlist for
    .model declarations but only the top-level components for refs,
    so subckt-internal .model lines were flagged as 'never used'.
    """
    cir = _write_cir(
        tmp_path, "subckt_model.cir",
        "* .model inside .subckt body must NOT be flagged as orphan\n"
        ".SUBCKT INV 1 2 3\n"
        "M1 1 2 0 0 NCH1 W=4u L=6u\n"
        "M2 1 2 3 3 PCH1 W=4u L=6u\n"
        ".model NCH1 NMOS Level=1 Kp=20u Vto=0.7\n"
        ".model PCH1 PMOS Level=1 Kp=10u Vto=-0.7\n"
        ".ENDS INV\n"
        "V1 3 0 5\n"
        "X1 OUT 2 3 INV\n"
        "RL OUT 0 100k\n"
        ".op\n"
        ".end\n"
    )
    r = run_cli("--check", str(cir))
    # No "never used" warning for NCH1/PCH1 (they are used inside the subckt body).
    assert "never used" not in r.stdout, \
        f"subckt-internal .model wrongly flagged: {r.stdout}"


def test_check_undefined_subckt_warning_when_truly_missing(tmp_path):
    """Positive test for the new subckt-undefined check.

    When an X invocation names a subckt that has no .subckt definition
    AND does not look like a standard library part, we should warn.
    """
    cir = _write_cir(
        tmp_path, "missing_subckt.cir",
        "* X invokes a subckt that is not defined inline\n"
        "V1 in 0 5\n"
        "X1 in out MY_CUSTOM_BLOCK\n"
        ".end\n"
    )
    r = run_cli("--check", str(cir))
    assert "MY_CUSTOM_BLOCK" in r.stdout, \
        f"expected missing subckt to be reported, got: {r.stdout}"


# =============================================================================
# E1: text-based check/info entry points (also used by the MCP server)
# =============================================================================


def test_check_text_clean_netlist():
    """check_text() returns no warnings for a well-formed netlist."""
    from ltspice_converter.cli import check_text
    netlist = (
        "* RC Lowpass\n"
        "V1 in 0 AC 1\n"
        "R1 in out 1k\n"
        "C1 out 0 1u\n"
        ".ac dec 20 1 100k\n"
        ".end\n"
    )
    info, warn = check_text(netlist, "cir")
    assert warn == [], f"unexpected warnings on clean netlist: {warn}"
    assert any("3 components" in m for m in info), info


def test_check_text_flags_duplicate_instance():
    """check_text() surfaces a duplicate-instance C5 lint warning."""
    from ltspice_converter.cli import check_text
    netlist = (
        "V1 in 0 AC 1\n"
        "R1 in out 1k\n"
        "R1 out 0 2k\n"   # duplicate
        ".end\n"
    )
    _info, warn = check_text(netlist, "cir")
    assert any("duplicate" in w.lower() for w in warn), warn


def test_check_text_voltage_controlled_switch_model_and_nodes():
    """Voltage-controlled S switches use four node pins plus a .model."""
    from ltspice_converter.cli import check_text
    netlist = (
        "* PWM buck switch seed\n"
        "VIN vin 0 DC 24\n"
        "VGATE gate 0 PULSE(0 10 0 20n 20n 2.3u 10u)\n"
        "SMAIN vin sw gate 0 SW_MAIN\n"
        "L1 sw out 120u\n"
        "RLOAD out 0 5\n"
        ".model SW_MAIN SW(Ron=50m Roff=10Meg Vt=4 Vh=0.5)\n"
        ".tran 0 1m\n"
        ".end\n"
    )
    _info, warn = check_text(netlist, "cir")
    assert not any("floating node" in w.lower() for w in warn), warn
    assert not any("sw_main" in w.lower() for w in warn), warn


def test_check_text_controlled_source_counts_control_pins():
    from ltspice_converter.cli import check_text
    netlist = (
        "* VCVS control pins are node pins\n"
        "VIN in 0 AC 1\n"
        "EGAIN out 0 in 0 10\n"
        "RLOAD out 0 10k\n"
        ".end\n"
    )
    _info, warn = check_text(netlist, "cir")
    assert not any("floating node" in w.lower() and "in" in w.lower() for w in warn), warn


def test_check_text_x_inline_params_use_last_bare_token_as_subckt():
    from ltspice_converter.cli import check_text
    netlist = (
        "* X inline params\n"
        "V1 in 0 DC 5\n"
        "X1 in out gain_block Avol=10\n"
        "RLOAD out 0 10k\n"
        ".subckt gain_block a b params: Avol=1\n"
        "E1 b 0 a 0 {Avol}\n"
        ".ends gain_block\n"
        ".end\n"
    )
    _info, warn = check_text(netlist, "cir")
    joined = "\n".join(warn).lower()
    assert "avol=10" not in joined, warn
    assert ".subckt declared but never used" not in joined, warn


def test_check_text_behavioral_source_v_reference_not_floating():
    from ltspice_converter.cli import check_text
    netlist = (
        "* Behavioral source references input node\n"
        "V1 in 0 SIN(0 1 1k)\n"
        "B1 out 0 V=limit(V(in)*2,-1.5,1.5)\n"
        "RLOAD out 0 1k\n"
        ".end\n"
    )
    _info, warn = check_text(netlist, "cir")
    assert not any("floating node" in w.lower() and "in" in w.lower() for w in warn), warn


def test_schemdraw_generated_script_preserves_raw_spice_lines():
    from ltspice_converter import netlist_to_schemdraw, schemdraw_to_netlist
    from ltspice_converter.topology import topology_equivalent

    netlist = (
        "* Buck switch raw-line preservation\n"
        "VIN vin 0 DC 24\n"
        "VGATE gate 0 PULSE(0 10 0 20n 20n 2.3u 10u)\n"
        "SMAIN vin sw gate 0 SW_MAIN\n"
        "L1 sw out 120u\n"
        "RLOAD out 0 5\n"
        ".model SW_MAIN SW(Ron=50m Roff=10Meg Vt=4 Vh=0.5)\n"
        ".end\n"
    )
    script = netlist_to_schemdraw(netlist, "buck_raw_preserve")
    recovered = schemdraw_to_netlist(script, "buck_raw_preserve")
    equivalent, info = topology_equivalent(netlist, recovered)
    assert equivalent, info
    assert "SMAIN vin sw gate 0 SW_MAIN" in recovered


def test_asc_roundtrip_preserves_unknown_subckt_raw_line():
    from ltspice_converter import asc_to_netlist, netlist_to_asc
    from ltspice_converter.topology import topology_equivalent

    netlist = (
        "* Unknown .asy X-line should stay simulation-correct\n"
        "V1 in 0 DC 5\n"
        "X1 in out gain_block Avol=10\n"
        "RLOAD out 0 10k\n"
        ".subckt gain_block a b params: Avol=1\n"
        "E1 b 0 a 0 {Avol}\n"
        ".ends gain_block\n"
        ".op\n"
        ".end\n"
    )
    asc = netlist_to_asc(netlist)
    assert "TEXT" in asc
    assert "!X1 in out gain_block Avol=10" in asc
    assert "SYMBOL gain_block" not in asc

    recovered = asc_to_netlist(asc, use_ltspice=False)
    equivalent, info = topology_equivalent(netlist, recovered)
    assert equivalent, info
    assert "X1 in out gain_block Avol=10" in recovered


def test_asc_layout_escapes_dense_source_lc_input_node():
    from ltspice_converter import asc_to_netlist, netlist_to_asc
    from ltspice_converter.topology import topology_equivalent

    netlist = (
        "* Dense source/L/C boost input should remain readable\n"
        "V1 in 0 DC 12\n"
        "L1 in sw 10u\n"
        "S1 sw 0 gate 0 SW\n"
        "D1 sw out D\n"
        "CIN in 0 22u\n"
        "RIN in cin_damp 1\n"
        "CDAMP cin_damp 0 2.2u\n"
        "COUT out 0 47u\n"
        "RLOAD out 0 24\n"
        "VCTRL gate 0 PULSE(0 10 0 20n 20n 2u 5u)\n"
        ".model SW SW(Ron=50m Roff=10Meg Vt=4 Vh=0.5)\n"
        ".model D D(Is=10n Rs=40m)\n"
        ".tran 0 2m\n"
        ".end\n"
    )
    asc = netlist_to_asc(netlist)
    recovered = asc_to_netlist(asc, use_ltspice=False)
    equivalent, info = topology_equivalent(netlist, recovered)
    assert equivalent, info

    symbols = {}
    current = None
    for line in asc.splitlines():
        match = re.match(r"^SYMBOL\s+\S+\s+(-?\d+)\s+(-?\d+)\s+\S+", line)
        if match:
            current = (int(match.group(1)), int(match.group(2)))
            continue
        match = re.match(r"^SYMATTR\s+InstName\s+(\S+)", line)
        if match and current is not None:
            symbols[match.group(1)] = current
            current = None

    assert {"V1", "L1", "CIN"} <= set(symbols)
    v1_l1_distance = math.hypot(
        symbols["V1"][0] - symbols["L1"][0],
        symbols["V1"][1] - symbols["L1"][1],
    )
    assert v1_l1_distance >= 64


def test_asc_conversion_preserves_local_semiconductor_model_as_raw_line():
    from ltspice_converter import asc_to_netlist, netlist_to_asc

    netlist = (
        "* Local semiconductor model should not trigger LTspice library clash\n"
        "VCC vcc 0 DC 12\n"
        "RB vcc base 100k\n"
        "RC vcc col 3.3k\n"
        "RE emit 0 680\n"
        "Q1 col base emit 2N3904\n"
        ".model 2N3904 NPN(Is=6.734f Bf=416.4 Vaf=74.03)\n"
        ".op\n"
        ".end\n"
    )
    asc = netlist_to_asc(netlist)
    assert "SYMBOL npn" not in asc
    assert "!Q1 col base emit 2N3904" in asc
    assert "!.model 2N3904 NPN(Is=6.734f Bf=416.4 Vaf=74.03)" in asc

    recovered = asc_to_netlist(asc, use_ltspice=False)
    assert "Q1 col base emit 2N3904" in recovered
    assert ".model 2N3904 NPN(Is=6.734f Bf=416.4 Vaf=74.03)" in recovered


def test_info_text_counts_by_class():
    """info_text() returns per-class component counts for a netlist."""
    from ltspice_converter.cli import info_text
    netlist = (
        "V1 in 0 AC 1\n"
        "R1 in m 1k\n"
        "R2 m out 1k\n"
        "C1 out 0 1u\n"
        ".end\n"
    )
    out = info_text(netlist, "cir")
    assert out["component_count"] == 4
    assert out["component_types"] == {"V": 1, "R": 2, "C": 1}


def test_mcp_tools_registered():
    """The mcp_server module exposes its documented public tool surface.

    Requires the optional ``[mcp]`` extra. CI runs with ``[test]`` only,
    so we skip if the ``mcp`` dependency isn't installed.
    """
    pytest.importorskip("mcp.server.fastmcp")
    from ltspice_converter import mcp_server
    names = {
        t.name if hasattr(t, "name") else getattr(t, "fn", t).__name__
        for t in mcp_server.mcp._tool_manager._tools.values()
    }
    expected = {
        "netlist_to_schemdraw", "schemdraw_to_netlist",
        "netlist_to_asc", "asc_to_netlist",
        "check_circuit", "info_circuit", "compare_topology",
        "circuit_knowledge", "buck_seed", "patentability_search_plan",
    }
    assert names == expected

    readme = (Path(__file__).parents[1] / "README.md").read_text(encoding="utf-8")
    assert "Exposes ten tools:" in readme
    for name in expected:
        assert f"`{name}" in readme
