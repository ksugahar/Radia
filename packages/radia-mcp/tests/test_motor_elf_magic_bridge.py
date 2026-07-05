# -*- coding: utf-8 -*-
"""Public-safe tests for the ELF/MAGIC bridge knowledge in radia-motor."""

import os
import sys


_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.motor.elf_magic_bridge_knowledge import get_elf_magic_bridge  # noqa: E402


def test_run_artifact_contract_keeps_mao_as_primary_result_log():
    doc = get_elf_magic_bridge("run_artifact_contract")

    assert ".mao" in doc
    assert "primary execution log" in doc
    assert "FLUM" in doc
    assert ".mag" in doc
    assert "field/result file" in doc
    assert ".mei" in doc
    assert "never treat it as the solver result" in doc


def test_bridge_all_includes_artifact_contract_without_private_paths():
    doc = get_elf_magic_bridge("all")

    assert "Source-tool run artifact contract" in doc
    for drive in ("S", "W", "C"):
        assert f"{drive}:" + "\\" not in doc
    assert "_cross" + "val" not in doc
