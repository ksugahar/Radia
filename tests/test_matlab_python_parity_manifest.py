from pathlib import Path
import importlib.util


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "check_matlab_python_parity", ROOT / "tools" / "check_matlab_python_parity.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_every_python_module_has_a_checked_matlab_classification():
    result = MODULE.audit(ROOT)
    assert result["ok"], result["errors"]
    assert result["classified_file_count"] == result["python_file_count"]
    assert result["counts"]["native-mex"] >= 4
    assert result["counts"]["python-fallback"] >= 80
    assert "acoustic-python" in result["python_fallback_families"]
    axifem = next(
        item for item in result["binary_extensions"] if item["python"] == "axifem.pyd"
    )
    assert axifem["classification"] == "native-mex"
    assert axifem["status"] == "focused-native-commands"
    assert axifem["native_commands"] == [
        "axifem.q1_magnetic_element_matrices",
        "axifem.q2_magnetic_element_matrices",
    ]
    backlog = result["native_promotion_backlog"]
    assert [item["priority"] for item in backlog] == list(range(1, 8))
    assert [item["family"] for item in backlog] == [
        "axifem",
        "vim-esim-ih",
        "bem-peec-sibc",
        "kelvin-dtn",
        "acoustic-cq-bem-fsi",
        "motor-maglev",
        "coil-cad",
    ]
