import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from radia_mcp.figure import tools


def test_overlap_check_preserves_hidden_text():
    fig, ax = plt.subplots()
    try:
        hidden = ax.text(0.1, 0.1, "hidden", visible=False)
        label = ax.text(0.5, 0.5, "label")
        tools.check_text_overlap(ax, label)
        assert not hidden.get_visible()
        assert label.get_visible()
    finally:
        plt.close(fig)


def test_overlap_check_restores_visibility_on_error(monkeypatch):
    fig, ax = plt.subplots()
    try:
        label = ax.text(0.5, 0.5, "label")

        def fail(*args):
            raise RuntimeError("drawing failed")

        monkeypatch.setattr(tools, "_occupied_points", fail)
        with pytest.raises(RuntimeError, match="drawing failed"):
            tools.check_text_overlap(ax, label)
        assert label.get_visible()
    finally:
        plt.close(fig)
