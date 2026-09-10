"""Keep every node an explicit renderer decision, including state-only nodes."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def test_renderer_dispatch_is_exhaustive():
    header = (ROOT / "src/equation_node.h").read_text(encoding="utf-8")
    enum = re.search(r"enum Tag\s*:\s*uint8_t\s*\{(.*?)\};", header, re.S)[1]
    tags = set(re.findall(r"\bk[A-Z]\w*", enum))
    render = (ROOT / "src/equation_render.cpp").read_text(encoding="utf-8")
    start = render.index("Layout layout_node(")
    opening = render.index("{", start)
    depth = 1
    end = opening + 1
    while depth:
        depth += (render[end] == "{") - (render[end] == "}")
        end += 1
    dispatch = render[opening:end]
    cases = set(re.findall(r"case Node::(k\w+)\s*:", dispatch))
    assert cases == tags, f"missing={tags - cases}, extra={cases - tags}"
    assert not re.search(r"\bdefault\s*:", dispatch)
    assert "layout_fallback(" not in render


if __name__ == "__main__":
    test_renderer_dispatch_is_exhaustive()
    print("PASS: node dispatch exhaustiveness")
