"""Native corpus/CLI/PPTX evidence. Run only in the isolated Eqnedit64 CI lane."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import types
import xml.etree.ElementTree as ET

M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
A = "http://schemas.openxmlformats.org/drawingml/2006/main"
OLD = "http://schemas.microsoft.com/office/2004/12/omml"
COMPONENT = Path(__file__).resolve().parents[1]
REPO = COMPONENT.parents[1]


def normalized(node):
    name = node.tag.split("}")[-1]
    if name in {"ctrlPr", "rPr", "baseJc"}:
        return None
    if name == "span":
        return None
    if name == "r":
        sty = next((n.attrib.get(f"{{{M}}}val", n.attrib.get(f"{{{OLD}}}val"))
                    for n in node.iter() if n.tag.split("}")[-1] == "sty"), "i")
        span = next((n for n in node if n.tag == "span"), None)
        if span is not None:
            colors = re.findall(r"(?:^|;)color:(#[0-9a-fA-F]{6}|black)", span.attrib.get("style", ""))
            color = colors[-1].lstrip("#").upper() if colors else "000000"
            if color == "BLACK": color = "000000"
            text = "".join(span.itertext())
        else:
            color = next((n.attrib["val"] for n in node.iter(f"{{{A}}}srgbClr")), "000000")
            text = "".join(n.text or "" for n in node.iter(f"{{{A}}}t"))
        return ["run", sty, color, text]
    children = []
    for child in node:
        item = normalized(child)
        if item is None: continue
        if children and item[0] == "run" and children[-1][:3] == item[:3]:
            children[-1][3] += item[3]
        else: children.append(item)
    attrs = sorted((k.split("}")[-1], v) for k, v in node.attrib.items())
    return [name, attrs, children]


def semantic_tree(root):
    """Compare delimiter content independent of its stretchable layout form."""
    root = ET.fromstring(ET.tostring(root))
    for parent in list(root.iter())[::-1]:
        for node in list(parent):
            if node.tag != f"{{{M}}}d":
                continue
            props = node.find(f"{{{M}}}dPr")
            tokens = [props.find(f"{{{M}}}{tag}").get(f"{{{M}}}val") for tag in ("begChr", "endChr")]
            def run(text):
                value = ET.Element(f"{{{M}}}r")
                sty = ET.SubElement(ET.SubElement(value, f"{{{M}}}rPr"), f"{{{M}}}sty")
                sty.set(f"{{{M}}}val", "p")
                ET.SubElement(value, f"{{{A}}}t").text = text
                return value
            index = list(parent).index(node)
            parent.remove(node)
            for offset, child in enumerate([run(tokens[0]), *list(node.find(f"{{{M}}}e")), run(tokens[1])]):
                parent.insert(index + offset, child)
    return root


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--app", type=Path, required=True)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if os.environ.get("EQNEDIT64_ISOLATED_TEST_SESSION") != "1" or os.environ.get("GITHUB_ACTIONS") != "true":
        raise RuntimeError("Native corpus export requires the isolated hosted CI lane")
    app, probe = args.app.resolve(), args.probe.resolve()
    if not app.is_file() or not probe.is_file():
        raise FileNotFoundError("Required native fixture binaries are missing")
    args.output.mkdir(parents=True, exist_ok=True)
    corpus = json.loads((COMPONENT / "tests/office_bit_corpus.json").read_text("utf-8"))
    report = {"schema": "eqnedit64.powerpoint-math.v1", "covered": 0, "rejected": 0,
              "failures": [], "powerpoint": "NOT RUN: separate approved desktop acceptance required",
              "release_acceptance": False,
              "schema_validation": "pending independent SDK gate",
              "semantic_comparison_scope": "structure/operators/text/style/colour; run grouping and dialect-specific control/vertical layout properties normalized",
              "app_sha256": hashlib.sha256(app.read_bytes()).hexdigest(),
              "source_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()}

    def call(binary, arguments):
        return subprocess.run([str(binary), *map(str, arguments)], capture_output=True,
                              text=True, encoding="utf-8", errors="replace", timeout=30)
    try:
        excluded = {f["id"] for f in corpus if f.get("coverage") == "excluded"}
        assert excluded == {"boxed", "binom", "cancel"} and len(corpus) == 64
        for fixture in corpus:
            source = args.output / (fixture["id"] + ".tex")
            output = args.output / (fixture["id"] + ".pptxmath")
            source.write_text(fixture["tex"], "utf-8")
            native = call(app, [source, output])
            if fixture["id"] in excluded:
                assert native.returncode != 0 and fixture["tex"] in native.stderr and not output.exists(), fixture["id"]
                report["rejected"] += 1
                continue
            assert native.returncode == 0, f'{fixture["id"]}: {native.stderr}'
            root = ET.fromstring(output.read_bytes())
            assert root.tag == f"{{{M}}}oMathPara"
            assert all(n.tag.startswith((f"{{{M}}}", f"{{{A}}}")) for n in root.iter())
            html = call(probe, [source])
            assert html.returncode == 0, html.stderr
            primary = html.stdout.split("<!--[if gte msEquation 12]>", 1)[1].split("<![endif]-->", 1)[0]
            semantic = semantic_tree(root)
            assert normalized(semantic_tree(ET.fromstring(primary))) == normalized(semantic), fixture["id"]
            if r"\left" in fixture["tex"] and fixture["id"] not in {"one-sided", "brace-cases"}:
                assert root.find(f".//{{{M}}}d/{{{M}}}dPr/{{{M}}}grow").get(f"{{{M}}}val") == "1", fixture["id"]
            for run in root.iter(f"{{{M}}}r"):
                text = "".join(run.itertext())
                cjk = any(0x3000 <= ord(c) <= 0x9FFF or 0xF900 <= ord(c) <= 0xFAFF or
                          0xFF00 <= ord(c) <= 0xFFEF or 0x20000 <= ord(c) <= 0x3134F for c in text)
                assert run.find(f"{{{A}}}rPr").get("lang") == ("ja-JP" if cjk else "en-US"), fixture["id"]
            for expected in fixture.get("run_styles", []):
                matching = "".join("".join(n.itertext()) for n in root.iter(f"{{{M}}}r")
                    if ((n.find(f"{{{M}}}rPr/{{{M}}}sty").get(f"{{{M}}}val")
                         if n.find(f"{{{M}}}rPr/{{{M}}}sty") is not None else None) == expected["sty"]))
                assert expected["text"] in matching, (fixture["id"], expected)
            if fixture.get("function_name"):
                names = []
                for node in root.iter(f"{{{M}}}fName"):
                    base = node.find(f"{{{M}}}limLow/{{{M}}}e")
                    names.append("".join((base if base is not None else node).itertext()))
                assert fixture["function_name"] in names, fixture["id"]
            if fixture.get("nary_tail_outside"):
                tail = fixture["nary_tail_outside"]
                parents = {child: parent for parent in semantic.iter() for child in parent}
                for nary in semantic.iter(f"{{{M}}}nary"):
                    siblings = list(parents[nary])
                    following = siblings[siblings.index(nary) + 1:]
                    assert any("".join(node.itertext()).endswith(tail) for node in following), fixture["id"]
            report["covered"] += 1
        assert report["covered"] == 61 and report["rejected"] == 3
        for i, tex in enumerate([r"\unknown{x}", r"\begin{unknown}x\end{unknown}",
                                 r"\frac{x}{y", r"\begin{aligned}x\end{matrix}"]):
            source = args.output / f"invalid-{i}.tex"
            source.write_text(tex, "utf-8")
            target = args.output / "preserved.pptxmath"
            target.write_bytes(b"previous output")
            result = call(app, [source, target])
            assert result.returncode != 0 and tex in result.stderr
            assert target.read_bytes() == b"previous output"
        api = load_module("eqnedit64_checked_api", REPO / "packages/eqnedit64/src/eqnedit64/api.py")
        api.backend_path = lambda: app
        sys.modules["eqnedit64"] = types.SimpleNamespace(__version__="3.1.4", export_powerpoint_math=api.export_powerpoint_math)
        assembly = load_module("native_markdown_equations", REPO / "packages/radia-mcp/src/radia_mcp/presentation/_markdown_equations.py")
        markdown = "# Native equations\n\nInline $\\mathbf{A}x$ and $\\sin(x+1)$.\n\n$$\\int_a^b x^2\\,dx$$\n\n> Explain $\\alpha$ and\n> $$\\frac{x}{y}$$\n"
        markdown += "\n# Structures\n\n$$\\sum_{k=1}^{n} k^2$$\n\n$$\\begin{matrix}a&b\\\\c&d\\end{matrix}$$\n\n$$\\left(\\frac{x}{y}\\right)$$\n"
        deck = args.output / "native-equations.pptx"
        report["deck"] = assembly.presentation_markdown_to_pptx(markdown, str(deck))
        # Do not store private runtime paths in public evidence.
        report["deck"].pop("output_path")
        report["deck_sha256"] = hashlib.sha256(deck.read_bytes()).hexdigest()
        report["passed"] = True
    except Exception as exc:
        report["failures"].append(f"{type(exc).__name__}: {exc}")
        report["passed"] = False
        raise
    finally:
        (args.output / "result.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), "utf-8")
    print(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    main()
