"""Reject retired product spellings in tracked paths and semantic text."""
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
FORBIDDEN = re.compile(r"eqn(?:edit|edt)32|eqn" + r"edt64|mt" + r"ef|math" + r"type|数式" + r"3\.0", re.I)


def notebook_text(document):
    for cell in document.get("cells", []):
        yield "".join(cell.get("source", []))
        for output in cell.get("outputs", []):
            yield "".join(output.get("text", []))
            for mime, value in output.get("data", {}).items():
                if mime.startswith("text/"):
                    yield "".join(value) if isinstance(value, list) else str(value)


def main():
    paths = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    bad = []
    for name in filter(None, paths):
        path = ROOT / name
        if FORBIDDEN.search(name):
            bad.append(name + ": path")
        if not path.is_file():
            continue
        data = path.read_bytes()
        if b"\0" in data:
            continue
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            continue
        texts = notebook_text(json.loads(text)) if path.suffix == ".ipynb" else [text]
        if any(FORBIDDEN.search(part) for part in texts):
            bad.append(name + ": text")
    if bad:
        raise SystemExit("\n".join(bad))
    print("PASS: public product names")


if __name__ == "__main__":
    main()
