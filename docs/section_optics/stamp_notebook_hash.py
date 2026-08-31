"""Stamp the executed notebook's hash into its result sidecar.

The notebook and its JSON are a synchronized pair, and the JSON records
``notebook_sha256`` so a reader can tell whether the numbers in front of
them came from the committed ``.ipynb``.  A notebook cannot hash itself
while it runs -- executing it is what changes the file -- so the stamp is
necessarily a step afterwards.  Keeping it as a script rather than a
manual edit means re-running the notebook cannot silently leave the pair
unsynchronized:

    jupyter nbconvert --to notebook --execute --inplace section_optics_design.ipynb
    python stamp_notebook_hash.py
"""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_NOTEBOOK = HERE / "section_optics_design.ipynb"
DEFAULT_SIDECAR = HERE / "section_optics_design_results.json"


def stamp(notebook=DEFAULT_NOTEBOOK, sidecar=DEFAULT_SIDECAR):
    notebook, sidecar = Path(notebook), Path(sidecar)
    if not notebook.is_file():
        raise SystemExit(f"no notebook at {notebook}")
    if not sidecar.is_file():
        raise SystemExit(
            f"no sidecar at {sidecar}: execute the notebook first, its last "
            "cell writes the results")
    payload = json.loads(sidecar.read_text(encoding="utf-8"))
    digest = hashlib.sha256(notebook.read_bytes()).hexdigest()
    payload["notebook"] = notebook.name
    payload["notebook_sha256"] = digest
    payload["notebook_bytes"] = notebook.stat().st_size
    ordered = {key: payload[key] for key in payload if key != "results"}
    if "results" in payload:
        ordered["results"] = payload["results"]
    sidecar.write_text(json.dumps(ordered, indent=2) + "\n", encoding="utf-8")
    return digest


if __name__ == "__main__":
    arguments = sys.argv[1:]
    print(f"notebook_sha256 {stamp(*arguments) if arguments else stamp()}")
