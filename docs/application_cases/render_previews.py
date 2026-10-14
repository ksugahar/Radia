"""Reproduce documentation previews without running a field solver.

The six gallery previews are redrawn by executing the JSON-only result readers.
Requires numpy, matplotlib, nbformat, nbclient and an ipykernel Python 3 kernel.
The historical coil-iteration image is an unchanged saved study output. The thermal curve is
the explicitly labelled analytical reference, not newly computed FEM data.
Run from any directory: python docs/application_cases/render_previews.py
"""
from pathlib import Path
import base64
import json
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


READERS = {
    "complex_coil": "complex_coil_geometry/complex_coil.ipynb",
    "motor": "electric_machine/cogging_skew_demo.ipynb",
    "winding": "stream_function/theory.ipynb",
    "lift": "maglev/maglev_showcase.ipynb",
    "particles": "gmsh_post/em_particle_orbits.ipynb",
    "heating": "esim_spatial/esim_spatial_demo.ipynb",
}


def render_readers(root, output):
    """Execute result readers (JSON/plots only) and extract their first figure."""
    import nbformat
    from nbclient import NotebookClient
    for case, relative in READERS.items():
        path = root / "docs" / relative
        notebook = nbformat.read(path, as_version=4)
        NotebookClient(notebook, timeout=120, kernel_name="python3", record_timing=False,
                       resources={"metadata": {"path": str(path.parent)}}).execute()
        for cell in notebook.cells:
            cell.metadata.pop("execution", None)
        nbformat.validate(notebook)
        images = [o["data"]["image/png"] for c in notebook.cells
                  if c.cell_type == "code" for o in c.get("outputs", [])
                  if "image/png" in o.get("data", {})]
        if not images:
            raise RuntimeError(f"No rendered figure for {case}")
        encoded = images[0]
        if isinstance(encoded, list):
            encoded = "".join(encoded)
        (output / f"{case}_rerun.png").write_bytes(base64.b64decode(encoded))
        nbformat.write(notebook, path)
        print(f"Rendered {case} from its committed JSON record")


def main():
    here = Path(__file__).absolute().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=here)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    root = here.parents[1]
    render_readers(root, output)
    notebook = root / "validation_test/showcase/studies/stream_function/theory.ipynb"
    cells = json.loads(notebook.read_text(encoding="utf-8"))["cells"]
    candidates = [c for c in cells if c["cell_type"] == "code"
                  and 'set_xlabel("Path-A iteration")' in "".join(c["source"])]
    if len(candidates) != 1:
        raise RuntimeError("Expected one saved Path-A convergence figure")
    images = [o["data"]["image/png"] for o in candidates[0].get("outputs", [])
              if "image/png" in o.get("data", {})]
    if len(images) != 1:
        raise RuntimeError("Expected one saved PNG; do not substitute a new solve")
    encoded = images[0]
    if isinstance(encoded, list):
        encoded = "".join(encoded)
    (output / "coil_iteration.png").write_bytes(base64.b64decode(encoded))

    # Exact solution declared in axisymmetric_p2_thermal.ipynb and its helper.
    radius, k, rho, cp, t0, a = .025, 46.6, 7800., 467., 293.15, 1e4
    radial = np.linspace(0., radius, 200)
    fig, ax = plt.subplots(figsize=(7.2, 4.3), constrained_layout=True)
    for seconds, color in [(0., "#64748b"), (5., "#0072b2"), (10., "#d55e00")]:
        temperature = t0 + a * radial**2 + 4 * k * a * seconds / (rho * cp)
        ax.plot(radial * 1000, temperature, label=f"t = {seconds:g} s", color=color)
    ax.set(xlabel="Radius r (mm)", ylabel="Temperature (K)",
           title="Boundary-heated cylinder: analytical reference")
    ax.legend(); ax.grid(alpha=.22)
    fig.savefig(output / "thermal_reference.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
