"""Reproduce documentation previews without running a field solver.

The coil image is an unchanged saved notebook output. The thermal curve is
the explicitly labelled analytical reference, not newly computed FEM data.
Run from any directory: python docs/application_cases/render_previews.py
"""
from pathlib import Path
import base64
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    here = Path(__file__).absolute().parent
    notebook = here.parent / "stream_function/theory.ipynb"
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
    (here / "coil_iteration.png").write_bytes(base64.b64decode(encoded))

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
    fig.savefig(here / "thermal_reference.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
