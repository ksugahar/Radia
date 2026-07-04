"""Result builder for the URN -> convolution-quadrature teaching notebook."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
RADIA_MCP_SRC = REPO_ROOT / "packages" / "radia-mcp" / "src"
if str(RADIA_MCP_SRC) not in sys.path:
    sys.path.insert(0, str(RADIA_MCP_SRC))

from radia_mcp.radia_ngsolve.cq_urn import make_cq_urn_bridge_artifact  # noqa: E402


def build_artifact(output_dir: str | Path | None = None) -> dict:
    """Run the deterministic URN/CQ bridge example and write result JSON."""

    out_dir = Path(output_dir) if output_dir is not None else Path(__file__).resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)
    artifact = make_cq_urn_bridge_artifact(n_steps=100, dt=0.01, hit_index=10)
    artifact["artifact_paths"] = {
        "notebook_source_path": "docs/universal_relaxation_network/cq_urn_bridge.ipynb",
        "result_json": "docs/universal_relaxation_network/cq_urn_bridge_results.json",
        "notebook_result_sidecar": "docs/universal_relaxation_network/cq_urn_bridge_result.json",
        "figure": "embedded in docs/universal_relaxation_network/cq_urn_bridge.ipynb",
    }
    artifact["notebook_source_artifact_id"] = (
        "docs.universal_relaxation_network.cq_urn_bridge.ipynb"
    )
    artifact["notebook_source_path"] = "docs/universal_relaxation_network/cq_urn_bridge.ipynb"
    result_path = out_dir / "cq_urn_bridge_results.json"
    result_path.write_text(json.dumps(artifact, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return artifact


def plot_artifact(artifact: dict, output_dir: str | Path | None = None) -> Path:
    """Write the waveform comparison figure used by the notebook.

    The default is ``C:/temp`` so the public repository does not need to track a
    binary PNG.  The executed notebook embeds the rendered image in its output.
    """

    import matplotlib.pyplot as plt

    out_dir = Path(output_dir) if output_dir is not None else Path("C:/temp")
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = artifact["timeseries"]
    t = ts["time_s"]
    u = ts["input_unit_step"]
    y_cq = ts["cq_response"]
    y_ifft = ts["ifft_periodic_response"]

    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    ax.plot(t, y_cq, label="CQ causal response", lw=2.2)
    ax.plot(t, y_ifft, label="periodic IFFT contrast", lw=1.5, ls="--")
    ax.plot(t, u, label="unit step hit", lw=1.2, color="0.25", alpha=0.8)
    ax.axvline(artifact["cq"]["hit_time_s"], color="0.2", lw=1.0, alpha=0.5)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("normalized response")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    ax.set_title("URN passive relaxation model used directly by BDF2 CQ")

    fig_path = out_dir / "cq_urn_bridge_waveforms.png"
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)
    return fig_path


if __name__ == "__main__":
    artifact = build_artifact()
    figure = plot_artifact(artifact)
    print(json.dumps({
        "pass": artifact["pass"],
        "result_json": "cq_urn_bridge_results.json",
        "figure": figure.name,
        "fit_relative_error": artifact["model"]["fit"]["relative_error"],
        "cq_prehit_max_abs": artifact["cq"]["prehit_max_abs"],
        "ifft_prehit_max_abs": artifact["ifft_periodic_contrast"]["prehit_max_abs"],
    }, indent=2))
