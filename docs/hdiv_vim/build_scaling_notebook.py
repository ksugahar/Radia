"""Build + execute docs/hdiv_vim/build_scaling_hdiv_vs_mmmm.ipynb -- the HDiv-VIM vs collocation-MMMM
matrix-BUILD scalability showcase.

The notebook LOADS the mdx-measured data (build_scaling_mdx_data.json) and plots it -- it does NOT
re-run any timing (per CLAUDE.md "計算時間の測定は mdx で MUST": timing lives on the quiet mdx host,
the notebook only renders the committed data).  Executed via nbconvert (LAB is fine for plotting).

    python docs/hdiv_vim/build_scaling_notebook.py
"""
import hashlib
import json
from pathlib import Path

import nbformat
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
NB = HERE / "build_scaling_hdiv_vs_mmmm.ipynb"
SIDE = HERE / "build_scaling_hdiv_vs_mmmm_result.json"
MD = nbformat.v4.new_markdown_cell
CODE = nbformat.v4.new_code_cell

CELLS = [
    MD("""\
# HDiv-VIM vs collocation MMMM — matrix-**BUILD** scalability (HACApK charge-Gram)

The HDiv(div) Galerkin VIM has one weakness vs the collocation moment method: its charge-Coulomb **Gram
matrix build** (singular face–face integrals) is heavy, where the moment method's matrix elements are
light point-to-face integrals.  The question that decides whether HDiv-VIM is a viable *main* solver at
scale: **how far does HACApK (the charge-Gram H-matrix) accelerate the HDiv build?**

This notebook renders the **mdx-measured** answer (`build_scaling_mdx_data.json`).  Per the lab
Benchmark Policy, *timing is measured on the quiet mdx host* — this notebook only loads + plots that
committed data (no timing is run here).

**Preliminary on TET**: released radia 4.95.5 HDiv-VIM is tet-only, so this is measured on tetrahedra;
the authoritative **hex-vs-hex** comparison (the paper's element) is pending the pure-hex HDiv-VIM
reaching mdx.  The HACApK build-scalability trend (compression → near-linear build) is an H-matrix
property and is expected to carry to hex."""),

    CODE("""\
import json, numpy as np, matplotlib
%matplotlib inline
import matplotlib.pyplot as plt
data = json.load(open("build_scaling_mdx_data.json"))
rows = data["rows"]
print("measured_on:", data["measured_on"])
print("radia", data["radia_version"], "| element:", data["element"], "|", data["num_threads"], "threads")
print(f"{'n_charge':>9} {'HDiv DoF':>9} {'HDiv build[s]':>13} {'MMMM build[s]':>13} {'H-compress':>11}")
for r in rows:
    mb = "-" if r["mmmm_build_s"] is None else f"{r['mmmm_build_s']:.3f}"
    print(f"{r['n_charge']:9d} {r['hdiv_ndof']:9d} {r['hdiv_build_s']:13.3f} {mb:>13} {r['hmat_compression']:11.3f}")
print("\\nHDiv build ~ n_charge^%.2f (near-linear => scalable)" % data["hdiv_build_exponent_vs_ncharge"])"""),

    MD("""\
## Build time vs problem size, and the H-matrix compression

Left: the isolated build time (HDiv charge-Gram H-matrix `build_time` vs the MMMM
`t_moment_system_build`).  The dashed line is the fitted `~ N^{1.23}` — the HACApK charge-Gram turns a
would-be `O(N^2)` dense Gram into a **near-linear** build.  Right: the H-matrix compression ratio
(fraction of the dense Gram avoided) grows strongly with N, which is *why* the build scales."""),

    CODE("""\
nc = np.array([r["n_charge"] for r in rows], float)
hb = np.array([r["hdiv_build_s"] for r in rows], float)
mb = np.array([r["mmmm_build_s"] if r["mmmm_build_s"] is not None else np.nan for r in rows], float)
comp = np.array([r["hmat_compression"] for r in rows], float)
p = data["hdiv_build_exponent_vs_ncharge"]
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
ax[0].loglog(nc, hb, "o-", label="HDiv charge-Gram build")
ax[0].loglog(nc, mb, "s-", label="MMMM moment build")
c = hb[-1] / nc[-1]**p
ax[0].loglog(nc, c * nc**p, "k--", lw=1, label=f"~N^{p:.2f} (near-linear)")
ax[0].set_xlabel("n_charge"); ax[0].set_ylabel("build time [s]"); ax[0].legend(fontsize=7)
ax[0].set_title("matrix BUILD time (mdx, tet)")
ax[1].semilogx(nc, comp, "o-", color="C2")
ax[1].set_xlabel("n_charge"); ax[1].set_ylabel("H-matrix compression"); ax[1].set_ylim(0, 1)
ax[1].set_title("charge-Gram compression grows with N")
plt.tight_layout(); plt.show()
# crossover
ratio = hb / mb
print("HDiv/MMMM build ratio:", {int(nc[i]): round(ratio[i], 1) for i in range(len(nc)) if not np.isnan(ratio[i])})
print(data["crossover"])"""),

    MD("""\
## Conclusion — a tiered role split, not a replacement

The HACApK charge-Gram **neutralises HDiv-VIM's one weakness** (the heavy singular-Gram build): the
build scales `~N^{1.2}` and reaches parity with the moment method's build by `~34k` DoF (below that,
MMMM's build is up to ~5× faster due to H-matrix overhead + weak small-N compression).

So the honest, measured positioning is **complementary / tiered**, not one method beating the other:

| regime | main | why |
|---|---|---|
| large-scale (≳34k DoF), accuracy, loop-free, distorted | **HDiv-VIM** | build scales `~N^1.2` (HACApK), on par with MMMM; loop-free + monotone convergence |
| small-scale (<10k DoF), optimisation inner loops, coarse fast passes | **MMMM** | ~5× faster build, DoF-economical |

**Authoritative next step**: repeat this on **hex** (the paper's element) once the pure-hex HDiv-VIM is
on mdx — timing on mdx, per policy."""),

    CODE("""\
import json, platform, sys
from datetime import datetime, timezone
out = {
    "source_data": "build_scaling_mdx_data.json (measured on mdx)",
    "hdiv_build_exponent": data["hdiv_build_exponent_vs_ncharge"],
    "parity_dof": "~34k HDiv DoF",
    "conclusion": data["conclusion"],
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "python_version": sys.version.split()[0], "platform": platform.platform(),
    "note": "notebook renders mdx-measured data; no timing run in-notebook (mdx-timing policy)",
}
json.dump(out, open("build_scaling_hdiv_vs_mmmm_result.json", "w"), indent=2)
print(json.dumps(out, indent=2))"""),
]


def main():
    nb = nbformat.v4.new_notebook(cells=CELLS)
    nb.metadata.update({"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"}})
    NotebookClient(nb, timeout=300, kernel_name="python3",
                   resources={"metadata": {"path": str(HERE)}}).execute()
    nbformat.write(nb, NB)
    sha = hashlib.sha256(NB.read_bytes()).hexdigest()
    d = json.loads(SIDE.read_text()); d["notebook_sha256"] = sha
    SIDE.write_text(json.dumps(d, indent=2))
    print("wrote", NB.name, "+ sidecar sha", sha[:16])


if __name__ == "__main__":
    main()
