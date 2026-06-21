# -*- coding: utf-8 -*-
"""Paper-quality figures for the DtN mesh-cutting (act2_11) and form-resolved p/h-convergence
(act2_12) datasheets, via the lab figure toolkit (radia_mcp.figure: IEEJ single-column profile +
no-in-figure-title + axes-area + legend-overlap quality gates).  Reads the COMMITTED JSON datasheets
next to it and writes <name>.pdf/.png beside them.  No in-figure title (lab rule).

Run: python plot_dtn_mesh_studies.py
"""
import os
import json
from radia_mcp.figure import paper_figure, emit_paper_figure
from radia_mcp.figure.tools import find_best_legend_loc

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with open(os.path.join(HERE, name)) as f:
        return json.load(f)


def place_legend(ax, **kw):
    loc, _ = find_best_legend_loc(ax)
    ax.legend(loc=loc, frameon=False, **kw)


d12 = load("act2_12_derham_dtn_pconv_hconv.json")
d11 = load("act2_11_kelvin_mesh_cutting_datasheet.json")

# ===========================================================================
# Figure A (act2_12): the DtN spectrum is FORM-dependent + the p-method
# ===========================================================================
figA, axesA = paper_figure(profile='ieej_single_column', nrows=1, ncols=2, panel_labels=True)
axLad, axP = axesA.flat

# (a) the two distinct de Rham ladders: DtN eigenvalue vs degree n
modes = d12["modes"]
axLad.plot(modes, d12["A_radial_ie"]["scalar_dtn"], 'o-', label=r'$H^1,\,H(\mathrm{div})_n$')
axLad.plot(modes, d12["A_radial_ie"]["toroidal_dtn"], 's--', label=r'$H(\mathrm{curl})_t$')
axLad.set_xlabel(r'multipole degree $n$')
axLad.set_ylabel(r'DtN eigenvalue $\lambda\,R$')
axLad.set_xticks(modes)
place_legend(axLad)

# (b) p-convergence (Kelvin is a p-method): defect vs element order
pc = d12["B_h1_fem"]["p_convergence"]
po = sorted(int(k) for k in pc)
axP.semilogy(po, [pc[str(p)]["d1"] for p in po], 'o-', label=r'$H^1$')
hc = d12["B_hcurl_fem"]
ho = sorted(int(k) for k in hc)
axP.semilogy(ho, [hc[str(o)]["rel"] for o in ho], 's--', label=r'$H(\mathrm{curl})$')
axP.set_xlabel(r'element order $p$')
axP.set_ylabel(r'rel. DtN defect')
axP.set_xticks([1, 2, 3, 4])
place_legend(axP)
emit_paper_figure(figA, os.path.join(HERE, 'act2_12_derham_dtn_pconv_hconv'),
                  profile='ieej_single_column', min_axes_fraction=0.72, on_fail='auto_tighten')
print("wrote act2_12_derham_dtn_pconv_hconv.{pdf,png}")

# ===========================================================================
# Figure B (act2_11): the min(p,k) law -- defect COLLAPSES onto min(p,k)
# ===========================================================================
ds = d11["datasheet"]


def mkey(dbm, h):                                   # JSON maxh key whose float equals h
    return next(kk for kk in dbm if abs(float(kk) - h) < 1e-12)


maxh_list = sorted(d11["maxh_list"])
nprobe = str(d11["n_probe"])
markers = ['o', 's', '^', 'D', 'v']

figB, axesB = paper_figure(profile='ieej_single_column', nrows=1, ncols=1)
axM = axesB.flat[0]
# finest-mesh defect at n=2 vs min(p,k): the three min=2 configs land together
for (cfg, mk) in zip(d11["configs"], markers):
    p, k = cfg
    dbm = ds[f"{p}_{k}"]["defect_by_maxh"]
    d2 = dbm[mkey(dbm, maxh_list[0])][nprobe]
    axM.semilogy([min(p, k)], [d2], mk, markersize=8, label=fr'$(p,k)=({p},{k})$')
axM.set_xlabel(r'$\min(p,k)$')
axM.set_ylabel(r'rel. DtN defect ($n{=}2$)')
axM.set_xticks([1, 2, 3])
place_legend(axM)
emit_paper_figure(figB, os.path.join(HERE, 'act2_11_kelvin_mesh_cutting_datasheet'),
                  profile='ieej_single_column', min_axes_fraction=0.72, on_fail='auto_tighten')
print("wrote act2_11_kelvin_mesh_cutting_datasheet.{pdf,png}")
