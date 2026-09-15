"""Prepare and execute the ESIM spatial notebook on a foreground mdx job."""
import argparse
import json
import os
from pathlib import Path
import sys
import nbformat as nbf


def prepare(root):
    cells = []
    def md(s): cells.append(nbf.v4.new_markdown_cell(s))
    def code(s): cells.append(nbf.v4.new_code_cell(s))
    md(r'''# Where does an ESIM workpiece heat up?

A coil heats a conducting steel cylinder unevenly. Radia's **PEEC coil +
nonlinear ESIM/BEM workpiece** route resolves surface-field-dependent impedance
without meshing every skin depth in the conductor. This notebook shows the
actual coil CAD, workpiece mesh, accepted tangential field, and surface-loss
distribution for **100 A at 10, 50 and 100 kHz**. The saved results come from
a fresh mdx2 calculation, not a prescribed-field illustration.

The **Radia MCP IH/ESIM usage tools** own current operating instructions,
contracts and limitations; Python here is the LLM-driven reproduction layer.
Simulink remains the formal human UI. This is a field-inspection notebook,
not an application workbench. See also [the ESIM benchmark](esim_showcase.ipynb).

## Model and equations

The coil is the bundled `ih_fem_kelvin_demo_coil.step`; the workpiece is the
bundled `ih_bem_sample_p1.vol`, with boundary label `sibc`. Conductivity is
$2\times10^6$ S/m; the supplied `em_sample_bh.txt` controls nonlinear steel
response. The cell half-thickness is 5 mm. These are demonstration material
parameters, not a measured specimen specification.

The original P1 volume file also contains coarse air and coil regions that
this PEEC/BEM route does **not** solve on. Their CAD discrepancies reach
8.45%, so the whole-file exploratory check uses a declared 10% tolerance;
the actual workpiece volume and `sibc` area must independently be within 1%
of their CAD references. This is not a volume-FEM quality certification.

For peak-amplitude phasors with $e^{j\omega t}$ convention, diffusion through
the surface-normal coordinate $s$ gives, for a linear material,

$$\frac{d^2 H_t}{ds^2}=j\omega\mu\sigma H_t,\qquad
\delta=\sqrt{\frac{2}{\omega\mu\sigma}},\qquad
Z_s\simeq\frac{1+j}{\sigma\delta}.$$

ESIM replaces the constant impedance with the nonlinear cell response
$Z_s(|H_t|,\omega)$ and updates it with the surface field
[@hollaus2026nonlinear]. Surface impedance removes the thin skin layer from
the outer spatial discretization, subject to the cell and curvature assumptions
[@yuferev2009surface]. Tangential electric field and cycle-averaged loss are

$$\mathbf E_t=Z_s(\mathbf n\times\mathbf H_t),\qquad
q_{\mathrm{local}}=\tfrac12\operatorname{Re}(Z_s)|H_t|^2,\qquad
P_{\mathrm{local}}=\int_\Gamma q_{\mathrm{local}}\,dS.$$

The nonlinear closure is iterated with relaxation and safeguarded Anderson
acceleration. Every displayed case must have `esim_converged=true`.
This is the **weak-coupled** route: do not infer strong coil back-reaction
accuracy, thermal evolution, or mesh convergence from solver convergence.

### Two different surface-heating quantities — do not conflate them

The local ESIM diagnostic above uses the saved accepted $|H_t|$ and
$\operatorname{Re}Z_s$, reconstructed in the solver's P1 vertex ordering.
Separately, the current solver writes a heat-transfer artifact with an incident
Biot–Savart pattern normalized to its BIE total power:

$$q_{\mathrm{transfer}}(x)=
P_{\mathrm{BIE}}\frac{|H_{t,\mathrm{incident}}(x)|^2}
{\int_\Gamma |H_{t,\mathrm{incident}}|^2\,dS}.$$

We display **both**, label them explicitly, and report their integrals.
Agreement of the transfer integral with total power is conservation by
construction, not an independent validation of the local heating pattern.
No volumetric eddy-current or temperature comparison is claimed here.
''')
    code('''import json, os
from pathlib import Path
import ngsolve as ng
from ngsolve.webgui import Draw
from IPython.display import display, Markdown
from esim_spatial_demo import load_fields
data = Path("spatial_demo_data")
samples = Path(os.environ.get("RADIA_ESIM_INPUTS", "../../src/radia/panels/samples"))
with ng.TaskManager():
    mesh, cases, check = load_fields(data, samples)
print("check-vol:", check["passed"], "vertices:", mesh.nv, "boundaries:", mesh.GetBoundaries())''')
    md('''## Coil and workpiece geometry

The first scene combines the source STEP coil with the cylindrical workpiece
geometry (radius 25 mm, height 25 mm, centred at the origin, as specified by
`ih_bem_sample.jou`). The cylinder bounds are checked against the actual mesh.
The second scene displays that solver workpiece mesh. Source coordinates are
in metres, with no deliberate axis exaggeration.
The workpiece scenes show only the surface: ESIM does not compute an interior
temperature or a volumetric loss field.''')
    code('''from netgen.occ import OCCGeometry, Cylinder, Pnt, Dir, Compound
from netgen.webgui import Draw as DrawGeometry
import numpy as np
coil_cad = OCCGeometry(str(samples / "ih_fem_kelvin_demo_coil.step"))
xyz = np.asarray([v.point for v in mesh.vertices])
assert np.allclose(xyz.min(axis=0), [-0.025,-0.025,-0.0125], atol=1e-6)
assert np.allclose(xyz.max(axis=0), [0.025,0.025,0.0125], atol=1e-6)
workpiece_cad = Cylinder(Pnt(0,0,-0.0125), Dir(0,0,1), r=0.025, h=0.025)
workpiece_cad.col = (0.45,0.55,0.7)
coil_cad.shape.col = (0.85,0.45,0.15)
coil_scene = DrawGeometry(Compound([coil_cad.shape, workpiece_cad]), width="100%", height="440px")''')
    code('''with ng.TaskManager():
    mesh_scene = Draw(mesh, name="ESIM_workpiece_mesh", draw_vol=False,
                      draw_surf=True, clipping=None,
                      width="100%", height="440px")''')
    md('''## Field and local surface-loss comparison

All three frequencies use the same current and geometry. Each quantity has a
**shared colour range across frequencies**; a brighter view therefore means
a larger physical value, not merely a rescaled colour bar. Magnetic field is
in A/m, surface resistance in ohms, and heating density in W/m².''')
    code('''hmax = max(max(c["result"]["esim_per_panel_H_t"]) for c in cases)
rmax = max(max(c["result"]["esim_per_panel_Z_s_real"]) for c in cases)
qmax = 0.5*rmax*hmax**2
scenes = []
for c in cases:
    display(Markdown(f"### {c['frequency_hz']/1000:g} kHz — accepted tangential field"))
    with ng.TaskManager():
        scenes.append(Draw(c["Ht"], mesh, name=f"Ht_A_per_m_{c['frequency_hz']}Hz",
                           draw_vol=False, draw_surf=True, autoscale=False,
                           min=0, max=hmax, clipping=None,
                           width="100%", height="440px"))''')
    code('''for c in cases:
    display(Markdown(f"### {c['frequency_hz']/1000:g} kHz — local ESIM diagnostic, not transfer artifact"))
    with ng.TaskManager():
        scenes.append(Draw(c["q_local"], mesh, name=f"q_local_W_per_m2_{c['frequency_hz']}Hz",
                           draw_vol=False, draw_surf=True, autoscale=False,
                           min=0, max=qmax, clipping=None,
                           width="100%", height="440px"))''')
    md('''## Why local impedance and the transfer artifact matter

At 50 kHz inspect the spatial resistance variation, then compare the
power-normalized transfer artifact with the local diagnostic above. The latter
does not have the same colour range: use the labels and integral table, not
colour alone, to compare these two definitions.''')
    code('''c = cases[1]
with ng.TaskManager():
    resistance_scene = Draw(c["resistance"], mesh, name="Re_Zs_ohm_50000Hz",
                            draw_vol=False, draw_surf=True, autoscale=False,
                            min=0, max=rmax, clipping=None,
                            width="100%", height="440px")
    transfer_scene = Draw(c["heat"], mesh, name="q_transfer_W_per_m2_50000Hz",
                          draw_vol=False, draw_surf=True, autoscale=True,
                          clipping=None, width="100%", height="440px")''')
    md('''## Power and acceptance evidence

These are distinct observables, not two independent solver validations.
Relative transfer-integral discrepancy is measured against the reported BIE
power. Full inputs, source text/hashes, commands, runtime and convergence
records accompany the notebook. The earlier analytical Bessel cross-check
and historical sweep remain in the main benchmark, with their own provenance.''')
    code('''rows = []
for c in cases:
    d = c["result"]
    rows.append(dict(frequency_hz=c["frequency_hz"], current_A=d["current_A"],
                     converged=d["esim_converged"], iterations=d["esim_iterations"],
                     P_BIE_W=d["P_wp_W"], P_local_W=c["P_local_W"],
                     P_transfer_W=c["P_heat_W"],
                     transfer_relative_gap=abs(c["P_heat_W"]/d["P_wp_W"]-1)))
display(Markdown("| f (kHz) | iterations | BIE (W) | local diagnostic (W) | transfer (W) | transfer gap |\\n"
                 "|---:|---:|---:|---:|---:|---:|\\n" + "\\n".join(
    f"| {r['frequency_hz']/1000:g} | {r['iterations']} | {r['P_BIE_W']:.6g} | {r['P_local_W']:.6g} | {r['P_transfer_W']:.6g} | {r['transfer_relative_gap']:.3e} |" for r in rows)))
(data / "spatial_metrics.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")''')
    # Discovery is visual-first; keep the derivation available below the scenes.
    opening, marker, theory = cells[0].source.partition('## Model and equations')
    cells[0].source = opening
    cells.append(nbf.v4.new_markdown_cell(marker + theory))
    nb=nbf.v4.new_notebook(cells=cells)
    nb.metadata.update(kernelspec=dict(name='python3',display_name='Python 3',language='python'),
        radia=dict(notebook_role='example',webgui_required=True,webgui_field_required=True,
                   citation_keys=['hollaus2026nonlinear','yuferev2009surface']))
    nbf.write(nb, root/'esim_spatial_demo.ipynb')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['prepare','execute'])
    p.add_argument('root',type=Path)
    args=p.parse_args()
    if args.mode == 'prepare': prepare(args.root)
    else:
        from nbclient import NotebookClient
        os.environ['JUPYTER_PATH']=str(Path(sys.prefix)/'share/jupyter')
        os.environ['RADIA_ESIM_INPUTS']=str(args.root.resolve())
        path=args.root/'esim_spatial_demo.ipynb'
        nb=nbf.read(path,as_version=4)
        NotebookClient(nb,timeout=180,kernel_name='python3',
                       resources={'metadata':{'path':str(args.root.resolve())}}).execute()
        nbf.write(nb,path)
