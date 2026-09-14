# See what Cubit Mesh Export makes possible

A curved surface does not have to become a dense collection of flat tetrahedra.
Cubit Mesh Export carries existing Cubit meshes into solver-facing `.vol`
files, including the curved HEX geometry demonstrated here.

Start with the [interactive mesh showcase](cubit_mesh_export_showcase.ipynb):

- Rotate the same 56-HEX sphere at geometry orders 1, 2 and 3.
- Compare the faceted silhouette with curved faces and edges.
- Inspect the actual material and boundary labels and saved `check-vol` reports.
- Look inside a manufactured Poisson solution on the cubic curved HEX mesh.

The notebook saves the scenes and numerical outputs. A WebGUI-capable notebook
viewer is needed for interaction; GitHub's static preview may not activate widgets.
These are reloaded, committed historical Cubit exports, not a newly certified
export from today's plugin build. Geometry checks and the scalar PDE example
do not certify every mixed-element mesh or electromagnetic FE space.

The [p-convergence notebook](netgen/p_convergence_demo.ipynb) provides a
complementary saved mesh-evaluation campaign.
The [mesh inputs](hex_sphere_highorder/) remain available for reproduction.

## Where to go next

The **`radia-cubit` MCP server is the operating manual**. Use `cubit_status`
and `netgen_workflow_guide` for current export routes, high-order geometry
constraints and validation. `radia-analysis` owns downstream FEM workflows.
Notebook Python is the LLM-driven reproduction layer.

Numerical checks and historical campaign summaries live under
[`validation_test/cubit_mesh_export/`](../../validation_test/cubit_mesh_export/).
The [export package](../../packages/cubit-mesh-export/) owns the implementation.
GMSH and Nastran files are downstream interchange artifacts; `.vol` is the
NGSolve mesh route.
