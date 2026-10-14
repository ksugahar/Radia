# Numerical studies and research evidence

Application introductions live in the [selected gallery](../docs/APPLICATION_GUIDE.md#browse-results).
This index collects detailed studies for numerical review and paper support.
Validation drivers, inputs, acceptance criteria and result records belong in
`validation_test/`. Public theory can remain in docs; a paper plan is not a validated result.

## Study inventory

These 32 existing notebooks were selected by their numerical or methodological
focus. Inclusion is an editorial classification, not a new acceptance result.
Most presentation notebooks retain their current location while their relative
imports, assets and consumers are reviewed. Follow each notebook to its original
driver and evidence; do not treat all entries as independent executable gates.

| Study | Current presentation | Placement status |
| --- | --- | --- |
| `analytical_formulas/analytical_formulas.ipynb` | [Notebook](../docs/analytical_formulas/analytical_formulas.ipynb) | Presentation retained; dependency review precedes relocation |
| `axifem/AXIFEM_ELEMENT_EVIDENCE.ipynb` | [Notebook](axifem/AXIFEM_ELEMENT_EVIDENCE.ipynb) | Moved beside validation evidence; saved results preserved |
| `axifem/README.ipynb` | [Notebook](../docs/axifem/README.ipynb) | Presentation retained; dependency review precedes relocation |
| `bem_extractor/bem_inductance_limitations.ipynb` | [Notebook](../docs/bem_extractor/bem_inductance_limitations.ipynb) | Presentation retained; dependency review precedes relocation |
| `clebsch_hodograph/edge_focusing_tracking.ipynb` | [Notebook](../docs/clebsch_hodograph/edge_focusing_tracking.ipynb) | Presentation retained; dependency review precedes relocation |
| `clebsch_hodograph/excitation_invariant_field.ipynb` | [Notebook](../docs/clebsch_hodograph/excitation_invariant_field.ipynb) | Presentation retained; dependency review precedes relocation |
| `clebsch_hodograph/hodograph_bending_sy.ipynb` | [Notebook](../docs/clebsch_hodograph/hodograph_bending_sy.ipynb) | Presentation retained; dependency review precedes relocation |
| `clebsch_hodograph/hodograph_feasibility_2d.ipynb` | [Notebook](../docs/clebsch_hodograph/hodograph_feasibility_2d.ipynb) | Presentation retained; dependency review precedes relocation |
| `clebsch_hodograph/public_demo.ipynb` | [Notebook](../docs/clebsch_hodograph/public_demo.ipynb) | Presentation retained; dependency review precedes relocation |
| `cohomology/loop_dof_cut_selection.ipynb` | [Notebook](../docs/cohomology/loop_dof_cut_selection.ipynb) | Presentation retained; dependency review precedes relocation |
| `cohomology/tomega_wire.ipynb` | [Notebook](../docs/cohomology/tomega_wire.ipynb) | Presentation retained; dependency review precedes relocation |
| `cubit_mesh_export/netgen/p_convergence_demo.ipynb` | [Notebook](../docs/cubit_mesh_export/netgen/p_convergence_demo.ipynb) | Presentation retained; dependency review precedes relocation |
| `electric_machine/em_reference_audit.ipynb` | [Notebook](../docs/electric_machine/em_reference_audit.ipynb) | Presentation retained; dependency review precedes relocation |
| `electric_machine/planar_vim_motor.ipynb` | [Notebook](../docs/electric_machine/planar_vim_motor.ipynb) | Presentation retained; dependency review precedes relocation |
| `electrostatics/electrostatics.ipynb` | [Notebook](../docs/electrostatics/electrostatics.ipynb) | Presentation retained; dependency review precedes relocation |
| `equivalence_source/demos.ipynb` | [Notebook](../docs/equivalence_source/demos.ipynb) | Presentation retained; dependency review precedes relocation |
| `hdiv_vim/c_type_three_formulation_mixed_omega.ipynb` | [Notebook](../docs/hdiv_vim/c_type_three_formulation_mixed_omega.ipynb) | Presentation retained; dependency review precedes relocation |
| `hdiv_vim/hdiv_curved_showcase.ipynb` | [Notebook](../docs/hdiv_vim/hdiv_curved_showcase.ipynb) | Presentation retained; dependency review precedes relocation |
| `hdiv_vim/isochronous_topopt.ipynb` | [Notebook](../docs/hdiv_vim/isochronous_topopt.ipynb) | Presentation retained; dependency review precedes relocation |
| `hdiv_vim/isochronous_topopt_shape_regen.ipynb` | [Notebook](../docs/hdiv_vim/isochronous_topopt_shape_regen.ipynb) | Presentation retained; dependency review precedes relocation |
| `hysteresis/hysteresis_validation.ipynb` | [Notebook](../docs/hysteresis/hysteresis_validation.ipynb) | Presentation retained; dependency review precedes relocation |
| `induction_heating/axisymmetric_p2_thermal.ipynb` | [Notebook](../docs/induction_heating/axisymmetric_p2_thermal.ipynb) | Presentation retained; dependency review precedes relocation |
| `induction_heating/induction_heating_demo_showcase.ipynb` | [Notebook](../docs/induction_heating/induction_heating_demo_showcase.ipynb) | Presentation retained; dependency review precedes relocation |
| `kelvin/kelvin_exterior_source_and_aphi.ipynb` | [Notebook](../docs/kelvin/kelvin_exterior_source_and_aphi.ipynb) | Presentation retained; dependency review precedes relocation |
| `kelvin/Supplement/cg_smoother_demo.ipynb` | [Notebook](../docs/kelvin/Supplement/cg_smoother_demo.ipynb) | Presentation retained; dependency review precedes relocation |
| `mixed_galerkin/mixed_galerkin_results.ipynb` | [Notebook](../docs/mixed_galerkin/mixed_galerkin_results.ipynb) | Presentation retained; dependency review precedes relocation |
| `nodal_force/nodal_force.ipynb` | [Notebook](../docs/nodal_force/nodal_force.ipynb) | Presentation retained; dependency review precedes relocation |
| `peec/dowell_surface_impedance_demo.ipynb` | [Notebook](../docs/peec/dowell_surface_impedance_demo.ipynb) | Presentation retained; dependency review precedes relocation |
| `stream_function/deformation.ipynb` | [Notebook](../docs/stream_function/deformation.ipynb) | Presentation retained; dependency review precedes relocation |
| `stream_function/regularization.ipynb` | [Notebook](../docs/stream_function/regularization.ipynb) | Presentation retained; dependency review precedes relocation |
| `universal_relaxation_network/cq_urn_bridge.ipynb` | [Notebook](../docs/universal_relaxation_network/cq_urn_bridge.ipynb) | Presentation retained; dependency review precedes relocation |
| `universal_relaxation_network/urn_showcase.ipynb` | [Notebook](../docs/universal_relaxation_network/urn_showcase.ipynb) | Presentation retained; dependency review precedes relocation |

## Use the evidence

Choose the relevant [validation lane](README.md#select-an-application-sample).
Record model conditions, reference, measured quantity, tolerance and execution
environment. A stored plot or notebook output does not establish that current
code passes. Heavy reruns follow the [compute-host policy](README.md#running).
