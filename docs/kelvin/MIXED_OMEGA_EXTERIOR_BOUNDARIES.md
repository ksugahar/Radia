# Finite exterior boundaries in mixed scalar potential solves

The default natural boundary specifies zero **total** normal flux.
`reduced_zero_normal_boundary` instead specifies zero normal **correction**
field, so that total normal H equals source normal H. These are different
physical conditions. The latter requires the fixed source-normal boundary
term in the weak right-hand side. Only named exterior faces of reduced
regions are accepted; internal interfaces are rejected.

A uniform-source control distinguishes these conditions analytically
(`tests/test_kelvin_mixed_omega.py::test_reduced_normal_boundary_is_not_total_normal_boundary`).
