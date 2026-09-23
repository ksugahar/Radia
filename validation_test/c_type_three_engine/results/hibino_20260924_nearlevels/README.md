# Near quadrature sensitivity across three mesh levels

HDiv-only diagnostic, not three-engine acceptance or an accuracy certificate.
Coarse and medium checked meshes each run default/intorder9/intorder13,
BDM1, mu_r1000, 8 threads, Gram eps1e-14, solve tolerance1e-10.
Fine fields are reused from ../hibino_20260924_nearquad/near_quad_comparison.json.
All six new native solves converged. JSON records residuals and timings.

The historical difference key ending `previous_fine` in these two level JSONs
actually uses the same-level default result; baseline_note explicitly states this.
increments.json applies the campaign median-plane axial projection/core mask,
normalizing both increments by each setting's fine vector norm.
Coarse-to-medium stays about 0.0163%; medium-to-fine is 1.00149% default,
1.18604% at intorder9, and 1.20047% at intorder13. Increased Gram quadrature
does not restore contraction. This does not identify a specific discretization
bug and is not a speed superiority claim. Next: investigate FE-order sensitivity.

Raw archive S:/Radia/validation_artifacts/hdiv_nearlevels_20260924/recovery.zip
SHA256 944d4a218803e93350d14fa1d97e61f1dc57c82e6d3e8828132d23106ec16dbe
matches LAB and hibino. Includes wheel, driver, checked meshes, command and logs.
Dedicated remote directory and ZIP to be removed after commit; shared mesh
family and system installation remain unchanged.
