# Two magnetic-surface families (prototype)

Diagnostic prototype, not an adopted validation case: no acceptance criterion
is defined here.

`two_families.py` reconstructs both Clebsch magnetic-surface families of the
nonaxisymmetric analytic map

    x = u + a sin(v) sin(s),  y = v + b sin(x) sin(2 s),  z = s,  a = b = 0.5,

whose field B = dX/ds = grad(u) x grad(v) has det(DX) = 1. The labels (u, v)
serve only as the oracle; the reconstruction sees the B evaluator and the inlet
section z = 0.

It reports three error sources separately:

1. projection of B into B_h = curl(A_h) (HCurl order k, structured hexahedra,
   simply connected box, natural boundary condition), with cell divergence and
   one-sided face normal flux of the HDiv copy;
2. the face-aware Gauss-Legendre tracker against an independent adaptive DOP853
   trace of the same B_h, and transported patch flux (closed A_h line integral);
3. surface reconstruction of both families under seed and along-line
   refinement, with signed and absolute normal flux and oracle surface distance;
   oracle-sampled surfaces with the exact B isolate pure reconstruction error.

Run: `python two_families.py --out result.json --configs 4:2,8:3`.
Gmsh is not used; any display export is post-processing only.
