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

Run: `python two_families.py --out result.json --configs 4:2,8:3` (a config is
`n:k` or `NXxNYxNZ:k`, e.g. `1x1x1:6,2x1x1:6` for one hex and two hexes sharing
the face x = 0).
Gmsh is not used; any display export is post-processing only.

## Derivation (`derivation.wls`)

`wolframscript -file derivation.wls` builds on the repository FEEC library
(`radia_mcp/mathematica/basis_functions/hcurl.wls`) and prints every check:

- S1: the analytic field identities used by the code. The z = pi map is the
  identity, so a patch outlet at z = pi has no curvature. The field line
  (u = c, v = 0) lies in the element face x = c over its whole length.
- S2: for the hex Nedelec/RT pair (p = 2..4), curl maps HCurl into HDiv and
  div curl = 0. A zero tangential trace forces a zero normal curl, so B_h.n is
  continuous. A_n and the tangential curl stay free, so they jump across faces.
  A_h.dl along a chord in a z = const plane has degree 2p+1, so 6-point Gauss
  per face-split piece is exact.
- S3: discrete Clebsch. With u, v in H1 (degree p), u grad v is H(curl)
  conforming and grad u x grad v lies in HDivHexT[2p+1] but not in [2p]. The
  wedge is conforming and exactly div-free, but a given B_h of order k is not
  in general a single wedge of order-k labels.
- S4: d/dz(B_z det J) = det J div B.
- S5: reconstruction leakage is O(dt^2) + O(dz^3) and geometry error is
  O(dt^2) + O(dz^4), from series on the analytic surfaces. The polygon flux
  error over one chord is h^3 (smooth), h^2 (kink), h^(3/2) (fold) or h (jump).
  A grazing crossing with tangential jump J shifts the outlet by J*2 sqrt(d/c).

A field line lying in a face plus the tangential jump J of B_h gives the discrete
flow map a steep layer. Its width is about the B_h error, and across it the
outlet moves by about J times the hugging length. Until the loop spacing h
resolves this layer, the transported patch flux error is the jump-type O(J h).

## Display (`display_surfaces.py`)

Display only. On the two-hex case (default `2x1x1:6`) the script tracks five
surfaces of each family from the inlet section through B_h = curl(A_h) and
writes `<out>.msh` (MSH 4.2 ASCII for Gmsh) plus `<out>.json`. The JSON
carries the same diagnostics as `two_families.py`: projection, signed patch
flux at z = 0, pi/2, pi, and per surface the tracked and oracle leakage
together with the crossings of the shared face. The `.msh` carries the
surfaces, their field lines, the element edges and the shared face x = 0.
Its views are the distance to the exact surface (nodes) and the family label
and local |B_h.n|/|B_h| (quad centres). The quantitative leakage is the
integrated value in the JSON, not the view.

    python display_surfaces.py --out <dir>/display --threads 1
