"""Shared FEM reference text for the Radia MCP topic dispatchers."""

CORE_TOPICS = {
    "index": (
        "ngsolve-mcp.core: NGSolve FEM spaces, solvers, and TaskManager.\n"
        "Topics: spaces, solvers, taskmanager, eddy_current."
    ),
    "spaces": (
        "NGSolve finite element spaces:\n"
        "- H1(mesh, order=p, dirichlet=...): scalar nodal potential (phi, T)\n"
        "- HCurl(mesh, order=p, dirichlet=..., nograds=...): vector potential A / E-field\n"
        "- HDiv(mesh, order=p): flux density B / current density J\n"
        "- L2(mesh, order=p-1): discontinuous scalar/vector fields"
    ),
    "solvers": (
        "NGSolve linear solvers & preconditioners:\n"
        "- Direct: sparsecholesky for suitable symmetric systems; check backend availability and matrix assumptions.\n"
        "- Iterative: CG / GMRES with bddc / multigrid / local preconditioners\n"
        "- External: sparsesolv (compact_ams / compact_amg)"
    ),
    "taskmanager": (
        "Wrap assembly and solves with `with TaskManager():` for shared-memory\n"
        "parallelism; call `SetNumThreads(n)` before entering TaskManager if needed."
    ),
    "eddy_current": (
        'See ngsolve-mcp.weak-form: weak_form_usage(topic="eddy_current") '
        'for the A-phi formulation and weak_form_usage(topic="ih_reference") '
        "for independent induction-heating FEM checks."
    ),
}


WEAK_FORM_TOPICS = {
    "overview": (
        "ngsolve-mcp.weak-form: NGSolve variational weak-form patterns.\n"
        "Topics: magnetostatics_a, eddy_current, sibc_robin, ih_reference, hdiv_flux, symbolic_bfi."
    ),
    "magnetostatics_a": (
        "Vector-potential A-formulation in HCurl:\n"
        "  a += nu * InnerProduct(curl(u), curl(v)) * dx\n"
        "  f += InnerProduct(J, v) * dx"
    ),
    "eddy_current": """# Frequency-domain A-phi weak form

Use exp(+i omega t) and scaled phi = electric potential / (i omega), so
E = -i omega (A + grad(phi)). Use HCurl for A and an H1 scalar space on
conductors for phi. Fix the scalar gauge on each independent conductor and
handle the HCurl gradient nullspace consistently with the port constraints.

With (A, phi) and (v, q) from the compound trial and test spaces, use the
current NGSolve integrator syntax (nu = inverse permeability, sigma = conductivity):

    a = BilinearForm(fes)
    a += nu * InnerProduct(curl(A), curl(v)) * dx
    a += 1j * omega * sigma * InnerProduct(A + grad(phi), v + grad(q)) * dx("conductor")

Use the supplied conductor region selector in place of the example label.
This bilinear form is complex symmetric with real reciprocal material
coefficients and compatible constraints. The trial/test proxies above are
not conjugated; in power or energy post-processing use a conjugated field
product without double conjugation by InnerProduct. The q test enforces
charge conservation; specify terminal and insulating boundary conditions
before choosing source terms. A-only variants require separate assumptions.
""",
    "sibc_robin": """# SIBC as a Robin term in FEM

In a frequency-domain FEM weak form, integrate the impedance contribution
over the selected boundary using tangential traces. For an A-only form:

    a += alpha * InnerProduct(A.Trace(), v.Trace()) * ds("sibc")

Here alpha is the coefficient after eliminating the exterior tangential
magnetic field with the chosen SIBC. Derive its sign and frequency factor
from the time convention, normal orientation and E-to-A relation; do not
substitute Zs itself. For A-phi, include the tangential scalar gradient where
the impedance law uses the total tangential electric field. Verify that the
chosen convention gives nonnegative dissipated power for a passive surface.

Example: let the FEM domain be the air/coil side, with its outward normal n
pointing into the conductor. With exp(+i omega t) and E_t = -i omega A_t,
the Leontovich condition using the conductor normal -n gives
n cross H = (i omega / Zs) A_t. For a good conductor with skin depth delta
and conductivity sigma, Zs = (1 + i) / (sigma * delta):

    a += (1j * omega / Zs) * InnerProduct(A.Trace(), v.Trace()) * ds("sibc")

Re(i omega / Zs) = omega * sigma * delta / 2 > 0 for this convention.
Derive the coefficient again if the normal, phasor convention or unknowns differ.""",
    "ih_reference": """# Independent NGSolve FEM reference for induction heating

For an axisymmetric reference, use the azimuthal A formulation with the
axisymmetric integration weight and the correct axis condition. For a fully
three-dimensional coil use HCurl A-phi with conductor and terminal constraints.
Resolve the skin layer for a volume reference; document the assumptions and
impedance coefficient if using the independent SIBC alternative.

Match geometry, frequency, material properties and excitation conventions
before comparing. Report complex coil impedance from V/I and time-average
workpiece power from the conductor Joule integral or the consistent surface
loss integral. Record whether amplitudes are peak or RMS. Check power balance,
gauge and mesh convergence independently. This topic is a reference workflow
description, not a numerical qualification or a production Radia route.""",
    "hdiv_flux": (
        "Mixed HDiv-L2 flux-potential formulation:\n"
        "  Enforces normal flux continuity across material interfaces."
    ),
    "symbolic_bfi": (
        "Legacy syntax reference. For current code, add the integrand times dx or ds\n"
        "directly to BilinearForm or LinearForm; see magnetostatics_a and eddy_current."
    ),
}


GAUGING_TOPICS = {
    "overview": (
        "ngsolve-mcp.gauging: Electromagnetic gauging reference and gates.\n"
        "Topics: tree_cotree, coil_av_tree, coulomb, open_boundary."
    ),
    "tree_cotree": (
        "Tree-Cotree gauging for HCurl magnetostatics/eddy-current formulations:\n"
        "Remove gradient-space nullspace (ker curl) by eliminating spanning-tree\n"
        "edge DoFs on the mesh graph before solving."
    ),
    "coil_av_tree": """# Tree gauge for a voltage-driven coil

Keep source and sink ports separate in the A-V space. Build the tree from
active free edges, excluding Dirichlet edges and periodic slaves. Contract
periodic copies consistently. Preserve tangential A on each connected SIBC
trace by contracting each trace component separately in the tree graph.

Choose the shortest admissible source-to-sink mesh-edge path off the SIBC
trace. Hold the line integral of A along that path fixed to enforce the
terminal voltage convention. Reject the setup if no admissible path exists. Do not use a
tree that short-circuits the ports or eliminates the prescribed terminal drop.
For disconnected free-edge graphs use one spanning tree per component.

Compare two admissible trees with distinct free-edge masks. Verify terminal
impedance, field energy, current-normalized inductance, and both coil and
workpiece losses. Check the true residual and the unnormalized terminal
current: a small residual can coexist with an incorrect port drive. Include
internal reactance when comparing with a BEM or surface-current air-coil
reference, and check loss convergence separately. The IH voltage scaling, production recipe and comparison criteria
remain in the IH application guidance.""",
    "coulomb": (
        "Coulomb gauge (div A = 0) and grad-div stabilization:\n"
        "Enforce div A = 0 via Lagrange multiplier or augmented grad-div term."
    ),
    "open_boundary": (
        "Open-boundary gauging with Kelvin transformation or BEM-FEM coupling:\n"
        "Ensure exterior potential decay and compatible gauge across interfaces."
    ),
}


KELVIN_TOPICS = {
    "overview": (
        "ngsolve-mcp.kelvin: Kelvin transformation for unbounded magnetostatic/eddy-current domains.\n"
        "Topics: inversion, jacobian, interface_coupling, verify_first."
    ),
    "inversion": (
        "Spherical/ellipsoidal Kelvin inversion maps exterior r > R0 to interior r' = R0^2 / r,\n"
        "preserving harmonicity with scaled material weight (R0 / r')^2."
    ),
    "jacobian": (
        "In NGSolve, apply the Jacobian-weighted bilinear form on the inverted region\n"
        "so that exterior energy and far-field decay match exact open-boundary conditions."
    ),
    "verify_first": """# Verify-first checks for a periodic Kelvin pair

These are checks 2 and 3 of the IH example. Check 1 (mesh and labels) stays
in the IH recipe.

2. Periodic slave DOFs: build the underlying HCurl space with the intended
   essential constraints, then wrap it with Periodic. Compare FreeDofs().NumSet()
   before and after. A nonempty Kelvin identification should remove free slave
   DOFs in this setup. Confirm expected labels, pairing, orientation and phase
   before proceeding; a lower count alone does not certify the pairing.

3. Interior/exterior ratio: set a nonzero tangential test field on kelvin_int in
   the periodic space and evaluate its squared trace norms on kelvin_int and
   kelvin_ext with Integrate(InnerProduct(g, g), mesh, BND, definedon=...).
   NGSolve InnerProduct conjugates complex coefficient functions in this
   evaluation; do not add Conj to the trace. Reject zero interior norm. Compare
   norm_ext / norm_int to one for the matching unit-phase Kelvin interface.
   The existing IH verify-first example uses abs(ratio - 1.0) < 1e-3; do not
   generalize this criterion to interfaces with different weighting or phase.""",
    "interface_coupling": (
        "Match tangential trace (HCurl) or normal flux (HDiv) continuously across r = R0\n"
        "between the physical interior mesh and the Kelvin-inverted exterior mesh."
    ),
}


DTN_TOPICS = {
    "overview": (
        "ngsolve-mcp.dtn: Dirichlet-to-Neumann (DtN) open-boundary condition.\n"
        "Topics: spherical_harmonics, coarse_mesh, bem_fem_comparison."
    ),
    "spherical_harmonics": (
        "Exact DtN map on a spherical artificial boundary r = R0 expands exterior\n"
        "harmonic modes in spherical harmonics Y_l^m with radial decay (R0/r)^(l+1)."
    ),
    "coarse_mesh": (
        "Use curved high-order boundary elements on r = R0 so DtN quadrature\n"
        "remains accurate even on coarse interior meshes."
    ),
    "bem_fem_comparison": (
        "DtN is diagonal in spherical modes (ideal for spherical truncation);\n"
        "BEM-FEM Schur or Kelvin transformation handles arbitrary boundary shapes."
    ),
}
