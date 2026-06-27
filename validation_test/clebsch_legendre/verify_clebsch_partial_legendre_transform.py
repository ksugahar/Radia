"""Symbolic verification of the PARTIAL Legendre (hodograph) transformation
of Clebsch-potential magnetostatics: swap only the pair z <-> psi.

Forward formulation (independent variables x, y, z):

    B = grad(phi) x grad(psi),   curl B = 0   (current-free region)

Transformed formulation (independent variables x, y, psi; unknowns
z(x,y,psi) and phi(x,y,psi); admissible where psi_z != 0 so the flux
surfaces psi = const are graphs z = z(x,y,psi)):

    J = z_psi
    grad(psi) = ( -z_x, -z_y, 1 ) / z_psi
    B = (1/z_psi) ( phi_y, -phi_x, z_x phi_y - z_y phi_x )

    covariant components (coordinates u = (x, y, psi)):
    b_x   = [ (1 + z_x^2) phi_y - z_x z_y phi_x ] / z_psi
    b_y   = [ z_x z_y phi_y - (1 + z_y^2) phi_x ] / z_psi
    b_psi = z_x phi_y - z_y phi_x

    (F1)  d(b_y)/dx   = d(b_x)/dy
    (F2)  d(b_psi)/dy = d(b_y)/dpsi
    (F3)  d(b_psi)/dx = d(b_x)/dpsi

Energy functional:

    W = 1/(2 mu0) Int [ phi_x^2 + phi_y^2 + (z_x phi_y - z_y phi_x)^2 ]
                      / z_psi  dx dy dpsi

whose Euler-Lagrange equations are the force-free projections

    EL_phi = -2 C_psi                         (C_psi = d(b_y)/dx - d(b_x)/dy)
    EL_z   = -2 (C_x phi_x + C_y phi_y) / z_psi

with C_x = d(b_psi)/dy - d(b_y)/dpsi, C_y = d(b_x)/dpsi - d(b_psi)/dx.

In the 2D limit (B_z = 0, no z-dependence) the Clebsch pair is
(phi, psi) = (A_z, z): the surface unknown becomes trivial (z = psi)
and phi reduces to the familiar magnetic vector potential component
A_z(x, y) — this formulation is its 3D generalization.

Run:  python verify_clebsch_partial_legendre_transform.py
"""

import sympy as sp


def check(name, expr):
    residual = sp.simplify(expr)
    ok = residual == 0
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {name}")
    if not ok:
        print("  residual:", residual)
    assert ok, name


def check_vec(name, vec):
    for i, comp in enumerate(vec):
        check(f"{name} [{'xyz'[i]}]", comp)


# ---------------------------------------------------------------------------
# Check 1: pointwise Jacobian algebra
#   Coordinates u = (x, y, psi), position r = (x, y, z(x,y,psi)),
#   phi = phi(x, y, psi).  Verify the B formula and that phi_psi drops out.
# ---------------------------------------------------------------------------
print("-- Check 1: B = (phi_y, -phi_x, z_x phi_y - z_y phi_x)/z_psi (algebra)")

zx, zy, zp = sp.symbols("z_x z_y z_psi")
fx, fy, fp = sp.symbols("phi_x phi_y phi_psi")

A = sp.Matrix([[1, 0, 0],
               [0, 1, 0],
               [zx, zy, zp]])       # A_ij = d x_i / d u_j
check("J = z_psi", A.det() - zp)

Ainv = A.inv()                       # rows = grad x, grad y, grad psi
grad_x = sp.Matrix(Ainv[0, :]).T
grad_y = sp.Matrix(Ainv[1, :]).T
grad_psi = sp.Matrix(Ainv[2, :]).T
check_vec("grad(psi) - (-z_x, -z_y, 1)/z_psi",
          grad_psi - sp.Matrix([-zx, -zy, 1]) / zp)

grad_phi = fx * grad_x + fy * grad_y + fp * grad_psi
B = grad_phi.cross(grad_psi)
B_formula = sp.Matrix([fy, -fx, zx * fy - zy * fx]) / zp
check_vec("grad(phi) x grad(psi) - B_formula", B - B_formula)
check_vec("B independent of phi_psi", sp.diff(B, fp))

# covariant components
r_x = sp.Matrix([1, 0, zx])
r_y = sp.Matrix([0, 1, zy])
r_p = sp.Matrix([0, 0, zp])
check("b_x", B_formula.dot(r_x) - ((1 + zx**2) * fy - zx * zy * fx) / zp)
check("b_y", B_formula.dot(r_y) - (zx * zy * fy - (1 + zy**2) * fx) / zp)
check("b_psi", B_formula.dot(r_p) - (zx * fy - zy * fx))

# ---------------------------------------------------------------------------
# Check 2: arbitrary z(x,y,psi), phi(x,y,psi)
# ---------------------------------------------------------------------------
print("-- Check 2: generic-function identities (arbitrary z, phi)")

X, Y, P = sp.symbols("x y psi")
z = sp.Function("z")(X, Y, P)
phi = sp.Function("phi")(X, Y, P)


def inverse_system(z, phi):
    """Return J, covariant components, and Jacobian-scaled curl components."""
    zx, zy, zp = sp.diff(z, X), sp.diff(z, Y), sp.diff(z, P)
    fx, fy = sp.diff(phi, X), sp.diff(phi, Y)
    J = zp
    b_x = ((1 + zx**2) * fy - zx * zy * fx) / zp
    b_y = (zx * zy * fy - (1 + zy**2) * fx) / zp
    b_p = zx * fy - zy * fx
    C_x = sp.diff(b_p, Y) - sp.diff(b_y, P)   # C^i = J (curl B . grad u^i)
    C_y = sp.diff(b_x, P) - sp.diff(b_p, X)
    C_p = sp.diff(b_y, X) - sp.diff(b_x, Y)
    return J, b_x, b_y, b_p, C_x, C_y, C_p


J, b_x, b_y, b_p, C_x, C_y, C_p = inverse_system(z, phi)

# 2a: differential identity (div curl = 0) -> only two of F1-F3 independent
check("d(C_x)/dx + d(C_y)/dy + d(C_psi)/dpsi = 0",
      sp.diff(C_x, X) + sp.diff(C_y, Y) + sp.diff(C_p, P))

# 2b: Euler-Lagrange equations of the transformed energy functional
D = sp.diff(z, X) * sp.diff(phi, Y) - sp.diff(z, Y) * sp.diff(phi, X)
L = (sp.diff(phi, X) ** 2 + sp.diff(phi, Y) ** 2 + D**2) / sp.diff(z, P)


def euler_lagrange(L, w):
    return sum(sp.diff(sp.diff(L, sp.diff(w, v)), v) for v in (X, Y, P))


EL_phi = euler_lagrange(L, phi)
EL_z = euler_lagrange(L, z)
check("EL_phi = -2 C_psi", EL_phi + 2 * C_p)
check("EL_z = -2 (C_x phi_x + C_y phi_y) / z_psi",
      EL_z + 2 * (C_x * sp.diff(phi, X) + C_y * sp.diff(phi, Y)) / sp.diff(z, P))

# ---------------------------------------------------------------------------
# Check 3a: uniform transverse field B = (B0, 0, 0)
#   phi = B0 y, psi = z  ->  z = psi, phi = B0 y
# ---------------------------------------------------------------------------
print("-- Check 3a: uniform field z = psi, phi = B0 y")

B0 = sp.symbols("B_0", positive=True)
_, _, _, b_p_u, C_x_u, C_y_u, C_p_u = inverse_system(P, B0 * Y)
check("C_x = 0", C_x_u)
check("C_y = 0", C_y_u)
check("C_psi = 0", C_p_u)

# ---------------------------------------------------------------------------
# Check 3b: hyperbolic field B = (y, x, B0)
#   psi = (B0/2)(y - x) exp(z/B0)  ->  z = B0 log( 2 psi / (B0 (y - x)) )
#   phi = (x + y) exp(-z/B0)       ->  phi = B0 (y^2 - x^2) / (2 psi)
# ---------------------------------------------------------------------------
print("-- Check 3b: hyperbolic field B = (y, x, B0)")

z_e = B0 * sp.log(2 * P / (B0 * (Y - X)))
phi_e = B0 * (Y**2 - X**2) / (2 * P)

# the Cartesian field reconstructed from the inverse unknowns
zx_e, zy_e, zp_e = sp.diff(z_e, X), sp.diff(z_e, Y), sp.diff(z_e, P)
fx_e, fy_e = sp.diff(phi_e, X), sp.diff(phi_e, Y)
B_rec = sp.Matrix([fy_e, -fx_e, zx_e * fy_e - zy_e * fx_e]) / zp_e
check_vec("B reconstruction - (y, x, B0)", B_rec - sp.Matrix([Y, X, B0]))

_, _, _, _, C_x_e, C_y_e, C_p_e = inverse_system(z_e, phi_e)
check("(F2) C_x = 0", C_x_e)
check("(F3) C_y = 0", C_y_e)
check("(F1) C_psi = 0", C_p_e)

print("All checks passed.")
