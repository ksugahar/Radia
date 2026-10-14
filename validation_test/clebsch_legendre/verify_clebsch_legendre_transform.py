"""Symbolic verification of the Legendre (hodograph) transformation of
Clebsch-potential magnetostatics.

Forward formulation (independent variables x, y, z):

    B = grad(phi) x grad(psi)        (div B = 0 identically)
    curl B = 0                       (current-free region)

Transformed formulation (independent variables phi, psi, z; unknowns
x(phi,psi,z), y(phi,psi,z)):

    J   = x_phi * y_psi - x_psi * y_phi          (= 1 / B_z)
    b_phi = (x_phi*x_z + y_phi*y_z) / J          (covariant components of B)
    b_psi = (x_psi*x_z + y_psi*y_z) / J
    b_z   = (1 + x_z^2 + y_z^2) / J

    (E1)  d(b_psi)/d(phi) = d(b_phi)/d(psi)
    (E2)  d(b_z)/d(phi)   = d(b_phi)/d(z)
    (E3)  d(b_z)/d(psi)   = d(b_psi)/d(z)

Checks performed (all symbolic, sympy):

  1. Pointwise algebra: grad(phi) x grad(psi) = r_z / J in the inverse
     variables (B is tangent to the phi,psi = const lines, B_z = 1/J).
  2. For ARBITRARY x(phi,psi,z), y(phi,psi,z):
     a. the three curl components satisfy the differential identity
        d(C_phi)/d(phi) + d(C_psi)/d(psi) + d(C_z)/d(z) = 0
        (so only two of E1-E3 are independent);
     b. the Euler-Lagrange equations of the transformed energy functional
        W = (1/2 mu0) Int (1 + x_z^2 + y_z^2)/J dphi dpsi dz
        equal the force-balance projections
        EL_x = 2 (C_phi y_phi + C_psi y_psi) / J
        EL_y = -2 (C_phi x_phi + C_psi x_psi) / J
        i.e. the variational equations are ((curl B) x B = 0, force-free);
        curl B = 0 (vacuum) is the curl-free subset.
  3. Exact vacuum solutions satisfy E1-E3 identically:
     a. uniform field B = (0, 0, B0):    x = phi/B0, y = psi;
     b. hyperbolic field B = (y, x, B0):
        phi = (x+y) exp(-z/B0), psi = -(B0/2)(x-y) exp(z/B0),
        inverse: x = (phi e^{z/B0} - (2 psi/B0) e^{-z/B0}) / 2,
                 y = (phi e^{z/B0} + (2 psi/B0) e^{-z/B0}) / 2.

Run:  python verify_clebsch_legendre_transform.py
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
#   With r(phi,psi,z) = (x, y, z), the forward gradients are rows of the
#   inverse Jacobian; B = grad(phi) x grad(psi) must equal r_z / J.
# ---------------------------------------------------------------------------
print("-- Check 1: B = grad(phi) x grad(psi) = r_z / J (pointwise algebra)")

xf, xs, xz, yf, ys, yz = sp.symbols("x_phi x_psi x_z y_phi y_psi y_z")
A = sp.Matrix([[xf, xs, xz],
               [yf, ys, yz],
               [0, 0, 1]])          # A_ij = d x_i / d u_j, u = (phi, psi, z)
J1 = A.det()
check("J = x_phi y_psi - x_psi y_phi", J1 - (xf * ys - xs * yf))

Ainv = A.inv()                       # (A^-1)_ij = d u_i / d x_j
grad_phi = sp.Matrix(Ainv[0, :]).T   # row 1 -> column vector
grad_psi = sp.Matrix(Ainv[1, :]).T
B_clebsch = grad_phi.cross(grad_psi)
r_z = sp.Matrix([xz, yz, 1])
check_vec("grad(phi) x grad(psi) - r_z/J", B_clebsch - r_z / J1)

# Covariant components b_i = B . r_i reproduce the formulas used below
r_phi = sp.Matrix([xf, yf, 0])
r_psi = sp.Matrix([xs, ys, 0])
check("b_phi", (r_z / J1).dot(r_phi) - (xf * xz + yf * yz) / J1)
check("b_psi", (r_z / J1).dot(r_psi) - (xs * xz + ys * yz) / J1)
check("b_z", (r_z / J1).dot(r_z) - (1 + xz**2 + yz**2) / J1)

# ---------------------------------------------------------------------------
# Check 2: arbitrary x(phi,psi,z), y(phi,psi,z)
# ---------------------------------------------------------------------------
print("-- Check 2: generic-function identities (arbitrary x, y)")

f, s, z = sp.symbols("phi psi z")
x = sp.Function("x")(f, s, z)
y = sp.Function("y")(f, s, z)


def inverse_system(x, y):
    """Return J, covariant components, and curl components C^i * J."""
    J = sp.diff(x, f) * sp.diff(y, s) - sp.diff(x, s) * sp.diff(y, f)
    b_f = (sp.diff(x, f) * sp.diff(x, z) + sp.diff(y, f) * sp.diff(y, z)) / J
    b_s = (sp.diff(x, s) * sp.diff(x, z) + sp.diff(y, s) * sp.diff(y, z)) / J
    b_z = (1 + sp.diff(x, z) ** 2 + sp.diff(y, z) ** 2) / J
    C_f = sp.diff(b_z, s) - sp.diff(b_s, z)   # C^i = J (curl B . grad u^i),
    C_s = sp.diff(b_f, z) - sp.diff(b_z, f)   # the contravariant curl components
    C_z = sp.diff(b_s, f) - sp.diff(b_f, s)   # scaled by the Jacobian J
    return J, b_f, b_s, b_z, C_f, C_s, C_z


J, b_f, b_s, b_z, C_f, C_s, C_z = inverse_system(x, y)

# 2a: differential identity (div curl = 0) -> only two of E1-E3 independent
check("d(C_phi)/dphi + d(C_psi)/dpsi + d(C_z)/dz = 0",
      sp.diff(C_f, f) + sp.diff(C_s, s) + sp.diff(C_z, z))

# 2b: Euler-Lagrange equations of W = Int (1 + x_z^2 + y_z^2)/(2 mu0 J)
L = (1 + sp.diff(x, z) ** 2 + sp.diff(y, z) ** 2) / J


def euler_lagrange(L, w):
    return sum(sp.diff(sp.diff(L, sp.diff(w, v)), v) for v in (f, s, z))


EL_x = euler_lagrange(L, x)
EL_y = euler_lagrange(L, y)
check("EL_x = 2 (C_phi y_phi + C_psi y_psi) / J",
      EL_x - 2 * (C_f * sp.diff(y, f) + C_s * sp.diff(y, s)) / J)
check("EL_y = -2 (C_phi x_phi + C_psi x_psi) / J",
      EL_y + 2 * (C_f * sp.diff(x, f) + C_s * sp.diff(x, s)) / J)

# ---------------------------------------------------------------------------
# Check 3a: uniform field B = (0, 0, B0)
# ---------------------------------------------------------------------------
print("-- Check 3a: uniform field x = phi/B0, y = psi")

B0 = sp.symbols("B_0", positive=True)
_, _, _, b_z_u, C_f_u, C_s_u, C_z_u = inverse_system(f / B0, s)
check("C_phi = 0", C_f_u)
check("C_psi = 0", C_s_u)
check("C_z = 0", C_z_u)
check("b_z = B0", b_z_u - B0)

# ---------------------------------------------------------------------------
# Check 3b: hyperbolic field B = (y, x, B0)
# ---------------------------------------------------------------------------
print("-- Check 3b: hyperbolic field B = (y, x, B0)")

X, Y, Z = sp.symbols("x y z")
phi_e = (X + Y) * sp.exp(-Z / B0)
psi_e = -(B0 / 2) * (X - Y) * sp.exp(Z / B0)

grad = lambda g: sp.Matrix([sp.diff(g, X), sp.diff(g, Y), sp.diff(g, Z)])
B_forward = grad(phi_e).cross(grad(psi_e))
check_vec("forward: grad(phi) x grad(psi) - (y, x, B0)",
          B_forward - sp.Matrix([Y, X, B0]))
curlB = sp.Matrix([
    sp.diff(B_forward[2], Y) - sp.diff(B_forward[1], Z),
    sp.diff(B_forward[0], Z) - sp.diff(B_forward[2], X),
    sp.diff(B_forward[1], X) - sp.diff(B_forward[0], Y),
])
check_vec("forward: curl B = 0", curlB)

# inverse map
x_e = (f * sp.exp(z / B0) - (2 * s / B0) * sp.exp(-z / B0)) / 2
y_e = (f * sp.exp(z / B0) + (2 * s / B0) * sp.exp(-z / B0)) / 2
check("composition phi(x(phi,psi,z), y(...), z) = phi",
      phi_e.subs([(X, x_e), (Y, y_e), (Z, z)]) - f)
check("composition psi(x(phi,psi,z), y(...), z) = psi",
      psi_e.subs([(X, x_e), (Y, y_e), (Z, z)]) - s)

J_e, b_f_e, b_s_e, b_z_e, C_f_e, C_s_e, C_z_e = inverse_system(x_e, y_e)
check("(E3) C_phi = 0", C_f_e)
check("(E2) C_psi = 0", C_s_e)
check("(E1) C_z = 0", C_z_e)
check("B_z = 1/J reproduces B0", 1 / J_e - B0)

print("All checks passed.")
