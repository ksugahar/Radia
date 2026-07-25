"""Chi-based (Legendre potential) saturable hodograph mode-sum solver validation.

The Legendre potential chi = H.r - Psi turns the hodograph coordinates into a
single scalar unknown whose FIRST derivatives are the physical coordinates.
This driver validates the numerical machinery (radial mode solver, B-radial
material input, Robin/MMF boundary, coordinate recovery, orientation monitor)
against closed-form solutions, rung by rung:
  A. radial mode solver vs closed forms (power-law mu=q^k, q-radial form)
  B. radial mode solver vs closed forms (Froehlich mu_s(B), B-radial form)
  C. 2D mode-sum assembly vs exact multi-mode solution (chi, xi, eta, J)
  D. Robin (MMF-type q*chi_q - chi = g) boundary, manufactured solution
  E. design-flavored demo: Froehlich material, constant-Psi (MMF) pole face,
     outputs the designed face contour + J monitoring (synthetic entrance data)

Equation (q-radial):  d/dq( mu q chi_q ) + ((mu q)'/q) chi_thth = 0
Equation (B-radial):  d/dB( mu_d B chi_B ) + (mu_s/B) chi_thth = 0
Coordinates:          e_H.r = chi_q = mu_d chi_B,  e_perp.r = chi_th/q = (mu_s/B) chi_th
Orientation:          J = det d(x,y)/d(q,th) < 0 everywhere (no folding when |J|>0)

Rung E is a MACHINERY demo with synthetic entrance data, not a physical design
scenario; the physically grounded end-to-end design check lives in
verify_chaplygin_bend_design.py.

Golden bands: closed-form agreement < 1e-5 at M=401, FD convergence order
>= 1.8, Robin residual < 1e-10, J single-signed.  Run:
    python verify_chi_modesum_solver.py
Writes results_chi_modesum_solver.json next to this file (committed).
"""
import os
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json
import math
import numpy as np

# ---------------------------------------------------------------- materials
K_POWER = -0.5  # power-law exponent, -1 < k <= 0 (saturating)

def mu_power(q):
    return q ** K_POWER

A_FR = 500.0   # Froehlich initial slope (dimensionless)
BS_FR = 2.0    # Froehlich saturation flux density

def mu_s_fr(B):
    return A_FR * (1.0 - B / BS_FR)

def mu_d_fr(B):
    return A_FR * (1.0 - B / BS_FR) ** 2

def q_of_B_fr(B):
    return B / mu_s_fr(B)

# ---------------------------------------------------------------- radial solver
def solve_radial(x, acoef, bcoef, nu, bc_left, bc_right):
    """Solve (a(x) R')' - nu^2 b(x) R = 0 with conservative 2nd-order FD.

    bc_* = ("D", value)                      Dirichlet R = value
         = ("R", (alpha, beta, gamma))       alpha*R' + beta*R = gamma (right end only)
    """
    M = len(x)
    h = x[1] - x[0]
    xm = 0.5 * (x[:-1] + x[1:])
    am = acoef(xm)
    A = np.zeros((M, M))
    rhs = np.zeros(M)
    for i in range(1, M - 1):
        A[i, i - 1] = am[i - 1] / h**2
        A[i, i] = -(am[i - 1] + am[i]) / h**2 - nu**2 * bcoef(x[i])
        A[i, i + 1] = am[i] / h**2
    kind, val = bc_left
    if kind != "D":
        raise ValueError("left end supports Dirichlet only in this prototype")
    A[0, 0] = 1.0
    rhs[0] = val
    kind, val = bc_right
    if kind == "D":
        A[-1, -1] = 1.0
        rhs[-1] = val
    elif kind == "R":
        alpha, beta, gamma = val
        # one-sided 2nd-order derivative at the right end
        A[-1, -1] = alpha * 3.0 / (2 * h) + beta
        A[-1, -2] = alpha * (-4.0) / (2 * h)
        A[-1, -3] = alpha * 1.0 / (2 * h)
        rhs[-1] = gamma
    else:
        raise ValueError(kind)
    return np.linalg.solve(A, rhs)

def dx_center(f, h):
    d = np.empty_like(f)
    d[1:-1] = (f[2:] - f[:-2]) / (2 * h)
    d[0] = (-3 * f[0] + 4 * f[1] - f[2]) / (2 * h)
    d[-1] = (3 * f[-1] - 4 * f[-2] + f[-3]) / (2 * h)
    return d

# ---------------------------------------------------------------- rung A
def s_exponent(nu, k, branch):
    disc = math.sqrt(k * k + 4 * (k + 1) * nu * nu)
    return 0.5 * (-k + disc) if branch > 0 else 0.5 * (-k - disc)

def rung_A(report):
    k = K_POWER
    qa, qb = 0.6, 1.8
    acoef = lambda q: mu_power(q) * q          # mu q
    bcoef = lambda q: (k + 1) * q ** (k - 1)   # (mu q)'/q for mu=q^k
    cases = []
    for nu in (0.0, 1.0, 2.0, 4.0):
        if nu == 0.0:
            exact = lambda q: 1.3 + 0.6 * q ** (-k)          # s=0 and s=-k
        else:
            sp = s_exponent(nu, k, +1)
            sm = s_exponent(nu, k, -1)
            exact = lambda q, sp=sp, sm=sm: 0.8 * q ** sp + 0.5 * q ** sm
        errs = {}
        for M in (101, 201, 401):
            x = np.linspace(qa, qb, M)
            R = solve_radial(x, acoef, bcoef, nu, ("D", exact(qa)), ("D", exact(qb)))
            errs[M] = float(np.max(np.abs(R - exact(x))) / np.max(np.abs(exact(x))))
        order = math.log2(errs[201] / errs[401])
        cases.append({"nu": nu, "rel_err": errs, "order_201_401": round(order, 2)})
        print(f"  [A] nu={nu:3.1f}  rel_err(M=401)={errs[401]:.3e}  conv_order={order:.2f}")
    report["rung_A_power_law_radial"] = cases

# ---------------------------------------------------------------- rung B
def rung_B(report):
    Ba, Bb = 0.3, 1.6
    acoef = lambda B: mu_d_fr(B) * B
    bcoef = lambda B: mu_s_fr(B) / B
    def exact0(B):  # (mu_d B R')' = 0 : R = c1 + c2 * F(B), F' = 1/(mu_d B)
        u = B / BS_FR
        F = (np.log(u) - np.log(1.0 - u) + 1.0 / (1.0 - u)) / A_FR
        return 1.1 + 0.9 * F
    exact1 = lambda B: q_of_B_fr(B)            # translation mode chi = q(B) cos(th)
    cases = []
    for nu, exact, label in ((0.0, exact0, "n=0 closed form"), (1.0, exact1, "n=1 translation q(B)")):
        errs = {}
        for M in (101, 201, 401):
            x = np.linspace(Ba, Bb, M)
            R = solve_radial(x, acoef, bcoef, nu, ("D", exact(x[0])), ("D", exact(x[-1])))
            errs[M] = float(np.max(np.abs(R - exact(x))) / np.max(np.abs(exact(x))))
        order = math.log2(errs[201] / errs[401])
        cases.append({"nu": nu, "check": label, "rel_err": errs, "order_201_401": round(order, 2)})
        print(f"  [B] nu={nu:3.1f} ({label})  rel_err(M=401)={errs[401]:.3e}  conv_order={order:.2f}")
    report["rung_B_froehlich_radial"] = cases

# ---------------------------------------------------------------- 2D assembly
class ChiSolution:
    """chi(x_r, th) = sum_n R_n(x_r) cos(nu_n th) on [xa,xb] x [0,Theta] (Neumann th-walls)."""
    def __init__(self, x, Theta, modes):
        self.x = x
        self.h = x[1] - x[0]
        self.Theta = Theta
        self.modes = modes  # list of (nu, R array)

    def eval_grid(self, th):
        X, T = np.meshgrid(self.x, th, indexing="ij")
        chi = np.zeros_like(X)
        chi_x = np.zeros_like(X)
        chi_t = np.zeros_like(X)
        for nu, R in self.modes:
            Rp = dx_center(R, self.h)
            c = np.cos(nu * th)[None, :]
            s = np.sin(nu * th)[None, :]
            chi += R[:, None] * c
            chi_x += Rp[:, None] * c
            chi_t += -nu * R[:, None] * s
        return chi, chi_x, chi_t

def coords_from_chi(x, th, chi_x, chi_t, xi_factor, eta_factor):
    """xi = xi_factor(x)*chi_x, eta = eta_factor(x)*chi_t; then rotate by th."""
    XI = xi_factor(x)[:, None] * chi_x
    ETA = eta_factor(x)[:, None] * chi_t
    C = np.cos(th)[None, :]
    S = np.sin(th)[None, :]
    return XI * C - ETA * S, XI * S + ETA * C, XI, ETA

def jacobian_grid(xx, yy, hx, ht):
    x_r = np.gradient(xx, hx, axis=0)
    x_t = np.gradient(xx, ht, axis=1)
    y_r = np.gradient(yy, hx, axis=0)
    y_t = np.gradient(yy, ht, axis=1)
    return x_r * y_t - x_t * y_r

def rung_C(report):
    k = K_POWER
    qa, qb, Theta = 0.6, 1.8, math.pi / 2
    M, NT = 401, 181
    x = np.linspace(qa, qb, M)
    th = np.linspace(0.0, Theta, NT)
    acoef = lambda q: mu_power(q) * q
    bcoef = lambda q: (k + 1) * q ** (k - 1)
    spec = [(0.0, (1.3, 0.0, 0.6, -k)),  # c1*q^0 + c2*q^{-k}
            (2.0, None), (4.0, None)]
    modes, exact_modes = [], []
    for nu, z in spec:
        if nu == 0.0:
            c1, s1, c2, s2 = z
            exact = lambda q, c1=c1, s1=s1, c2=c2, s2=s2: c1 * q ** s1 + c2 * q ** s2
            exactp = lambda q, c1=c1, s1=s1, c2=c2, s2=s2: c1 * s1 * q ** (s1 - 1) + c2 * s2 * q ** (s2 - 1)
        else:
            sp, sm = s_exponent(nu, k, +1), s_exponent(nu, k, -1)
            exact = lambda q, sp=sp, sm=sm: 0.8 * q ** sp + 0.5 * q ** sm
            exactp = lambda q, sp=sp, sm=sm: 0.8 * sp * q ** (sp - 1) + 0.5 * sm * q ** (sm - 1)
        R = solve_radial(x, acoef, bcoef, nu, ("D", exact(qa)), ("D", exact(qb)))
        modes.append((nu, R))
        exact_modes.append((nu, exact, exactp))
    sol = ChiSolution(x, Theta, modes)
    chi, chi_x, chi_t = sol.eval_grid(th)
    chi_ex = np.zeros_like(chi); chix_ex = np.zeros_like(chi); chit_ex = np.zeros_like(chi)
    for nu, exact, exactp in exact_modes:
        c = np.cos(nu * th)[None, :]; s = np.sin(nu * th)[None, :]
        chi_ex += exact(x)[:, None] * c
        chix_ex += exactp(x)[:, None] * c
        chit_ex += -nu * exact(x)[:, None] * s
    xi_f = lambda q: np.ones_like(q)
    eta_f = lambda q: 1.0 / q
    xx, yy, XI, ETA = coords_from_chi(x, th, chi_x, chi_t, xi_f, eta_f)
    xx_e, yy_e, _, _ = coords_from_chi(x, th, chix_ex, chit_ex, xi_f, eta_f)
    scale = np.max(np.abs(chi_ex))
    cscale = max(np.max(np.abs(xx_e)), np.max(np.abs(yy_e)))
    e_chi = float(np.max(np.abs(chi - chi_ex)) / scale)
    e_xy = float(max(np.max(np.abs(xx - xx_e)), np.max(np.abs(yy - yy_e))) / cscale)
    J = jacobian_grid(xx, yy, sol.h, th[1] - th[0])
    ok_orient = bool(np.all(J < 0.0))
    report["rung_C_2d_assembly"] = {
        "rel_err_chi": e_chi, "rel_err_coords": e_xy,
        "J_all_negative": ok_orient, "min_absJ": float(np.min(np.abs(J))),
    }
    print(f"  [C] chi rel_err={e_chi:.3e}  coord rel_err={e_xy:.3e}  J<0 everywhere={ok_orient}  min|J|={np.min(np.abs(J)):.3e}")

def rung_D(report):
    k = K_POWER
    qa, qb, Theta = 0.6, 1.8, math.pi / 2
    M = 401
    x = np.linspace(qa, qb, M)
    acoef = lambda q: mu_power(q) * q
    bcoef = lambda q: (k + 1) * q ** (k - 1)
    errs = []
    for nu in (0.0, 2.0, 4.0):
        if nu == 0.0:
            exact = lambda q: 1.3 + 0.6 * q ** (-k)
            exactp = lambda q: 0.6 * (-k) * q ** (-k - 1)
        else:
            sp, sm = s_exponent(nu, k, +1), s_exponent(nu, k, -1)
            exact = lambda q, sp=sp, sm=sm: 0.8 * q ** sp + 0.5 * q ** sm
            exactp = lambda q, sp=sp, sm=sm: 0.8 * sp * q ** (sp - 1) + 0.5 * sm * q ** (sm - 1)
        gamma = qb * exactp(qb) - exact(qb)          # MMF-type Robin data
        R = solve_radial(x, acoef, bcoef, nu, ("D", exact(qa)), ("R", (qb, -1.0, gamma)))
        errs.append(float(np.max(np.abs(R - exact(x))) / np.max(np.abs(exact(x)))))
        print(f"  [D] nu={nu:3.1f}  Robin(q chi_q - chi) manufactured  rel_err={errs[-1]:.3e}")
    report["rung_D_robin_manufactured"] = {"rel_err": errs}

def rung_E(report):
    """Design-flavored demo (SYNTHETIC entrance data, machinery demo only):
    Froehlich material, B in [B_gap, B_face], th in [0, pi/2].
    - B = B_face: MMF Robin  q(B) mu_d chi_B - chi = Psi_p (constant)
    - B = B_gap: Dirichlet from a translated-uniform entrance state
    - th walls: Neumann (cos series)
    Output: designed pole-face contour (x,y) on the Robin edge + J monitor."""
    Ba, Bb, Theta = 0.3, 1.6, math.pi / 2
    M, NT, NMODE = 401, 181, 24
    Psi_p = -2.0e-3
    x = np.linspace(Ba, Bb, M)
    th = np.linspace(0.0, Theta, NT)
    acoef = lambda B: mu_d_fr(B) * B
    bcoef = lambda B: mu_s_fr(B) / B
    # entrance data: chi(B_gap, th) = q(B_gap) * (x0 cos th + small shaping)
    x0 = 1.0e-3
    f_ent = q_of_B_fr(Ba) * (x0 * np.cos(th) + 0.15 * x0 * np.cos(2 * th))
    # cosine projection (nu_n = 2n on quarter symmetry would need even modes only;
    # keep general nu_n = n*pi/Theta = 2n here)
    modes = []
    thq = np.linspace(0.0, Theta, 4001)
    f_ent_q = q_of_B_fr(Ba) * (x0 * np.cos(thq) + 0.15 * x0 * np.cos(2 * thq))
    alpha_R = q_of_B_fr(Bb) * mu_d_fr(Bb)
    for n in range(NMODE):
        nu = n * math.pi / Theta
        w = np.cos(nu * thq)
        cn = np.trapezoid(f_ent_q * w, thq) * (1.0 if n == 0 else 2.0) / Theta
        gamma = Psi_p if n == 0 else 0.0
        R = solve_radial(x, acoef, bcoef, nu, ("D", cn), ("R", (alpha_R, -1.0, gamma)))
        modes.append((nu, R))
    sol = ChiSolution(x, Theta, modes)
    chi, chi_B, chi_t = sol.eval_grid(th)
    xi_f = lambda B: mu_d_fr(B)
    eta_f = lambda B: mu_s_fr(B) / B
    xx, yy, XI, ETA = coords_from_chi(x, th, chi_B, chi_t, xi_f, eta_f)
    J = jacobian_grid(xx, yy, sol.h, th[1] - th[0])
    ok_orient = bool(np.all(J < 0.0))
    # residual check of the Robin condition on the face
    face_res = alpha_R * chi_B[-1, :] - chi[-1, :]
    face_dev = float(np.max(np.abs(face_res - Psi_p)) / abs(Psi_p))
    contour = {"x_mm": (1e3 * xx[-1, :]).round(6).tolist(),
               "y_mm": (1e3 * yy[-1, :]).round(6).tolist()}
    report["rung_E_design_demo"] = {
        "material": {"model": "Froehlich", "a": A_FR, "Bs": BS_FR},
        "domain": {"B_gap": Ba, "B_face": Bb, "Theta": Theta},
        "Psi_p": Psi_p, "n_modes": NMODE,
        "face_robin_rel_dev": face_dev,
        "J_all_negative": ok_orient,
        "min_absJ": float(np.min(np.abs(J))),
        "face_contour_mm": contour,
        "note": "synthetic entrance data; machinery demo, not a physical design scenario",
    }
    print(f"  [E] Froehlich design demo: face Robin residual dev={face_dev:.2e}  "
          f"J<0 everywhere={ok_orient}  min|J|={np.min(np.abs(J)):.3e}")
    print(f"      face contour x range [mm]: {1e3*xx[-1,:].min():.4f} .. {1e3*xx[-1,:].max():.4f}")
    print(f"      face contour y range [mm]: {1e3*yy[-1,:].min():.4f} .. {1e3*yy[-1,:].max():.4f}")

def main():
    report = {}
    print("rung A: power-law radial solver vs closed forms")
    rung_A(report)
    print("rung B: Froehlich B-radial solver vs closed forms")
    rung_B(report)
    print("rung C: 2D assembly vs exact multi-mode solution")
    rung_C(report)
    print("rung D: MMF-type Robin boundary, manufactured solution")
    rung_D(report)
    print("rung E: design-flavored Froehlich demo")
    rung_E(report)

    # ---------------- golden bands (fail loud) ----------------
    failures = []
    for case in report["rung_A_power_law_radial"]:
        if case["rel_err"][401] > 1e-5:
            failures.append(f"rung A nu={case['nu']}: err {case['rel_err'][401]:.2e}")
        if case["order_201_401"] < 1.8:
            failures.append(f"rung A nu={case['nu']}: order {case['order_201_401']}")
    for case in report["rung_B_froehlich_radial"]:
        if case["rel_err"][401] > 1e-5:
            failures.append(f"rung B nu={case['nu']}: err {case['rel_err'][401]:.2e}")
        if case["order_201_401"] < 1.8:
            failures.append(f"rung B nu={case['nu']}: order {case['order_201_401']}")
    rc = report["rung_C_2d_assembly"]
    if rc["rel_err_chi"] > 1e-5:
        failures.append(f"rung C chi err {rc['rel_err_chi']:.2e}")
    if rc["rel_err_coords"] > 5e-4:
        failures.append(f"rung C coord err {rc['rel_err_coords']:.2e}")
    if not rc["J_all_negative"]:
        failures.append("rung C: J changed sign")
    for e in report["rung_D_robin_manufactured"]["rel_err"]:
        if e > 1e-5:
            failures.append(f"rung D Robin err {e:.2e}")
    re_ = report["rung_E_design_demo"]
    if re_["face_robin_rel_dev"] > 1e-10:
        failures.append(f"rung E Robin dev {re_['face_robin_rel_dev']:.2e}")
    if not re_["J_all_negative"]:
        failures.append("rung E: J changed sign")
    report["golden"] = {"passed": not failures, "failures": failures}

    import datetime
    import platform
    report["meta"] = {
        "generated_at_utc": datetime.datetime.now(datetime.timezone.utc)
        .isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "python_version": platform.python_version(),
        "purpose": "correctness validation only (no timing claims)",
    }
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results_chi_modesum_solver.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"results -> {out}")
    if failures:
        for f_ in failures:
            print("GOLDEN FAIL:", f_)
        raise SystemExit(1)
    print("all golden bands passed")

if __name__ == "__main__":
    main()
