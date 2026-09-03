"""Lubich convolution-quadrature (CQ) time-domain acoustic BEM (sound-soft sphere).

The retarded single-layer time convolution V(d/dt) q = g is turned into N decoupled
Laplace-domain solves at s = delta(zeta)/dt (BDF2) on the rho-circle; each is a
Helmholtz single layer at COMPLEX wavenumber kappa = i s / c (ngsolve.bem accepts
complex kappa via the HelmholtzSL(integrand, kappa) form), and an FFT recovers the
time history.  Complement NGSolve: the per-frequency BEM is ngsolve.bem, the CQ
transform is the thin Python wrapper.

Validated (validation_test/acoustics): each CQ node's frequency-domain BEM scattered
field equals the analytic soft sphere at that COMPLEX wavenumber to ~3e-4
(soft_sphere_scattering_complex_k) -- the rigorous core check, independent of the
time-domain reconstruction; the recovered time signal is real and causal.

Observation points must lie in the x-z plane (y=0): the potential is evaluated on a
flat y=0 screen (the exterior of a general point cloud is not an ngsolve.bem target).
"""
import numpy as np
from scipy.special import eval_legendre, hankel1, jv


def bdf_delta(zeta, method="BDF2"):
    """BDF generating function delta(zeta) for the CQ Laplace nodes s = delta/dt."""
    values = np.asarray(zeta, dtype=complex)
    if not np.all(np.isfinite(values)):
        raise ValueError("zeta must be finite")
    if method == "BDF1":
        result = 1.0 - values
    elif method == "BDF2":
        result = 1.5 - 2.0 * values + 0.5 * values**2
    else:
        raise ValueError("method must be 'BDF1' or 'BDF2'")
    return result.item() if values.ndim == 0 else result


def _spherical_j(order, value):
    return np.sqrt(np.pi / (2.0 * value)) * jv(order + 0.5, value)


def _spherical_hankel(order, value):
    return np.sqrt(np.pi / (2.0 * value)) * hankel1(order + 0.5, value)


def soft_sphere_scattering_complex_k(k, radius, points, terms=28):
    """Sound-soft sphere partial-wave scattered field at a (possibly COMPLEX) k.

    This readable SciPy series is independent of the numerical BEM solve. It
    serves as the validation reference at CQ Laplace nodes ``kappa = i s / c``.
    """
    k = complex(k)
    radius = float(radius)
    terms = int(terms)
    if not np.isfinite(k.real) or not np.isfinite(k.imag) or abs(k) == 0.0:
        raise ValueError("wavenumber must be finite and nonzero")
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("radius must be positive and finite")
    if terms < 0 or terms > 512:
        raise ValueError("terms must lie between 0 and 512")
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    if len(pts) == 0 or not np.all(np.isfinite(pts)):
        raise ValueError("points must be a nonempty N-by-3 array of finite values")
    distance = np.linalg.norm(pts, axis=1)
    if np.any(distance < radius * (1.0 - 1.0e-9)):
        raise ValueError(
            "evaluation points must lie on or outside the sphere r >= R"
        )
    cosine = pts[:, 2] / distance
    scattered = np.zeros(len(pts), dtype=complex)
    for order in range(terms + 1):
        coefficient = (
            -(1j**order)
            * (2 * order + 1)
            * _spherical_j(order, k * radius)
            / _spherical_hankel(order, k * radius)
        )
        scattered += (
            coefficient
            * _spherical_hankel(order, k * distance)
            * eval_legendre(order, cosine)
        )
    return scattered


def _cq_grid(sample_count, time_step, sound_speed, method):
    sample_count = int(sample_count)
    time_step = float(time_step)
    sound_speed = float(sound_speed)
    if sample_count <= 0 or sample_count > 1_048_576:
        raise ValueError("sample_count must be a positive practical integer")
    if not np.isfinite(time_step) or time_step <= 0.0:
        raise ValueError("time_step must be positive and finite")
    if not np.isfinite(sound_speed) or sound_speed <= 0.0:
        raise ValueError("sound_speed must be positive and finite")
    index = np.arange(sample_count)
    radius = np.finfo(float).eps ** (1.0 / (2.0 * sample_count))
    zeta = radius * np.exp(-2j * np.pi * index / sample_count)
    nodes = bdf_delta(zeta, method) / time_step
    return {
        "backend": "numpy",
        "cq_radius": float(radius),
        "zeta": zeta,
        "cq_nodes": nodes,
        "cq_wavenumbers": 1j * nodes / sound_speed,
    }


def _build_sphere_screen(radius, maxh, order, screen_extent):
    from netgen.occ import WorkPlane, Axes, Sphere, Fuse, Compound, X, Y
    from ngsolve import SurfaceL2, Compress
    screen = WorkPlane(Axes((0, 0, 0), Y, X)).RectangleC(screen_extent, screen_extent).Face()
    sphere = Sphere((0, 0, 0), radius)
    screen = screen - sphere
    sp = Fuse(sphere.faces)
    screen.faces.name = "screen"
    sp.faces.name = "sphere"
    mesh = Compound([screen, sp]).GenerateMesh(maxh=maxh).Curve(order)
    fes = Compress(SurfaceL2(mesh, order=order, complex=True,
                             definedon=mesh.Boundaries("sphere")))
    return mesh, fes


def _frequency_scattered(mesh, fes, pre, kappa, ghat_cf, obs):
    """Solve the Laplace-domain single layer V(kappa) q = <ghat, v> and evaluate the
    single-layer potential S(kappa) q at obs (sound-soft direct single layer)."""
    from ngsolve import GridFunction, LinearForm, ds, solvers
    from ngsolve.bem import HelmholtzSL
    V = HelmholtzSL(fes.TrialFunction() * ds("sphere"), complex(kappa)) * fes.TestFunction() * ds
    rhs = LinearForm(ghat_cf * fes.TestFunction() * ds("sphere")).Assemble()
    q = GridFunction(fes)
    q.vec.data = solvers.GMRes(A=V.mat, b=rhs.vec, pre=pre, tol=1e-10, maxsteps=800, printrates=False)
    pot = HelmholtzSL(fes.TrialFunction() * ds("sphere"), complex(kappa))
    pot_cf = pot(q, mesh.Boundaries("screen"))
    return np.array([complex(pot_cf(mesh(px, py, pz))) for px, py, pz in obs])


def cq_soft_sphere_scattering(obs, radius=1.0, num_time=24, time_step=0.28,
                              sound_speed=1.0, pulse_center=None, pulse_width=None,
                              method="BDF2", order=3, maxh=0.4, screen_extent=7.0):
    """CQ time-domain sound-soft sphere scattering of an incident plane-wave pulse.

    Incident p_inc(x,t) = ricker((t - pulse_center - z/c)/pulse_width) travelling +z.
    ``obs`` (N x 3) must lie in the x-z plane.  Returns a dict with ``time``,
    ``scattered`` (num_time x nobs, real time-domain scattered pressure), the
    per-frequency ``pressure_hat`` / ``cq_wavenumbers`` (kappa) / ``pulse_transform``
    (A_l), and diagnostics (``max_imag`` of the pre-real signal).
    """
    from ngsolve import BilinearForm, TaskManager, ds, z as Zc, exp as ngexp

    c = float(sound_speed)
    R = float(radius)
    N = int(num_time)
    dt = float(time_step)
    obs = np.asarray(obs, float).reshape(-1, 3)
    t0 = 8 * dt if pulse_center is None else float(pulse_center)
    w = 2.5 * dt if pulse_width is None else float(pulse_width)

    n = np.arange(N)
    t = n * dt
    grid = _cq_grid(N, dt, c, str(method))
    rho = float(grid["cq_radius"])
    s = np.asarray(grid["cq_nodes"])
    kappa = np.asarray(grid["cq_wavenumbers"])

    finc = (1 - 2 * ((t - t0) / w) ** 2) * np.exp(-((t - t0) / w) ** 2)   # pulse at z=0
    A = np.fft.fft(rho**n * finc)                                        # CQ transform of the pulse

    mesh, fes = _build_sphere_screen(R, maxh, order, screen_extent)
    u, v = fes.TnT()
    pressure_hat = np.zeros((N, len(obs)), complex)
    with TaskManager():
        pre = BilinearForm(u * v * ds("sphere"), diagonal=True).Assemble().mat.Inverse()
        for l in range(N):
            ghat = -complex(A[l]) * ngexp(-complex(s[l]) / c * Zc)       # -A_l exp(i k_l z)
            pressure_hat[l, :] = _frequency_scattered(mesh, fes, pre, kappa[l], ghat, obs)

    pressure_complex = (rho ** (-n))[:, None] * np.fft.ifft(pressure_hat, axis=0)
    return {
        "kind": "cq_soft_sphere_time_domain_scattering",
        "time": t, "time_step": dt, "method": method, "cq_radius": rho,
        "scattered": pressure_complex.real,
        "cq_nodes": s, "cq_wavenumbers": kappa, "pulse_transform": A,
        "pressure_hat": pressure_hat,
        "max_imag": float(np.max(np.abs(pressure_complex.imag))),
        "obs": obs, "radius": R, "sound_speed": c,
    }
