"""radia.planar_charges -- SHARED 2D planar exterior field + Maxwell torque.

Method-agnostic postprocessing for planar (per-unit-length) magnetostatics: given ANY 2D mesh and a
per-element magnetization ``M_elem`` (nEl, 2), evaluate the exterior field and the Maxwell-stress
torque.  Used identically by the collocation MMMM (``radia.mmmm2d``) and the HDiv-VIM
(``radia.vim._vim2d``) -- both produce a per-element M, so the exterior field is the field of that
magnetization = its M.n equivalent bound charge on the element edges (log kernel -ln(r)/(2 pi)).

The hot loops (charge-cloud field, circle torque integral) are C++
(radia._radia_pybind.PlanarChargeField / PlanarMaxwellTorqueCircle); this module builds the M.n
cloud and calls them.

    from radia.planar_charges import maxwell_torque, exterior_field
    T = maxwell_torque(mesh, M_elem, Rc=0.05, H_ext=(H0x, H0y))   # reluctance torque, unit length
    H = exterior_field(mesh, M_elem, P)                            # H at points P (n,2), in air
"""
from __future__ import annotations

import numpy as np
import ngsolve as ng

import radia._radia_pybind as _rp

MU0 = 4e-7 * np.pi


def _gauss01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w          # nodes on [0,1], weights sum to 1


def mn_edge_cloud(mesh, M_elem, ngauss=4):
    """Build the M.n equivalent bound-charge cloud of a per-element magnetization.

    For each element (CCW), each edge carries a uniform line charge density lambda = M_elem . n_out;
    sampled at ``ngauss`` Gauss points -> point charges Q_a = lambda * L * w_a at X_a.  Internal
    edges get contributions from both adjacent elements (M_k - M_k').n, handled by superposition.

    Returns (Xq (nq,2), Q (nq,)).
    """
    if mesh.dim != 2:
        raise ValueError("planar_charges: mesh.dim must be 2 (got %d)" % mesh.dim)
    M_elem = np.asarray(M_elem)                    # keep dtype (real, or complex phasor for eddy)
    pts = np.array([list(mesh[v].point)[:2] for v in mesh.vertices])
    tg, wg = _gauss01(ngauss)
    Xs, Qs = [], []
    for el in mesh.Elements(ng.VOL):
        k = el.nr
        V = pts[[v.nr for v in el.vertices]]
        x, y = V[:, 0], V[:, 1]
        A2 = np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
        if A2 < 0:
            V = V[::-1]
        c = V.mean(axis=0)                    # vertex centroid is enough for the outward-normal test
        Mk = M_elem[k]
        nv = len(V)
        for i in range(nv):
            P0 = V[i]; P1 = V[(i + 1) % nv]
            t = P1 - P0; L = np.hypot(t[0], t[1]); th = t / L
            nout = np.array([th[1], -th[0]])
            mid = 0.5 * (P0 + P1)
            if nout @ (mid - c) < 0:
                nout = -nout
            lam = Mk @ nout                   # uniform line-charge density (real, or complex phasor)
            X = P0[None, :] + tg[:, None] * (P1 - P0)[None, :]     # (ngauss, 2)
            Xs.append(X)
            Qs.append(lam * L * wg)
    return np.vstack(Xs), np.concatenate(Qs)


def charge_field(Xq, Q, P):
    """H at points P (n,2) from a point-charge cloud (Xq (nq,2), Q (nq,)).  C++ hot loop."""
    P = np.ascontiguousarray(np.asarray(P, float).reshape(-1, 2))
    return _rp.PlanarChargeField(np.ascontiguousarray(Xq, float),
                                 np.ascontiguousarray(Q, float), P)


def exterior_field(mesh, M_elem, P, ngauss=4):
    """Exterior H (in air) at points P (n,2) of the magnetization M_elem on ``mesh``."""
    Xq, Q = mn_edge_cloud(mesh, M_elem, ngauss)
    return charge_field(Xq, Q, P)


def maxwell_torque_cloud(Xq, Q, Rc, H_ext=(0.0, 0.0), center=(0.0, 0.0), n=1440):
    """Maxwell torque per unit length on a circle Rc in air from a cloud + uniform applied H_ext."""
    return _rp.PlanarMaxwellTorqueCircle(np.ascontiguousarray(Xq, float),
                                         np.ascontiguousarray(Q, float),
                                         float(Rc), float(center[0]), float(center[1]), int(n),
                                         float(H_ext[0]), float(H_ext[1]))


def maxwell_torque(mesh, M_elem, Rc, H_ext=(0.0, 0.0), center=(0.0, 0.0), n=1440, ngauss=4):
    """Maxwell-stress torque per unit length about ``center`` on a circle of radius Rc in air, for
    the magnetization M_elem on ``mesh`` in a UNIFORM applied field H_ext.  The reluctance/alignment
    torque = mu0 A (M x H0); Rc must enclose the body and lie in air."""
    Xq, Q = mn_edge_cloud(mesh, M_elem, ngauss)
    return maxwell_torque_cloud(Xq, Q, Rc, H_ext=H_ext, center=center, n=n)


# ---- Maxwell-stress FORCE (maglev levitation / actuator; shared by MMMM and HDiv-VIM) ------------

def maxwell_force_cloud(Xq, Q, Rc, H_ext=(0.0, 0.0), center=(0.0, 0.0), n=1440):
    """Force per unit length (Fx,Fy) on a circle Rc in air from a cloud + uniform applied H_ext.
    A uniform field gives ~0 net force; CONCATENATE several bodies' clouds (with the circle
    enclosing ONE) for the inter-body force."""
    return np.asarray(_rp.PlanarMaxwellForceCircle(
        np.ascontiguousarray(Xq, float), np.ascontiguousarray(Q, float),
        float(Rc), float(center[0]), float(center[1]), int(n),
        float(H_ext[0]), float(H_ext[1])), float)


def maxwell_force(mesh, M_elem, Rc, H_ext=(0.0, 0.0), center=(0.0, 0.0), n=1440, ngauss=4):
    """Maxwell-stress force per unit length (Fx,Fy) on a circle Rc (air, enclosing the body) for the
    magnetization M_elem in a UNIFORM applied field.  NOTE a uniform field exerts ZERO net force
    (only torque) -- use force_between(...) for the maglev inter-body / field-gradient force."""
    Xq, Q = mn_edge_cloud(mesh, M_elem, ngauss)
    return maxwell_force_cloud(Xq, Q, Rc, H_ext=H_ext, center=center, n=n)


def force_between(bodies, Rc, center, n=1440, ngauss=4):
    """Inter-body force per unit length on the body enclosed by the circle (Rc, center) in air.
    ``bodies`` = list of (mesh, M_elem); ALL bodies' M.n clouds are summed (the circle must enclose
    exactly ONE of them).  The maglev / actuator force: F on body A in the field of body B."""
    Xs, Qs = [], []
    for mesh, M in bodies:
        X, Q = mn_edge_cloud(mesh, M, ngauss)
        Xs.append(X); Qs.append(Q)
    return maxwell_force_cloud(np.vstack(Xs), np.concatenate(Qs), Rc, center=center, n=n)


def maxwell_force_complex(mesh, M_elem, Rc, H_ext=(0.0, 0.0), center=(0.0, 0.0), n=1440, ngauss=4):
    """TIME-AVERAGED Maxwell force per unit length (Fx,Fy) for COMPLEX phasor magnetization M_elem
    and complex uniform H_ext (maglev levitation from induced eddy currents):
        <F_i> = mu0 Rc (2pi/n) sum 0.5 Re[ H_r conj(H_i) - 1/2 (H . conj(H)) n_i ],  H = H_body + H_ext."""
    Xq, Qc = mn_edge_cloud(mesh, M_elem, ngauss)
    phi = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    P = np.stack([center[0] + Rc * np.cos(phi), center[1] + Rc * np.sin(phi)], axis=1)
    H = field_complex(Xq, Qc, P) + np.asarray(H_ext, complex)[None, :]
    c, s = np.cos(phi), np.sin(phi)
    Hr = H[:, 0] * c + H[:, 1] * s
    H2 = H[:, 0] * np.conj(H[:, 0]) + H[:, 1] * np.conj(H[:, 1])
    Fx = np.real(Hr * np.conj(H[:, 0]) - 0.5 * H2 * c).sum()
    Fy = np.real(Hr * np.conj(H[:, 1]) - 0.5 * H2 * s).sum()
    return MU0 * Rc * (2 * np.pi / n) * 0.5 * np.array([Fx, Fy])


# ---- out-of-plane vector potential A_z (for the reduced-FEM eddy / maglev coupling) --------------

def charge_az(Xq, Q, P):
    """A_z at points P (n,2) from a point-charge cloud: A_z = mu0/(2pi) sum Q atan2(dy,dx).  C++."""
    P = np.ascontiguousarray(np.asarray(P, float).reshape(-1, 2))
    return _rp.PlanarChargeAz(np.ascontiguousarray(Xq, float), np.ascontiguousarray(Q, float), P)


def vector_potential_az(mesh, M_elem, P, ngauss=4):
    """Out-of-plane vector potential A_z (in air) at points P (n,2) of the magnetization M_elem.
    BRANCH-CUT CAVEAT (atan2): valid where the eval set sees the body from one side."""
    Xq, Q = mn_edge_cloud(mesh, M_elem, ngauss)
    return charge_az(Xq, Q, P)


# ---- COMPLEX (eddy-current phasor) field + TIME-AVERAGED torque (IM / maglev) --------------------

def field_complex(Xq, Qc, P):
    """Complex H at P from a COMPLEX point-charge cloud (phasors): the log-grad kernel is linear, so
    H = Field(Re Q) + 1j Field(Im Q).  Two real C++ calls."""
    Qc = np.asarray(Qc)
    Hr = charge_field(Xq, np.ascontiguousarray(Qc.real, float), P)
    Hi = charge_field(Xq, np.ascontiguousarray(Qc.imag, float), P)
    return Hr + 1j * Hi


def maxwell_torque_complex(mesh, M_elem, Rc, H_ext=(0.0, 0.0), center=(0.0, 0.0), n=1440, ngauss=4):
    """TIME-AVERAGED Maxwell torque per unit length for COMPLEX phasor magnetization M_elem (nEl,2)
    and complex uniform applied H_ext, on a circle Rc in air:
        <T> = mu0 Rc^2 (2 pi / n) sum_i (1/2) Re( H_r conj(H_phi) ),  H = H_body + H_ext.
    (The C++ real MaxwellTorqueCircle is the M_elem-real, H_ext-real special case.)"""
    Xq, Qc = mn_edge_cloud(mesh, M_elem, ngauss)                  # complex M -> complex Q
    phi = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    P = np.stack([center[0] + Rc * np.cos(phi), center[1] + Rc * np.sin(phi)], axis=1)
    H = field_complex(Xq, Qc, P) + np.asarray(H_ext, complex)[None, :]
    er = np.stack([np.cos(phi), np.sin(phi)], axis=1)
    et = np.stack([-np.sin(phi), np.cos(phi)], axis=1)
    Hr = (H * er).sum(axis=1)
    Ht = (H * et).sum(axis=1)
    return MU0 * Rc * Rc * (2 * np.pi / n) * 0.5 * float(np.real(Hr * np.conj(Ht)).sum())
