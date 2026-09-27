"""TEAM-13-geometry model on a Kelvin open boundary (NGSolve + Radia Kelvin helpers).

Geometry (m): TEAM 13 description / Hanser TU Wien geometry.py (lab copy in
W:/00_CAE/NGSolve/矢野/2026_03_30_TEAM_benchmark/Problem13).
  coil   : rounded square loop, outer 200x200 (corner radius 50), inner 150x150
           (corner radius 25), corner centres (+-50,+-50) mm, z in [-50, 50] mm
  steel  : centre plate  x in [-1.6, 1.6], y in [-25, 25], z in [-63.2, 63.2] mm
           right channel x in [2.1, 125.3], y in [15, 65], web/flange 3.2 mm
           left channel  mirror (x -> -x, y -> -y)
B-H    : TEAM 13 table (lab notebook), energy form via BSpline integral.

The transient parameters (SIGMA, TAU, NI) are PLACEHOLDERS for the skeleton;
they are replaced by the TEAM 10 original values once available.
"""
import numpy as np
from netgen.occ import Box, Pnt, Sphere, Glue, OCCGeometry, X, Y, Z
from ngsolve import (Mesh, HCurl, H1, Periodic, BilinearForm, LinearForm, GridFunction, BSpline,
                     curl, grad, dx, sqrt, x, y, z, IfPos, CF, TaskManager, Integrate, InnerProduct)
from radia.kelvin_geometry import add_kelvin_exterior_domain
from radia.kelvin_material import make_kelvin_nu_cf, NU_0, MU_0

import os
R_K = 0.25
OFFSET = (0.6, 0.0, 0.0)
COIL_AREA = 0.025 * 0.100          # m^2
# TEAM_CASE=team13 (default): TEAM 13 geometry (channels shifted +-40 mm in y), sigma placeholder 7.5e6.
# TEAM_CASE=team10: channels aligned with the centre plate (y in [-25, 25] mm), sigma = 7.505e6 S/m
#   (Yan & Jin, PIER 153, 33-55, 2015; Yan et al., PIER 153, 69-91, 2015).  Excitation 162 turns,
#   5.64 A (1 - exp(-t/0.05 s)).  B-H: the TEAM 13 table below coincides with the measured points
#   plotted in Yan & Jin Fig. 2(a); not yet confirmed against the TEAM 10 original.
CASE = os.environ.get('TEAM_CASE', 'team13')
if CASE not in ('team13', 'team10'): raise ValueError(f'TEAM_CASE={CASE!r}')
SIGMA_PLACEHOLDER = 7.505e6 if CASE == 'team10' else 7.5e6   # S/m (name kept for existing callers)
CHANNEL_Y = (-.025, .025) if CASE == 'team10' else None       # None: TEAM 13 shifted channels
# TEAM_QUARTER=1 (TEAM 10 only): additionally cut at y = 0.  The coil current crosses y = 0
# normally, so A is normal there: tangential A = 0 (Dirichlet 'sym_y'); B_y = 0 on the plane.
QUARTER = os.environ.get('TEAM_QUARTER', '0') == '1'
if QUARTER and CASE != 'team10': raise ValueError('TEAM_QUARTER needs TEAM_CASE=team10 (TEAM 13 is not y-symmetric)')
# TEAM_COIL2=1 (TEAM 10 only): a second, independently driven rounded-rectangular coil around the web
# of the right channel (a changing excitation PATTERN).  Loop axis z, centre (123.7 mm, 0), inner
# 14 x 64 mm (corner radius 4 mm), 5 mm radial build, z in [25, 35] mm (mirrored below z = 0 by the
# symmetry); it crosses y = 0 normally and keeps both symmetry planes of the quarter model.
COIL2 = os.environ.get('TEAM_COIL2', '0') == '1'
if COIL2 and CASE != 'team10': raise ValueError('TEAM_COIL2 needs TEAM_CASE=team10')
COIL2_CENTER, COIL2_CORNER = (0.1237, 0.0), (0.003, 0.028)     # corner-centre offsets from the loop centre
COIL2_R_IN, COIL2_BUILD, COIL2_Z = 0.004, 0.005, (0.025, 0.035)
COIL2_AREA = COIL2_BUILD * (COIL2_Z[1] - COIL2_Z[0])

H_KL = [0.0, 16, 30, 54, 93, 143, 191, 210, 222, 233, 247, 258, 272, 289, 313, 342, 377, 433, 509, 648,
        933, 1228, 1934, 2913, 4993, 7189, 9423]
B_KL = [0.0, .0025, .005, .0125, .025, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0, 1.1, 1.2, 1.3, 1.4,
        1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8]
A_SAT, B_SAT, C_SAT, MS = -2.381e-10, 2.327e-5, 1.590, 2.16    # TEAM 13 Eq. (1)


def bh_table():
    """TEAM 13 B-H: table to 1.8 T, then Eq. (1): B = mu0 H + (aH^2+bH+c) up to 2.22 T,
    B = mu0 H + Ms beyond.  Returned as monotone (B, H) samples."""
    H = list(H_KL); B = list(B_KL)
    for h in np.linspace(9423, 48000, 40)[1:]:
        b = MU_0 * h + A_SAT * h * h + B_SAT * h + C_SAT
        if b >= 2.22: break
        H.append(float(h)); B.append(float(b))
    for h in np.geomspace(H[-1] * 1.2, 3e7, 40):
        H.append(float(h)); B.append(float(MU_0 * h + MS))
    order = np.argsort(B); B = np.array(B)[order]; H = np.array(H)[order]
    keep = np.r_[True, np.diff(B) > 1e-9]
    return B[keep], H[keep]


def build_mesh(maxh_steel=0.004, maxh_coil=0.015, maxh_air=0.03, curve=2):
    def zedges(s):
        return [e for e in s.edges if abs(e.start.x - e.end.x) < 1e-9 and abs(e.start.y - e.end.y) < 1e-9]
    outer = Box(Pnt(-.1, -.1, -.05), Pnt(.1, .1, .05)); outer = outer.MakeFillet(zedges(outer), .05)
    inner = Box(Pnt(-.075, -.075, -.05), Pnt(.075, .075, .05)); inner = inner.MakeFillet(zedges(inner), .025)
    coil = outer - inner; coil.mat('coil'); coil.maxh = maxh_coil
    t = .0032; hz = .0632
    centre = Box(Pnt(-.0016, -.025, -hz), Pnt(.0016, .025, hz))
    right = Box(Pnt(.0021, .015, -hz), Pnt(.0021 + .12 + t, .065, hz)) - Box(Pnt(0, .014, -hz + t), Pnt(.0021 + .12, .066, hz - t))
    left = Box(Pnt(-(.0021 + .12 + t), -.065, -hz), Pnt(-.0021, -.015, hz)) - Box(Pnt(-(.0021 + .12), -.066, -hz + t), Pnt(0, -.014, hz - t))
    steel = []
    for s, name in [(centre, 'steel_c'), (right, 'steel_r'), (left, 'steel_l')]:
        s.mat('steel'); s.maxh = maxh_steel
        for f in s.faces: f.name = 'steel_surf'
        steel.append(s)
    air = Sphere(Pnt(0, 0, 0), R_K)
    for f in air.faces: f.name = 'kelvin_int'
    air = air - coil - steel[0] - steel[1] - steel[2]; air.mat('air'); air.maxh = maxh_air
    geo, info = add_kelvin_exterior_domain([air, coil] + steel, OFFSET, R_K, inner_maxh=maxh_air, outer_maxh_factor=2.0)
    mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh_air * 2))
    if curve: mesh.Curve(curve)
    return mesh


def build_mesh_half(maxh_steel=0.004, maxh_coil=0.015, maxh_air=0.03, curve=2):
    """z >= 0 half model.  The Kelvin centre lies on z=0, so the inversion maps the
    upper half-ball onto the upper half-ball: the two hemispherical faces are
    identified (periodic, as in the full two-sphere convention) and every z=0 face
    is the symmetry plane 'sym_z' (B normal -> natural BC for the A-formulation)."""
    from netgen.occ import IdentificationType, Vertex
    def zedges(s):
        return [e for e in s.edges if abs(e.start.x - e.end.x) < 1e-9 and abs(e.start.y - e.end.y) < 1e-9]
    def upper(s, extent=2.0):
        box = Box(Pnt(-extent, 0 if QUARTER else -extent, 0), Pnt(extent, extent, extent))
        for f in box.faces:
            c = f.center
            f.name = 'sym_y' if (QUARTER and abs(c.y) < 1e-9) else 'sym_z'
        return s * box
    outer = Box(Pnt(-.1, -.1, -.05), Pnt(.1, .1, .05)); outer = outer.MakeFillet(zedges(outer), .05)
    inner = Box(Pnt(-.075, -.075, -.05), Pnt(.075, .075, .05)); inner = inner.MakeFillet(zedges(inner), .025)
    coil = upper(outer - inner); coil.mat('coil'); coil.maxh = maxh_coil
    t = .0032; hz = .0632
    (r0, r1), (l0, l1) = ((CHANNEL_Y, CHANNEL_Y) if CHANNEL_Y else ((.015, .065), (-.065, -.015)))
    parts = [Box(Pnt(-.0016, -.025, -hz), Pnt(.0016, .025, hz)),
             Box(Pnt(.0021, r0, -hz), Pnt(.0021 + .12 + t, r1, hz)) - Box(Pnt(0, r0 - .001, -hz + t), Pnt(.0021 + .12, r1 + .001, hz - t)),
             Box(Pnt(-(.0021 + .12 + t), l0, -hz), Pnt(-.0021, l1, hz)) - Box(Pnt(-(.0021 + .12), l0 - .001, -hz + t), Pnt(0, l1 + .001, hz - t))]
    steel = []
    for s in parts:
        for f in s.faces: f.name = 'steel_surf'
        s = upper(s); s.mat('steel'); s.maxh = maxh_steel; steel.append(s)
    extra = []
    if COIL2:
        (xc, yc), (ax, ay), ri, rb, (z0, z1) = COIL2_CENTER, COIL2_CORNER, COIL2_R_IN, COIL2_BUILD, COIL2_Z
        o2 = Box(Pnt(xc - ax - ri - rb, yc - ay - ri - rb, z0), Pnt(xc + ax + ri + rb, yc + ay + ri + rb, z1)); o2 = o2.MakeFillet(zedges(o2), ri + rb)
        i2 = Box(Pnt(xc - ax - ri, yc - ay - ri, z0), Pnt(xc + ax + ri, yc + ay + ri, z1)); i2 = i2.MakeFillet(zedges(i2), ri)
        coil2 = upper(o2 - i2); coil2.mat('coil2'); coil2.maxh = 0.004; extra.append(coil2)
    ball = Sphere(Pnt(0, 0, 0), R_K)
    for f in ball.faces: f.name = 'kelvin_int'
    air = upper(ball) - coil - steel[0] - steel[1] - steel[2]
    for c2 in extra: air = air - c2
    air.mat('air'); air.maxh = maxh_air
    kb = Sphere(Pnt(*OFFSET), R_K)
    for f in kb.faces: f.name = 'kelvin_ext'
    kel = upper(kb); kel.mat('kelvin'); kel.maxh = 2 * maxh_air
    fi = [f for f in air.faces if f.name == 'kelvin_int']; fo = [f for f in kel.faces if f.name == 'kelvin_ext']
    if len(fi) != 1 or len(fo) != 1: raise RuntimeError(('hemisphere faces', len(fi), len(fo)))
    fi[0].Identify(fo[0], 'kelvin_periodic', IdentificationType.PERIODIC)
    gnd = Vertex(Pnt(*OFFSET)); gnd.name = 'GND'
    geo = Glue([air, coil] + extra + steel + [kel, gnd])
    mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh_air * 2))
    if curve: mesh.Curve(curve)
    return mesh


def coil_direction(center=(0.0, 0.0), corner=(.05, .05)):
    """Unit tangent of a rounded-rectangle loop around the z axis through `center`, whose corner
    arcs are centred at center +- corner (clockwise seen from +z, as in the lab TEAM 13
    notebook: +x on the y>0 side)."""
    X0, Y0 = x - center[0], y - center[1]; ax, ay = corner
    cx = IfPos(X0 - ax, X0 - ax, IfPos(-ax - X0, X0 + ax, 0))
    cy = IfPos(Y0 - ay, Y0 - ay, IfPos(-ay - Y0, Y0 + ay, 0))
    r = sqrt(cx * cx + cy * cy + 1e-30)
    return CF((cy / r, -cx / r, 0))


def source_linear_form(fes, mesh, NI, order, region='coil', area=COIL_AREA, direction=None):
    """Coil current NI/area along the loop, Hodge-projected in the coil so that
    f(grad psi) = 0 for every H1 psi (discretely divergence-free source)."""
    J = NI / area * (coil_direction() if direction is None else direction)
    # quarter model: the current leaves the coil through y = 0, so psi = 0 there
    h1 = H1(mesh, order=order + 1, definedon=region, dirichlet='sym_y' if 'sym_y' in mesh.GetBoundaries() else '')
    p, q = h1.TnT()
    a = BilinearForm(grad(p) * grad(q) * dx(region) + 1e-8 * p * q * dx(region)).Assemble()
    rhs = LinearForm(J * grad(q) * dx(region)).Assemble()
    phi = GridFunction(h1); phi.vec.data = a.mat.Inverse(h1.FreeDofs(), inverse='sparsecholesky') * rhs.vec
    v = fes.TestFunction()
    f = LinearForm((J - grad(phi)) * v * dx(region, bonus_intorder=2)).Assemble()
    return f, J - grad(phi)


def coil2_source(fes, mesh, NI, order):
    return source_linear_form(fes, mesh, NI, order, region='coil2', area=COIL2_AREA,
                              direction=coil_direction(COIL2_CENTER, COIL2_CORNER))


def spaces(mesh, order=1):
    sym_y = 'sym_y' in mesh.GetBoundaries()
    return Periodic(HCurl(mesh, order=order, dirichlet='sym_y' if sym_y else '', dirichlet_bbnd='GND'))


def energy_density():
    Bt, Ht = bh_table()
    bs = BSpline(2, [0.0] + list(Bt), list(Ht))
    return bs.Integrate(), bs


# --- TEAM 10 (digitised by eye, +-0.03 T, from Yan & Jin PIER 153 (2015) Fig. 3(c)(d);
#     original data: Nakata, Takahashi, Fujiwara, COMPEL 14(2/3) 103-112, 1995) -------------
TEAM10_TURNS, TEAM10_IM, TEAM10_TAU = 162, 5.64, 0.05
TEAM10_MEASURED_DIGITISED = dict(
    t_ms=[10, 20, 30, 40, 50, 60, 80, 100, 120, 150],
    S1=[0.21, 0.52, 0.86, 1.13, 1.34, 1.49, 1.59, 1.62, 1.64, 1.65],
    S2=[0.10, 0.27, 0.43, 0.57, 0.69, 0.77, 0.84, 0.87, 0.89, 0.92],
    S3=[0.10, 0.28, 0.48, 0.66, 0.79, 0.89, 1.00, 1.05, 1.08, 1.11])


def team10_sections(B, mesh, n=(4, 12), s2_x=(0.0397, 0.0418)):
    """Average flux density over the TEAM 10 cross-sections (3.2 x 50 mm):
    S1 centre plate at mid-height (Bz), S3 channel web at mid-height (Bz),
    S2 upper flange (Bx) at x = 39.7 mm and at x = 2.1 + 39.7 mm (the reference of
    the 39.7 mm dimension in Yan & Jin Fig. 3(a) is ambiguous; both are reported)."""
    s = 1e-3; z0 = 1e-4
    # B_x, B_z are even in y, so the quarter model averages over y >= 0 only
    yi = (np.linspace(0.5, 24.5, n[1]) if QUARTER else np.linspace(-24.5, 24.5, n[1])) * s
    def avg(comp, pts): return float(np.mean([B(mesh(*p))[comp] for p in pts]))
    out = dict(S1=abs(avg(2, [(a, b, z0) for a in np.linspace(-1.5, 1.5, n[0]) * s for b in yi])),
               S3=abs(avg(2, [(a, b, z0) for a in np.linspace(122.2, 125.2, n[0]) * s for b in yi])))
    zi = np.linspace(60.1, 63.1, n[0]) * s
    for xx in s2_x:
        out[f'S2_x{xx*1000:.1f}'] = abs(avg(0, [(xx, b, c) for b in yi for c in zi]))
    return out


# --- TEAM 13 measured average flux densities (1000 AT), positions 1-25 ------
TEAM13_MEASURED_1000AT = [1.33, 1.329, 1.286, 1.225, 1.129, 0.985, 0.655,
                          0.259, 0.453, 0.554, 0.637, 0.698, 0.755, 0.809, 0.901, 0.945, 0.954, 0.956,
                          0.960, 0.965, 0.970, 0.974, 0.981, 0.984, 0.985]


def average_B_positions(B, mesh, n=(4, 12)):
    """Average flux-density component over the TEAM 13 cross-sections 1-25
    (same positions/components as the lab measurement.py, mm -> m)."""
    s = 1e-3; out = []
    def avg(comp, pts):
        vals = []
        for p in pts:
            vals.append(B(mesh(*p))[comp])
        return float(np.mean(vals))
    xi = np.linspace(.05, 1.55, n[0]) * s; yi = np.linspace(-24.5, 24.5, n[1]) * s
    for zz in [0, 10, 20, 30, 40, 50, 60]:
        zz = max(zz, .1)
        out.append(abs(avg(2, [(a, b, zz * s) for a in xi for b in yi])))
    yi = np.linspace(15.5, 64.5, n[1]) * s; zi = np.linspace(60.1, 63.1, n[0]) * s
    for xx in [2.1, 10, 20, 30, 40, 50, 60, 80, 100, 110, 122.1]:
        xx = max(xx, 2.2)
        out.append(abs(avg(0, [(xx * s, b, c) for b in yi for c in zi])))
    xi = np.linspace(122.2, 125.2, n[0]) * s
    for zz in [60, 50, 40, 30, 20, 10, 0]:
        zz = max(zz, .1)
        out.append(abs(avg(2, [(a, b, zz * s) for a in xi for b in yi])))
    return out
