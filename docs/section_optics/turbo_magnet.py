"""TURBO-regime curved combined-function iron electromagnet.

Parameter regime taken from the PUBLIC paper

  A. F. Steinberg, R. B. Appleby, J. S. L. Yap, S. L. Sheehy,
  "Design of a large energy acceptance beamline using fixed field
  accelerator optics", Phys. Rev. Accel. Beams 27, 071601 (2024),
  arXiv:2402.01120  (open access)

Section IV of that paper specifies, for the TURBO demonstrator arc:
100 mm magnets, drifts of 75 / 112.5 mm, a 30 deg total bend over four
cells of three magnets, protons of 0.5-3.0 MeV (rigidity 0.1-0.25 T m),
momentum acceptance +-42 %, a working-point scan over field index
k = 40-120 with B0F = 0.30-0.60 T, and an Enge fringe with lambda = 7 mm.
Its own magnets are permanent-magnet Halbach arrays; those are NOT used
here.  This example is a conventional iron-dominated electromagnet
because the design variable of interest is the POLE SHAPE: the coil sets
the excitation, the pole contour sets the field profile.

Every number below is either from that paper (cited inline) or an
explicit modelling choice stated as such -- nothing comes from any
non-public source.

Why this regime and not the earlier C-magnet toy:

  quantity                  C-magnet toy    here        ratio
  frame curvature h         0.123 /m        2.0 /m      16x
  h*x at x = 10 mm          0.0012          0.020       16x
  on-orbit gradient         0.56 T/m        15 T/m      27x
  fringe / magnet length    small           ~30 %       --
  momentum spread           single orbit    +-42 %      multi-orbit

so the non-paraxial exact-sqrt formulation, p_s >= 2, the graded fringe
subdivision and the multi-orbit battery all become load-bearing rather
than decorative.
"""
import numpy as np
import ngsolve as ng
from netgen.occ import Axis, Face, OCCGeometry, Pnt, Segment, Vec, Wire, Z

# ---- regime constants (paper) ------------------------------------------
RIGIDITY_T_M = 0.25          # top of the 0.1-0.25 T m range, Sec. IV
ARC_LENGTH_M = 0.100         # "each magnet is 10 cm long", Sec. IV
FRINGE_LAMBDA_M = 0.007      # Enge lambda = 7 mm, Sec. III B
# ---- modelling choices (stated, not from the paper) --------------------
B0_T = 0.5                   # inside the 0.30-0.60 T scan window
FIELD_INDEX = 15.0           # below the paper's k = 40-120 scan so the
#                              pole taper stays machinable; the regime
#                              (strong combined function) is preserved
HALF_GAP_M = 0.010           # 20 mm full gap -> Enge lambda ~ 7 mm
POLE_INNER_M = 0.480         # pole spans r0 -+ 20 mm
POLE_OUTER_M = 0.520
YOKE_INNER_M = 0.450
LEG_INNER_M = 0.560
LEG_OUTER_M = 0.590
POLE_ROOT_Z_M = 0.060        # underside of the top plate
YOKE_TOP_Z_M = 0.090
COIL_Z_M = 0.035
COIL_INNER_M = 0.470
COIL_OUTER_M = 0.530
CONTOUR_POINTS = 24

REFERENCE_RADIUS_M = RIGIDITY_T_M / B0_T          # 0.5 m
CURVATURE_PER_M = 1.0 / REFERENCE_RADIUS_M        # 2.0 /m
SECTOR_ANGLE_RAD = ARC_LENGTH_M / REFERENCE_RADIUS_M
GRADIENT_T_PER_M = FIELD_INDEX * B0_T / REFERENCE_RADIUS_M
AMPERE_TURNS = 2.0 * HALF_GAP_M * B0_T / (4.0e-7 * np.pi)


def half_gap_at(radius):
    """Scaling-FFA gap law: B ~ (r/r0)^k  =>  g ~ (r/r0)^-k."""
    return HALF_GAP_M * (np.asarray(radius, dtype=float)
                         / REFERENCE_RADIUS_M) ** (-FIELD_INDEX)


def build_upper_half_yoke(sector_angle=None):
    """Upper (z >= 0) half of the curved C-yoke, revolved about +z."""
    angle = SECTOR_ANGLE_RAD if sector_angle is None else float(sector_angle)
    radii = np.linspace(POLE_INNER_M, POLE_OUTER_M, CONTOUR_POINTS)
    contour = [(float(r), float(z)) for r, z in zip(radii, half_gap_at(radii))]
    profile = ([(YOKE_INNER_M, POLE_ROOT_Z_M)]
               + [(POLE_INNER_M, POLE_ROOT_Z_M)]
               + contour
               + [(POLE_OUTER_M, POLE_ROOT_Z_M),
                  (LEG_INNER_M, POLE_ROOT_Z_M),
                  (LEG_INNER_M, 0.0),
                  (LEG_OUTER_M, 0.0),
                  (LEG_OUTER_M, YOKE_TOP_Z_M),
                  (YOKE_INNER_M, YOKE_TOP_Z_M)])
    segments = [Segment(Pnt(profile[i][0], 0.0, profile[i][1]),
                        Pnt(profile[(i + 1) % len(profile)][0], 0.0,
                            profile[(i + 1) % len(profile)][1]))
                for i in range(len(profile))]
    face = Face(Wire(segments))
    solid = face.Revolve(Axis(Pnt(0, 0, 0), Z), np.degrees(angle))
    # Place the sector so the design orbit passes through the origin
    # heading +x.  The native tracker stops on a plane x = const crossed
    # upward, which needs x to be monotone along the orbit; a sector
    # centred on the +x axis instead turns back and would trip the exit
    # test at the magnet's midpoint.
    solid = solid.Rotate(Axis(Pnt(0, 0, 0), Z),
                         np.degrees(-0.5 * np.pi - 0.5 * angle))
    solid = solid.Move(Vec(0.0, REFERENCE_RADIUS_M, 0.0))
    return solid.mat("iron")


def coil_filaments(turns_points=48):
    """Two pancake loops (z = +-COIL_Z_M) encircling the pole.

    The pair is symmetric under z -> -z with the same circulation, so
    B_z is even and B_x, B_y are odd -- the parity the z-mirror image
    contract requires.
    """
    angles = np.linspace(-0.5 * SECTOR_ANGLE_RAD - 0.10,
                         0.5 * SECTOR_ANGLE_RAD + 0.10, turns_points)
    loops = []
    shift = -0.5 * np.pi
    for z_value in (COIL_Z_M, -COIL_Z_M):
        outward = [[COIL_INNER_M * np.cos(a + shift),
                    COIL_INNER_M * np.sin(a + shift) + REFERENCE_RADIUS_M,
                    z_value] for a in angles]
        back = [[COIL_OUTER_M * np.cos(a + shift),
                 COIL_OUTER_M * np.sin(a + shift) + REFERENCE_RADIUS_M,
                 z_value] for a in angles[::-1]]
        loop = outward + back
        loop.append(loop[0])
        loops.append(np.asarray(loop))
    return loops


def orbit_entrance():
    """Entry point and heading of the reference orbit.

    The magnet is centred on azimuth 0 and the design orbit is the circle
    r = r0, so the particle enters at the sector's start face heading
    tangentially.  The entrance is backed off by one gap so the tracker
    starts in the fringe-free region.
    """
    # The design circle is followed only INSIDE the magnet; the entrance
    # drift is field free and the particle travels STRAIGHT along it.  So
    # the start must sit on the TANGENT LINE at the magnet's entrance
    # face, not on the circle itself: starting on the circle with the
    # local tangent leaves the particle d^2/(2R) off the design orbit by
    # the time it arrives (35 mm for a 190 mm drift here, which threw the
    # first tracked orbit clean off the 40 mm wide pole).
    drift = 0.060
    phi_in = -0.5 * np.pi - 0.5 * SECTOR_ANGLE_RAD
    entrance = np.array([REFERENCE_RADIUS_M * np.cos(phi_in),
                         REFERENCE_RADIUS_M * np.sin(phi_in)
                         + REFERENCE_RADIUS_M, 0.0])
    heading = np.array([-np.sin(phi_in), np.cos(phi_in), 0.0])
    point = entrance - drift * heading
    phi_out = -0.5 * np.pi + 0.5 * SECTOR_ANGLE_RAD
    exit_x = REFERENCE_RADIUS_M * np.cos(phi_out) + 0.5 * drift
    return point, heading, float(exit_x)


if __name__ == "__main__":
    import time

    print("TURBO-regime curved combined-function electromagnet")
    print(f"  reference radius   {REFERENCE_RADIUS_M*1e3:8.1f} mm "
          f"(curvature {CURVATURE_PER_M:.2f} /m)")
    print(f"  sector angle       {np.degrees(SECTOR_ANGLE_RAD):8.2f} deg "
          f"(arc {ARC_LENGTH_M*1e3:.0f} mm)")
    print(f"  central field      {B0_T:8.2f} T   rigidity "
          f"{RIGIDITY_T_M:.3f} T m")
    print(f"  field index k      {FIELD_INDEX:8.1f}     gradient "
          f"{GRADIENT_T_PER_M:.1f} T/m")
    print(f"  half gap at r0     {HALF_GAP_M*1e3:8.1f} mm  -> at pole edges "
          f"{half_gap_at(POLE_INNER_M)*1e3:.1f} / "
          f"{half_gap_at(POLE_OUTER_M)*1e3:.1f} mm")
    print(f"  h*x at x = 10 mm   {CURVATURE_PER_M*0.010:8.3f}  "
          f"(C-magnet toy: 0.001)")
    print(f"  ampere-turns       {AMPERE_TURNS:8.0f} A per gap")
    started = time.perf_counter()
    with ng.TaskManager():
        solid = build_upper_half_yoke()
        mesh = ng.Mesh(OCCGeometry(solid).GenerateMesh(maxh=0.012))
        volumes = np.asarray(ng.Integrate(ng.CoefficientFunction(1.0), mesh,
                                          element_wise=True))
    print(f"\nmesh {mesh.ne} elements in {time.perf_counter()-started:.1f} s; "
          f"iron volume {volumes.sum()*1e6:.0f} cm^3")
    loops = coil_filaments()
    print(f"coil: {len(loops)} pancake loops, "
          f"{loops[0].shape[0]} points each, "
          f"{AMPERE_TURNS/2:.0f} A-turns per loop")
    entry, heading, exit_x = orbit_entrance()
    print(f"orbit entry {np.round(1e3*entry, 1)} mm heading "
          f"{np.round(heading, 3)}; exit plane x = {1e3*exit_x:.1f} mm "
          f"(the orbit passes the origin heading +x)")
    centre = np.asarray(ng.Integrate(
        ng.CoefficientFunction((ng.x, ng.y, ng.z)), mesh, element_wise=False))
    print(f"iron centroid {np.round(1e3*centre/volumes.sum(), 1)} mm")
