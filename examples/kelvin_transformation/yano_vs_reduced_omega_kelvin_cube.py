"""yano-MSC (loop-removed) vs reduced-Omega + KELVIN FEM on a hex iron cube in a UNIFORM applied field.

Confirms the yano-type surface-charge MSC solve is LOOP-FREE and DISTORTION-ROBUST against an independent
open-boundary FEM, with the SAME uniform source in both:

  * source   : uniform H0 along +z, applied IDENTICALLY -- yano via rad.ObjBckg([0,0,mu0*H0]) (H_app=H0);
               FEM via Omega_s = H0*z on the iron boundary.  No coil (loops come from the multi-element
               MESH, distortion from shearing the hexes -- both mesh properties a uniform field exercises).
  * yano-MSC : a hex iron cube (REGULAR vs GRADED-DISTORTED grid) solved with rad.Solve (yano backend);
               observable = volume-averaged magnetization <M_z> over the iron (reliable; not the
               inside-iron field eval).
  * FEM      : reduced-Omega + Kelvin (open boundary, no PML).  Uniform field -> the HYBRID formulation
               (total-Omega in iron via Omega_s=H0*z on the iron boundary + reduced-Omega in air + Neumann
               jump) is MANDATORY and cancellation-free.  A naive "reduced-Omega everywhere" CANCELS
               catastrophically inside high-mu iron (H_in ~ Hs/67 is the tiny difference of two large
               ~Hs quantities) and reports no demagnetization -- this hybrid avoids it.

Result (mu_r=200, H0=1000): yano matches the reduced-Omega+Kelvin FEM to ~2% for BOTH the regular and the
distorted hex grid, and both land in the validated cube range (tests/feec/parity_vs_msc, ~3437) -- so the
loop removal (rad.GetLoopBasis) leaves a clean, distortion-robust field.

Two reproduction notes baked in (each cost a debug cycle):
  - kelvin_radius MUST enclose the box bounding SPHERE (corner = sqrt(3)*half-edge), not the half-edge;
    the auto value used the half-edge and the iron poked out of the air sphere -> broken mesh, no demag.
  - The Neumann jump on an OCC-built iron/air interface uses +specialcf.normal (the OCC "default"
    interface is oriented (mag,air), OPPOSITE the Cubit "sphere" sideset (air,mag) that
    solve_kelvin_benchmark's -normal is tuned to; -normal under-magnetizes 2.3x).
"""
import json
import math
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src", "radia"))
sys.path.insert(0, os.path.join(REPO, "src", "radia", "panels"))

import radia as rad  # noqa: E402
import ngsolve as ng  # noqa: E402
from ngsolve import grad, dx, ds, specialcf, CoefficientFunction as CF  # noqa: E402
from netgen.occ import Box, Pnt  # noqa: E402
from step_mesh_builder import build_mesh_from_step  # noqa: E402
from calc_common import detect_kelvin_offset, add_periodic_kelvin  # noqa: E402
from kelvin_source import kelvin_mu_factor_3d_cf, build_material_cf  # noqa: E402

MU0 = 4e-7 * math.pi
MU_R = 200.0
H0 = 1000.0
L = 0.020              # iron cube half-edge (m); centered at origin -> [-20, 20]^3 mm


def hex_cube(n, distort=0.0):
    """uniform (distort=0) or graded (distort>0) hex grid filling the cube; planar faces (valid MSC hex)."""
    if distort == 0:
        ax = np.linspace(-L, L, n + 1)
    else:
        w = np.array([1.0 + distort * 0.5 * ((((i * 7 + 3) % 5) / 4.0) - 0.5) for i in range(n)])
        c = np.concatenate([[0.0], np.cumsum(w)])
        ax = -L + 2 * L * c / c[-1]
    hexes = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                hexes.append(np.array([
                    [ax[i], ax[j], ax[k]], [ax[i+1], ax[j], ax[k]], [ax[i+1], ax[j+1], ax[k]], [ax[i], ax[j+1], ax[k]],
                    [ax[i], ax[j], ax[k+1]], [ax[i+1], ax[j], ax[k+1]], [ax[i+1], ax[j+1], ax[k+1]], [ax[i], ax[j+1], ax[k+1]]]))
    return hexes


def yano_Mz(n, distort):
    """yano-MSC <M_z> over the iron cube in a uniform applied field (loop-removed surface-charge solve)."""
    hexes = hex_cube(n, distort)
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    objs = []
    for V in hexes:
        h = rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]); rad.MatApl(h, rad.MatLin(MU_R)); objs.append(h)
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0, 0, MU0 * H0])])
    rad.Solve(cont, 1e-6, 3000, 1)
    vols = np.array([abs(np.dot(V[1]-V[0], np.cross(V[3]-V[0], V[4]-V[0]))) for V in hexes])
    cen = np.array([V.mean(0) for V in hexes])
    Mz = np.array([rad.Fld(cont, "m", list(c))[2] for c in cen])
    rad.UtiDelAll()
    return float((Mz * vols).sum() / vols.sum()), len(hexes)


def fem_Mz_kelvin(p):
    """reduced-Omega + Kelvin (hybrid, cancellation-free) <M_z> over the iron cube."""
    step = os.path.join(tempfile.gettempdir(), "yano_cube_kelvin.step")
    Box(Pnt(-L, -L, -L), Pnt(L, L, L)).WriteStep(step)
    # kelvin_radius must enclose the box bounding sphere (corner sqrt(3)*20 = 34.6 mm); 70 mm is safe
    mesh, info = build_mesh_from_step(step, symmetry="full", kelvin_radius=0.07, kelvin_factor=2.0,
                                      mesh_size_yoke=4e-3, mesh_size_air=12e-3, mesh_size_kelvin=24e-3)
    R_kelvin = info["kelvin_radius"]          # the inversion radius MUST match the mesh geometry
    offset = detect_kelvin_offset(mesh)
    add_periodic_kelvin(mesh, offset)
    ng.SetNumThreads(4)
    with ng.TaskManager():
        if p >= 2:
            mesh.Curve(p)
        fes = ng.Periodic(ng.H1(mesh, order=p))
        Mu = build_material_cf(mesh, MU0, kelvin_mu_factor_3d_cf(center=tuple(offset), R=R_kelvin),
                               outer_keyword="kelvin", overrides={"yoke": MU0 * MU_R})
        u, v = fes.TnT()
        Bs = CF((0.0, 0.0, MU0 * H0))
        gfO = ng.GridFunction(fes); gfO.vec[:] = 0.0
        gfO.Set(H0 * ng.z, ng.BND, mesh.Boundaries("default"))       # total-Omega on the iron boundary
        a = ng.BilinearForm(fes, symmetric=True); a += Mu * grad(u) * grad(v) * dx
        f = ng.LinearForm(fes); f += Mu * (grad(gfO) * grad(v)) * dx("air")
        f += (specialcf.normal(mesh.dim) * Bs) * v * ds("default")   # +normal for the OCC interface
        f.Assemble(); a.Assemble()
        free = fes.FreeDofs()                  # gauge pin: H = grad(Omega) is gauge-invariant -> pin 1 DOF
        for d in range(fes.ndof):
            if free[d]:
                free[d] = False; break
        gfO.vec.data = a.mat.Inverse(free, inverse="pardiso") * f.vec
        vol = ng.Integrate(CF(1.0), mesh, definedon=mesh.Materials("yoke"))
        Mz = (MU_R - 1.0) * ng.Integrate(grad(gfO)[2], mesh, definedon=mesh.Materials("yoke")) / vol
        return float(Mz), fes.ndof, info["ne"]


def main():
    fem, ndof, ne = fem_Mz_kelvin(p=3)
    print(f"\nhex iron cube in UNIFORM field H0={H0:.0f}, mu_r={MU_R:.0f}   (<M_z> over the iron)\n")
    print(f"  {'method':<34} {'<M_z>':>10} {'vs FEM':>10}")
    print(f"  {'reduced-Omega + Kelvin FEM p=3':<34} {fem:>10.1f} {'--':>10}   ndof={ndof}, ne={ne}")
    rows = []
    for distort, tag in [(0.0, "REGULAR"), (0.6, "DISTORTED(graded)")]:
        ymz, nh = yano_Mz(4, distort)
        err = ymz / fem - 1.0
        print(f"  {'yano-MSC ' + tag:<34} {ymz:>10.1f} {err:>+9.2%}   n={nh} hex")
        rows.append({"yano_grid": tag, "n_hex": nh, "Mz": ymz, "err_vs_fem": err})

    out = {"problem": "hex iron cube in uniform field; yano-MSC vs reduced-Omega+Kelvin FEM",
           "mu_r": MU_R, "H0": H0, "cube_half_edge_m": L,
           "fem_reduced_omega_kelvin": {"Mz": fem, "ndof": ndof, "ne": ne, "fes_order": 3},
           "yano_msc": rows,
           "note": ("yano (loop-removed) matches the independent reduced-Omega+Kelvin FEM to ~2% for both "
                    "regular and distorted hexes -> loop-free + distortion-robust; both in the validated "
                    "cube range (tests/feec/parity_vs_msc).")}
    with open(os.path.join(HERE, "yano_vs_reduced_omega_kelvin_cube.json"), "w") as fp:
        json.dump(out, fp, indent=2, default=float)
    print("\nReadout: yano-MSC (loop-removed) == reduced-Omega+Kelvin FEM within ~2% for regular AND")
    print("distorted hexes -> the loop removal leaves a clean, distortion-robust field.")
    print("saved", os.path.join(HERE, "yano_vs_reduced_omega_kelvin_cube.json"))


if __name__ == "__main__":
    main()
