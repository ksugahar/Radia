"""Independent HCurl scattered-A reference for a prescribed four-patch Zs.

Default invocation writes a pending comparison plan only. --solve is explicit;
it does not authorize compute admission or certify nonlinear ESIM physics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
import sys
import platform
import importlib.metadata

import numpy as np

HERE = Path(__file__).resolve().parent
MU0 = 4e-7 * np.pi


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_case(path):
    case = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if case.get("schema") != "radia.ih-element-impedance-case.v1":
        raise ValueError("Unsupported element-impedance case")
    vertices = np.asarray(case["vertices_m"], dtype=float)
    raw_triangles = np.asarray(case["triangles_one_based"], dtype=float)
    if not np.isfinite(raw_triangles).all() or np.any(raw_triangles != raw_triangles.astype(int)):
        raise ValueError("Integer one-based triangle indices required")
    triangles = raw_triangles.astype(int) - 1
    zs = np.asarray(case["patch_zs_ohm_real"]) + 1j * np.asarray(case["patch_zs_ohm_imag"])
    if vertices.shape != (4, 3) or triangles.shape != (4, 3) or zs.shape != (4,):
        raise ValueError("This independently authored case requires four tetrahedron patches")
    if not np.isfinite(vertices).all() or not np.isfinite(zs).all() or np.any(zs.real <= 0):
        raise ValueError("Nonfinite geometry or nonpassive impedance")
    if np.any(triangles < 0) or np.any(triangles >= 4):
        raise ValueError("Invalid one-based triangle indices")
    if len({tuple(sorted(t)) for t in triangles}) != 4 or any(len(set(t)) != 3 for t in triangles):
        raise ValueError("Four distinct tetrahedron faces required")
    for tri in triangles:
        a, b, c = vertices[tri]
        if np.dot(np.cross(b-a, c-a), (a+b+c)/3-vertices.mean(axis=0)) <= 0:
            raise ValueError("Outward oriented nondegenerate patches required")
    if not np.isfinite([case["frequency_Hz"], case["current_peak_A"]]).all() or case["frequency_Hz"] <= 0 or case["current_peak_A"] <= 0:
        raise ValueError("Positive frequency and peak current required")
    if case["phasor"] != "exp(+i omega t), peak":
        raise ValueError("Unexpected phasor convention")
    return case, vertices, triangles, zs


def plan(case_path, run_index, order, offset):
    record = dict(schema="radia.ih-element-impedance-comparison-plan.v1",
                status="not-performed", case_sha256=digest(case_path),
                source_sha256=digest(__file__), run_index=run_index,
                filament_gauss_order=order, surface_probe_offset_m=offset,
                reference="independently meshed HCurl scattered-A, sparsecholesky",
                boundary_law="four prescribed complex patches, never averaged",
                gates=[
                    dict(quantity=q, norm=n, region=r, normalization=normalization,
                         absolute_floor=floor, tolerance=tol, rationale=reason)
                    for q, n, r, normalization, floor, tol, reason in (
                        ("power", "absolute scalar difference", "body", "FEM body power", 1e-18, .02, "0.5% FEM mesh + 0.5% exterior + 1% BEM budget"),
                        ("reaction_impedance", "complex magnitude difference", "filament port", "FEM body reaction magnitude", 1e-12, .03, "body contribution, not supplied coil impedance"),
                        ("q_surf", "area-weighted L2 difference", "all patches", "FEM area-weighted L2", 1e-18, .10, "separate local derivative/loss error budget"),
                        ("tangential_H", "complex vector area-weighted L2 difference", "all patches", "FEM area-weighted vector L2", 1e-12, .05, "phase-preserving local derivative budget"),
                        ("linear_residual", "Euclidean", "free HCurl DOFs", "RHS norm", 1e-30, 1e-6, "mandatory linear true-residual contract"),
                        ("power_identity", "absolute scalar difference", "body and filament port", "FEM body power", 1e-18, .01, "independent work versus surface dissipation"))],
                additional_controls="PLAN.md: mesh, exterior, BEM refinement, quadrature and probe-offset convergence are also required",
                nonlinear_volumetric_validation="not-performed",
                experimental_validation="not-performed")
    for label, tolerance in (("fem_mesh", .005), ("fem_exterior", .005),
                             ("bem_refinement", .01), ("filament_quadrature", .001)):
        for quantity, region, floor in (("power", "body", 1e-18), ("reaction_impedance", "filament port", 1e-12)):
            record["gates"].append(dict(quantity=label+"_"+quantity, norm="complex/scalar magnitude difference",
                region=region, normalization="declared finer reference magnitude", absolute_floor=floor,
                tolerance=tolerance, rationale="Independent numerical error-budget component fixed in PLAN.md"))
    for quantity, floor in (("tangential_H", 1e-12), ("q_surf", 1e-18)):
        record["gates"].append(dict(quantity="probe_offset_"+quantity, norm="area-weighted vector/scalar L2 difference",
            region="all patches", normalization="half-offset reference area-weighted L2", absolute_floor=floor,
            tolerance=.001, rationale="Air-side trace evaluation must not dominate spatial discretization error"))
    record.update(case=json.loads(Path(case_path).read_text()), command=[sys.executable, *sys.argv],
                  runtime=dict(host=platform.node(), python=sys.executable, python_version=platform.python_version()))
    return record


def read_samples(path, case_path, vertices, triangles):
    samples = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if samples["case_file_sha256"] != digest(case_path):
        raise ValueError("BEM samples and FEM must use the identical case file")
    points = np.asarray(samples["points_m"], dtype=float)
    patches = np.asarray(samples["patch_one_based"], dtype=float).reshape(-1)
    weights = np.asarray(samples["area_weights_m2"], dtype=float).reshape(-1)
    if points.shape != (len(patches), 3) or weights.shape != patches.shape:
        raise ValueError("Invalid common quadrature array shapes")
    if not np.isfinite(points).all() or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("Finite points and positive area weights required")
    if not np.isfinite(patches).all() or np.any(patches != patches.astype(int)) or np.any((patches < 1) | (patches > 4)):
        raise ValueError("Explicit integer one-based patch IDs required")
    for index, tri in enumerate(triangles):
        a, b, c = vertices[tri]
        selected = points[patches == index+1]
        normal = np.cross(b-a, c-a)
        area = np.linalg.norm(normal)/2
        if abs(weights[patches == index+1].sum()/area-1) > 1e-10:
            raise ValueError("Common quadrature must cover each original patch area")
        for point in selected:
            if abs(np.dot(point-a, normal)/np.linalg.norm(normal)) > 1e-10:
                raise ValueError("Comparison point is off its declared patch plane")
            bary = np.linalg.lstsq(np.column_stack((b-a, c-a)), point-a, rcond=None)[0]
            if min(bary) < -1e-10 or bary.sum() > 1+1e-10:
                raise ValueError("Comparison point is outside its declared patch")
    samples["patch_one_based"] = patches.astype(int).tolist()
    samples["area_weights_m2"] = weights.tolist()
    return samples


def make_mesh(case, vertices, triangles, run, work, reuse_result, case_path):
    import ngsolve as ng
    from netgen.occ import Face, Wire, Segment, Pnt, Sew, Solid, Sphere, OCCGeometry
    vol = work / "air.vol"
    if reuse_result:
        prior = json.loads(reuse_result.read_text())
        original = reuse_result.parent / "air.vol"
        if (prior["run"] != run or prior["plan"]["case_sha256"] != digest(case_path)
                or prior["plan"]["source_sha256"] != digest(__file__)
                or prior["ngsolve"] != ng.__version__ or prior["mesh_sha256"] != digest(original)):
            raise ValueError("Reused mesh must match the case, run, source, runtime and original hash")
        shutil.copyfile(original, vol)
    else:
        faces = []
        for index, tri in enumerate(triangles):
            p = [Pnt(*vertices[j]) for j in tri]
            face = Face(Wire([Segment(p[j], p[(j+1) % 3]) for j in range(3)]))
            face.name = f"wp_{index+1}"
            face.maxh = run["surface_maxh_m"]
            faces.append(face)
        body = Solid(Sew(faces, tolerance=1e-10))
        outer = Sphere(Pnt(0, 0, 0), run["outer_radius_m"])
        outer.faces.name = "outer"
        air = outer-body
        air.mat("air")
        mesh = ng.Mesh(OCCGeometry(air).GenerateMesh(maxh=run["air_maxh_m"]))
        mesh.Curve(case["fem_geometry_order"])
        mesh.ngmesh.Save(str(vol))
    from cubit_mesh_export.check import check_consistency
    contract = dict(schema="radia.vol-label-contract.v1", application="element_impedance_crosscheck",
                    required=dict(materials=["air"], boundaries=["outer", "wp_1", "wp_2", "wp_3", "wp_4"]))
    checked = check_consistency(str(vol), contract=contract)
    if not checked["passed"]:
        raise RuntimeError("Independent FEM mesh failed check-vol")
    if reuse_result:
        mesh = ng.Mesh(str(vol))  # Preserve saved curvature; do not re-Curve a checked input.
    for index, tri in enumerate(triangles):
        a, b, c = vertices[tri]
        expected = np.linalg.norm(np.cross(b-a, c-a)) / 2
        actual = ng.Integrate(1, mesh, ng.BND, definedon=mesh.Boundaries(f"wp_{index+1}"))
        if abs(actual/expected-1) > 1e-8:
            raise RuntimeError("FEM patch geometry/area differs from the prescribed planar case")
    return mesh, dict(mesh_sha256=digest(vol), label_check_passed=True)


def incident(case, order):
    """Independent Gauss filament expressions, not a Radia BEM field call."""
    import ngsolve as ng
    gauss, weights = np.polynomial.legendre.leggauss(order)
    coordinates = (ng.x, ng.y, ng.z)
    potential = [0, 0, 0]
    magnetic = [0, 0, 0]
    coil = np.asarray(case["coil_vertices_m"], dtype=float)
    for start, end in zip(coil, np.roll(coil, -1, axis=0)):
        dl = end-start
        for x, w in zip(gauss, weights):
            point = start + (x+1)*dl/2
            r = [coordinates[k]-float(point[k]) for k in range(3)]
            distance = ng.sqrt(sum(x*x for x in r))
            for k in range(3):
                potential[k] += MU0/(4*np.pi)*case["current_peak_A"]*w/2*float(dl[k])/distance
                j, l = (k+1) % 3, (k+2) % 3
                magnetic[k] += case["current_peak_A"]/(4*np.pi)*w/2*(float(dl[j])*r[l]-float(dl[l])*r[j])/distance**3
    return ng.CF(tuple(potential)), ng.CF(tuple(magnetic))


def solve(case, vertices, triangles, zs, args, samples):
    """Helpers follow the caller-owned TaskManager policy."""
    import ngsolve as ng
    run = case["fem_runs"][args.run_index]
    mesh, mesh_record = make_mesh(case, vertices, triangles, run, args.work, args.reuse_result, args.case)
    ai, hi = incident(case, args.quadrature_order)
    omega = 2*np.pi*case["frequency_Hz"]
    space = ng.HCurl(mesh, order=case["fem_hcurl_order"], complex=True, dirichlet="outer")
    u, v = space.TnT()
    normal = ng.specialcf.normal(3)
    tangential = lambda x: x-normal*(x*normal)
    matrix = ng.BilinearForm(space, symmetric=True)
    matrix += (ng.curl(u)*ng.curl(v) + case["fem_gauge_per_m2"]*u*v)/MU0*ng.dx
    rhs = ng.LinearForm(space)
    for index, impedance in enumerate(zs):
        boundary = ng.ds(f"wp_{index+1}")
        robin = 1j*omega/complex(impedance)
        matrix += robin*u.Trace()*v.Trace()*boundary
        rhs += (ng.Cross(normal, hi)-robin*tangential(ai))*v.Trace()*boundary
    matrix.Assemble()
    rhs.Assemble()
    inverse = matrix.mat.Inverse(space.FreeDofs(), inverse="sparsecholesky")
    solution = ng.GridFunction(space)
    solution.vec.data = inverse*rhs.vec
    residual = rhs.vec.CreateVector()
    free = np.asarray(list(space.FreeDofs()), dtype=bool)
    scale = max(np.linalg.norm(rhs.vec.FV().NumPy()[free]), 1e-30)
    for correction in range(4):
        residual.data = rhs.vec-matrix.mat*solution.vec
        relative = float(np.linalg.norm(residual.FV().NumPy()[free])/scale)
        if np.isfinite(relative) and relative <= 1e-6:
            break
        if correction == 3:
            raise RuntimeError(f"FEM true residual rejected: {relative}")
        solution.vec.data += inverse*residual
    total_trace = solution.Trace()+tangential(ai)
    powers = []
    for index, impedance in enumerate(zs):
        ht = -1j*omega/complex(impedance)*ng.Cross(normal, total_trace)
        q = .5*float(impedance.real)*(ht*ng.Conj(ht)).real
        powers.append(float(ng.Integrate(q, mesh, ng.BND, definedon=mesh.Boundaries(f"wp_{index+1}"), order=12)))
    # Composite integration avoids treating a many-element FE trace as one polynomial.
    gauss, weights = np.polynomial.legendre.leggauss(args.quadrature_order)
    coil = np.asarray(case["coil_vertices_m"], dtype=float)
    linkage = 0j
    for start, end in zip(coil, np.roll(coil, -1, axis=0)):
        count = int(np.ceil(np.linalg.norm(end-start)/case["reaction_segment_maxh_m"]))
        dl = (end-start)/count
        for sub in range(count):
            for x, w in zip(gauss, weights):
                point = start+(sub+(x+1)/2)*dl
                linkage += np.dot(np.asarray(solution(mesh(*point))), dl)*w/2
    reaction = 1j*omega/case["current_peak_A"]*linkage
    power = sum(powers)
    out = dict(status="computed-not-yet-crosschecked", **mesh_record, ndof=space.ndof,
               run=run, linear_true_relative_residual=relative, residual_corrections=correction,
               patch_power_W=powers, total_power_W=power,
               reaction_impedance_ohm_real=float(reaction.real), reaction_impedance_ohm_imag=float(reaction.imag),
               power_identity_relative=float(abs(power-.5*reaction.real*case["current_peak_A"]**2)/max(power, 1e-18)))
    if samples is not None:
        fields = []
        for point, patch in zip(samples["points_m"], samples["patch_one_based"]):
            tri = triangles[patch-1]
            a, b, c = vertices[tri]
            n = np.cross(b-a, c-a); n /= np.linalg.norm(n)
            point = np.asarray(point, dtype=float)
            if abs(np.dot(point-a, n)) > 1e-10:
                raise ValueError("A comparison point is not on its declared original patch")
            query = point + args.probe_offset*n
            at = np.asarray(solution(mesh(*query)))+np.asarray(ai(mesh(*query)))
            fields.append(1j*omega/zs[patch-1]*np.cross(n, at))
        fields = np.asarray(fields)
        out.update(H_t_real=fields.real.tolist(), H_t_imag=fields.imag.tolist(),
                   points_m=samples["points_m"], patch_one_based=samples["patch_one_based"],
                   area_weights_m2=samples["area_weights_m2"])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=HERE/"case.json")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--run-index", type=int, choices=range(3), default=0)
    parser.add_argument("--quadrature-order", type=int, choices=(16, 24), default=16)
    parser.add_argument("--probe-offset", type=float, choices=(1e-9, 5e-10), default=1e-9)
    parser.add_argument("--samples", type=Path)
    parser.add_argument("--reuse-result", type=Path, help="Reuse that run's checked air.vol for quadrature/trace controls")
    parser.add_argument("--threads", type=int, choices=range(1, 17), default=2)
    parser.add_argument("--solve", action="store_true")
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=False)
    case, vertices, triangles, zs = read_case(args.case)
    report = plan(args.case, args.run_index, args.quadrature_order, args.probe_offset)
    report["source_commit"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    report["source_dirty"] = bool(subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=all", "--", str(HERE)], cwd=HERE, text=True).strip())
    (args.work/"comparison-plan.json").write_text(json.dumps(report, indent=2)+"\n")
    if not args.solve:
        print("Pending comparison plan written; no mesh or solve performed")
        return
    if report["source_dirty"]:
        raise RuntimeError("Commit and review the exact comparison source before a numerical run")
    samples = read_samples(args.samples, args.case, vertices, triangles) if args.samples else None
    try:
        import ngsolve as ng
        ng.SetNumThreads(args.threads)
        with ng.TaskManager():
            result = solve(case, vertices, triangles, zs, args, samples)
        result.update(plan=report, ngsolve=ng.__version__, threads=args.threads,
                      numpy=np.__version__, netgen=importlib.metadata.version("netgen-mesher"),
                      mesh_checker=importlib.metadata.version("cubit-mesh-export"),
                      nonlinear_volumetric_validation="not-performed", experimental_validation="not-performed")
    except Exception as error:
        result = dict(status="failed", error=str(error), plan=report)
        raise
    finally:
        (args.work/"result.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
