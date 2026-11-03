"""Device nonlinear HDiv route against the CPU route on one checked mesh; JSON record.

Material energy, gradient and tangent action at a state of known saturation,
with central finite differences of the device energy and gradient, a near-root
state, then the whole Newton solve against ``radia.vim.Solve``. The two routes
share the charge-Gram near/far profile and quadrature, so this is operator and
implementation parity, not independent physical accuracy: for that use
``validate_nonlinear_sphere.py``. The limits are those of the retired host-route
checks and are proposed, not adopted.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import numpy as np

LIMITS = dict(material=1e-8, finite_difference=1e-5, coefficients=1e-8, mean=1e-8, nonlinear_residual=1e-10)
H_SAMPLES = np.array([0., 1., 10., 100., 1000., 1e4, 1e5])
MS, CHI = 1e6, 1000.  # Synthetic monotone saturation law; not measured steel.


def relative(a, b):
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b))/max(np.linalg.norm(b), np.finfo(float).tiny))


def bh_table():
    return np.column_stack([H_SAMPLES, 4e-7*np.pi*(H_SAMPLES + MS*CHI*H_SAMPLES/(MS + CHI*H_SAMPLES))])


def cpu_material(fes, table, state, direction):
    import ngsolve as ng
    from radia.vim._nonlinear import _EnergyMaterialQuadrature
    from scipy.sparse import coo_matrix
    material = _EnergyMaterialQuadrature(fes, table)
    energy = material.update(state)
    u, v = fes.TnT()
    load = ng.LinearForm(fes)
    load += material.field*v*material.measure
    load.Assemble()
    form = ng.BilinearForm(fes)
    form += material.bilinear_integrator(ng.InnerProduct(material.tangent*u, v))
    form.Assemble()
    rows, columns, values = form.mat.COO()
    tangent = coo_matrix((np.array(values), (np.array(rows), np.array(columns))), shape=(fes.ndof, fes.ndof)).tocsr()
    return dict(energy=float(energy), gradient=load.vec.FV().NumPy().copy(), tangent_action=tangent @ direction)


def saturating_state(fes, table, seed=20261103):
    """A random state scaled so its largest point magnetization sits on the saturating branch."""
    from radia.vim._nonlinear import _EnergyMaterialQuadrature
    rng = np.random.default_rng(seed)
    material = _EnergyMaterialQuadrature(fes, table)
    state = rng.standard_normal(fes.ndof)
    material.m_field.vec.FV().NumPy()[:] = state
    material.samples.Interpolate(material.m_field)
    points = material.samples.vec.FV().NumPy().reshape(3, material.npoints).T
    direction = rng.standard_normal(fes.ndof)
    return state*.3*MS/np.max(np.linalg.norm(points, axis=1)), direction/np.linalg.norm(direction)


def compare_material(device, cpu):
    return dict(energy=abs(device["energy"] - cpu["energy"])/max(abs(cpu["energy"]), np.finfo(float).tiny),
                gradient=relative(device["gradient"], cpu["gradient"]),
                tangent=relative(device["tangent_action"], cpu["tangent_action"]))


def finite_differences(material, state, direction, at):
    step = 1e-6*max(np.linalg.norm(state), 1.)
    empty = np.zeros(0)
    plus = material.evaluate(np.ascontiguousarray(state + step*direction), empty)
    minus = material.evaluate(np.ascontiguousarray(state - step*direction), empty)
    slope = float(at["gradient"] @ direction)
    return dict(energy_gradient=abs((plus["energy"] - minus["energy"])/(2*step) - slope)/max(abs(slope), np.finfo(float).tiny),
                gradient_tangent=relative((plus["gradient"] - minus["gradient"])/(2*step), at["tangent_action"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--mesh-check-receipt", type=Path, required=True)
    parser.add_argument("--mesh-check-provenance", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--order", type=int, choices=[1, 2], default=2)
    parser.add_argument("--curve-order", type=int, choices=[1, 2], default=1)
    parser.add_argument("--field", type=float, default=1e5)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    import ngsolve as ng
    from radia.hdiv_gpu import PreparedHDivNonlinear
    from radia.vim import Solve
    ng.SetNumThreads(args.threads)
    table = bh_table()
    result = dict(schema="radia.hdiv-device-nonlinear.v1", accepted=False, units="SI",
                  comparison_kind="operator and implementation parity on a matched near/far profile",
                  adoption="proposed", proposed_limits=LIMITS, order=args.order, curve_order=args.curve_order,
                  field_A_per_m=args.field, threads=args.threads, ngsolve_version=ng.__version__,
                  machine=platform.machine(), source_sha=os.environ.get("RADIA_SOURCE_SHA"),
                  mesh_sha256=hashlib.sha256(args.mesh.read_bytes()).hexdigest())
    try:
        receipt = dict(label_contract=str(args.contract), mesh_check_receipt=str(args.mesh_check_receipt),
                       mesh_check_receipt_sha256=hashlib.sha256(args.mesh_check_receipt.read_bytes()).hexdigest(),
                       mesh_check_provenance=json.loads(args.mesh_check_provenance.read_text(encoding="utf-8")))
        with ng.TaskManager():
            mesh = ng.Mesh(str(args.mesh))
            if args.curve_order == 2:
                mesh.Curve(2)
            fes = ng.HDiv(mesh, order=args.order)
            state, direction = saturating_state(fes, table)
            near_root = state*1e-6
            cpu = cpu_material(fes, table, state, direction)
            cpu_near = cpu_material(fes, table, near_root, direction)
            answer = Solve(mesh, bh_table=table, H_ext=ng.CoefficientFunction((0, 0, args.field)), order=args.order,
                           nl_tol=1e-10, nl_maxit=50, newton_inner_tol="fixed", newton_warmstart="linear",
                           preconditioner="mass-riesz", curve_gauss=8)
        with PreparedHDivNonlinear.from_bh_table(args.mesh, table, backend="gpu", order=args.order,
                                                 curve_order=args.curve_order, **receipt) as prepared:
            material = prepared._material
            device = material.evaluate(np.ascontiguousarray(state), np.ascontiguousarray(direction))
            device_near = material.evaluate(np.ascontiguousarray(near_root), np.ascontiguousarray(direction))
            result["material"] = compare_material(device, cpu)
            result["near_root_diagnostic"] = compare_material(device_near, cpu_near)
            result["finite_difference"] = finite_differences(material, state, direction, device)
            solution = prepared.solve((0., 0., args.field))
        coefficients = np.asarray(answer["_m_coefficients"])
        result["solve"] = dict(
            coefficients=relative(solution["x"], coefficients),
            mean=relative(solution["mean_magnetization"], np.asarray(answer["M_avg"])),
            nonlinear_residual=float(solution["relative_residual"]),
            newton_iterations=int(solution["newton_iterations"]), dofs=fes.ndof)
        result["accepted"] = bool(
            max(result["material"].values()) <= LIMITS["material"]
            and max(result["finite_difference"].values()) <= LIMITS["finite_difference"]
            and result["solve"]["coefficients"] <= LIMITS["coefficients"]
            and result["solve"]["mean"] <= LIMITS["mean"]
            and result["solve"]["nonlinear_residual"] <= LIMITS["nonlinear_residual"])
    except BaseException as exc:
        result["error"] = repr(exc)
        raise
    finally:
        (args.out / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main()
