"""Independent original-operator CPU gate for a retained-geometry material sweep.

MAT exchange is validation-only. GPU production needs no Python adapter.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time
import re
import ntpath


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(args):
    result = dict(schema="radia.hcurl-geometry-independent-cpu.v1", accepted=False,
                  timing_accepted=False, accepted_for_adoption=False, cases=[],
                  scope="Independent original-operator validation, not GPU timing")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.exists():
        raise FileExistsError("Never overwrite previous validation evidence")
    started = time.perf_counter()
    try:
        report = json.loads(args.mat.with_suffix(args.mat.suffix + ".json").read_text())
        module = json.loads(args.module.read_text(encoding="utf-8-sig"))
        receipt = json.loads(args.receipt.read_text())
        binaries = [value for key, value in receipt["sha256"].items()
                    if key.endswith("/radia_cuda_ams.mexw64")]
        if len(binaries) != 1 or module["sha256"] != binaries[0] or receipt["source_dirty"]:
            raise RuntimeError("Actual loaded MEX/build identity mismatch")
        if not re.fullmatch('[0-9a-f]{40}', receipt['source_commit']):
            raise RuntimeError('Full build source SHA required')
        if ntpath.normcase(ntpath.normpath(module['path'])) != ntpath.normcase(ntpath.normpath(report['loaded_binary'])):
            raise RuntimeError('Reported MEX path differs from actually loaded module')
        if report["geometry"]["mesh_sha256"] != sha(args.vol):
            raise RuntimeError("Independent CPU mesh differs from geometry snapshot")
        if not report["internal_parity_accepted"] or not report["registry_restored"]:
            raise RuntimeError("Internal geometry prerequisites failed")
        lifecycle = getattr(args, 'lifecycle_evidence', False)
        eligible = (receipt.get('distribution_eligible') is True
                    and receipt.get('geometry_lifecycle_test_build') is False)
        if not eligible and not (lifecycle and receipt.get('geometry_lifecycle_test_build') is True
                                 and receipt.get('distribution_eligible') is False):
            raise RuntimeError('Production distribution eligibility required; test receipt needs --lifecycle-evidence')
        dirichlet = report.get('dirichlet')
        labels = report.get('boundary_labels')
        if not isinstance(dirichlet, str) or not dirichlet or not isinstance(labels, list) or not labels or any(not isinstance(x, str) or not x for x in labels) or len(set(labels)) != len(labels):
            raise RuntimeError('Invalid declared Dirichlet/boundary labels contract')
        import numpy as np
        from scipy.io import loadmat
        result.update(mex_sha256=binaries[0], source_commit=receipt["source_commit"],
                      module_receipt_sha256=sha(args.module), build_receipt_sha256=sha(args.receipt),
                      input_mat_sha256=sha(args.mat), mesh_sha256=sha(args.vol),
                      build_classification=receipt.get('build_classification','Unclassified historical receipt'),
                      numpy=np.__version__)
        payload = loadmat(args.mat, simplify_cells=True)["solutions"]
        cases = [payload] if isinstance(payload, dict) else list(np.asarray(payload, dtype=object).ravel())
        if len(cases) != 3:
            raise RuntimeError("Expected three declared material/source cases")
        declared = np.asarray(report['parameters'], dtype=float)
        sources = np.asarray(report['current_sources'], dtype=float)
        if declared.shape != (3,2) or sources.shape != (3,3) or not np.isfinite(declared).all() or not np.isfinite(sources).all():
            raise RuntimeError('Invalid declared material/source contract')
        for i, case in enumerate(cases):
            if not np.array_equal([case['reluctivity'],case['mass']],declared[i]) or not np.array_equal(np.asarray(case['current']).ravel(),sources[i]):
                raise RuntimeError('MAT solutions differ from declared material/RHS identity')
        import ngsolve as ng
        result['ngsolve']=ng.__version__
        ng.SetNumThreads(1)
        # The validation driver owns the CPU region; no GPU workflow calls it.
        with ng.TaskManager():
            mesh = ng.Mesh(str(args.vol))
            if set(mesh.GetBoundaries()) != set(labels):
                raise RuntimeError('Independent mesh boundary labels differ from geometry report')
            result.update(dirichlet=dirichlet, boundary_labels=labels,
                          lifecycle_evidence=lifecycle)
            fes = ng.HCurl(mesh, order=1, nograds=True, dirichlet=dirichlet)
            free = np.asarray(list(fes.FreeDofs()), dtype=bool)
            for case in cases:
                u, v = fes.TnT()
                form = ng.BilinearForm(fes)
                form += (float(case["reluctivity"]) * ng.curl(u) * ng.curl(v)
                         + float(case["mass"]) * u * v) * ng.dx
                rhs = ng.LinearForm(fes)
                rhs += ng.CoefficientFunction(tuple(np.asarray(case["current"]).ravel())) * v * ng.dx
                form.Assemble(); rhs.Assemble()
                rhs.vec.FV().NumPy()[~free] = 0
                inverse = form.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
                reference = rhs.vec.CreateVector(); reference.data = inverse * rhs.vec
                truth = reference.FV().NumPy().copy()
                checks = []
                for key, scale in (("value", 1.0), ("scaled_value", 2.0)):
                    value = np.asarray(case[key], dtype=float).ravel()
                    if value.size != fes.ndof or not np.isfinite(value).all():
                        raise RuntimeError("Invalid solution shape/values")
                    trial = rhs.vec.CreateVector(); trial.FV().NumPy()[:] = value
                    applied = rhs.vec.CreateVector(); applied.data = form.mat * trial
                    error = float(np.linalg.norm(value - scale * truth) / np.linalg.norm(scale * truth))
                    residual = float(np.linalg.norm((applied.FV().NumPy() - scale * rhs.vec.FV().NumPy())[free])
                                     / np.linalg.norm(scale * rhs.vec.FV().NumPy()[free]))
                    checks.append(dict(scale=scale, solution_relative=error, true_relative_residual=residual))
                    if not np.isfinite([error, residual]).all() or error > 1e-8 or residual > 1e-6:
                        raise RuntimeError(f"Independent CPU/true residual failed: {checks[-1]}")
                result["cases"].append(dict(reluctivity=float(case["reluctivity"]),
                                            mass=float(case["mass"]), checks=checks))
                del inverse
        result["accepted"] = True
        result['accepted_for_adoption'] = eligible and not lifecycle
    except BaseException as error:
        result["error"] = repr(error)
        raise
    finally:
        result["validation_seconds"] = time.perf_counter() - started
        args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lifecycle-evidence', action='store_true',
                        help='Allow test-only receipt; never qualifies adoption')
    for name in ("vol", "mat", "module", "receipt", "out"):
        parser.add_argument("--" + name, type=Path, required=True)
    main(parser.parse_args())
