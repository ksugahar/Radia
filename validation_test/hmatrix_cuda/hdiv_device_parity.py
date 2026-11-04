"""Same-kernel pybind11 producer for the device HDiv frontend parity (PYBIND_MEX_BITWISE.md).

``make-input`` writes the canonical input HDF5 of one case: the NGSolve
device-construction tables of a 6:1:1 bar (``device-flat-tet`` a Netgen mesh, BDM2;
``device-mixed`` a connected tet, pyramid and hex grid), the production inverse-BH
payload of a fixed soft-iron table, every input vector and every solver control. ``pybind`` reads only that file, runs ``_hdiv_cuda`` and writes
the pybind11 record. ``produceHDivDeviceMexRecord.m`` runs the plain MEX on the
same file. Both frontends link one nvcc kernel object; the comparator checks it.
This is frontend parity, not CPU/GPU accuracy or performance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

CASES = ("device-flat-tet", "device-mixed")
MIXED_SHAPE = (8, 2, 2)
OUTPUTS = ("apply", "linear_solution", "nonlinear_solution", "nonlinear_residual_history")
CHARGE_FIELDS = ("cell_vertices", "face_vertices", "charge_host", "charge_kind", "charge_exponents", "tet_rule",
                 "tri_rule", "tet_low_rule", "tri_low_rule", "cell_class", "cell_dofs", "cell_tables",
                 "cell_exponents", "face_class", "face_dofs", "face_tables", "face_exponents")
MIXED_GRAM_FIELDS = ("hex_rules", "hex_bands", "host_types")
MIXED_CHARGE_FIELDS = ("host_class", "class_info", "host_dofs", "class_tables")
EVALUATION_FIELDS = ("rule", "reference_shapes", "vertices", "element_class", "element_dofs", "element_points")
# A 6:1:1 bar gives well-separated clusters, so device ACA builds low-rank leaves (the unit cube did not).
ACA = dict(eps=1e-10, eta=2.0, leaf_size=32, max_rank=200, rank_cap="dense")
BAR, MAXH = (6., 1., 1.), .4
# Soft-iron [H A/m, B T] rows; any strictly increasing origin-anchored table serves.
BH_TABLE = ((0., 0.), (100., .5), (300., 1.), (1000., 1.4), (5000., 1.7), (20000., 1.9), (100000., 2.1))
LINEAR = dict(inv_chi=.1, tolerance=1e-11, maximum_iterations=2000, preconditioner="mass-diagonal")
NEWTON = dict(nonlinear_tolerance=1e-10, linear_tolerance=5e-14, maximum_iterations=40,
              maximum_linear_iterations=10000, preconditioner="reference-diagonal", degree=32, lower_fraction=.01)
APPLIED_FIELD = (0., 0., 1e4)
COMMON_RECEIPT = ("source_sha", "kernel_object_sha256", "cuda_compiler", "python_binary_sha256",
                  "mex_binary_sha256")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def device_identity() -> dict:
    row = subprocess.check_output(["nvidia-smi", "--query-gpu=uuid,compute_cap,name", "--format=csv,noheader"],
                                  text=True).strip().splitlines()
    if len(row) != 1:
        raise RuntimeError("Parity requires exactly one visible GPU, got %d" % len(row))
    uuid, cap, name = (part.strip() for part in row[0].split(","))
    return dict(device_uuid=uuid, architecture="sm_" + cap.replace(".", ""), device_name=name)


def make_input(out: Path, order: int, maxh: float = MAXH, leaf_size: int = ACA["leaf_size"], length: float = BAR[0],
               case: str = CASES[0], shape: tuple = MIXED_SHAPE) -> None:
    import h5py
    import ngsolve as ng
    import numpy as np
    from netgen.csg import CSGeometry, OrthoBrick, Pnt
    from radia import hdiv_gpu
    from radia.gpu import _native
    from radia.vim._nonlinear import _bh_inverse_payload, _validate_bh_table
    mesh_path = out.with_suffix(".vol")
    if case == "device-mixed":
        from hdiv_mixed_meshes import mixed_grid
        mixed_grid(shape, lengths=(length, BAR[1], BAR[2])).ngmesh.Save(str(mesh_path))
    else:
        geometry = CSGeometry()
        geometry.Add(OrthoBrick(Pnt(0, 0, 0), Pnt(length, BAR[1], BAR[2])))
        geometry.GenerateMesh(maxh=maxh).Save(str(mesh_path))
    with ng.TaskManager():
        setup = hdiv_gpu._prepare_device_components(mesh_path, order, _native("_hdiv_cuda"), aca=dict(ACA, leaf_size=leaf_size),
                                                    compression="aca")
    table = np.asarray(BH_TABLE)
    _validate_bh_table(table, context="device parity")
    law = _bh_inverse_payload(table[:, 0], table[:, 1])
    nd = setup["dofs"]
    rhs = setup["rhs_operator"] @ np.asarray(APPLIED_FIELD)
    with h5py.File(out, "x") as f:
        f.attrs["schema"] = "radia.hdiv-device-parity-input.v1"
        f.attrs["case"], f.attrs["order"], f.attrs["dofs"] = case, order, nd
        f.attrs["mesh_sha256"] = sha256(mesh_path)
        gram = setup["charge_gram"]
        charges = f.create_group("charges")
        for name in CHARGE_FIELDS:
            charges[name] = np.ascontiguousarray(gram[name])
        if case == "device-mixed":
            for name in MIXED_GRAM_FIELDS:
                charges[name] = np.ascontiguousarray(gram[name])
            for name in MIXED_CHARGE_FIELDS:
                charges[name] = np.ascontiguousarray(gram["mixed_charges"][name])
        charges.attrs["far_factor"], charges.attrs["compression"] = float(gram["far_factor"]), gram["compression"]
        charges.attrs["family"] = setup["family"]
        for key, value in gram["aca"].items():
            charges.attrs["aca_" + key] = value
        for index, tables in enumerate(setup["evaluation_tables"]["tables"], 1):
            evaluation = f.create_group("evaluation/group%d" % index)
            for name in EVALUATION_FIELDS:
                evaluation[name] = np.ascontiguousarray(tables[name])
        group = f.create_group("law")
        group["knots"] = law["knots"]
        group["coefficients"] = np.asarray(law["coefficients"]).reshape(-1, 4)
        group["integral_base"] = law["integral_base"]
        for key in ("cap", "barrier_slope", "last_field", "initial_reluctivity"):
            group.attrs[key] = float(law[key])
        vectors = f.create_group("vectors")
        vectors["apply_input"] = np.sin(np.arange(1, nd + 1, dtype=np.float64))
        vectors["rhs"] = np.ascontiguousarray(rhs)
        vectors["initial"] = np.zeros(nd)
        for name, controls in (("linear", LINEAR), ("newton", NEWTON)):
            group = f.create_group(name)
            for key, value in controls.items():
                group.attrs[key] = value
    print(json.dumps(dict(input=str(out), input_sha256=sha256(out), dofs=nd)))


def text(value) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def run_pybind(input_path: Path, out: Path, receipt: dict) -> None:
    import h5py
    import numpy as np
    from radia.gpu import _native
    cuda = _native("_hdiv_cuda")
    digest = sha256(input_path)
    with h5py.File(input_path, "r") as f:
        case = text(f.attrs["case"])
        if text(f.attrs["schema"]) != "radia.hdiv-device-parity-input.v1" or case not in CASES:
            raise ValueError("Not a device parity input")
        nd = int(f.attrs["dofs"])
        c = f["charges"]
        charges = {name: c[name][...] for name in CHARGE_FIELDS}
        if case == "device-mixed":
            charges.update({name: c[name][...] for name in MIXED_GRAM_FIELDS},
                           mixed_charges={name: c[name][...] for name in MIXED_CHARGE_FIELDS})
        aca = dict(eps=float(c.attrs["aca_eps"]), eta=float(c.attrs["aca_eta"]), leaf_size=int(c.attrs["aca_leaf_size"]),
                   max_rank=int(c.attrs["aca_max_rank"]), rank_cap=text(c.attrs["aca_rank_cap"]))
        evaluation = dict(tables=[{name: f["evaluation"]["group%d" % k][name][...] for name in EVALUATION_FIELDS}
                                  for k in range(1, len(f["evaluation"]) + 1)])
        law = dict(knots=f["law/knots"][...], coefficients=np.ascontiguousarray(f["law/coefficients"][...].ravel()),
                   integral_base=f["law/integral_base"][...], **{k: float(v) for k, v in f["law"].attrs.items()})
        vectors = {name: f["vectors"][name][...] for name in ("apply_input", "rhs", "initial")}
        linear, newton = dict(f["linear"].attrs), dict(f["newton"].attrs)
        compression, far, family = text(c.attrs["compression"]), float(c.attrs["far_factor"]), text(c.attrs["family"])
    operator = cuda._ResidentHDiv.from_flat_charges(**charges, far_factor=far, dofs=nd, compression=compression,
                                                    aca=aca, backend="gpu", family=family)
    mass = cuda._ResidentMass.from_flat_tables(**evaluation, dofs=nd, backend="gpu")
    material = cuda._ResidentBHMaterial.from_flat_tables(**evaluation, dofs=nd, law=law, backend="gpu")
    try:
        counters = operator.assembly()["hmatrix"]
        applied = operator.apply(vectors["apply_input"])
        lin = operator.solve_linear(mass, vectors["rhs"], float(linear["inv_chi"]), float(linear["tolerance"]),
                                    int(linear["maximum_iterations"]), text(linear["preconditioner"]))
        new = operator.solve_nonlinear_material(material, vectors["rhs"], vectors["initial"],
                                                float(newton["nonlinear_tolerance"]), float(newton["linear_tolerance"]),
                                                int(newton["maximum_iterations"]), int(newton["maximum_linear_iterations"]),
                                                text(newton["preconditioner"]), int(newton["degree"]),
                                                float(newton["lower_fraction"]))
    finally:
        material.close()
        mass.close()
        operator.close()
    if int(counters.get("lowrank_leaves", 0)) < 1:
        raise RuntimeError("Parity input did not exercise device ACA low-rank leaves: %r" % (counters,))
    outputs = dict(apply=applied, linear_solution=lin["x"], nonlinear_solution=new["x"],
                   nonlinear_residual_history=np.array([row["relative_residual"] for row in new["history"]]))
    write_record(out, dict(receipt, frontend="pybind11", case=case, input_sha256=digest, output_names=list(OUTPUTS),
                           artifact_sha256=sha256(Path(cuda.__file__)), hmatrix=dict(counters),
                           linear_relative_residual=lin["relative_residual"],
                           nonlinear_relative_residual=new["relative_residual"], **device_identity()), outputs)


def write_record(out: Path, receipt: dict, outputs: dict) -> None:
    import h5py
    import numpy as np
    with h5py.File(out, "x") as f:
        f.attrs["receipt_json"] = json.dumps(receipt, sort_keys=True)
        group = f.create_group("outputs")
        for name in OUTPUTS:
            group[name] = np.ascontiguousarray(outputs[name], dtype=np.float64)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    make = sub.add_parser("make-input")
    make.add_argument("--out", type=Path, required=True)
    make.add_argument("--order", type=int, choices=(1, 2), default=2)
    make.add_argument("--maxh", type=float, default=MAXH)
    make.add_argument("--leaf-size", type=int, default=ACA["leaf_size"])
    make.add_argument("--length", type=float, default=BAR[0])
    make.add_argument("--case", choices=CASES, default=CASES[0])
    make.add_argument("--shape", type=int, nargs=3, default=MIXED_SHAPE)
    run = sub.add_parser("pybind")
    run.add_argument("--input", type=Path, required=True)
    run.add_argument("--out", type=Path, required=True)
    run.add_argument("--common-receipt", type=Path, required=True,
                     help="JSON shared by both producers: " + ", ".join(COMMON_RECEIPT))
    args = parser.parse_args()
    if args.command == "make-input":
        make_input(args.out, args.order, args.maxh, args.leaf_size, args.length, args.case, tuple(args.shape))
    else:
        run_pybind(args.input, args.out, common_receipt(args.common_receipt))


def common_receipt(path: Path) -> dict:
    """Build identity written once from the pyd and MEX build receipts of one source SHA."""
    receipt = json.loads(path.read_text(encoding="utf-8"))
    if set(receipt) != set(COMMON_RECEIPT) or not all(receipt.values()):
        raise ValueError("Common receipt requires exactly %s" % (COMMON_RECEIPT,))
    return receipt


if __name__ == "__main__":
    main()
