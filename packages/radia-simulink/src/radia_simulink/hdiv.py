"""Thin MATLAB publisher of the Radia-owned HDiv device setup.

Staging loads the CUDA-linked native module, so MATLAB, which holds its own CUDA
runtime, runs this module as a separate process and reads the setup back from HDF5.
"""

_CHARGE_FIELDS = ("cell_vertices", "face_vertices", "charge_host", "charge_kind", "charge_exponents", "tet_rule",
                  "tri_rule", "tet_low_rule", "tri_low_rule", "cell_class", "cell_dofs", "cell_tables",
                  "cell_exponents", "face_class", "face_dofs", "face_tables", "face_exponents")
_EVALUATION_FIELDS = ("rule", "reference_shapes", "vertices", "element_class", "element_dofs", "element_points")
# Present only for the curved P2 family (``curve_order=2``) or the HEX family; the MEX rejects them otherwise.
_CURVED_CHARGE_FIELDS = ("cell_numbers", "face_numbers", "curve_gl", "curve_gw", "cell_mass_rule", "face_mass_rule")
_HEX_CHARGE_FIELDS = ("hex_rules", "hex_bands")
_FAMILY_FIELDS = {"affine-tet": (), "curved-tet": _CURVED_CHARGE_FIELDS, "hex": _HEX_CHARGE_FIELDS}


def prepare_hdiv_device_setup(mesh_path, *, backend, label_contract, mesh_check_receipt,
                              mesh_check_receipt_sha256, mesh_check_provenance, order=1, gram_aca=None,
                              gram_compression="aca", curve_order=1):
    """Publish ``radia.hdiv-device-setup.v1`` for the plain MEX ``radia_hdiv_device``.

    Every table crosses as an owned flat C-order column (float64 or int32), the
    layout the MEX reads; dimensions are implied by the counts and checked
    natively. Only primitive tables and the NG-assembled source/mean vectors
    cross; the device builds every operator. Caller-owned NGSolve TaskManager.
    """
    import numpy as np
    from radia.hdiv_gpu import prepare_hdiv_device_setup as prepare
    data = prepare(mesh_path, backend=backend, label_contract=label_contract, mesh_check_receipt=mesh_check_receipt,
                   mesh_check_receipt_sha256=mesh_check_receipt_sha256, mesh_check_provenance=mesh_check_provenance,
                   order=order, gram_aca=gram_aca, gram_compression=gram_compression, curve_order=curve_order)

    def flat(value, dtype):
        return np.array(value, dtype=dtype, order="C").reshape(-1, 1)

    def ints(name, source):
        value = np.asarray(source[name])
        if value.dtype != np.int32:
            raise ValueError("Device setup integer table %s must be int32" % name)
        return flat(value, np.int32)

    gram = data["charge_gram"]
    fields = _CHARGE_FIELDS + _FAMILY_FIELDS[data["family"]]
    charges = {name: ints(name, gram) if np.asarray(gram[name]).dtype.kind in "iu" else flat(gram[name], np.float64)
               for name in fields}
    aca, compression = gram["aca"], str(gram["compression"])
    charges.update(far_factor=float(gram["far_factor"]), compression=compression, dofs=float(data["dofs"]),
                   aca=np.zeros((0, 0)) if compression == "dense" else dict(eps=float(aca["eps"]), eta=float(aca["eta"]),
                       leaf_size=float(aca["leaf_size"]), max_rank=float(aca["max_rank"]), rank_cap=str(aca["rank_cap"])))
    tables = data["evaluation_tables"]
    evaluation = {name: ints(name, tables) if np.asarray(tables[name]).dtype.kind in "iu" else flat(tables[name], np.float64)
                  for name in _EVALUATION_FIELDS}
    return dict(schema=data["schema"], backend="gpu", family=data["family"], order=float(data["order"]),
                dofs=float(data["dofs"]), mesh_sha256=data["mesh_sha256"],
                label_contract_sha256=data["label_contract_sha256"], charges=charges, evaluation=evaluation,
                rhs_operator=np.array(data["rhs_operator"], dtype=np.float64, order="C"),
                mean_operator=np.array(data["mean_operator"], dtype=np.float64, order="C"))


def write_hdiv_device_setup(request, output):
    """Write the published setup to a new HDF5 file: SI units, 0-based indices, row-major tables."""
    import h5py
    expected = {"mesh_path", "backend", "label_contract", "mesh_check_receipt", "mesh_check_receipt_sha256",
                "mesh_check_provenance", "order", "curve_order", "gram_compression"}
    if not expected <= set(request) <= expected | {"gram_aca"}:
        raise ValueError("Invalid HDiv device setup request fields: %s" % sorted(request))
    keywords = dict(request)
    mesh_path = keywords.pop("mesh_path")
    if "gram_aca" in keywords:
        kinds = dict(eps=float, eta=float, leaf_size=int, max_rank=int, rank_cap=str)
        keywords["gram_aca"] = {key: kinds[key](value) for key, value in keywords["gram_aca"].items()}
    setup = prepare_hdiv_device_setup(mesh_path, **keywords)
    with h5py.File(output, "x") as f:
        f.attrs.update(units="SI", index_base=0, schema=setup["schema"], backend=setup["backend"], family=setup["family"],
                       order=setup["order"], dofs=setup["dofs"], mesh_sha256=setup["mesh_sha256"],
                       label_contract_sha256=setup["label_contract_sha256"])
        f["rhs_operator"], f["mean_operator"] = setup["rhs_operator"], setup["mean_operator"]
        for name in ("charges", "evaluation"):
            group = f.create_group(name)
            for key, value in setup[name].items():
                if isinstance(value, dict):
                    group.create_group(key).attrs.update(value)
                elif key == "aca":
                    continue  # the dense Gram carries no ACA parameters
                elif hasattr(value, "ndim"):
                    group[key] = value.ravel()
                else:
                    group.attrs[key] = value


def main(argv=None):
    """``python -m radia_simulink.hdiv <request.json> <setup.h5>``; this process owns its TaskManager."""
    import argparse
    import json
    import ngsolve
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("request")
    parser.add_argument("output")
    args = parser.parse_args(argv)
    with open(args.request, encoding="utf-8") as stream:
        request = json.load(stream)
    with ngsolve.TaskManager():
        write_hdiv_device_setup(request, args.output)


if __name__ == "__main__":
    main()
