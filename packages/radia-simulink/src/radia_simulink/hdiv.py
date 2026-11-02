"""Thin MATLAB publisher of the Radia-owned HDiv device setup."""

_CHARGE_FIELDS = ("cell_vertices", "face_vertices", "charge_host", "charge_kind", "charge_exponents", "tet_rule",
                  "tri_rule", "tet_low_rule", "tri_low_rule", "cell_class", "cell_dofs", "cell_tables",
                  "cell_exponents", "face_class", "face_dofs", "face_tables", "face_exponents")
_EVALUATION_FIELDS = ("rule", "reference_shapes", "vertices", "element_class", "element_dofs", "element_points")


def prepare_hdiv_device_setup(mesh_path, *, backend, label_contract, mesh_check_receipt,
                              mesh_check_receipt_sha256, mesh_check_provenance, order=1, gram_aca=None):
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
                   order=order, gram_aca=gram_aca)

    def flat(value, dtype):
        return np.array(value, dtype=dtype, order="C").reshape(-1, 1)

    def ints(name, source):
        value = np.asarray(source[name])
        if value.dtype != np.int32:
            raise ValueError("Device setup integer table %s must be int32" % name)
        return flat(value, np.int32)

    gram = data["charge_gram"]
    charges = {name: ints(name, gram) if np.asarray(gram[name]).dtype.kind in "iu" else flat(gram[name], np.float64)
               for name in _CHARGE_FIELDS}
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
