"""Thin MATLAB publisher of the Radia-owned HDiv host setup."""

def prepare_hdiv_host_setup(mesh_path, *, backend, label_contract, mesh_check_receipt,
                            mesh_check_receipt_sha256, mesh_check_provenance, order=1):
    """Prepare owned numeric FE data for explicit MATLAB shared-GPU initialization.

    Caller-owned NGSolve TaskManager; no CPU solve, Engine or resident GPU
    allocation. Returns arrays/sparse matrices in memory, not a required file
    exchange. MATLAB selects Assembly='ngsolve-host-normalized' explicitly.
    Unsupported charge geometry still raises in the existing Radia builder.
    """
    from radia.hdiv_gpu import prepare_hdiv_host_setup as prepare
    data = prepare(mesh_path, backend=backend, label_contract=label_contract, mesh_check_receipt=mesh_check_receipt,
                   mesh_check_receipt_sha256=mesh_check_receipt_sha256,
                   mesh_check_provenance=mesh_check_provenance, order=order, maximum_dofs=4096)
    mesh_hash = data["mesh_sha256"]
    contract_hash = data["label_contract_sha256"]
    import numpy as np
    snapshot = data["snapshot"]
    if data["dofs"] > 4096 or snapshot["dimension"] > 4608 or len(data["weights"]) > 65536:
        raise ValueError("MATLAB shared host setup exceeds bounded scope")
    leaves = []
    for leaf in snapshot["leaves"]:
        leaves.append(dict(row=float(leaf["row"]+1), column=float(leaf["column"]+1),
            rows=float(leaf["rows"]), columns=float(leaf["columns"]), rank=float(leaf["rank"]),
            a1=np.asarray(leaf["a1"], dtype=np.float64).reshape(-1, 1).copy(),
            a2=np.asarray(leaf["a2"], dtype=np.float64).reshape(-1, 1).copy()))
    operator = dict(schema_version=2., assembly="normalized_leaves", dimension=float(snapshot["dimension"]),
        permutation=(np.asarray(snapshot["permutation"], dtype=np.float64)+1).reshape(-1, 1),
        scale=np.asarray(snapshot["scale"], dtype=np.float64).reshape(-1, 1).copy(),
        leaves=leaves, operator_kind="demag", charge_map=data["charge_map"].tocsc(copy=True))
    return dict(schema="radia.hdiv-host-setup.v1", backend="gpu", order=float(order),
        mesh_sha256=mesh_hash, label_contract_sha256=contract_hash, operator=operator,
        mass=data["mass"].tocsc(), rhs_operator=data["rhs_operator"], mean_operator=data["mean_operator"],
        evaluation=data["evaluation"].tocsc(), weights=data["weights"].reshape(-1, 1),
        quadrature_order=float(data["quadrature_order"]))
