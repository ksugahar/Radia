import hashlib,json,os,re,shutil
from collections import Counter
from pathlib import Path
from .optuna_boundary import matlab_optuna_mcp_route
from .optuna_quality import matlab_optuna_health
HERE=Path(__file__).parent; EXT=HERE/"extensions/radia-matlab-tools.json"; MATLAB=HERE/"matlab"
PROFILES={"existing":("existing","desktop"),"auto_nodesktop":("auto","nodesktop"),"new_nodesktop":("new","nodesktop")}
def matlab_extension_path():
    if not EXT.is_file(): raise FileNotFoundError(EXT)
    return EXT
def matlab_extension_contract():
    raw=matlab_extension_path().read_bytes(); data=json.loads(raw); tools=data.get("tools",[]); sig=data.get("signatures",{}); files=list((MATLAB/"+radia_mcp_matlab").glob("*.m")); names={p.stem for p in files}; errors=[]
    if len(tools)!=43 or len(sig)!=43: errors.append("expected 43 tools/signatures")
    if len(files)!=86: errors.append("expected 86 MATLAB functions")
    if any(not t["name"].startswith("matlab_") for t in tools): errors.append("invalid tool prefix")
    for p in files:
        text=p.read_text(encoding="utf-8")
        if "acoustic" in text.lower() or "fembem" in text.lower(): errors.append(f"{p.name}: namespace leak")
        errors += [f"{p.name}: missing {d}" for d in re.findall(r"radia_mcp_matlab\.([A-Za-z]\w*)",text) if d not in names]
    return {"schema":"radia-mcp.matlab-extension/v1","ok":not errors,"status":"ok" if not errors else "error","runtime_owner":"MathWorks MATLAB MCP Server","extension_file":str(EXT),"matlab_root":str(MATLAB),"sha256":hashlib.sha256(raw).hexdigest(),"tool_count":len(tools),"tool_names":[t["name"] for t in tools],"signature_count":len(sig),"matlab_function_count":len(files),"errors":errors}
def matlab_official_server_config(profile="existing",*,include_generic_extension=False):
    if profile not in PROFILES: raise ValueError(profile)
    session,display=PROFILES[profile]; args=[f"--matlab-session-mode={session}"]; files=[]; setup=""
    if include_generic_extension:
        c=matlab_extension_contract()
        if not c["ok"]: raise RuntimeError(c["errors"])
        files=[c["extension_file"]]; args += [f"--extension-file={c['extension_file']}"]; setup=f"addpath('{c['matlab_root']}');"
    if display=="nodesktop": args.append("--matlab-display-mode=nodesktop")
    candidates=[Path.home()/".matlab/agentic-toolkits/bin/matlab-mcp-server-windows-x64.exe"]
    candidates += [Path(x) for x in [shutil.which("matlab-mcp-server"),shutil.which("matlab-mcp-server-windows-x64.exe"),shutil.which("matlab-mcp-core-server"),shutil.which("matlab-mcp-core-server-win64.exe")] if x]
    cmd=os.getenv("RADIA_MATLAB_MCP_SERVER") or next((str(x) for x in candidates if x.is_file()),"matlab-mcp-server")
    return {"schema":"radia-mcp.matlab-server-config/v1","status":"ok","runtime_owner":"MathWorks MATLAB MCP Server","integration_owner":"radia-mcp.matlab","command_id":"matlab-mcp-server","command":cmd,"profile":profile,"args":args,"extension_files":files,"matlab_setup_code":setup}
def matlab_radia_acoustic_interface_contract(): return {"runtime_owner":"MathWorks MATLAB MCP Server","matlab_workflow_owner":"MathWorks MATLAB Agentic Toolkit","simulink_workflow_owner":"MathWorks Simulink Agentic Toolkit","generic_operations_owner":"radia_mcp.matlab","generic_matlab_package":"radia_mcp_matlab","production_owner":"radia.acoustics","education_solver_owner":"radia_mcp.acoustic_fembem","radia_extension_scope":["Radia MATLAB/MEX domain APIs","LTspice orchestration and result import","MATLAB Optuna 4.9.0 differential subset"]}


def _radia_repo_root():
    configured = os.environ.get("RADIA_REPO_ROOT", "").strip()
    candidates = [Path(configured)] if configured else []
    candidates.extend(HERE.resolve().parents)
    for candidate in candidates:
        source = candidate / "src" / "matlab" / "radia_mex.cpp"
        if source.is_file():
            return candidate, source
    return None, None


def _mex_command_branches():
    """Read the separated Radia and Optuna command lists from shared source."""
    _root, source = _radia_repo_root()
    if source is None:
        return [], []
    text = source.read_text(encoding="utf-8")
    match = re.search(
        r"mxArray\* Commands\(\)\s*\{(?P<body>.*?)constexpr std::size_t count",
        text,
        flags=re.DOTALL,
    )
    if match is None:
        return [], []
    body = match.group("body")
    conditional = re.search(
        r"#ifdef\s+RADIA_OPTUNA_MEX_ONLY(?P<optuna>.*?)"
        r"#else(?P<radia>.*?)#endif",
        body,
        flags=re.DOTALL,
    )
    if conditional is None:
        return re.findall(r'"([^"\\]*(?:\\.[^"\\]*)*)"', body), []
    extract = lambda text: re.findall(r'"([^"\\]*(?:\\.[^"\\]*)*)"', text)
    return extract(conditional.group("radia")), extract(conditional.group("optuna"))


def _radia_mex_commands():
    """Read commands compiled into the general Radia MEX gateway."""
    return _mex_command_branches()[0]


def _optuna_mex_commands():
    """Read commands compiled only into the lightweight Optuna gateway."""
    return _mex_command_branches()[1]


def _optuna_native_kernels(optuna_commands):
    """Return the optimizer kernels, excluding the api.* gateway commands.

    Selected by namespace rather than by a positional slice: a slice silently
    reports the wrong kernel list the moment another non-kernel command is
    added ahead of them.
    """
    return [name for name in optuna_commands if name.startswith("optuna.")]


def _pybind_public_names():
    root, _source = _radia_repo_root()
    if root is None:
        return []
    source = root / "src" / "lib" / "radia_pybind.cpp"
    if not source.is_file():
        return []
    text = source.read_text(encoding="utf-8")
    return [
        name for name in re.findall(r'\bm\.def\s*\(\s*"([^"]+)"', text)
        if not name.startswith("_")
    ]


def _pybind_all_top_level_names():
    root, _source = _radia_repo_root()
    if root is None:
        return []
    source = root / "src" / "lib" / "radia_pybind.cpp"
    if not source.is_file():
        return []
    text = source.read_text(encoding="utf-8")
    return re.findall(r'(?m)^\s{4}m\.def\s*\(\s*"([^"]+)"', text)


def _pybind_class_surface():
    """Return class/member names declared by the Radia pybind translation unit."""
    root, _source = _radia_repo_root()
    if root is None:
        return []
    source = root / "src" / "lib" / "radia_pybind.cpp"
    if not source.is_file():
        return []
    text = source.read_text(encoding="utf-8")
    starts = list(re.finditer(r'(?m)^\s{4}(?=py::class_|m\.def\s*\()', text))
    surface = []
    for index, start in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        chunk = text[start.start():end]
        if not chunk.lstrip().startswith("py::class_"):
            continue
        match = re.search(r'\bm\s*,\s*"([^"]+)"', chunk)
        if match is None:
            continue
        class_name = match.group(1)
        surface.append(f"{class_name}.__binding__")
        if re.search(r'\.def\s*\(\s*py::init', chunk):
            surface.append(f"{class_name}.__init__")
        for pattern in (
            r'\.def\s*\(\s*"([^"]+)"',
            r'\.def_static\s*\(\s*"([^"]+)"',
            r'\.def_property(?:_readonly)?\s*\(\s*"([^"]+)"',
            r'\.def_read(?:only|write)?\s*\(\s*"([^"]+)"',
        ):
            surface.extend(
                f"{class_name}.{member}" for member in re.findall(pattern, chunk)
            )
    return sorted(set(surface))


_PYBIND_INTERNAL_NUMERICAL_COMMANDS = {
    "_AcousticSoftSphere": ("acoustic.soft_sphere",),
    "_AcousticRigidSphere": ("acoustic.rigid_sphere",),
    "_AcousticFluidSphere": ("acoustic.fluid_sphere",),
    "_AcousticElasticSphere": ("acoustic.elastic_sphere",),
    "_AcousticSoftSphereComplexK": ("acoustic.soft_sphere_complex_k",),
    "_AcousticBDFDelta": ("acoustic.bdf_delta",),
    "_AcousticCQGrid": ("acoustic.cq_grid",),
    "_EVRSTMethodAlgebra": ("evrs.tmethod",),
    "_HybridVIMSchurComplement": ("hybrid_vim.schur",),
    "_HybridVIMSolve": ("hybrid_vim.solve",),
    "_SkinImpedance": ("hybrid_vim.skin_impedance",),
    "_SIBCAdmittanceTail": ("hybrid_vim.sibc_admittance_tail",),
    "_SIBCSchurTerminationImpedance": ("hybrid_vim.sibc_termination_impedance",),
    "_SIBCSchurTerminationAdmittance": ("hybrid_vim.sibc_termination_admittance",),
    "_TetHCurlReducedGram": ("hcurl.tet_reduced_gram",),
    "_AffineCellSelfEnergyShapeDerivative": ("hdiv.affine_cell_self_energy_shape_derivative",),
    "_AssembleSLDL_Galerkin": ("bem.assemble_sldl",),
    "_AssembleSLDL_Galerkin_P2": ("bem.assemble_sldl_p2",),
    "_HFromSegmentsComplex": ("biot_savart.h_segments_complex",),
    "_AFromSegmentsComplex": ("biot_savart.a_segments_complex",),
    "_AFromTrianglesComplex": ("biot_savart.a_triangles_complex",),
    "_BFromTrianglesComplex": ("biot_savart.b_triangles_complex",),
    "_EquivalenceSourceStaticH": ("equivalence.static_h",),
    "_EquivalenceSourceHarmonic": ("equivalence.harmonic",),
    "_average_B_in_box": ("radia.AverageBInBox",),
    "_average_demag_tensor": ("radia.AverageDemagTensor",),
    "_evaluate_ngsolve_coefficient_at_point": (
        "ngsolve.coefficient_function.evaluate",
    ),
    "_stream_aca_tsvd": ("stream.aca_tsvd",),
}

_PYBIND_INTERNAL_EXCLUSIONS = {
    "_compute_ho_bnd_nodes": "NGSolve-owned mesh/high-order boundary plumbing",
    "_volume_element_vertex_counts": "NGSolve-owned mesh inspection plumbing",
    "_TestPEECHACApKSanity": "Python regression-test helper, not a user API",
}

_RETIRED_UNSAFE_CONSTRUCTORS = (
    "radia.ObjMltExtPgn",
    "radia.ObjMltExtRtg",
    "radia.ObjMltExtTri",
)

_RETIRED_UNSAFE_C_ABI_SYMBOLS = (
    "RadObjMltExtPgn",
    "RadObjMltExtRtg",
    "RadObjMltExtTri",
)

_PYBIND_CLASS_EXCLUSIONS = {
    "_PEECBuilderInternal.__binding__": "implementation holder",
    "_HACApKPEECManagerInternal.__binding__": "implementation holder",
    "_HACApKBEMManagerInternal.__binding__": "implementation holder",
    "_HDivVectorPotentialCoefficient.__binding__": (
        "private radia.vim adapter that materializes an NGSolve-owned Python "
        "CoefficientFunction; MATLAB uses the declared vim-public Python "
        "boundary for construction and the projected HCurl .sol beam MEX "
        "contract for native tracking"),
    "_HDivVectorPotentialCoefficient.__init__": (
        "private radia.vim adapter that materializes an NGSolve-owned Python "
        "CoefficientFunction; MATLAB uses the declared vim-public Python "
        "boundary for construction and the projected HCurl .sol beam MEX "
        "contract for native tracking"),
    "_HDivVectorPotentialCoefficient.source_count": (
        "diagnostic on the private NGSolve CoefficientFunction adapter; source "
        "construction remains in the declared vim-public Python boundary"),
    "_HDivExactVectorPotentialCoefficient.__binding__": (
        "private radia.vim adapter for the analytic equivalent-current A "
        "source (construction=\"exact\"); it materializes an NGSolve-owned "
        "Python CoefficientFunction, so MATLAB uses the declared vim-public "
        "Python boundary for construction and the projected HCurl .sol beam "
        "MEX contract for native tracking"),
    "_HDivExactVectorPotentialCoefficient.__init__": (
        "private radia.vim adapter for the analytic equivalent-current A "
        "source; same vim-public Python / projected-HCurl MEX boundary as "
        "_HDivVectorPotentialCoefficient"),
    "_HDivExactVectorPotentialCoefficient.tetrahedron_count": (
        "diagnostic on the private analytic equivalent-current adapter; the "
        "source statistics are already reported through the vim-public "
        "vector_potential_coefficient_stats result entry"),
    "_HDivExactVectorPotentialCoefficient.triangle_count": (
        "diagnostic on the private analytic equivalent-current adapter; the "
        "source statistics are already reported through the vim-public "
        "vector_potential_coefficient_stats result entry"),
    "_HDivFieldCoefficient.reflection_normal": (
        "introspection-only mirror of the FieldCoefficientFromSolution "
        "reflection_normal argument on the private NGSolve adapter; the "
        "physical field MATLAB consumes flows through the projected .sol "
        "beam MEX contract"),
    "_HDivFieldCoefficient.reflection_symmetrized": (
        "introspection-only flag on the private NGSolve adapter; same "
        "projected-.sol MEX boundary as reflection_normal"),
    "_ChargeGramHMatrix.charge_sigma": (
        "sigma-normalization diagnostic for the roundoff-amplification "
        "regression tests (empty before BuildHMatrix; every public apply "
        "already wraps S back, so MATLAB sees the physical Gram)"),
    "_ChargeGramHMatrix._reduce_configured_candidate_directional_schur": (
        "private proposal-only acceleration used by the Python topology "
        "optimizer; it is not a public pybind capability, and MATLAB uses "
        "the public reduce_configured_candidate_schur MEX contract"),
}


def _class_commands(prefix, members, info_command):
    return {f"{prefix}.{member}": (info_command,) for member in members}


_PYBIND_CLASS_COMMANDS = {
    "_EnergyStopMaterial.__binding__": ("energy_stop.create",),
    "_EnergyStopMaterial.__init__": ("energy_stop.create",),
    "_EnergyStopMaterial.state0": ("energy_stop.state0",),
    "_EnergyStopMaterial.forward_batch": ("energy_stop.forward",),
    "_EnergyStopMaterial.commit_batch": ("energy_stop.commit",),
    "_EnergyStopMaterial.stored_energy_batch": ("energy_stop.stored_energy",),
    **_class_commands("_EnergyStopMaterial", (
        "branch_count", "state_size", "alpha", "b_max", "nu_bound", "eta", "gamma"),
        "energy_stop.info"),
    "_HDivFieldEvaluator.__binding__": ("hdiv.field_evaluator.from_tet",),
    "_HDivFieldEvaluator.from_tet": ("hdiv.field_evaluator.from_tet",),
    "_HDivFieldEvaluator.from_cloud": ("hdiv.field_evaluator.from_cloud",),
    "_HDivFieldEvaluator.field": ("hdiv.field_evaluator.field",),
    "_HDivFieldEvaluator.field_gradient": ("hdiv.field_evaluator.field_gradient",),
    "_HDivFieldEvaluator.candidate_algorithm_for": ("hdiv.field_evaluator.candidate_algorithm",),
    "_HDivFieldEvaluator.last_algorithm": ("hdiv.field_evaluator.last_algorithm",),
    "_HDivFieldEvaluator.stats": ("hdiv.field_evaluator.stats",),
    "_PlanarHDivFieldEvaluator.__binding__": ("hdiv.planar_evaluator.create",),
    "_PlanarHDivFieldEvaluator.field": ("hdiv.planar_evaluator.field",),
    "_PlanarHDivFieldEvaluator.az": ("hdiv.planar_evaluator.az",),
    "_PlanarHDivFieldEvaluator.stats": ("hdiv.planar_evaluator.stats",),
    "_HDivFieldCoefficient.__binding__": ("hdiv.field_evaluator.as_coefficient",),
    "_HDivFieldCoefficient.__init__": ("hdiv.field_evaluator.as_coefficient",),
    "_HDivFieldCoefficient.algorithm": ("ngsolve.coefficient_function.info",),
    "_PlanarHDivFieldCoefficient.__binding__": ("hdiv.planar_evaluator.as_coefficient",),
    "_PlanarHDivFieldCoefficient.__init__": ("hdiv.planar_evaluator.as_coefficient",),
    "_PlanarHDivFieldCoefficient.source_angle": ("ngsolve.coefficient_function.info",),
    "_PlanarHDivFieldCoefficient.target_angle": ("ngsolve.coefficient_function.info",),
    "_HDivDemagMatrix.__binding__": ("hacapk.charge_gram.demag_matrix",),
    "_ComplexDiagonalInverseMatrix.__binding__": ("ngsolve.matrix.diagonal_preconditioner",),
    "_ProjectedBaseMatrix.__binding__": ("ngsolve.matrix.projected_create",),
    "_ProjectedBaseMatrix.__init__": ("ngsolve.matrix.projected_create",),
    "_ReducedBlockMatrix.__binding__": ("ngsolve.matrix.reduced_block_create",),
    "_ReducedBlockMatrix.__init__": ("ngsolve.matrix.reduced_block_create",),
    "_ReducedBlockMatrix.diagonal_preconditioner": ("ngsolve.matrix.diagonal_preconditioner",),
    "_ReducedBlockMatrix.term_count": ("ngsolve.matrix.term_count",),
    "_ChargeGramHMatrix.__binding__": ("hacapk.charge_gram.create_monopole",),
    "_ChargeGramHMatrix.__init__": (
        "hacapk.charge_gram.create_monopole", "hacapk.charge_gram.create_analytic_tet",
        "hacapk.charge_gram.create_analytic_polytope", "hacapk.charge_gram.create_high_order_tet",
        "hacapk.charge_gram.create_curved_high_order_tet", "hacapk.charge_gram.create_curved_polytope",
        "hacapk.charge_gram.create_hex", "hacapk.charge_gram.create_wedge",
        "hacapk.charge_gram.create_planar_2d"),
    "_ChargeGramHMatrix.ndof": ("hacapk.charge_gram.info",),
    "_ChargeGramHMatrix.set_image_rotations": ("hacapk.charge_gram.set_image_rotations",),
    "_ChargeGramHMatrix.build_hmatrix": ("hacapk.charge_gram.build",),
    "_ChargeGramHMatrix.matvec": ("hacapk.charge_gram.matvec",),
    "_ChargeGramHMatrix.matvec_transpose": ("hacapk.charge_gram.matvec_transpose",),
    "_ChargeGramHMatrix.matvec_sym": ("hacapk.charge_gram.matvec_sym",),
    "_ChargeGramHMatrix.entry": ("hacapk.charge_gram.entry",),
    "_ChargeGramHMatrix.hex_volume_self_block_directional_derivative": ("hacapk.charge_gram.hex_volume_self_block_directional_derivative",),
    "_ChargeGramHMatrix.hex_face_self_block_directional_derivative": ("hacapk.charge_gram.hex_face_self_block_directional_derivative",),
    "_ChargeGramHMatrix.hex_charge_gram_directional_derivative": ("hacapk.charge_gram.hex_directional_derivative",),
    "_ChargeGramHMatrix.tet_volume_self_block_directional_derivative": ("hacapk.charge_gram.tet_volume_self_block_directional_derivative",),
    "_ChargeGramHMatrix.tet_face_self_block_directional_derivative": ("hacapk.charge_gram.tet_face_self_block_directional_derivative",),
    "_ChargeGramHMatrix.tet_charge_gram_directional_derivative": ("hacapk.charge_gram.tet_directional_derivative",),
    "_ChargeGramHMatrix.tet_charge_map_row_directional_rates": ("hacapk.charge_gram.tet_charge_map_row_directional_rates",),
    "_ChargeGramHMatrix.wedge_volume_self_block_directional_derivative": ("hacapk.charge_gram.wedge_volume_self_block_directional_derivative",),
    "_ChargeGramHMatrix.wedge_face_self_block_directional_derivative": ("hacapk.charge_gram.wedge_face_self_block_directional_derivative",),
    "_ChargeGramHMatrix.wedge_charge_gram_directional_derivative": ("hacapk.charge_gram.wedge_directional_derivative",),
    "_ChargeGramHMatrix.directional_derivative_operator": ("hacapk.charge_gram.directional_derivative_operator",),
    "_ChargeGramHMatrix.directional_derivative_contractions": ("hacapk.charge_gram.directional_derivative_contractions",),
    "_ChargeGramHMatrix.directional_derivative_contractions_many": ("hacapk.charge_gram.directional_derivative_contractions_many",),
    "_ChargeGramHMatrix.hex_state_check": ("hacapk.charge_gram.hex_state_check",),
    "_ChargeGramHMatrix.hex_stored_nodes": ("hacapk.charge_gram.hex_stored_nodes",),
    "_ChargeGramHMatrix.hex_state_breakdown": ("hacapk.charge_gram.hex_state_breakdown",),
    "_ChargeGramHMatrix.configure_charge_map": ("hacapk.charge_gram.configure_charge_map",),
    "_ChargeGramHMatrix.configure_vector_charge_map": ("hacapk.charge_gram.configure_vector_charge_map",),
    "_ChargeGramHMatrix.configure_mass_matrix": ("hacapk.charge_gram.configure_mass_matrix",),
    "_ChargeGramHMatrix.configure_mass_matrix_ngsolve": ("hacapk.charge_gram.configure_mass_matrix_ngsolve",),
    "_ChargeGramHMatrix.configure_geometry_mass_matrix": ("hacapk.charge_gram.configure_geometry_mass_matrix",),
    "_ChargeGramHMatrix.configure_geometry_mass_matrix_ngsolve": ("hacapk.charge_gram.configure_geometry_mass_matrix_ngsolve",),
    "_ChargeGramHMatrix.restore_geometry_mass_matrix": ("hacapk.charge_gram.restore_geometry_mass_matrix",),
    "_ChargeGramHMatrix.set_configured_constraints": ("hacapk.charge_gram.set_configured_constraints",),
    "_ChargeGramHMatrix.configured_active_hmatrix_stats": ("hacapk.charge_gram.configured_active_hmatrix_stats",),
    "_ChargeGramHMatrix.demag_matrix": ("hacapk.charge_gram.demag_matrix",),
    "_ChargeGramHMatrix.apply_configured_demag": ("hacapk.charge_gram.demag_apply",),
    "_ChargeGramHMatrix.apply_configured_geometry_mass": ("hacapk.charge_gram.geometry_mass_apply",),
    "_ChargeGramHMatrix.apply_configured_linear_material_operator": ("hacapk.charge_gram.apply_configured_linear_material_operator",),
    "_ChargeGramHMatrix.apply_configured_linear_material_operator_many": ("hacapk.charge_gram.apply_configured_linear_material_operator_many",),
    "_ChargeGramHMatrix.configured_linear_material_element_blocks": ("hacapk.charge_gram.configured_linear_material_element_blocks",),
    "_ChargeGramHMatrix.configured_linear_material_candidate_clusters": ("hacapk.charge_gram.configured_linear_material_candidate_clusters",),
    "_ChargeGramHMatrix.reduce_configured_candidate_schur": ("hacapk.charge_gram.reduce_configured_candidate_schur",),
    "_ChargeGramHMatrix._reduce_configured_candidate_directional_schur": ("hacapk.charge_gram.reduce_configured_candidate_directional_schur",),
    "_ChargeGramHMatrix.apply_configured_mass_riesz": ("hacapk.charge_gram.mass_riesz",),
    "_ChargeGramHMatrix.solve_configured_linear_material_mass_riesz": ("hacapk.charge_gram.solve_configured_linear_material",),
    "_ChargeGramHMatrix.solve_configured_linear_material_auto_prec": ("hacapk.charge_gram.solve_configured_linear_material_auto_prec",),
    "_ChargeGramHMatrix.solve_configured_linear_material_auto_prec_many": ("hacapk.charge_gram.solve_configured_linear_material_auto_prec_many",),
    "_ChargeGramHMatrix.configured_field_functional_rows": ("hacapk.charge_gram.configured_field_functional_rows",),
    "_ChargeGramHMatrix.configured_field_functional_rows_directional_derivative": ("hacapk.charge_gram.configured_field_functional_rows_directional_derivative",),
    "_ChargeGramHMatrix.configured_field_values_shape_derivative": ("hacapk.charge_gram.configured_field_values_shape_derivative",),
    "_ChargeGramHMatrix.create_field_evaluator": ("hacapk.charge_gram.create_field_evaluator",),
    "_ChargeGramHMatrix.create_planar_field_evaluator": ("hacapk.charge_gram.create_planar_field_evaluator",),
    "_ChargeGramHMatrix.stats": ("hacapk.charge_gram.stats",),
    "_ChargeGramHMatrix.from_sampled_laplace": ("hacapk.charge_gram.create_sampled_laplace",),
    "_ChargeGramHMatrix.from_sampled_planar_log": ("hacapk.charge_gram.create_sampled_planar_log",),
    "_ChargeGramHMatrix.from_local_polynomials": ("hacapk.charge_gram.create_local_polynomials",),
    "_ChargeGramHMatrix.operator_configured": ("hacapk.charge_gram.operator_info",),
    "_ChargeGramHMatrix.constraint_count": ("hacapk.charge_gram.operator_info",),
    "HACApKChargeGramDerivative.__binding__": ("hacapk.charge_gram.directional_derivative_operator",),
    "HACApKChargeGramDerivative.ndof": ("hacapk.charge_gram_derivative.info",),
    "HACApKChargeGramDerivative.entry": ("hacapk.charge_gram_derivative.entry",),
    "HACApKChargeGramDerivative.matvec_sym": ("hacapk.charge_gram_derivative.matvec_sym",),
    "HACApKChargeGramDerivative.stats": ("hacapk.charge_gram_derivative.info",),
    "RadiaField.__binding__": ("ngsolve.radia_field.create",),
    "RadiaField.__init__": ("ngsolve.radia_field.create",),
    "RadiaField.PrepareCache": ("ngsolve.radia_field.prepare_cache",),
    "RadiaField.ClearCache": ("ngsolve.radia_field.clear_cache",),
    "RadiaField.GetCacheStats": ("ngsolve.radia_field.cache_stats",),
    "RadiaField.as_voxel_cf": ("ngsolve.radia_field.as_voxel_coefficient",),
    **_class_commands("RadiaField", (
        "radia_obj", "field_type", "use_transform", "precision"),
        "ngsolve.radia_field.info"),
    "HACApKPEECManager.__binding__": ("hacapk.peec.create",),
    "HACApKPEECManager.__init__": ("hacapk.peec.create",),
    "HACApKPEECManager.BuildHMatrix": ("hacapk.peec.build",),
    "HACApKPEECManager.MatVec": ("hacapk.peec.matvec",),
    **_class_commands("HACApKPEECManager", ("GetNDOF", "IsValid", "GetStats"),
        "hacapk.peec.info"),
    "HACApKBEMManager.__binding__": ("hacapk.bem.create",),
    "HACApKBEMManager.__init__": ("hacapk.bem.create",),
    "HACApKBEMManager.BuildHMatrix": ("hacapk.bem.build",),
    "HACApKBEMManager.MatVec": ("hacapk.bem.matvec",),
    **_class_commands("HACApKBEMManager", ("GetNDOF", "IsValid", "GetStats"),
        "hacapk.bem.info"),
}


def _pybind_command_name(name):
    # The beam Lie-map / reference-orbit surface groups its MEX commands
    # under beam.* rather than the default radia.<name>, so the audit needs
    # the mapping explicitly. Each pair below shares one C++ kernel
    # (rad_lie::LieMapTensorsFromSpolyArrays, DragtFinnFactorizeFourthOrder,
    # ApplyDragtFinnMapBatch, rad_orbit::TrackReferenceOrbit3D, and
    # rad_orbit::TrackReferenceOrbit3DToPlane), which is what makes the
    # coverage claim true rather than name-level.
    beam_names = {
        "lie_map_tensors_spoly": "beam.lie.map_tensors_spoly",
        "lie_dragt_finn_factorize": "beam.lie.dragt_finn_factorize",
        "lie_apply_dragt_finn_batch": "beam.lie.apply_dragt_finn_batch",
        "track_reference_orbit_native": "beam.orbit.track_reference_3d",
        "track_reference_orbit_to_plane_native": "beam.orbit.track_reference_to_plane",
    }
    if name in beam_names:
        return beam_names[name]
    hlu_names = {
        "HLUSetTruncTol": "hlu.set_trunc_tol",
        "HLUGetTruncTol": "hlu.get_trunc_tol",
        "HLULastTimings": "hlu.last_timings",
        "HLUMaterializeStats": "hlu.materialize_stats",
        "HLUSetParallel": "hlu.set_parallel",
        "HLUGetParallel": "hlu.get_parallel",
        "HLUSetParCutoff": "hlu.set_par_cutoff",
        "HLUMaxThreads": "hlu.max_threads",
        "HLUSetAccumCap": "hlu.set_accum_cap",
        "HLUGetAccumCap": "hlu.get_accum_cap",
        "HLUMixedBreakdown": "hlu.mixed_breakdown",
        "HLUSelfTest": "hlu.self_test",
        "HLUSelfTestRk": "hlu.self_test_rk",
        "HLUSelfTestAddmulRkRk": "hlu.self_test_addmul_rkrk",
        "HLUSelfTestRadiaExactWithMatrix": "hlu.self_test_radia_exact_with_matrix",
        "HLUSelfTestRadiaExactDiag": "hlu.self_test_radia_exact_diag",
        "HLUSelfTestRadiaExact": "hlu.self_test_radia_exact",
        "HLUSelfTestDepth3Asymmetric": "hlu.self_test_depth3_asymmetric",
        "HLUSelfTestMixedSiblingViaConversion": "hlu.self_test_mixed_sibling_via_conversion",
        "HLUSelfTestMixedSiblingNonUniform": "hlu.self_test_mixed_sibling_nonuniform",
        "HLUSelfTestMixedSibling": "hlu.self_test_mixed_sibling",
        "HLUSelfTestRkDeep": "hlu.self_test_rk_deep",
    }
    return hlu_names.get(name, f"radia.{name}")


def _radia_mex_group(command):
    prefix = command.split(".", 1)[0]
    if prefix == "radia":
        return "radia-core"
    if prefix == "hlu":
        return "hacapk-hlu"
    return prefix


def matlab_radia_mex_contract(topic="all"):
    """Return the shared Python/NGSolve/MATLAB MEX capability contract."""
    topic = str(topic).lower().strip()
    allowed = {"all", "mex", "ngsolve", "optimization", "simulink", "limitations"}
    if topic not in allowed:
        raise ValueError(f"unknown topic: {topic}; available: {', '.join(sorted(allowed))}")

    repo_root, source = _radia_repo_root()
    commands = _radia_mex_commands()
    optuna_commands = _optuna_mex_commands()
    pybind_names = _pybind_public_names()
    pybind_command_names = [_pybind_command_name(name) for name in pybind_names]
    pybind_missing = sorted(set(pybind_command_names) - set(commands))
    all_top_level_names = _pybind_all_top_level_names()
    internal_names = sorted(name for name in all_top_level_names if name.startswith("_"))
    internal_unclassified = sorted(
        set(internal_names)
        - set(_PYBIND_INTERNAL_NUMERICAL_COMMANDS)
        - set(_PYBIND_INTERNAL_EXCLUSIONS)
    )
    internal_missing = sorted({
        command
        for name in internal_names
        for command in _PYBIND_INTERNAL_NUMERICAL_COMMANDS.get(name, ())
        if command not in commands
    })
    class_surface = _pybind_class_surface()
    class_relevant = sorted(
        name for name in class_surface if name not in _PYBIND_CLASS_EXCLUSIONS
    )
    class_unmapped = sorted(
        name for name in class_relevant if name not in _PYBIND_CLASS_COMMANDS
    )
    class_missing = sorted({
        command
        for name in class_relevant
        for command in _PYBIND_CLASS_COMMANDS.get(name, ())
        if command not in commands
    })
    class_covered = sum(
        name in _PYBIND_CLASS_COMMANDS
        and all(command in commands for command in _PYBIND_CLASS_COMMANDS[name])
        for name in class_relevant
    )
    retired_unsafe_constructor_leaks = sorted(
        name for name in _RETIRED_UNSAFE_CONSTRUCTORS
        if name in commands or name.removeprefix("radia.") in pybind_names
    )
    retired_unsafe_c_abi_leaks = []
    if repo_root is not None:
        for relative_path in (
            Path("src/lib/radentry.h"),
            Path("src/lib/radentry.cpp"),
            Path("src/lib/raddll.def"),
        ):
            candidate = repo_root / relative_path
            if not candidate.is_file():
                continue
            candidate_text = candidate.read_text(encoding="utf-8", errors="ignore")
            retired_unsafe_c_abi_leaks.extend(
                f"{relative_path.as_posix()}:{symbol}"
                for symbol in _RETIRED_UNSAFE_C_ABI_SYMBOLS
                if re.search(rf"\b{re.escape(symbol)}\b", candidate_text)
            )
    retired_unsafe_c_abi_leaks.sort()
    parity_complete = not (
        pybind_missing or internal_missing or internal_unclassified
        or class_missing or class_unmapped or retired_unsafe_constructor_leaks
        or retired_unsafe_c_abi_leaks
    )
    radia_matlab_root = (repo_root / "matlab") if repo_root is not None else None
    radia_package = (radia_matlab_root / "+radia") if radia_matlab_root is not None else None
    matlab_functions = sorted(path.stem for path in radia_package.rglob("*.m")) if radia_package is not None and radia_package.is_dir() else []
    optuna_package = (radia_package / "+optuna") if radia_package is not None else None
    optuna_files = sorted(optuna_package.glob("*.m")) if optuna_package is not None and optuna_package.is_dir() else []
    optuna_classes = sorted(
        path.stem for path in optuna_files
        if "classdef" in path.read_text(encoding="utf-8", errors="ignore")
    )
    optuna_functions = sorted(path.stem for path in optuna_files if path.stem not in optuna_classes)
    optuna_health = matlab_optuna_health(str(repo_root) if repo_root is not None else "")
    groups = Counter(_radia_mex_group(command) for command in commands)
    mex_ready = bool(commands and source is not None)
    base = {
        "schema": "radia-mcp.matlab-radia-mex/v1",
        "status": "ready" if mex_ready else "source_unavailable",
        "parity_status": (
            "complete_for_radia_pybind_numerics"
            if parity_complete else "incomplete"
        ),
        "source_of_truth": "src/matlab/radia_mex.cpp: Commands()",
        "mex_entrypoint": "radia_mex",
        "command_count": len(commands),
        "optuna_mex_entrypoint": "optuna_mex",
        "optuna_mex_command_count": len(optuna_commands),
        "optuna_mex_command_names": optuna_commands,
        "command_groups": dict(sorted(groups.items())),
        "command_names": commands,
        "pybind_public_count": len(pybind_names),
        "pybind_covered_count": len(pybind_command_names) - len(pybind_missing),
        "pybind_missing": pybind_missing,
        "retired_unsafe_constructors": list(_RETIRED_UNSAFE_CONSTRUCTORS),
        "retired_unsafe_constructor_leaks": retired_unsafe_constructor_leaks,
        "retired_unsafe_c_abi_symbols": list(_RETIRED_UNSAFE_C_ABI_SYMBOLS),
        "retired_unsafe_c_abi_leaks": retired_unsafe_c_abi_leaks,
        "pybind_internal_numerical_count": len(_PYBIND_INTERNAL_NUMERICAL_COMMANDS),
        "pybind_internal_missing": internal_missing,
        "pybind_internal_unclassified": internal_unclassified,
        "pybind_internal_exclusions": dict(sorted(_PYBIND_INTERNAL_EXCLUSIONS.items())),
        "pybind_class_surface_count": len(class_relevant),
        "pybind_class_covered_count": class_covered,
        "pybind_class_missing_commands": class_missing,
        "pybind_class_unmapped": class_unmapped,
        "pybind_class_exclusions": dict(sorted(_PYBIND_CLASS_EXCLUSIONS.items())),
        "matlab_wrapper_count": len(matlab_functions),
        "matlab_optuna_class_count": len(optuna_classes),
        "matlab_optuna_function_count": len(optuna_functions),
        "matlab_optuna_classes": optuna_classes,
        "matlab_optuna_file_count": optuna_health.get("distribution", {}).get("matlab_file_count"),
        "matlab_optuna_expected_file_count": optuna_health.get("distribution", {}).get("expected_matlab_file_count"),
        "matlab_optuna_public_api": optuna_health.get("public_api"),
        "matlab_optuna_distribution_health": optuna_health,
        "interfaces": {
            "python": "radia._radia_pybind and radia Python modules",
            "matlab": "radia_mex owns Radia/NGSolve operations; optuna_mex independently owns optimization kernels and has no NGSolve, MKL, Radia-core, or Python dependency",
            "ngsolve": "NGSolve owns meshes/spaces/transforms; MATLAB uses numeric/struct contracts",
            "mcp": ["matlab_radia_mex_contract", "matlab_optuna_health", "matlab_optuna_oracle_plan", "matlab_optuna_benchmark_plan", "matlab_optuna_release_gate", "matlab_optuna_simulink_contract", "radia_usage", "ngsolve_usage"],
        },
        "ngsolve_boundary": {
            "project_bridge_commands": [name for name in commands if name.startswith("ngsolve.") or name.startswith("hcurl.") or name.startswith("hdiv.")],
            "numeric_field_and_matrix_boundary": ["ngsolve.space_info", "ngsolve.matrix_dump", "ngsolve.mesh.*", "ngsolve.fespace.*", "ngsolve.bilinear_form.*", "ngsolve.matrix.*", "ngsolve.grid_function.from_fespace", "ngsolve.linear_form.*", "ngsolve.coefficient_function.*", "ngsolve.grid_function.*", "ngsolve.vector.*", "hdiv.field_evaluator.*", "hdiv.planar_evaluator.*"],
            "canonical_matlab_names": ["radia.ngsolve.space_info", "radia.ngsolve.matrix_dump", "radia.ngsolve.Mesh", "radia.ngsolve.FESpace", "radia.ngsolve.BilinearForm", "radia.ngsolve.Matrix", "radia.ngsolve.LinearForm", "radia.ngsolve.CoefficientFunction", "radia.ngsolve.GridFunction", "radia.ngsolve.Vector", "radia.ngsolve.hcurl_eddy_cln_native_basis", "radia.ngsolve.hcurl_eddy_cln_model", "radia.hcurl.tet_reduced_gram"],
            "full_external_ngsolve_api": "not duplicated; NGSolve remains the FE object owner",
            "topology_deformation": "radia.ngsolve.Mesh setDeformation/unsetDeformation/trafoQuality use native VectorH1 and GetTrafo; MATLAB receives per-element Jacobian ratios and spectral condition numbers",
        },
        "verified_contract": {
            "matlab_gate": "runtests('tests/matlab')",
            "python_numerical_parity_gate": "runtests('tests/matlab/test_radia_ngsolve_parity.m')",
            "acoustic_python_mex_parity_gate": "runtests('tests/matlab/test_acoustic_mex.m')",
            "axifem_python_mex_parity_gate": "runtests('tests/matlab/test_axifem_mex.m')",
            "hcurl_topology_python_mex_parity_gate": "runtests('tests/matlab/test_hcurl_topology_optimization.m')",
            "topology_two_level_gate": "runtests('tests/matlab/test_topology_optimization.m')",
            "optuna_native_kernel_gate": "runtests('tests/matlab/test_radia_mex.m', Name={'test_radia_mex/testOptunaGatewaySeparation','test_radia_mex/testOptunaParetoRankCrowding','test_radia_mex/testOptunaParzenLogPdfKernels'})",
            "optuna_table_gate": "runtests('tests/matlab/test_optuna_table.m')",
            "optuna_simulink_block_gate": "runtests('tests/matlab/test_optuna_simulink_block.m')",
            "optuna_distribution_health_gate": "matlab_optuna_health",
            "optuna_differential_oracle_gate": "matlab_optuna_oracle_plan",
            "optuna_performance_gate": "matlab_optuna_benchmark_plan -> matlab_optuna_release_gate",
            "optuna_native_kernel_benchmark": "validation_test/optimization/results_matlab_optuna_mex_benchmark_20260806.json",
            "optuna49_performance_benchmark": "validation_test/optimization/results_matlab_optuna49_performance_20260825.json",
            "runtime_probe": "radia.quickCheck()",
            "native_build": "pwsh -ExecutionPolicy Bypass -File .\\Build.ps1 -MatlabMexOnly -Verbose",
            "optuna_native_build": "pwsh -ExecutionPolicy Bypass -File .\\Build.ps1 -OptunaMexOnly",
            "native_motor_family_artifact": "validation_test/radia_mcp/artifacts/annular_motor_dual_lane_v1/native_motor_angle_family.json",
            "native_motor_family_gate": "run('validation_test/radia_mcp/generate_motor_angle_family_mex_artifact.m')",
            "openmp_runtime_policy": "radia.setup excludes foreign libiomp5md.dll directories from the MATLAB process PATH",
        },
        "limitations": [
            "Python object identity is not shared with MATLAB; handles and numeric/struct snapshots are the boundary.",
            "ObjMltExtPgn, ObjMltExtRtg, and ObjMltExtTri were deleted from Python, MATLAB, and the legacy C ABI; no compatibility shims remain.",
            "Mesh-backed high-order NGSolve spaces remain NGSolve-owned; MATLAB receives assembled data or diagnostics. The native HCurl CLN builder now projects M, curl-curl, and ports in C++ and returns a local diffusion model, but it does not replace topology-aware VIM/BEM/SIBC assembly.",
            "Python convenience modules are tracked separately by matlab/python_api_parity_manifest.json; native MEX coverage alone is not module-level parity.",
            "MATLAB Optuna shared behavior is judged only against pinned optuna==4.9.0 fixtures; MATLAB-only table/MAT, Simulink, parallel, and MEX behavior is extension evidence, not upstream-parity evidence.",
            "The native MEX target does not embed or launch Python; Python callback objects are rejected at the MEX boundary and numeric/handle equivalents are used instead. Full Python-DLL independence requires an NGSolve/Netgen build without Python support.",
            "beam.orbit.track_reference_3d shares the rad_orbit::TrackReferenceOrbit3D kernel with the pybind route but drives Radia-object sources only; the HDiv iron evaluator stays a pybind-owned handle until an evaluator handle exists in the MEX registry.",
        ],
    }
    sections = {
        "mex": {"command_count": len(commands), "command_groups": dict(sorted(groups.items())), "command_names": commands, "optuna_mex_command_count": len(optuna_commands), "optuna_mex_command_names": optuna_commands, "pybind_missing": pybind_missing, "pybind_internal_missing": internal_missing, "pybind_class_missing_commands": class_missing, "pybind_class_unmapped": class_unmapped},
        "ngsolve": {"owner": "NGSolve", "matlab_boundary": "ngsolve.space_info, ngsolve.matrix_dump, persistent Mesh/FESpace/BilinearForm/Matrix/LinearForm handles, native CoefficientFunction/GridFunction handles, native HCurl response projection, and numeric HDiv field evaluators", "policy": "Use NGSolve for FE plumbing; exchange typed native handles, numeric matrices, vectors, fields, and metadata. Persistent forms expose built-in real/complex volume integrators, scalar CoefficientFunction-weighted volume and trace matrices, native CoefficientFunction volume and boundary RHS assembly in real/complex H1/HCurl/HDiv, and native sparse matvec/inverse operations; arbitrary callbacks and tensor-valued forms remain explicit gaps."},
        "optimization": {"package": "radia.optuna", "native_gateway": "optuna_mex", "native_gateway_required": True, "native_command_count": len(optuna_commands), "classes": optuna_classes, "factory_functions": optuna_functions, "storage": "MAT-file containing readable MATLAB tables and CAE trial/failure artifacts", "native_kernels": _optuna_native_kernels(optuna_commands), "mcp_tool": "matlab_optuna_simulink_contract"},
        "simulink": {"class": "radia.optuna.SimulinkRunner", "workflow": "SimulationInput -> sim/parsim -> score/constraints/validation/artifacts -> Study.tell or typed failure", "blocks": ["Radia/Applications/Induction Heating: readable Level-2 MATLAB Eddy and Thermal S-Functions backed by radia_mex handles", "radia.simulink.buildTeam28CLNModel", "radia.simulink.buildHCurlEddyCLNModel Block=radia-mex", "radia.simulink.buildMotorAngleFamilyModel", "Optimization/Optuna Optimization: start/cancel, failure telemetry, and automatic optimizer MEX kernels"], "native_state_space_commands": ["simulink.state_space.create", "simulink.state_space.info", "simulink.state_space.output", "simulink.state_space.update", "simulink.state_space.step", "simulink.state_space.snapshot", "simulink.state_space.restore", "simulink.state_space.reset", "simulink.state_space.destroy"], "native_state_space_overloads": {"static": "create(A,B,C,D,x0); output(handle,u); update(handle,u)", "periodic_motor_family": "create(grid,period,A,B,C,D,Q,R,S,x0); output(handle,mechanical_angle,u) -> [linear_outputs; torque]; update(handle,mechanical_angle,u)", "standalone_debug": "step is the atomic output-plus-update probe", "sim_state": "snapshot/restore preserve native state for Simulink CustomSimState"}, "mcp_tool": "matlab_optuna_simulink_contract"},
        "reinforcement_learning": {"package": "radia.rl", "workflow": "reset -> MEX/Radia step -> reward -> next observation", "adapter": "rlFunctionEnv when Reinforcement Learning Toolbox is installed"},
        "limitations": {"items": base["limitations"]},
    }
    base["topic"] = topic
    base["topic_data"] = {"sections": sections} if topic == "all" else sections[topic]
    return base


def matlab_optuna_simulink_contract():
    """Return the MATLAB-native Optuna and Simulink difference contract."""
    health = matlab_optuna_health()
    from .optuna_oracle import matlab_optuna_compatibility_contract

    upstream = matlab_optuna_compatibility_contract()
    return {
        "schema": "radia-mcp.matlab-optuna-simulink/v3",
        "status": "ready",
        "package": "radia.optuna",
        "distribution": "radia-optuna",
        "upstream_oracle": upstream,
        "upstream_oracle_version": "optuna==4.9.0",
        "mcp_ownership": matlab_optuna_mcp_route(),
        "distribution_health": health,
        "seed_semantics": {
            "explicit": "equal explicit uint32 seeds consume the checked NumPy RandomState-compatible stream",
            "unseeded": "seed=None draws fresh private entropy per sampler without mutating MATLAB's global RNG; exact proposal sequences are intentionally nondeterministic",
        },
        "classes": {
            "Study": "ask/tell, scalar or vector objectives, normalized ObjectiveTable, MAT persistence, bestTrial/bestValue/bestParams/bestSolution or paretoFront",
            "Trial": "suggestFloat/suggestInteger/suggestCategorical, report, shouldPrune, user attributes",
            "FixedTrial": "objective evaluation against supplied parameters with upstream-oracled suggestions, warnings, attributes, and compatibility errors",
            "TrialPruned": "throw(radia.optuna.TrialPruned()) is caught by Study.optimize and matches upstream PRUNED state, last intermediate value, and callback behavior",
            "RandomSampler": "explicitly seeded proposals match the checked Optuna 4.9.0 oracle",
            "TPESampler": "seeded scalar, mixed, grouped multivariate, constrained, callable gamma/weights, and categorical-distance Parzen proposals plus independent-fallback warnings match the checked Optuna 4.9.0 oracle",
            "MOTPESampler": "seeded multi-objective TPE and constrained Pareto behavior match the checked Optuna 4.9.0 oracle",
            "CmaEsSampler": "seeded numeric proposals and independently seeded fallback-sampler controls match Optuna 4.9.0; restart, margin, separable, source-trial warm-start, and learning-rate modes remain explicit gaps",
            "GPSampler": "Backend='upstream-python' delegates startup, constrained PyTorch/SciPy LogEI acquisition, and deterministic persisted-history replay to pinned Optuna 4.9.0; Backend='matlab-native' is integration-only",
            "NSGAIISampler": "seeded population behavior and all six built-in crossover proposal sequences match Optuna 4.9.0",
            "NSGAIIISampler": "seeded reference-line population behavior matches Optuna 4.9.0",
            "QMCSampler": "seeded scrambled and deterministic unscrambled Sobol/Halton proposals match Optuna/SciPy within the documented dimension boundary",
            "GridSampler": "finite Cartesian-grid exhaustion matches Optuna 4.9.0",
            "BruteForceSampler": "fixed and conditional define-by-run tree exhaustion matches Optuna 4.9.0",
            "PartialFixedSampler": "fixed and delegated proposal sequences match Optuna 4.9.0",
            "Pruners": "Percentile, threshold, patient, no-op, successive-halving, Hyperband, and Wilcoxon decisions are upstream-oracled",
            "Importance": "get_param_importances delegates fANOVA, mean-decrease-impurity, and PED-ANOVA evaluation to pinned Optuna 4.9.0",
            "Termination": "MaxTrialsCallback, BestValueStagnationEvaluator, Terminator, and TerminatorCallback provide upstream-oracled stopping behavior",
            "LiveMonitor": "trial progress, objective history, duration, parameter history, and live two-objective Pareto display",
            "SimulinkRunner": "CAE-aware SimulationInput configuration, adaptive-batch sim/parsim, c <= 0 constraints, validation, typed failures, artifact manifests, model SHA-256, and four-part timing",
            "LTspiceRunner": "serial or parfeval LTspice trials with isolated output directories, complex RAW, and client-owned Study updates",
            "SheetMetalRunner": "Optuna outer trials over the native HEX-sheet or HCurl sheet-metal two-level drivers, with isolated NGSolve/Cubit work directories and result.json provenance",
        },
        "simulink_blocks": [
            "Radia/Applications/Induction Heating uses readable Level-2 MATLAB S-Functions backed by separate native Eddy and Thermal object handles for the distributed-field kernels; Eddy receives current, angle, and temperature distribution, while Thermal receives heat distribution, ambient temperature, and angle.",
            "makeIHNativeConfig plus validateIHNativeConfig enforce native matrix dimensions, shared thermal CSR sparsity, positive cell weights, preassembled FEM/PEEC/BEM-A/BIM method metadata, and .vol preflight before initialization.",
            "IH Thermal publishes its previous accepted field in the Level-2 Outputs callback, feeds it directly to Eddy, and advances only in Update. Linear current changes use current-squared heat scaling; temperature-dependent operator data trigger a new field solve. Nonlinear BH fails fast in the preview. LUT and lumped IH builders are removed.",
            "buildTeam28CLNModel uses the validated 50 Hz six-stage CLN force LUT; makeTeam28CoilBuilderPlant and radia_team28_coilbuilder_dynamic.slx advance only the slower mechanical motion, so the cycle-averaged lane must not be called a full electromagnetic transient.",
            "ExportHCurlEddyCLNJSON plus loadHCurlEddyCLNModel maps numeric R/L/P reduced HCurl-VIM data into a passive discrete state-space block; solveHCurlEddyCLNHarmonic uses hybrid_vim.solve. radia.ngsolve.hcurl_eddy_cln_model additionally assembles a Python-free local HCurl diffusion projection directly from a VOL mesh.",
            "buildHCurlEddyCLNModel(..., Block=\"radia-mex\") uses the Python-free native state-space handle for fixed reduced HCurl models; makeMotorAngleFamily plus buildMotorAngleFamilyModel use the same persistent MEX handle for periodic angle interpolation and quadratic torque output.",
            "ExportHCurlEddyCLNFamilyJSON plus loadHCurlEddyCLNFamily/interpolateHCurlEddyCLNFamily provides a common-basis height-indexed family; buildHCurlEddyCLNFamilyModel exposes height, derivative-current, and coil-current ports.",
            "Optimization/Optuna Optimization advances one trial per block sample, accepts independent start/cancel signals, continues after recorded CAE failures, exports attempted/failed/failure-code telemetry, and uses the required independent optuna_mex gateway without Python per trial.",
            "Optimization/Sheet Metal Optimization evaluates a radia.optuna.SheetMetalRunner and uses the same best-update and Pareto monitor contract.",
        ],
        "team28": {
            "package": "radia.simulink",
            "workflow": "makeTeam28CLNLUT -> evaluateTeam28CLNForce -> buildTeam28CLNModel",
            "coilbuilder_dynamic_workflow": "makeTeam28CoilBuilderLUT -> makeTeam28CoilBuilderPlant -> radia_team28_coilbuilder_dynamic.slx",
            "frequency_hz": 50,
            "force_convention": "force_N is the physical signed force; upward_lift_N is positive upward.",
            "high_fidelity_boundary": "numeric p=6 HCurl-VIM and common-basis height-family export are supported; compute-host generation and comparison against the 25-point reference remain the validation gate.",
            "validated_dynamic_scope": "cycle_averaged_mechanical_motion",
            "electromagnetic_model_class": "fixed_frequency_cycle_averaged_force_height_lut",
            "height_coupling": "quasi_steady_interpolation",
            "electromagnetic_state_transient_included": False,
            "motional_emf_included": False,
            "damping_identified_from_measurement": False,
            "artifact_gate": "mcp-server-maglev:team28_cycle_averaged_motion_gate",
            "unsupported_claims": [
                "full_electromagnetic_transient",
                "carrier_resolved_electromagnetic_waveform",
                "motion_induced_emf",
                "experimentally_identified_damping",
            ],
            "full_transient_next_step": "Advance position-dependent electromagnetic states and the mechanical state together, including motion derivative terms of the R/L/P family, then close current-force-energy-motion histories.",
        },
        "hcurl_eddy_cln": {
            "package": "radia.simulink",
            "workflow": "NGSolve vim.ExportHCurlEddyCLNJSON -> MATLAB loadHCurlEddyCLNModel -> solveHCurlEddyCLNHarmonic/buildHCurlEddyCLNModel",
            "equation": "(R+sL)c=-sPi, with u=-di/dt and y=P'c in the time domain.",
            "mex_kernel": "hybrid_vim.solve",
            "force": "Optional K(k,a,b) force operator is evaluated by evaluateHCurlEddyCLNForce.",
            "moving_family": "ExportHCurlEddyCLNFamilyJSON -> loadHCurlEddyCLNFamily -> makeMotorAngleFamily -> buildMotorAngleFamilyModel; C++ performs periodic interpolation, state advance, linear outputs, and quadratic torque without Python or MATLAB matrix algebra per step.",
            "native_motor_angle_family": {
                "matlab_factory": "radia.simulink.makeMotorAngleFamily",
                "simulink_builder": "radia.simulink.buildMotorAngleFamilyModel",
                "s_function": "radia_motor_angle_family_mex_sfunction",
                "state": "persistent C++ handle with separated output/update, CustomSimState snapshot/restore, reset, and destroy lifecycle",
                "interpolation": "periodic linear interpolation over one mechanical-angle period",
                "torque": "0.5*x'*Q*x + x'*R*u + 0.5*u'*S*u",
                "validation_artifact": "validation_test/radia_mcp/artifacts/annular_motor_dual_lane_v1/native_motor_angle_family.json",
                "verified_tests": 74,
            },
            "sibc": "SIBC termination must be rationalized before state-space export; the current numeric bridge records this boundary explicitly.",
        },
        "tables": ["TrialTable", "ObjectiveTable", "ParamTable", "IntermediateTable", "UserAttrTable", "ConstraintTable", "SamplerStateTable"],
        "native_acceleration": {
            "status_api": "radia.optuna.nativeStatus",
            "gateway": "optuna_mex",
            "command_count": health.get("native", {}).get("command_count"),
            "required": True,
            "commands": [
                "optuna.pareto.rank_crowding",
                "optuna.parzen.log_pdf_numerical",
                "optuna.parzen.log_pdf_categorical",
                "optuna.sobol.points",
            ],
            "policy": "Require the exact independent 21-command optuna_mex gateway for checked optimizer kernels; keep Study/Trial, persistence, callbacks, and unsupported fused-option orchestration in readable MATLAB. A missing or incompatible optuna_mex is an error and never redirects through radia_mex or a silent substitute.",
            "missing_mex_fallback": False,
            "full_optimizer_in_cpp": False,
            "python_per_trial": "false for MATLAB-native samplers; true only for GPSampler Backend='upstream-python'",
            "upstream_python_gp_python_per_trial": True,
        },
        "cae_trial_contract": {
            "success_schema": "radia.optuna.cae-trial.v1",
            "failure_schema": "radia.optuna.cae-failure.v1",
            "identity": ["model path", "model SHA-256", "MATLAB version", "Simulink version", "geometry/mesh/material/excitation context"],
            "constraints": "c <= 0",
            "timing": ["configuration", "simulation", "postprocess", "total"],
            "failure_policy": "Classify mesh, convergence, license/resource, timeout, configuration, and observable failures; preserve the failed trial and continue when requested.",
        },
        "multi_objective": {
            "directions": "one minimize/maximize entry per returned objective",
            "storage": "ObjectiveTable has one row per trial and objective index; TrialTable.Value retains objective 1 for compatibility",
            "selection": "bestTrial, bestValue, bestParams, and bestSolution are single-objective only; paretoFront returns non-dominated trials",
            "visualization": "LiveMonitor plots objective history for one objective and a live Pareto scatter for two or more objectives",
            "samplers": ["RandomSampler", "MOTPESampler", "NSGAIISampler"],
        },
        "sampler_quality": {
            "intended_scale": "expensive CAE trials where field-solver time dominates MATLAB table-backed ask/tell overhead",
            "tpe_boundary": "seeded scalar, mixed, multivariate, constrained, and callable gamma/weights behavior is checked against Optuna 4.9.0; group decomposition and categorical-distance hooks remain gaps",
            "correlated_continuous": "CmaEsSampler has checked seeded numeric and independent-fallback parity; GPSampler Backend='upstream-python' owns exact checked LogEI behavior",
            "multi_objective": "use MOTPESampler or NSGAIISampler and inspect front error plus coverage, not Pareto point count alone",
            "validation": "validation_test/optimization/validate_matlab_optuna_quality.m",
            "python_parity_claim": "Only behavior mapped to pinned optuna==4.9.0 fixtures is parity evidence; API coverage must close before a complete-compatibility claim.",
            "simulink_auto": "CmaEsSampler for one objective and NSGAIISampler for multiple objectives; the mask also exposes every sampler explicitly",
        },
        "parallel_trials": {
            "simulink": "SimulinkRunner.optimizeParallel asks adaptive client-side batches and evaluates each SimulationInput batch with parsim so completed trials influence later proposals",
            "ltspice": "LTspiceRunner.optimizeParallel suggests parameters on the client and runs isolated LTspice jobs with parfeval",
            "state_safety": "Study, sampler, and table mutations stay on the client; workers return numerical results only",
            "requirement": "Parallel Computing Toolbox; absence is an explicit error, never a silent serial fallback",
        },
        "ltspice_integrated_workflow": {
            "raw": "readRaw/RawRead normalize ASCII or binary, real or complex results under radia.ltspice.raw.v2 while preserving source schema",
            "parallel": "LTspiceRunner.optimizeParallel uses isolated parfeval runs and returns numerical results to the client-owned Study",
            "pareto": "vector ScoreFcn -> ObjectiveTable -> MOTPE/NSGA-II -> LiveMonitor Pareto scatter",
        },
        "cad_topology_optimization": {
            "cad_owner": "Cubit design cells, element IDs, material blocks, and final CAD reconstruction",
            "physics": "Radia-VIM analytic system linearization A*dm_i=db_i-dA_i*m",
            "sensitivity_policy": "No cell-wise finite differences: factor/solve the VIM tangent system with all design perturbations as multiple right-hand sides",
            "subproblem": "bounded linear programming with volume, field inequality, and move-limit constraints",
            "matlab_api": ["radia.topopt.linearizeVIM", "radia.topopt.solveLPUpdate", "radia.topopt.optimizeVIMLP", "radia.topopt.writeCubitJournal"],
            "abe_element_fill": {
                "status": "native-mex-ready",
                "mex_command": "topopt.abe_element_fill_plan",
                "matlab_api": [
                    "radia.topopt.solveAbeElementFillPlan",
                    "radia.topopt.contractHDivElementFillResponse",
                    "radia.topopt.composeSpecificationFillResponse",
                    "radia.topopt.binElementFillToInterfaceHeight",
                    "radia.topopt.blendedInterfaceDisplacement",
                ],
                "variables": "one signed fill fraction per design element; existing iron [-1,0], addable air [0,1]",
                "factorization": "HACApK ACA+ followed by QR-TSVD, factored once and reused through bounded residual corrections",
                "ngsolve_boundary": "local/global HDiv DOF transforms, GetTrafo deformation, complete field re-solves, and exact transfer-map acceptance remain caller-owned",
            },
            "adjoint_optimization": {
                "status": "ready",
                "contract": "Evaluator returns scalar objective, design gradient, c <= 0 inequalities with design-by-constraint Jacobian, and optional SQP equalities/Jacobian",
                "matlab_api": [
                    "radia.topopt.checkAdjointGradient",
                    "radia.topopt.optimizeAdjoint",
                    "radia.topopt.optimizeHCurlActivationAdjoint",
                ],
                "mma": "Finite-bound moving-asymptote iterations with analytic convex-subproblem Jacobians; inequalities only",
                "sqp": "fmincon SQP with supplied objective and constraint gradients; inequalities and equalities",
                "finite_difference_policy": "Directional finite differences are explicit QA only and are never an optimizer fallback",
                "hybrid_role": "TPE/CMA-ES outer global or discrete search followed by MMA/SQP continuous field-sensitive refinement",
            },
            "mcp_tool": "matlab_cad_topology_build",
            "implemented_surface": "mesh-cell material optimization and deterministic Cubit solid/void block journal",
            "cad_reconstruction_route": "density field -> level set -> Coreform Cubit Sculpt -> Exodus mesh import",
            "cubit_validation_gates": [
                "cubit_ato_levelset_sculpt_source_replay_gate",
                "cubit_levelset_sculpt_hex_validation_gate",
            ],
            "boundary": "The generic journal classifies an existing hex mesh; density-to-level-set/Sculpt reconstruction is design-domain-specific and must pass the Cubit gates before being called CAD-complete",
            "sheet_metal": {
                "variables": ["neutral-surface normal displacement", "thickness", "material activation"],
                "constraints": ["linearized volume", "thickness bounds", "discrete curvature", "move limits", "user field inequalities"],
                "routing": "Jacobian determinant, condition number, relative displacement, and topology change select ngsolve_deform, ngsolve_refine, or cubit_rebuild",
                "python_executor": "radia.sheet_metal_optimization.apply_ngsolve_mesh_route uses Mesh.SetDeformation or SetRefinementFlag+Refine; Cubit rebuild remains an explicit external gated action",
                "matlab_api": ["radia.topopt.localTrustRegion", "radia.topopt.solveSheetMetalLP", "radia.topopt.backtrackTrafoDeformation", "radia.topopt.optimizeHexSheetTopology", "radia.topopt.CubitHexRemeshBackend", "radia.topopt.routeMeshUpdate"],
                "optuna_runner": "radia.optuna.SheetMetalRunner",
                "simulink_block": "Optimization/Sheet Metal Optimization",
                "simulink_builder": "radia.simulink.buildSheetMetalOptimizationBlock",
                "trial_artifact": "radia.optuna.sheet-metal-trial.v1 result.json in an isolated trial directory",
                "mcp_tool": "matlab_sheet_metal_topology_build",
                "two_level_loop": {"inner": "5-20 NGSolve GetTrafo iterations with continuous activation and no H-matrix rebuild", "outer": "0.35/0.65 activation hysteresis accumulates pending topology changes until their fraction or age reaches the Cubit batch threshold", "post_cubit": "rebuild the H-matrix exactly once, then start the next inner batch"},
                "continuity_policy": "Activation stays continuous through the complete inner batch; threshold classification is committed only at the outer Cubit boundary",
                "analytic_shape_tangent": {
                    "operator": "A=inv_chi*M+B^T*G*B differentiated by the full product rule",
                    "laplace_pairs": "distance and quadrature-weight derivatives are analytic; no optimization finite differences",
                    "matlab_api": ["radia.topopt.linearizeLaplacePairGram", "radia.topopt.linearizeVIMOperator", "radia.topopt.linearizeVIM"],
                    "self_panel_boundary": "singular self-panel values and their shape derivatives must come from the topology-specific analytic moment kernel; the generic pair helper intentionally leaves its diagonal zero",
                },
                "hcurl_eddy_bubble": {
                    "status": "native-mex-ready",
                    "numerical_owner": "C++ MEX with NGSolve-owned HCurl assembly and HACApK scalar Gram",
                    "matlab_api": [
                        "radia.topopt.HCurlTopologyOperator",
                        "radia.topopt.assembleHCurlResistanceShapeTangents",
                        "radia.topopt.assembleHCurlCellCurlGrams",
                        "radia.topopt.sampleHCurlSubtetVelocities",
                        "radia.topopt.linearizeHCurlMultifrequencyJoule",
                        "radia.topopt.linearizeHCurlActivationMultifrequencyJoule",
                        "radia.topopt.optimizeAdjoint",
                        "radia.topopt.optimizeHCurlActivationAdjoint",
                        "radia.topopt.linearizeAndSolveHCurlSheetJouleLP",
                        "radia.topopt.linearizeAndSolveHCurlActivationSheetJouleLP",
                        "radia.topopt.optimizeHCurlEddyBubbleHexSheet",
                        "radia.topopt.optimizeHCurlEddyBubbleActivationHexSheet",
                    ],
                    "mex_commands": [
                        "hcurl.topopt.operator.*",
                        "hcurl.topopt.resistance_shape_tangents",
                        "hcurl.topopt.cell_curl_grams",
                        "hcurl.topopt.multifrequency_joule",
                        "hcurl.topopt.activation_multifrequency_joule",
                    ],
                    "shape_gradient": "analytic HCurl Piola dR plus HACApK dG and dB contractions; no directional dL matrix is materialized",
                    "activation_gradient": "cell-local NGSolve curl Grams plus SIMP conductivity and matrix-free HACApK activation contractions",
                    "load_cases": "weighted complex adjoints over multiple frequencies and excitations in one MEX call",
                    "python_boundary": "none in the MATLAB optimization loop",
                },
            },
        },
        "persistence": "MAT-file via Study.save or AutoSave; tables remain inspectable in MATLAB, while SimulinkRunner records versioned CAE trial/failure structs and artifact manifests. bestSolution derives a reloadable best snapshot from the tables.",
        "simulink_workflow": [
            "Create a radia.optuna.Study.",
            "Configure each Trial into Simulink.SimulationInput.",
            "Run sim/parsim and return finite objectives, c <= 0 constraints, validation, artifacts, and timing through the configured callbacks.",
            "Record COMPLETE, PRUNED, or typed FAIL state and persist both readable tables and versioned CAE evidence.",
        ],
        "reinforcement_learning_workflow": [
            "Create radia.rl.Environment for a native numerical step contract, or place a Simulink RL Agent around a production application block.",
            "Call reset and step around Radia/NGSolve MEX numerical kernels.",
            "Optionally adapt the environment to rlFunctionEnv.",
            "For IH, train against the native distributed Eddy/Thermal block; do not substitute a LUT or lumped thermal environment.",
        ],
        "python_relation": "The official optuna/optuna-mcp server owns every shared operation in its live tools/list. radia-mcp owns only MATLAB/Simulink differences and never performs a Python-to-MATLAB-Engine call per trial; exact GPSampler Backend='upstream-python' intentionally uses pinned Optuna 4.9.0 in-process.",
        "mcp_route": "matlab_optuna_mcp_route",
    }


def matlab_simulink_library_contract():
    """Return the installation and compatibility contract for Radia blocks."""
    return {
        "schema": "radia-mcp.matlab-simulink-library/v3",
        "status": "ready",
        "library": "radia_simulink_library",
        "browser_name": "Radia",
        "registration_code": [
            "radia.simulink.buildLibrary",
            "sl_refresh_customizations",
        ],
        "registration_timing": "Run once after installation or after replacing the packaged library.",
        "blocks": [
            "Applications/Electromagnet",
            "Applications/PCB PEEC",
            "Applications/Motor",
            "Applications/Stream Function",
            "Applications/Induction Heating",
            "Applications/Magnetic Levitation",
            "Applications/Field Study",
            "Material Models/Material Dictionary",
            "Coupling/Winding Dictionary",
            "Coupling/Field Study Configuration",
            "LTspice/LTspice Circuit",
            "LTspice/Hysteretic LTspice Plant",
            "Optimization/Optuna Optimization",
            "Optimization/Sheet Metal Optimization",
            "Optimization/Optuna Monitor",
        ],
        "applications": {
            "runner": "radia.simulink.application",
            "matlab_entry_point": "radia.simulink.runApplication",
            "config_writer": "radia.simulink.writeApplicationConfig",
            "config_schema": "radia.simulink.application_config.v1",
            "result_schema": "radia.simulink.application_run.v1",
            "trigger": "one batch solve on a Boolean rising edge",
            "initial_backend": "python-headless-cli",
            "per_step_python": "forbidden",
            "mex_policy": (
                "optional only after numerical parity, error propagation, "
                "state lifecycle, repeated-run, and long-run stability tests"
            ),
            "notebook_policy": (
                "notebook workbenches are retired for every application, "
                "including IH; docs notebooks are result-bearing examples "
                "whose field scenes use Draw(field, mesh, name=..., ...)"
            ),
            "spatial_artifact": {
                "format": "GMSH .msh v4.1",
                "location": "application run directory",
                "required_when": "the selected mode computes a spatial field",
                "result_index": "radia_result.artifacts.gmsh",
                "scalar_only": "not-applicable",
                "purpose": "visualization and post-processing only",
            },
            "mesh_preflight": {
                "input_format": "Netgen .vol",
                "checker": "check-vol",
                "timing": "after export and before solver initialization",
                "required_when": "the selected mode consumes a .vol mesh",
                "label_contract_schema": "radia.vol-label-contract.v1",
                "report_schema": "cubit-mesh-export.vol-check.v1",
                "material_constants": (
                    "validated by DesignSpec/configuration and never inferred "
                    "from labels"
                ),
            },
        },
        "material_dictionary": {
            "setup_type": "MATLAB dictionary with scalar-cell material structs",
            "material_constructor": "radia.simulink.makeMaterialSpec",
            "compiler": "radia.simulink.compileMaterialDictionary",
            "mesh_inventory": "radia.simulink.inspectVolMaterials",
            "mesh_format": "Netgen .vol",
            "region_policy": "exact named coverage; no guessed material constants",
            "runtime_bus": "RadiaMaterialBus",
            "runtime_shape": "fixed-width numeric and logical fields",
            "per_step_dictionary_lookup": False,
            "per_step_strings": False,
            "per_step_python": False,
            "supported_properties": [
                "linear and tabulated nonlinear magnetic",
                "conductivity",
                "relative permittivity",
                "remanence vector",
                "density, specific heat, and thermal conductivity",
                "B-input energy hysteresis parameters",
            ],
        },
        "field_study": {
            "application_block": "Applications/Field Study",
            "configuration_block": "Coupling/Field Study Configuration",
            "study_constructor": "radia.simulink.makeFieldStudySpec",
            "study_compiler": "radia.simulink.compileFieldStudy",
            "request_writer": "radia.simulink.writeFieldStudyRequest",
            "runtime_bus": "RadiaStudyBus",
            "mesh_format": "Netgen .vol",
            "physics": [
                "electrostatic",
                "current_flow (DC or harmonic lossy dielectric)",
                "steady_heat",
                "harmonic_eddy",
            ],
            "formulations": ["planar", "axisymmetric"],
            "harmonic_eddy_operator": "(K + j*omega*M_sigma) a = S i",
            "harmonic_eddy_gates": [
                "frequency_hz > 0",
                "linear non-hysteretic permeability",
                "one nonzero current per compiled winding identity",
                "branch real power closes 0.5*omega^2*a^H*M_sigma*a",
            ],
            "scalar_gates": [
                "complete material-region coverage",
                "Dirichlet or positive thermal Robin constraint",
                "terminal reaction and energy/power balance",
            ],
            "execution": "one owned Python/NGSolve batch worker per rising trigger",
            "per_step_python": False,
            "artifacts": ["versioned JSON", "timing breakdown", "Gmsh .msh v4.1", ".geo/.opt companions"],
            "retirement_claim": (
                "FEMM physics are exposed through the default Simulink interface; "
                "full retirement still requires frozen live solver artifacts for every lane"
            ),
        },
        "circuit_field_coupling": {
            "native_backend": "exact-ZOH reduced field/circuit MEX S-Function",
            "detailed_circuit_backend": "LTspice interval coupling",
            "shared_identity": [
                "winding terminal names",
                "turn count and polarity",
                "series/parallel connection",
                "voltage/current sign convention",
                "mechanical position and speed state",
            ],
            "purpose": "control, power circuit, motor, loss, and thermal co-simulation",
            "winding_constructor": "radia.simulink.makeWindingSpec",
            "winding_compiler": "radia.simulink.compileWindingDictionary",
            "winding_bus": "RadiaWindingBus",
            "coil_side_sign": "one RegionPolarity value per .vol winding region",
            "command_bus": "RadiaMachineCommandBus",
            "response_bus": "RadiaMachineResponseBus",
            "mechanical_owner": "Simulink or Simscape",
            "motion_coordinates": [
                "rotor angle and angular speed",
                "3-D translation position and velocity",
                "load torque and force",
            ],
            "field_outputs": [
                "terminal current",
                "flux linkage",
                "back EMF",
                "electromagnetic torque and force",
                "copper, iron, and eddy-current losses",
            ],
        },
        "ltspice": {
            "supported_distribution": "Current Analog Devices LTspice only",
            "executable": "LTspice.exe",
            "simulink_runtime": "Level-2 MATLAB S-Function: radia_ltspice_sfun",
            "native_boundary": "Standalone MEX function ABI only when a measured hot path requires it",
            "automatic_locations": [
                r"C:\Program Files\ADI\LTspice\LTspice.exe",
                r"%LOCALAPPDATA%\Programs\ADI\LTspice\LTspice.exe",
            ],
            "explicit_override": "Pass Executable=... or set the block executable parameter.",
            "legacy_ltc_versions": "not supported",
            "browser_dependency": "Registration does not require LTspice; circuit execution does.",
        },
        "kicad_ltspice": {
            "source_of_truth": "KiCad .kicad_sch remains the editable PCB-design source; generated CIR/ASC files are analysis artifacts.",
            "workflow": [
                "radia.kicad.exportSpiceNetlist",
                "radia.kicad.prepareLTspice",
                "radia.kicad.buildLTspiceBlock",
            ],
            "pipeline": ".kicad_sch -> kicad-cli SPICE CIR -> editable LTspice ASC and/or Simulink LTspice Circuit block",
            "reverse_sync": "LTspice ASC edits are not written back to KiCad schematic or PCB metadata.",
        },
        "optimization": {
            "execution": "one complete Optuna trial per Simulink sample after a rising start trigger",
            "monitor": "standard Simulink Scope and XY Graph; no web browser",
            "best_update": "one-sample pulse only when the primary-objective incumbent changes",
            "pareto": "fixed-size X/Y arrays, active point count, and revision counter",
            "sheet_metal_runner": "radia.optuna.SheetMetalRunner over native MATLAB/MEX plus NGSolve/Cubit drivers",
            "quality_validation": "validation_test/optimization/validate_matlab_optuna_quality.m",
            "sampler_mask": "auto, random, tpe, cmaes, motpe, or nsgaii; auto selects CMA-ES for one objective and NSGA-II for multiple objectives",
        },
        "execute_with": "official MATLAB MCP evaluate_matlab_code",
    }
