"""Fast structural checks for the public, result-bearing notebook corpus."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_docs_notebooks_are_parseable_and_free_of_replacement_glyphs():
    notebooks = sorted((ROOT / "docs").rglob("*.ipynb"))
    assert notebooks

    for path in notebooks:
        data = json.loads(path.read_text(encoding="utf-8"))
        assert data.get("nbformat", 0) >= 4, path
        assert data.get("cells"), path
        assert "\ufffd" not in json.dumps(data, ensure_ascii=False), path


def test_public_example_notebooks_keep_their_saved_scene_contract():
    for path in sorted((ROOT / "docs").rglob("*.ipynb")):
        data = json.loads(path.read_text(encoding="utf-8"))
        radia = data.get("metadata", {}).get("radia", {})
        if radia.get("notebook_role") != "example":
            continue

        assert radia.get("webgui_required") is True, path
        if radia.get("webgui_field_required"):
            assert radia.get("webgui_required") is True, path
        assert any(
            cell.get("outputs")
            for cell in data.get("cells", [])
            if cell.get("cell_type") == "code"
        ), path


def test_docs_do_not_own_validation_or_benchmark_json():
    offenders = [
        path.relative_to(ROOT)
        for path in (ROOT / "docs").rglob("*.json")
        if "validation" in path.name.lower() or "benchmark" in path.name.lower()
    ]
    assert not offenders, (
        "Move validation and benchmark JSON records to validation_test/: "
        + ", ".join(map(str, offenders))
    )


def test_urn_docs_use_the_packaged_implementation():
    urn_docs = ROOT / "docs" / "universal_relaxation_network"
    duplicate_modules = {
        urn_docs / "relaxation_basis_library.py",
        urn_docs / "universal_relaxation_network.py",
    }
    assert not any(path.exists() for path in duplicate_modules)

    stale_imports = (
        "from universal_relaxation_network import",
        "from relaxation_basis_library import",
    )
    offenders = [
        path.relative_to(ROOT)
        for path in urn_docs.rglob("*.py")
        if any(token in path.read_text(encoding="utf-8") for token in stale_imports)
    ]
    assert not offenders, (
        "URN docs must import the canonical radia.urn package: "
        + ", ".join(map(str, offenders))
    )


def test_urn_docs_do_not_own_generated_validation_artifacts():
    urn_docs = ROOT / "docs" / "universal_relaxation_network"
    offenders = [
        path.relative_to(ROOT)
        for directory in (urn_docs / "results", urn_docs / "validation_output")
        if directory.exists()
        for path in directory.rglob("*")
        if path.is_file()
    ]
    assert not offenders, (
        "URN runtime and validation artifacts belong under validation_test/: "
        + ", ".join(map(str, offenders))
    )


def test_docs_do_not_track_notebook_checksum_sidecars():
    offenders = []
    for path in (ROOT / "docs").rglob("*_result.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if payload.get("schema") == "radia.notebook_result.v1":
            offenders.append(path.relative_to(ROOT))
    assert not offenders, (
        "Executed notebooks own their saved outputs; remove checksum sidecars: "
        + ", ".join(map(str, offenders))
    )


def test_docs_indexes_use_demo_or_workflow_names():
    offenders = sorted(
        path.relative_to(ROOT)
        for name in ("examples.md", "EXAMPLES.md")
        for path in (ROOT / "docs").rglob(name)
    )
    assert not offenders, (
        "The examples tier is retired; name public indexes demos or workflows: "
        + ", ".join(map(str, offenders))
    )


def test_docs_do_not_restore_completed_migration_ledgers():
    retired = {
        ROOT / "docs" / "clebsch_hodograph": {
            "examples_catalog",
            "examples_migration",
        },
        ROOT / "docs" / "induction_heating": {
            "induction_heating_examples_catalog",
            "public_demo",
        },
        ROOT / "docs" / "kelvin": {
            "kelvin_classic_demos",
            "kelvin_examples_migration",
        },
        ROOT / "docs" / "peec_integration": {
            "cleanup_routing",
            "examples_catalog",
            "post_examples_migration",
            "public_demo",
            "verification_migration",
        },
        ROOT / "docs" / "stream_function": {
            "demo_gallery",
            "examples_catalog",
        },
    }
    offenders = [
        path.relative_to(ROOT)
        for directory, retired_stems in retired.items()
        for path in directory.iterdir()
        if path.is_file()
        and path.suffix.lower() in {".ipynb", ".json"}
        and path.stem.removesuffix("_results") in retired_stems
    ]
    assert not offenders, (
        "Git history owns completed migration bookkeeping; remove docs ledgers: "
        + ", ".join(map(str, offenders))
    )


def test_docs_do_not_restore_retired_result_bookkeeping():
    retired = {
        ROOT / "docs" / "axifem" / "axifem_element_evidence.json",
        ROOT / "docs" / "cubit_mesh_export" / "cubit_mesh_export_showcase_results.json",
        ROOT / "docs" / "cubit_mesh_export" / "netgen" / "p_convergence_demo_results.json",
        ROOT / "docs" / "electric_machine" / "cogging_skew_demo_results.json",
        ROOT / "docs" / "clebsch_hodograph" / "excitation_invariant_field.json",
        ROOT / "docs" / "clebsch_hodograph" / "hodograph_bending_sy.json",
        ROOT / "docs" / "clebsch_hodograph" / "hodograph_feasibility_2d.json",
        ROOT / "docs" / "clebsch_hodograph" / "edge_focusing_tracking.json",
        ROOT / "docs" / "clebsch_hodograph" / "edge_focusing_fem_results.json",
        ROOT / "docs" / "gmsh_post" / "em_fieldlines_results.json",
        ROOT / "docs" / "gmsh_post" / "em_particle_orbits_results.json",
        ROOT / "docs" / "gmsh_post" / "em_post_gallery_results.json",
        ROOT / "docs" / "nodal_force" / "nodal_force_results.json",
        *(ROOT / "docs" / "mesh_fusion" / f"results_phase{i}.json" for i in range(1, 6)),
        ROOT / "docs" / "gmsh_animation" / "gmsh_animation_export_results.json",
        ROOT / "docs" / "gmsh_animation" / "gmsh_animation_results.json",
        ROOT / "docs" / "kelvin" / "kelvin_exterior_source_and_aphi_results.json",
        ROOT / "docs" / "hdiv_vim" / "compare_curved_vs_radia_field.json",
        ROOT / "docs" / "hdiv_vim" / "hdiv_curved_nonlinear_field.json",
        ROOT / "docs" / "hdiv_vim" / "hdiv_demag_curved.json",
        ROOT / "docs" / "hdiv_vim" / "c_type_three_formulation_tosca_mixed_iron.vol-check.json",
        ROOT / "docs" / "section_optics" / "section_optics_design_results.json",
        ROOT / "docs" / "section_optics" / "stamp_notebook_hash.py",
        ROOT / "docs" / "universal_relaxation_network" / "cq_urn_bridge_results.json",
        ROOT / "docs" / "universal_relaxation_network" / "urn_vs_vf_comparison.json",
    }
    retired.update((ROOT / "docs" / "clebsch_hodograph" / "demos").glob("*.json"))
    retired.update((ROOT / "docs" / "maglev" / "demos").rglob("*.json"))
    offenders = sorted(path.relative_to(ROOT) for path in retired if path.exists())
    assert not offenders, (
        "Docs notebooks embed display results; validation_test owns checked evidence: "
        + ", ".join(map(str, offenders))
    )


def test_urn_validation_drivers_stay_out_of_docs():
    docs_dir = ROOT / "docs" / "universal_relaxation_network"
    validation_only = {
        "ablation_study.py",
        "benchmark_urn_vs_skrf_vf.py",
        "run_ltspice_verification.py",
        "urn_benchmark_improved.py",
        "urn_benchmark_suite.py",
        "urn_vs_vf_comparison.py",
        "urn_vs_vf_quick.py",
        "validate_all_datasets.py",
        "validate_real_data.py",
        "validate_real_world.py",
        "validate_tdk_ferrite.py",
        "validate_urn_vs_vf.py",
        "verify_timedomain_stability.py",
        "vf_benchmark_same_data.py",
    }
    offenders = sorted(path.name for path in docs_dir.iterdir() if path.name in validation_only)
    assert not offenders, (
        "Move URN validation and benchmark drivers to validation_test: "
        + ", ".join(offenders)
    )
