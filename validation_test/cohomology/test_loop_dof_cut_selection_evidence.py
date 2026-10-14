import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "loop_dof_cut_selection_data.json"


def test_loop_dof_cut_selection_evidence_is_topologically_consistent():
    data = json.loads(EVIDENCE.read_text(encoding="utf-8"))
    mesh = data["mesh"]
    expected_checks = {
        "torus_genus_one": mesh["chi"] == 0 and mesh["genus"] == 1,
        "first_betti_number_two": data["b1"] == 2,
        "qr_basis_is_class_mixed": data["qr_basis_is_class_mixed"],
        "pure_toroidal_candidate_found": data["n_pure_toroidal_candidates"] > 0,
        "selected_cut_is_toroidal": data["selected_cut_winding"] == [1, 0],
    }

    assert data["schema"] == "radia.validation.cohomology.loop-dof-cut-selection.v1"
    assert data["checks"] == expected_checks
    assert all(expected_checks.values())
    assert data["n_fundamental_cycles"] > data["b1"]
    assert data["selected_cut_length_mm"] > 0.0
    assert data["selected_cut_n_vertices"] >= 3
