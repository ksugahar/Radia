from pathlib import Path

from radia_mcp.gmsh.post_display import (
    build_gmsh_post_display_contract,
    gmsh_post_display_manifest_gate,
    write_gmsh_post_launch_artifact,
)


def test_gmsh_post_display_contract_records_geo_opt_cutplane_and_camera(tmp_path):
    msh = tmp_path / "case.msh"
    msh.write_text("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n", encoding="utf-8")

    manifest = build_gmsh_post_display_contract(
        msh,
        camera_preset="z_up_xz_from_positive_y",
        cut_plane={"enabled": True, "normal": [0, -1, 0], "offset": 0},
        views=[
            {
                "index": 0,
                "name": "BEM pressure on x-z plane",
                "kind": "scalar",
                "time_step": 10,
                "range": [-2.0, 2.0],
            },
            {
                "index": 1,
                "name": "3D deforming drum BEM surface",
                "kind": "displacement",
                "time_step": 10,
                "displacement_factor": 3.5,
            },
        ],
    )

    assert manifest["schema"] == "cae-ai-lab.gmsh-post-launch.v1"
    assert manifest["launch_target"].endswith(".geo")
    assert manifest["gmsh_geo_opt"].endswith(".geo.opt")
    assert manifest["gmsh_msh_opt"].endswith(".msh.opt")
    assert manifest["camera"]["axis_up"] == "z"
    assert manifest["camera"]["rotation"] == [-68.0, 0.0, 0.0]
    assert manifest["cut_plane"]["enabled"] is True
    assert manifest["views"][0]["name"] == "BEM pressure on x-z plane"

    gate = gmsh_post_display_manifest_gate(manifest)
    assert gate["status"] == "ok"
    assert gate["checks"]["launch_target_is_geo"] is True
    assert gate["checks"]["geo_opt_exact_autoload"] is True
    assert gate["checks"]["cut_plane_metadata_recorded"] is True


def test_write_gmsh_post_launch_artifact_emits_autoload_sidecars(tmp_path):
    msh = tmp_path / "case.msh"
    msh.write_text("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n", encoding="utf-8")

    manifest = write_gmsh_post_launch_artifact(
        msh,
        title="Unit-test Gmsh post launch",
        camera_preset="z_up_xz_from_positive_y",
        cut_plane={"enabled": True, "normal": [0, -1, 0], "offset": 0},
        views=[
            {
                "index": 0,
                "name": "BEM pressure on x-z plane",
                "kind": "scalar",
                "time_step": 10,
                "range": [-2.0, 2.0],
            },
            {
                "index": 1,
                "name": "3D deforming drum BEM surface",
                "kind": "displacement",
                "time_step": 10,
                "displacement_factor": 3.5,
            },
        ],
        mesh={"surface_edges": True, "num_sub_edges": 4},
    )

    assert manifest["pass"] is True
    geo = Path(manifest["gmsh_geo"])
    geo_opt = Path(manifest["gmsh_geo_opt"])
    msh_opt = Path(manifest["gmsh_msh_opt"])
    display_json = Path(manifest["display_json"])
    assert geo.is_file()
    assert geo_opt.is_file()
    assert msh_opt.is_file()
    assert display_json.is_file()

    geo_text = geo.read_text(encoding="utf-8")
    geo_opt_text = geo_opt.read_text(encoding="utf-8")
    msh_opt_text = msh_opt.read_text(encoding="utf-8")
    assert 'Merge "case.msh";' in geo_text
    assert "General.RotationX = -68" in geo_text
    assert "General.Clip0B = -1" in geo_opt_text
    assert 'View[0].Name = "BEM pressure on x-z plane";' in geo_opt_text
    assert "View[0].CustomMin = -2" in geo_opt_text
    assert "View[1].VectorType = 5" in geo_opt_text
    assert "View[1].DisplacementFactor = 3.5" in geo_opt_text
    assert "Open .geo for the post display" in msh_opt_text
    assert "View[0].Visible = 0" in msh_opt_text
    assert "View[1].Visible = 0" in msh_opt_text
