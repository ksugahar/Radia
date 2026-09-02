from radia_mcp.radia_ngsolve.knowledge.install_deploy import (
    INSTALL_DEPLOY,
    get_install_deploy_documentation,
)


def test_install_deploy_topics_match_current_machine_roles():
    overview = get_install_deploy_documentation("overview")
    ci_compute = get_install_deploy_documentation("ci_compute")

    assert "LAB is not a CI runner" in overview
    assert "mdx" in overview and "CI priority" in overview
    assert "impact-selected" in ci_compute
    assert "validation_test/" in ci_compute


def test_install_deploy_does_not_publish_retired_binary_copy_routes():
    assert "push_pyds_to_mdx" not in INSTALL_DEPLOY
    assert "download_binaries" not in INSTALL_DEPLOY
    assert "manually dropping" in INSTALL_DEPLOY
    assert "cubit-plugin-install" in get_install_deploy_documentation("cubit")
