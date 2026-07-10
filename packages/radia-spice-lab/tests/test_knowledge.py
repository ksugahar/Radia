def test_new_public_import_alias_exports_conversion_and_knowledge():
    import spice_circuit_lab as scl

    assert callable(scl.netlist_to_asc)
    assert callable(scl.buck_seed)
    assert callable(scl.circuit_knowledge)
    assert callable(scl.patentability_search_plan)


def test_buck_seed_generates_clean_switching_netlist():
    from ltspice_converter.cli import check_text
    from ltspice_converter.knowledge import buck_seed

    seed = buck_seed(24, 5, 1, fsw_hz=100_000)
    netlist = seed.to_netlist()
    info, warn = check_text(netlist, "cir")

    assert seed.to_dict()["duty"] > 0.2
    assert any("9 components" in item for item in info), info
    assert warn == []


def test_circuit_knowledge_returns_public_rules():
    from ltspice_converter import circuit_knowledge

    out = circuit_knowledge("buck converter")
    assert "rules" in out
    assert any("feedback" in rule.lower() for rule in out["rules"])


def test_patentability_search_plan_expands_circuit_queries():
    from ltspice_converter import patentability_search_plan

    plan = patentability_search_plan(
        title="snubber-assisted boost converter",
        features=["boost converter", "switch-node RC snubber", "soft start"],
        effects=["reduced ringing", "lower overshoot"],
        domains=["power electronics", "circuit"],
    )

    assert "google_scholar" in plan
    assert "google_patents" in plan
    assert plan["query_sets"]["google_scholar"] == plan["google_scholar"]
    assert plan["query_sets"]["web"] == plan["web"]
    assert "jplatpat_keywords_ja" in plan
    assert any("boost converter" in query for query in plan["google_scholar"])
    assert "昇圧" in plan["jplatpat_keywords_ja"]
    assert "スナバ" in plan["jplatpat_keywords_ja"]
    assert "not a legal opinion" in plan["disclaimer"]
