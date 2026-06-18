"""Coverage for the synthesized Coreform Cubit YouTube knowledge corpus.

Locks in the topic set, alias resolution, and the index/all dispatch so the
caption-derived corpus (Tutorials + Getting-started + Flex/IGA + Clips) stays
queryable via cubit_docs("coreform_<topic>").
"""
import pytest

from radia_mcp.cubit.knowledge.coreform_webinars import (
    get_coreform_webinar_documentation as get_doc,
    _TOPICS,
    _ALIASES,
)

EXPECTED_TOPICS = {
    "python_automation", "meshing_strategy", "solver_workflows",
    "getting_started", "advanced_meshing", "ml_and_gui",
    "neutronics_fusion", "domain_applications",
    "getting_started_2025", "flex_iga", "clips_highlights",
}


def test_topic_set():
    assert set(_TOPICS) == EXPECTED_TOPICS


@pytest.mark.parametrize("topic", sorted(EXPECTED_TOPICS))
def test_each_topic_is_substantive(topic):
    doc = get_doc(topic)
    assert doc.lstrip().startswith("##"), topic
    assert "### Sources" in doc, topic
    assert len(doc) > 1500, topic
    assert '"""' not in doc, topic  # must stay r-string-safe


def test_index_lists_every_topic():
    idx = get_doc("index")
    for topic in EXPECTED_TOPICS:
        assert topic in idx, topic


def test_all_concatenates_everything():
    all_doc = get_doc("all")
    # at least the sum of the individual topic lengths' worth of content
    assert len(all_doc) > sum(len(_TOPICS[k]) for k in _TOPICS) * 0.9
    assert "U-spline" in all_doc       # from flex_iga
    assert "Enceladus" in all_doc      # from clips_highlights
    assert "Power Tools" in all_doc    # from getting_started_2025


@pytest.mark.parametrize("alias,expected", [
    ("iga", "flex_iga"),
    ("u_splines", "flex_iga"),
    ("flex", "flex_iga"),
    ("bezier", "flex_iga"),
    ("nurbs", "flex_iga"),
    ("power_tools", "getting_started_2025"),
    ("journaling", "getting_started_2025"),
    ("associate_edition", "getting_started_2025"),
    ("trelis", "getting_started_2025"),
    ("geodynamics", "clips_highlights"),
    ("enceladus", "clips_highlights"),
    ("clips", "clips_highlights"),
    # pre-existing aliases still resolve
    ("hex", "meshing_strategy"),
    ("dagmc", "neutronics_fusion"),
    ("moose", "solver_workflows"),
])
def test_alias_resolution(alias, expected):
    assert _ALIASES[alias] == expected
    assert get_doc(alias) == _TOPICS[expected]


def test_unknown_topic_is_graceful():
    out = get_doc("does_not_exist_xyz")
    assert "Unknown coreform topic" in out


def test_new_topic_content_markers():
    flex = get_doc("flex_iga")
    assert "Bezier extraction" in flex
    assert "U-spline" in flex
    gs = get_doc("getting_started_2025")
    assert "Journal" in gs
    assert "Associate" in gs or "Cubit Learn" in gs
    clips = get_doc("clips_highlights")
    assert "geodynamic" in clips.lower()
