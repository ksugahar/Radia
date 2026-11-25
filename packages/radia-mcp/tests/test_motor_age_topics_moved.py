"""AGE method notes moved from ngsolve-mcp.age into radia-mcp.motor."""

import subprocess
import sys


def test_age_topics_same_text_on_legacy_and_new_route():
    from ngsolve_mcp.age.server import age_usage
    from radia_mcp.motor.age import AGE_TOPICS
    from radia_mcp.motor.age_quality_knowledge import get_age_quality_report

    for topic, text in AGE_TOPICS.items():
        assert get_age_quality_report(topic) == text
        assert age_usage(topic) == text
    assert "motor_age_quality" in age_usage("overview")


def test_age_package_does_not_import_radia():
    import os
    import radia_mcp

    src = os.path.dirname(os.path.dirname(radia_mcp.__file__))
    code = ("import sys; sys.path.insert(0, sys.argv[1]); import radia_mcp.motor.age; "
            "sys.exit(1 if any(m == 'radia' or m.startswith('radia.') for m in sys.modules) else 0)")
    assert subprocess.run([sys.executable, "-I", "-c", code, src]).returncode == 0
