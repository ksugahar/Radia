"""AGE guidance and import boundary in the Radia motor package."""

import subprocess
import sys


def test_age_topics_are_served_by_motor_guidance():
    from radia_mcp.motor.age import AGE_TOPICS
    from radia_mcp.motor.age_quality_knowledge import get_age_quality_report

    for topic, text in AGE_TOPICS.items():
        assert get_age_quality_report(topic) == text


def test_age_package_does_not_import_radia():
    import os
    import radia_mcp

    src = os.path.dirname(os.path.dirname(radia_mcp.__file__))
    code = ("import sys; sys.path.insert(0, sys.argv[1]); import radia_mcp.motor.age.lanes; "
            "sys.exit(1 if any(m == 'radia' or m.startswith('radia.') for m in sys.modules) else 0)")
    assert subprocess.run([sys.executable, "-I", "-c", code, src]).returncode == 0
