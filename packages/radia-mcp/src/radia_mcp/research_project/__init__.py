"""research_project sub-module — cross-cutting project-level checks.

While each sub-skill (grant_writing / paper_writing / presentation /
poster / bibliography / pdf / doc_convert) covers its own artifact, a
research project usually has 5-10 of these artifacts in one directory.
``research_project`` orchestrates them.

Tools (prefix ``research_project_``):

    Tier 1 — Discovery:
        research_project_scan(directory, recursive=True)
        research_project_health_dashboard(directory, recursive=False)

    Tier 2 — Cross-artifact:
        research_project_deadline_gantt(plan_path, reference_date=None)
        research_project_consistency_check(directory)
"""
from . import tools as _tools


def register(mcp) -> int:
    count = 0
    for name in dir(_tools):
        if name.startswith("research_project_") and callable(getattr(_tools, name)):
            mcp.tool()(getattr(_tools, name))
            count += 1
    return count
