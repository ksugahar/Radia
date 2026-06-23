"""Load knowledge as external .md files via importlib.resources.

Learned from wjc9011/COMSOL_Multiphysics_MCP — they keep
`src/knowledge/prompts/{mph_api,physics_guide,workflow}.md` as
separate markdown files and load them as plain text. This avoids
the `r\"\"\"...\"\"\"` raw-string parser collision we hit twice
(see memory/feedback_*.md if it exists, or just remember: any
Python triple-quoted docstring inside a raw-string r\"\"\"..\"\"\"
breaks because `\"\"\"` cannot be escaped inside raw strings).

Convention adopted for radia_mcp:
  packages/radia-mcp/src/radia_mcp/<subpackage>/prompts/<name>.md

Usage:
    from radia_mcp.common import load_prompt
    text = load_prompt(\"ih\", \"karl_overview\")   # → str

Migration order recommended (highest pain first):
  1. ih  — historical r-string bug example
  2. ih.sibc — also hit
  3. peec.hoibc — also hit
  4. bayesian_opt / evolutionary — preventive

Migration mechanics:
  - Extract one constant `FOO_BAR = r\"\"\"...\"\"\"` into
    `<subpackage>/prompts/foo_bar.md`
  - In knowledge.py: `FOO_BAR = load_prompt(\"<subpackage>\", \"foo_bar\")`
  - Keep dispatcher unchanged.
"""

from __future__ import annotations
from importlib import resources
from typing import List


def load_prompt(subpackage: str, name: str) -> str:
    """Load a markdown prompt from radia_mcp.<subpackage>.prompts.

    Args:
        subpackage: subpackage name (e.g. "ih", "peec.hoibc"
                    — dots OK for nested).
        name: prompt file basename WITHOUT the .md extension
              (e.g. "hayaho_pm_motor" → loads hayaho_pm_motor.md).

    Returns:
        The full markdown text as a string.

    Raises:
        FileNotFoundError: if the prompt file does not exist (with a
                           clear list of what IS available in that
                           subpackage, to help debug).
    """
    package = f"radia_mcp.{subpackage}.prompts"
    filename = f"{name}.md"
    try:
        return resources.files(package).joinpath(filename).read_text(
            encoding="utf-8")
    except FileNotFoundError:
        # Helpful error: list what IS available
        try:
            available = list_prompts(subpackage)
            msg = (f"Prompt '{name}.md' not found in {package}. "
                   f"Available: {', '.join(available) if available else '(none)'}")
        except Exception:
            msg = f"Prompt '{name}.md' not found in {package}."
        raise FileNotFoundError(msg)


def list_prompts(subpackage: str) -> List[str]:
    """List all .md prompts in radia_mcp.<subpackage>.prompts.

    Returns:
        List of prompt basenames (without .md extension), sorted.

    Raises:
        ModuleNotFoundError: if the subpackage has no prompts/ directory.
    """
    package = f"radia_mcp.{subpackage}.prompts"
    names = []
    for entry in resources.files(package).iterdir():
        if entry.is_file() and entry.name.endswith(".md"):
            names.append(entry.name[:-3])
    return sorted(names)
