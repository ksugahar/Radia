"""Post-processing examples on existing MSH 4.2 results."""

_EXAMPLES = {
    "display": """
# Review an existing result
Use gmsh_usage("workflow") for the complete sequence and runtime contract.
Inspect and validate field.msh, then create its launch companions with
gmsh_write_post_launch_artifact(msh_path=..., output_base=...). Open field.geo
in Gmsh, or render it with gmsh_render(path=..., view="B", time_step=0).
Keep axes equal, check SI coordinate units, and use a shared color range when
comparing cases. STEP overlays are display-only geometry.
""",
    "fields": """
# Sample and compare result fields
Select the view by its inspected name. gmsh_probe gives values and found status;
gmsh_line_profile samples a declared segment; gmsh_export_csv writes the
existing node/element values. gmsh_compare_fields compares two result fields
with explicit interpolation/domain settings. Samples outside the mesh are
missing data, never a zero field. Keep physical units and vector components.
""",
    "animation": """
# Export existing time-dependent views
Inspect the actual recorded step times. Use gmsh_export_animation with
link_views=True, a fixed camera and a shared range; playback delay is not the
physical timestep. For moving geometry, use the displacement view with its
declared scale. gmsh_usage("magnetic_movies") explains field-line seeding and
time-resolution acceptance. Retain data/companion hashes and runtime identity.
""",
}


def get_gmsh_examples(topic: str = "all") -> str:
    """Return display, field extraction, or animation examples."""
    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(_EXAMPLES.values())
    if topic in _EXAMPLES:
        return _EXAMPLES[topic]
    return f"Unknown topic: '{topic}'. Available: all, {', '.join(_EXAMPLES)}"
