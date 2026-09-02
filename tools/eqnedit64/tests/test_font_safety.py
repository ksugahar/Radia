"""Static guardrails for the Windows session-safe font path."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    render = (ROOT / "src/equation_render.cpp").read_text(encoding="utf-8")
    workflow = (
        ROOT.parents[1] / ".github/workflows/eqnedit64.yml"
    ).read_text(encoding="utf-8")

    failures: list[str] = []
    if "AddFontMemResourceEx(" in render:
        failures.append("unsafe memory-font registration returned")
    for required in (
        "AddFontResourceExW(path.c_str(), FR_PRIVATE | FR_NOT_ENUM",
        "cache_embedded_math_font",
        "file_matches_bytes(target, bytes, size)",
    ):
        if required not in render:
            failures.append(f"missing file-backed font invariant: {required}")
    if "test_font_session.ps1" not in workflow:
        failures.append("Eqnedit64 CI does not run the font-session stress gate")
    font_session = (ROOT / "build/test_font_session.ps1").read_text(
        encoding="utf-8"
    )
    ui_fuzz = (ROOT / "build/test_ui_fuzz.ps1").read_text(encoding="utf-8")
    for name, script in (("font lifecycle", font_session),
                         ("UI fuzz", ui_fuzz)):
        for required in ("INCONCLUSIVE", "pre-test control", "FromMinutes(10)"):
            if required not in script:
                failures.append(
                    f"{name} guard cannot distinguish ambient font-host churn: "
                    f"missing {required!r}"
                )

    if failures:
        for failure in failures:
            print(f"FAIL  {failure}")
        return 1
    print("PASS: file-backed private font and CI session-health contract")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
