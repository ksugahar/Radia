"""Fast guards for the current kernel-specific HACApK architecture."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "src" / "core"


def test_retired_compact_magnetostatic_manager_stays_removed():
    sources = (
        CORE / "rad_hacapk.h",
        CORE / "rad_hacapk.cpp",
        CORE / "rad_interaction.h",
        CORE / "rad_interaction.cpp",
    )
    for path in sources:
        assert "RadHACApKMagnetostaticManager" not in path.read_text(
            encoding="utf-8"
        )


def test_hacapk_base_is_owned_by_current_kernel_managers():
    base_header = (CORE / "rad_hacapk.h").read_text(encoding="utf-8")
    assert "class RadHACApKBase" in base_header

    expected_managers = {
        "rad_hacapk_hdiv.h": "class RadHACApKChargeGram : public RadHACApKBase",
        "rad_hacapk_peec.h": "class RadHACApKPEECManager : public RadHACApKBase",
        "rad_hacapk_bem.h": "class RadHACApKBEMManager : public RadHACApKBase",
    }
    for filename, declaration in expected_managers.items():
        assert declaration in (CORE / filename).read_text(encoding="utf-8")


def test_callback_state_contains_no_retired_kernel_payload():
    callback_files = (
        (CORE / "rad_hacapk.h").read_text(encoding="utf-8"),
        (CORE / "rad_hacapk.cpp").read_text(encoding="utf-8"),
    )
    retired_names = (
        "SetInteraction",
        "SetInvChi",
        "GetInvChi",
        "SetLod",
        "ClearLod",
        "g_hacapk_generation",
    )
    for source in callback_files:
        for name in retired_names:
            assert name not in source
