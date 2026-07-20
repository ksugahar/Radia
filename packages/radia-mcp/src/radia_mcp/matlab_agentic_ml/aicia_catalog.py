"""Public-metadata catalog gate for AIcia Solid CAE-AI lesson promotion."""

from __future__ import annotations

from collections.abc import Mapping

CHANNEL_ID = "UC2lJYodMaAfFeFQrGUwhlaQ"
REQUIRED_PROMOTION_GATES = {
    "seed_or_determinism",
    "units_and_schema",
    "provenance",
    "independent_forward_solver_verification",
}


def validate_aicia_catalog(catalog: Mapping[str, object]) -> dict[str, object]:
    """Reject incomplete channel snapshots and unsafe CAE lesson promotion.

    Titles and public channel metadata may nominate a lesson, but never establish
    a numerical result. Every promoted candidate must retain the independent
    forward-solver gate.
    """
    counts = catalog.get("counts") if isinstance(catalog.get("counts"), Mapping) else {}
    policy = catalog.get("policy") if isinstance(catalog.get("policy"), Mapping) else {}
    items = catalog.get("items") if isinstance(catalog.get("items"), list) else []
    pending = [item for item in items if isinstance(item, Mapping) and not item.get("previously_verified")]
    ids = [item.get("id") for item in items if isinstance(item, Mapping)]
    dispositions = {"candidate", "review", "not_promoted"}
    candidates = [item for item in pending if item.get("disposition") == "candidate"]
    checks = {
        "channel_is_bound": catalog.get("channel_id") == CHANNEL_ID,
        "full_channel_scope_is_397": counts.get("videos") == 245 and counts.get("shorts") == 2 and counts.get("streams") == 150 and counts.get("total") == 397,
        "ids_are_complete_and_unique": len(items) == 397 and len(ids) == len(set(ids)) == counts.get("unique_ids"),
        "prior_and_new_counts_close": counts.get("previously_verified") == 244 and counts.get("processed_now") == len(pending) == 153,
        "every_new_item_is_dispositioned": all(item.get("disposition") in dispositions for item in pending),
        "disposition_counts_close": sum(counts.get(name, -1) for name in dispositions) == 153,
        "metadata_only_policy_is_enforced": policy.get("metadata_only") is True and policy.get("transcripts_downloaded") is False and policy.get("media_downloaded") is False,
        "generated_candidates_are_not_ground_truth": policy.get("generated_candidate_is_ground_truth") is False and policy.get("promotion_requires_forward_solver_verification") is True,
        "candidate_promotions_keep_solver_gates": all(REQUIRED_PROMOTION_GATES <= set(item.get("promotion_requirements", [])) for item in candidates),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "aicia_public_metadata_cae_promotion_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "processed_now": len(pending),
        "candidate_count": len(candidates),
    }
