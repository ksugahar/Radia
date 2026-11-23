"""Geometric and label-based comparison of STEP readbacks."""


def _boolean_result(operation):
    from build123d import Shape

    if not operation.IsDone():
        raise ValueError("STEP comparison Boolean operation failed")
    result = operation.Shape()
    return None if result.IsNull() else Shape.cast(result)


def _union(shape):
    first, *others = shape.solids()
    if not others:
        return first
    result = first.fuse(*others)
    if result is None or not result.is_valid or not result.solids():
        raise ValueError("STEP solid union did not produce valid solids")
    return result


def _volume(shape):
    return 0.0 if shape is None else sum(solid.volume for solid in shape.solids())


def _measure(shape, occupied):
    bounds = shape.bounding_box()
    return {
        "volume_mm3": _volume(occupied),
        "bbox_min_mm": list(bounds.min),
        "bbox_max_mm": list(bounds.max),
        "topology": {
            "solids": len(shape.solids()),
            "faces": len(shape.faces()),
            "edges": len(shape.edges()),
            "vertices": len(shape.vertices()),
        },
    }


def _compare(before, after, occupied_before, occupied_after):
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut

    old, new = _measure(before, occupied_before), _measure(after, occupied_after)
    delta = {
        "volume_mm3": new["volume_mm3"] - old["volume_mm3"],
        "topology": {key: new["topology"][key] - count
                     for key, count in old["topology"].items()},
    }
    for key in ("bbox_min_mm", "bbox_max_mm"):
        delta[key] = [b - a for a, b in zip(old[key], new[key], strict=True)]
    symmetric_difference = sum(
        _volume(_boolean_result(BRepAlgoAPI_Cut(a.wrapped, b.wrapped)))
        for a, b in ((occupied_before, occupied_after), (occupied_after, occupied_before))
    )
    return {"before": old, "after": new, "delta": delta,
            "symmetric_difference_volume_mm3": symmetric_difference}


def _by_label(pairs):
    by_label = {}
    for solid_id, (label, solid) in enumerate(pairs, 1):
        by_label.setdefault(label.casefold(), []).append((solid_id, label, solid))
    return by_label


def compare_steps(before, after, before_pairs, after_pairs):
    old, new = _by_label(before_pairs), _by_label(after_pairs)
    parts, unresolved = [], []
    for key in sorted(old.keys() | new.keys()):
        left, right = old.get(key, []), new.get(key, [])
        if key and len(left) == len(right) == 1:
            (old_id, old_label, a), (new_id, new_label, b) = left[0], right[0]
            parts.append({
                "before": {"id": old_id, "label": old_label},
                "after": {"id": new_id, "label": new_label},
                "comparison": _compare(a, b, a, b),
            })
        else:
            issues = []
            if not key:
                issues.append("missing_label")
            else:
                if not left or not right:
                    issues.append("missing_match")
                if len(left) > 1 or len(right) > 1:
                    issues.append("duplicate_label")
            unresolved.append({
                "label_key": key,
                "issues": issues,
                "before": [{"id": i, "label": label} for i, label, _ in left],
                "after": [{"id": i, "label": label} for i, label, _ in right],
            })
    return {
        "assembly": _compare(before, after, _union(before), _union(after)),
        "matched_parts": parts,
        "unresolved": unresolved,
    }
