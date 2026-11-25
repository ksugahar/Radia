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


def _distance(a, b):
    return sum((u - v) ** 2 for u, v in zip(a.center(), b.center(), strict=True)) ** 0.5


def _match_by_placement(left, right, tolerance_mm=1e-6):
    """Pair same-label instances by centroid distance, or return None when ambiguous.

    A pair is accepted only when its distance is below half of the next-best
    alternative for both instances, so an equal shift of a pattern stays unresolved.
    """
    rows = [[_distance(a, b) for _, _, b in right] for _, _, a in left]
    order = sorted((d, i, j) for i, row in enumerate(rows) for j, d in enumerate(row))
    used_i, used_j, chosen = set(), set(), []
    for d, i, j in order:
        if i not in used_i and j not in used_j:
            used_i.add(i)
            used_j.add(j)
            chosen.append((i, j, d))
    for i, j, d in chosen:
        others = [rows[i][k] for k in range(len(right)) if k != j]
        others += [rows[k][j] for k in range(len(left)) if k != i]
        runner_up = min(others)
        if not (runner_up > 2 * d and runner_up - d > tolerance_mm):
            return None, rows
    return sorted(chosen), rows


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
                "match": "label",
                "comparison": _compare(a, b, a, b),
            })
            continue
        chosen = rows = None
        if key and len(left) == len(right) > 1:
            chosen, rows = _match_by_placement(left, right)
        if chosen:
            for i, j, d in chosen:
                (old_id, old_label, a), (new_id, new_label, b) = left[i], right[j]
                parts.append({
                    "before": {"id": old_id, "label": old_label},
                    "after": {"id": new_id, "label": new_label},
                    "match": "placement",
                    "centroid_distance_mm": d,
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
                if rows is not None:
                    issues.append("ambiguous_placement")
            entry = {
                "label_key": key,
                "issues": issues,
                "before": [{"id": i, "label": label} for i, label, _ in left],
                "after": [{"id": i, "label": label} for i, label, _ in right],
            }
            if rows is not None:
                entry["centroid_distance_mm"] = [
                    {"before_id": left[i][0], "after_id": right[j][0], "distance_mm": d}
                    for i, row in enumerate(rows) for j, d in enumerate(row)
                ]
            unresolved.append(entry)
    return {
        "assembly": _compare(before, after, _union(before), _union(after)),
        "matched_parts": parts,
        "unresolved": unresolved,
    }
