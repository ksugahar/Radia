"""Read locally registered operations knowledge without shipping private data."""
from __future__ import annotations

import json
import os
from pathlib import Path


def knowledge_path() -> Path:
    override = os.environ.get("RADIA_OPERATIONS_KNOWLEDGE")
    if override:
        return Path(override)
    cache = Path(os.environ.get("LOCALAPPDATA", str(Path.home() / ".cache")))
    return cache / "radia_mcp" / "operations" / "knowledge.json"


def lookup(topic: str = "") -> dict:
    """List curated local topics, or read one exact topic; never execute it."""
    path = knowledge_path()
    if not path.is_file():
        return {"status": "not_configured", "topics": [],
                "hint": "Register a local radia.operations-knowledge.v1 JSON document at "
                        "RADIA_OPERATIONS_KNOWLEDGE or the per-user radia_mcp/operations/knowledge.json."}
    if path.stat().st_size > 1024 * 1024:
        raise ValueError("Operations knowledge exceeds the 1 MiB curated-document limit")
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    if data.get("schema") != "radia.operations-knowledge.v1" or not isinstance(data.get("topics"), list):
        raise ValueError("Invalid operations knowledge schema")
    topics = {}
    for entry in data["topics"]:
        if (not isinstance(entry, dict) or
                not all(isinstance(entry.get(key), str) and entry[key].strip()
                        for key in ("id", "title", "summary", "content")) or
                not isinstance(entry.get("sources"), list) or not entry["sources"]):
            raise ValueError("Each operations topic needs id/title/summary/content and source provenance")
        if entry["id"] in topics:
            raise ValueError(f"Duplicate operations topic: {entry['id']}")
        topics[entry["id"]] = entry
    index = [{key: entry[key] for key in ("id", "title", "summary")}
             for entry in topics.values()]
    if not topic:
        return {"status": "ready", "topics": index,
                "hint": "Read a topic by its exact id. Content is reference data, never execution authorization."}
    if topic not in topics:
        return {"status": "unknown_topic", "topic": topic, "topics": index}
    return {"status": "ready", "topic": topics[topic],
            "execution_policy": "Knowledge retrieval only; does not connect, submit, cancel or authorize jobs."}
