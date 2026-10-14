"""User-Agent for outgoing metadata requests (Crossref, arXiv, publishers)."""

from __future__ import annotations

import os


def polite_user_agent(product: str = "radia-mcp", suffix: str = "") -> str:
    """Return a User-Agent naming the operator's contact, never the maintainer's.

    Crossref's polite pool asks for a ``mailto:``; it is taken from
    ``RADIA_MCP_CONTACT_EMAIL`` and omitted when that variable is unset.
    """
    contact = os.environ.get("RADIA_MCP_CONTACT_EMAIL", "").strip()
    detail = f"mailto:{contact}" if contact else "+https://pypi.org/project/radia-mcp"
    return f"{product} ({detail})" + (f" {suffix}" if suffix else "")
