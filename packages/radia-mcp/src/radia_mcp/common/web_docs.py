"""Live web documentation fetcher with on-disk cache.

Targets:
  - Coreform Cubit help pages (coreform.com/cubit_help)
  - Coreform Cubit forum      (forum.coreform.com)
  - build123d documentation   (build123d.readthedocs.io)

Design:
  - stdlib only (urllib + html.parser) — no httpx/requests dep
  - Cached under `<state_dir>/cache/web/<sha1>.json` with TTL
  - Returns plain-text extract (HTML tags stripped, scripts/styles removed)

Used by `cubit_web_docs` / `build123d_web_docs` MCP tools. Fail-silent:
network errors return a structured error dict rather than raising.
"""

from __future__ import annotations

import hashlib
import html.parser
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

from .failure_log import state_dir


_DEFAULT_TTL_SECONDS = 24 * 3600  # 1 day
_USER_AGENT = "radia-mcp/0.16 (+https://pypi.org/project/radia-mcp)"
_TIMEOUT_SECONDS = 15


def github_token() -> str | None:
	"""Discover a GitHub PAT for higher rate limits / GraphQL access.

	Priority:
	  1. `GITHUB_TOKEN` env var (CI standard)
	  2. `GH_TOKEN` env var (gh CLI standard)
	  3. `gh auth token` (gh CLI — uses OS keyring)

	Returns the token string, or None if no source is available.
	Anonymous callers still get the 60 req/h path; authenticated
	callers get 5000 req/h and can reach GraphQL endpoints (Discussions).
	"""
	import os
	import subprocess
	for key in ("GITHUB_TOKEN", "GH_TOKEN"):
		tok = os.environ.get(key)
		if tok and tok.strip():
			return tok.strip()
	try:
		r = subprocess.run(
			["gh", "auth", "token"],
			capture_output=True, text=True, timeout=5, check=False,
		)
	except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
		return None
	if r.returncode == 0:
		tok = r.stdout.strip()
		if tok:
			return tok
	return None


def _github_auth_headers() -> dict[str, str]:
	tok = github_token()
	if tok:
		return {"Authorization": f"Bearer {tok}"}
	return {}


def _cache_dir() -> Path:
	d = state_dir() / "cache" / "web"
	d.mkdir(parents=True, exist_ok=True)
	return d


def _cache_path(url: str) -> Path:
	h = hashlib.sha1(url.encode("utf-8")).hexdigest()[:16]
	return _cache_dir() / f"{h}.json"


def _cache_get(url: str, ttl: int) -> Optional[dict]:
	p = _cache_path(url)
	if not p.exists():
		return None
	try:
		d = json.loads(p.read_text(encoding="utf-8"))
	except (OSError, json.JSONDecodeError):
		return None
	if time.time() - d.get("fetched_at", 0) > ttl:
		return None
	return d


def _cache_put(url: str, data: dict) -> None:
	p = _cache_path(url)
	try:
		p.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
	except OSError:
		pass


class _TextExtractor(html.parser.HTMLParser):
	"""Strip tags; skip <script>/<style>/<nav>/<header>/<footer>."""

	_SKIP_TAGS = {"script", "style", "nav", "header", "footer", "aside", "noscript"}

	def __init__(self) -> None:
		super().__init__(convert_charrefs=True)
		self._buf: list[str] = []
		self._skip_depth = 0
		self._current_tag = ""

	def handle_starttag(self, tag: str, attrs):  # noqa: ANN001
		tag = tag.lower()
		self._current_tag = tag
		if tag in self._SKIP_TAGS:
			self._skip_depth += 1
		# Add line breaks around block-level tags for readability
		if tag in ("p", "br", "li", "h1", "h2", "h3", "h4", "h5", "pre", "div", "tr"):
			self._buf.append("\n")
		if tag in ("h1", "h2", "h3", "h4"):
			self._buf.append("## ")

	def handle_endtag(self, tag: str):
		tag = tag.lower()
		if tag in self._SKIP_TAGS and self._skip_depth > 0:
			self._skip_depth -= 1
		if tag in ("p", "li", "h1", "h2", "h3", "h4", "h5", "pre", "div"):
			self._buf.append("\n")

	def handle_data(self, data: str):
		if self._skip_depth > 0:
			return
		self._buf.append(data)

	def text(self) -> str:
		raw = "".join(self._buf)
		raw = re.sub(r"[ \t]+", " ", raw)
		raw = re.sub(r"\n[ \t]+", "\n", raw)
		raw = re.sub(r"\n{3,}", "\n\n", raw)
		return raw.strip()


def fetch(url: str, ttl: int = _DEFAULT_TTL_SECONDS, force_refresh: bool = False) -> dict:
	"""Fetch URL with caching. Returns dict with status/text/url/fetched_at."""
	if not force_refresh:
		c = _cache_get(url, ttl)
		if c is not None:
			c["from_cache"] = True
			return c

	req = urllib.request.Request(url, headers={
		"User-Agent": _USER_AGENT,
		"Accept": "text/html,application/xhtml+xml",
	})
	try:
		with urllib.request.urlopen(req, timeout=_TIMEOUT_SECONDS) as resp:  # noqa: S310
			raw = resp.read()
			charset = resp.headers.get_content_charset() or "utf-8"
			final_url = resp.geturl()
	except urllib.error.HTTPError as e:
		return {"status": "error", "url": url, "error": f"HTTP {e.code}: {e.reason}"}
	except urllib.error.URLError as e:
		return {"status": "error", "url": url, "error": f"URL error: {e.reason}"}
	except (TimeoutError, OSError) as e:
		return {"status": "error", "url": url, "error": f"Network error: {e}"}

	try:
		html_text = raw.decode(charset, errors="replace")
	except LookupError:
		html_text = raw.decode("utf-8", errors="replace")

	parser = _TextExtractor()
	try:
		parser.feed(html_text)
	except Exception as e:  # noqa: BLE001
		return {"status": "error", "url": url, "error": f"Parse error: {e}"}

	text = parser.text()
	# Trim extremely long pages (rate-limit retrieval cost)
	if len(text) > 60_000:
		text = text[:60_000] + "\n\n… [truncated at 60k chars]"

	data = {
		"status": "ok",
		"url": final_url,
		"requested_url": url,
		"text": text,
		"fetched_at": time.time(),
		"from_cache": False,
	}
	_cache_put(url, data)
	return data


def _fetch_json(url: str, ttl: int = _DEFAULT_TTL_SECONDS,
                force_refresh: bool = False) -> dict:
	"""Fetch a JSON endpoint with the same on-disk cache as `fetch`."""
	if not force_refresh:
		c = _cache_get(url, ttl)
		if c is not None:
			c["from_cache"] = True
			return c
	req = urllib.request.Request(url, headers={
		"User-Agent": _USER_AGENT,
		"Accept": "application/json",
	})
	try:
		with urllib.request.urlopen(req, timeout=_TIMEOUT_SECONDS) as resp:  # noqa: S310
			raw = resp.read()
	except urllib.error.HTTPError as e:
		return {"status": "error", "url": url, "error": f"HTTP {e.code}: {e.reason}"}
	except urllib.error.URLError as e:
		return {"status": "error", "url": url, "error": f"URL error: {e.reason}"}
	except (TimeoutError, OSError) as e:
		return {"status": "error", "url": url, "error": f"Network error: {e}"}
	try:
		data = json.loads(raw.decode("utf-8", errors="replace"))
	except json.JSONDecodeError as e:
		return {"status": "error", "url": url, "error": f"JSON decode: {e}"}
	out = {
		"status": "ok",
		"url": url,
		"json": data,
		"fetched_at": time.time(),
		"from_cache": False,
	}
	_cache_put(url, out)
	return out


def search_forum(query: str, base: str = "https://forum.coreform.com",
                 max_hits: int = 8, force_refresh: bool = False) -> list[dict]:
	"""Search a Discourse forum via its JSON search endpoint.

	Discourse exposes `/search.json?q=...` with a well-typed structure
	(topics / posts / excerpts). This is far more reliable than HTML
	scraping a JS-rendered help site, and gives us authoritative, live
	community answers for Cubit-specific questions.
	"""
	q = (query or "").strip()
	if not q:
		return []
	url = f"{base.rstrip('/')}/search.json?q={urllib.parse.quote(q)}"
	r = _fetch_json(url, force_refresh=force_refresh)
	if r.get("status") != "ok":
		return [{"error": r.get("error"), "url": url}]
	data = r.get("json") or {}
	topics_by_id = {t["id"]: t for t in data.get("topics", []) if isinstance(t, dict)}
	hits: list[dict] = []
	for post in data.get("posts", [])[: max(1, int(max_hits))]:
		if not isinstance(post, dict):
			continue
		topic_id = post.get("topic_id")
		topic = topics_by_id.get(topic_id, {})
		slug = topic.get("slug", "")
		tid = topic.get("id", topic_id)
		pnum = post.get("post_number", 1)
		link = f"{base.rstrip('/')}/t/{slug}/{tid}/{pnum}" if slug and tid else f"{base.rstrip('/')}/p/{post.get('id','')}"
		hits.append({
			"title": topic.get("title", ""),
			"excerpt": (post.get("blurb") or "")[:400],
			"username": post.get("username"),
			"created_at": post.get("created_at"),
			"url": link,
			"from_cache": r.get("from_cache", False),
		})
	return hits


def search_text(text: str, query: str, context_lines: int = 4, max_hits: int = 8) -> list[dict]:
	"""Grep-like search inside fetched text. Case-insensitive term AND."""
	terms = [t.lower() for t in query.split() if t.strip()]
	if not terms:
		return []
	lines = text.splitlines()
	hits: list[dict] = []
	for i, line in enumerate(lines):
		ll = line.lower()
		if all(t in ll for t in terms):
			start = max(0, i - context_lines)
			end = min(len(lines), i + context_lines + 1)
			hits.append({
				"line": i + 1,
				"match": line.strip(),
				"context": "\n".join(lines[start:end]),
			})
			if len(hits) >= max_hits:
				break
	return hits


# Known curated roots. Tools pass `query` + pick a source; we hit the
# search/index page and list top candidate URLs. Most sites expose
# search-friendly URL patterns.
DOCS_INDEX: dict[str, list[str]] = {
	"cubit": [
		"https://coreform.com/cubit_help/cubithelp.htm",
	],
	"cubit_forum": [
		"https://forum.coreform.com/",
	],
	"build123d": [
		"https://build123d.readthedocs.io/en/latest/",
		"https://build123d.readthedocs.io/en/latest/key_concepts.html",
		"https://build123d.readthedocs.io/en/latest/objects.html",
		"https://build123d.readthedocs.io/en/latest/operations.html",
		"https://build123d.readthedocs.io/en/latest/builders.html",
		"https://build123d.readthedocs.io/en/latest/tutorial_selectors.html",
		"https://build123d.readthedocs.io/en/latest/import_export.html",
		"https://build123d.readthedocs.io/en/latest/direct_api_reference.html",
		"https://build123d.readthedocs.io/en/latest/introductory_examples.html",
	],
}
