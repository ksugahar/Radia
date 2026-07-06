"""Scraped example library for build123d + Cubit.

Two sources, two layouts:

* **build123d** — GitHub `gumyr/build123d/examples/` (~65 curated Python
  files).  Listing via `api.github.com/repos/.../contents/examples`,
  bodies via `raw.githubusercontent.com`.  One cache file per example.

* **cubit** — Coreform forum (`forum.coreform.com`, Discourse).  Posts
  tagged / matching `tutorial example mesh` are fetched (JSON API),
  markdown code fences (```...```) are extracted as individual snippets.
  Each post becomes one "example" (metadata + concatenated code).

Both layouts serialize to `<state_dir>/examples/{source}/`:
  - `index.json` — list of `{name, title, url, code_path, tokens, token_set}`
  - `<name>.py` / `<name>.jou` — the actual code

The tf-idf retrieval machinery (`search_examples`) is shared with the
lookup tools — same scoring / heading-boost approach.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from .failure_log import state_dir
from .web_docs import _USER_AGENT, _TIMEOUT_SECONDS


_WORD_RE = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*")
_CODE_FENCE_RE = re.compile(r"```[a-z]*\n?([\s\S]*?)```", re.MULTILINE)
_REFRESH_TTL_SECONDS = 7 * 24 * 3600  # one week


def _examples_dir(source: str) -> Path:
	d = state_dir() / "examples" / source
	d.mkdir(parents=True, exist_ok=True)
	return d


def _index_path(source: str) -> Path:
	return _examples_dir(source) / "index.json"


def _tokenize(s: str) -> list[str]:
	"""Lower-case identifier tokens + their underscore-split sub-tokens.

	`roller_coaster` indexes as both `roller_coaster` and `roller` +
	`coaster`, so both the exact identifier and free-form queries match.
	"""
	out: list[str] = []
	for m in _WORD_RE.finditer(s):
		w = m.group(0).lower()
		out.append(w)
		if "_" in w:
			out.extend(p for p in w.split("_") if p)
	return out


def _http_get(url: str, accept: str = "text/plain") -> bytes | None:
	headers = {
		"User-Agent": _USER_AGENT,
		"Accept": accept,
	}
	# Authenticate GitHub API calls when a token is available: 60 → 5000
	# req/h, and allows GraphQL / Discussions access. Safe no-op otherwise.
	if "api.github.com" in url or "raw.githubusercontent.com" in url:
		from . import web_docs as _wd
		headers.update(_wd._github_auth_headers())
	req = urllib.request.Request(url, headers=headers)
	try:
		with urllib.request.urlopen(req, timeout=_TIMEOUT_SECONDS) as resp:  # noqa: S310
			return resp.read()
	except (urllib.error.HTTPError, urllib.error.URLError,
	        TimeoutError, OSError):
		return None


def _safe_name(s: str) -> str:
	"""Sanitize a filename component."""
	return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)[:96] or "unnamed"


# ---------------------------------------------------------------------------
# build123d: GitHub `gumyr/build123d/examples/*.py`
# ---------------------------------------------------------------------------

_BUILD123D_LIST_URL = ("https://api.github.com/repos/gumyr/build123d/"
                       "contents/examples")

_BD_WAREHOUSE_LIST_URL = ("https://api.github.com/repos/gumyr/bd_warehouse/"
                           "contents/src/bd_warehouse")
_BD_WAREHOUSE_EXAMPLES_URL = ("https://api.github.com/repos/gumyr/bd_warehouse/"
                               "contents/examples")

# Default local roots walked by refresh_cubit_local_examples.  The retired
# Radia examples/ tree is intentionally absent; durable in-repo material lives
# in docs/, validation_test/, and panel samples.
_DEFAULT_LOCAL_CUBIT_ROOTS = [
	r"public-safe curated corpus",
	r"repo:/docs",
	r"repo:/validation_test",
	r"repo:/src/radia/panels/samples",
]

# Source families — `search_examples(family, ...)` unions all sub-sources
# so a single query sees forum + local + durable Radia docs/validation assets,
# GitHub examples + bd_warehouse for build123d, etc.
FAMILIES: dict[str, list[str]] = {
	"cubit": ["cubit", "cubit_local",
	          "cubit_youtube", "cubit_jou_github"],
	"build123d": ["build123d", "bd_warehouse",
	              "build123d_discussions",
	              "build123d_github_discussions",
	              "build123d_youtube"],
	"gmsh": ["gmsh_issues", "gmsh_stackoverflow",
	         "gmsh_youtube"],
}


def refresh_build123d_examples(limit: int = 0) -> dict:
	"""Fetch the build123d examples folder from GitHub and populate the cache.

	Args:
	    limit: cap on number of files to fetch (0 = all).

	Returns stats dict: {fetched, failed, total, index_path}.
	"""
	raw = _http_get(_BUILD123D_LIST_URL, accept="application/vnd.github+json")
	if not raw:
		return {"status": "error", "error": "GitHub listing fetch failed"}
	try:
		listing = json.loads(raw)
	except json.JSONDecodeError as e:
		return {"status": "error", "error": f"listing JSON decode: {e}"}
	if not isinstance(listing, list):
		return {"status": "error",
		        "error": f"unexpected listing shape: {str(listing)[:200]}"}

	py_files = [f for f in listing
	            if isinstance(f, dict) and f.get("name", "").endswith(".py")]
	if limit > 0:
		py_files = py_files[:limit]

	d = _examples_dir("build123d")
	index: list[dict[str, Any]] = []
	fetched = 0
	failed: list[str] = []

	for f in py_files:
		name = f.get("name", "")
		dl = f.get("download_url")
		if not (name and dl):
			failed.append(name)
			continue
		body_raw = _http_get(dl, accept="text/plain")
		if body_raw is None:
			failed.append(name)
			continue
		body = body_raw.decode("utf-8", errors="replace")

		code_path = d / _safe_name(name)
		try:
			code_path.write_text(body, encoding="utf-8")
		except OSError:
			failed.append(name)
			continue

		# Extract title from docstring/comments
		title = _guess_python_title(body, fallback=name)
		tokens = _tokenize(body)
		index.append({
			"name": name,
			"title": title,
			"url": f.get("html_url")
			       or f"https://github.com/gumyr/build123d/blob/dev/examples/{name}",
			"code_path": str(code_path),
			"size": f.get("size", 0),
			"tokens": tokens,
			"token_set": sorted(set(tokens)),
		})
		fetched += 1

	index_data = {
		"fetched_at": time.time(),
		"source": "build123d",
		"items": index,
	}
	_index_path("build123d").write_text(
		json.dumps(index_data, ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"fetched": fetched,
		"failed": failed,
		"total": len(py_files),
		"index_path": str(_index_path("build123d")),
	}


def _guess_python_title(body: str, fallback: str) -> str:
	"""Pick a 1-line title from a Python file.

	Priority: module docstring first line → first comment line → fallback.
	"""
	# Module docstring
	m = re.match(r'\s*(?:[ru]?["\']{3})([\s\S]*?)["\']{3}', body)
	if m:
		first = m.group(1).strip().splitlines()
		if first:
			return first[0][:120]
	# First `# ...` comment line
	for line in body.splitlines()[:20]:
		ls = line.strip()
		if ls.startswith("#") and len(ls) > 1:
			return ls.lstrip("# ").strip()[:120]
	return fallback


# ---------------------------------------------------------------------------
# build123d_discussions: GitHub Issues on `gumyr/build123d` (de-facto forum)
# ---------------------------------------------------------------------------

_B3D_ISSUES_URL = "https://api.github.com/repos/gumyr/build123d/issues"
_B3D_COMMENTS_URL = (
	"https://api.github.com/repos/gumyr/build123d/issues/{num}/comments"
)


def refresh_build123d_discussions(max_issues: int = 60,
                                   include_comments: bool = True,
                                   min_discussion_signals: int = 1) -> dict:
	"""Fetch recent issues (with comments) from `gumyr/build123d`.

	GitHub Discussions require GraphQL + PAT; Issues on this repo are
	the pragmatic equivalent — questions, bug reports, and multi-turn
	Q&A all live there, REST-accessible anonymously.

	Walks `GET /repos/gumyr/build123d/issues?state=all` pagewise until
	`max_issues` bodies are indexed. Pull requests are skipped (body
	quality there is code-diff, not discussion). Issues with
	`comments >= min_discussion_signals` also have their comment
	threads folded in.

	Anonymous GitHub API rate is 60 req/h — tune `max_issues` to stay
	under it (default 60 issues × 1 comment fetch each = up to 61 req).

	Args:
	    max_issues: cap on issues indexed (default 60).
	    include_comments: fetch comment threads for issues that have them.
	    min_discussion_signals: minimum number of comments to consider
	        "discussion-like" (default 1).
	"""
	d = _examples_dir("build123d_discussions")
	index: list[dict[str, Any]] = []
	fetched = 0
	failed: list[str] = []

	page = 1
	while fetched < max_issues:
		url = f"{_B3D_ISSUES_URL}?state=all&per_page=100&page={page}&sort=updated"
		raw = _http_get(url, accept="application/vnd.github+json")
		if raw is None:
			failed.append(f"issues page {page}")
			break
		try:
			issues = json.loads(raw)
		except json.JSONDecodeError:
			failed.append(f"issues page {page} JSON")
			break
		if not isinstance(issues, list) or not issues:
			break
		for issue in issues:
			if fetched >= max_issues:
				break
			if not isinstance(issue, dict):
				continue
			# Skip PRs (they have a `pull_request` field)
			if "pull_request" in issue:
				continue
			num = issue.get("number")
			title = issue.get("title") or f"issue {num}"
			body = issue.get("body") or ""
			n_comments = int(issue.get("comments", 0))
			html_url = issue.get("html_url") or \
				f"https://github.com/gumyr/build123d/issues/{num}"
			user = (issue.get("user") or {}).get("login", "?")
			state = issue.get("state", "?")

			# Compose the indexed blob (prose + all comments)
			parts = [
				f"# {title}",
				f"# source: {html_url}",
				f"# state: {state}  author: {user}  comments: {n_comments}",
				"",
				"# --- issue body ---",
				body.strip(),
			]

			if include_comments and n_comments >= min_discussion_signals:
				c_url = _B3D_COMMENTS_URL.format(num=num) + "?per_page=100"
				c_raw = _http_get(c_url, accept="application/vnd.github+json")
				if c_raw is not None:
					try:
						comments = json.loads(c_raw)
					except json.JSONDecodeError:
						comments = []
					if isinstance(comments, list):
						for c in comments:
							cu = (c.get("user") or {}).get("login", "?")
							cb = c.get("body") or ""
							parts.append("")
							parts.append(f"# --- comment by {cu} ---")
							parts.append(cb.strip())
				else:
					failed.append(f"comments for #{num}")

			blob = "\n".join(parts) + "\n"
			name = f"issue_{num:04d}_{_safe_name(title)[:40]}.md"
			code_path = d / name
			try:
				code_path.write_text(blob, encoding="utf-8")
			except OSError:
				failed.append(name)
				continue
			tokens = _tokenize(blob)
			index.append({
				"name": name,
				"title": title,
				"url": html_url,
				"code_path": str(code_path),
				"issue_number": num,
				"state": state,
				"comments": n_comments,
				"tokens": tokens,
				"token_set": sorted(set(tokens)),
			})
			fetched += 1
		page += 1
		# Safety cap in case GitHub gives unexpected pagination
		if page > 20:
			break

	_index_path("build123d_discussions").write_text(
		json.dumps({"fetched_at": time.time(),
		            "source": "build123d_discussions",
		            "items": index},
		           ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"fetched": fetched,
		"failed": failed[:10],
		"index_path": str(_index_path("build123d_discussions")),
	}


# ---------------------------------------------------------------------------
# YouTube transcripts (shared scraper for cubit / build123d / gmsh)
# ---------------------------------------------------------------------------

_YT_SEARCH_URL = "https://www.youtube.com/results?search_query={q}"
_YT_VIDEO_URL = "https://www.youtube.com/watch?v={vid}"
_YT_VID_RE = re.compile(r'"videoId":"([a-zA-Z0-9_-]{11})"')
_YT_TITLE_RE = re.compile(r'"title":\{"runs":\[\{"text":"([^"]{5,120})"')


def _youtube_search(query: str, max_videos: int) -> list[tuple[str, str]]:
	"""Search YouTube via the public results page, return [(id, title)]."""
	q = urllib.parse.quote(query)
	req = urllib.request.Request(
		_YT_SEARCH_URL.format(q=q),
		headers={"User-Agent": "Mozilla/5.0 (radia-mcp)"},
	)
	try:
		with urllib.request.urlopen(req, timeout=_TIMEOUT_SECONDS) as resp:  # noqa: S310
			html = resp.read().decode("utf-8", errors="replace")
	except (urllib.error.HTTPError, urllib.error.URLError,
	        TimeoutError, OSError):
		return []
	ids = _YT_VID_RE.findall(html)
	titles = _YT_TITLE_RE.findall(html)
	out: list[tuple[str, str]] = []
	seen: set[str] = set()
	for i, t in zip(ids, titles):
		if i in seen:
			continue
		seen.add(i)
		out.append((i, t))
		if len(out) >= max_videos:
			break
	return out


def _fetch_youtube_transcript(video_id: str) -> tuple[str, list[str]]:
	"""Fetch English transcript for a video. Returns (text, languages_tried).

	Returns ('', []) when the video has no transcript / subtitles disabled
	/ youtube-transcript-api not installed.
	"""
	try:
		from youtube_transcript_api import YouTubeTranscriptApi
	except ImportError:
		return "", []
	try:
		api = YouTubeTranscriptApi()
		# Returns a FetchedTranscript or list of snippets — handle both
		fetched = api.fetch(video_id, languages=("en", "en-US", "en-GB", "ja"))
	except Exception:
		return "", ["en", "en-US", "en-GB", "ja"]
	# Snippets can be dicts or attribute objects depending on version
	parts: list[str] = []
	for snip in fetched:
		t = getattr(snip, "text", None) if not isinstance(snip, dict) \
			else snip.get("text")
		if t:
			parts.append(str(t))
	return "\n".join(parts), []


def refresh_youtube_transcripts(query: str, source_name: str,
                                  max_videos: int = 12) -> dict:
	"""Generic: search YouTube for `query`, fetch top-N transcripts,
	store as the `source_name` example sub-source.

	Used by `refresh_cubit_youtube` / `refresh_build123d_youtube` /
	`refresh_gmsh_youtube`. Optional dep: `youtube-transcript-api`
	(install via `pip install radia-mcp[youtube]`).

	Args:
	    query: search string passed to youtube.com/results.
	    source_name: index sub-source name (e.g., "cubit_youtube").
	    max_videos: cap on videos to enumerate (transcripts may be
	        unavailable on some — those are skipped).
	"""
	d = _examples_dir(source_name)
	results = _youtube_search(query, max_videos)
	if not results:
		return {"status": "error",
		        "error": "youtube search returned 0 videos",
		        "query": query}

	index: list[dict[str, Any]] = []
	fetched = 0
	skipped: list[str] = []
	for vid, title in results:
		text, _langs = _fetch_youtube_transcript(vid)
		if not text:
			skipped.append(vid)
			continue
		url = _YT_VIDEO_URL.format(vid=vid)
		blob = (
			f"# {title}\n"
			f"# source: {url}\n"
			f"# video_id: {vid}\n\n"
			f"# --- transcript ---\n"
			f"{text}\n"
		)
		name = f"{vid}_{_safe_name(title)[:50]}.md"
		code_path = d / name
		try:
			code_path.write_text(blob, encoding="utf-8")
		except OSError:
			skipped.append(vid)
			continue
		tokens = _tokenize(blob)
		index.append({
			"name": name,
			"title": title,
			"url": url,
			"code_path": str(code_path),
			"video_id": vid,
			"tokens": tokens,
			"token_set": sorted(set(tokens)),
		})
		fetched += 1

	_index_path(source_name).write_text(
		json.dumps({"fetched_at": time.time(),
		             "source": source_name,
		             "items": index,
		             "search_query": query},
		            ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"fetched": fetched,
		"skipped_no_transcript": len(skipped),
		"query": query,
		"index_path": str(_index_path(source_name)),
	}


def refresh_cubit_youtube(max_videos: int = 12) -> dict:
	return refresh_youtube_transcripts(
		query="Coreform Cubit meshing tutorial",
		source_name="cubit_youtube",
		max_videos=max_videos,
	)


def refresh_build123d_youtube(max_videos: int = 12) -> dict:
	return refresh_youtube_transcripts(
		query="build123d python CAD tutorial",
		source_name="build123d_youtube",
		max_videos=max_videos,
	)


def refresh_gmsh_youtube(max_videos: int = 12) -> dict:
	return refresh_youtube_transcripts(
		query="gmsh meshing tutorial",
		source_name="gmsh_youtube",
		max_videos=max_videos,
	)


# ---------------------------------------------------------------------------
# Coreform training .zip — bulk .jou + .sat sample download
# ---------------------------------------------------------------------------

_CUBIT_TRAINING_ZIP_URL = (
	"https://coreform.com/downloads/cubit-training/examples_only.zip"
)


def refresh_coreform_training_zip(force: bool = False) -> dict:
	"""Download Coreform's `examples_only.zip` (Cubit training pack)
	and extract under `<state_dir>/coreform_training/`. Each `.jou`
	is rewalked into the existing `cubit_local` index so it shows up
	in `cubit_examples` / `cubit_ask` automatically.
	"""
	import zipfile as _zip
	target_dir = state_dir() / "coreform_training"
	target_dir.mkdir(parents=True, exist_ok=True)
	zip_local = target_dir / "examples_only.zip"

	if force or not zip_local.exists():
		raw = _http_get(_CUBIT_TRAINING_ZIP_URL,
		                 accept="application/zip")
		if raw is None:
			return {"status": "error",
			        "error": "download failed",
			        "url": _CUBIT_TRAINING_ZIP_URL}
		try:
			zip_local.write_bytes(raw)
		except OSError as e:
			return {"status": "error", "stage": "write_zip",
			        "error": str(e)}

	# Extract
	extracted = 0
	try:
		with _zip.ZipFile(zip_local) as zf:
			members = zf.namelist()
			zf.extractall(target_dir)
			extracted = len(members)
	except _zip.BadZipFile as e:
		return {"status": "error", "stage": "unzip", "error": str(e)}

	# Re-walk into cubit_local with this dir added
	jou_count = sum(1 for p in target_dir.rglob("*.jou"))
	stat = refresh_cubit_local_examples(
		roots=[*_DEFAULT_LOCAL_CUBIT_ROOTS, str(target_dir)],
	)
	return {
		"status": "ok",
		"download_url": _CUBIT_TRAINING_ZIP_URL,
		"local_zip": str(zip_local),
		"zip_size": zip_local.stat().st_size,
		"extracted_members": extracted,
		"jou_in_pack": jou_count,
		"cubit_local_after": stat.get("indexed"),
	}


# ---------------------------------------------------------------------------
# GitHub-wide `.jou` file code search (PAT-gated)
# ---------------------------------------------------------------------------

_GH_CODE_SEARCH = ("https://api.github.com/search/code?q=extension:jou"
                    "+{q}&per_page=30")


def refresh_jou_github_code_search(extra_query: str = "cubit",
                                    max_files: int = 30) -> dict:
	"""GitHub-wide search for `.jou` files (Cubit journals) and pull
	their content. Requires a PAT (gh auth token / GITHUB_TOKEN) —
	GitHub's code search endpoint is auth-only.

	Indexed under `cubit_jou_github` sub-source; joined into the
	`cubit` family.
	"""
	from . import web_docs as _wd
	tok = _wd.github_token()
	if not tok:
		return {"status": "skipped",
		        "reason": "no GitHub PAT (GITHUB_TOKEN / GH_TOKEN / "
		                  "`gh auth login`) — code search requires auth"}
	d = _examples_dir("cubit_jou_github")
	url = _GH_CODE_SEARCH.format(q=urllib.parse.quote(extra_query))
	req = urllib.request.Request(url, headers={
		"User-Agent": _USER_AGENT,
		"Accept": "application/vnd.github+json",
		"Authorization": f"Bearer {tok}",
	})
	try:
		with urllib.request.urlopen(req, timeout=_TIMEOUT_SECONDS) as resp:  # noqa: S310
			body = json.loads(resp.read().decode("utf-8", errors="replace"))
	except (urllib.error.HTTPError, urllib.error.URLError,
	        TimeoutError, OSError, json.JSONDecodeError) as e:
		return {"status": "error", "stage": "search",
		        "error": f"{type(e).__name__}: {e}"}
	items = body.get("items", [])[: max_files]
	index: list[dict[str, Any]] = []
	failed: list[str] = []
	for it in items:
		repo = it.get("repository", {}).get("full_name", "?")
		path = it.get("path", "")
		html_url = it.get("html_url") or ""
		# Convert blob URL → raw URL for content
		raw_url = it.get("url")
		if raw_url:
			# /repos/owner/repo/contents/path?ref=sha
			c_raw = _http_get(raw_url, accept="application/vnd.github+json")
			if c_raw is None:
				failed.append(html_url)
				continue
			try:
				meta = json.loads(c_raw)
				import base64 as _b64
				content_b64 = meta.get("content", "")
				body_bytes = _b64.b64decode(content_b64)
				body_text = body_bytes.decode("utf-8", errors="replace")
			except (json.JSONDecodeError, ValueError):
				failed.append(html_url)
				continue
		else:
			continue
		blob_text = (
			f"# {repo}: {path}\n"
			f"# source: {html_url}\n\n"
			f"{body_text}\n"
		)
		name = _safe_name(f"{repo.replace('/', '_')}__{path.replace('/', '_')}")[:120]
		code_path = d / name
		try:
			code_path.write_text(blob_text, encoding="utf-8")
		except OSError:
			failed.append(html_url)
			continue
		tokens = _tokenize(blob_text)
		index.append({
			"name": name,
			"title": f"{repo}: {path}",
			"url": html_url,
			"code_path": str(code_path),
			"repo": repo,
			"path": path,
			"tokens": tokens,
			"token_set": sorted(set(tokens)),
		})
	_index_path("cubit_jou_github").write_text(
		json.dumps({"fetched_at": time.time(),
		             "source": "cubit_jou_github",
		             "items": index,
		             "query": extra_query},
		            ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"fetched": len(index),
		"failed": failed[:10],
		"index_path": str(_index_path("cubit_jou_github")),
	}


# ---------------------------------------------------------------------------
# gmsh family: GitLab issues (gmsh/gmsh) + StackOverflow [gmsh]
# ---------------------------------------------------------------------------

_GMSH_GITLAB_ISSUES = ("https://gitlab.onelab.info/api/v4/projects/3/issues"
                       "?state=all&per_page=100&order_by=updated_at&page={page}")
_GMSH_GITLAB_NOTES = ("https://gitlab.onelab.info/api/v4/projects/3/"
                       "issues/{iid}/notes?per_page=100&sort=asc")


def refresh_gmsh_issues(max_issues: int = 60,
                         include_notes: bool = True) -> dict:
	"""Fetch recent gmsh GitLab issues + comment threads (anonymous).

	Walks `https://gitlab.onelab.info/api/v4/projects/3/issues`
	(`gmsh/gmsh`, 3000+ issues total as of 2026-04). Each indexed
	entry carries: title, state, labels, issue body, comment thread.

	Anonymous GitLab API allows unauthenticated reads with moderate
	rate limits; default `max_issues=60` = 1 listing request + 60
	note requests, well within budget.

	Args:
	    max_issues: cap on issues indexed (default 60).
	    include_notes: also fetch discussion notes per issue.
	"""
	d = _examples_dir("gmsh_issues")
	index: list[dict[str, Any]] = []
	fetched = 0
	failed: list[str] = []

	page = 1
	while fetched < max_issues:
		url = _GMSH_GITLAB_ISSUES.format(page=page)
		raw = _http_get(url, accept="application/json")
		if raw is None:
			failed.append(f"issues page {page}")
			break
		try:
			issues = json.loads(raw)
		except json.JSONDecodeError:
			failed.append(f"issues page {page} json")
			break
		if not isinstance(issues, list) or not issues:
			break
		for issue in issues:
			if fetched >= max_issues:
				break
			iid = issue.get("iid")
			title = issue.get("title") or f"issue {iid}"
			description = issue.get("description") or ""
			state = issue.get("state", "?")
			html_url = issue.get("web_url") or \
				f"https://gitlab.onelab.info/gmsh/gmsh/-/issues/{iid}"
			n_notes = int(issue.get("user_notes_count", 0) or 0)
			author = (issue.get("author") or {}).get("username", "?")
			labels = issue.get("labels") or []

			parts = [
				f"# {title}",
				f"# source: {html_url}",
				f"# state: {state}  author: {author}  notes: {n_notes}",
				f"# labels: {', '.join(labels) if labels else '(none)'}",
				"",
				"# --- issue body ---",
				description.strip(),
			]

			if include_notes and n_notes > 0:
				notes_url = _GMSH_GITLAB_NOTES.format(iid=iid)
				nraw = _http_get(notes_url, accept="application/json")
				if nraw is not None:
					try:
						notes = json.loads(nraw)
					except json.JSONDecodeError:
						notes = []
					if isinstance(notes, list):
						for n in notes:
							if n.get("system"):
								continue  # skip auto label/state changes
							nu = (n.get("author") or {}).get("username", "?")
							nb = n.get("body") or ""
							parts.append("")
							parts.append(f"# --- note by {nu} ---")
							parts.append(nb.strip())
				else:
					failed.append(f"notes for #{iid}")

			blob = "\n".join(parts) + "\n"
			name = f"issue_{iid:04d}_{_safe_name(title)[:40]}.md"
			code_path = d / name
			try:
				code_path.write_text(blob, encoding="utf-8")
			except OSError:
				failed.append(name)
				continue
			tokens = _tokenize(blob)
			index.append({
				"name": name,
				"title": title,
				"url": html_url,
				"code_path": str(code_path),
				"iid": iid,
				"state": state,
				"labels": labels,
				"tokens": tokens,
				"token_set": sorted(set(tokens)),
			})
			fetched += 1
		page += 1
		if page > 20:
			break

	_index_path("gmsh_issues").write_text(
		json.dumps({"fetched_at": time.time(),
		             "source": "gmsh_issues",
		             "items": index},
		            ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"fetched": fetched,
		"failed": failed[:10],
		"index_path": str(_index_path("gmsh_issues")),
	}


_GMSH_SE_QUESTIONS = (
	"https://api.stackexchange.com/2.3/questions"
	"?tagged=gmsh&site={site}&pagesize={ps}&order=desc&sort=votes"
	"&filter=withbody"
)
_GMSH_SE_ANSWERS = (
	"https://api.stackexchange.com/2.3/questions/{qid}/answers"
	"?site={site}&filter=withbody&order=desc&sort=votes"
)


def _fetch_se_gzip(url: str) -> dict | None:
	"""StackExchange gzip-encoded JSON fetcher."""
	import gzip as _gz
	req = urllib.request.Request(url, headers={
		"User-Agent": _USER_AGENT if False else "radia-mcp",
		"Accept": "application/json",
	})
	try:
		with urllib.request.urlopen(req, timeout=_TIMEOUT_SECONDS) as resp:  # noqa: S310
			raw = resp.read()
			if resp.headers.get("Content-Encoding") == "gzip":
				raw = _gz.decompress(raw)
	except (urllib.error.HTTPError, urllib.error.URLError,
	        TimeoutError, OSError):
		return None
	try:
		return json.loads(raw.decode("utf-8", errors="replace"))
	except json.JSONDecodeError:
		return None


def refresh_gmsh_stackoverflow(max_questions: int = 30,
                                 include_scicomp: bool = True) -> dict:
	"""Fetch top-voted StackOverflow (and optionally SciComp.SE)
	questions tagged `gmsh` + their accepted/top answers.

	The StackExchange API is anonymous (300 req/day shared), and
	returns HTML bodies — we strip tags for the indexed blob.
	"""
	import re as _re
	import html as _html
	d = _examples_dir("gmsh_stackoverflow")
	tag_re = _re.compile(r"<[^>]+>")

	def _strip(s: str) -> str:
		return _html.unescape(tag_re.sub("", s or "")).strip()

	sites = ["stackoverflow"]
	if include_scicomp:
		sites.append("scicomp")

	index: list[dict[str, Any]] = []
	fetched = 0
	failed: list[str] = []

	for site in sites:
		url = _GMSH_SE_QUESTIONS.format(site=site,
		                                  ps=min(max_questions, 30))
		data = _fetch_se_gzip(url)
		if not data or "items" not in data:
			failed.append(f"questions on {site}")
			continue
		for item in data["items"][: max_questions]:
			if fetched >= max_questions * len(sites):
				break
			qid = item.get("question_id")
			title = item.get("title") or f"question {qid}"
			body = _strip(item.get("body") or "")
			score = int(item.get("score", 0) or 0)
			n_ans = int(item.get("answer_count", 0) or 0)
			has_accepted = bool(item.get("accepted_answer_id"))
			link = item.get("link") or f"https://{site}.com/q/{qid}"

			parts = [
				f"# {title}",
				f"# source: {link}",
				f"# site: {site}  score: {score}  answers: {n_ans}  "
				f"accepted: {has_accepted}",
				"",
				"# --- question ---",
				body,
			]

			# Fetch answers (fold top into body)
			a_url = _GMSH_SE_ANSWERS.format(qid=qid, site=site)
			adata = _fetch_se_gzip(a_url)
			if adata and "items" in adata:
				for a in adata["items"][:5]:
					ab = _strip(a.get("body") or "")
					ascore = int(a.get("score", 0) or 0)
					accepted = bool(a.get("is_accepted"))
					tag = "ACCEPTED ANSWER" if accepted else "answer"
					parts.append("")
					parts.append(f"# --- {tag} (score {ascore}) ---")
					parts.append(ab)

			blob = "\n".join(parts) + "\n"
			name = f"so_{qid}_{_safe_name(title)[:40]}.md"
			code_path = d / name
			try:
				code_path.write_text(blob, encoding="utf-8")
			except OSError:
				failed.append(name)
				continue
			tokens = _tokenize(blob)
			index.append({
				"name": name,
				"title": title,
				"url": link,
				"code_path": str(code_path),
				"site": site,
				"score": score,
				"answer_count": n_ans,
				"tokens": tokens,
				"token_set": sorted(set(tokens)),
			})
			fetched += 1

	_index_path("gmsh_stackoverflow").write_text(
		json.dumps({"fetched_at": time.time(),
		             "source": "gmsh_stackoverflow",
		             "items": index},
		            ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"fetched": fetched,
		"failed": failed[:10],
		"index_path": str(_index_path("gmsh_stackoverflow")),
	}


_B3D_GRAPHQL_DISCUSSIONS = """
query($cursor: String) {
  repository(owner: "gumyr", name: "build123d") {
    discussions(first: 50, after: $cursor,
                orderBy: {field: UPDATED_AT, direction: DESC}) {
      pageInfo { hasNextPage endCursor }
      nodes {
        number
        title
        bodyText
        url
        updatedAt
        category { name }
        answer { bodyText author { login } }
        comments(first: 30) {
          nodes { bodyText author { login } }
        }
      }
    }
  }
}
"""


def refresh_build123d_github_discussions(max_discussions: int = 100) -> dict:
	"""Fetch the actual GitHub Discussions (not just Issues) from
	`gumyr/build123d` via GraphQL. Requires a PAT — discovered via
	`github_token()` helper. Skipped (returns `status='skipped'`) if
	no PAT is available.
	"""
	import json as _json
	from . import web_docs as _wd
	tok = _wd.github_token()
	if not tok:
		return {"status": "skipped",
		        "reason": "no GitHub PAT (set GITHUB_TOKEN or run `gh auth login`)"}

	d = _examples_dir("build123d_github_discussions")
	index: list[dict[str, Any]] = []
	fetched = 0
	failed: list[str] = []
	cursor = None
	while fetched < max_discussions:
		payload = _json.dumps({
			"query": _B3D_GRAPHQL_DISCUSSIONS,
			"variables": {"cursor": cursor},
		}).encode("utf-8")
		req = urllib.request.Request(
			"https://api.github.com/graphql",
			data=payload,
			headers={
				"User-Agent": "radia-mcp/0.16",
				"Authorization": f"Bearer {tok}",
				"Content-Type": "application/json",
				"Accept": "application/vnd.github+json",
			},
			method="POST",
		)
		try:
			with urllib.request.urlopen(req, timeout=20) as resp:  # noqa: S310
				body = _json.loads(resp.read().decode("utf-8"))
		except (urllib.error.HTTPError, urllib.error.URLError,
		        TimeoutError, OSError, _json.JSONDecodeError) as e:
			return {"status": "error", "error": f"GraphQL: {e}"}
		if "errors" in body:
			return {"status": "error",
			        "error": str(body["errors"])[:400]}
		conn = body.get("data", {}).get("repository", {}).get("discussions", {})
		nodes = conn.get("nodes", [])
		if not nodes:
			break
		for disc in nodes:
			if fetched >= max_discussions:
				break
			num = disc.get("number")
			title = disc.get("title") or f"discussion {num}"
			url = disc.get("url") or f"https://github.com/gumyr/build123d/discussions/{num}"
			body_text = disc.get("bodyText") or ""
			cat = (disc.get("category") or {}).get("name", "?")
			answer = disc.get("answer") or {}
			comments = (disc.get("comments") or {}).get("nodes", [])
			parts = [
				f"# {title}",
				f"# source: {url}",
				f"# category: {cat}",
				"",
				"# --- body ---",
				body_text.strip(),
			]
			if answer and answer.get("bodyText"):
				au = (answer.get("author") or {}).get("login", "?")
				parts.append("")
				parts.append(f"# --- ACCEPTED ANSWER by {au} ---")
				parts.append(answer["bodyText"].strip())
			for c in comments:
				cb = c.get("bodyText") or ""
				cu = (c.get("author") or {}).get("login", "?")
				parts.append("")
				parts.append(f"# --- comment by {cu} ---")
				parts.append(cb.strip())
			blob = "\n".join(parts) + "\n"
			name = f"discussion_{num:04d}_{_safe_name(title)[:40]}.md"
			code_path = d / name
			try:
				code_path.write_text(blob, encoding="utf-8")
			except OSError:
				failed.append(name)
				continue
			tokens = _tokenize(blob)
			index.append({
				"name": name,
				"title": title,
				"url": url,
				"code_path": str(code_path),
				"discussion_number": num,
				"category": cat,
				"has_answer": bool(answer),
				"comments": len(comments),
				"tokens": tokens,
				"token_set": sorted(set(tokens)),
			})
			fetched += 1
		info = conn.get("pageInfo", {})
		if not info.get("hasNextPage"):
			break
		cursor = info.get("endCursor")

	_index_path("build123d_github_discussions").write_text(
		_json.dumps({"fetched_at": time.time(),
		             "source": "build123d_github_discussions",
		             "items": index},
		            ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"fetched": fetched,
		"failed": failed[:10],
		"index_path": str(_index_path("build123d_github_discussions")),
	}


# ---------------------------------------------------------------------------
# bd_warehouse: GitHub `gumyr/bd_warehouse` (industrial-parts library)
# ---------------------------------------------------------------------------

def refresh_bd_warehouse_examples(include_examples_dir: bool = True) -> dict:
	"""Fetch bd_warehouse source modules + examples from GitHub.

	bd_warehouse is gumyr's industrial-parts companion to build123d
	(bearings, fasteners, flanges, gears, open_builds extrusions,
	pipes, sprockets). Each module is a single file with parametric
	classes. Scraped the same way as build123d examples.

	Args:
	    include_examples_dir: also index `examples/` subdir if present.
	"""
	d = _examples_dir("bd_warehouse")
	index: list[dict[str, Any]] = []
	fetched = 0
	failed: list[str] = []

	endpoints: list[tuple[str, str]] = [("src/bd_warehouse", _BD_WAREHOUSE_LIST_URL)]
	if include_examples_dir:
		endpoints.append(("examples", _BD_WAREHOUSE_EXAMPLES_URL))

	for subdir, url in endpoints:
		raw = _http_get(url, accept="application/vnd.github+json")
		if not raw:
			failed.append(f"listing: {subdir}")
			continue
		try:
			listing = json.loads(raw)
		except json.JSONDecodeError:
			failed.append(f"listing json: {subdir}")
			continue
		if not isinstance(listing, list):
			# Some sub-dirs (e.g. examples) might not exist — skip
			continue
		for f in listing:
			if not isinstance(f, dict):
				continue
			name = f.get("name", "")
			if not name.endswith(".py"):
				continue
			dl = f.get("download_url")
			if not dl:
				failed.append(name)
				continue
			body_raw = _http_get(dl, accept="text/plain")
			if body_raw is None:
				failed.append(name)
				continue
			body = body_raw.decode("utf-8", errors="replace")

			stored_name = f"{subdir.replace('/', '_')}__{_safe_name(name)}"
			code_path = d / stored_name
			try:
				code_path.write_text(body, encoding="utf-8")
			except OSError:
				failed.append(name)
				continue
			title = _guess_python_title(body, fallback=name)
			tokens = _tokenize(body)
			index.append({
				"name": name,
				"stored": stored_name,
				"title": title,
				"url": f.get("html_url")
				       or f"https://github.com/gumyr/bd_warehouse/blob/main/{subdir}/{name}",
				"code_path": str(code_path),
				"size": f.get("size", 0),
				"subdir": subdir,
				"tokens": tokens,
				"token_set": sorted(set(tokens)),
			})
			fetched += 1

	_index_path("bd_warehouse").write_text(
		json.dumps({"fetched_at": time.time(), "source": "bd_warehouse",
		            "items": index}, ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"fetched": fetched,
		"failed": failed,
		"index_path": str(_index_path("bd_warehouse")),
	}


# ---------------------------------------------------------------------------
# cubit_local: walks on-disk directories (curated corpus + Radia durable lanes)
# ---------------------------------------------------------------------------

_LOCAL_EXTS = (".jou", ".py")
_LOCAL_SKIP_DIRS = frozenset({
	".git", "__pycache__", "build", "dist", ".venv", "venv", "node_modules",
	"packages",  # don't re-index radia-mcp's own source
})
_LOCAL_MAX_FILE_BYTES = 200_000  # skip huge files (data dumps, generated meshes)


def _resolve_local_root(root_s: str) -> Path:
	"""Resolve `repo:/...` roots from an editable Radia checkout."""
	if root_s.startswith("repo:/"):
		rel = root_s[len("repo:/"):].lstrip("/\\")
		here = Path(__file__).resolve()
		candidates = []
		if len(here.parents) > 5:
			candidates.append(here.parents[5] / rel)
		candidates.append(Path.cwd() / rel)
		for candidate in candidates:
			if candidate.exists():
				return candidate
	return Path(root_s)


def _guess_local_title(body: str, path: Path) -> str:
	"""Pull a 1-line title from a .jou or .py file."""
	for raw_line in body.splitlines()[:30]:
		line = raw_line.strip()
		if not line:
			continue
		# Python docstring: opening triple-quote then text
		if line.startswith(("'''", '"""')) and len(line) > 3:
			return line.strip("'\" ")[:120]
		# Python comment / Cubit journal comment (both start with #)
		if line.startswith("#"):
			t = line.lstrip("# \t").strip()
			if t and not t.startswith("!"):
				return t[:120]
		# First meaningful non-comment line (useful when files have no header)
		return line[:120]
	return path.name


def refresh_cubit_local_examples(roots: list[str] | None = None,
                                 extra_skip_dirs: list[str] | None = None) -> dict:
	"""Walk local directories and index .jou / .py files as Cubit examples.

	Default roots: `public-safe curated corpus` (the lab's years of curated
	Cubit projects, ~145 files), plus Radia's durable docs/,
	validation_test/, and panel-sample lanes. Users can pass `roots=[...]`
	to override.

	Each file:
	  - title from first non-empty comment / docstring line
	  - tokens from full body (tf-idf indexable)
	  - path kept as absolute — NOT copied to state_dir (these are
	    already on a shared drive; duplicating wastes space and drifts
	    as files evolve)

	Files larger than 200 kB are skipped to keep the index sane.
	"""
	roots = roots or _DEFAULT_LOCAL_CUBIT_ROOTS
	skip_dirs = set(_LOCAL_SKIP_DIRS)
	if extra_skip_dirs:
		skip_dirs.update(extra_skip_dirs)

	index: list[dict[str, Any]] = []
	indexed_paths: set[str] = set()
	skipped_big = 0
	read_errors = 0
	roots_walked: list[str] = []

	for root_s in roots:
		root = _resolve_local_root(root_s)
		if not root.exists():
			continue
		roots_walked.append(str(root))
		for path in root.rglob("*"):
			if not path.is_file():
				continue
			if path.suffix.lower() not in _LOCAL_EXTS:
				continue
			# Skip any path with a component in skip_dirs
			if any(part in skip_dirs for part in path.parts):
				continue
			try:
				size = path.stat().st_size
			except OSError:
				continue
			if size > _LOCAL_MAX_FILE_BYTES or size == 0:
				if size > _LOCAL_MAX_FILE_BYTES:
					skipped_big += 1
				continue
			p_abs = str(path.resolve())
			if p_abs in indexed_paths:
				continue
			try:
				body = path.read_text(encoding="utf-8", errors="replace")
			except OSError:
				read_errors += 1
				continue
			tokens = _tokenize(body)
			if not tokens:
				continue
			rel = str(path).replace("\\", "/")
			title = _guess_local_title(body, path)
			index.append({
				"name": path.name,
				"title": title,
				"url": f"file:///{rel}",
				"code_path": p_abs,
				"size": size,
				"root": str(root),
				"parent_dir": path.parent.name,
				"ext": path.suffix.lower(),
				"tokens": tokens,
				"token_set": sorted(set(tokens)),
			})
			indexed_paths.add(p_abs)

	_index_path("cubit_local").write_text(
		json.dumps({"fetched_at": time.time(), "source": "cubit_local",
		            "items": index, "roots_walked": roots_walked},
		           ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"indexed": len(index),
		"roots_walked": roots_walked,
		"skipped_too_big": skipped_big,
		"read_errors": read_errors,
		"index_path": str(_index_path("cubit_local")),
	}


# ---------------------------------------------------------------------------
# cubit: Discourse forum (forum.coreform.com)
# ---------------------------------------------------------------------------

_CUBIT_FORUM_BASE = "https://forum.coreform.com"
_CUBIT_QUERIES = [
	"tutorial example mesh",
	"hex meshing tutorial",
	"sweep scheme example",
	"journal file example",
	"export example",
	"boundary layer mesh",
	"thin shell mesh",
	"mesh quality metric",
	"mesh refinement",
	"abaqus export",
	"high order curving",
	"polyhedron scheme",
	"multi sweep source target",
	"imprint merge",
	"webcut journal",
]


def _walk_forum_latest(max_pages: int = 30) -> list[int]:
	"""Walk `/latest.json` pagewise, collecting all topic IDs.

	Returns topic IDs in recency order. Stops at an empty page or at
	`max_pages` (safety cap — Coreform currently has ~263 topics, so
	30 pages × 30 topics/page is comfortably enough headroom).
	"""
	ids: list[int] = []
	for page in range(max_pages):
		url = (f"{_CUBIT_FORUM_BASE}/latest.json"
		       f"?no_definitions=true&page={page}")
		raw = _http_get(url, accept="application/json")
		if raw is None:
			break
		try:
			data = json.loads(raw)
		except json.JSONDecodeError:
			break
		topics = data.get("topic_list", {}).get("topics", [])
		if not topics:
			break
		for t in topics:
			if isinstance(t, dict):
				tid = t.get("id")
				if isinstance(tid, int):
					ids.append(tid)
	# dedupe preserving order
	seen: set[int] = set()
	out: list[int] = []
	for tid in ids:
		if tid not in seen:
			seen.add(tid)
			out.append(tid)
	return out


def _fetch_topic_posts(topic_id: int) -> list[dict]:
	"""Fetch a topic's full post stream via Discourse's t/<id>.json."""
	url = f"{_CUBIT_FORUM_BASE}/t/{topic_id}.json?include_raw=1"
	raw = _http_get(url, accept="application/json")
	if raw is None:
		return []
	try:
		data = json.loads(raw)
	except json.JSONDecodeError:
		return []
	posts = data.get("post_stream", {}).get("posts", [])
	# Attach topic metadata onto each post entry for convenience
	meta = {
		"title": data.get("title", ""),
		"slug": data.get("slug", ""),
		"topic_id": topic_id,
	}
	out = []
	for p in posts:
		if isinstance(p, dict):
			p = dict(p)
			p["_topic"] = meta
			out.append(p)
	return out


def _refresh_cubit_forum_full(max_topics: int = 300,
                               concurrency: int = 8) -> dict:
	"""Walk /latest.json across all pages, fetch every topic's posts in
	parallel, extract code fences, and index as cubit_examples.

	Preferred over `refresh_cubit_examples` (seed-query based) when you
	want full-archive coverage. Coreform forum is small (~263 topics as
	of 2026-04), so walking the whole thing with a small thread pool
	takes <30 s.
	"""
	import concurrent.futures as _cf

	d = _examples_dir("cubit")
	tids = _walk_forum_latest()
	if max_topics > 0:
		tids = tids[:max_topics]

	# Fetch all topics in parallel — each is one HTTP round-trip
	topic_posts: dict[int, list[dict]] = {}
	failed: list[str] = []
	with _cf.ThreadPoolExecutor(max_workers=max(1, concurrency)) as pool:
		futures = {pool.submit(_fetch_topic_posts, tid): tid for tid in tids}
		for fut in _cf.as_completed(futures):
			tid = futures[fut]
			try:
				posts = fut.result()
			except Exception:
				failed.append(f"topic {tid}")
				continue
			if posts:
				topic_posts[tid] = posts
			else:
				failed.append(f"topic {tid}")

	index: list[dict[str, Any]] = []
	fetched = 0
	seen_post_ids: set[int] = set()
	# Preserve original ordering (most-recent first)
	for tid in tids:
		posts = topic_posts.get(tid)
		if not posts:
			continue
		for post in posts:
			pid = post.get("id")
			if not pid or pid in seen_post_ids:
				continue
			seen_post_ids.add(pid)
			body = post.get("raw") or post.get("cooked") or ""
			fences = _CODE_FENCE_RE.findall(body)
			if not fences:
				continue
			code = "\n\n# ---- next snippet ----\n\n".join(
				f.strip() for f in fences if f.strip()
			)
			if not code:
				continue
			meta = post.get("_topic", {})
			title = meta.get("title", "") or f"forum post {pid}"
			slug = meta.get("slug", "")
			pnum = post.get("post_number", 1)
			html_url = (f"{_CUBIT_FORUM_BASE}/t/{slug}/{tid}/{pnum}"
			            if slug else f"{_CUBIT_FORUM_BASE}/p/{pid}")
			name = f"{_safe_name(title)[:60]}__p{pid}.jou"
			code_path = d / name
			prose_head = body.splitlines()[:60]
			blob = (f"# {title}\n# source: {html_url}\n# post id: {pid}\n\n"
			        f"# --- prose ---\n"
			        + "\n".join("# " + ln for ln in prose_head)
			        + "\n\n# --- code ---\n" + code + "\n")
			try:
				code_path.write_text(blob, encoding="utf-8")
			except OSError:
				failed.append(name)
				continue
			tokens = _tokenize(blob)
			index.append({
				"name": name,
				"title": title,
				"url": html_url,
				"code_path": str(code_path),
				"post_id": pid,
				"tokens": tokens,
				"token_set": sorted(set(tokens)),
			})
			fetched += 1
	_index_path("cubit").write_text(
		json.dumps({"fetched_at": time.time(), "source": "cubit",
		            "items": index, "topics_walked": len(tids)},
		           ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"fetched": fetched,
		"failed": failed[:10],
		"topics_walked": len(tids),
		"posts_seen": len(seen_post_ids),
		"concurrency": concurrency,
		"index_path": str(_index_path("cubit")),
	}


def refresh_cubit_examples(limit_per_query: int = 5,
                           full_walk: bool = True,
                           max_topics: int = 300) -> dict:
	"""Fetch Cubit example snippets from the Coreform forum.

	Default (`full_walk=True`): walks the entire forum via
	/latest.json pagination, covering every topic. Best recall for a
	small forum (~263 topics as of 2026-04). Falls back to the
	legacy seed-query behavior if `full_walk=False`.

	Args:
	    limit_per_query: only used when full_walk=False (seed queries).
	    full_walk: if True, walk /latest.json across pages.
	    max_topics: cap on topics walked (safety).

	Returns stats dict.
	"""
	if full_walk:
		return _refresh_cubit_forum_full(max_topics=max_topics)
	d = _examples_dir("cubit")
	seen_post_ids: set[int] = set()
	index: list[dict[str, Any]] = []
	fetched = 0
	failed: list[str] = []

	for q in _CUBIT_QUERIES:
		url = (f"{_CUBIT_FORUM_BASE}/search.json?q="
		       + urllib.parse.quote(q))
		raw = _http_get(url, accept="application/json")
		if raw is None:
			failed.append(f"search: {q}")
			continue
		try:
			search_json = json.loads(raw)
		except json.JSONDecodeError:
			failed.append(f"search json: {q}")
			continue
		topics_by_id = {t["id"]: t for t in search_json.get("topics", [])
		                if isinstance(t, dict)}
		posts = search_json.get("posts", [])[:limit_per_query]
		for post in posts:
			pid = post.get("id")
			if not pid or pid in seen_post_ids:
				continue
			seen_post_ids.add(pid)
			post_url = f"{_CUBIT_FORUM_BASE}/posts/{pid}.json"
			praw = _http_get(post_url, accept="application/json")
			if praw is None:
				failed.append(f"post {pid}")
				continue
			try:
				pdata = json.loads(praw)
			except json.JSONDecodeError:
				failed.append(f"post json {pid}")
				continue
			body = pdata.get("raw") or ""
			fences = _CODE_FENCE_RE.findall(body)
			if not fences:
				continue
			code = "\n\n# ---- next snippet ----\n\n".join(
				f.strip() for f in fences if f.strip()
			)
			if not code:
				continue

			topic = topics_by_id.get(post.get("topic_id"), {})
			title = topic.get("title", "") or f"forum post {pid}"
			slug = topic.get("slug", "")
			tid = topic.get("id", post.get("topic_id"))
			pnum = post.get("post_number", 1)
			html_url = (f"{_CUBIT_FORUM_BASE}/t/{slug}/{tid}/{pnum}"
			            if slug and tid else f"{_CUBIT_FORUM_BASE}/p/{pid}")

			name = f"{_safe_name(title)[:60]}__p{pid}.jou"
			code_path = d / name
			# Store both prose and code so retrieval indexes both
			blob = (f"# {title}\n# source: {html_url}\n# "
			        f"post id: {pid}\n\n# --- prose ---\n"
			        + "\n".join("# " + line for line in body.splitlines()[:60])
			        + "\n\n# --- code ---\n" + code + "\n")
			try:
				code_path.write_text(blob, encoding="utf-8")
			except OSError:
				failed.append(name)
				continue

			tokens = _tokenize(blob)
			index.append({
				"name": name,
				"title": title,
				"url": html_url,
				"code_path": str(code_path),
				"post_id": pid,
				"tokens": tokens,
				"token_set": sorted(set(tokens)),
			})
			fetched += 1

	index_data = {
		"fetched_at": time.time(),
		"source": "cubit",
		"items": index,
	}
	_index_path("cubit").write_text(
		json.dumps(index_data, ensure_ascii=False, indent=2),
		encoding="utf-8",
	)
	return {
		"status": "ok",
		"fetched": fetched,
		"failed": failed,
		"posts_seen": len(seen_post_ids),
		"index_path": str(_index_path("cubit")),
	}


# ---------------------------------------------------------------------------
# Retrieval
# ---------------------------------------------------------------------------

def _load_index(source: str) -> dict:
	p = _index_path(source)
	if not p.exists():
		return {"fetched_at": 0, "items": []}
	try:
		return json.loads(p.read_text(encoding="utf-8"))
	except (OSError, json.JSONDecodeError):
		return {"fetched_at": 0, "items": []}


def _is_stale(index: dict, ttl: int = _REFRESH_TTL_SECONDS) -> bool:
	return time.time() - float(index.get("fetched_at", 0)) > ttl


REFRESH_FUNCS = {
	"build123d": refresh_build123d_examples,
	"bd_warehouse": refresh_bd_warehouse_examples,
	"build123d_discussions": refresh_build123d_discussions,
	"build123d_github_discussions": refresh_build123d_github_discussions,
	"build123d_youtube": refresh_build123d_youtube,
	"cubit": refresh_cubit_examples,
	"cubit_local": refresh_cubit_local_examples,
	"cubit_youtube": refresh_cubit_youtube,
	"cubit_jou_github": refresh_jou_github_code_search,
	"gmsh_issues": refresh_gmsh_issues,
	"gmsh_stackoverflow": refresh_gmsh_stackoverflow,
	"gmsh_youtube": refresh_gmsh_youtube,
}


def _resolve_family(name: str) -> list[str]:
	"""Return the list of sub-sources for a family, or [name] if not a family."""
	return FAMILIES.get(name, [name])


def _ensure_fresh(source: str, force: bool) -> None:
	if force or _is_stale(_load_index(source)):
		fn = REFRESH_FUNCS.get(source)
		if fn is not None:
			try:
				fn()
			except Exception:
				# Refresh failure should not break search of the remaining sources
				pass


def search_examples(source: str, query: str, limit: int = 5,
                    force_refresh: bool = False,
                    auto_refresh_if_empty: bool = True) -> dict:
	"""tf-idf search across cached examples for `source`.

	`source` may be a single concrete source (`build123d`, `cubit`,
	`bd_warehouse`, `cubit_local`) or a family name (`build123d` resolves
	to `[build123d, bd_warehouse]`; `cubit` resolves to `[cubit,
	cubit_local]`).  For families the indexes are unioned and ranked
	together so one query sees GitHub + local + forum hits at once.

	Args:
	    source: concrete source name or family name.
	    query: keywords (AND-scored, stop-words stripped).
	    limit: cap on returned examples.
	    force_refresh: re-scrape ignoring cache.
	    auto_refresh_if_empty: if any sub-source's cache is empty or
	        stale, refresh that one before searching. Default True.

	Returns dict with `results` (each hit: source / name / title / url /
	score / excerpt / path) + per-sub-source `stats`.
	"""
	sub_sources = _resolve_family(source)
	if not sub_sources:
		return {"status": "error", "error": f"unknown source {source!r}",
		        "known": list(REFRESH_FUNCS.keys())}

	for sub in sub_sources:
		if force_refresh or auto_refresh_if_empty:
			_ensure_fresh(sub, force=force_refresh)

	combined: list[tuple[dict, str]] = []  # (item, owning_sub_source)
	stats: dict[str, dict] = {}
	for sub in sub_sources:
		idx = _load_index(sub)
		items = idx.get("items", [])
		stats[sub] = {"indexed": len(items),
		              "fetched_at": idx.get("fetched_at")}
		for it in items:
			combined.append((it, sub))

	if not combined:
		return {"status": "ok", "source": source, "query": query,
		        "sub_sources": sub_sources, "stats": stats,
		        "results": [], "note": "empty index"}

	q_terms = [t for t in _tokenize(query or "")]
	if not q_terms:
		return {"status": "error", "error": "empty query"}

	n = max(1, len(combined))
	idf: dict[str, float] = {}
	for t in q_terms:
		df = sum(1 for it, _ in combined if t in (it.get("token_set") or []))
		idf[t] = math.log((n + 1) / (df + 1)) + 1.0

	scored: list[tuple[float, dict, str]] = []
	for it, sub in combined:
		tokens = it.get("tokens") or []
		if not tokens:
			continue
		tf: dict[str, int] = {}
		for tok in tokens:
			tf[tok] = tf.get(tok, 0) + 1
		inv_len = 1.0 / (1.0 + math.log(1 + len(tokens)))
		score = 0.0
		hits = 0
		for term in q_terms:
			f = tf.get(term, 0)
			if f == 0:
				continue
			hits += 1
			score += (1.0 + math.log(f)) * idf.get(term, 1.0) * inv_len
		if hits == 0:
			continue
		title_lc = (it.get("title") or "").lower()
		head_hits = sum(1 for t in q_terms if t in title_lc)
		if head_hits:
			score *= 1.0 + 0.75 * head_hits
		score *= 1.0 + 0.25 * (hits / max(1, len(q_terms)))
		# Mild local-source boost: lab-curated examples are usually more
		# relevant to Radia users than random forum posts / general libs.
		if sub in ("cubit_local",):
			score *= 1.1
		scored.append((score, it, sub))

	scored.sort(key=lambda x: (-x[0], x[2], x[1].get("name", "")))
	top = scored[: max(1, int(limit))]

	results: list[dict] = []
	for score, it, sub in top:
		excerpt = ""
		p = Path(it.get("code_path", ""))
		if p.exists():
			try:
				excerpt = p.read_text(encoding="utf-8", errors="replace")[:1600]
			except OSError:
				pass
		results.append({
			"source": sub,
			"name": it.get("name"),
			"title": it.get("title"),
			"url": it.get("url"),
			"score": round(score, 3),
			"parent_dir": it.get("parent_dir"),  # local only
			"excerpt": excerpt,
			"path": it.get("code_path"),
		})

	return {
		"status": "ok",
		"source": source,
		"sub_sources": sub_sources,
		"query": query,
		"stats": stats,
		"total_scored": len(scored),
		"results": results,
	}
