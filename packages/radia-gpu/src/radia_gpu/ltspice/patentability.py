"""Public prior-art search helpers for circuit invention triage.

These helpers do not perform a legal patentability opinion.  They create a
repeatable search plan that an engineer or agent can use with Google Scholar,
Google Patents, J-PlatPat, and general web searches before consulting a patent
professional.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable


JPLATPAT_ENTRY = "https://www.j-platpat.inpit.go.jp/?uri=%2Fs0000%2Fen"
JPO_JPLATPAT_HELP = "https://www.jpo.go.jp/e/support/j_platpat/patent_search.html"
GOOGLE_PATENTS_ADVANCED = "https://patents.google.com/advanced"
GOOGLE_SCHOLAR = "https://scholar.google.com/"

_JA_CIRCUIT_TERMS = {
    "circuit": ["回路", "電気回路"],
    "converter": ["コンバータ", "電力変換"],
    "boost": ["昇圧", "ブースト"],
    "buck": ["降圧", "バック"],
    "snubber": ["スナバ", "スナバー"],
    "switching": ["スイッチング"],
    "ripple": ["リップル"],
    "damping": ["減衰", "ダンピング"],
    "soft start": ["ソフトスタート", "起動"],
}


@dataclass(frozen=True)
class PriorArtSearchPlan:
    """Search query bundle for non-legal patentability triage."""

    title: str
    features: list[str]
    effects: list[str]
    domains: list[str]
    google_scholar: list[str]
    google_patents: list[str]
    jplatpat_keywords_en: list[str]
    jplatpat_keywords_ja: list[str]
    web: list[str]
    report_questions: list[str]
    source_urls: dict[str, str]
    disclaimer: str

    def to_dict(self) -> dict:
        data = asdict(self)
        data["query_sets"] = {
            "google_scholar": self.google_scholar,
            "google_patents": self.google_patents,
            "jplatpat_keywords_en": self.jplatpat_keywords_en,
            "jplatpat_keywords_ja": self.jplatpat_keywords_ja,
            "web": self.web,
        }
        return data


def patentability_search_plan(
    title: str,
    features: Iterable[str],
    effects: Iterable[str] | None = None,
    domains: Iterable[str] | None = None,
    include_japanese: bool = True,
) -> dict:
    """Build public search queries for novelty / inventive-step triage.

    Args:
        title: Short invention name.
        features: Essential structural or functional features.
        effects: Technical effects, problems solved, or advantages.
        domains: Technical fields such as "power electronics" or "sensor".
        include_japanese: Include Japanese J-PlatPat keyword expansions.

    Returns:
        Dict containing query strings and report prompts for use with Google
        Scholar, Google Patents, J-PlatPat, and general web search.
    """
    feature_list = _clean_list(features)
    effect_list = _clean_list(effects or [])
    domain_list = _clean_list(domains or [])
    title = title.strip()
    if not feature_list:
        raise ValueError("at least one essential feature is required")

    feature_blob = " ".join(_quote(item) for item in feature_list[:4])
    effect_blob = " ".join(_quote(item) for item in effect_list[:3])
    domain_blob = " ".join(domain_list[:3])
    pairs = [
        f"{_quote(left)} {_quote(right)}"
        for idx, left in enumerate(feature_list)
        for right in feature_list[idx + 1:]
    ]

    google_scholar = _dedupe([
        f"{feature_blob} {domain_blob}",
        f"{feature_blob} {effect_blob}",
        *[f"{pair} {domain_blob}" for pair in pairs[:8]],
    ])
    google_patents = _dedupe([
        f"({' '.join(feature_list[:4])}) ({' '.join(effect_list[:3])})",
        _quote(title) if title else feature_blob,
        *pairs[:8],
    ])
    jplatpat_en = _dedupe([
        " ".join(feature_list[:4]),
        " ".join(feature_list[:3] + effect_list[:3]),
    ])
    jplatpat_ja = (
        _dedupe(_expand_japanese_terms(feature_list + effect_list))
        if include_japanese else []
    )
    web = _dedupe([
        f"{feature_blob} filetype:pdf",
        f"{feature_blob} datasheet OR manual OR application note",
        f"{feature_blob} github OR repository OR thesis",
        f"{feature_blob} standard OR specification",
    ])

    plan = PriorArtSearchPlan(
        title=title,
        features=feature_list,
        effects=effect_list,
        domains=domain_list,
        google_scholar=[q.strip() for q in google_scholar if q.strip()],
        google_patents=[q.strip() for q in google_patents if q.strip()],
        jplatpat_keywords_en=[q.strip() for q in jplatpat_en if q.strip()],
        jplatpat_keywords_ja=jplatpat_ja,
        web=[q.strip() for q in web if q.strip()],
        report_questions=[
            "Does one reference disclose every essential feature?",
            "What is the closest prior art and what feature is missing?",
            "Would the missing feature be a routine substitution or predictable optimization?",
            "What technical effect should be proven experimentally before filing?",
            "Which claim elements must be narrowed to avoid the closest reference?",
        ],
        source_urls={
            "google_scholar": GOOGLE_SCHOLAR,
            "google_patents_advanced": GOOGLE_PATENTS_ADVANCED,
            "jplatpat": JPLATPAT_ENTRY,
            "jpo_jplatpat_help": JPO_JPLATPAT_HELP,
        },
        disclaimer=(
            "This is a search-plan aid for prior-art triage, not a legal "
            "opinion or filing clearance."
        ),
    )
    return plan.to_dict()


def _clean_list(items: Iterable[str]) -> list[str]:
    out: list[str] = []
    for item in items:
        text = str(item).strip()
        if text:
            out.append(text)
    return _dedupe(out)


def _quote(text: str) -> str:
    text = text.strip()
    if not text:
        return ""
    if text.startswith('"') and text.endswith('"'):
        return text
    return f'"{text}"'


def _dedupe(items: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        key = item.strip()
        if key and key not in seen:
            seen.add(key)
            out.append(key)
    return out


def _expand_japanese_terms(items: Iterable[str]) -> list[str]:
    terms: list[str] = []
    for item in items:
        lower = item.lower()
        for key, values in _JA_CIRCUIT_TERMS.items():
            if key in lower:
                terms.extend(values)
    return terms
