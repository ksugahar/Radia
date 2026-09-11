"""Lossless normalization of bare and resolver-URL DOI identifiers."""
import re
import urllib.parse


def normalize_doi(doi: str) -> str:
    doi = re.sub(r"^doi\s*:\s*", "", doi.strip(), flags=re.IGNORECASE)
    prefix = re.match(r"^(?:https?://)?(?:dx\.)?doi\.org/", doi, re.IGNORECASE)
    return urllib.parse.unquote(doi[prefix.end():]) if prefix else doi
