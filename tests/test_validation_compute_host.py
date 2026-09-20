"""Compute-host marker recognition must include numbered mdx/hibino nodes."""

import platform

import pytest

from validation_test import conftest as validation_conftest


@pytest.mark.parametrize(
    "hostname",
    ["mdx", "mdx1", "mdx2", "mdx-worker", "hibino", "hibino2", "hibino-worker"],
)
def test_numbered_compute_hosts_are_recognized(monkeypatch, hostname):
    monkeypatch.setattr(platform, "node", lambda: hostname)
    assert validation_conftest._is_compute_host()


@pytest.mark.parametrize("hostname", ["LAB", "100", "developer-workstation"])
def test_non_compute_hosts_are_rejected(monkeypatch, hostname):
    monkeypatch.setattr(platform, "node", lambda: hostname)
    assert not validation_conftest._is_compute_host()
