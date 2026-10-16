"""Isolate mocked controller/Engine contracts from fixed system scratch."""
import tempfile
import pytest


@pytest.fixture(autouse=True)
def isolated_controller_scratch(monkeypatch, tmp_path):
    temporary_directory = tempfile.TemporaryDirectory

    class IsolatedTemporaryDirectory(temporary_directory):
        def __init__(self, *args, **kwargs):
            if str(kwargs.get("dir", "")).replace("/", "\\").casefold() == r"c:\temp":
                kwargs["dir"] = tmp_path
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(tempfile, "TemporaryDirectory", IsolatedTemporaryDirectory)

    def owned_file_factory(factory):
        def create(*args, **kwargs):
            if str(kwargs.get("dir", "")).replace("/", "\\").casefold() == r"c:\temp":
                kwargs["dir"] = tmp_path
            return factory(*args, **kwargs)
        return create

    monkeypatch.setattr(tempfile, "TemporaryFile", owned_file_factory(tempfile.TemporaryFile))
    monkeypatch.setattr(tempfile, "NamedTemporaryFile", owned_file_factory(tempfile.NamedTemporaryFile))
