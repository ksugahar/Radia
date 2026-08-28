"""Long crash-isolated transformed-container field regression."""

import pytest

from tests.test_radiafield_transformed_container import (
    TestBatchFldTransformedContainer,
)

pytestmark = pytest.mark.slow


def test_batch_fld_no_crash_and_matches_reference():
    TestBatchFldTransformedContainer()._validate_batch_fld_no_crash_and_matches_reference()
