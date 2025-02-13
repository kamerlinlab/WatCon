"""
Unit and regression test for the WatCon package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import WatCon


def test_WatCon_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "WatCon" in sys.modules
