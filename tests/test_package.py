import re

import endocleave


def test_public_api_is_present():
    assert hasattr(endocleave, "EndoproteaseModel")
    assert hasattr(endocleave, "simulate")
    assert hasattr(endocleave, "set_log_level")
    assert callable(endocleave.simulate)
    assert callable(endocleave.set_log_level)


def test_version_is_a_real_looking_version_string():
    # PEP 440-ish: starts with digits, e.g. "0.1.0.dev0" -- just a sanity check that the
    # importlib.metadata wiring in __init__.py actually resolved to something real, rather
    # than silently falling back to "0.0.0+unknown" (which would mean the installed
    # package's metadata isn't being found)
    assert re.match(r"^\d+\.\d+\.\d+", endocleave.__version__)
    assert endocleave.__version__ != "0.0.0+unknown"
