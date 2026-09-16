import copy
import json
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
LACTOFERRIN_JSON = REPO_ROOT / "data" / "lactoferrin.json"


def load_jsonc(path: Path) -> dict:
    """Parse a JSON-with-comments file the same way the C++ loader does (// stripped)."""
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"//[^\n]*", "", text)
    return json.loads(text)


@pytest.fixture
def lactoferrin_config() -> dict:
    """The real sample configuration, scaled down (quantity=5) and seeded for determinism.

    A fresh dict is returned on every call (via the module-level cache below plus a deep
    copy), so tests can freely mutate it without affecting other tests.
    """
    config = copy.deepcopy(_lactoferrin_config_template())
    config["proteins"][0]["quantity"] = 5
    config["parameters"]["randomSeed"] = 42
    config["parameters"]["maxDH"] = 0.05
    return config


@pytest.fixture
def tiny_deterministic_config() -> dict:
    """A minimal, fully deterministic config.

    "aaffaa" has exactly one cuttable bond (the "f"-"f" one, the only pair with a
    probability > 0), with probability 1.0 -- i.e. it always cuts, and every other bond
    never does. Combined with the default maxAttemptsPerTime=10 (more than this 6-residue
    protein's length), the single possible cut is guaranteed to happen on the very first
    simulation iteration, regardless of the random seed: nothing here is actually left to
    chance, so this config is for tests that must not be able to flake.
    """
    return {
        "parameters": {"maxDH": 1.0, "maxAttempts": 50},
        "proteins": [{"name": "toy", "quantity": 1, "sequence": "aaffaa"}],
        "cuts": {"f": {"f": 1.0}},
    }


@pytest.fixture
def tmp_json_file(tmp_path):
    """Write a given dict as a JSON file under tmp_path and return its path (as a str)."""

    def _write(config: dict, name: str = "config.json") -> str:
        path = tmp_path / name
        path.write_text(json.dumps(config), encoding="utf-8")
        return str(path)

    return _write


_lactoferrin_template_cache = None


def _lactoferrin_config_template() -> dict:
    global _lactoferrin_template_cache
    if _lactoferrin_template_cache is None:
        _lactoferrin_template_cache = load_jsonc(LACTOFERRIN_JSON)
    return _lactoferrin_template_cache
