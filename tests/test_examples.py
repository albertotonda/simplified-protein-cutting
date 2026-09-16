"""Tests for the bundled example-data helpers (example_config, lactoferrin_protein,
pepsin_cuts). These read data/lactoferrin.json back out of the installed package (bundled as
package data via cpp/CMakeLists.txt's install(FILES ...) rule) -- these tests are as much a
check that the packaging actually shipped the file as they are a check of the functions
themselves.
"""

import pandas as pd
import pytest

import seqcleave


def test_example_config_has_the_expected_shape():
    config = seqcleave.example_config()

    assert set(config.keys()) == {"parameters", "proteins", "cuts", "alterations", "terminalAlterations"}
    assert len(config["proteins"]) == 1
    assert config["proteins"][0]["quantity"] == 500


def test_example_config_runs():
    config = seqcleave.example_config()
    config["proteins"][0]["quantity"] = 5
    config["parameters"]["maxDH"] = 0.02

    df = seqcleave.simulate(config)

    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0


def test_lactoferrin_protein_matches_the_readme():
    # the README states these facts about the bundled sample protein -- if this ever
    # changes, the README needs updating to match, not just this test
    protein = seqcleave.lactoferrin_protein()

    assert set(protein.keys()) == {"name", "quantity", "sequence", "disulfideBonds"}
    assert protein["name"] == "lactoferrin"
    assert protein["quantity"] == 500
    assert len(protein["sequence"]) == 689
    assert len(protein["disulfideBonds"]) == 36


def test_pepsin_cuts_has_exactly_the_enzyme_behavior_keys():
    cuts = seqcleave.pepsin_cuts()

    assert set(cuts.keys()) == {"cuts", "alterations", "terminalAlterations"}
    assert len(cuts["cuts"]) > 0
    assert len(cuts["alterations"]) > 0


@pytest.mark.parametrize("loader", [seqcleave.example_config, seqcleave.lactoferrin_protein, seqcleave.pepsin_cuts])
def test_each_loader_returns_an_independent_copy_every_call(loader):
    first = loader()
    # mutate deeply, not just a top-level key, to make sure this is a real deep copy
    if "parameters" in first:
        first["parameters"]["maxDH"] = -1
    elif "cuts" in first:
        first["cuts"].clear()
    else:
        first["quantity"] = -1

    second = loader()

    assert second != first, "mutating one call's result affected a later call -- not a real copy"


def test_lactoferrin_protein_and_pepsin_cuts_compose_to_the_same_result_as_example_config():
    # the whole point of splitting these into three functions instead of one: combining the
    # two granular ones by hand must be *exactly* equivalent to the convenience function,
    # since example_config() is implemented as nothing more than that combination
    config = seqcleave.example_config()
    config["parameters"]["randomSeed"] = 42
    config["proteins"][0]["quantity"] = 5
    config["parameters"]["maxDH"] = 0.05

    manual = {
        "parameters": dict(config["parameters"]),
        "proteins": [seqcleave.lactoferrin_protein()],
        **seqcleave.pepsin_cuts(),
    }
    manual["proteins"][0]["quantity"] = 5

    df_from_example_config = seqcleave.simulate(config)
    df_from_manual_composition = seqcleave.simulate(manual)

    pd.testing.assert_frame_equal(df_from_example_config, df_from_manual_composition)


def test_pepsin_cuts_can_be_combined_with_a_different_protein():
    # the actual selling point of splitting the enzyme data out on its own: it should work
    # against a protein that has nothing to do with lactoferrin
    config = {
        "parameters": {"maxDH": 0.05, "randomSeed": 7},
        "proteins": [{"name": "toy", "quantity": 5, "sequence": "aprknvrwctisqpewfkcrrwqwrmkklga"}],
        **seqcleave.pepsin_cuts(),
    }

    df = seqcleave.simulate(config)

    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
