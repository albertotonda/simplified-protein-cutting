import pandas as pd
import pytest

import endocleave


def test_simulate_with_a_dict_returns_a_dataframe(tiny_deterministic_config):
    df = endocleave.simulate(tiny_deterministic_config)

    assert isinstance(df, pd.DataFrame)
    assert list(df.columns[:3]) == ["time", "time2", "enzyme"]
    assert len(df) > 0


def test_simulate_with_a_file_path(tiny_deterministic_config, tmp_json_file):
    path = tmp_json_file(tiny_deterministic_config)
    df = endocleave.simulate(path)

    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0


def test_dict_and_equivalent_file_give_identical_results(lactoferrin_config, tmp_json_file):
    # same config, same fixed seed, two different ways in -- must match exactly
    path = tmp_json_file(lactoferrin_config)

    df_from_dict = endocleave.simulate(lactoferrin_config)
    df_from_file = endocleave.simulate(path)

    pd.testing.assert_frame_equal(df_from_dict, df_from_file)


def test_output_parameter_writes_a_csv_matching_the_returned_dataframe(tiny_deterministic_config, tmp_path):
    output_path = tmp_path / "statistics.csv"

    df = endocleave.simulate(tiny_deterministic_config, output=str(output_path))

    assert output_path.exists()
    df_from_csv = pd.read_csv(output_path)
    # check_dtype=False: "enzyme" is a C++ double, but the CSV writer (like any plain text
    # format) doesn't carry a distinct int/float type -- when every value in the column
    # happens to be a whole number (e.g. always 1.0), it's written as "1", and pandas then
    # infers int64 reading it back instead of float64. Values matching is what matters here.
    pd.testing.assert_frame_equal(df, df_from_csv, check_dtype=False)


def test_simulate_rejects_a_config_missing_required_fields():
    with pytest.raises(ValueError):
        # no "proteins" key at all
        endocleave.simulate({"parameters": {"maxDH": 0.1}, "cuts": {}})


def test_simulate_rejects_a_missing_file():
    with pytest.raises(ValueError):
        endocleave.simulate("this-file-does-not-exist.json")


def test_period_parameter_changes_row_count(lactoferrin_config):
    df_coarse = endocleave.simulate(lactoferrin_config, period=50)
    df_fine = endocleave.simulate(lactoferrin_config, period=10)

    # a larger sampling period means fewer (or equal) rows -- same simulation, coarser sampling
    assert len(df_coarse) <= len(df_fine)
