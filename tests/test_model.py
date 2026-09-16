import os

import pytest

from seqcleave import EndoproteaseModel


def test_default_parameters_match_the_cpp_constructor():
    # these mirror EndoproteaseModel's C++ constructor defaults (cpp/EndoproteaseModel.cpp)
    # -- if this test breaks, either the binding or the constructor changed, and the other
    # side needs to catch up
    model = EndoproteaseModel()
    assert model.current_enzyme == pytest.approx(1.0)
    assert model.max_attempts == 1000
    assert model.max_attempts_per_time == 10
    assert model.max_dh == pytest.approx(100.0)
    assert model.max_time == 10000
    assert model.enzyme_always_dying is False
    assert model.enzyme_dying_ratio == pytest.approx(1.0)
    assert model.random_seed == 0
    assert model.t == 0


def test_parameters_are_settable_before_run(tiny_deterministic_config):
    model = EndoproteaseModel()
    assert model.read_config(tiny_deterministic_config) == 0

    model.max_attempts = 5
    model.random_seed = 123

    assert model.max_attempts == 5
    assert model.random_seed == 123


def test_read_json_reports_failure_for_a_missing_file():
    model = EndoproteaseModel()
    # the low-level binding mirrors the C++ return code (0 success, non-zero failure)
    # rather than raising -- that translation happens one level up, in simulate()
    assert model.read_json("this-file-does-not-exist.json") != 0


def test_run_actually_advances_the_simulation(tiny_deterministic_config):
    model = EndoproteaseModel()
    assert model.read_config(tiny_deterministic_config) == 0
    assert model.t == 0

    model.run()

    assert model.t > 0


def test_the_one_guaranteed_cut_happens(tiny_deterministic_config):
    # see tiny_deterministic_config's docstring: this config has exactly one possible cut,
    # with probability 1.0, guaranteed to happen on the very first iteration -- so time2
    # (the cut counter) must reach exactly 1, deterministically, no matter the random seed.
    # period=1 here: the default period=10 only samples every 10th cut, so a single cut
    # would never show up in the sampled series at all (see test_regression.py's note on
    # what "period" means) -- period=1 samples every cut instead.
    model = EndoproteaseModel()
    model.read_config(tiny_deterministic_config)
    model.run()

    series = model.compute_time_series(period=1)
    assert series["time2"][-1] == 1


def test_compute_time_series_shape(tiny_deterministic_config):
    model = EndoproteaseModel()
    model.read_config(tiny_deterministic_config)
    model.run()

    series = model.compute_time_series()

    assert set(["time", "time2", "enzyme"]).issubset(series.keys())
    # every column must be the same length (one entry per sampled row)
    lengths = {len(values) for values in series.values()}
    assert len(lengths) == 1


def test_write_log_produces_a_csv_with_the_expected_header(tiny_deterministic_config, tmp_path):
    model = EndoproteaseModel()
    model.read_config(tiny_deterministic_config)
    model.run()

    output_path = tmp_path / "out.csv"
    model.write_log(str(output_path))

    assert output_path.exists()
    header = output_path.read_text(encoding="utf-8").splitlines()[0]
    assert header.startswith('"time","time2","enzyme"')
