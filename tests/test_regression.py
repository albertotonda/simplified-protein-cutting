"""Regression tests.

These formalize the manual "run it, diff the CSV" verification that was done by hand,
repeatedly, throughout this project's C++ -> JSON -> logging -> RNG -> Python-bindings
history. A fixed seed is meaningful here specifically because the random engine is
std::mt19937 (a standardized, fully portable algorithm) -- unlike the old libc rand() this
project used to use, which was not standardized and would have made a byte-for-byte golden
test like this fragile across platforms/compilers.
"""

import pandas as pd

import endocleave


def test_mass_balance_invariant_holds_at_every_sampled_row(lactoferrin_config):
    # a cut always takes one peptide and produces two: net +1 peptide per cut, regardless
    # of *which* peptide gets cut or into what. So at any point in time, the total count of
    # all peptides must equal the original protein quantity plus the number of cuts so far
    # -- true for any config, any seed, independent of the exact stochastic path taken.
    df = endocleave.simulate(lactoferrin_config)

    original_quantity = lactoferrin_config["proteins"][0]["quantity"]
    peptide_columns = df.columns[3:]

    total_peptides = df[peptide_columns].sum(axis=1)
    expected = original_quantity + df["time2"]

    pd.testing.assert_series_equal(total_peptides, expected, check_names=False)


def test_time_and_time2_are_non_decreasing(lactoferrin_config):
    df = endocleave.simulate(lactoferrin_config)

    assert df["time"].is_monotonic_increasing
    assert df["time2"].is_monotonic_increasing


def test_fixed_seed_reproduces_known_values(lactoferrin_config):
    # captured from an actual run of this exact fixture (quantity=5, randomSeed=42,
    # maxDH=0.05) -- if this ever changes, either something genuinely regressed, or the
    # change was intentional and these constants need updating to match.
    #
    # NOTE: this is only meaningful *because* EndoproteaseModel::randomUnit()/randomIndex()
    # deliberately avoid std::uniform_real_distribution/uniform_int_distribution (see their
    # comments in EndoproteaseModel.cpp) -- those are built on a portable engine
    # (std::mt19937) but are themselves only implementation-defined, so libstdc++ (Linux,
    # MinGW) and libc++ (macOS/Clang) silently produced different results from the exact
    # same seed. Found via this very test failing in CI on macOS but not Linux/Windows.
    df = endocleave.simulate(lactoferrin_config)

    assert df.shape == (18, 327)

    last_row = df.iloc[-1]
    assert last_row["time"] == 1927
    assert last_row["time2"] == 170
    assert last_row["enzyme"] == 1.0
