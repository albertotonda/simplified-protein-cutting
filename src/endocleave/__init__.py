"""endocleave: stochastic simulation of protein hydrolysis by endoproteases.

Simulates an endoprotease (pepsin, by default -- but the model works for any enzyme, given
per-position cleavage-frequency data) cutting one or more copies of a protein, and reports how
the resulting peptide population evolves over the course of the simulated reaction.

Two ways to use this package:

- The high-level :func:`simulate` function: point it at a JSON configuration file (or an
  equivalent Python dict) and get back a :class:`pandas.DataFrame`.
- The lower-level :class:`EndoproteaseModel` class (a near-direct binding of the C++ core) for
  full control -- read the configuration, tweak parameters, run the simulation, and pull out
  results as a CSV file or an in-memory dict.

New to this and don't have a configuration handy? :func:`example_config` returns one that's
ready to run (``endocleave.simulate(endocleave.example_config())``); :func:`lactoferrin_protein`
and :func:`pepsin_cuts` give you its two halves separately, e.g. to try pepsin's published
cleavage data against a protein of your own.

Logging goes through the standard library's ``logging`` module, under the name
``"endocleave"``. Verbosity is controlled with :func:`set_log_level` rather than
``logging.getLogger("endocleave").setLevel(...)`` directly -- see that function's docstring
for why.
"""

import copy
import functools
import json
import os
import re
from importlib import resources
from typing import Any, Dict, Mapping, Optional, Union

import pandas as pd

from ._core import EndoproteaseModel, set_log_level

try:
    from importlib.metadata import PackageNotFoundError, version as _version

    __version__ = _version("endocleave")
except PackageNotFoundError:
    # package is being imported without being installed (e.g. running straight from a
    # source checkout without "pip install -e ."): fall back to a clearly-not-real version
    __version__ = "0.0.0+unknown"

__all__ = [
    "EndoproteaseModel",
    "set_log_level",
    "simulate",
    "example_config",
    "lactoferrin_protein",
    "pepsin_cuts",
    "__version__",
]


@functools.lru_cache(maxsize=1)
def _load_example_config() -> Dict[str, Any]:
    """Parse the bundled sample configuration once, cached for the life of the process.

    This is the same file as the repository's data/lactoferrin.json (bundled as package data
    at build time, see cpp/CMakeLists.txt) -- there's exactly one copy of this data, this just
    reads it back out of wherever it landed inside the installed package.
    """
    text = resources.files(__package__).joinpath("data", "lactoferrin.json").read_text(encoding="utf-8")
    # strip "//" comments the same way the C++ loader does, before handing to json.loads
    return json.loads(re.sub(r"//[^\n]*", "", text))


def simulate(
    config: Union[str, "os.PathLike[str]", Mapping],
    output: Optional[Union[str, "os.PathLike[str]"]] = None,
    period: int = 10,
) -> pd.DataFrame:
    """Run a full simulation and return the results as a DataFrame.

    Parameters
    ----------
    config:
        Either a path to a JSON configuration file, or a dict with the same schema (see the
        project README's "Configuration format" section for the schema itself).
    output:
        If given, also write the full simulation history to this path as a CSV file (the
        same format the CLI's ``--output`` produces). Optional -- by default, nothing is
        written to disk.
    period:
        Sample the simulation every ``period`` cuts (rather than every single simulation
        iteration) when building the returned DataFrame. Matches the CLI's own default.

    Returns
    -------
    A DataFrame with columns ``time`` (iteration count), ``time2`` (cuts so far), ``enzyme``
    (current enzyme quantity), and one column per distinct peptide produced during the
    simulation, in alphabetical order -- one row per sampled point in time.
    """
    model = EndoproteaseModel()

    if isinstance(config, Mapping):
        status = model.read_config(dict(config))
    else:
        status = model.read_json(os.fspath(config))

    if status != 0:
        raise ValueError(
            "Failed to load the simulation configuration -- see the log "
            '(set_log_level("debug") for detail) for the specific error.'
        )

    model.run()

    if output is not None:
        model.write_log(os.fspath(output))

    return pd.DataFrame(model.compute_time_series(period))


def lactoferrin_protein() -> Dict[str, Any]:
    """Return the bovine lactoferrin protein entry from the original 2017 paper's case study.

    A plain dict matching one entry of the configuration schema's ``"proteins"`` list (keys:
    ``name``, ``quantity``, ``sequence``, ``disulfideBonds``) -- pass it inside a ``"proteins"``
    list of your own configuration, e.g. to try a *different* enzyme's cleavage data against
    this same protein. Combined with :func:`pepsin_cuts`, this is exactly what
    :func:`example_config` returns.

    A fresh dict is returned on every call, safe to mutate without affecting later calls.
    """
    return copy.deepcopy(_load_example_config()["proteins"][0])


def pepsin_cuts() -> Dict[str, Any]:
    """Return pepsin's published cleavage-probability data, on its own.

    From Hamuro et al., 2008 and Powers et al., 1977 (see the README's citation section) -- a
    dict with ``"cuts"``, ``"alterations"``, and ``"terminalAlterations"`` keys, ready to merge
    into a configuration for *any* protein, e.g.::

        config = {
            "parameters": {"maxDH": 0.1},
            "proteins": [{"sequence": "your own sequence here", "quantity": 100}],
            **endocleave.pepsin_cuts(),
        }

    A fresh dict is returned on every call, safe to mutate without affecting later calls.
    """
    config = _load_example_config()
    return {
        "cuts": copy.deepcopy(config["cuts"]),
        "alterations": copy.deepcopy(config["alterations"]),
        "terminalAlterations": copy.deepcopy(config["terminalAlterations"]),
    }


def example_config() -> Dict[str, Any]:
    """Return a complete, ready-to-run configuration: bovine lactoferrin cut by pepsin, the
    case study from the original paper (500 copies of the protein; see the README's
    "Configuration format" section for what each parameter means).

    Equivalent to combining :func:`lactoferrin_protein` and :func:`pepsin_cuts` with the same
    ``parameters`` used in the bundled sample file -- if you want to change quantity, maxDH,
    or the random seed, get this dict and edit it directly rather than building one by hand::

        config = endocleave.example_config()
        config["proteins"][0]["quantity"] = 5
        config["parameters"]["maxDH"] = 0.02
        df = endocleave.simulate(config)

    A fresh dict is returned on every call, safe to mutate without affecting later calls.
    """
    return copy.deepcopy(_load_example_config())
