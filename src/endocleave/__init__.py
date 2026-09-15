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

Logging goes through the standard library's ``logging`` module, under the name
``"endocleave"``. Verbosity is controlled with :func:`set_log_level` rather than
``logging.getLogger("endocleave").setLevel(...)`` directly -- see that function's docstring
for why.
"""

import os
from typing import Mapping, Optional, Union

import pandas as pd

from ._core import EndoproteaseModel, set_log_level

try:
    from importlib.metadata import PackageNotFoundError, version as _version

    __version__ = _version("endocleave")
except PackageNotFoundError:
    # package is being imported without being installed (e.g. running straight from a
    # source checkout without "pip install -e ."): fall back to a clearly-not-real version
    __version__ = "0.0.0+unknown"

__all__ = ["EndoproteaseModel", "set_log_level", "simulate", "__version__"]


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
