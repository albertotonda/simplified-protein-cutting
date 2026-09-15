"""endocleave: stochastic simulation of protein hydrolysis by endoproteases.

This package will wrap the C++ simulation core (see the ``cpp/`` directory at the
repository root) via pybind11 bindings. The bindings themselves are not implemented
yet -- this is currently just the package skeleton, put in place while reorganizing
the repository for packaging (see the project README's "Project status / roadmap"
section).

Once the compiled extension exists, this module will expose:
- ``EndoproteaseModel``: a near 1:1 binding of the C++ class, for full control.
- ``simulate(config, ...)``: a high-level convenience function taking either a path
  to a JSON configuration file or a plain Python dict, returning a pandas DataFrame.
"""

__version__ = "0.1.0.dev0"
