# Simplified in-silico protein cutting

[![CI](https://github.com/albertotonda/simplified-protein-cutting/actions/workflows/ci.yml/badge.svg)](https://github.com/albertotonda/simplified-protein-cutting/actions/workflows/ci.yml)

Modeling protein hydrolysis and release of peptides by endoproteases requires complex simulations, typically taking into account the 3D structure of both the enzymes and the target protein. Such structures can sometimes be difficult to predict starting from the protein's acido-aminic sequence.

This repository contains the code for an alternative approach, published in [Tonda et al. (2017), _In silico modeling of protein hydrolysis by endoproteases: a case study on pepsin digestion of bovine lactoferrin_, Food & Function, Vol. 8, Issue 12, DOI: 10.1039/C7FO00830A](https://pubs.rsc.org/fo/article-abstract/8/12/4404/566726/In-silico-modeling-of-protein-hydrolysis-by). The idea is to just consider the linear sequence of amino-acids, and then simulate the behavior of an enzyme starting from the frequency of cuts observed during previous experiments. The final peptides obtained by the simulation are qualitatively coherent with real-world experiments, even though the exact absolute quantities might be different.

If you use this software in your publications, please cite the paper (see [Citation](#citation) below).

## What is this software?

In a nutshell, this software simulates the action of an enzyme on several copies of a protein (protein structure given in input). The model is not specific to any one enzyme: it works for any endoprotease, as long as cleavage frequency data is available for it. The enzyme's behavior is considered stochastic, and during the simulation it will cut bonds with a certain probability, depending on the amino-acids to the left and right of a bond (positions P4, P3, P2, P1, P1', P2', P3', P4'). The sample configuration shipped with this repository models pepsin, using probabilities computed from analyses performed by [Hamuro et al., 2008](https://pubmed.ncbi.nlm.nih.gov/18327892/) and [Powers et al., 1977](https://link.springer.com/chapter/10.1007/978-1-4757-0719-9_9).

## Repository structure

```
data/                sample input (lactoferrin.json), a small test script, and a
                      scaled-down variant used for quick smoke tests / CI
scripts/              utility scripts (e.g. the old XML -> JSON converter)
cpp/                  C++ source code: the core, the CLI, and the pybind11 bindings
cpp/thirdparty/       vendored dependencies (nlohmann/json, spdlog, pybind11_json)
src/endocleave/       the Python package (pure-Python wrapper; the compiled
                      extension lands here too once built, see "Python package" below)
tests/                pytest suite for the Python package (see "Running the tests" below)
pyproject.toml        Python packaging config (scikit-build-core)
```

The original code is in C++, and is contained in the `cpp/` subfolder. A Python package (`endocleave`) wrapping it via pybind11 lives in `src/endocleave/` (a "src-layout" Python package, following the convention expected by Python's packaging tools — not to be confused with `cpp/`, which holds the C++ sources).

## Building the C++ code

You will need [CMake](https://cmake.org/) (3.15+) and a C++17 compiler. The code has no external dependencies to install: [nlohmann/json](https://github.com/nlohmann/json) (JSON parsing) and [spdlog](https://github.com/gabime/spdlog) (logging) are vendored, header-only, directly in `cpp/thirdparty/`, so no network access or package manager is required at build time.

```sh
cd cpp
mkdir build && cd build
cmake ..
cmake --build .
```

This produces an executable called `protein-cutting` (`protein-cutting.exe` on Windows) in the `build` directory. Any generator CMake supports should work (Unix Makefiles, Ninja, Visual Studio, ...); on Windows with MinGW/Ninja, for example: `cmake -G Ninja ..`.

The code should be cross-compiling on any platform with ISO C++ and CMake support (tested on Ubuntu 14.04/16.04 originally, and on Windows with MinGW-w64).

## Running a simulation

To run a simulation, you need a JSON file describing the protein(s) and the cut probabilities (see [Configuration format](#configuration-format-json) below). A sample file, `data/lactoferrin.json`, is provided: it contains the structure of bovine lactoferrin and simulates cutting 500 copies of the protein, using probabilities taken from Hamuro et al., 2008 and Powers et al., 1977 (see above for the DOIs).

```sh
./protein-cutting --input lactoferrin.json
```

**This can take a while** (minutes to tens of minutes, depending on your hardware and on the `quantity` of proteins simulated), and writing the final result to disk can itself take a long time — don't quit at that point.

| Argument | Description |
|---|---|
| `--input <file.json>` | **Required.** The JSON configuration file describing the simulation. |
| `--output <file.csv>` | Where to write the results. Default: `statistics.csv`. |
| `--verbose` | Prints debug-level detail to the console (shorthand for `--log-level debug`). |
| `--log-level <level>` | Sets the console log level explicitly: `trace`, `debug`, `info`, `warn`, `error`, or `off`. Default: `info`. |
| `--log-file <file.log>` | Additionally writes a full trace-level debug log to the given file. **No file is ever written unless this is explicitly given.** |

By default the program only prints progress/warning/error messages to the console (`info` level) — no per-position debug trace, which used to flood the console/log files in earlier versions of this code.

## Configuration format (JSON)

Simulations are configured entirely from a JSON file — no need to modify the source code to change the protein(s), probabilities, or simulation parameters. Comments (`//` and `/* */`) are supported by the loader and stripped before parsing, so configuration files can be annotated just like code; `data/lactoferrin.json` is heavily commented and is the best starting point for writing your own.

The file has four top-level sections:

- **`parameters`**: simulation-wide settings — `randomSeed` (`null` for a time-based seed), `maxTime` (max iterations), `maxDH` (stop once this degree of hydrolysis is reached), `maxAttemptsPerTime`, `maxAttempts` (stop after this many consecutive failed cut attempts), and the experimental `initialEnzyme` / `enzymeAlwaysDying` / `enzymeDyingRatio` (enzyme activity decaying over time).
- **`proteins`**: an array of proteins to simulate, each with a `name`, a `quantity` (number of copies), a `sequence` (the amino-acid chain), and `disulfideBonds` (1-indexed positions the enzyme finds harder to cut — not all of them are strictly disulfide bonds, some are glycosylations).
- **`cuts`**: base probability of cutting a bond, keyed by the amino-acid to the left (P1) and right (P1') of the bond, e.g. `"cuts": { "f": { "y": 0.65, "f": 0.85, ... }, ... }`. Bonds not listed default to probability 0.
- **`alterations`** / **`terminalAlterations`**: position-dependent adjustments to the base probability for amino-acids found further away from the bond (P2-P4 / P2'-P4'), and multipliers applied near either end of a peptide chain.

Until 2026, configuration files were XML, parsed with the [tinyxml](http://www.grinninglizard.com/tinyxml/) library; the format was switched to JSON (parsed with nlohmann/json) for easier editing and future Python bindings. The original sample file, `data/lactoferrin.xml`, is kept in the repository as a historical reference — it is no longer read by the code. The one-off script used for the conversion, `scripts/xml_to_json.py`, is kept for reference in case other old XML configuration files need migrating.

## Output format

The program produces a CSV file (`statistics.csv` by default) tracking the quantity of each peptide over the course of the simulation. Columns are `time` (iteration count), `time2` (number of cuts so far), `enzyme` (currently always 1.0, reserved for future developments), followed by one column per distinct peptide produced during the simulation, in alphabetical order. Each row gives the count of each peptide at that point in the simulation. The resulting file is usually large (~70 MB for the lactoferrin example); extracting meaningful information from it typically requires a separate analysis script rather than manual inspection.

## Python package

The C++ core is also available from Python, as the `endocleave` package: [pybind11](https://github.com/pybind/pybind11) bindings over the same core used by the CLI, built with [scikit-build-core](https://github.com/scikit-build/scikit-build-core) so a normal `pip install` compiles everything automatically — no separate C++ build step, and no dependency on CMake or a compiler once installed.

```sh
pip install .
```

(run from the repository root; not yet published on PyPI).

```python
import endocleave

# high-level: JSON file or dict in, a pandas.DataFrame out (same columns as the CSV above)
df = endocleave.simulate("data/lactoferrin.json")

# optionally, also write the CSV file, same as the CLI's --output
df = endocleave.simulate("data/lactoferrin.json", output="statistics.csv")
```

Don't have a configuration handy? The bovine lactoferrin / pepsin case study from the paper ships with the package, in three forms — a ready-to-run combination, and its two halves separately (useful, for example, to try pepsin's published cleavage data against a protein of your own):

```python
df = endocleave.simulate(endocleave.example_config())   # ready to run as-is

endocleave.lactoferrin_protein()   # -> just the "proteins" entry (sequence, disulfideBonds, quantity)
endocleave.pepsin_cuts()           # -> just the enzyme data ("cuts", "alterations", "terminalAlterations")

# e.g. pepsin's cleavage data against a different protein:
config = {
    "parameters": {"maxDH": 0.1},
    "proteins": [{"sequence": "your own sequence here", "quantity": 100}],
    **endocleave.pepsin_cuts(),
}
```

For full control, use the lower-level `EndoproteaseModel` class directly — a near 1-to-1 binding of the C++ class, with plain read/write attributes for every simulation parameter:

```python
from endocleave import EndoproteaseModel
import pandas as pd

model = EndoproteaseModel()
model.read_json("data/lactoferrin.json")   # or model.read_config({...}) for a native dict
model.max_dh = 0.05
model.run()
df = pd.DataFrame(model.compute_time_series())
```

Logging goes through the standard `logging` module, under the name `"endocleave"`. Verbosity is controlled with `endocleave.set_log_level(...)` rather than `logging.getLogger("endocleave").setLevel(...)` directly: the native core uses its own log level as a performance gate (deciding whether to even format a message), so the two have to stay in sync, and `set_log_level()` does that in one call.

```python
endocleave.set_log_level("debug")
```

### Running the tests

```sh
pip install -e ".[test]"
pytest
```

The suite (`tests/`) covers the Python bindings and the `simulate()`/`EndoproteaseModel` API: a fixed-seed regression check (meaningful and portable across platforms, since the random engine is `std::mt19937`, a standardized algorithm), a mass-balance invariant (every cut turns one peptide into two, so the total peptide count must always equal the original quantity plus the number of cuts so far — true for any config or seed), and error handling for malformed input.

## Project status / roadmap

- ✅ Configuration format switched from XML to JSON.
- ✅ Logging rewritten (leveled, quiet by default, opt-in file output) in preparation for reuse from other languages.
- ✅ Model and parameter names generalized (`EndoproteaseModel`, `enzyme*` fields) — the simulation was never pepsin-specific, and now neither is its naming.
- ✅ Repository reorganized (`cpp/` for the C++ core, `src/endocleave/` for the Python package).
- ✅ pybind11 bindings, a `simulate()` convenience API, and a working `pip install .` (via scikit-build-core).
- ✅ A `pytest` suite (`tests/`) covering the Python API, a fixed-seed regression check, and a seed-independent structural invariant.
- ✅ CI (`.github/workflows/ci.yml`): builds the CLI and runs the pytest suite on Linux, macOS, and Windows on every push/PR.
- ⏳ Planned: publish `endocleave` on PyPI, with prebuilt wheels (via `cibuildwheel`) for the common platforms.

## Citation

If you use this software in your publications, please cite:

> Tonda, Alberto and Grosvenor, Anita J and Clerens, Stefan and Le Feunteun, Steven,
> "In silico modeling of protein hydrolysis by endoproteases: a case study on pepsin digestion of bovine lactoferrin",
> Food & Function, 2017, DOI: 10.1039/C7FO00830A

<details open>
<summary>BibTeX</summary>

```bibtex
@article{tonda2017insilico,
    author = {Tonda, Alberto and Grosvenor, Anita and Clerens, Stefan and Le Feunteun, Steven},
    title = {In silico modeling of protein hydrolysis by endoproteases: a case study on pepsin digestion of bovine lactoferrin},
    journal = {Food \& Function},
    volume = {8},
    number = {12},
    pages = {4404-4413},
    year = {2017},
    month = {12},
    abstract = {This paper presents a novel model of protein hydrolysis and release of peptides by endoproteases. It requires the amino-acid sequence of the protein substrate to run, and makes use of simple Monte-Carlo in silico simulations to qualitatively and quantitatively predict the peptides that are likely to be produced during the course of the proteolytic reaction. In the present study, the model is applied to the case of pepsin, the gastric protease. Unlike pancreatic proteases, pepsin has a low substrate specificity and therefore displays a stochastic behavior that is particularly challenging to model and predict. Two versions of the model are studied and compared with peptidomic data obtained during pepsin hydrolysis of bovine lactoferrin. The first version of the model takes into account cleavage probabilities according to the amino acids in position P1–P1′ only, whereas the second version also accounts for the influence of neighbor amino acids (P4, P3, P2, P2′, P3′, P4′) and peptide terminal ends. The second version of the model was able to reproduce many real-world features of the reported behavior of pepsin, such as the peptide size distribution, or the quantity of free amino-acids. More remarkably, 50\% of the experimentally monitored peptides (44/87) lay within the 120 most abundant simulated peptides. The presented methodology has the advantage of being applicable not only to different proteins, but to different enzymes as well, as long as cleavage frequency data are available.},
    issn = {2042-6496},
    doi = {10.1039/c7fo00830a},
    url = {https://doi.org/10.1039/c7fo00830a},
    eprint = {https://pubs.rsc.org/fo/article-pdf/8/12/4404/5512608/c7fo00830a.pdf},
}
```

</details>

## License

Copyright (c) 2017, Alberto Tonda \<alberto.tonda@gmail.com\>

Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted, provided that the above copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

`nlohmann/json` (vendored in `cpp/thirdparty/nlohmann/`) is authored by Niels Lohmann and distributed under the MIT license. `spdlog` (vendored in `cpp/thirdparty/spdlog/`) is authored by Gabi Melman and distributed under the MIT license. `pybind11_json` (vendored in `cpp/thirdparty/pybind11_json/`) is authored by Martin Renou and distributed under the BSD 3-Clause license. [pybind11](https://github.com/pybind/pybind11) itself and [scikit-build-core](https://github.com/scikit-build/scikit-build-core) are build-time-only dependencies (not vendored, resolved automatically by `pip` from `pyproject.toml`), both distributed under permissive licenses (BSD-style and Apache 2.0, respectively). The original tinyxml library (no longer used, kept out of the repository) was authored by Lee Thomason, Yves Berquin, and Andrew Ellerton.

In case you need help, advice, or you notice a bug, please contact Alberto Tonda \<alberto.tonda@gmail.com\>.
