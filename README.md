# Simplified in-silico protein cutting

[![CI](https://github.com/albertotonda/simplified-protein-cutting/actions/workflows/ci.yml/badge.svg)](https://github.com/albertotonda/simplified-protein-cutting/actions/workflows/ci.yml)

Modeling protein hydrolysis and release of peptides by endoproteases requires complex, compute-intensive simulations, typically taking into account the 3D structure of both the enzymes and the target protein. Furthermore, despite recent advancements in folding, such structures can still be difficult to predict starting from the protein's acido-aminic sequence.

This repository contains the code for an alternative approach, originally published in [Tonda et al. (2017), _In silico modeling of protein hydrolysis by endoproteases: a case study on pepsin digestion of bovine lactoferrin_, Food & Function, Vol. 8, Issue 12, DOI: 10.1039/C7FO00830A](https://pubs.rsc.org/fo/article-abstract/8/12/4404/566726/In-silico-modeling-of-protein-hydrolysis-by). The idea is to just consider the linear sequence of amino-acids, and then simulate the behavior of an enzyme starting from the frequency of cuts observed during previous experiments, plus other information, for example the presence of disuflide bonds. The final peptides obtained by the simulation aim to be qualitatively coherent with real-world experiments, even though the exact absolute quantities might be different.

If you use this software in your publications, please cite the original paper (see [Citing this software](#citing-this-software) below).

## What is this software?

In a nutshell, this software simulates the action of an enzyme on several copies of a protein (protein structure given in input). The model is not specific to any one enzyme: it works for any endoprotease, as long as cleavage frequency data is available for it. The enzyme's behavior is assumed to be **stochastic**, and during the simulation it will cut bonds with a certain probability, depending on the four amino-acids to the left (positions traditionally labeled as P4, P3, P2, P1) and right (P1', P2', P3', P4') of a bond. Here below is the flowchart of the algorithm.

<p align="center">
  <a href="https://raw.githubusercontent.com/albertotonda/simplified-protein-cutting/main/figures/flowchart.jpg">
    <img src="https://raw.githubusercontent.com/albertotonda/simplified-protein-cutting/main/figures/flowchart.jpg" alt="Flow chart of the simulation" width="360">
  </a>
</p>

The sample configuration shipped with this package models **pepsin**, using probabilities computed from analyses performed by [Hamuro et al., 2008](https://pubmed.ncbi.nlm.nih.gov/18327892/) and [Powers et al., 1977](https://link.springer.com/chapter/10.1007/978-1-4757-0719-9_9), which was the case study for the original [Tonda et al., 2017](https://pubs.rsc.org/fo/article-abstract/8/12/4404/566726/In-silico-modeling-of-protein-hydrolysis-by) paper.

Most people will want the **Python package**, described next. Under the Python hood, there is a core written in C++, for speed. The C++ core, CLI, and file formats behind it are documented further down, in [C++ core](#c-core).

## Python package

You can install this software as the `seqcleave` package, with:

```sh
pip install seqcleave
```

To install from a checkout of this repository instead (e.g. to test local changes), run `pip install .` from the repository root; that will also compile the extension automatically, via `scikit-build-core`.

You can check whether everything works fine by running the case study, loading the configuration for **pepsin** enzyme cleaving the protein **bovine lactoferrin**:

```python
import seqcleave

# load complete configuration as a dictionary
config = seqcleave.example_config()
# set the random seed for reproducibility
config["randomSeed"] = 42
# reduce initial quantity of the first protein (pepsin) to speed up the computation
config["proteins"][0]["quantity"] = 10

# launch the simulation (it might take a few seconds)
df = seqcleave.simulate(config)
# the DataFrame contains the quantity of different peptides at each instant of time 
print(df)
# if an output filename is specified, the result is also saved as a CSV file
#df = seqcleave.simulate(config, output="peptides.csv")
```
The output should be:
```
    time  time2  enzyme  a  aca  acaf  acafltr  ad  aedvgdva  ...  
0      0      0     1.0  0    0     0        0   0         0  ...     
1    114     10     1.0  0    0     0        0   0         0  ...      
2    165     20     1.0  0    0     0        0   ..   ...    ...                            ...                 ...                    ...         
68  9951    680     1.0  5    4     2        0   2         1  ...      

[69 rows x 1015 columns]
```

`seqcleave.simulate()` also accepts a JSON file in input. The configuration includes the protein sequence(s), the cutting probability for each bond, and several other parameters used in the simulation:

```python
print(seqcleave.example_config())
```
Output:
```python
{
 'parameters': {'randomSeed': None,
                'maxTime': 1200000,
                'maxDH': 0.1,
                'maxAttemptsPerTime': 1,
                'maxAttempts': 1000,
                'initialEnzyme': 1.0,
                'enzymeAlwaysDying': False,
                'enzymeDyingRatio': 1.0},
 'proteins': [{'name': 'lactoferrin',
               'quantity': 500,
               'disulfideBonds': [9,
                                  19,
                                  ...
                                  545],
               'sequence': 'aprknvrwctisqpewfkcrrwqwrmkklga...'}],
 'cuts': {'f': {'y': 0.65,
                'f': 0.85,
                    ...
                'g': 0.28},
          'l': {'y': 0.68,
                'f': 0.84,
                    ...
                'g': 0.07},
          ...
          'k': {'i': 0.02, 'g': 0.02},
          'h': {'p': 0.05}},
 'alterations': {
    'f': {'left': {'2': 0.068, '3': 0.065, '4': 0.178},
          'right': {'2': 0.076, '3': 0.138, '4': 0.138}},
                ...
    'p': {'left': {'2': 0.002, '3': 0.137, '4': 0.194},
          'right': {'2': 0.013, '3': 0.027, '4': 0.207}}},
 'terminalAlterations': {'left': {'2': 0.3333, '3': 0.3333},
                         'right': {'2': 0.3333, '3': 0.3333}}
}
```

The information related to the enzyme and the protein are also accessible separately:
```python
seqcleave.lactoferrin_protein()   # -> just the "proteins" entry (sequence, disulfideBonds, quantity)
seqcleave.pepsin_cuts()           # -> just the enzyme data ("cuts", "alterations", "terminalAlterations")

# e.g. pepsin's cleavage data against a different protein:
config = {
    "parameters": {"maxDH": 0.1},
    "proteins": [{"sequence": "your own sequence here", "quantity": 100}],
    **seqcleave.pepsin_cuts(),
}
```

See [Configuration format (dictionary and JSON)](#configuration-format-dictionary-and-json) below for the full schema `simulate()` accepts (as a Python dictionary or a JSON file).

For full control, it is also possible to use the lower-level `EndoproteaseModel` class directly — a near 1-to-1 binding of the C++ class, with plain read/write attributes for every simulation parameter:

```python
from seqcleave import EndoproteaseModel
import pandas as pd

model = EndoproteaseModel()
model.read_json("data/lactoferrin.json")   # or model.read_config({...}) for a native dict
model.max_dh = 0.05
model.run()
df = pd.DataFrame(model.compute_time_series())
```

The simulation returns a `pandas.DataFrame` (named `df` in the examples above), tracking the quantity of each type of peptide generated by the simulation over time. See [Output format](#output-format) below for more detailed information.

## Configuration format (dictionary and JSON)

Simulations — whether run through the Python package or the CLI — are configured entirely from a JSON file or Python dictionary, no need to modify the source code to change the protein(s), probabilities, or simulation parameters. Comments (`//` and `/* */`) are supported by the loader and stripped before parsing, so configuration files can be annotated just like code; [`lactoferrin.json` in the original repository](https://github.com/albertotonda/simplified-protein-cutting/blob/main/data/lactoferrin.json) is heavily commented and is the best starting point for writing your own JSON.

The file has four top-level sections:

- **`parameters`**: simulation-wide settings   
  - `randomSeed` is the seed for the random number generation; fixing it to an integer value (e.g. `42`) will ensure repeatable result; `null` is the option for a time-based seed
  - `maxTime` is the maximum number of cutting iterations (approximately corresponding to instants of time), once reached the simulation will stop
  - `maxDH` is the maximum degree of hydrolysis, once reached the simulation will stop
  - `maxAttemptsPerTime` controls the amount of cutting attempts that the enzyme will perform per iteration before going to the next iteration; sometimes, especially when only small peptides with strong bonds are available, the enzyme will not manage to cut
  - `maxAttempts` once this number of failed attempts is reached, the simulation will stop; the enzyme did not find anything to cut
  - `initialEnzyme` / `enzymeAlwaysDying` / `enzymeDyingRatio` are all experimental parameters to regulate enzyme activity decaying over time; we recommend leaving them at the default values (`initialEnzyme : 1.0, enzymeAlwaysDying : False, enzymeDyingRatio : 1.0`)
- **`proteins`**: an array of proteins to simulate, each with a `name`, a `quantity` (number of copies), a `sequence` (the amino-acid chain), and `disulfideBonds` (1-indexed positions the enzyme finds harder to cut — despite the name, not all of them are strictly disulfide bonds, some can be glycosylations).
- **`cuts`**: base probability of cutting a bond, keyed by the amino-acid to the left (P1) and right (P1') of the bond, e.g. `"cuts": { "f": { "y": 0.65, "f": 0.85, ... }, ... }` means that the probability of the enzyme cutting a `f-y` bond is `0.65` (65%), while a `f-f` bond will be cut with a `0.85` (85%) probability. Bonds not listed default to probability `0.0` and therefore will never be cut.
- **`alterations`** / **`terminalAlterations`**: position-dependent adjustments to the base probability for amino-acids found further away from the bond (P2-P4 / P2'-P4'), and multipliers applied near either end of a peptide chain (typically more difficult to cut for an enzyme).

In the original code, configuration files were XML; the format was switched to JSON (parsed with nlohmann/json) for easier editing. If you are curious, [the original sample file, `data/lactoferrin.xml`](https://github.com/albertotonda/simplified-protein-cutting/blob/main/data/lactoferrin.xml), is kept in the repository as a reference.

## Output format

The result of a `seqcleave.simulate()` is a `pandas.DataFrame` table (and possibly a CSV file), tracking the quantity of each peptide over the course of the simulation. Columns are `time` (iteration count), `time2` (number of cuts so far), `enzyme` (currently always 1.0, reserved for future developments), followed by one column per distinct peptide produced during the simulation, in alphabetical order. Each row gives the count of each peptide at that point in the simulation. The full CSV is usually large (~70 MB for the default lactoferrin example). We recommend using a separate analysis script rather than manual inspection, to extract meaningful information from it.

## C++ core

The Python package wraps a C++ core, which can also be built and run standalone as a CLI. This section is for building from source or working on the C++ code directly — most users won't need it.

### Repository structure

```
data/                sample input (lactoferrin.json), a small test script, and a
                      scaled-down variant used for quick smoke tests / CI
scripts/              utility scripts (e.g. the old XML -> JSON converter)
cpp/                  C++ source code: the core, the CLI, and the pybind11 bindings
cpp/thirdparty/       vendored dependencies (nlohmann/json, spdlog, pybind11_json)
src/seqcleave/       the Python package (pure-Python wrapper; the compiled
                      extension lands here too once built, see "Python package" above)
tests/                pytest suite for the Python package (see "Running the tests" above)
pyproject.toml        Python packaging config (scikit-build-core)
```

The original code is in C++, and is contained in the `cpp/` subfolder. A Python package (`seqcleave`) wrapping it via pybind11 lives in `src/seqcleave/` (a "src-layout" Python package, following the convention expected by Python's packaging tools — not to be confused with `cpp/`, which holds the C++ sources).

### Building the C++ code

You will need [CMake](https://cmake.org/) (3.15+) and a C++17 compiler. The code has no external dependencies to install: [nlohmann/json](https://github.com/nlohmann/json) (JSON parsing) and [spdlog](https://github.com/gabime/spdlog) (logging) are vendored, header-only, directly in `cpp/thirdparty/`, so no network access or package manager is required at build time.

```sh
cd cpp
mkdir build && cd build
cmake ..
cmake --build .
```

This produces an executable called `protein-cutting` (`protein-cutting.exe` on Windows) in the `build` directory. Any generator CMake supports should work (Unix Makefiles, Ninja, Visual Studio, ...); on Windows with MinGW/Ninja, for example: `cmake -G Ninja ..`.

The code should be cross-compiling on any platform with ISO C++ and CMake support (tested on Ubuntu 14.04/16.04 originally, and on Windows with MinGW-w64).

### Running a simulation (CLI)

To run a simulation, you need a JSON file describing the protein(s) and the cut probabilities (see [Configuration format](#configuration-format-json) above). A sample file, `data/lactoferrin.json`, is provided: it contains the structure of bovine lactoferrin and simulates cutting 500 copies of the protein, using probabilities taken from Hamuro et al., 2008 and Powers et al., 1977 (see above for the DOIs).

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

## Project status / roadmap

- ✅ Configuration format switched from XML to JSON.
- ✅ Logging rewritten (leveled, quiet by default, opt-in file output) in preparation for reuse from other languages.
- ✅ Model and parameter names generalized (`EndoproteaseModel`, `enzyme*` fields) — the simulation was never pepsin-specific, and now neither is its naming.
- ✅ Repository reorganized (`cpp/` for the C++ core, `src/seqcleave/` for the Python package).
- ✅ pybind11 bindings, a `simulate()` convenience API, and a working `pip install .` (via scikit-build-core).
- ✅ A `pytest` suite (`tests/`) covering the Python API, a fixed-seed regression check, and a seed-independent structural invariant.
- ✅ CI (`.github/workflows/ci.yml`): builds the CLI and runs the pytest suite on Linux, macOS, and Windows on every push/PR.
- ✅ `seqcleave` [published on PyPI](https://pypi.org/project/seqcleave/), with prebuilt wheels (via `cibuildwheel`) for Linux, Windows, and macOS (Intel + Apple Silicon).

## Citing this software

If you use this software in your publications, please cite:

> Tonda, Alberto and Grosvenor, Anita J. and Clerens, Stefan and Le Feunteun, Steven,
> "In silico modeling of protein hydrolysis by endoproteases: a case study on pepsin digestion of bovine lactoferrin",
> Food & Function, 2017, DOI: 10.1039/C7FO00830A

<details open>
<summary>BibTeX entry</summary>

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

Copyright (c) 2017-2026, Alberto Tonda \<alberto.tonda@gmail.com\>

Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted, provided that the above copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

`nlohmann/json` (vendored in `cpp/thirdparty/nlohmann/`) is authored by Niels Lohmann and distributed under the MIT license. `spdlog` (vendored in `cpp/thirdparty/spdlog/`) is authored by Gabi Melman and distributed under the MIT license. `pybind11_json` (vendored in `cpp/thirdparty/pybind11_json/`) is authored by Martin Renou and distributed under the BSD 3-Clause license. [pybind11](https://github.com/pybind/pybind11) itself and [scikit-build-core](https://github.com/scikit-build/scikit-build-core) are build-time-only dependencies (not vendored, resolved automatically by `pip` from `pyproject.toml`), both distributed under permissive licenses (BSD-style and Apache 2.0, respectively). The original tinyxml library (no longer used, kept out of the repository) was authored by Lee Thomason, Yves Berquin, and Andrew Ellerton.

In case you need help, advice, or you notice a bug, please contact **Alberto Tonda** \<alberto.tonda@gmail.com\> or open a Pull Request on GitHub.
