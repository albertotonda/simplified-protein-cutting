README
======
In order to compile and run this software, you will need CMAKE installed. The software has been tested under Ubuntu (14.04, then 16.04), but in theory the code is written in ISO C++ and should thus be cross-compiling using CMAKE, even under Windows or Mac OS.

In case you need help, advice, or you notice a bug, please contact Alberto Tonda <alberto.tonda@gmail.com>

If you use this software in your publications, please cite:

Tonda, Alberto and Grosvenor, Anita J and Clerens, Stefan and Le Feunteun, Steven, 
"In silico modeling of protein hydrolysis by endoproteases: a case study on pepsin digestion of bovine lactoferrin", 
Food & Function, 2017, DOI: 10.1039/C7FO00830A

What is this software?
----------------------
In a nutshell, this software simulates the action of the pepsin enzyme on several copies of a protein (protein structure given in input). The behavior of pepsin is considered stochastic, and during the simulation it will cut bonds with a certain probability, depending on the amino-acids to the left and right of a bond (positions P4, P3, P2, P1, P1', P2', P3', P4'). The probabilities are computed on the basis of analyses performed by Hamuro et al., 2008 (DOI:10.1002/rcm.3467) and Powers et al., 1977 (DOI: 10.1007/978-1-4757-0719-9_9).

For more information, refer to: Tonda, Alberto and Grosvenor, Anita J and Clerens, Stefan and Le Feunteun, Steven, "In silico modeling of protein hydrolysis by endoproteases: a case study on pepsin digestion of bovine lactoferrin", Food & Function, 2017, DOI: 10.1039/C7FO00830A

In order to parse JSON configuration files, the nlohmann/json library is vendored (as a single header) in src/thirdparty/nlohmann/. We did not write it, this excellent library is authored by Niels Lohmann; for more information, see https://github.com/nlohmann/json
(Until 2026, configuration files were XML, parsed with the tinyxml library; the format was switched to JSON, see below.)

Compiling instructions for Linux
--------------------------------
Install CMAKE (under Ubuntu, 'sudo apt-get install cmake'). 

Once you have CMAKE installed, you have two options: you can either run './test.sh', a shell script that will perform all the steps below, compile the program and run it; or go step-by-step.

If you decided to go step-by-step, run 'cmake .' in the main directory of the project. Then, run 'make'. The result of this procedure should be an executable called 'protein-cutting'.

In order to run a simulation, you will need a .JSON file that describes the structure of the protein, and the probabilities that the simulated pepsin enzyme will use to cut the bonds. In the data directory a sample file is provided, "lactoferrin.json": this file contains the structure of bovine lactoferrin, and it will simulate cutting 500 copies of the protein, using probabilities taken from Hamuro et al., 2008 and Powers et al., 1977 (see above for the DOIs).

To run a simulation, call './protein-cutting --input lactoferrin.json'

THE SIMULATION WILL TAKE SOME TIME TO COMPLETE (probably minutes to tens of minutes, depending on your hardware). ALSO, AT THE END, THE PROGRAM WILL TAKE A LONG TIME TO WRITE THE RESULT TO DISK. DON'T QUIT AT THAT POINT.

In the end, the program will produce a file called 'statistics.csv', that marks the quantity of each peptide during 'time'. On the columns, you have 'time' (number of iterations until that moment), 'time2' (number of cuts so far), 'pepsin' (it's always 1.0, but it's there for future developments), then all peptides produced during the simulations, in alphabetical order. In each row, there's the number of each type of peptide for that iteration. The resulting file is usually large (~73 MB for the example), and you should probably write another script to extract the meaningful information from it (analyzing it manually would take a long time).

The program has a lot of output that is used mainly for debugging purposes (yes, I am guilty of 'printf debugging'), but the good news is that the output to screen can be used to understand what it is doing at each iteration.

Using the software for your work
--------------------------------
In order to use this software for your purpose, you don't need to modify the source code: the behavior of the simulations can be altered by simply editing the JSON configuration file. The sample file 'lactoferrin.json' is heavily commented (the parser accepts "//" and "/* */" comments), I would suggest starting from that and modifying it for your purpose, changing the protein(s) or the probabilities.

License
-------
Copyright (c) 2017, Alberto Tonda <alberto.tonda@gmail.com>

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

If you use this software in your publications, please cite:

Tonda, Alberto and Grosvenor, Anita J and Clerens, Stefan and Le Feunteun, Steven, 
"In silico modeling of protein hydrolysis by endoproteases: a case study on pepsin digestion of bovine lactoferrin", 
Food & Function, 2017, DOI: 10.1039/C7FO00830A

