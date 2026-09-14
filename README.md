# Simplified in-silico protein cutting
Modeling protein hydrolysis and release of peptides by endoproteases requires complex simulations, typically taking into account the 3D structure of both the enzymes and the target protein. Such structures can sometimes be difficult to predict starting from the protein's acido-aminic sequence.

This repository contains the code for an alternative approach, published in (Tonda et al. (2017) _In silico modeling of protein hydrolysis by endoproteases: a case study on pepsin digestion of bovine lactoferrin_, Food & Function, Vol. 8, Issue 12, DOI: 10.1039/C7FO00830A)[https://pubs.rsc.org/fo/article-abstract/8/12/4404/566726/In-silico-modeling-of-protein-hydrolysis-by]. The idea is to just consider from the linear sequence of amino-acids, and then simulate the behavior of an enzyme starting from the frequency of cuts observed during previous experiments. The final peptides obtained by the simulation are qualitatively coherent with real-world experiments, even though the exact absolute quantities might be different.

## Installation instructions
The original code is in C++, and is contained in the `/src` subfolder.