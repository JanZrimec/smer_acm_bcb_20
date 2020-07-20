# DNA structural representations

Structural representations of DNA regulatory substrates can enhance sequence-based algorithms by associating functional sequence variants

<img src=https://github.com/JanZrimec/smer_acm_bcb_20/blob/master/docs/smir_fig_smers.png alt="drawing" width="400">
Figure. Schematic depiction of the (A) construction and (B) usage of structural representations. In a structural representation of a given DNA sequence, each central nucleotide position and its neighboring regions define a k-mer from 3 to 9 bp in length, and are encoded as an s-mer with n structural dimensions (S. dim.) that can be defined as a sequence of s-mer cluster centroids.

---------------

Supplementary Table S1 and Figure S1 are available [here](https://github.com/JanZrimec/smer_acm_bcb_20/blob/master/docs/2020_Zrimec-supplements.pdf).

This repository contains scripts to reproduce the figures and analysis. The data is available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3891576.svg)](https://doi.org/10.5281/zenodo.3891576).

Requires Matlab v2017b (Bioinformatics Toolbox v4.9), 
Python v3.6 (numpy v1.15.4, pandas v0.24.2, pyfaidx v0.5.5.2, seaborn v0.9.0, scipy v1.1.0, biopython v1.72) and 
R v3.5 (DNAshapeR v1.10.0), or higher. 

See also the repositories [DNA_structural_variables predictor](https://github.com/JanZrimec/DNA_structural_variables) for prediction of DNA structural properties and
[Non-parametric_multivariate ANOVA_with bootstraps](https://github.com/JanZrimec/NP_MANOVA_bootstrap).

