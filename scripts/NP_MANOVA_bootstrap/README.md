# Non-parametric MANOVA (Permanova) with bootstraps

Used in [Zrimec & Lapanje 2018: DNA structure at the plasmid origin-of-transfer indicates its potential transfer range](https://www.nature.com/articles/s41598-018-20157-y)

## Description

<img src="https://github.com/JanZrimec/NP_MANOVA_bootstrap/blob/master/Figure_1.png" width="160">

The Permanova method (Anderson 2001) implemented in Matlab. Enables any distance metric to be used as a measure of variance.

Both the F-test and t-test (pairwise F-test) are included for custom distance matrix input, or sequence input where p-distance (Hamming) is used. Bootstraps of input sequences are generated for the background distribution.

## Usage
See script Script_example_run.m.

```[F, Sw, St] = Anova_F_ratio(seqs,y);```

where:
* seqs .. character array of input sequnces
* y .. grouping variable, array of equal length as seqs
* F .. F-stastic, p value is obtained from bootstrap F ditribution with functions make_bootstraps.m, pvalue_boots.m
* Sw, St .. within group and total sum of squares used to calculate R2 with function get_R2.m
