# XgeneR

This repository contains the code written by Maria Carilli to run analysis of _cis_ and _trans_ regulatory differences between homozygous strains using RNA-seq data from parents of both strains and their hybrid crosses, as well as code to visualize the results. It requires that replicates are available. 

- Source code is in the `R/` directory.  
- Example vignettes for single-condition or multiple-condition regulation are available in the `vignettes/` directory.

A Python port of `XgeneR` is available: see [`XgenePy`](https://github.com/pachterlab/XgenePy).

## Installation

`XgeneR` can be installed in R with the following command:

```r
devtools::install_github("pachterlab/XgeneR")
```

## Reference

The `XgeneR` method is described in 

Ingileif B. Hallgrímsdóttir, Maria Carilli,  Lior Pachter, [Estimating cis and trans contributions to differences in gene regulation](https://www.biorxiv.org/content/10.1101/2024.07.13.603403v2), bioRxiv, 2024, doi.org/10.1101/2024.07.13.603403.
