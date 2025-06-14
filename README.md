# XgeneR

This repository contains the code written by Maria Carilli to run analysis of _cis_ and _trans_ regulatory differences between homozygous strains using RNA-seq data from parents of both strains and their hybrid crosses, as well as code to visualize the results. It requires that replicates are available. 

- Source code is in the `R/` directory.  
- Example vignettes for single-condition or multiple-condition regulation are available in the `vignettes/` directory.

The preprint **"Estimating cis and trans contributions to differences in gene regulation"** (Hallgrímsdóttir et al., 2024) describes the theory behind the model, statistical tests, and visualization.

## Installation

XgeneR can be installed in R with the following command:

```r
devtools::install_github("pachterlab/XgeneR")
```
