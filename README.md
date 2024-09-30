# powerPLS
 
This package provides a tool to perform **power analysis** in Partial Least Squares (PLS) for classification when two classes are analyzed. 

## Installation

You can install the released version of `powerPLS` with:

``` r
devtools::install_github("angeella/powerPLS")
```

## Quick overview

The main functions are 
- `computeSampleSize()` which estimated the power considering several values of sample size and number of score components.

- `computePower()` which estimated the power considering a fixed sample size and several number of score components.


``` r
datas <- simulatePilotData(nvar = 10, clus.size = c(5,5),m = 6,nvar_rel = 5,A = 2)
out <- computePower(X = datas$X, Y = datas$Y, A = 3, n = 20, test = "R2")
out <- computeSampleSize(X = datas$X, Y = datas$Y, A = 2, A = 3, n = 20, test = "R2")
```

## References

Andreella, A., Finos, L., Scarpa, B. and Stocchero, M. "Towards a power analysis for PLS-based methods" 	arXiv:2403.10289 stat.ME. **link**: https://arxiv.org/abs/2403.10289

## Did you find some bugs?

Please write to angela.andreella[\at]unive[\dot]it or insert a reproducible example using reprex on my issue github page.
