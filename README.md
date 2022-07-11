
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Penalized Decomposition Using Residuals (PeDecURe)

PeDecURe can be used for feature extraction that with built-in
adjustment for nuisance variables, such as confounders. Our method
identifies sources of variation in that data that are shared between
features (e.g., measurements derived from neuroimaging scans) and an
outcome of interest (e.g., diagnosis), while substantially reducing
overlap with information about nuisance variables (e.g., age or sex).

Researchers often use methods like principal components analysis (PCA)
or partial least squares (PLS) for feature extraction, but these methods
do not adjust for nuisance variables, such as confounders, which may
hinder generalizability or interpretability of downstream analyses that
incorporate those features. In our
[manuscript](https://www.biorxiv.org/content/10.1101/2022.01.27.477859v1),
we introduce the intuition behind our method and illustrate that
features extracted using PeDecURe are predictive of an outcome of
interest, have low correlations with nuisance variables, and show
promise for out-of-sample generalizability.

<!-- badges: start -->
<!-- badges: end -->

## Installation

``` r
# install.packages("devtools")
devtools::install_github("smweinst/PeDecURe")
```

## Example

``` r
library(PeDecURe)

# get residuals:
resid.dat = get.resid(X,Y,A)
X.star = resid.dat$X.star
X.tilde = resid.dat$X.tilde

# tune lambda:
lambda.tune = pedecure.tune(X.orig = X,
                            X.max = X.star,
                            X.penalize = X.penalize,
                            lambdas = seq(0,10,by=0.1),
                            A = A,
                            Y = Y,
                            nPC = 3)
best.lambda = lambda.tune$lambda_tune

# run pedecure:
pedecure.out = pedecure(X = X.star,
                        X.penalize = X.tilde,
                        A = A,
                        Y = Y,
                        lambda = best.lambda,
                        nPC = 3)
```
