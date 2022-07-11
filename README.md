
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Penalized Decomposition Using Residuals (PeDecURe)

PeDecURe provides feature extraction with built-in adjustment for
nuisance variables. Our method identifies sources of variation in that
data that are shared between features (e.g., measurements derived from
neuroimaging scans) and an outcome of interest (e.g., diagnosis), while
substantially reducing overlap with information about nuisance variables
(e.g., age or sex).

In our
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

*Notation:*

-   `X` : feature matrix
    ![(n \times p)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28n%20%5Ctimes%20p%29 "(n \times p)")
-   `A` : matrix of nuisance variables
    ![(n \times q)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28n%20%5Ctimes%20q%29 "(n \times q)")
-   `Y` : vector (length
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n"))
    with outcome labels (e.g., disease group)

*Implementation using `PeDecURe` R package:*

``` r
library(PeDecURe)

# get residuals:
resid.dat = get.resid(X,Y,A)
X.star = resid.dat$X.star
X.tilde = resid.dat$X.tilde

# tune lambda:
lambda.tune = pedecure.tune(X.orig = X,
                            X.max = X.star,
                            X.penalize = X.tilde,
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

# PC scores - these are our new features that can be used for an association study, predictive model, etc.
## note: X should be centered by column
PC.scores = X%*%pedecure.out$vectors

# Look at correlations between the first few PC scores and the nuisance variables (A1, A2, Y)
cor.scores = partial.cor(PC.scores, A, Y)
scores.partial.cor = cor.scores$partial$estimates
scores.marginal.cor = cor.scores$marginal$estimates
```

### Applying PeDecURe output in a new sample

*PC scores in new sample: multiply new feature matrix by PC loadings
from above.*

-   `X.test` : feature matrix in a different sample
    ![(m \times p)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28m%20%5Ctimes%20p%29 "(m \times p)")

``` r
PC.scores.test = X.test%*%pedecure.out$vectors
# note: pedecure.out$vectors was the output from applying PeDecURe in training sample above
```

*If
![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
and
![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
are observed in the new sample, we can also look at their correlations
with the PC scores in the test sample:*

-   `A.test` : matrix of nuisance variables in the new sample
    ![(m \times q)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28m%20%5Ctimes%20q%29 "(m \times q)")
    (if observed)
-   `Y.test` : vector (length
    ![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m"))
    with outcome labels (if observed)

``` r
cor.scores.test = partial.cor(PC.scores.test, A.test, Y.test)
scores.partial.cor.test = cor.scores.test$partial$estimates
scores.marginal.cor.test = cor.scores.test$marginal$estimates
```
