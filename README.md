## Instruction

We proposed a multi-kernel framework with boosted microbiome distances
for classification with microbiome data.

#### The source codes and examples are available [**Here**](https://github.com/HXu06/Multi_kernel-microbiome).

## Requirements

The following R packages are required:

-   `tibble`, `pRoc`, `Rcpp`,`MiSPU`

You could download them directly in CRAN through the following commands
in your R console.

    install.packages(c('tibble', 'pRoc','Rcpp','MiSPU'))

## Setup

Before running for the first time, you need to compile the C extensions.
To finish the compiling, you should install `Rcpp` package first, then
run the following code:

    library(Rcpp)
    sourceCpp("GUniFrac.cpp")

## Usage

-   `mklpre.R` is the main program to perform prediction.
-   `allfunctions.R` contains all functions needed for prediction.
-   `example.RData` is a sample data set.

### Input

-   subOTUtable: a large N by p matrix with each row representing an
    individual; each column representing an OTU. OTU counts data of both
    training samples and testing samples are included.

<!-- -->

    subOTUtable[1:5,1:10]

    ##   OTU1883 OTU3114 OTU1483 OTU2576 OTU4408 OTU2013 OTU2859 OTU1671 OTU73 OTU2510
    ## 1       0       0       0      79       2       0       0       7     0       0
    ## 2       0       0       0       0       0       0       0       5     7       0
    ## 3       0       6       0       0       0       0       0       5     0       0
    ## 4       0       1       0       0       0       0       0       0     0       0
    ## 5       0       6       0       0       0       0       0       1     0       0

-   subtree: a phylogenetic tree that captures evolutionary
    relationships among species. OTUs that are close on a phylogenetic
    tree is refereed as phylogenetically-related OTUs, which are usually
    also functionally related.

<!-- -->

    attributes(subtree)

    ## $names
    ## [1] "edge"        "Nnode"       "tip.label"   "edge.length"
    ## 
    ## $class
    ## [1] "phylo"

-   y: a binary outcome for training data.

-   Lcov: kernel list for covariates for both training and testing data

-   tran\_ind: index of training data in subOTUtable and Lcov.

-   test\_ind: index of testin data in subOTUtable and Lcov.

-   rho\_ls: a grid of rho\_ls, default is data-driven, or you can give
    a grid by yourself.

-   grid: (grid+1) values are given to choose *ρ*, default is 10.

-   kfold: kfold cross-validation is used to choose *ρ*, default is 5.
    \### Output

-   The results are output as a list named as follows:

<!-- -->

    resl <- mklpre(subOTUtable = subOTUtable, subtree = subtree, y=y,  Lcov = Lcov,tran_ind = tran_ind, test_ind= test_ind)
    names(resl)

    ## [1] "weight"     "outcome_pr" "rho"

-   Kernel weights, prediction probabilities of membership and value of
    *ρ* determined by CV are reported. Given the probabilities and the
    true status of testing data, we can calulate the Area under the ROC
    Curve (AUC).

<!-- -->

    resl$weight

    ##           [,1]
    ## [1,] 0.1973882
    ## [2,] 0.1950701
    ## [3,] 0.2010016
    ## [4,] 0.1955624
    ## [5,] 0.2109777

    resl$outcome_pr[1:10]

    ##  [1] 0.81457366 0.02710860 0.90090740 0.15995356 0.09322434 0.99144924
    ##  [7] 0.15930572 0.91354733 0.35377634 0.52018647

    resl$rho

    ## [1] 52623.8

    pROC::roc(y_test, resl$outcome_pr, levels = c(0,1), direction = '<')$auc

    ## Area under the curve: 0.7523
