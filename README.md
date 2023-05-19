# transferUQ
An R package for uncertainty quantification of transfer performance for prediction models

## Overview
This R package implements the forecast intervals and confidence intervals on transfer errors of prediction models using observations from multiple domains, along with the comparison between two prediction models in terms of worst-case- and everywhere-dominance. All methods are introduced in our paper: [The Transfer Performance of Economics Models](https://www.dropbox.com/s/f4elounoyqu7tty/Theory_Transfer.pdf?dl=0). It also includes the dataset `CertaintyEquivalents` analyzed in the paper, which includes unnormalized and normalized transfer error matrices on certainty equivalents for binary lotteries from 44 domains. The descriptions of main functions can be found [here]()

## Installation         

```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("lihualei71/transferUQ")
```
To access the vignette, run the following code to build it. 
```
devtools::install_github("lihualei71/transferUQ", build_vignettes = TRUE)
```

The package requires [lpSolve](https://cran.r-project.org/web/packages/lpSolve/index.html) to be installed. 
