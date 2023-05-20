# transferUQ
An R package for uncertainty quantification of transfer performance for prediction models

## Overview
This R package implements the forecast intervals and confidence intervals on transfer errors of prediction models using observations from multiple domains, along with the comparison between two prediction models in terms of worst-case- and everywhere-dominance. All methods are introduced in our paper: [The Transfer Performance of Economics Models](https://lihualei71.github.io/Theory_Transfer.pdf). It also includes the dataset `CertaintyEquivalents` analyzed in the paper, which includes unnormalized and normalized transfer error matrices on certainty equivalents for binary lotteries from 44 domains. A full description of all functions can be found in the manual or the [package site] ()

## Installation         

```
if (!require("devtools")){
    install.packages("devtools")
}
devtools::install_github("lihualei71/transferUQ")
```
The package requires [lpSolve](https://cran.r-project.org/web/packages/lpSolve/index.html) to be installed. 

To access the vignette, run the following code to build it. 
```
devtools::install_github("lihualei71/transferUQ", build_vignettes = TRUE)
```
Make sure that [knitr](https://cran.r-project.org/web/packages/knitr/index.html), [rmarkdown](https://rmarkdown.rstudio.com/), [latex2exp](https://cran.r-project.org/web/packages/latex2exp/index.html), and [tidyverse](https://www.tidyverse.org/) are installed.

## Usage Examples
We illustrate the usage of transferUQ package using the `CertaintyEquivalents` dataset. For details please read the vignette (`vignette("transferUQ_demo", package = "transferUQ")`) and the manual.

```
library("transferUQ")

# Extract the unnormalized transfor error matrix for the CPT model and Random Forests
err_CPT <- CertaintyEquivalents$unnormalized[["abdg"]]
err_RF <- CertaintyEquivalents$unnormalized[["RF"]]

# Compute 80% two-sided forecast intervals
forecast_interval(err_CPT, 0.8)
forecast_interval(err_RF, 0.8)

# Compute 80% two-sided confidence intervals for the median transfer errors
qte_interval(err_CPT, 0.8, betas = 0.5)
qte_interval(err_RF, 0.8, betas = 0.5)

# Check if CPT worst-case-upper-dominates RF
check_worstcase_dominance(err_CPT, err_RF)

# Check if CPT everywhere-upper-dominates RF
check_everywhere_dominance(err_CPT, err_RF)
```


