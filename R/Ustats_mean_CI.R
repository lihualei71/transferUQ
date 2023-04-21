#' Upper/Lower confidence bound of mean of U-statistics of order k
#' via Hoeffding-Bentkus-Maurer inequalities
#'
#' @param U the value of U-statistic
#' @param n sample size
#' @param k order of U-statistic
#' @param alpha confidence level
#' 
HBM_mean_upper <- function(U, n, k, alpha){
    U <- pmin(pmax(U, 1e-10), 1-1e-10)
    inv_ustat_HBM_upper(U, n, k, alpha, type = "mu")
}

HBM_mean_lower <- function(U, n, k, alpha){
    1 - HBM_mean_upper(1 - U, n, k, alpha)
}

#' Two-sided confidence interval for the expectation of
#' transfer errors 
#'
#' @param errmat a square matrix with \code{errmat[i, j]} being the transfer error using the i-th domain for training and j-th domain for testing
#' @param alpha confidence level
#' @param bd the upper bound of the loss function
#'
#' @return a 2-d vector giving the left and right ends of the confidence interval
#' 
#' @export
ETE_interval <- function(errmat, alpha, bd = 1){
    n <- nrow(errmat)
    diag(errmat) <- NA
    errvec <- as.vector(errmat) / bd
    U <- mean(errvec, na.rm = TRUE)
    lb <- HBM_mean_lower(U, n, 2, alpha / 2)
    ub <- HBM_mean_upper(U, n, 2, alpha / 2)    
    c(lb, ub) * bd
}
