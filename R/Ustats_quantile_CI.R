#' Upper/Lower confidence bound for uhat defined in the note
#' via Hoeffding-Bentkus-Maurer inequalities
#'
#' @param beta quantile
#' @param n sample size
#' @param k order of U-statistic
#' @param alpha confidence level
HBM_quantile_lower <- function(beta, n, k, alpha){
    inv_ustat_HBM_upper(beta, n, k, alpha, type = "x")
}

HBM_quantile_upper <- function(beta, n, k, alpha){
    1 - inv_ustat_HBM_upper(1 - beta, n, k, alpha, type = "x")
}

#' Two-sided confidence interval for the quantile of
#' transfer errors 
#'
#' @param errmat a square matrix with \code{errmat[i, j]} being the transfer error using the i-th domain for training and j-th domain for testing
#' @param alpha confidence level
#' @param betas a vector of quantiles at which confidence intervals are computed
#'
#' @return a data.frame with three columns: \code{beta} for the quantile, \code{lb} for the lower confidence bound, and \code{ub} for the upper confidence bound.
#' 
#' @export
QTE_interval <- function(errmat, alpha,
                         betas = seq(0.2, 0.8, 0.01)){
    n <- nrow(errmat)
    res <- sapply(betas, function(beta){
        beta_plus <- HBM_quantile_upper(beta, n, 2, alpha / 2)
        beta_minus <- HBM_quantile_lower(beta, n, 2, alpha / 2)
        diag(errmat) <- NA
        errvec <- as.vector(errmat)
        errvec <- errvec[!is.na(errvec)]
        lb <- quantile(errvec, probs = beta_minus, type = 1)
        ub <- quantile(errvec, probs = beta_plus, type = 1)
        c(beta, lb, ub)
    })
    res <- data.frame(t(res))
    names(res) <- c("beta", "lb", "ub")
    return(res)
}
