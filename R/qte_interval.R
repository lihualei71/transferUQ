#' Lower confidence bound for uhat defined in the note
#' via Hoeffding-Bentkus-Maurer inequalities
#'
#' @param beta quantile
#' @param n sample size
#' @param k order of U-statistic
#' @param alpha confidence level
#'
#' @noRd
#' @keywords internal
HBM_quantile_lower <- function(beta, n, k, alpha){
    inv_ustat_HBM_upper(beta, n, k, alpha, type = "x")
}

#' Upper confidence bound for uhat defined in the note
#' via Hoeffding-Bentkus-Maurer inequalities
#'
#' @param beta quantile
#' @param n sample size
#' @param k order of U-statistic
#' @param alpha confidence level
#'
#' @noRd
#' @keywords internal
HBM_quantile_upper <- function(beta, n, k, alpha){
    1 - inv_ustat_HBM_upper(1 - beta, n, k, alpha, type = "x")
}

#' Confidence interval for quantile transfer errors
#'
#' \code{qte_interval} generates one- or two-sided forecast
#' intervals for quantile transfer errors (Section 5.2) with
#' a given coverage level. 
#'
#' @details \code{err} should be in one of the following forms.
#' \itemize{
#' \item A square asymmetric matrix when the number of training domains \eqn{r = 1}. In this case, \code{err[i, j]} records the transfer error using the i-th domain for training and j-th domain for testing.
#' \item A data.frame with \eqn{(r+2)} columns. The last column records the transfer errors, the second last column records the indices of test domains, and the first r columns record the indices of training domains.
#' \item (Not recommended) a data.frame with two columns. The second column records the transfer errors and the first column records the indices of test domains. 
#' }
#' For this version, it is recommended to include the transfer errors for all n-choose-(r+1) combinations of training and test domains.
#' Otherwise, \code{ete_interval} can still
#' output an interval, though the theoretical guarantee is unclear.
#'
#' The input \code{r} should be specified only when \code{err} is a data.frame with two columns (i.e., without training domain indices).
#' 
#' @param err error object; see Details
#' @param coverage target coverage level
#' @param side "two" for two-sided, "left" for intervals of form \eqn{[a, \infty)}, "right" for intervals of form \eqn{(-\infty, a]}
#' @param betas a vector of quantiles at which confidence intervals are computed; \code{seq(0.2, 0.8, 0.01)} by default
#' @param r number of training domains (optional); see Details
#' 
#' @return a data.frame with three columns: \code{beta} for the quantile, \code{lb} for the lower confidence bound, and \code{ub} for the upper confidence bound
#'
#' @examples
#' \donttest{# Generate an error matrix
#' n <- 100
#' set.seed(1)
#' errmat <- matrix(runif(n^2), nrow = n)
#' diag(errmat) <- NA
#'
#' # Two-sided confidence intervals for 
#' qte_interval(errmat, 0.9, "two")
#' 
#' # One-sided confidence intervals
#' qte_interval(errmat, 0.9, "left")
#' qte_interval(errmat, 0.9, "right")
#'
#' # Generate an error data.frame with r+2=3 columns
#' idx <- expand.grid(train = 1:n, test = 1:n)
#' errdf <- cbind(idx, data.frame(val = as.numeric(errmat)))
#' errdf <- errdf[!is.na(errdf$val), ]
#' qte_interval(errdf, 0.9, "two")
#' qte_interval(errmat, 0.9, "two")
#'
#' # Generate an error data.frame with only two columns (not recommended)
#' errdf2 <- errdf[, 2:3]
#' qte_interval(errdf2, 0.9, "two", r = 1)
#'
#' # Generate an error data.frame with r>1
#' n <- 40
#' set.seed(1)
#' idx <- expand.grid(train1 = 1:n, train2 = 1:n, test = 1:n)
#' idx <- idx[idx$train1 != idx$train2 & idx$train1 != idx$test & idx$train2 != idx$test, ]
#' err <- cbind(idx, data.frame(val = runif(nrow(idx))))
#'
#' # Two-sided confidence interval
#' qte_interval(err, 0.9, "two")
#' 
#' # One-sided confidence intervals
#' qte_interval(err, 0.9, "left")
#' qte_interval(err, 0.9, "right")
#' }
#' @export
qte_interval <- function(err, coverage,
                         side = c("two", "left", "right"),
                         betas = seq(0.2, 0.8, 0.01),
                         r = NULL){
    side <- side[1]
    err_obj <- standardize_err(err)
    err <- err_obj$err
    n <- err_obj$n    
    if (is.null(r)){
        r <- err_obj$r
        if (is.na(r)){
            stop("Missing r. Add training domain indices into err or add r as an input.")
        }
    } else if (!is.na(err_obj$r) && err_obj$r != r){
        stop("The input r does not match the number of training domains in err.")
    }
    
    res <- sapply(betas, function(beta){
        if (side == "two"){
            alpha <- (1 - coverage) / 2
            beta_minus <- HBM_quantile_lower(beta, n, r, alpha)
            beta_plus <- HBM_quantile_upper(beta, n, r, alpha)
            lb <- quantile(err$val, probs = beta_minus, type = 1)
            ub <- quantile(err$val, probs = beta_plus, type = 1)
        } else if (side == "left"){
            alpha <- 1 - coverage
            beta_minus <- HBM_quantile_lower(beta, n, r, alpha)
            lb <- quantile(err$val, probs = beta_minus, type = 1)
            ub <- Inf
        } else if (side == "right"){
            alpha <- 1 - coverage
            beta_plus <- HBM_quantile_upper(beta, n, r, alpha)
            lb <- -Inf
            ub <- quantile(err$val, probs = beta_plus, type = 1)
        }
        c(beta, lb, ub)
    })
    res <- data.frame(t(res))
    names(res) <- c("beta", "lb", "ub")
    return(res)
}
