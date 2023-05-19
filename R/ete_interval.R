#' Lower confidence bound of mean of U-statistics of order k
#' via Hoeffding-Bentkus-Maurer inequalities
#'
#' @param U the value of U-statistic
#' @param n sample size
#' @param k order of U-statistic
#' @param alpha confidence level
#'
#' @noRd
#' @keywords internal
HBM_mean_upper <- function(U, n, k, alpha){
    U <- pmin(pmax(U, 1e-10), 1-1e-10)
    inv_ustat_HBM_upper(U, n, k, alpha, type = "mu")
}

#' Upper confidence bound of mean of U-statistics of order k
#' via Hoeffding-Bentkus-Maurer inequalities
#'
#' @param U the value of U-statistic
#' @param n sample size
#' @param k order of U-statistic
#' @param alpha confidence level
#'
#' @noRd
#' @keywords internal
HBM_mean_lower <- function(U, n, k, alpha){
    1 - HBM_mean_upper(1 - U, n, k, alpha)
}

#' Confidence interval for expected transfer errors 
#'
#' \code{ete_interval} generates one- or two-sided forecast
#' intervals for expected transfer errors (Section 5.3) with
#' a given coverage level. This function requires the transfer
#' error to be bounded by [0, \code{bd}]. 
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
#' @param bd upper bound of the loss function; 1 by default
#' @param r number of training domains (optional); see Details 
#'
#' @return a 2-d vector giving the left and right ends of the confidence interval
#'
#' @examples
#' \donttest{# Generate an error matrix
#' n <- 100
#' set.seed(1)
#' errmat <- matrix(runif(n^2), nrow = n)
#' diag(errmat) <- NA
#'
#' # Two-sided confidence interval
#' ete_interval(errmat, 0.9, "two", 1)
#' 
#' # One-sided confidence intervals
#' ete_interval(errmat, 0.9, "left", 1)
#' ete_interval(errmat, 0.9, "right", 1)
#'
#' # Generate an error data.frame with r+2=3 columns
#' idx <- expand.grid(train = 1:n, test = 1:n)
#' errdf <- cbind(idx, data.frame(val = as.numeric(errmat)))
#' errdf <- errdf[!is.na(errdf$val), ]
#' ete_interval(errdf, 0.9, "two", 1)
#' ete_interval(errmat, 0.9, "two", 1)
#'
#' # Generate an error data.frame with only two columns (not recommended)
#' errdf2 <- errdf[, 2:3]
#' ete_interval(errdf2, 0.9, "two", 1, r = 1)
#'
#' # Generate an error data.frame with r>1
#' n <- 40
#' set.seed(1)
#' idx <- expand.grid(train1 = 1:n, train2 = 1:n, test = 1:n)
#' idx <- idx[idx$train1 != idx$train2 & idx$train1 != idx$test & idx$train2 != idx$test, ]
#' err <- cbind(idx, data.frame(val = runif(nrow(idx))))
#'
#' # Two-sided confidence interval
#' ete_interval(err, 0.9, "two", 1)
#' 
#' # One-sided confidence intervals
#' ete_interval(err, 0.9, "left", 1)
#' ete_interval(err, 0.9, "right", 1)
#' }
#' @export
ete_interval <- function(err, coverage,
                         side = c("two", "left", "right"),
                         bd = 1,
                         r = NULL){
    side <- side[1]
    err_obj <- standardize_err(err)
    err <- err_obj$err / bd
    n <- err_obj$n
    if (is.null(r)){
        r <- err_obj$r
        if (is.na(r)){
            stop("Missing r. Add training domain indices into err or add r as an input.")
        }
    } else if (!is.na(err_obj$r) && err_obj$r != r){
        print(err_obj$r)
        print(r)
        stop("The input r does not match the number of training domains in err.")
    }
    U <- mean(err$val)

    if (side == "two"){
        alpha <- (1 - coverage) / 2
        CI_left <- HBM_mean_lower(U, n, r, alpha)
        CI_right <- HBM_mean_upper(U, n, r, alpha)
    } else if (side == "left"){
        alpha <- 1 - coverage
        CI_left <- HBM_mean_lower(U, n, r, alpha)
        CI_right <- Inf
    } else if (side == "right"){
        alpha <- 1 - coverage
        CI_left <- -Inf
        CI_right <- HBM_mean_upper(U, n, r, alpha)
    }
    return(c(CI_left, CI_right) * bd)
}
