#' \eqn{Q_j(\Gamma)} function defined in Theorem Q.1
#'
#' @noRd
#' @keywords internal
Qj_Gamma <- function(Psi, j, Gamma){
    n <- length(Psi)
    Psibar <- cumsum(sort(Psi)) / (1:n)
    min((Psibar - j/n) / (1 + n / (Gamma^2 - 1) / (1:n))) + j/n
}

#' The worst-case transfer error \eqn{\bar{e}_{\tau}(\Gamma)} as
#' a function of \eqn{\tau} for a given \eqn{\Gamma}.
#'
#' @param err the error object for method 1; see Details
#' @param Gamma the Gamma value for selection bias
#'
#' @return a data.frame with two columns: \code{prob} for all turning points of \eqn{\tau} where \eqn{\bar{e}_{\tau}(\Gamma)} changes and \code{val} for the corresponding bounds
#'
#' @noRd
#' @keywords internal
Gamma_curve <- function(err, Gamma){
    err_obj <- standardize_err(err)
    err <- err_obj$err
    n <- err_obj$n    
    m <- n * (n - 1)
    Qvec <- rep(m)
    for (j in 1:m){
        Psi <- as.numeric(table(as.factor(err$test)[1:j]))
        Qvec[j] <- Qj_Gamma(Psi, j, Gamma)
    }
    data.frame(prob = Qvec / (n - 1), val = err$val)
}

#' Check if method 1 worst-case-upper-dominates method 2
#'
#' \code{check_worstcase_dominance} evaluates whether method 1
#' worst-case-upper-dominates method 2 at the tau-th quantile over
#' \eqn{w\in [\Gamma^{-1}, \Gamma]}. The definition in Appendix Q.1
#' corresponds to \eqn{\Gamma = \infty}. 
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
#' @param err1 error object for method 1; see Details
#' @param err2 error object for method 2; see Details
#' @param taulist a vector of quantiles at which the methods are compared
#' @param Gamma a scalar. \code{Inf} by default
#'
#' @return a data.frame with two columns: the first column records the values of \eqn{\tau} and the second column records whether method 1 worst-case-upper-dominates method 2 at each quantile over \eqn{w\in [\Gamma^{-1}, \Gamma]}
#'
#' @examples
#' \donttest{# Generate three error matrices such that
#' # errmat1 and errmat2 are similar 
#' # and errmat1 is smaller than errmat3
#' n <- 40
#' set.seed(1)
#' errmat1 <- matrix(runif(n^2), nrow = n)
#' errmat2 <- matrix(runif(n^2), nrow = n)
#' errmat3 <- errmat1 + 0.1
#' diag(errmat1) <- diag(errmat2) <- diag(errmat3) <- NA
#' 
#' check_worstcase_dominance(errmat1, errmat2)
#' check_worstcase_dominance(errmat1, errmat3)
#' 
#' # Transform error matrices into data.frames with r+2=3 columns
#' idx <- expand.grid(train = 1:n, test = 1:n)
#' errdf1 <- cbind(idx, data.frame(val = as.numeric(errmat1)))
#' errdf1 <- errdf1[!is.na(errdf1$val), ]
#' errdf2 <- cbind(idx, data.frame(val = as.numeric(errmat2)))
#' errdf2 <- errdf2[!is.na(errdf2$val), ]
#' errdf3 <- cbind(idx, data.frame(val = as.numeric(errmat3)))
#' errdf3 <- errdf3[!is.na(errdf3$val), ]
#'
#' check_worstcase_dominance(errdf1, errdf2)
#' check_worstcase_dominance(errdf1, errdf3)
#'
#' # Transform error matrices into data.frames with only two columns (not recommended)
#' errdf1_2 <- errdf1[, 2:3]
#' errdf2_2 <- errdf2[, 2:3]
#' errdf3_2 <- errdf3[, 2:3]
#'
#' check_worstcase_dominance(errdf1_2, errdf2_2)
#' check_worstcase_dominance(errdf1_2, errdf3_2)
#'
#' # Generate an error data.frame with r>1
#' n <- 20
#' set.seed(1)
#' idx <- expand.grid(train1 = 1:n, train2 = 1:n, test = 1:n)
#' idx <- idx[idx$train1 != idx$train2 & idx$train1 != idx$test & idx$train2 != idx$test, ]
#' err1 <- cbind(idx, data.frame(val = runif(nrow(idx))))
#' err2 <- cbind(idx, data.frame(val = runif(nrow(idx))))
#' err3 <- cbind(idx, data.frame(val = err1$val + 0.1))
#'
#' check_worstcase_dominance(err1, err2)
#' check_worstcase_dominance(err1, err3)
#' }
#' @export
check_worstcase_dominance <- function(err1, err2,
                                      taulist = seq(0.5, 1, 0.01),
                                      Gamma = Inf){
    Gamma_curve1 <- Gamma_curve(err1, Gamma)
    Gamma_curve2 <- Gamma_curve(err2, Gamma)
    res <- sapply(taulist, function(tau){
        id1 <- max(which(Gamma_curve1$prob < tau))
        val1 <- Gamma_curve1$val[id1]
        id2 <- max(which(Gamma_curve2$prob < tau))
        val2 <- Gamma_curve2$val[id2]
        return(val1 < val2)
    })
    return(data.frame(tau = taulist, res = res))
}
