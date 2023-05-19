#' The LP subproblem for everywhere dominance
#' (Eq. (Q.4) in Appendix Q.2)
#'
#' @noRd
#' @keywords internal
lp_subproblem <- function(Psi1, Psi2, tau){
    n <- length(Psi1)
    objective <- Psi1
    constraint <- rbind(Psi2, diag(rep(1, n)), rep(1, n))
    dir <- c(rep(">=", n+1), "=")
    rhs <- c(tau, rep(0, n), 1)
    lpSolve::lp("min", objective, constraint, dir, rhs)$objval / (n - 1)
}

#' Check if method 1 everywhere-upper-dominates method 2
#' 
#' \code{check_everywhere_dominance} evaluates whether method 1
#' everywhere-upper-dominates method 2 at the tau-th quantile
#' (Appendix Q.2).
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
#'
#' @return a data.frame with two columns: the first column records the values of \eqn{\tau} and the second column records whether method 1 everywhere-upper-dominates method 2 at each quantile
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
#' check_everywhere_dominance(errmat1, errmat2)
#' check_everywhere_dominance(errmat1, errmat3)
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
#' check_everywhere_dominance(errdf1, errdf2)
#' check_everywhere_dominance(errdf1, errdf3)
#'
#' # Transform error matrices into data.frames with only two columns (not recommended)
#' errdf1_2 <- errdf1[, 2:3]
#' errdf2_2 <- errdf2[, 2:3]
#' errdf3_2 <- errdf3[, 2:3]
#'
#' check_everywhere_dominance(errdf1_2, errdf2_2)
#' check_everywhere_dominance(errdf1_2, errdf3_2)
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
#' check_everywhere_dominance(err1, err2)
#' check_everywhere_dominance(err1, err3)
#' }
#' @export
check_everywhere_dominance <- function(err1, err2,
                                       taulist = seq(0.5, 1, 0.01)){
    err_obj1 <- standardize_err(err1)
    err1 <- err_obj1$err
    n <- err_obj1$n
    
    err_obj2 <- standardize_err(err2)
    err2 <- err_obj2$err
    if (err_obj2$n != n){
        stop("The number of domains do not match!")
    }

    res <- sapply(taulist, function(tau){
        jstart <- ceiling(tau * n * (n - 1))
        flag <- TRUE
        for (j in jstart:(n*(n-1))){
            val <- err2$val[j]
            idx <- which(err1$val >= val)
            if (length(idx) == 0){
                ## All realized errors of method 1 are smaller
                return(TRUE)
            }
            idx <- min(idx) - 1
            Psi2 <- as.numeric(table(as.factor(err2$test)[1:j]))
            Psi1 <- as.numeric(table(as.factor(err1$test)[1:idx]))
            if (all(Psi2 <= Psi1)){
                next
            }
            if (min(Psi1) / (n - 1) > tau){
                next
            } else if (lp_subproblem(Psi1, Psi2, tau) <= tau){
                flag <- FALSE
                break
            }
        }
        return(flag)
    })
    return(data.frame(tau = taulist, res = res))    
}
