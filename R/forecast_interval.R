sort_2d_offdiag <- function(mat){
    n <- nrow(mat)
    diag(mat) <- NA
    vec <- as.vector(mat)
    res <- sort(vec, na.last = TRUE, index.return = TRUE)
    colid <- ceiling(res$ix / n)
    rowid <- res$ix - (colid - 1) * n
    res <- data.frame(row = rowid, col = colid, val = res$x)
    res[!is.na(res$val), ]
}

Qj_Gamma <- function(Psi, j, Gamma){
    n <- length(Psi)
    Psibar <- cumsum(sort(Psi)) / (1:n)
    min((Psibar - j/n) / (1 + n / (Gamma^2 - 1) / (1:n))) + j/n
}

etau_Gamma <- function(errmat, tau, Gamma){
    n <- nrow(errmat)
    sort_df <- sort_2d_offdiag(errmat)
    jstart <- ceiling(tau * n * (n - 1))
    for (j in jstart:(n*(n-1))){
        Psi <- as.numeric(table(as.factor(sort_df$col)[1:j]))
        Qj <- Qj_Gamma(Psi, j, Gamma)
        if (Qj >= tau * n){
            break
        }
    }
    return(sort_df$val[j])
}

#' The worst-case transfer error \eqn{\bar{e}_{\tau}(\Gamma)} as
#' a function of \eqn{\tau} for a given \eqn{\Gamma}.
#'
#' @param errmat a square matrix with \code{errmat[i, j]} being the transfer error using the i-th domain for training and j-th domain for testing
#' @param Gamma the Gamma value for selection bias
#'
#' @return a data.frame with two columns: \code{prob} for all turning points of \eqn{\tau} where \eqn{\bar{e}_{\tau}(\Gamma)} changes and \code{val} for the corresponding bounds
#'
#' @export
Gamma_curve <- function(errmat, Gamma){
    n <- nrow(errmat)
    m <- n * (n - 1)
    sort_df <- sort_2d_offdiag(errmat)
    Qvec <- rep(m)
    for (j in 1:m){
        Psi <- as.numeric(table(as.factor(sort_df$col)[1:j]))
        Qvec[j] <- Qj_Gamma(Psi, j, Gamma)
    }
    data.frame(prob = Qvec / (n - 1), val = sort_df$val)
}

#' Check if method 1 everywhere-upper-dominates method 2
#' at the tau-th quantile
#'
#' @param errmat1 a square matrix with \code{errmat[i, j]} being the transfer error using the i-th domain for training and j-th domain for testing for method 1
#' @param errmat2 the error matrix for method 2
#' @param tau the quantile at which the methods are compared
#'
#' @return a logical indicating whether method 1 everywhere-upper-dominates method 2 at the tau-th quantile
#' 
#' @export
check_everywhere_dominance <- function(errmat1, errmat2, tau){
    n <- nrow(errmat1)
    err1 <- sort_2d_offdiag(errmat1)
    err2 <- sort_2d_offdiag(errmat2)
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
        Psi2 <- as.numeric(table(as.factor(err2$col)[1:j]))
        Psi1 <- as.numeric(table(as.factor(err1$col)[1:idx]))
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
}

lp_subproblem <- function(Psi1, Psi2, tau){
    n <- length(Psi1)
    objective <- Psi1
    constraint <- rbind(Psi2, diag(rep(1, n)), rep(1, n))
    dir <- c(rep(">=", n+1), "=")
    rhs <- c(tau, rep(0, n), 1)
    lpSolve::lp("min", objective, constraint, dir, rhs)$objval / (n - 1)
}
