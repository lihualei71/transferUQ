## Log upper-tail inequalities of U statistics
h1 <- function(y, mu){
    y * log(y / mu) + (1 - y) * log((1 - y) / (1 - mu))
}

ustat_hoeffding_upper <- function(mu, x, n, k){
    n <- floor(n / k)
    -n * h1(pmin(mu, x), mu)
}

ustat_bentkus_upper <- function(mu, x, n, k){
    n <- floor(n / k)
    log(pbinom(ceiling(n * x), n, mu, lower.tail = TRUE)) + 1
}

ustat_maurer_upper <- function(mu, x, n, k){
    lambda_fun <- function(lambda){
        Glambda <- exp(lambda) - lambda - 1
        if (is.nan(Glambda)){
            Glambda <- Inf
        }
        lambda * (x - mu / (k * Glambda / lambda + 1))
    }
    left <- 1e-5
    right <- k
    while (TRUE){
        obj <- optimize(lambda_fun, c(1e-5, k))
        if (obj$minimum >= k){
            right <- 10 * right
            left <- right
            obj <- optimize(lambda_fun, c(left, right))
        } else {
            break
        }
    }
    logprob <- n * obj$objective / k
    logprob
}

ustat_HBM_upper <- function(mu, x, n, k){
    hoeffding_logtail <- ustat_hoeffding_upper(mu, x, n, k)
    bentkus_logtail <- ustat_bentkus_upper(mu, x, n, k)
    maurer_logtail <- ustat_maurer_upper(mu, x, n, k)
    min(hoeffding_logtail, bentkus_logtail, maurer_logtail)
}

#' Inverse of \eqn{B_{n, k}(x; \mu)} at \eqn{alpha} given \eqn{x} or
#' \eqn{\mu} val measures "x" if type = "mu" and "mu" if
#' type = "x". Setting type = "x" yields \eqn{\inf{y: B_{n, k}(y; \mu) \ge \alpha}} and setting type = "mu" yields
#' \eqn{\inf{\nu: B_{n, k}(x; \nu) \ge \alpha}}
#' 
#' @param val the value for either "x" or "mu"
#' @param n sample size
#' @param k order of the U statistic
#' @param alpha confidence level
#' @param type see Details
#'
#' @return the inverse of \eqn{B_{n, k}(x; \mu)}
inv_ustat_HBM_upper <- function(val, n, k, alpha,
                                type = c("x", "mu")){
    type <- type[1]
    if (type == "x"){
        tailprob_x <- function(x){
            ustat_HBM_upper(val, x, n, k) - log(alpha)
        }
        if (tailprob_x(1e-4) > 0){
            0
        } else {
            uniroot(tailprob_x, c(1e-4, val - 1e-10))$root
        }
    } else if (type == "mu"){
        tailprob_mu <- function(mu){
            ustat_HBM_upper(mu, val, n, k) - log(alpha)
        }
        if (tailprob_mu(1 - 1e-4) > 0){
            1
        } else {
            uniroot(tailprob_mu, c(val + 1e-10, 1 - 1e-4))$root
        }
    }
}
