sort_2d_offdiag <- function(mat){
    n <- nrow(mat)
    diag(mat) <- NA
    vec <- as.vector(mat)
    res <- sort(vec, na.last = TRUE, index.return = TRUE)
    testid <- ceiling(res$ix / n)
    res <- data.frame(test = testid, val = res$x)
    res[!is.na(res$val), ]
}

standardize_err <- function(err){
    if (is.matrix(err)){
        n <- nrow(err)
        r <- 1
        err <- sort_2d_offdiag(err)
    } else if (is.data.frame(err)){
        d <- ncol(err)
        if (d < 2){
            stop("err must have at least two columns with the last column recording the transfer errors, the second last column recording the indices of test domains.")
        } else if (d == 2){
            r <- NA
        } else {
            r <- d - 2
        }
        err <- err[, (d-1):d]
        names(err) <- c("test", "val")
        err <- err[order(err$val), ]
        n <- length(unique(err$test))
    } else {
        stop("err should be a matrix (when r = 1) or a data.frame with the last column recording the transfer errors, the second last column recording the indices of test domains, and, optionally, first r columns recording the indices of training domains")
    }
    return(list(err = err, n = n, r = r))
}
