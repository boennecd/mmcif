#' Computes the Log Cholesky Decomposition and the Inverse
#'
#' @param x A positive definite matrix or a vector with a log Cholesky
#' decomposition.
#'
#' @examples
#' set.seed(1)
#' S <- drop(rWishart(1, 10, diag(10)))
#' stopifnot(isTRUE(all.equal(S, log_chol_inv(log_chol(S)))))
#'
#' @export
log_chol <- function(x){
  x <- chol(x)
  diag(x) <- log(diag(x))
  x[upper.tri(x, TRUE)]
}

#' @rdname log_chol
#' @export
log_chol_inv <- function(x){
  dim <- (sqrt(8 * length(x) + 1) - 1) / 2
  out <- matrix(0, dim, dim)
  out[upper.tri(out, TRUE)] <- x
  diag(out) <- exp(diag(out))
  crossprod(out)
}
