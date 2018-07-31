#' Modified KernelMatrix function to produce less errors
#' @param kernel The kernel use from kernlab. The available kernels are
#'   "rbfdot", "polydot", "vanilladot", "tanhdot", "laplacedot", "besseldot" or
#'   "anovadot".
#' @param x The data.
#' @param y Currently not used.
#' @importFrom kernlab as.kernelMatrix
#' @export
kernelMat <- function(kernel, x, y = NULL) {
  if (is(x, "vector"))
    x <- as.matrix(x)
  if (is(y, "vector"))
    y <- as.matrix(y)
  if (!is(x, "matrix"))
    stop("x must be a matrix")
  if (!is(y, "matrix") && !is.null(y))
    stop("y must be a matrix")
  n <- nrow(x)
  res <- matrix(rep(0, n * n), ncol = n)
  if (is.null(y)) {
    for (i in 1:n) {
      for (j in i:n) {
        res[i, j] <- as.numeric(kernel(x[i, ], x[j, ]))
      }
    }
    res <- res + t(res)
    diag(res) <- diag(res)/2
  }
  if (is(y, "matrix")) {
    m <- dim(y)[1]
    res <- matrix(0, dim(x)[1], dim(y)[1])
    for (i in 1:n) {
      for (j in 1:m) {
        res[i, j] <- kernel(x[i, ], y[j, ])
      }
    }
  }
  return(as.kernelMatrix(res))
}
