#' Fits a quantile regression SVM based on the Pinball Loss
#'
#' @param x An n X m matrix containing the predictors (n= number of observatiosn, m = number of predictors).
#' @param y The Response onto which the qrsvm shall be fitted.
#' @param kernel a string giving the type of kernels from package \link{\code{kernlab}} to use f.e.
#'   "rbfdot" for Radial Basis Function Kernel. All Kernels except "stringdot" supported.
#' @param cost The Cost parameter see f.e. package \link{\code{e1071}} and \link{\code{kernlab}}
#' @param tau The Quantile that shall be estimated. 0<=tau<=1.
#' @param sigma,degree,scale,offset,order A possible tuning parameter for specific Kernelfunctions, see \link{\code{kernlab}}.
#' @details There is no preimplemented scaling of the input variables which should be considered beforehand.
#'   Also optimization is based on "quadprog:solve.QP" function which can be considerably slow compared to
#'   other SVM implementations.
#' @return An object of class "qrsvm".
#' @references "Nonparametric Quantile Regression" by I.Takeuchi, Q.V. Le, T. Sears, A.J. Smola (2004)
#' @importFrom kernlab rbfdot polydot vanilladot tanhdot laplacedot besseldot anovadot
#' @importFrom Matrix nearPD
#' @importFrom quadprog solve.QP
#' @examples
#' # data generation
#' n <- 200
#' x <- as.matrix(seq(-1.5, 1.5, length.out = n))
#' y <- rnorm(n) * (0.3 + abs(sin(x)))
#'
#' # fit models
#' mod005 <- qrsvm(x, y, tau = 0.05)
#' mod050 <- qrsvm(x, y, tau = 0.5)
#' mod095 <- qrsvm(x, y, tau = 0.95)
#'
#' # methods
#' print(mod050)
#' summary(mod050)
#' fittedDf <- data.frame(x, fitted05 = fitted(mod005), fitted50 = fitted(mod050),
#'                        fitted95 = fitted(mod095))
#' predict(mod050, c(-1, 0, 1))
#'
#' # graph
#' library(ggplot2)
#' ggplot(data.frame(x, y), aes(x, y)) + geom_point() +
#'   geom_line(aes(x, fitted05, colour = "P05"), fittedDf) +
#'   geom_line(aes(x, fitted50, colour = "P50"), fittedDf) +
#'   geom_line(aes(x, fitted95, colour = "P95"), fittedDf) +
#'   labs(colour = expression(tau))
#' @export
qrsvm <- function(x, y, kernel = "rbfdot", cost = 1,    tau = 0.95,
                  sigma = 5, degree = 2, scale = 1, offset = 1, order = 1) {
  if (class(kernel) == "character") {
    if (kernel == "rbfdot") {
      kern <- rbfdot(sigma = sigma)
      kernmat <- kernelMat(kern, x)
    } else if (kernel == "polydot") {
      kern <- polydot(degree = degree, scale = scale, offset = offset)
      kernmat <- kernelMat(kern, x)
    } else if (kernel == "vanilladot") {
      kern <- vanilladot()
      kernmat <- kernelMat(kern, x)
    } else if (kernel == "tanhdot") {
      kern <- tanhdot(scale = scale, offset = offset)
      kernmat <- kernelMat(kern, x)
    } else if (kernel == "laplacedot") {
      kern <- laplacedot(sigma = sigma)
      kernmat <- kernelMat(kern, x)
    } else if (kernel == "besseldot") {
      kern <- besseldot(sigma = sigma, order = order, degree = degree)
      kernmat <- kernelMat(kern, x)
    } else if (kernel == "anovadot") {
      kern <- anovadot(sigma = sigma, degree = degree)
      kernmat <- kernelMat(kern, x)
    } else {
      stop("kernelMat not valid! Check if valid kernel type stated!")
    }
    pdmat <- as.matrix(nearPD(kernmat)$mat)
  } else if (class(kernel) == "matrix") {
    pdmat <- kernel
    kernmat <- kernel
    if (nrow(kernel) != ncol(kernel))
      stop("Given kernelMat has different col- rownumbers!! Check!")
  }

  Amat <- rbind(rep(1, nrow(x)), diag(x = 1, nrow(x), nrow(x)),
                diag(x = -1, nrow(x), nrow(x)))
  b0 <- c(0, rep((cost * (tau - 1)), nrow(x)), rep(-cost * tau, nrow(x)))
  erg <- solve.QP(Dmat = pdmat, dvec = y, Amat = t(Amat),
                  bvec = b0, meq = 1, factorized = FALSE)
  alpha <- erg$solution
  f <- alpha %*% kernmat
  offshift <- which.min((round(alpha, 3) - (cost * tau))^2 + (round(alpha, 3) - (cost * (tau - 1)))^2)
  b <- y[offshift] - f[offshift]
  fnew <- alpha %*% kernmat + b

  if (class(kernel) == "character") {
    model <- list(alpha = alpha, xtrain = x, kernel = kern,
                  sigma = sigma, cost = cost, b0 = b, fitted = as.numeric(fnew),
                  tau = tau, scale = scale, offset = offset, order = order,
                  kernstring = kernel, y = y)
  } else if (class(kernel) == "matrix") {
    model <- list(alpha = alpha, xtrain = x, kernel = 0,
                  sigma = sigma, cost = cost, b0 = b, fitted = as.numeric(fnew),
                  tau = tau, scale = scale, offset = offset, order = order,
                  kernstring = kernel, y = y)
  }
  class(model) <- "qrsvm"
  return(model)
}
