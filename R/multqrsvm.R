#' @importFrom utils globalVariables data
utils::globalVariables(c("tt")) # avoid the variable checking

#' Fits multiple Quantile Regression SVM
#'
#' @param x An n X m matrix containing the predictors (n = number of observatiosn, m = number of predictors) .
#' @param y The Response onto which the qrsvm shall be fitted.
#' @param kernel A string giving the type of kernels from \link[kernlab]{kernelMatrix}.
#'   Default value is "rbfdot" for Radial Basis Function Kernel. All Kernels except \link[kernlab]{stringdot} supported.
#' @param cost The cost parameter see \link[e1071]{svm} and \link[kernlab]{kernelMatrix}.
#' @param tau The quantiles that shall be estimated. 0<=tau<=1.
#' @param sigma,degree,scale,offset,order A possible tuning parameter for specific Kernelfunctions,
#'   see \link[kernlab]{rbfdot}, \link[kernlab]{polydot}, \link[kernlab]{vanilladot}, \link[kernlab]{tanhdot},
#'   \link[kernlab]{laplacedot}, \link[kernlab]{besseldot} or \link[kernlab]{anovadot}.
#' @param doPar Should a parallel backend be used. Logical.
#' @param clustnum The number of parallel tasks to use given doPar==TRUE. Default = 2.
#' @details There is no preimplemented scaling of the input variables which should be considered beforehand.
#'  Also optimization is based on "quadprog:solve.QP" function which can be considerably slow compared to
#'  other SVM implementations.
#' @return An object of class "multqrsvm"
#' @references "Nonparametric Quantile Regression" by I.Takeuchi, Q.V. Le, T. Sears, A.J. Smola (2004)
#' @importFrom kernlab rbfdot polydot vanilladot tanhdot laplacedot besseldot anovadot
#' @importFrom Matrix nearPD
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom utils modifyList
#' @examples
#' # data generation
#' n <- 200
#' x <- as.matrix(seq(-2, 2, length.out = n))
#' y <- rnorm(n)*(0.3 + abs(sin(x)))
#'
#' # fit models
#' quant <- c(0.01, 0.25, 0.5, 0.75, 0.99)
#' models <- multqrsvm(x, y, tau = quant, doPar = FALSE, sigma = 1)
#'
#' # methods
#' print(models)
#' fittedDf <- data.frame(cbind(x, fitted(models)))
#' names(fittedDf) <- c("x", sprintf("fitted_%02i", quant * 100))
#' predict(models, c(-1, 0, 1))
#'
#' # graph
#' library(ggplot2)
#' g <- ggplot(data.frame(x, y), aes(x, y)) + geom_point()
#' for (i in seq_along(models)) {
#'   mapping <- aes_string("x", names(fittedDf)[i+1])
#'   mapping$colour <- sprintf("P%02i", quant[i] * 100)
#'   g <- g + geom_line(mapping, fittedDf)
#' }
#' g + labs(colour = expression(tau))
#' @export
multqrsvm <- function(x, y, kernel = "rbfdot", cost = 1, tau = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      sigma = 5, degree = 2, scale = 1, offset = 1, order = 1, doPar = FALSE, clustnum = 2){
  kernelInfo <- getKernFunc(x, kernel, sigma, degree, scale, offset, order)
  `%foreach_op%` <- `%do%`
  if (doPar) {
    message("Running parallelly with ", clustnum, " clusters.")
    registerDoParallel(clustnum)
    `%foreach_op%` <- `%dopar%`
  }
  models <- foreach(tt = tau, .packages = "qrsvm") %foreach_op% {
    model <- qrsvm(x = x, y = y, kernel = kernelInfo$pdmat, cost = cost,
                   tau = tt, sigma = sigma, degree = degree,
                   scale = scale, offset = offset, order = order)
    model <- modifyList(model, list(kernel = kernelInfo$kernel))
    return(model)
  }

  class(models) <- "multqrsvm"
  return(models)
}
