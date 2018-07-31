#' Fits multiple Quantile Regression SVM
#'
#' @param x An n X m matrix containing the predictors (n= number of observatiosn, m = number of predictors) .
#' @param y The Response onto which the qrsvm shall be fitted.
#' @param kernel a string giving the type of kernels from package \link{\code{kernlab}} to use f.e.
#'   "rbfdot" for Radial Basis Function Kernel. All Kernels except "stringdot" supported.
#' @param cost The Cost parameter see f.e. package \link{\code{e1071}} and \link{\code{kernlab}}
#' @param tau The Quantile that shall be estimated. 0<=tau<=1.
#' @param sigma,degree,scale,offset,order A possible tuning parameter for specific Kernelfunctions, see \link{\code{kernlab}}.
#' @param doPar Should a parallel backend be used. Logical.
#' @param clustnum The number of parallel tasks to use given doPar==TRUE. Default = 2.
#' @details There is no preimplemented scaling of the input variables which should be considered beforehand.
#'  Also optimization is based on "quadprog:solve.QP" function which can be considerably slow compared to
#'  other SVM implementations.
#' @return An object of class "qrsvm"
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

  `%foreach_op%` <- `%do%`
  if (doPar) {
    message("Running parallelly with ", clustnum, " clusters.")
    registerDoParallel(clustnum)
    `%foreach_op%` <- `%dopar%`
  }
  models <- foreach(tt = tau, .packages = "qrsvm") %foreach_op% {
    model <- qrsvm(x = x, y = y, kernel = pdmat, cost = cost,
                   tau = tt, sigma = sigma, degree = degree,
                   scale = scale, offset = offset, order = order)
    modifyList(model, list(kernel = kern))
    return(model)
  }

  class(models) <- "multqrsvm"
  return(models)
}
