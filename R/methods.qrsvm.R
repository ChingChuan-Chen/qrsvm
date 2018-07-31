#' @method print qrsvm
#' @export
print.qrsvm <- function(x, ...) {
  cat(paste0("Quantile Regression SVM. Cost is ", x$cost), fill = TRUE)
  cat(paste0("Estimated Quantile is ", x$tau), fill = TRUE)
}

#' @method summary qrsvm
#' @export
summary.qrsvm <- function(object, ...) {
  idx <- c(4L:6L, 8L:12L)
  mapply(function(n, x){
    cat(paste(n, " is ", x), fill = TRUE)
  }, names(object[idx]), object[idx])
  invisible(NULL)
}

#' Fitted values of class "qrsvm"
#'
#' @param object An object of class \link{qrsvm}.
#' @param ... other arguments.
#' @return A numeric vector of fitted values.
#' @method fitted qrsvm
#' @importFrom stats fitted
#' @export
fitted.qrsvm <- function(object, ...) {
  return(object$fitted)
}

#' Predict an Object of class "qrsvm"
#'
#' @param object An object of class \link{qrsvm}.
#' @param newdata The predictors of the predictable data in an n X m Matrix.
#' @param ... other arguments.
#' @return A numeric vector of predicted values.
#' @method predict qrsvm
#' @importFrom stats predict
#' @importFrom kernlab kernelMult
#' @export
predict.qrsvm <- function(object, newdata, ...) {
  if (!is.matrix(newdata))
    newdata <- as.matrix(newdata)
  xold <- object$xtrain
  if (ncol(newdata) != ncol(xold))
    stop("Newdata has different number of columns than xtrain please check consistency!")
  return(kernelMult(object$kernel, newdata, xold, object$alpha) + object$b0)
}
