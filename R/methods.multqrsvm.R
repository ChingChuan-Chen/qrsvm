#' @method print multqrsvm
#' @export
print.multqrsvm <- function(x, ...){
  for (i in seq_along(x)) {
    print(x[[i]])
  }
}

#' Fitted values of class "multqrsvm"
#'
#' @param object An object of class \link{multqrsvm}.
#' @param ... other arguments.
#' @return A numeric matrix of fitted values.
#' @method fitted multqrsvm
#' @importFrom stats fitted
#' @export
fitted.multqrsvm <- function(object, ...) {
  return(sapply(object, fitted))
}

#' Predict an Object of class "multqrsvm"
#'
#' @param object An object of class \link{multqrsvm}.
#' @param newdata The predictors of the predictable data in an n X m Matrix.
#' @param ... other arguments.
#' @return A list of predicted values.
#' @method predict multqrsvm
#' @importFrom stats predict
#' @importFrom kernlab kernelMult
#' @export
predict.multqrsvm <- function(object, newdata, ...) {
  if (!is.matrix(newdata))
    newdata <- as.matrix(newdata)
  if (ncol(newdata) != ncol(object[[1]]$xtrain))
    stop("Newdata has different number of columns than xtrain please check consistency!")
  return(sapply(object, predict, newdata = newdata))
}
