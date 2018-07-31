#' @method print multqrsvm
#' @export
print.multqrsvm <- function(x, ...){
  for (i in seq_along(x)) {
    print(x[[i]])
  }
}

#' Fitted values of class "multqrsvm"
#'
#' @param object An object of class "multqrsvm"
#' @param ... other arguments.
#' @return A numeric vector of predicted values
#' @method fitted multqrsvm
#' @importFrom stats fitted
#' @export
fitted.multqrsvm <- function(object, ...) {
  return(sapply(object, fitted))
}

#' Predict an Object of class "multqrsvm"
#'
#' @param object An object of class "multqrsvm".
#' @param newdata The predictors of the predictable data in an n X m Matrix.
#' @param ... other arguments.
#' @return A list of predicted values
#' @method predict multqrsvm
#' @importFrom kernlab kernelMult
#' @export
predict.multqrsvm <- function(object, newdata, ...) {
  prediction <- vector("list", length(object))
  if (ncol(newdata) != ncol(object[[1]]$xtrain)) {
    stop("Newdata has different number of columns than xtrain please check consistency!")
  }
  prediction <- lapply(object, function(o){
    kernelMult(o$kernel, newdata, o$xtrain, o$alpha) + o$b0
  })
  return(prediction)
}
