#' @title The Function for Constructing the Covariance Matrix in \code{GPM} Package
#' @description Builds the covariance matrix given two datasets and the type and parameters of the covariance function.
#' @param X1,X2 Matrices containing the data points. The rows and columns of both \code{X1} and \code{X2} denote individual observation settings and dimension, respectively.
#' @param CovType The covariance function of the GP model. Choices include \code{'G'} (default), \code{'PE'}, \code{'LBG'}, and \code{'LB'}. See the \code{references} for the details.
#' @param w The vector storing all the (hyper)parameters of the covariance function. The length of \code{w} depends on the \code{CovType}. See \code{reference 1}.
#'
#' @return R The covariance matrix with size \code{nrow(X1)}-by-\code{nrow(X2)}. See \href{https://en.wikipedia.org/wiki/Covariance_matrix}{here}.
#'
#' @note This function is \strong{NOT} exported once the GPM package is loaded.
#' @references
#' \enumerate{
#' \item T. Kearney, S. Tao, R. Bostanabad, D.W. Apley, and W. Chen (2017). GPM: An R Package for Gaussian Process Modeling of Deterministic Computer Simulations. The Journal of Statistical Software.
#' \item M. Plumlee, D.W. Apley (2016). Lifted Brownian kriging models, Technometrics.
#' }
#' @seealso
#' \code{\link[GPM]{Fit}} to see how a GP model can be fitted to a training dataset.\cr
#' \code{\link[GPM]{Predict}} to use the fitted GP model for prediction.\cr
#' \code{\link[GPM]{Draw}} to plot the response via the fitted model.
#' @examples
#' # see the examples in \code{\link[GPM]{Fit}}

CovMat <-  function(X1, X2, CovType, w){
k = nrow(X1)
p = ncol(X1)
m = nrow(X2)
R = matrix(0, k, m)
w = as.vector(w)
if (CovType == 'G'){
  if (p != length(w)){
    stop(paste('In a Gaussian covariance function, there should be', toString(p), 'parameters!'))
  }
  if (k >= m){
    for (i in 1: m) {
      #R[, i]=exp(-apply((X1 - t(replicate(k, X2[i, ])))**2*t(replicate(k, 10^w)), 1, sum))
      #R[ , i] = exp(-apply(t(t(t(t(X1) - X2[i, ])**2)*(10^w)), 1, sum))
      #R[i:k, i] = exp(-colSums((t(X1[i:k, , drop = FALSE] - X2[i, ]))^2*(10^w-1)))
      R[, i] = colSums((t(X1) - X2[i, ])^2*(10^w))
    }
    R = exp(-R)
  } else{
    for (i in 1: k) {
      R[i , ] = colSums((t(X2) - X1[i, ])^2*(10^w))
    }
    R = exp(-R)
  }
} else if (CovType == 'PE'){
  if (p != length(w) - 1){
    stop(paste('In a PE covariance function, there should be', toString(p+1), 'parameters!'))
  }
  if (k >= m){
    for (i in 1: m) {
      R[ , i] = exp(-colSums(abs((t(X1) - X2[i, ]))^w[p+1]*(10^w[-(p+1)])))
    }
  } else{
    for (i in 1: k) {
      R[i , ] = exp(-colSums(abs((t(X2) - X1[i, ]))^w[p+1]*(10^w[-(p+1)])))
    }
  }
} else if (CovType == 'LBG'){
  if (p != length(w) - 1){
    stop(paste('In an LBG covariance function, there should be', toString(p+1), 'parameters!'))
  }
  A = 10^w[1:p]
  Beta = w[p+1]
  if (k >= m){
    for (i in 1: m) {
      R[, i] = (1 + colSums(t(X1**2)*A))**Beta + (1 + sum(X2[i, ]**2*A))**Beta -
        (1 + colSums((t(X1) - X2[i, ])**2*A))**Beta - 1
    }
  } else{
    for (i in 1: k) {
      R[i, ] = (1 + colSums(t(X2**2)*A))**Beta + (1 + sum(X1[i, ]**2*A))**Beta -
        (1 + colSums((t(X2) - X1[i, ])**2*A))**Beta - 1
    }
  }
} else if (CovType == 'LB'){
  if (p != length(w) - 2){
    stop(paste('In an LB covariance function, there should be', toString(p+2), 'parameters!'))
  }
  A = 10^w[1:p]
  Beta = w[p+1]
  Gamma = w[p+2]
  if (k >= m){
    for (i in 1: m) {
      R[, i] = (1 + colSums(t(X1**2)*A)**Gamma)**Beta + (1 + sum(X2[i, ]**2*A)**Gamma)**Beta -
        (1 + colSums((t(X1) - X2[i, ])**2*A)**Gamma)**Beta - 1
    }
  } else{
    for (i in 1: k) {
      R[i, ] = (1 + colSums(t(X2**2)*A)**Gamma)**Beta + (1 + sum(X1[i, ]**2*A)**Gamma)**Beta -
        (1 + colSums((t(X2) - X1[i, ])**2*A)**Gamma)**Beta - 1
    }
  }
} else {
  stop('The type of the covariance is not supported!')
}
R[R < (.Machine$double.eps)] <- 0

return(R)
}
